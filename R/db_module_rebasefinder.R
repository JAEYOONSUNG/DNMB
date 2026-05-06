.dnmb_rebasefinder_module_name <- function() {
  "rebasefinder"
}

.dnmb_rebasefinder_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_rebasefinder_empty_status <- function() {
  .dnmb_rebasefinder_status_row(character(), character(), character())
}

.dnmb_rebasefinder_cache_dir <- function(cache_root = NULL, create = TRUE) {
  root <- .dnmb_db_cache_root(cache_root = cache_root, create = create)
  path <- file.path(root, "rebasefinder", "cache")
  if (isTRUE(create)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

# ── REBASE bairoch metadata: methylation type + position lookup ──
# Downloads once from rebase.neb.com/rebase/link_bairoch, parses the
# bairoch-format flat file, and builds a lookup table keyed by enzyme
# name. Each entry includes the recognition sequence, methylation
# type(s) (m6A / m5C / Nm4C), and position(s) within the rec_seq.
#
# This is the ONLY authoritative source for which base in the
# recognition sequence is methylated and how — motif-based inference
# is a heuristic fallback at best.

.dnmb_rebasefinder_bairoch_cache_path <- function(cache_root = NULL) {
  file.path(.dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE), "rebase_bairoch_lookup.rds")
}

.dnmb_rebasefinder_download_bairoch <- function(cache_root = NULL, force = FALSE) {
  cache_path <- .dnmb_rebasefinder_bairoch_cache_path(cache_root)
  if (!isTRUE(force) && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  legacy_path <- file.path(
    .dnmb_db_cache_root(cache_root = cache_root, create = FALSE),
    "db_modules", "rebasefinder", "cache", "rebase_bairoch_lookup.rds"
  )
  if (!isTRUE(force) && file.exists(legacy_path)) {
    lookup <- readRDS(legacy_path)
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    tryCatch(saveRDS(lookup, cache_path), error = function(e) NULL)
    return(lookup)
  }
  url <- "http://rebase.neb.com/rebase/link_bairoch"
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  tryCatch(utils::download.file(url, tmp, quiet = TRUE, method = "libcurl"),
           error = function(e) NULL)
  if (!file.exists(tmp) || file.info(tmp)$size < 1000) {
    return(NULL)
  }
  lookup <- .dnmb_rebasefinder_parse_bairoch(tmp)
  if (!is.null(lookup) && nrow(lookup) > 0) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(lookup, cache_path)
  }
  lookup
}

.dnmb_rebasefinder_parse_bairoch <- function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- list()
  current <- list(id = NA_character_, rs = NA_character_, ms = NA_character_)
  for (line in lines) {
    if (startsWith(line, "ID   ")) {
      current$id <- trimws(sub("^ID   ", "", line))
    } else if (startsWith(line, "RS   ")) {
      rs_raw <- trimws(sub("^RS   ", "", line))
      rs_parts <- strsplit(rs_raw, ",")[[1]]
      current$rs <- trimws(gsub(";", "", rs_parts[1]))
    } else if (startsWith(line, "MS   ")) {
      current$ms <- trimws(sub("^MS   ", "", line))
    } else if (line == "//") {
      if (!is.na(current$id) && !is.na(current$ms)) {
        entries[[length(entries) + 1L]] <- current
      }
      current <- list(id = NA_character_, rs = NA_character_, ms = NA_character_)
    }
  }
  if (!length(entries)) return(NULL)
  # Parse MS field: "2(m6A),-4(m6A)" → list of {pos, type}
  rows <- lapply(entries, function(e) {
    ms_parts <- strsplit(gsub(";$", "", e$ms), ",")[[1]]
    ms_parsed <- lapply(trimws(ms_parts), function(p) {
      m <- regmatches(p, regexec("^(-?[0-9?]+)\\(([^)]+)\\)", p))[[1]]
      if (length(m) == 3L) list(pos = m[2], type = m[3]) else NULL
    })
    ms_parsed <- Filter(Negate(is.null), ms_parsed)
    # Primary methylation: first entry
    primary_pos <- if (length(ms_parsed)) ms_parsed[[1]]$pos else NA_character_
    primary_type <- if (length(ms_parsed)) ms_parsed[[1]]$type else NA_character_
    data.frame(
      enzyme_name = e$id,
      rec_seq = if (is.na(e$rs)) NA_character_ else e$rs,
      meth_type = primary_type,
      meth_pos = primary_pos,
      meth_all = e$ms,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.dnmb_rebasefinder_lookup_methylation <- function(enzyme_name, cache_root = NULL) {
  lookup <- .dnmb_rebasefinder_download_bairoch(cache_root = cache_root)
  if (is.null(lookup) || !nrow(lookup)) {
    return(list(meth_type = NA_character_, meth_pos = NA_character_,
                rec_seq = NA_character_, meth_all = NA_character_))
  }
  # Match by enzyme name (exact or prototype)
  idx <- match(enzyme_name, lookup$enzyme_name)
  if (is.na(idx)) {
    return(list(meth_type = NA_character_, meth_pos = NA_character_,
                rec_seq = NA_character_, meth_all = NA_character_))
  }
  as.list(lookup[idx, ])
}

.dnmb_rebasefinder_role_from_hit <- function(hit_name) {
  hit_name <- base::as.character(hit_name)
  base::ifelse(base::grepl("^M[0-9]?\\.", hit_name), "M",
    base::ifelse(base::grepl("^R[0-9]?\\.", hit_name), "R",
      base::ifelse(base::grepl("^S[0-9]?\\.", hit_name), "S", NA_character_)
    )
  )
}

.dnmb_rebasefinder_clean_text <- function(...) {
  parts <- base::lapply(list(...), base::as.character)
  parts <- base::unlist(parts, use.names = FALSE)
  parts <- parts[!base::is.na(parts) & base::nzchar(parts)]
  if (!base::length(parts)) return(NA_character_)
  base::tolower(base::paste(parts, collapse = " | "))
}

.dnmb_rebasefinder_gene_annotation_text <- function(genes) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  text_cols <- base::grep(
    "^(product|gene|note)$|description|accession|pfam|interpro|cdd|gene3d|superfamily|ncbi|smart|prosite|eggNOG|PADLOC|DefenseFinder",
    base::names(genes),
    ignore.case = TRUE,
    value = TRUE
  )
  if (!base::length(text_cols)) return(base::rep(NA_character_, base::nrow(genes)))
  base::vapply(base::seq_len(base::nrow(genes)), function(i) {
    .dnmb_rebasefinder_clean_text(genes[i, text_cols, drop = TRUE])
  }, character(1))
}

.dnmb_rebasefinder_partial_threshold <- function(family_id, enzyme_role) {
  family_id <- base::as.character(family_id)[1]
  enzyme_role <- base::as.character(enzyme_role)[1]
  if (base::is.na(family_id)) family_id <- ""
  if (base::is.na(enzyme_role)) enzyme_role <- ""

  if (identical(family_id, "Type I") && identical(enzyme_role, "R")) return(700L)
  if (identical(family_id, "Type I") && identical(enzyme_role, "M")) return(350L)
  if (identical(family_id, "Type I") && identical(enzyme_role, "S")) return(300L)
  if (identical(family_id, "Type III") && identical(enzyme_role, "R")) return(650L)
  if (identical(family_id, "Type III") && identical(enzyme_role, "M")) return(300L)
  if (identical(family_id, "Type IV") && identical(enzyme_role, "R")) return(250L)
  if (identical(enzyme_role, "M")) return(200L)
  if (identical(enzyme_role, "R")) return(150L)
  if (identical(enzyme_role, "S")) return(250L)
  100L
}

.dnmb_rebasefinder_sequence_partial_table <- function(tbl,
                                                       family_col = "REBASEfinder_family_id",
                                                       role_col = "REBASEfinder_enzyme_role") {
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  n <- base::nrow(tbl)
  seq <- if ("translation" %in% base::names(tbl)) {
    .dnmb_normalize_translation(tbl$translation)
  } else {
    base::rep(NA_character_, n)
  }
  aa_len <- base::nchar(seq)
  aa_len[base::is.na(seq)] <- NA_integer_
  family <- if (family_col %in% base::names(tbl)) base::as.character(tbl[[family_col]]) else base::rep(NA_character_, n)
  role <- if (role_col %in% base::names(tbl)) base::as.character(tbl[[role_col]]) else base::rep(NA_character_, n)
  expected_min <- base::mapply(.dnmb_rebasefinder_partial_threshold, family, role, USE.NAMES = FALSE)

  text_cols <- base::intersect(
    c("product", "note", "support", "REBASEfinder_support"),
    base::names(tbl)
  )
  anno <- if (base::length(text_cols)) {
    base::vapply(base::seq_len(n), function(i) {
      base::tolower(base::paste(stats::na.omit(base::as.character(tbl[i, text_cols, drop = TRUE])), collapse = " "))
    }, character(1))
  } else {
    base::rep("", n)
  }
  annotated_partial <- base::grepl(
    "partial|truncat|fragment|incomplete|missing[ -]?(n|c)[ -]?terminus|frameshift|pseudogene",
    anno,
    perl = TRUE
  )
  short <- !base::is.na(aa_len) & aa_len < expected_min
  status <- base::ifelse(short | annotated_partial, "partial_or_short", "full_length")
  reason <- base::vapply(base::seq_len(n), function(i) {
    parts <- c(
      if (short[[i]]) base::sprintf("aa_len=%s below expected_min=%s", aa_len[[i]], expected_min[[i]]) else NULL,
      if (annotated_partial[[i]]) "annotation suggests partial/truncated sequence" else NULL
    )
    if (!base::length(parts)) return(NA_character_)
    base::paste(parts, collapse = "; ")
  }, character(1))

  base::data.frame(
    aa_len = aa_len,
    expected_min_aa = expected_min,
    partial_status = status,
    partial_reason = reason,
    stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_typeiii_context_role <- function(text, known_rm_type = NA_character_, known_role = NA_character_) {
  text <- base::tolower(base::as.character(text))
  known_rm_type <- base::as.character(known_rm_type)
  known_role <- base::as.character(known_role)
  if (!base::is.na(known_rm_type) && identical(known_rm_type, "Type III") &&
      !base::is.na(known_role) && known_role %in% c("M", "R", "S")) {
    return(known_role)
  }
  if (base::is.na(text) || !base::nzchar(text)) return(NA_character_)

  typeiii_token <- base::grepl("type[ _-]?iii|\\biii\\b|res[ _-]?(iii|3)|mod[ _-]?(iii|3)", text, perl = TRUE)
  mod_token <- base::grepl(
    "type[ _-]?iii.*(mod|methyl|methyltransferase)|\\bmod\\b|modification methylase|site-specific dna-methyltransferase|dna methyltransferase|sam-dependent methyltransferase|methyltransferase type 11",
    text,
    perl = TRUE
  )
  res_token <- base::grepl(
    "type[ _-]?iii.*(res|restriction)|\\bres\\b|res[ _-]?(iii|3)|restriction endonuclease|dead/deah|dexh|superfamily 2 helicase|helicase|atp-dependent|atpase|pd-\\(?d/e\\)?xk|pdexk|pd-dexk",
    text,
    perl = TRUE
  )

  if (typeiii_token && mod_token) return("M")
  if (typeiii_token && res_token) return("R")
  if (base::grepl("specificity subunit|hsds|target recognition", text, perl = TRUE)) return("S")
  if (base::grepl("res subunit|restriction endonuclease subunit r|type[ _-]?iii restriction enzyme", text, perl = TRUE)) return("R")
  if (base::grepl("mod subunit|type[ _-]?iii methyltransferase", text, perl = TRUE)) return("M")
  NA_character_
}

.dnmb_rebasefinder_typeiii_sequence_signals <- function(translation) {
  seq <- base::toupper(base::gsub("[^A-Za-z]", "", base::as.character(translation)))
  if (base::is.na(seq) || !base::nzchar(seq)) {
    return(list(
      role = NA_character_,
      signals = character(),
      has_sam = FALSE,
      has_mtase = FALSE,
      res_signal_count = 0L
    ))
  }
  hits <- c(
    SAM = base::grepl("[FYW].G.[GA]", seq, perl = TRUE),
    MTase = base::grepl("[DNS]PP[YFW]", seq, perl = TRUE),
    `Helicase-DEAD` = base::grepl("DE[AHCFQ][DHQ]", seq, perl = TRUE),
    `ResIII-WA` = base::grepl("[AG].{4}GK[ST]", seq, perl = TRUE),
    `ResIII-WB` = base::grepl("DE.H", seq, perl = TRUE),
    `ResIII-PD` = base::grepl("PD.{1,25}[DE].K", seq, perl = TRUE)
  )
  res_signal_count <- base::sum(hits[c("ResIII-WA", "ResIII-WB", "ResIII-PD")])
  role <- if (res_signal_count >= 2L) {
    "R"
  } else if (hits[["Helicase-DEAD"]] && base::nchar(seq) >= 600L && !hits[["MTase"]]) {
    "R"
  } else if (hits[["SAM"]] && hits[["MTase"]]) {
    "M"
  } else {
    NA_character_
  }
  list(
    role = role,
    signals = base::names(hits)[hits],
    has_sam = hits[["SAM"]],
    has_mtase = hits[["MTase"]],
    has_helicase = hits[["Helicase-DEAD"]],
    res_signal_count = res_signal_count
  )
}

.dnmb_rebasefinder_typeiii_sequence_role <- function(translation) {
  .dnmb_rebasefinder_typeiii_sequence_signals(translation)$role
}

.dnmb_rebasefinder_strip_typeiii_support <- function(x) {
  base::vapply(base::as.character(x), function(one) {
    if (base::is.na(one) || !base::nzchar(one)) return(NA_character_)
    parts <- base::trimws(base::strsplit(one, ";", fixed = TRUE)[[1]])
    parts <- parts[base::nzchar(parts)]
    parts <- parts[!base::grepl("^(typeIII_context=|partners=)", parts)]
    if (!base::length(parts)) NA_character_ else base::paste(parts, collapse = "; ")
  }, character(1))
}

.dnmb_rebasefinder_typei_context_role <- function(text, known_rm_type = NA_character_, known_role = NA_character_) {
  text <- base::tolower(base::as.character(text))
  known_rm_type <- base::as.character(known_rm_type)
  known_role <- base::as.character(known_role)
  if (!base::is.na(known_rm_type) && identical(known_rm_type, "Type I") &&
      !base::is.na(known_role) && known_role %in% c("M", "R", "S")) {
    return(known_role)
  }
  if (base::is.na(text) || !base::nzchar(text)) return(NA_character_)

  if (base::grepl("hsds|specificity subunit|target recognition|restriction endonuclease subunit s|\\bsubunit s\\b", text, perl = TRUE)) {
    return("S")
  }
  if (base::grepl("hsdr|type[ _-]?i restriction-modification system endonuclease|restriction-modification system endonuclease|restriction endonuclease subunit r|\\bsubunit r\\b", text, perl = TRUE)) {
    return("R")
  }
  if (base::grepl("dead/deah.*helicase|helicase family protein|superfamily 2 helicase|atp-dependent helicase", text, perl = TRUE)) {
    return("R")
  }
  if (base::grepl("hsdm|class[ _-]?i .*methyl|type[ _-]?i .*methyl|n-6 dna methylase|sam-dependent dna methyltransferase|dna methyltransferase", text, perl = TRUE)) {
    return("M")
  }
  NA_character_
}

.dnmb_rebasefinder_typei_sequence_signals <- function(translation) {
  seq <- base::toupper(base::gsub("[^A-Za-z]", "", base::as.character(translation)))
  if (base::is.na(seq) || !base::nzchar(seq)) {
    return(list(role = NA_character_, signals = character()))
  }
  hits <- c(
    SAM = base::grepl("[FYW].G.[GA]", seq, perl = TRUE),
    MTase = base::grepl("[DNS]PP[YFW]", seq, perl = TRUE),
    Ploop = base::grepl("[AG].{4}GK[ST]", seq, perl = TRUE),
    DEAD = base::grepl("DE[AHCFQ][DHQ]", seq, perl = TRUE),
    PD = base::grepl("PD.{1,25}[DE].K", seq, perl = TRUE)
  )
  role <- if ((hits[["Ploop"]] && hits[["DEAD"]]) || (hits[["PD"]] && hits[["DEAD"]])) {
    "R"
  } else if (hits[["MTase"]]) {
    "M"
  } else {
    NA_character_
  }
  list(role = role, signals = base::names(hits)[hits])
}

.dnmb_rebasefinder_strip_typei_support <- function(x) {
  base::vapply(base::as.character(x), function(one) {
    if (base::is.na(one) || !base::nzchar(one)) return(NA_character_)
    parts <- base::trimws(base::strsplit(one, ";", fixed = TRUE)[[1]])
    parts <- parts[base::nzchar(parts)]
    parts <- parts[!base::grepl("^(typeI_context=|partners=)", parts)]
    if (!base::length(parts)) NA_character_ else base::paste(parts, collapse = "; ")
  }, character(1))
}

.dnmb_rebasefinder_gene_order <- function(genes) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  for (col in c("start", "end")) {
    if (col %in% base::names(genes)) genes[[col]] <- suppressWarnings(base::as.numeric(genes[[col]]))
  }
  if (!"contig" %in% base::names(genes)) genes$contig <- "contig"
  if (!"direction" %in% base::names(genes)) genes$direction <- NA_character_
  genes$.dnmb_row <- base::seq_len(base::nrow(genes))
  genes <- genes[!base::is.na(genes$locus_tag) & base::nzchar(genes$locus_tag), , drop = FALSE]
  genes <- genes[base::order(genes$contig, genes$start, genes$end, genes$.dnmb_row), , drop = FALSE]
  genes$.dnmb_order <- base::seq_len(base::nrow(genes))
  rownames(genes) <- NULL
  genes
}

.dnmb_rebasefinder_intergenic_gap <- function(a, b) {
  if (base::is.na(a$start) || base::is.na(a$end) || base::is.na(b$start) || base::is.na(b$end)) return(NA_real_)
  base::max(0, base::max(a$start, b$start) - base::min(a$end, b$end))
}

.dnmb_rebasefinder_structure_path <- function(output_dir, structure_validation = NULL) {
  candidates <- c(
    structure_validation,
    file.path(output_dir, "foldseek_results.tsv"),
    file.path(output_dir, "foldseek_validation.tsv"),
    file.path(output_dir, "structure_validation.tsv"),
    file.path(output_dir, "dnmb_module_rebasefinder", "foldseek_results.tsv"),
    file.path(output_dir, "dnmb_module_rebasefinder", "structure_validation.tsv")
  )
  candidates <- candidates[!base::is.na(candidates) & base::nzchar(candidates)]
  candidates <- candidates[base::file.exists(candidates)]
  if (base::length(candidates)) base::normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE) else NULL
}

.dnmb_rebasefinder_structure_reference_manifest <- function(path = NULL) {
  candidates <- c(
    path,
    system.file("extdata", "rebasefinder_structure_refs.tsv", package = "DNMB"),
    file.path("inst", "extdata", "rebasefinder_structure_refs.tsv")
  )
  candidates <- candidates[!base::is.na(candidates) & base::nzchar(candidates) & base::file.exists(candidates)]
  if (!base::length(candidates)) return(NULL)
  refs <- tryCatch(utils::read.delim(candidates[[1]], stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) NULL)
  if (base::is.null(refs) || !base::nrow(refs) ||
      !all(c("reference_id", "pdb_id", "rm_type", "enzyme_role") %in% base::names(refs))) {
    return(NULL)
  }
  refs$reference_id <- base::as.character(refs$reference_id)
  refs$pdb_id <- base::toupper(base::as.character(refs$pdb_id))
  refs
}

.dnmb_rebasefinder_structure_reference_id <- function(target) {
  target <- base::as.character(target)
  out <- base::sub("__.*$", "", target)
  out <- base::sub("\\.(cif|pdb|mmcif)$", "", out, ignore.case = TRUE)
  out[!base::nzchar(out)] <- NA_character_
  out
}

.dnmb_rebasefinder_structure_chain_id <- function(target) {
  target <- base::basename(base::as.character(target))
  has_chain <- base::grepl("__[A-Za-z0-9]{4}_[^_]+$", target)
  out <- base::rep(NA_character_, base::length(target))
  out[has_chain] <- base::sub("^.*__[A-Za-z0-9]{4}_([^_]+)$", "\\1", target[has_chain])
  out <- base::sub("\\.(cif|pdb|mmcif)$", "", out, ignore.case = TRUE)
  out
}

.dnmb_rebasefinder_chain_role <- function(chain, chain_role_map) {
  chain <- base::as.character(chain)
  chain_role_map <- base::as.character(chain_role_map)
  base::vapply(base::seq_along(chain), function(i) {
    if (base::is.na(chain[i]) || !base::nzchar(chain[i]) ||
        base::is.na(chain_role_map[i]) || !base::nzchar(chain_role_map[i])) {
      return(NA_character_)
    }
    parts <- base::strsplit(chain_role_map[i], ",", fixed = TRUE)[[1]]
    kv <- base::strsplit(parts, ":", fixed = TRUE)
    keys <- base::vapply(kv, function(x) if (base::length(x) >= 1L) base::trimws(x[[1]]) else NA_character_, character(1))
    vals <- base::vapply(kv, function(x) if (base::length(x) >= 2L) base::trimws(x[[2]]) else NA_character_, character(1))
    idx <- base::match(chain[i], keys)
    if (base::is.na(idx)) NA_character_ else vals[[idx]]
  }, character(1))
}

.dnmb_rebasefinder_read_structure_validation <- function(path,
                                                         max_evalue = 1e-3,
                                                         min_probability = 0.50,
                                                         min_tmscore = 0.45,
                                                         reference_manifest = NULL) {
  if (base::is.null(path) || !base::file.exists(path)) return(NULL)
  lines <- base::readLines(path, warn = FALSE)
  lines <- lines[base::nzchar(lines) & !base::grepl("^#", lines)]
  if (!base::length(lines)) return(NULL)
  first <- base::strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
  first_l <- base::tolower(first)
  has_header <- base::any(first_l %in% c("query", "qseqid", "locus_tag", "target", "tseqid", "sseqid", "prob", "probability", "evalue", "bits", "bitscore", "alntmscore", "tmscore"))
  tbl <- utils::read.delim(
    text = base::paste(lines, collapse = "\n"),
    header = has_header,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!has_header) {
    nms <- c("query", "target", "fident", "alnlen", "mismatch", "gapopen",
             "qstart", "qend", "tstart", "tend", "evalue", "bits")
    base::names(tbl) <- nms[seq_len(base::ncol(tbl))]
  }
  names_l <- base::tolower(base::names(tbl))
  pick <- function(candidates) {
    idx <- base::match(candidates, names_l)
    idx <- idx[!base::is.na(idx)]
    if (base::length(idx)) base::names(tbl)[idx[[1]]] else NA_character_
  }
  q_col <- pick(c("query", "qseqid", "locus_tag", "protein", "protein_id"))
  t_col <- pick(c("target", "tseqid", "sseqid", "subject", "hit", "template"))
  if (base::is.na(q_col) || base::is.na(t_col)) return(NULL)
  e_col <- pick(c("evalue", "eval", "expect"))
  b_col <- pick(c("bits", "bitscore", "score"))
  p_col <- pick(c("prob", "probability"))
  tm_col <- pick(c("alntmscore", "tmscore", "tm_score", "qtm", "ttm", "qtmscore", "ttmscore"))
  rmsd_col <- pick(c("rmsd"))

  out <- data.frame(
    query = .dnmb_module_clean_annotation_key(tbl[[q_col]]),
    structure_hit = base::as.character(tbl[[t_col]]),
    structure_evalue = if (!base::is.na(e_col)) suppressWarnings(base::as.numeric(tbl[[e_col]])) else NA_real_,
    structure_bitscore = if (!base::is.na(b_col)) suppressWarnings(base::as.numeric(tbl[[b_col]])) else NA_real_,
    structure_probability = if (!base::is.na(p_col)) suppressWarnings(base::as.numeric(tbl[[p_col]])) else NA_real_,
    structure_tmscore = if (!base::is.na(tm_col)) suppressWarnings(base::as.numeric(tbl[[tm_col]])) else NA_real_,
    structure_rmsd = if (!base::is.na(rmsd_col)) suppressWarnings(base::as.numeric(tbl[[rmsd_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out[!base::is.na(out$query) & base::nzchar(out$query), , drop = FALSE]
  if (!base::nrow(out)) return(NULL)
  out$structure_reference_id <- .dnmb_rebasefinder_structure_reference_id(out$structure_hit)
  out$structure_chain <- .dnmb_rebasefinder_structure_chain_id(out$structure_hit)
  refs <- .dnmb_rebasefinder_structure_reference_manifest(reference_manifest)
  if (!base::is.null(refs)) {
    ref_idx <- base::match(out$structure_reference_id, refs$reference_id)
    out$structure_family <- refs$rm_type[ref_idx]
    out$structure_role <- refs$enzyme_role[ref_idx]
    out$structure_class <- refs$structural_class[ref_idx]
    out$structure_chain_role <- if ("chain_role_map" %in% base::names(refs)) {
      .dnmb_rebasefinder_chain_role(out$structure_chain, refs$chain_role_map[ref_idx])
    } else {
      NA_character_
    }
  } else {
    out$structure_family <- NA_character_
    out$structure_role <- NA_character_
    out$structure_class <- NA_character_
    out$structure_chain_role <- NA_character_
  }
  missing_family <- base::is.na(out$structure_family) | !base::nzchar(out$structure_family)
  out$structure_family <- base::ifelse(base::grepl("type[ _-]?i(\\b|[^iv])", out$structure_hit, ignore.case = TRUE), "Type I",
    base::ifelse(base::grepl("type[ _-]?ii(\\b|[^i])", out$structure_hit, ignore.case = TRUE), "Type II",
      base::ifelse(base::grepl("type[ _-]?iii|resiii|modiii", out$structure_hit, ignore.case = TRUE), "Type III",
        base::ifelse(base::grepl("type[ _-]?iv", out$structure_hit, ignore.case = TRUE), "Type IV", NA_character_)
      )
    )
  )
  out$structure_family[!missing_family] <- refs$rm_type[base::match(out$structure_reference_id[!missing_family], refs$reference_id)]
  has_chain_role <- !base::is.na(out$structure_chain_role) & base::nzchar(out$structure_chain_role)
  out$structure_role[has_chain_role] <- out$structure_chain_role[has_chain_role]
  inferred_role <- .dnmb_rebasefinder_role_from_hit(out$structure_hit)
  missing_role <- base::is.na(out$structure_role) | !base::nzchar(out$structure_role)
  out$structure_role[missing_role] <- inferred_role[missing_role]
  missing_role <- base::is.na(out$structure_role) | !base::nzchar(out$structure_role)
  out$structure_role[missing_role] <- base::vapply(out$structure_hit[missing_role], function(x) {
    .dnmb_rebasefinder_typeiii_context_role(x)
  }, character(1))
  excluded_chain <- base::tolower(out$structure_role) %in% c("exclude", "excluded", "non_rm", "ocr", "anti_restriction")
  out$structure_pass <- !excluded_chain & (
    (!base::is.na(out$structure_evalue) & out$structure_evalue <= max_evalue) |
      (!base::is.na(out$structure_probability) & out$structure_probability >= min_probability) |
      (!base::is.na(out$structure_tmscore) & out$structure_tmscore >= min_tmscore)
  )
  out$structure_status <- base::ifelse(excluded_chain, "structure_excluded",
                                base::ifelse(out$structure_pass, "structure_supported", "structure_weak"))
  out <- out[base::order(
    out$query,
    !out$structure_pass,
    base::ifelse(base::is.na(out$structure_evalue), Inf, out$structure_evalue),
    -base::ifelse(base::is.na(out$structure_probability), -Inf, out$structure_probability),
    -base::ifelse(base::is.na(out$structure_tmscore), -Inf, out$structure_tmscore),
    -base::ifelse(base::is.na(out$structure_bitscore), -Inf, out$structure_bitscore)
  ), , drop = FALSE]
  out <- out[!base::duplicated(out$query), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_rebasefinder_structure_support <- function(row) {
  parts <- c(
    if (!base::is.na(row$structure_hit) && base::nzchar(row$structure_hit)) base::paste0("fold=", row$structure_hit) else NULL,
    if (!base::is.na(row$structure_probability)) base::paste0("prob=", base::round(row$structure_probability, 3)) else NULL,
    if (!base::is.na(row$structure_tmscore)) base::paste0("TM=", base::round(row$structure_tmscore, 3)) else NULL,
    if (!base::is.na(row$structure_evalue)) base::paste0("fold_e=", format(row$structure_evalue, scientific = TRUE, digits = 3)) else NULL
  )
  if (!base::length(parts)) NA_character_ else base::paste(parts, collapse = "; ")
}

.dnmb_rebasefinder_structure_family_compatible <- function(candidate_family, structure_family) {
  candidate_family <- base::as.character(candidate_family)
  structure_family <- base::as.character(structure_family)
  base::vapply(base::seq_along(candidate_family), function(i) {
    cf <- candidate_family[[i]]
    sf <- structure_family[[i]]
    if (base::is.na(cf) || !base::nzchar(cf) || base::is.na(sf) || !base::nzchar(sf)) return(TRUE)
    if (identical(cf, sf)) return(TRUE)
    # Treat Type IIS/IIG/IIP references as Type II-compatible.
    if (identical(cf, "Type II") && base::grepl("^Type II", sf)) return(TRUE)
    if (identical(sf, "Type II") && base::grepl("^Type II", cf)) return(TRUE)
    FALSE
  }, logical(1))
}

.dnmb_rebasefinder_structure_role_compatible <- function(candidate_role, structure_role, structure_chain_role = NA_character_) {
  candidate_role <- base::as.character(candidate_role)
  structure_role <- base::as.character(structure_role)
  structure_chain_role <- base::as.character(structure_chain_role)
  base::vapply(base::seq_along(candidate_role), function(i) {
    cr <- candidate_role[[i]]
    if (base::is.na(cr) || !base::nzchar(cr)) return(TRUE)
	    roles <- stats::na.omit(c(structure_chain_role[[i]], structure_role[[i]]))
	    roles <- roles[base::nzchar(roles)]
	    if (!base::length(roles)) return(TRUE)
	    role_tokens <- base::unlist(base::strsplit(roles, "[/,+;[:space:]]+", perl = TRUE), use.names = FALSE)
	    role_tokens <- role_tokens[base::nzchar(role_tokens)]
	    cr %in% role_tokens
	  }, logical(1))
}

.dnmb_rebasefinder_structure_candidate_consistency <- function(candidate_family,
                                                              candidate_role,
                                                              structure_family,
                                                              structure_role,
                                                              structure_chain_role = NA_character_) {
  family_ok <- .dnmb_rebasefinder_structure_family_compatible(candidate_family, structure_family)
  role_ok <- .dnmb_rebasefinder_structure_role_compatible(candidate_role, structure_role, structure_chain_role)
  base::data.frame(
    structure_family_consistent = family_ok,
    structure_role_consistent = role_ok,
    structure_candidate_consistent = family_ok & role_ok,
    stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_merge_structure_validation <- function(hits, genes, structure_tbl) {
  if (base::is.null(structure_tbl) || !base::is.data.frame(structure_tbl) || !base::nrow(structure_tbl)) return(hits)
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)

  extra_cols <- base::setdiff(base::names(structure_tbl), "query")
  for (col in extra_cols) {
    if (!col %in% base::names(hits)) hits[[col]] <- .dnmb_na_vector_like(structure_tbl[[col]], base::nrow(hits))
  }
  idx <- base::match(hits$query, structure_tbl$query)
  include <- !base::is.na(idx)
  for (col in extra_cols) hits[[col]][include] <- structure_tbl[[col]][idx[include]]

  if ("structure_pass" %in% base::names(hits)) {
    cons <- .dnmb_rebasefinder_structure_candidate_consistency(
      hits$family_id,
      hits$enzyme_role,
      if ("structure_family" %in% base::names(hits)) hits$structure_family else NA_character_,
      if ("structure_role" %in% base::names(hits)) hits$structure_role else NA_character_,
      if ("structure_chain_role" %in% base::names(hits)) hits$structure_chain_role else NA_character_
    )
    hits$structure_family_consistent <- cons$structure_family_consistent
    hits$structure_role_consistent <- cons$structure_role_consistent
    hits$structure_candidate_consistent <- cons$structure_candidate_consistent
    raw_pass <- !base::is.na(hits$structure_pass) & hits$structure_pass
    pass <- raw_pass & hits$structure_candidate_consistent
    mismatch <- raw_pass & !hits$structure_candidate_consistent
    if (base::any(mismatch) && "structure_status" %in% base::names(hits)) {
      hits$structure_status[mismatch & !hits$structure_family_consistent] <- "structure_family_mismatch"
      hits$structure_status[mismatch & hits$structure_family_consistent & !hits$structure_role_consistent] <- "structure_role_mismatch"
    }
    hits$evidence_mode[pass & hits$evidence_mode %in% c("low_confidence", "annotation_only")] <- "structure_supported"
    hits$typing_eligible[pass] <- TRUE
    missing_family <- pass & (base::is.na(hits$family_id) | !base::nzchar(hits$family_id)) &
      !base::is.na(hits$structure_family) & base::nzchar(hits$structure_family)
    hits$family_id[missing_family] <- hits$structure_family[missing_family]
    missing_role <- pass & (base::is.na(hits$enzyme_role) | !base::nzchar(hits$enzyme_role)) &
      !base::is.na(hits$structure_role) & base::nzchar(hits$structure_role)
    hits$enzyme_role[missing_role] <- hits$structure_role[missing_role]
    support_add <- base::vapply(base::seq_len(base::nrow(hits)), function(i) {
      if (!isTRUE(pass[i])) return(NA_character_)
      .dnmb_rebasefinder_structure_support(hits[i, , drop = FALSE])
    }, character(1))
    add <- !base::is.na(support_add) & base::nzchar(support_add)
    hits$support[add] <- base::ifelse(
      base::is.na(hits$support[add]) | !base::nzchar(hits$support[add]),
      support_add[add],
      base::paste(hits$support[add], support_add[add], sep = "; ")
    )
  }

  new_struct <- structure_tbl[structure_tbl$structure_pass & !structure_tbl$query %in% hits$query, , drop = FALSE]
  if (base::nrow(new_struct)) {
    gene_idx <- base::match(new_struct$query, genes$locus_tag)
    new_struct <- new_struct[!base::is.na(gene_idx), , drop = FALSE]
    gene_idx <- gene_idx[!base::is.na(gene_idx)]
  }
  if (base::nrow(new_struct)) {
    add <- data.frame(
      query = new_struct$query,
      source = "rebasefinder",
      family_system = "REBASEfinder",
      family_id = base::ifelse(!base::is.na(new_struct$structure_family) & base::nzchar(new_struct$structure_family),
                               new_struct$structure_family, "Structure-supported R-M"),
      hit_label = new_struct$structure_hit,
      enzyme_role = new_struct$structure_role,
      evidence_mode = "structure_only",
      substrate_label = NA_character_,
      support = base::vapply(base::seq_len(base::nrow(new_struct)), function(i) {
        .dnmb_rebasefinder_structure_support(new_struct[i, , drop = FALSE])
      }, character(1)),
      typing_eligible = TRUE,
      stringsAsFactors = FALSE
    )
    for (col in extra_cols) add[[col]] <- new_struct[[col]]
    add_cons <- .dnmb_rebasefinder_structure_candidate_consistency(
      add$family_id,
      add$enzyme_role,
      if ("structure_family" %in% base::names(add)) add$structure_family else NA_character_,
      if ("structure_role" %in% base::names(add)) add$structure_role else NA_character_,
      if ("structure_chain_role" %in% base::names(add)) add$structure_chain_role else NA_character_
    )
    add$structure_family_consistent <- add_cons$structure_family_consistent
    add$structure_role_consistent <- add_cons$structure_role_consistent
    add$structure_candidate_consistent <- add_cons$structure_candidate_consistent
    for (col in base::setdiff(base::names(hits), base::names(add))) {
      add[[col]] <- .dnmb_na_vector_like(hits[[col]], base::nrow(add))
    }
    hits <- base::rbind(hits, add[, base::names(hits), drop = FALSE])
  }
  hits
}

.dnmb_rebasefinder_add_typeiii_context <- function(hits,
                                                   genes,
                                                   max_operon_gap = 5000,
                                                   max_intervening = 2L,
                                                   max_neighbors = 3L) {
  if (base::is.null(hits) || !base::is.data.frame(hits)) return(hits)
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  genes_ord <- .dnmb_rebasefinder_gene_order(genes)
  if (!base::nrow(genes_ord)) return(hits)

  anno_text <- .dnmb_rebasefinder_gene_annotation_text(genes_ord)
  hit_idx <- base::match(genes_ord$locus_tag, hits$query)
  genes_ord$.rm_type <- NA_character_
  genes_ord$.rm_role <- NA_character_
  genes_ord$.hit_label <- NA_character_
  have_hit <- !base::is.na(hit_idx)
  genes_ord$.rm_type[have_hit] <- hits$family_id[hit_idx[have_hit]]
  genes_ord$.rm_role[have_hit] <- hits$enzyme_role[hit_idx[have_hit]]
  genes_ord$.hit_label[have_hit] <- hits$hit_label[hit_idx[have_hit]]
  context_role <- base::mapply(
    .dnmb_rebasefinder_typeiii_context_role,
    text = anno_text,
    known_rm_type = genes_ord$.rm_type,
    known_role = genes_ord$.rm_role,
    USE.NAMES = FALSE
  )
  translations <- if ("translation" %in% base::names(genes_ord)) genes_ord$translation else base::rep(NA_character_, base::nrow(genes_ord))
  seq_signals <- base::lapply(translations, .dnmb_rebasefinder_typeiii_sequence_signals)
  seq_role <- base::vapply(seq_signals, `[[`, character(1), "role")
  seq_signal_label <- base::vapply(seq_signals, function(x) {
    sig <- x$signals
    if (!base::length(sig)) NA_character_ else base::paste(sig, collapse = "+")
  }, character(1))
  genes_ord$.typeiii_role <- context_role
  use_seq_role <- base::is.na(genes_ord$.typeiii_role) & !base::is.na(seq_role) & seq_role %in% c("M", "R")
  genes_ord$.typeiii_role[use_seq_role] <- seq_role[use_seq_role]
  evidence_mode <- if ("evidence_mode" %in% base::names(hits)) base::as.character(hits$evidence_mode) else base::rep(NA_character_, base::nrow(hits))
  genes_ord$.evidence_mode <- NA_character_
  genes_ord$.evidence_mode[have_hit] <- evidence_mode[hit_idx[have_hit]]
  operon_context_hit <- !base::is.na(genes_ord$.evidence_mode) & genes_ord$.evidence_mode == "operon_context"
  known_typeiii_role <- !operon_context_hit &
    !base::is.na(genes_ord$.rm_type) & genes_ord$.rm_type == "Type III" &
    !base::is.na(genes_ord$.rm_role) & genes_ord$.rm_role %in% c("M", "R", "S")
  operon_source <- base::ifelse(
    operon_context_hit & base::grepl(":sequence_motif", genes_ord$.hit_label, fixed = TRUE),
    "sequence_motif",
    base::ifelse(
      operon_context_hit & base::grepl(":annotation", genes_ord$.hit_label, fixed = TRUE),
      "annotation",
      base::ifelse(operon_context_hit, "operon_context", NA_character_)
    )
  )
  genes_ord$.typeiii_role_source <- base::ifelse(
    known_typeiii_role,
    "rebase_hit",
    base::ifelse(!base::is.na(operon_source), operon_source,
                 base::ifelse(!base::is.na(context_role), "annotation",
                 base::ifelse(use_seq_role, "sequence_motif", NA_character_))
    )
  )
  genes_ord$.typeiii_sequence_signals <- seq_signal_label

  is_typeiii_seed <- !base::is.na(genes_ord$.rm_type) & genes_ord$.rm_type == "Type III"
  if (!base::any(is_typeiii_seed)) {
    for (col in c("typeiii_context_status", "typeiii_context_roles", "typeiii_context_partners",
                  "typeiii_context_sources", "typeiii_operon_supported", "typeiii_context_window_bp")) {
      if (!col %in% base::names(hits)) hits[[col]] <- NA_character_
    }
    return(hits)
  }

  context_rows <- list()
  for (seed_i in base::which(is_typeiii_seed)) {
    seed <- genes_ord[seed_i, , drop = FALSE]
    same_contig <- genes_ord$contig == seed$contig
    idx <- base::which(same_contig & base::abs(genes_ord$.dnmb_order - seed$.dnmb_order) <= max_neighbors)
    idx <- idx[idx != seed_i]
    if (base::length(idx)) {
      distances <- base::vapply(idx, function(j) {
        .dnmb_rebasefinder_intergenic_gap(seed, genes_ord[j, , drop = FALSE])
      }, numeric(1))
      idx <- idx[base::is.na(distances) | distances <= max_operon_gap]
    }
    ctx_idx <- c(seed_i, idx)
    ctx <- genes_ord[ctx_idx, , drop = FALSE]
    ctx$distance_bp <- base::vapply(ctx_idx, function(j) {
      if (j == seed_i) return(0)
      .dnmb_rebasefinder_intergenic_gap(seed, genes_ord[j, , drop = FALSE])
    }, numeric(1))
    ctx$relative_gene <- ctx$.dnmb_order - seed$.dnmb_order
    ctx$same_direction <- base::is.na(seed$direction) | base::is.na(ctx$direction) | seed$direction == ctx$direction
    ctx$operon_link <- ctx$same_direction & (base::is.na(ctx$distance_bp) | ctx$distance_bp <= max_operon_gap) &
      base::abs(ctx$relative_gene) <= max_intervening + 1L

    roles_present <- base::sort(base::unique(stats::na.omit(ctx$.typeiii_role[ctx$operon_link | ctx$relative_gene == 0])))
    source_rows <- ctx[!base::is.na(ctx$.typeiii_role) & (ctx$operon_link | ctx$relative_gene == 0), , drop = FALSE]
    source_label <- if (base::nrow(source_rows)) {
      base::paste(
        base::sort(base::unique(base::paste0(source_rows$.typeiii_role, "=", source_rows$.typeiii_role_source))),
        collapse = ";"
      )
    } else {
      NA_character_
    }
    has_m <- "M" %in% roles_present
    has_r <- "R" %in% roles_present
    status <- if (has_m && has_r) "complete_mod_res"
      else if (identical(seed$.typeiii_role, "M")) "missing_res_check_neighbors"
      else if (identical(seed$.typeiii_role, "R")) "missing_mod_check_neighbors"
      else "typeiii_context_weak"
    partners <- ctx[ctx$relative_gene != 0 & !base::is.na(ctx$.typeiii_role) & ctx$operon_link, , drop = FALSE]
    partner_label <- if (base::nrow(partners)) {
      base::paste(
        base::sprintf("%s:%s:%s:%+dg:%sbp", partners$locus_tag, partners$.typeiii_role,
                      partners$.typeiii_role_source,
                      partners$relative_gene, base::ifelse(base::is.na(partners$distance_bp), "NA", partners$distance_bp)),
        collapse = " | "
      )
    } else {
      NA_character_
    }
    context_rows[[base::length(context_rows) + 1L]] <- data.frame(
      query = seed$locus_tag,
      typeiii_context_status = status,
      typeiii_context_roles = base::paste(roles_present, collapse = "/"),
      typeiii_context_partners = partner_label,
      typeiii_context_sources = source_label,
      typeiii_operon_supported = has_m && has_r,
      typeiii_context_window_bp = base::max(ctx$end, na.rm = TRUE) - base::min(ctx$start, na.rm = TRUE) + 1,
      stringsAsFactors = FALSE
    )
  }
  context_tbl <- do.call(rbind, context_rows)

  for (col in base::setdiff(base::names(context_tbl), "query")) {
    if (!col %in% base::names(hits)) hits[[col]] <- .dnmb_na_vector_like(context_tbl[[col]], base::nrow(hits))
  }
  hidx <- base::match(hits$query, context_tbl$query)
  include <- !base::is.na(hidx)
  for (col in base::setdiff(base::names(context_tbl), "query")) hits[[col]][include] <- context_tbl[[col]][hidx[include]]

  support_add <- base::ifelse(
    include,
    base::paste0("typeIII_context=", hits$typeiii_context_status, "; partners=", hits$typeiii_context_partners),
    NA_character_
  )
  add <- !base::is.na(support_add) & base::nzchar(support_add)
  if (base::any(add)) {
    hits$support[add] <- .dnmb_rebasefinder_strip_typeiii_support(hits$support[add])
  }
  hits$support[add] <- base::ifelse(
    base::is.na(hits$support[add]) | !base::nzchar(hits$support[add]),
    support_add[add],
    base::paste(hits$support[add], support_add[add], sep = "; ")
  )
  hits$typing_eligible[add & hits$typeiii_operon_supported %in% TRUE] <- TRUE

  # Add likely Type III Res/Mod partners that were absent from REBASE BLAST hits.
  seed_idx <- base::which(is_typeiii_seed)
  candidate_idx <- base::which(
    !genes_ord$locus_tag %in% hits$query &
      !base::is.na(genes_ord$.typeiii_role) &
      genes_ord$.typeiii_role %in% c("M", "R") &
      base::vapply(base::seq_len(base::nrow(genes_ord)), function(i) {
        any(base::vapply(seed_idx, function(si) {
          genes_ord$contig[i] == genes_ord$contig[si] &&
            base::abs(genes_ord$.dnmb_order[i] - genes_ord$.dnmb_order[si]) <= max_intervening + 2L &&
            {
              gap <- .dnmb_rebasefinder_intergenic_gap(genes_ord[i, , drop = FALSE], genes_ord[si, , drop = FALSE])
              base::is.na(gap) || gap <= max_operon_gap
            }
        }, logical(1)))
      }, logical(1))
  )
  if (base::length(candidate_idx)) {
    add_genes <- genes_ord[candidate_idx, , drop = FALSE]
    add_support <- base::ifelse(
      add_genes$.typeiii_role_source == "sequence_motif",
      base::paste0(
        "Type III operon-context candidate absent from primary REBASE BLAST hit table; role inferred from protein motifs=",
        base::ifelse(base::is.na(add_genes$.typeiii_sequence_signals), "NA", add_genes$.typeiii_sequence_signals)
      ),
      "Type III operon-context candidate absent from primary REBASE BLAST hit table"
    )
    add <- data.frame(
      query = add_genes$locus_tag,
      source = "rebasefinder",
      family_system = "REBASEfinder",
      family_id = "Type III",
      hit_label = base::paste0("typeIII_context:", add_genes$.typeiii_role, ":", add_genes$.typeiii_role_source),
      enzyme_role = add_genes$.typeiii_role,
      evidence_mode = "operon_context",
      substrate_label = NA_character_,
      support = add_support,
      typing_eligible = FALSE,
      stringsAsFactors = FALSE
    )
    for (col in base::setdiff(base::names(hits), base::names(add))) {
      add[[col]] <- .dnmb_na_vector_like(hits[[col]], base::nrow(add))
    }
    hits <- base::rbind(hits, add[, base::names(hits), drop = FALSE])
  }

  hits
}

.dnmb_rebasefinder_add_typei_context <- function(hits,
                                                 genes,
                                                 max_operon_gap = 5000,
                                                 max_intervening = 1L,
                                                 max_neighbors = 4L) {
  if (base::is.null(hits) || !base::is.data.frame(hits)) return(hits)
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  genes_ord <- .dnmb_rebasefinder_gene_order(genes)
  if (!base::nrow(genes_ord)) return(hits)

  anno_text <- .dnmb_rebasefinder_gene_annotation_text(genes_ord)
  hit_idx <- base::match(genes_ord$locus_tag, hits$query)
  genes_ord$.rm_type <- NA_character_
  genes_ord$.rm_role <- NA_character_
  genes_ord$.hit_label <- NA_character_
  have_hit <- !base::is.na(hit_idx)
  genes_ord$.rm_type[have_hit] <- hits$family_id[hit_idx[have_hit]]
  genes_ord$.rm_role[have_hit] <- hits$enzyme_role[hit_idx[have_hit]]
  genes_ord$.hit_label[have_hit] <- hits$hit_label[hit_idx[have_hit]]
  context_role <- base::mapply(
    .dnmb_rebasefinder_typei_context_role,
    text = anno_text,
    known_rm_type = genes_ord$.rm_type,
    known_role = genes_ord$.rm_role,
    USE.NAMES = FALSE
  )
  translations <- if ("translation" %in% base::names(genes_ord)) genes_ord$translation else base::rep(NA_character_, base::nrow(genes_ord))
  seq_signals <- base::lapply(translations, .dnmb_rebasefinder_typei_sequence_signals)
  seq_role <- base::vapply(seq_signals, `[[`, character(1), "role")
  seq_signal_label <- base::vapply(seq_signals, function(x) {
    sig <- x$signals
    if (!base::length(sig)) NA_character_ else base::paste(sig, collapse = "+")
  }, character(1))
  genes_ord$.typei_role <- context_role
  use_seq_role <- base::is.na(genes_ord$.typei_role) & !base::is.na(seq_role) & seq_role %in% c("M", "R")
  genes_ord$.typei_role[use_seq_role] <- seq_role[use_seq_role]
  evidence_mode <- if ("evidence_mode" %in% base::names(hits)) base::as.character(hits$evidence_mode) else base::rep(NA_character_, base::nrow(hits))
  genes_ord$.evidence_mode <- NA_character_
  genes_ord$.evidence_mode[have_hit] <- evidence_mode[hit_idx[have_hit]]
  operon_context_hit <- !base::is.na(genes_ord$.evidence_mode) & genes_ord$.evidence_mode == "operon_context"
  known_typei_role <- !operon_context_hit &
    !base::is.na(genes_ord$.rm_type) & genes_ord$.rm_type == "Type I" &
    !base::is.na(genes_ord$.rm_role) & genes_ord$.rm_role %in% c("M", "R", "S")
  operon_source <- base::ifelse(
    operon_context_hit & base::grepl(":sequence_motif", genes_ord$.hit_label, fixed = TRUE),
    "sequence_motif",
    base::ifelse(
      operon_context_hit & base::grepl(":annotation", genes_ord$.hit_label, fixed = TRUE),
      "annotation",
      base::ifelse(operon_context_hit, "operon_context", NA_character_)
    )
  )
  genes_ord$.typei_role_source <- base::ifelse(
    known_typei_role,
    "rebase_hit",
    base::ifelse(!base::is.na(operon_source), operon_source,
                 base::ifelse(!base::is.na(context_role), "annotation",
                              base::ifelse(use_seq_role, "sequence_motif", NA_character_))
    )
  )
  genes_ord$.typei_sequence_signals <- seq_signal_label

  is_typei_seed <- !base::is.na(genes_ord$.rm_type) & genes_ord$.rm_type == "Type I"
  if (!base::any(is_typei_seed)) {
    for (col in c("typei_context_status", "typei_context_roles", "typei_context_partners",
                  "typei_context_sources", "typei_operon_supported", "typei_context_window_bp")) {
      if (!col %in% base::names(hits)) hits[[col]] <- NA_character_
    }
    return(hits)
  }

  context_rows <- list()
  for (seed_i in base::which(is_typei_seed)) {
    seed <- genes_ord[seed_i, , drop = FALSE]
    same_contig <- genes_ord$contig == seed$contig
    idx <- base::which(same_contig & base::abs(genes_ord$.dnmb_order - seed$.dnmb_order) <= max_neighbors)
    idx <- idx[idx != seed_i]
    if (base::length(idx)) {
      distances <- base::vapply(idx, function(j) {
        .dnmb_rebasefinder_intergenic_gap(seed, genes_ord[j, , drop = FALSE])
      }, numeric(1))
      idx <- idx[base::is.na(distances) | distances <= max_operon_gap]
    }
    ctx_idx <- c(seed_i, idx)
    ctx <- genes_ord[ctx_idx, , drop = FALSE]
    ctx$distance_bp <- base::vapply(ctx_idx, function(j) {
      if (j == seed_i) return(0)
      .dnmb_rebasefinder_intergenic_gap(seed, genes_ord[j, , drop = FALSE])
    }, numeric(1))
    ctx$relative_gene <- ctx$.dnmb_order - seed$.dnmb_order
    ctx$same_direction <- base::is.na(seed$direction) | base::is.na(ctx$direction) | seed$direction == ctx$direction
    ctx$operon_link <- ctx$same_direction &
      (base::is.na(ctx$distance_bp) | ctx$distance_bp <= max_operon_gap) &
      base::abs(ctx$relative_gene) <= max_intervening + 2L

    roles_present <- base::sort(base::unique(stats::na.omit(ctx$.typei_role[ctx$operon_link | ctx$relative_gene == 0])))
    source_rows <- ctx[!base::is.na(ctx$.typei_role) & (ctx$operon_link | ctx$relative_gene == 0), , drop = FALSE]
    source_label <- if (base::nrow(source_rows)) {
      base::paste(
        base::sort(base::unique(base::paste0(source_rows$.typei_role, "=", source_rows$.typei_role_source))),
        collapse = ";"
      )
    } else {
      NA_character_
    }
    has_m <- "M" %in% roles_present
    has_r <- "R" %in% roles_present
    has_s <- "S" %in% roles_present
    status <- if (has_m && has_r && has_s) "complete_mrs"
      else if (has_m && has_r) "missing_specificity_check_neighbors"
      else if (has_m && has_s) "missing_restriction_check_neighbors"
      else if (has_r && has_s) "missing_methylase_check_neighbors"
      else "typei_context_weak"
    partners <- ctx[ctx$relative_gene != 0 & !base::is.na(ctx$.typei_role) & ctx$operon_link, , drop = FALSE]
    partner_label <- if (base::nrow(partners)) {
      base::paste(
        base::sprintf("%s:%s:%s:%+dg:%sbp", partners$locus_tag, partners$.typei_role,
                      partners$.typei_role_source,
                      partners$relative_gene, base::ifelse(base::is.na(partners$distance_bp), "NA", partners$distance_bp)),
        collapse = " | "
      )
    } else {
      NA_character_
    }
    context_rows[[base::length(context_rows) + 1L]] <- data.frame(
      query = seed$locus_tag,
      typei_context_status = status,
      typei_context_roles = base::paste(roles_present, collapse = "/"),
      typei_context_partners = partner_label,
      typei_context_sources = source_label,
      typei_operon_supported = has_m && has_r && has_s,
      typei_context_window_bp = base::max(ctx$end, na.rm = TRUE) - base::min(ctx$start, na.rm = TRUE) + 1,
      stringsAsFactors = FALSE
    )
  }
  context_tbl <- do.call(rbind, context_rows)

  for (col in base::setdiff(base::names(context_tbl), "query")) {
    if (!col %in% base::names(hits)) hits[[col]] <- .dnmb_na_vector_like(context_tbl[[col]], base::nrow(hits))
  }
  hidx <- base::match(hits$query, context_tbl$query)
  include <- !base::is.na(hidx)
  for (col in base::setdiff(base::names(context_tbl), "query")) hits[[col]][include] <- context_tbl[[col]][hidx[include]]

  support_add <- base::ifelse(
    include,
    base::paste0("typeI_context=", hits$typei_context_status, "; partners=", hits$typei_context_partners),
    NA_character_
  )
  add <- !base::is.na(support_add) & base::nzchar(support_add)
  if (base::any(add)) {
    hits$support[add] <- .dnmb_rebasefinder_strip_typei_support(hits$support[add])
  }
  hits$support[add] <- base::ifelse(
    base::is.na(hits$support[add]) | !base::nzchar(hits$support[add]),
    support_add[add],
    base::paste(hits$support[add], support_add[add], sep = "; ")
  )
  hits$typing_eligible[add & hits$typei_operon_supported %in% TRUE] <- TRUE

  seed_idx <- base::which(is_typei_seed)
  candidate_idx <- base::which(
    !genes_ord$locus_tag %in% hits$query &
      !base::is.na(genes_ord$.typei_role) &
      genes_ord$.typei_role %in% c("M", "R", "S") &
      base::vapply(base::seq_len(base::nrow(genes_ord)), function(i) {
        any(base::vapply(seed_idx, function(si) {
          same_direction <- base::is.na(genes_ord$direction[i]) | base::is.na(genes_ord$direction[si]) |
            genes_ord$direction[i] == genes_ord$direction[si]
          genes_ord$contig[i] == genes_ord$contig[si] &&
            same_direction &&
            base::abs(genes_ord$.dnmb_order[i] - genes_ord$.dnmb_order[si]) <= max_neighbors &&
            {
              gap <- .dnmb_rebasefinder_intergenic_gap(genes_ord[i, , drop = FALSE], genes_ord[si, , drop = FALSE])
              base::is.na(gap) || gap <= max_operon_gap
            }
        }, logical(1)))
      }, logical(1))
  )
  if (base::length(candidate_idx)) {
    add_genes <- genes_ord[candidate_idx, , drop = FALSE]
    add_support <- base::ifelse(
      add_genes$.typei_role_source == "sequence_motif",
      base::paste0(
        "Type I operon-context candidate absent from primary REBASE BLAST hit table; role inferred from protein motifs=",
        base::ifelse(base::is.na(add_genes$.typei_sequence_signals), "NA", add_genes$.typei_sequence_signals)
      ),
      "Type I operon-context candidate absent from primary REBASE BLAST hit table"
    )
    add <- data.frame(
      query = add_genes$locus_tag,
      source = "rebasefinder",
      family_system = "REBASEfinder",
      family_id = "Type I",
      hit_label = base::paste0("typeI_context:", add_genes$.typei_role, ":", add_genes$.typei_role_source),
      enzyme_role = add_genes$.typei_role,
      evidence_mode = "operon_context",
      substrate_label = NA_character_,
      support = add_support,
      typing_eligible = FALSE,
      stringsAsFactors = FALSE
    )
    for (col in base::setdiff(base::names(hits), base::names(add))) {
      add[[col]] <- .dnmb_na_vector_like(hits[[col]], base::nrow(add))
    }
    hits <- base::rbind(hits, add[, base::names(hits), drop = FALSE])
  }

  hits
}

.dnmb_rebasefinder_typeiv_candidate_mask <- function(genes) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(genes)) return(logical())
  text <- .dnmb_rebasefinder_gene_annotation_text(genes)
  seqs <- if ("translation" %in% base::names(genes)) {
    base::toupper(base::gsub("[^A-Za-z]", "", base::as.character(genes$translation)))
  } else {
    base::rep(NA_character_, base::nrow(genes))
  }
  mrr_like <- !base::is.na(seqs) & base::nzchar(seqs) &
    base::grepl("D.{8,15}[EQ].[KR].{20,60}[DE].{0,5}[KR]", seqs, perl = TRUE)
  explicit_typeiv <- base::grepl(
    "type[ _-]?iv restriction|\\bmrr\\b|\\bmcr[abc]?\\b|modified[ -]?dna.*restriction|methylated[ -]?dna.*restriction|restriction.*modified[ -]?dna|modification-dependent restriction",
    text,
    perl = TRUE
  )
  generic_rease <- base::grepl("\\brestriction endonuclease(-like)?( protein)?\\b", text, perl = TRUE)
  not_rm_subunit <- !base::grepl("subunit s|type[ _-]?i restriction|type[ _-]?ii|type[ _-]?iii|crispr|cas[0-9]?|tnpb|uma2|mutl|muts|iscb|endonuclease iii|endonuclease q", text, perl = TRUE)
  (explicit_typeiv | (generic_rease & not_rm_subunit & mrr_like))
}

.dnmb_rebasefinder_add_typeiv_candidates <- function(hits, genes) {
  if (base::is.null(hits) || !base::is.data.frame(hits)) return(hits)
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes_ord <- .dnmb_rebasefinder_gene_order(genes)
  if (!base::nrow(genes_ord)) return(hits)
  mask <- .dnmb_rebasefinder_typeiv_candidate_mask(genes_ord)
  if (!base::any(mask, na.rm = TRUE)) return(hits)
  add_genes <- genes_ord[mask & !genes_ord$locus_tag %in% hits$query, , drop = FALSE]
  if (!base::nrow(add_genes)) return(hits)
  add <- data.frame(
    query = add_genes$locus_tag,
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = "Type IV",
    hit_label = "typeIV_context:Mrr_or_modified_DNA_REase",
    enzyme_role = "R",
    evidence_mode = "annotation_motif_candidate",
    substrate_label = NA_character_,
    support = "Type IV restriction candidate from modified-DNA/Mrr-like annotation or Mrr-like motif plus restriction-endonuclease annotation",
    typing_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  for (col in base::setdiff(base::names(hits), base::names(add))) {
    add[[col]] <- .dnmb_na_vector_like(hits[[col]], base::nrow(add))
  }
  hits <- base::rbind(hits, add[, base::names(hits), drop = FALSE])
  hits
}

.dnmb_rebasefinder_add_sequence_partial_status <- function(hits, genes) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  if (!"locus_tag" %in% base::names(genes)) {
    hits$aa_len <- NA_integer_
    hits$expected_min_aa <- NA_integer_
    hits$partial_status <- NA_character_
    hits$partial_reason <- NA_character_
    return(hits)
  }
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  idx <- base::match(.dnmb_module_clean_annotation_key(hits$query), genes$locus_tag)
  gene_tbl <- genes[idx, , drop = FALSE]
  gene_tbl$family_id <- hits$family_id
  gene_tbl$enzyme_role <- hits$enzyme_role
  gene_tbl$support <- hits$support
  partial <- .dnmb_rebasefinder_sequence_partial_table(
    gene_tbl,
    family_col = "family_id",
    role_col = "enzyme_role"
  )
  hits$aa_len <- partial$aa_len
  hits$expected_min_aa <- partial$expected_min_aa
  hits$partial_status <- partial$partial_status
  hits$partial_reason <- partial$partial_reason
  hits
}

.dnmb_rebasefinder_supplemental_query_mask <- function(hits) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(logical())
  family <- if ("family_id" %in% base::names(hits)) base::as.character(hits$family_id) else base::rep(NA_character_, base::nrow(hits))
  role <- if ("enzyme_role" %in% base::names(hits)) base::as.character(hits$enzyme_role) else base::rep(NA_character_, base::nrow(hits))
  evidence <- if ("evidence_mode" %in% base::names(hits)) base::as.character(hits$evidence_mode) else base::rep(NA_character_, base::nrow(hits))
  blast_missing <- if ("blast_identity" %in% base::names(hits)) base::is.na(hits$blast_identity) else base::rep(TRUE, base::nrow(hits))
  family %in% c("Type I", "Type III", "Type IV") &
    role %in% c("R", "S") &
    blast_missing &
    evidence %in% c("operon_context", "annotation_motif_candidate", "structure_supported", "structure_only")
}

.dnmb_rebasefinder_cached_rebase_data <- function(cache_root = NULL) {
  global_db <- base::get0("rmscan_rebase_db", envir = .GlobalEnv, inherits = FALSE)
  if (base::is.data.frame(global_db) && base::nrow(global_db)) return(global_db)

  cache_dir <- .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = FALSE)
  rds <- base::file.path(cache_dir, "rebase_data.rds")
  if (base::file.exists(rds)) {
    out <- tryCatch(base::readRDS(rds), error = function(e) NULL)
    if (base::is.data.frame(out) && base::nrow(out)) return(out)
  }
  NULL
}

.dnmb_rebasefinder_set_defenseviz_cache <- function(cache_root = NULL, verbose = FALSE) {
  cache_dir <- .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE)
  if (!requireNamespace("DefenseViz", quietly = TRUE)) return(cache_dir)
  tryCatch({
    cache_path <- cache_dir
    override_fn <- local({
      path <- cache_path
      function() {
        if (!base::dir.exists(path)) base::dir.create(path, recursive = TRUE, showWarnings = FALSE)
        path
      }
    })
    utils::assignInNamespace("get_rebase_cache_dir", override_fn, ns = "DefenseViz")
    if (isTRUE(verbose)) message("[REBASEfinder] REBASE cache: ", cache_dir)
  }, error = function(e) {
    if (isTRUE(verbose)) message("[REBASEfinder] Could not override cache dir: ", conditionMessage(e))
  })
  cache_dir
}

.dnmb_rebasefinder_rebase_fasta <- function(cache_root = NULL) {
  cache_dir <- .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE)
  candidates <- c(base::file.path(cache_dir, "rebase_db.fasta"))
  for (fasta in candidates) {
    if (base::file.exists(fasta) && base::file.info(fasta)$size > 1e6) return(fasta)
  }
  NA_character_
}

.dnmb_rebasefinder_run_blastp <- function(query_fasta,
                                          rebase_fasta,
                                          output_file,
                                          evalue = 1e-5,
                                          num_threads = 1L,
                                          max_target_seqs = 10L,
                                          verbose = TRUE) {
  blastp <- Sys.which("blastp")
  makeblastdb <- Sys.which("makeblastdb")
  if (!nzchar(blastp)) stop("blastp not found", call. = FALSE)
  if (!nzchar(makeblastdb)) stop("makeblastdb not found", call. = FALSE)
  blast_work <- tempfile("dnmb_rebase_blast_")
  base::dir.create(blast_work, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(blast_work, recursive = TRUE, force = TRUE), add = TRUE)
  safe_query <- base::file.path(blast_work, "query.faa")
  safe_rebase <- base::file.path(blast_work, "rebase_db.fasta")
  safe_db <- base::file.path(blast_work, "rebase_db")
  safe_out <- base::file.path(blast_work, "blast.tsv")
  if (!base::file.copy(query_fasta, safe_query, overwrite = TRUE)) {
    stop("Could not stage supplemental query FASTA", call. = FALSE)
  }
  linked <- suppressWarnings(base::file.symlink(rebase_fasta, safe_rebase))
  if (!isTRUE(linked) && !base::file.copy(rebase_fasta, safe_rebase, overwrite = TRUE)) {
    stop("Could not stage REBASE FASTA for supplemental BLAST", call. = FALSE)
  }

  if (isTRUE(verbose)) message("  [1/2] Creating BLAST database from REBASE FASTA...")
  mk <- system2(
    makeblastdb,
    args = c("-in", safe_rebase, "-dbtype", "prot", "-out", safe_db, "-title", "REBASE"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(mk, "status")
  if (is.null(status)) status <- 0L
  db_files <- paste0(safe_db, c(".phr", ".pin", ".psq"))
  if (!identical(as.integer(status), 0L) || !all(base::file.exists(db_files))) {
    stop("makeblastdb failed for supplemental REBASE BLAST: ", paste(mk, collapse = " "), call. = FALSE)
  }

  if (isTRUE(verbose)) {
    query_count <- length(grep("^>", readLines(safe_query, warn = FALSE)))
    message("  [2/2] Running supplemental BLASTP: ", query_count, " queries vs REBASE database...")
  }
  args <- c(
    "-query", safe_query,
    "-db", safe_db,
    "-out", safe_out,
    "-evalue", as.character(evalue),
    "-num_threads", as.character(max(1L, as.integer(num_threads))),
    "-max_target_seqs", as.character(max(1L, as.integer(max_target_seqs))),
    "-outfmt", shQuote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
  )
  run <- system2(blastp, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(run, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    stop("blastp failed for supplemental REBASE BLAST: ", paste(run, collapse = " "), call. = FALSE)
  }
  if (base::file.exists(safe_out) && base::file.info(safe_out)$size > 0) {
    base::dir.create(base::dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (!base::file.copy(safe_out, output_file, overwrite = TRUE)) {
      stop("Could not copy supplemental BLAST output to final path", call. = FALSE)
    }
    output_file
  } else {
    NULL
  }
}

.dnmb_rebasefinder_best_from_blast <- function(blast_tbl, rebase_data = NULL, min_identity = 0.10) {
  if (base::is.null(blast_tbl) || !base::is.data.frame(blast_tbl) || !base::nrow(blast_tbl)) {
    return(base::data.frame())
  }
  blast_tbl <- base::as.data.frame(blast_tbl, stringsAsFactors = FALSE)
  if (!"pct_identity" %in% base::names(blast_tbl) && "pident" %in% base::names(blast_tbl)) {
    blast_tbl$pct_identity <- blast_tbl$pident / 100
  }
  blast_tbl <- blast_tbl[!base::is.na(blast_tbl$pct_identity) & blast_tbl$pct_identity >= min_identity, , drop = FALSE]
  if (!base::nrow(blast_tbl)) return(base::data.frame())

  if (base::is.data.frame(rebase_data) && base::nrow(rebase_data)) {
    best <- tryCatch(
      DefenseViz::get_best_rebase_match(blast_tbl, rebase_data, min_identity = min_identity),
      error = function(e) NULL
    )
    if (base::is.data.frame(best) && base::nrow(best)) return(base::as.data.frame(best, stringsAsFactors = FALSE))
  }

  blast_tbl <- blast_tbl[base::order(blast_tbl$query_id, -blast_tbl$pct_identity, -blast_tbl$bitscore), , drop = FALSE]
  blast_tbl <- blast_tbl[!base::duplicated(blast_tbl$query_id), , drop = FALSE]
  subject <- if ("rebase_enzyme" %in% base::names(blast_tbl)) {
    base::as.character(blast_tbl$rebase_enzyme)
  } else {
    base::as.character(blast_tbl$subject_id)
  }
  base::data.frame(
    query_id = blast_tbl$query_id,
    best_rebase_match = subject,
    best_match_identity = blast_tbl$pct_identity,
    rm_type = NA_character_,
    subunit = .dnmb_rebasefinder_role_from_hit(subject),
    rec_seq = NA_character_,
    evalue = blast_tbl$evalue,
    bitscore = blast_tbl$bitscore,
    length = blast_tbl$length,
    stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_supplemental_rebase_blast <- function(hits,
                                                          genes,
                                                          output_dir,
                                                          blast_min_identity = 0.10,
                                                          blast_min_length = 50,
                                                          cache_root = NULL,
                                                          cpu = 1L,
                                                          download_rebase_if_missing = TRUE,
                                                          verbose = TRUE) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  empty_status <- .dnmb_rebasefinder_status_row(character(), character(), character())
  if (!base::nrow(hits) || !"translation" %in% base::names(genes)) {
    return(list(hits = hits, status = empty_status))
  }

  mask <- .dnmb_rebasefinder_supplemental_query_mask(hits)
  if (!base::any(mask)) {
    return(list(hits = hits, status = .dnmb_rebasefinder_status_row(
      "supplemental_rebase_blast", "skipped", "No R/S/Type IV context candidates missing REBASE BLAST"
    )))
  }

  .dnmb_rebasefinder_set_defenseviz_cache(cache_root = cache_root, verbose = FALSE)
  rebase_fasta <- .dnmb_rebasefinder_rebase_fasta(cache_root = cache_root)
  rebase_data <- .dnmb_rebasefinder_cached_rebase_data(cache_root = cache_root)
  if ((base::is.na(rebase_fasta) || !base::file.exists(rebase_fasta)) &&
      !base::is.data.frame(rebase_data) &&
      isTRUE(download_rebase_if_missing) &&
      requireNamespace("DefenseViz", quietly = TRUE)) {
    rebase_data <- tryCatch(
      DefenseViz::get_rebase_data(force_download = FALSE, verbose = verbose),
      error = function(e) {
        if (isTRUE(verbose)) message("[REBASEfinder] Could not load/download REBASE data for supplemental BLAST: ", conditionMessage(e))
        NULL
      }
    )
    rebase_fasta <- .dnmb_rebasefinder_rebase_fasta(cache_root = cache_root)
  }
  if (base::is.na(rebase_fasta) && base::is.data.frame(rebase_data) && base::nrow(rebase_data)) {
    rebase_fasta <- base::file.path(.dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE), "rebase_db.fasta")
    rebase_for_fasta <- rebase_data
    if (!"locus_tag" %in% base::names(rebase_for_fasta) && "enzyme_name" %in% base::names(rebase_for_fasta)) {
      rebase_for_fasta$locus_tag <- rebase_for_fasta$enzyme_name
    }
    if (!"translation" %in% base::names(rebase_for_fasta) && "sequence" %in% base::names(rebase_for_fasta)) {
      rebase_for_fasta$translation <- rebase_for_fasta$sequence
    }
    write_fasta_for_blast <- base::get("write_fasta_for_blast", envir = base::asNamespace("DefenseViz"))
    tryCatch(
      write_fasta_for_blast(rebase_for_fasta, rebase_fasta, id_col = "locus_tag", seq_col = "translation"),
      error = function(e) {
        if (isTRUE(verbose)) message("[REBASEfinder] Could not write supplemental REBASE FASTA: ", conditionMessage(e))
      }
    )
  }
  if (base::is.na(rebase_fasta) || !base::file.exists(rebase_fasta)) {
    return(list(hits = hits, status = .dnmb_rebasefinder_status_row(
      "supplemental_rebase_blast", "skipped",
      "REBASE FASTA/database not available and could not be created; supplemental R/S BLAST was not run"
    )))
  }
  if (!nzchar(Sys.which("blastp"))) {
    return(list(hits = hits, status = .dnmb_rebasefinder_status_row(
      "supplemental_rebase_blast", "skipped", "blastp not found"
    )))
  }

  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  query_ids <- unique(.dnmb_module_clean_annotation_key(hits$query[mask]))
  q_idx <- match(query_ids, genes$locus_tag)
  query_tbl <- genes[q_idx[!is.na(q_idx)], , drop = FALSE]
  query_tbl$translation <- .dnmb_normalize_translation(query_tbl$translation)
  query_tbl <- query_tbl[!is.na(query_tbl$translation) & nchar(query_tbl$translation) >= blast_min_length, , drop = FALSE]
  if (!base::nrow(query_tbl)) {
    return(list(hits = hits, status = .dnmb_rebasefinder_status_row(
      "supplemental_rebase_blast", "skipped", "No supplemental query sequence passed length filter"
    )))
  }

  supp_dir <- base::file.path(output_dir, "supplemental_rebase_blast")
  base::dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  query_fasta <- base::file.path(supp_dir, "supplemental_queries.faa")
  .dnmb_write_protein_fasta(
    data.frame(protein_label = query_tbl$locus_tag, protein_seq = query_tbl$translation, stringsAsFactors = FALSE),
    query_fasta
  )
  raw_file <- base::file.path(supp_dir, "supplemental_blast_results.txt")
  best_file <- base::file.path(supp_dir, "supplemental_best_hits.tsv")

  parse_blast_results <- base::get("parse_blast_results", envir = base::asNamespace("DefenseViz"))
  filter_blast_results <- base::get("filter_blast_results", envir = base::asNamespace("DefenseViz"))

  blast_file <- tryCatch(
    .dnmb_rebasefinder_run_blastp(
      query_fasta,
      rebase_fasta,
      raw_file,
      evalue = 1e-5,
      num_threads = max(1L, as.integer(cpu)),
      max_target_seqs = 10,
      verbose = verbose
    ),
    error = function(e) {
      if (isTRUE(verbose)) message("[REBASEfinder] Supplemental REBASE BLAST failed: ", conditionMessage(e))
      NULL
    }
  )
  if (base::is.null(blast_file) || !base::file.exists(blast_file)) {
    return(list(hits = hits, status = .dnmb_rebasefinder_status_row(
      "supplemental_rebase_blast", "failed", "Supplemental BLAST produced no result file"
    )))
  }
  blast_tbl <- tryCatch(parse_blast_results(blast_file, verbose = FALSE), error = function(e) data.frame())
  blast_tbl <- tryCatch(
    filter_blast_results(
      blast_tbl,
      min_identity = blast_min_identity,
      min_length = blast_min_length,
      max_evalue = 0.001,
      verbose = FALSE
    ),
    error = function(e) data.frame()
  )
  best <- .dnmb_rebasefinder_best_from_blast(blast_tbl, rebase_data, min_identity = blast_min_identity)
  if (!base::nrow(best)) {
    return(list(hits = hits, status = .dnmb_rebasefinder_status_row(
      "supplemental_rebase_blast", "ok",
      base::sprintf("Supplemental BLAST ran for %d candidates but found no hits passing filters", nrow(query_tbl))
    )))
  }
  utils::write.table(best, best_file, sep = "\t", quote = FALSE, row.names = FALSE)

  for (col in c("supplemental_blast_match", "supplemental_blast_identity", "supplemental_blast_evalue",
                "supplemental_blast_bitscore", "supplemental_blast_length")) {
    if (!col %in% names(hits)) hits[[col]] <- .dnmb_na_vector_like(if (grepl("identity|evalue|bitscore|length", col)) numeric() else character(), nrow(hits))
  }
  h_idx <- match(best$query_id, hits$query)
  h_idx <- h_idx[!is.na(h_idx)]
  if (length(h_idx)) {
    b <- best[match(hits$query[h_idx], best$query_id), , drop = FALSE]
    hits$supplemental_blast_match[h_idx] <- b$best_rebase_match
    hits$supplemental_blast_identity[h_idx] <- b$best_match_identity
    hits$supplemental_blast_evalue[h_idx] <- if ("evalue" %in% names(b)) b$evalue else NA_real_
    hits$supplemental_blast_bitscore[h_idx] <- if ("bitscore" %in% names(b)) b$bitscore else NA_real_
    hits$supplemental_blast_length[h_idx] <- if ("length" %in% names(b)) b$length else NA_real_

    replace_blast <- is.na(hits$blast_identity[h_idx])
    replace_idx <- h_idx[replace_blast]
    bb <- b[replace_blast, , drop = FALSE]
    hits$blast_identity[replace_idx] <- bb$best_match_identity
    hits$blast_evalue[replace_idx] <- if ("evalue" %in% names(bb)) bb$evalue else NA_real_
    hits$blast_bitscore[replace_idx] <- if ("bitscore" %in% names(bb)) bb$bitscore else NA_real_
    hits$blast_length[replace_idx] <- if ("length" %in% names(bb)) bb$length else NA_real_
    trusted_override <- !is.na(bb$best_match_identity) & bb$best_match_identity >= 0.50
    hits$hit_label[replace_idx[trusted_override]] <- bb$best_rebase_match[trusted_override]
    known_family <- trusted_override & !is.na(bb$rm_type) & nzchar(bb$rm_type)
    hits$family_id[replace_idx[known_family]] <- bb$rm_type[known_family]
    known_role <- trusted_override & !is.na(bb$subunit) & nzchar(bb$subunit)
    hits$enzyme_role[replace_idx[known_role]] <- bb$subunit[known_role]
    if ("rec_seq" %in% names(bb)) {
      known_rec <- trusted_override & !is.na(bb$rec_seq) & nzchar(bb$rec_seq)
      hits$rec_seq[replace_idx[known_rec]] <- bb$rec_seq[known_rec]
    }
    hits$evidence_mode[replace_idx] <- "operon_context_rebase_blast"
    hits$typing_eligible[replace_idx] <- !is.na(hits$blast_identity[replace_idx]) &
      hits$blast_identity[replace_idx] >= 0.5
    add_support <- base::sprintf(
      "supplemental_REBASE_BLAST=%s; identity=%s%%",
      bb$best_rebase_match,
      round(bb$best_match_identity * 100, 1)
    )
    hits$support[replace_idx] <- base::ifelse(
      is.na(hits$support[replace_idx]) | !nzchar(hits$support[replace_idx]),
      add_support,
      base::paste(hits$support[replace_idx], add_support, sep = "; ")
    )
  }

  list(hits = hits, status = .dnmb_rebasefinder_status_row(
    "supplemental_rebase_blast", "ok",
    base::sprintf("Supplemental REBASE BLAST matched %d/%d R/S/Type IV context candidates", nrow(best), nrow(query_tbl))
  ))
}

.dnmb_rebasefinder_write_augmented_results <- function(output_dir, hits) {
  out_dir <- output_dir
  base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tsv <- base::file.path(out_dir, "DNMB_REBASEfinder_augmented_hits.tsv")
  utils::write.table(hits, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  xlsx <- base::file.path(out_dir, "DNMB_REBASEfinder_augmented_hits.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(openxlsx::write.xlsx(hits, xlsx, overwrite = TRUE), error = function(e) NULL)
  }
  invisible(list(tsv = tsv, xlsx = if (base::file.exists(xlsx)) xlsx else NA_character_))
}

.dnmb_rebasefinder_structure_query_keep <- function(hits) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(logical())
  base::rep(TRUE, base::nrow(hits))
}

.dnmb_rebasefinder_structure_query_priority <- function(hits) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(character())
  has_blast <- if ("blast_identity" %in% base::names(hits)) !base::is.na(hits$blast_identity) else base::rep(FALSE, base::nrow(hits))
  evidence <- if ("evidence_mode" %in% base::names(hits)) {
    base::as.character(hits$evidence_mode)
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  family <- if ("family_id" %in% base::names(hits)) {
    base::as.character(hits$family_id)
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  typeiii_context <- if ("typeiii_context_status" %in% base::names(hits)) {
    !base::is.na(hits$typeiii_context_status) & base::nzchar(base::as.character(hits$typeiii_context_status))
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  structure_evidence <- evidence %in% c("structure_supported", "structure_only")
  operon_context <- evidence %in% "operon_context"
  base::ifelse(
    structure_evidence, "foldseek_supported",
    base::ifelse(
      operon_context, "operon_context",
      base::ifelse(
        has_blast, "rebase_blast",
        base::ifelse(family %in% "Type III" | typeiii_context, "typeIII_context", "annotation_candidate")
      )
    )
  )
}

.dnmb_rebasefinder_write_structure_queries <- function(output_dir, genes, hits) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(genes) || !base::nrow(hits) || !"translation" %in% base::names(genes)) {
    return(NA_character_)
  }
  hits <- hits[.dnmb_rebasefinder_structure_query_keep(hits), , drop = FALSE]
  if (!base::nrow(hits)) return(NA_character_)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  idx <- base::match(hits$query, genes$locus_tag)
  keep <- !base::is.na(idx)
  if (!base::any(keep)) return(NA_character_)
  out_tbl <- genes[idx[keep], , drop = FALSE]
  hit_tbl <- hits[keep, , drop = FALSE]
  seqs <- .dnmb_normalize_translation(out_tbl$translation)
  keep_seq <- !base::is.na(seqs) & base::nzchar(seqs)
  if (!base::any(keep_seq)) return(NA_character_)
  out_tbl <- out_tbl[keep_seq, , drop = FALSE]
  hit_tbl <- hit_tbl[keep_seq, , drop = FALSE]
  seqs <- seqs[keep_seq]
  priority <- .dnmb_rebasefinder_structure_query_priority(hit_tbl)

  fasta_path <- base::file.path(output_dir, "DNMB_REBASEfinder_structure_queries.faa")
  con <- base::file(fasta_path, "w")
  on.exit(base::close(con), add = TRUE)
  for (i in base::seq_along(seqs)) {
    label <- base::paste(
      c(
        out_tbl$locus_tag[[i]],
        base::paste0("rm_type=", hit_tbl$family_id[[i]]),
        base::paste0("role=", hit_tbl$enzyme_role[[i]]),
        base::paste0("priority=", priority[[i]]),
        base::paste0("hit=", hit_tbl$hit_label[[i]])
      ),
      collapse = " "
    )
    label <- base::gsub("[\r\n\t]+", " ", label)
    base::writeLines(base::paste0(">", label), con)
    base::writeLines(base::gsub("(.{1,80})", "\\1\n", seqs[[i]], perl = TRUE), con)
  }
  base::normalizePath(fasta_path, winslash = "/", mustWork = FALSE)
}

.dnmb_rebasefinder_read_fasta_records <- function(path) {
  if (base::is.null(path) || base::is.na(path) || !base::file.exists(path)) {
    return(base::data.frame(query = character(), header = character(), sequence = character(), stringsAsFactors = FALSE))
  }
  lines <- base::readLines(path, warn = FALSE)
  if (!base::length(lines)) {
    return(base::data.frame(query = character(), header = character(), sequence = character(), stringsAsFactors = FALSE))
  }
  records <- list()
  header <- NULL
  seq_parts <- character()
  flush <- function() {
    if (base::is.null(header)) return(NULL)
    h <- base::sub("^>", "", header)
    q <- base::sub("\\s.*$", "", h)
    base::data.frame(query = q, header = h, sequence = base::paste(seq_parts, collapse = ""), stringsAsFactors = FALSE)
  }
  for (line in lines) {
    if (base::grepl("^>", line)) {
      rec <- flush()
      if (!base::is.null(rec)) records[[base::length(records) + 1L]] <- rec
      header <- line
      seq_parts <- character()
    } else if (base::nzchar(base::trimws(line))) {
      seq_parts <- c(seq_parts, base::trimws(line))
    }
  }
  rec <- flush()
  if (!base::is.null(rec)) records[[base::length(records) + 1L]] <- rec
  if (!base::length(records)) {
    return(base::data.frame(query = character(), header = character(), sequence = character(), stringsAsFactors = FALSE))
  }
  out <- base::do.call(base::rbind, records)
  out$query <- .dnmb_module_clean_annotation_key(out$query)
  out
}

.dnmb_rebasefinder_write_fasta_records <- function(records, path) {
  records <- base::as.data.frame(records, stringsAsFactors = FALSE)
  if (!base::nrow(records)) return(NA_character_)
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- base::file(path, "w")
  on.exit(base::close(con), add = TRUE)
  for (i in base::seq_len(base::nrow(records))) {
    header <- if ("header" %in% base::names(records) && !base::is.na(records$header[[i]]) && base::nzchar(records$header[[i]])) {
      records$header[[i]]
    } else {
      records$query[[i]]
    }
    seq <- base::gsub("\\s+", "", base::as.character(records$sequence[[i]]), perl = TRUE)
    if (base::is.na(seq) || !base::nzchar(seq)) next
    base::writeLines(base::paste0(">", header), con)
    base::writeLines(base::gsub("(.{1,80})", "\\1\n", seq, perl = TRUE), con)
  }
  base::normalizePath(path, winslash = "/", mustWork = FALSE)
}

.dnmb_rebasefinder_structure_dirs <- function(output_dir, structure_validation_path = NULL) {
  dirs <- c(
    base::file.path(output_dir, "query_structures"),
    base::file.path(output_dir, "alphafold_query_structures"),
    base::file.path(output_dir, "rebasefinder_query_structures"),
    base::file.path(output_dir, "dnmb_module_rebasefinder", "query_structures"),
    base::file.path(output_dir, "dnmb_module_rebasefinder", "alphafold_query_structures")
  )
  if (!base::is.null(structure_validation_path) && !base::is.na(structure_validation_path) && base::nzchar(structure_validation_path)) {
    root <- base::dirname(structure_validation_path)
    dirs <- c(
      dirs,
      base::file.path(root, "query_structures"),
      base::file.path(root, "alphafold_query_structures"),
      base::file.path(root, "rebasefinder_query_structures")
    )
  }
  base::unique(dirs[base::dir.exists(dirs)])
}

.dnmb_rebasefinder_structure_file_for_query <- function(query, dirs) {
  if (!base::length(dirs) || base::is.na(query) || !base::nzchar(query)) return(NA_character_)
  safe_query <- base::gsub("[^A-Za-z0-9_.-]+", "_", query)
  candidates <- base::as.vector(base::outer(
    dirs,
    base::paste0(base::rep(c(query, safe_query), each = 2L), c(".pdb", ".cif")),
    base::file.path
  ))
  candidates <- base::unique(candidates)
  hit <- candidates[base::file.exists(candidates)]
  if (!base::length(hit)) return(NA_character_)
  base::normalizePath(hit[[1]], winslash = "/", mustWork = FALSE)
}

.dnmb_rebasefinder_foldseek_hit_queries <- function(path) {
  if (base::is.null(path) || base::is.na(path) || !base::file.exists(path)) return(character())
  tbl <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
  if (base::is.null(tbl) || !base::nrow(tbl)) return(character())
  q_col <- base::intersect(c("query", "qseqid", "Query", "query_id"), base::names(tbl))[1]
  if (base::is.na(q_col)) q_col <- base::names(tbl)[[1]]
  base::unique(.dnmb_module_clean_annotation_key(base::as.character(tbl[[q_col]])))
}

.dnmb_rebasefinder_write_structure_coverage <- function(output_dir,
                                                        genes,
                                                        hits,
                                                        fasta_path = NULL,
                                                        structure_validation_path = NULL,
                                                        structure_tbl = NULL) {
  if (base::is.null(fasta_path) || base::is.na(fasta_path) || !base::file.exists(fasta_path)) {
    fasta_path <- .dnmb_rebasefinder_write_structure_queries(output_dir, genes, hits)
  }
  records <- .dnmb_rebasefinder_read_fasta_records(fasta_path)
  if (!base::nrow(records)) return(NULL)

  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if ("query" %in% base::names(hits)) hits$query <- .dnmb_module_clean_annotation_key(hits$query)
  if ("locus_tag" %in% base::names(genes)) genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)

  hidx <- if ("query" %in% base::names(hits)) base::match(records$query, hits$query) else base::rep(NA_integer_, base::nrow(records))
  gidx <- if ("locus_tag" %in% base::names(genes)) base::match(records$query, genes$locus_tag) else base::rep(NA_integer_, base::nrow(records))
  hit_rows <- hits[hidx, , drop = FALSE]
  gene_rows <- genes[gidx, , drop = FALSE]
  priority <- .dnmb_rebasefinder_structure_query_priority(hit_rows)
  if (base::length(priority) != base::nrow(records)) priority <- base::rep(NA_character_, base::nrow(records))

  if (base::is.null(structure_validation_path) || base::is.na(structure_validation_path)) {
    structure_validation_path <- .dnmb_rebasefinder_structure_path(output_dir)
  }
  structure_dirs <- .dnmb_rebasefinder_structure_dirs(output_dir, structure_validation_path)
  structure_file <- base::vapply(records$query, .dnmb_rebasefinder_structure_file_for_query, character(1), dirs = structure_dirs)
  structure_file_exists <- !base::is.na(structure_file) & base::nzchar(structure_file) & base::file.exists(structure_file)
  foldseek_hit_queries <- .dnmb_rebasefinder_foldseek_hit_queries(structure_validation_path)
  foldseek_hit_present <- records$query %in% foldseek_hit_queries

	  raw_supported <- base::rep(FALSE, base::nrow(records))
	  supported <- base::rep(FALSE, base::nrow(records))
	  best_status <- base::rep(NA_character_, base::nrow(records))
	  best_hit <- base::rep(NA_character_, base::nrow(records))
	  best_family <- base::rep(NA_character_, base::nrow(records))
	  best_role <- base::rep(NA_character_, base::nrow(records))
	  best_chain_role <- base::rep(NA_character_, base::nrow(records))
	  if (!base::is.null(structure_tbl) && base::is.data.frame(structure_tbl) && base::nrow(structure_tbl) && "query" %in% base::names(structure_tbl)) {
	    structure_tbl$query <- .dnmb_module_clean_annotation_key(structure_tbl$query)
	    sidx <- base::match(records$query, structure_tbl$query)
	    ok <- !base::is.na(sidx)
	    if (base::any(ok)) {
	      if ("structure_pass" %in% base::names(structure_tbl)) {
	        raw_supported[ok] <- !base::is.na(structure_tbl$structure_pass[sidx[ok]]) & structure_tbl$structure_pass[sidx[ok]]
	      }
	      if ("structure_status" %in% base::names(structure_tbl)) {
	        best_status[ok] <- base::as.character(structure_tbl$structure_status[sidx[ok]])
	      }
	      if ("structure_hit" %in% base::names(structure_tbl)) {
	        best_hit[ok] <- base::as.character(structure_tbl$structure_hit[sidx[ok]])
	      }
	      if ("structure_family" %in% base::names(structure_tbl)) {
	        best_family[ok] <- base::as.character(structure_tbl$structure_family[sidx[ok]])
	      }
	      if ("structure_role" %in% base::names(structure_tbl)) {
	        best_role[ok] <- base::as.character(structure_tbl$structure_role[sidx[ok]])
	      }
	      if ("structure_chain_role" %in% base::names(structure_tbl)) {
	        best_chain_role[ok] <- base::as.character(structure_tbl$structure_chain_role[sidx[ok]])
	      }
	    }
	  }

  value_or_na <- function(tbl, col, idx) {
    if (!col %in% base::names(tbl)) return(base::rep(NA_character_, base::length(idx)))
    base::as.character(tbl[[col]][idx])
  }
	  num_value_or_na <- function(tbl, col, idx) {
	    if (!col %in% base::names(tbl)) return(base::rep(NA_real_, base::length(idx)))
	    suppressWarnings(base::as.numeric(tbl[[col]][idx]))
	  }

	  structure_consistency <- .dnmb_rebasefinder_structure_candidate_consistency(
	    value_or_na(hits, "family_id", hidx),
	    value_or_na(hits, "enzyme_role", hidx),
	    best_family,
	    best_role,
	    best_chain_role
	  )
	  supported <- raw_supported & structure_consistency$structure_candidate_consistent
	  mismatch <- foldseek_hit_present & raw_supported & !supported
	  best_status[mismatch & !structure_consistency$structure_family_consistent] <- "structure_family_mismatch"
	  best_status[mismatch & structure_consistency$structure_family_consistent &
	                !structure_consistency$structure_role_consistent] <- "structure_role_mismatch"

	  coverage_status <- base::ifelse(
	    supported, "foldseek_supported",
	    base::ifelse(
	      mismatch & !structure_consistency$structure_family_consistent, "foldseek_family_mismatch",
	      base::ifelse(
	        mismatch & !structure_consistency$structure_role_consistent, "foldseek_role_mismatch",
	        base::ifelse(
	          foldseek_hit_present, "foldseek_checked_no_pass",
	          base::ifelse(
	            structure_file_exists, "structure_available_no_foldseek_hit",
	            "structure_missing"
	          )
	        )
	      )
	    )
	  )

  out <- base::data.frame(
    query = records$query,
    family_id = value_or_na(hits, "family_id", hidx),
    enzyme_role = value_or_na(hits, "enzyme_role", hidx),
    evidence_mode = value_or_na(hits, "evidence_mode", hidx),
    hit_label = value_or_na(hits, "hit_label", hidx),
    priority = priority,
    aa_len = base::nchar(records$sequence),
    gene_start = num_value_or_na(genes, "start", gidx),
    gene_end = num_value_or_na(genes, "end", gidx),
    structure_file = structure_file,
    structure_file_exists = structure_file_exists,
	    foldseek_hit_present = foldseek_hit_present,
	    structure_supported = supported,
	    structure_status = best_status,
	    structure_best_hit = best_hit,
	    structure_family_consistent = structure_consistency$structure_family_consistent,
	    structure_role_consistent = structure_consistency$structure_role_consistent,
	    structure_candidate_consistent = structure_consistency$structure_candidate_consistent,
	    coverage_status = coverage_status,
	    stringsAsFactors = FALSE
	  )
  out <- out[base::order(out$coverage_status, out$family_id, out$gene_start, out$query), , drop = FALSE]

  tsv <- base::file.path(output_dir, "DNMB_REBASEfinder_structure_coverage.tsv")
  utils::write.table(out, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  xlsx <- base::file.path(output_dir, "DNMB_REBASEfinder_structure_coverage.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(openxlsx::write.xlsx(out, xlsx, overwrite = TRUE), error = function(e) NULL)
  }

  missing <- records[records$query %in% out$query[out$coverage_status == "structure_missing"], , drop = FALSE]
  unchecked <- records[records$query %in% out$query[!out$foldseek_hit_present], , drop = FALSE]
  missing_faa <- .dnmb_rebasefinder_write_fasta_records(
    missing,
    base::file.path(output_dir, "DNMB_REBASEfinder_structure_missing_queries.faa")
  )
  unchecked_faa <- .dnmb_rebasefinder_write_fasta_records(
    unchecked,
    base::file.path(output_dir, "DNMB_REBASEfinder_structure_unchecked_queries.faa")
  )

	  list(
	    tsv = base::normalizePath(tsv, winslash = "/", mustWork = FALSE),
	    xlsx = if (base::file.exists(xlsx)) base::normalizePath(xlsx, winslash = "/", mustWork = FALSE) else NA_character_,
	    missing_faa = missing_faa,
	    unchecked_faa = unchecked_faa,
	    table = out,
	    n_queries = base::nrow(out),
	    n_structure_files = base::sum(out$structure_file_exists, na.rm = TRUE),
	    n_foldseek_hits = base::sum(out$foldseek_hit_present, na.rm = TRUE),
    n_supported = base::sum(out$structure_supported, na.rm = TRUE),
    n_missing = base::sum(out$coverage_status == "structure_missing", na.rm = TRUE),
    n_unchecked = base::sum(!out$foldseek_hit_present, na.rm = TRUE)
	  )
	}

.dnmb_rebasefinder_apply_structure_coverage <- function(hits, coverage_tbl) {
  if (base::is.null(coverage_tbl) || !base::is.data.frame(coverage_tbl) ||
      !base::nrow(coverage_tbl) || !"query" %in% base::names(coverage_tbl)) {
    return(hits)
  }
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits) || !"query" %in% base::names(hits)) return(hits)
  hits$query <- .dnmb_module_clean_annotation_key(hits$query)
  coverage_tbl$query <- .dnmb_module_clean_annotation_key(coverage_tbl$query)
  idx <- base::match(hits$query, coverage_tbl$query)
  ok <- !base::is.na(idx)
  if (!base::any(ok)) return(hits)

  ensure_col <- function(col, value) {
    if (!col %in% base::names(hits)) hits[[col]] <<- base::rep(value, base::nrow(hits))
  }
  ensure_col("structure_status", NA_character_)
  ensure_col("structure_pass", FALSE)
  ensure_col("structure_file", NA_character_)
  ensure_col("structure_file_exists", FALSE)
  ensure_col("foldseek_hit_present", FALSE)
  ensure_col("structure_supported", FALSE)
  ensure_col("structure_coverage_status", NA_character_)

  cov <- coverage_tbl[idx[ok], , drop = FALSE]
  cov_status <- base::as.character(cov$coverage_status)
	  mapped_status <- base::ifelse(
	    cov_status == "foldseek_supported", "structure_supported",
	    base::ifelse(
	      cov_status == "foldseek_role_mismatch", "structure_role_mismatch",
	      base::ifelse(
	        cov_status == "foldseek_family_mismatch", "structure_family_mismatch",
	        base::ifelse(
	      cov_status == "foldseek_checked_no_pass", "structure_checked_no_pass",
	      base::ifelse(
	        cov_status == "structure_available_no_foldseek_hit", "structure_available_no_foldseek_hit",
	        "structure_missing"
	      )
	        )
	      )
	    )
	  )

	  status_missing <- base::is.na(hits$structure_status[ok]) | !base::nzchar(base::as.character(hits$structure_status[ok]))
	  repl_idx <- base::which(ok)[status_missing]
	  hits$structure_status[repl_idx] <- mapped_status[status_missing]

	  pass <- cov_status == "foldseek_supported"
	  pass_idx <- base::which(ok)
	  if (base::any(pass)) hits$structure_status[pass_idx[pass]] <- "structure_supported"
	  not_supported <- !pass
	  if (base::any(not_supported)) hits$structure_status[pass_idx[not_supported]] <- mapped_status[not_supported]
	  hits$structure_pass[pass_idx] <- pass
	  hits$structure_supported[pass_idx] <- pass
	  hits$structure_coverage_status[pass_idx] <- cov_status

	  if ("structure_file" %in% base::names(cov)) hits$structure_file[pass_idx] <- base::as.character(cov$structure_file)
	  if ("structure_file_exists" %in% base::names(cov)) hits$structure_file_exists[pass_idx] <- cov$structure_file_exists
	  if ("foldseek_hit_present" %in% base::names(cov)) hits$foldseek_hit_present[pass_idx] <- cov$foldseek_hit_present
	  if ("structure_family_consistent" %in% base::names(cov)) hits$structure_family_consistent[pass_idx] <- cov$structure_family_consistent
	  if ("structure_role_consistent" %in% base::names(cov)) hits$structure_role_consistent[pass_idx] <- cov$structure_role_consistent
	  if ("structure_candidate_consistent" %in% base::names(cov)) hits$structure_candidate_consistent[pass_idx] <- cov$structure_candidate_consistent

	  hits
	}

.dnmb_rebasefinder_check_package <- function() {
  if (!requireNamespace("DefenseViz", quietly = TRUE)) {
    return(FALSE)
  }
  TRUE
}

.dnmb_rebasefinder_install_package <- function(verbose = TRUE) {
  if (.dnmb_rebasefinder_check_package()) {
    return(TRUE)
  }
  if (verbose) message("[REBASEfinder] Installing DefenseViz from GitHub...")
  tryCatch({
    if (!requireNamespace("devtools", quietly = TRUE)) {
      utils::install.packages("devtools", repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    devtools::install_github("JAEYOONSUNG/DefenseViz", quiet = TRUE, upgrade = "never")
    .dnmb_rebasefinder_check_package()
  }, error = function(e) {
    if (verbose) message("[REBASEfinder] DefenseViz install failed: ", conditionMessage(e))
    FALSE
  })
}

.dnmb_rebasefinder_prepare_input <- function(genes, output_dir) {
  genes_df <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  required <- c("locus_tag", "start", "end", "direction", "translation", "product")
  missing <- base::setdiff(required, base::names(genes_df))
  if (base::length(missing)) {
    base::stop(
      "genbank_table is missing columns required by REBASEfinder: ",
      base::paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  tsv_path <- base::file.path(output_dir, "rebasefinder_input.tsv")
  utils::write.table(genes_df, file = tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
  tsv_path
}

.dnmb_rebasefinder_normalize_hits <- function(rm_comprehensive, id_col = "locus_tag") {
  if (base::is.null(rm_comprehensive) || !base::is.data.frame(rm_comprehensive) || !base::nrow(rm_comprehensive)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  tbl <- base::as.data.frame(rm_comprehensive, stringsAsFactors = FALSE)

  if (!id_col %in% base::names(tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  tbl$query <- .dnmb_module_clean_annotation_key(tbl[[id_col]])

  tbl$source <- "rebasefinder"
  tbl$family_system <- "REBASEfinder"

  tbl$family_id <- base::ifelse(
    !base::is.na(tbl[["rm_type"]]) & base::nzchar(base::as.character(tbl[["rm_type"]])),
    base::as.character(tbl[["rm_type"]]),
    NA_character_
  )

  tbl$hit_label <- base::ifelse(
    !base::is.na(tbl[["blast_match"]]) & base::nzchar(base::as.character(tbl[["blast_match"]])),
    base::as.character(tbl[["blast_match"]]),
    NA_character_
  )

  tbl$enzyme_role <- base::ifelse(
    !base::is.na(tbl[["subunit"]]) & base::nzchar(base::as.character(tbl[["subunit"]])),
    base::as.character(tbl[["subunit"]]),
    NA_character_
  )
  # Override role from REBASE hit name when subunit is missing or contradicts
  # the naming convention: M./M1./M2. = methyltransferase, not restriction
  if ("blast_match" %in% base::names(tbl)) {
    hit_name <- base::as.character(tbl[["blast_match"]])
    inferred <- .dnmb_rebasefinder_role_from_hit(hit_name)
    needs_fix <- (!base::is.na(inferred)) &
      (base::is.na(tbl$enzyme_role) | tbl$enzyme_role != inferred)
    tbl$enzyme_role[needs_fix] <- inferred[needs_fix]
  }

  tbl$evidence_mode <- base::ifelse(
    !base::is.na(tbl[["blast_identity"]]) & tbl[["blast_identity"]] >= 0.5,
    "high_confidence",
    base::ifelse(
      !base::is.na(tbl[["blast_identity"]]),
      "low_confidence",
      "annotation_only"
    )
  )

  tbl$substrate_label <- base::ifelse(
    !base::is.na(tbl[["rec_seq"]]) & base::nzchar(base::as.character(tbl[["rec_seq"]])) & tbl[["rec_seq"]] != "?",
    base::as.character(tbl[["rec_seq"]]),
    NA_character_
  )

  support_parts <- base::vapply(base::seq_len(base::nrow(tbl)), function(i) {
    parts <- c(
      if (!base::is.na(tbl[["rm_type"]][[i]]) && base::nzchar(tbl[["rm_type"]][[i]])) base::paste0("rm_type=", tbl[["rm_type"]][[i]]) else NULL,
      if (!base::is.na(tbl[["subunit"]][[i]]) && base::nzchar(tbl[["subunit"]][[i]])) base::paste0("subunit=", tbl[["subunit"]][[i]]) else NULL,
      if (!base::is.na(tbl[["blast_match"]][[i]]) && base::nzchar(tbl[["blast_match"]][[i]])) base::paste0("rebase_match=", tbl[["blast_match"]][[i]]) else NULL,
      if (!base::is.na(tbl[["blast_identity"]][[i]])) base::paste0("identity=", base::round(tbl[["blast_identity"]][[i]] * 100, 1), "%") else NULL,
      if (!base::is.na(tbl[["rec_seq"]][[i]]) && base::nzchar(tbl[["rec_seq"]][[i]]) && tbl[["rec_seq"]][[i]] != "?") base::paste0("rec_seq=", tbl[["rec_seq"]][[i]]) else NULL,
      if (!base::is.na(tbl[["operon_id"]][[i]]) && base::nzchar(tbl[["operon_id"]][[i]])) base::paste0("operon=", tbl[["operon_id"]][[i]]) else NULL
    )
    if (!base::length(parts)) return(NA_character_)
    base::paste(parts, collapse = "; ")
  }, character(1))
  tbl$support <- support_parts

  tbl$typing_eligible <- !base::is.na(tbl[["blast_identity"]]) & tbl[["blast_identity"]] >= 0.5

  extra_cols <- c("blast_identity", "blast_evalue", "blast_bitscore", "blast_length",
                  "rec_seq", "operon_id", "partner_locus_tag")
  extra_cols <- base::intersect(extra_cols, base::names(tbl))

  out <- base::data.frame(
    query = tbl$query,
    source = tbl$source,
    family_system = tbl$family_system,
    family_id = tbl$family_id,
    hit_label = tbl$hit_label,
    enzyme_role = tbl$enzyme_role,
    evidence_mode = tbl$evidence_mode,
    substrate_label = tbl$substrate_label,
    support = tbl$support,
    typing_eligible = tbl$typing_eligible,
    stringsAsFactors = FALSE
  )

  for (col in extra_cols) {
    out[[col]] <- tbl[[col]]
  }

  out <- out[!base::is.na(out$query), , drop = FALSE]
  out <- out[base::order(out$query, -base::ifelse(base::is.na(out$blast_identity), -Inf, out$blast_identity)), , drop = FALSE]
  out <- out[!base::duplicated(out$query), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_rebasefinder_output_table <- function(genes, hits) {
  selected_hits <- if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    .dnmb_module_empty_optional_long_table()
  } else {
    hits
  }
  .dnmb_module_output_table(genes = genes, hits = selected_hits)
}

dnmb_run_rebasefinder_module <- function(genes,
                                         output_dir,
                                         genbank = NULL,
                                         compare_rebase = TRUE,
                                         search_motifs = TRUE,
                                         blast_min_identity = 0.10,
                                         blast_min_length = 50,
                                         max_operon_gap = 5000,
                                         max_intervening = 1,
                                         typeiii_context_neighbors = 3L,
                                         structure_validation = NULL,
                                         structure_max_evalue = 1e-3,
                                         structure_min_probability = 0.50,
                                         structure_min_tmscore = 0.45,
                                         install = TRUE,
                                         cache_root = NULL,
                                         cpu = 1L,
                                         verbose = TRUE) {

  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_rebasefinder_empty_status()

  if (!.dnmb_rebasefinder_check_package()) {
    if (base::isTRUE(install)) {
      installed <- .dnmb_rebasefinder_install_package(verbose = verbose)
      if (!installed) {
        status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
          "rebasefinder_install", "failed",
          "DefenseViz package not available. Install with: devtools::install_github('JAEYOONSUNG/DefenseViz')"
        ))
        return(list(
          ok = FALSE,
          status = status,
          hits = .dnmb_module_empty_optional_long_table(),
          output_table = base::data.frame()
        ))
      }
    } else {
      status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
        "rebasefinder_install", "skipped",
        "DefenseViz not installed and install=FALSE"
      ))
      return(list(
        ok = FALSE,
        status = status,
        hits = .dnmb_module_empty_optional_long_table(),
        output_table = base::data.frame()
      ))
    }
  }

  status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
    "rebasefinder_install", "ok", "DefenseViz available"
  ))

  # Override DefenseViz cache dir -> DNMB's db_modules/rebasefinder/cache.
  # This works regardless of DefenseViz version (no env var dependency).
  .dnmb_rebasefinder_set_defenseviz_cache(cache_root = cache_root, verbose = verbose)

  input_path <- .dnmb_rebasefinder_prepare_input(genes, output_dir)

  # Patch detect_methyltransferases: deduplicate after union to prevent
  # row explosion (keyword fallback can produce duplicates via left_join)
  base::tryCatch({
    orig_detect <- DefenseViz:::detect_methyltransferases
    patched_detect <- function(dnmb_data, pfam_col = NULL, product_col = "product") {
      result <- orig_detect(dnmb_data, pfam_col = pfam_col, product_col = product_col)
      id_col <- base::intersect(c("locus_tag", "protein_id"), base::names(result))[1]
      if (!base::is.null(id_col) && base::nrow(result) > 0) {
        result <- result[!base::duplicated(result[[id_col]]), , drop = FALSE]
      }
      # Cap keyword-only candidates to prevent OOM in operon analysis
      if ("detection_method" %in% base::names(result) && base::nrow(result) > 500L) {
        keep <- result$detection_method != "keyword_only" | result$mtase_confidence %in% c("high", "medium")
        if (base::sum(keep) > 0) result <- result[keep, , drop = FALSE]
      }
      result
    }
    utils::assignInNamespace("detect_methyltransferases", patched_detect, ns = "DefenseViz")
    if (isTRUE(verbose)) message("[REBASEfinder] Patched detect_methyltransferases (dedup + cap)")
  }, error = function(e) {
    if (isTRUE(verbose)) message("[REBASEfinder] Could not patch detect_methyltransferases: ", conditionMessage(e))
  })

  pipeline_result <- base::tryCatch({
    DefenseViz::rmscan_pipeline(
      input_file = input_path,
      output_dir = output_dir,
      compare_rebase = compare_rebase,
      search_motifs = search_motifs,
      blast_min_identity = blast_min_identity,
      blast_min_length = blast_min_length,
      max_operon_gap = max_operon_gap,
      max_intervening = max_intervening,
      verbose = verbose,
      save_intermediates = TRUE
    )
  }, error = function(e) {
    list(.error = conditionMessage(e))
  })

  if (!base::is.null(pipeline_result$.error)) {
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "rebasefinder_run", "failed", pipeline_result$.error
    ))
    return(list(
      ok = FALSE,
      status = status,
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = base::data.frame()
    ))
  }

  rm_comprehensive <- pipeline_result$rm_comprehensive
  if (base::is.null(rm_comprehensive) || !base::is.data.frame(rm_comprehensive)) {
    rm_comprehensive <- base::data.frame()
  }

  hits <- .dnmb_rebasefinder_normalize_hits(rm_comprehensive)
  structure_path <- .dnmb_rebasefinder_structure_path(output_dir, structure_validation = structure_validation)
  structure_tbl <- .dnmb_rebasefinder_read_structure_validation(
    structure_path,
    max_evalue = structure_max_evalue,
    min_probability = structure_min_probability,
    min_tmscore = structure_min_tmscore
  )
  if (!base::is.null(structure_tbl) && base::nrow(structure_tbl)) {
    hits <- .dnmb_rebasefinder_merge_structure_validation(hits, genes, structure_tbl)
    n_struct <- base::sum(structure_tbl$structure_pass, na.rm = TRUE)
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "structure_validation", "ok",
      base::sprintf("%d Foldseek/structure-supported candidates from %s", n_struct, basename(structure_path))
    ))
  } else if (!base::is.null(structure_validation)) {
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "structure_validation", "missing",
      "No readable Foldseek/structure validation TSV found"
    ))
  } else {
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "structure_validation", "not_run",
      "No Foldseek/structure TSV found; pass rebasefinder_structure_tsv or place foldseek_results.tsv in the REBASEfinder output directory"
    ))
  }

  hits <- .dnmb_rebasefinder_add_typeiii_context(
    hits,
    genes,
    max_operon_gap = max_operon_gap,
    max_intervening = max_intervening,
    max_neighbors = typeiii_context_neighbors
  )
  hits <- .dnmb_rebasefinder_add_typei_context(
    hits,
    genes,
    max_operon_gap = max_operon_gap,
    max_intervening = max_intervening,
    max_neighbors = max(4L, as.integer(typeiii_context_neighbors))
  )
  hits <- .dnmb_rebasefinder_add_typeiv_candidates(hits, genes)
  supplemental <- .dnmb_rebasefinder_supplemental_rebase_blast(
    hits,
    genes,
    output_dir = output_dir,
    blast_min_identity = blast_min_identity,
    blast_min_length = blast_min_length,
    cache_root = cache_root,
    cpu = cpu,
    verbose = verbose
  )
  hits <- supplemental$hits
  if (base::nrow(supplemental$status)) {
    status <- dplyr::bind_rows(status, supplemental$status)
	  }
	  hits <- .dnmb_rebasefinder_add_sequence_partial_status(hits, genes)
	  augmented_files <- list()
	  augmented_files$structure_queries_faa <- .dnmb_rebasefinder_write_structure_queries(output_dir, genes, hits)
	  structure_coverage <- .dnmb_rebasefinder_write_structure_coverage(
	    output_dir,
	    genes,
    hits,
    fasta_path = augmented_files$structure_queries_faa,
	    structure_validation_path = structure_path,
	    structure_tbl = structure_tbl
	  )
	  if (!base::is.null(structure_coverage)) {
	    hits <- .dnmb_rebasefinder_apply_structure_coverage(hits, structure_coverage$table)
	    augmented_files$structure_coverage_tsv <- structure_coverage$tsv
	    augmented_files$structure_coverage_xlsx <- structure_coverage$xlsx
	    augmented_files$structure_missing_queries_faa <- structure_coverage$missing_faa
    augmented_files$structure_unchecked_queries_faa <- structure_coverage$unchecked_faa
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "structure_coverage",
      if (structure_coverage$n_unchecked > 0L || structure_coverage$n_missing > 0L) "partial" else "ok",
      base::sprintf(
        "Foldseek coverage: %d/%d queries have structure files; %d/%d have Foldseek hits; %d supported. Missing FASTA: %s",
        structure_coverage$n_structure_files,
        structure_coverage$n_queries,
        structure_coverage$n_foldseek_hits,
        structure_coverage$n_queries,
        structure_coverage$n_supported,
        if (!base::is.na(structure_coverage$missing_faa)) base::basename(structure_coverage$missing_faa) else "none"
	      )
	    ))
	  }
	  augmented_result_files <- .dnmb_rebasefinder_write_augmented_results(output_dir, hits)
	  augmented_files <- c(augmented_result_files, augmented_files)
	  if (!base::is.na(augmented_files$structure_queries_faa) &&
	      base::file.exists(augmented_files$structure_queries_faa)) {
    n_structure_queries <- base::length(base::grep(
      "^>",
      base::readLines(augmented_files$structure_queries_faa, warn = FALSE)
    ))
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "structure_queries", "ok",
      base::sprintf(
        "Candidate FASTA for structure prediction: %s (%d REBASEfinder candidates with translations)",
        base::basename(augmented_files$structure_queries_faa),
        n_structure_queries
      )
    ))
  }
  output_table <- .dnmb_rebasefinder_output_table(genes = genes, hits = hits)

  n_hits <- base::nrow(hits)
  n_high <- base::sum(hits$typing_eligible, na.rm = TRUE)
  n_typeiii <- base::sum(hits$family_id == "Type III", na.rm = TRUE)
  n_typeiii_context <- if ("typeiii_operon_supported" %in% base::names(hits)) {
    base::sum(hits$typeiii_operon_supported %in% TRUE, na.rm = TRUE)
  } else 0L
  n_typei <- base::sum(hits$family_id == "Type I", na.rm = TRUE)
  n_typei_context <- if ("typei_operon_supported" %in% base::names(hits)) {
    base::sum(hits$typei_operon_supported %in% TRUE, na.rm = TRUE)
  } else 0L
  n_typeiv <- base::sum(hits$family_id == "Type IV", na.rm = TRUE)

  status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
    "rebasefinder_run", "ok",
    base::sprintf(
      "%d R-M candidates (%d high-confidence; %d Type I/%d Type I operon-supported; %d Type III/%d Type III operon-supported; %d Type IV)",
      n_hits, n_high, n_typei, n_typei_context, n_typeiii, n_typeiii_context, n_typeiv
    )
  ))

  list(
    ok = n_hits > 0L || base::nrow(rm_comprehensive) == 0L,
    status = status,
    files = augmented_files,
    hits = hits,
    output_table = output_table,
    pipeline_result = pipeline_result
  )
}
