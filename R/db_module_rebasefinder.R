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

.dnmb_rebasefinder_download_bairoch <- function(cache_root = NULL,
                                                force = FALSE,
                                                max_age_days = 31L) {
  cache_path <- .dnmb_rebasefinder_bairoch_cache_path(cache_root)
  stale_lookup <- NULL
  if (!isTRUE(force) && file.exists(cache_path)) {
    stale_lookup <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    cache_age <- base::as.numeric(base::difftime(
      base::Sys.time(), base::file.info(cache_path)$mtime, units = "days"
    ))
    if (base::is.data.frame(stale_lookup) &&
        identical(base::attr(stale_lookup, "dnmb_bairoch_schema"), 2L) &&
        !base::is.na(cache_age) && cache_age <= max_age_days) {
      return(stale_lookup)
    }
  }
  legacy_path <- file.path(
    .dnmb_db_cache_root(cache_root = cache_root, create = FALSE),
    "db_modules", "rebasefinder", "cache", "rebase_bairoch_lookup.rds"
  )
  if (!isTRUE(force) && file.exists(legacy_path)) {
    lookup <- readRDS(legacy_path)
    legacy_age <- base::as.numeric(base::difftime(
      base::Sys.time(), base::file.info(legacy_path)$mtime, units = "days"
    ))
    if (identical(base::attr(lookup, "dnmb_bairoch_schema"), 2L) &&
        !base::is.na(legacy_age) && legacy_age <= max_age_days) {
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      tryCatch(saveRDS(lookup, cache_path), error = function(e) NULL)
      return(lookup)
    }
    if (base::is.null(stale_lookup)) stale_lookup <- lookup
  }
  url <- "https://rebase.neb.com/rebase/link_bairoch"
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  tryCatch(utils::download.file(url, tmp, quiet = TRUE, method = "libcurl"),
           error = function(e) NULL)
  if (!file.exists(tmp) || file.info(tmp)$size < 1000) {
    return(stale_lookup)
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
      if (!is.na(current$id) && (!is.na(current$rs) || !is.na(current$ms))) {
        entries[[length(entries) + 1L]] <- current
      }
      current <- list(id = NA_character_, rs = NA_character_, ms = NA_character_)
    }
  }
  if (!length(entries)) return(NULL)
  # Parse MS field: "2(m6A),-4(m6A)" → list of {pos, type}
  rows <- lapply(entries, function(e) {
    ms_parsed <- list()
    if (!is.na(e$ms) && nzchar(e$ms)) {
      ms_parts <- strsplit(gsub(";$", "", e$ms), ",")[[1]]
      ms_parsed <- lapply(trimws(ms_parts), function(p) {
        m <- regmatches(p, regexec("^(-?[0-9?]+)\\(([^)]+)\\)", p))[[1]]
        if (length(m) == 3L) list(pos = m[2], type = m[3]) else NULL
      })
      ms_parsed <- Filter(Negate(is.null), ms_parsed)
    }
    # Primary methylation: first entry
    primary_pos <- if (length(ms_parsed)) ms_parsed[[1]]$pos else NA_character_
    primary_type <- if (length(ms_parsed)) ms_parsed[[1]]$type else NA_character_
    data.frame(
      enzyme_name = e$id,
      rec_seq = if (is.na(e$rs)) NA_character_ else e$rs,
      meth_type = primary_type,
      meth_pos = primary_pos,
      meth_all = if (is.na(e$ms) || !nzchar(e$ms)) NA_character_ else e$ms,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  base::attr(out, "dnmb_bairoch_schema") <- 2L
  out
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

.dnmb_rebasefinder_valid_recognition <- function(x) {
  x <- base::trimws(base::as.character(x))
  !base::is.na(x) & base::nzchar(x) & !x %in% c("?", "NA", "-")
}

.dnmb_rebasefinder_clean_rebase_subject <- function(x) {
  x <- base::as.character(x)
  base::sub("_[0-9]+$", "", x, perl = TRUE)
}

.dnmb_rebasefinder_system_key <- function(x) {
  x <- .dnmb_rebasefinder_clean_rebase_subject(x)
  x <- base::sub("^[MRS][0-9]*\\.", "", x, perl = TRUE)
  base::tolower(x)
}

.dnmb_rebasefinder_recognition_lookup <- function(rebase_data = NULL, bairoch = NULL) {
  rows <- list()
  if (base::is.data.frame(rebase_data) && base::nrow(rebase_data) &&
      base::all(c("enzyme_name", "rec_seq") %in% base::names(rebase_data))) {
    tmp <- base::data.frame(
      enzyme_name = base::as.character(rebase_data$enzyme_name),
      rec_seq = base::as.character(rebase_data$rec_seq),
      recognition_source = "rebase_protein_header",
      is_gold_standard = if ("is_gold_standard" %in% base::names(rebase_data)) {
        base::as.logical(rebase_data$is_gold_standard)
      } else {
        FALSE
      },
      stringsAsFactors = FALSE
    )
    rows[[base::length(rows) + 1L]] <- tmp
  }
  if (base::is.data.frame(bairoch) && base::nrow(bairoch) &&
      base::all(c("enzyme_name", "rec_seq") %in% base::names(bairoch))) {
    rows[[base::length(rows) + 1L]] <- base::data.frame(
      enzyme_name = base::as.character(bairoch$enzyme_name),
      rec_seq = base::as.character(bairoch$rec_seq),
      recognition_source = "rebase_bairoch",
      is_gold_standard = TRUE,
      stringsAsFactors = FALSE
    )
  }
  if (!base::length(rows)) {
    return(base::data.frame(
      enzyme_name = character(), rec_seq = character(), recognition_source = character(),
      is_gold_standard = logical(), system_key = character(), stringsAsFactors = FALSE
    ))
  }
  out <- base::do.call(base::rbind, rows)
  out <- out[.dnmb_rebasefinder_valid_recognition(out$rec_seq) &
               !base::is.na(out$enzyme_name) & base::nzchar(out$enzyme_name), , drop = FALSE]
  if (!base::nrow(out)) return(out)
  out$system_key <- .dnmb_rebasefinder_system_key(out$enzyme_name)
  out <- out[base::order(
    out$enzyme_name,
    -base::as.integer(out$is_gold_standard %in% TRUE),
    -base::nchar(out$rec_seq)
  ), , drop = FALSE]
  out <- out[!base::duplicated(base::paste(out$enzyme_name, out$rec_seq, sep = "\r")), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_rebasefinder_enrich_recognition <- function(hits, rebase_data = NULL, bairoch = NULL) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  if (!"rec_seq" %in% base::names(hits)) hits$rec_seq <- NA_character_
  if (!"substrate_label" %in% base::names(hits)) hits$substrate_label <- NA_character_
  if (!"reference_rec_seq" %in% base::names(hits)) hits$reference_rec_seq <- NA_character_
  if (!"recognition_source" %in% base::names(hits)) hits$recognition_source <- NA_character_
  if (!"recognition_match" %in% base::names(hits)) hits$recognition_match <- NA_character_
  if (!"recognition_match_subject_id" %in% base::names(hits)) hits$recognition_match_subject_id <- NA_character_
  if (!"recognition_donor" %in% base::names(hits)) hits$recognition_donor <- NA_character_

  lookup <- .dnmb_rebasefinder_recognition_lookup(rebase_data, bairoch)
  raw_hit_name <- if ("hit_label" %in% base::names(hits)) {
    base::as.character(hits$hit_label)
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  hit_name <- .dnmb_rebasefinder_clean_rebase_subject(raw_hit_name)
  for (i in base::seq_len(base::nrow(hits))) {
    if (.dnmb_rebasefinder_valid_recognition(hits$rec_seq[[i]])) {
      if (base::is.na(hits$recognition_source[[i]]) || !base::nzchar(hits$recognition_source[[i]])) {
        hits$recognition_source[[i]] <- "pipeline_exact"
      }
      if (base::is.na(hits$recognition_match[[i]]) || !base::nzchar(hits$recognition_match[[i]])) {
        hits$recognition_match[[i]] <- hit_name[[i]]
      }
      if (base::is.na(hits$recognition_match_subject_id[[i]]) ||
          !base::nzchar(hits$recognition_match_subject_id[[i]])) {
        hits$recognition_match_subject_id[[i]] <- raw_hit_name[[i]]
      }
      if (!.dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq[[i]])) {
        hits$reference_rec_seq[[i]] <- hits$rec_seq[[i]]
      }
      if (base::is.na(hits$recognition_donor[[i]]) || !base::nzchar(hits$recognition_donor[[i]])) {
        hits$recognition_donor[[i]] <- hit_name[[i]]
      }
      hits$substrate_label[[i]] <- hits$rec_seq[[i]]
      next
    }
    if (.dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq[[i]])) {
      if (base::is.na(hits$recognition_source[[i]]) || !base::nzchar(hits$recognition_source[[i]])) {
        hits$recognition_source[[i]] <- "pipeline_reference_exact"
      }
      if (base::is.na(hits$recognition_match[[i]]) || !base::nzchar(hits$recognition_match[[i]])) {
        hits$recognition_match[[i]] <- hit_name[[i]]
      }
      if (base::is.na(hits$recognition_match_subject_id[[i]]) ||
          !base::nzchar(hits$recognition_match_subject_id[[i]])) {
        hits$recognition_match_subject_id[[i]] <- raw_hit_name[[i]]
      }
      if (base::is.na(hits$recognition_donor[[i]]) || !base::nzchar(hits$recognition_donor[[i]])) {
        hits$recognition_donor[[i]] <- hit_name[[i]]
      }
      next
    }
    if (!base::nrow(lookup) || base::is.na(hit_name[[i]]) || !base::nzchar(hit_name[[i]])) next
    exact <- base::which(lookup$enzyme_name == hit_name[[i]])
    relation <- "exact"
    candidate <- exact
    if (!base::length(candidate)) {
      key <- .dnmb_rebasefinder_system_key(hit_name[[i]])
      candidate <- base::which(lookup$system_key == key)
      relation <- "cognate_system"
    }
    if (!base::length(candidate)) next
    donor <- candidate[[1]]
    hits$reference_rec_seq[[i]] <- lookup$rec_seq[[donor]]
    hits$recognition_match[[i]] <- hit_name[[i]]
    hits$recognition_match_subject_id[[i]] <- raw_hit_name[[i]]
    hits$recognition_donor[[i]] <- lookup$enzyme_name[[donor]]
    hits$recognition_source[[i]] <- base::paste0(lookup$recognition_source[[donor]], "_", relation)
  }
  hits
}

.dnmb_rebasefinder_accession_key <- function(x) {
  x <- base::toupper(base::trimws(base::as.character(x)))
  x <- base::sub("^NZ_", "", x)
  base::sub("\\.[0-9]+$", "", x)
}

.dnmb_rebasefinder_html_text <- function(x) {
  x <- base::gsub("<[^>]+>", "", x, perl = TRUE)
  x <- base::gsub("&nbsp;", " ", x, fixed = TRUE)
  x <- base::gsub("&amp;", "&", x, fixed = TRUE)
  x <- base::gsub("&#39;|&apos;", "'", x, perl = TRUE)
  x <- base::gsub("&quot;", '"', x, fixed = TRUE)
  base::trimws(base::gsub("[[:space:]]+", " ", x, perl = TRUE))
}

.dnmb_rebasefinder_role_from_hit <- function(hit_name) {
  hit_name <- base::as.character(hit_name)
  base::ifelse(base::grepl("^M[0-9]*\\.", hit_name), "M",
    base::ifelse(base::grepl("^R[0-9]*\\.", hit_name), "R",
      base::ifelse(base::grepl("^S[0-9]*\\.", hit_name), "S", NA_character_)
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
  if (identical(family_id, "Type IV") && identical(enzyme_role, "R")) return(100L)
  if (identical(enzyme_role, "RM")) return(250L)
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
  status <- base::ifelse(
    short | annotated_partial,
    "partial_or_short",
    base::ifelse(base::is.na(aa_len), "unknown_length", "full_length")
  )
  reason <- base::vapply(base::seq_len(n), function(i) {
    parts <- c(
      if (short[[i]]) base::sprintf("aa_len=%s below expected_min=%s", aa_len[[i]], expected_min[[i]]) else NULL,
      if (annotated_partial[[i]]) "annotation suggests partial/truncated sequence" else NULL,
      if (base::is.na(aa_len[[i]])) "translation unavailable" else NULL
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

  explicit_res <- base::grepl(
    "type[ _-]?iii[^|;]{0,60}(restriction (endonuclease|enzyme)|endonuclease|res subunit)|\\bres[ _-]?(iii|3)\\b|\\bres subunit\\b",
    text,
    perl = TRUE
  )
  explicit_mod <- base::grepl(
    "type[ _-]?iii[^|;]{0,60}(mod subunit|modification methylase|methylase|methyltransferase)|\\bmod[ _-]?(iii|3)\\b|\\bmod subunit\\b|modification methylase|pf01555|ipr002295",
    text,
    perl = TRUE
  )

  if (base::grepl("specificity subunit|hsds|target recognition", text, perl = TRUE)) return("S")
  if (explicit_res) return("R")
  if (explicit_mod) return("M")
  if (base::grepl("res subunit|restriction endonuclease subunit r|type[ _-]?iii restriction enzyme", text, perl = TRUE)) return("R")
  if (base::grepl("mod subunit|type[ _-]?iii methyltransferase", text, perl = TRUE)) return("M")
  NA_character_
}

.dnmb_rebasefinder_typeiii_sequence_signals <- function(translation) {
  seq <- .dnmb_rebasefinder_normalize_protein(translation)
  if (base::is.na(seq)) {
    return(list(
      role = NA_character_,
      signals = character(),
      has_sam = FALSE,
      has_mtase = FALSE,
      res_signal_count = 0L
    ))
  }
  mtase <- .dnmb_rebasefinder_mtase_sequence_signals(seq)
  motif_positions <- function(name) {
    .dnmb_rebasefinder_motif_all_positions(seq, .dnmb_rebasefinder_motif_pattern(name))
  }
  motor <- .dnmb_rebasefinder_motif_ordered(
    list(motif_positions("ResIII-WA"), motif_positions("ResIII-WB"), motif_positions("ResIII-MIII")),
    max_span = 750L
  )
  hits <- c(
    SAM = !base::is.na(mtase$sam_match),
    MTase = mtase$verified,
    `ResIII-WA` = base::length(motif_positions("ResIII-WA")) > 0L,
    `ResIII-WB` = base::length(motif_positions("ResIII-WB")) > 0L,
    `ResIII-MIII` = base::length(motif_positions("ResIII-MIII")) > 0L,
    `ResIII-PD` = base::length(motif_positions("ResIII-PD")) > 0L
  )
  res_signal_count <- base::sum(hits[c("ResIII-WA", "ResIII-WB", "ResIII-MIII", "ResIII-PD")])
  role <- if (base::nchar(seq) >= 600L && motor && hits[["ResIII-PD"]]) {
    "R"
  } else if (isTRUE(mtase$verified)) {
    "M"
  } else {
    NA_character_
  }
  list(
    role = role,
    signals = base::names(hits)[hits],
    has_sam = hits[["SAM"]],
    has_mtase = hits[["MTase"]],
    has_helicase = motor,
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
  if (base::grepl("hsdm|class[ _-]?i .*methyl|type[ _-]?i .*methyl|n-6 dna methylase|sam-dependent dna methyltransferase|dna methyltransferase", text, perl = TRUE)) {
    return("M")
  }
  NA_character_
}

.dnmb_rebasefinder_typei_sequence_signals <- function(translation) {
  seq <- .dnmb_rebasefinder_normalize_protein(translation)
  if (base::is.na(seq)) {
    return(list(role = NA_character_, signals = character()))
  }
  mtase <- .dnmb_rebasefinder_mtase_sequence_signals(seq)
  motif_positions <- function(name) {
    .dnmb_rebasefinder_motif_all_positions(seq, .dnmb_rebasefinder_motif_pattern(name))
  }
  motor <- .dnmb_rebasefinder_motif_ordered(
    list(motif_positions("P-loop"), motif_positions("HsdR-WB"), motif_positions("HsdR-MIII")),
    max_span = 750L
  )
  hits <- c(
    SAM = !base::is.na(mtase$sam_match),
    MTase = mtase$verified,
    Ploop = base::length(motif_positions("P-loop")) > 0L,
    `HsdR-WB` = base::length(motif_positions("HsdR-WB")) > 0L,
    `HsdR-MIII` = base::length(motif_positions("HsdR-MIII")) > 0L,
    PD = base::length(motif_positions("HsdR-PD")) > 0L
  )
  role <- if (base::nchar(seq) >= 600L && motor && hits[["PD"]]) {
    "R"
  } else if (isTRUE(mtase$verified)) {
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

.dnmb_rebasefinder_strip_typeii_support <- function(x) {
  base::vapply(base::as.character(x), function(one) {
    if (base::is.na(one) || !base::nzchar(one)) return(NA_character_)
    parts <- base::trimws(base::strsplit(one, ";", fixed = TRUE)[[1]])
    parts <- parts[base::nzchar(parts)]
    parts <- parts[!base::grepl("^(Type II operon_context=|anchor=|component=|motifs=)", parts)]
    if (!base::length(parts)) NA_character_ else base::paste(parts, collapse = "; ")
  }, character(1))
}

.dnmb_rebasefinder_reset_operon_context <- function(hits, family) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  family <- base::match.arg(family, c("typei", "typeii", "typeiii"))
  specs <- switch(
    family,
    typei = list(
      label = "typeI", supported = "typei_operon_supported",
      columns = list(
        typei_context_status = NA_character_, typei_context_roles = NA_character_,
        typei_context_partners = NA_character_, typei_context_sources = NA_character_,
        typei_context_anchor = NA_character_, typei_operon_supported = NA,
        typei_context_window_bp = NA_real_
      ),
      strip = .dnmb_rebasefinder_strip_typei_support
    ),
    typeii = list(
      label = "typeII", supported = "typeii_operon_supported",
      columns = list(
        typeii_context_status = NA_character_, typeii_context_roles = NA_character_,
        typeii_context_partners = NA_character_, typeii_context_sources = NA_character_,
        typeii_operon_candidate = NA, typeii_operon_supported = NA,
        typeii_context_window_bp = NA_real_, operon_component = NA_character_
      ),
      strip = .dnmb_rebasefinder_strip_typeii_support
    ),
    typeiii = list(
      label = "typeIII", supported = "typeiii_operon_supported",
      columns = list(
        typeiii_context_status = NA_character_, typeiii_context_roles = NA_character_,
        typeiii_context_partners = NA_character_, typeiii_context_sources = NA_character_,
        typeiii_context_anchor = NA_character_, typeiii_operon_supported = NA,
        typeiii_context_window_bp = NA_real_
      ),
      strip = .dnmb_rebasefinder_strip_typeiii_support
    )
  )
  if (!base::nrow(hits)) {
    for (col in base::names(specs$columns)) hits[[col]] <- base::rep(specs$columns[[col]], 0L)
    return(hits)
  }

  previous_supported <- if (specs$supported %in% base::names(hits)) {
    hits[[specs$supported]] %in% TRUE
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  if ("typing_eligible" %in% base::names(hits) && base::any(previous_supported)) {
    baseline <- base::rep(FALSE, base::nrow(hits))
    if ("raw_typing_eligible" %in% base::names(hits)) baseline <- baseline | hits$raw_typing_eligible %in% TRUE
    if ("curated_blast_promoted" %in% base::names(hits)) baseline <- baseline | hits$curated_blast_promoted %in% TRUE
    if ("structure_supported" %in% base::names(hits)) baseline <- baseline | hits$structure_supported %in% TRUE
    hits$typing_eligible[previous_supported] <- baseline[previous_supported]
  }

  hit_label <- if ("hit_label" %in% base::names(hits)) base::as.character(hits$hit_label) else NA_character_
  evidence_mode <- if ("evidence_mode" %in% base::names(hits)) base::as.character(hits$evidence_mode) else NA_character_
  promoted <- if ("curated_blast_promoted" %in% base::names(hits)) hits$curated_blast_promoted %in% TRUE else FALSE
  context_only <- !base::is.na(hit_label) &
    base::grepl(base::paste0("^", specs$label, "_context:"), hit_label) &
    !base::is.na(evidence_mode) & base::grepl("operon_context", evidence_mode, fixed = TRUE) &
    !promoted
  hits <- hits[!context_only, , drop = FALSE]

  for (col in base::names(specs$columns)) {
    hits[[col]] <- base::rep(specs$columns[[col]], base::nrow(hits))
  }
  if ("support" %in% base::names(hits)) hits$support <- specs$strip(hits$support)
  if ("operon_id" %in% base::names(hits)) {
    stale_operon <- !base::is.na(hits$operon_id) &
      base::grepl(
        base::paste0("^DNMB_", specs$label, "_"),
        base::as.character(hits$operon_id),
        ignore.case = TRUE
      )
    replacement <- if ("raw_operon_id" %in% base::names(hits)) {
      base::as.character(hits$raw_operon_id[stale_operon])
    } else {
      base::rep(NA_character_, base::sum(stale_operon))
    }
    hits$operon_id[stale_operon] <- replacement
  }
  base::rownames(hits) <- NULL
  hits
}

.dnmb_rebasefinder_gene_order <- function(genes) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  for (col in c("start", "end")) {
    if (col %in% base::names(genes)) genes[[col]] <- suppressWarnings(base::as.numeric(genes[[col]]))
  }
  if (!"contig" %in% base::names(genes)) genes$contig <- "contig"
  if ("contig_number" %in% base::names(genes)) {
    contig_number <- base::trimws(base::as.character(genes$contig_number))
    has_number <- !base::is.na(contig_number) & base::nzchar(contig_number)
    genes$.dnmb_contig_label <- genes$contig
    genes$contig[has_number] <- base::paste0(
      base::as.character(genes$contig[has_number]),
      "::record_",
      contig_number[has_number]
    )
  }
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

.dnmb_rebasefinder_coherent_gene_path <- function(genes_ord, i, j) {
  if (base::is.na(i) || base::is.na(j) || i < 1L || j < 1L ||
      i > base::nrow(genes_ord) || j > base::nrow(genes_ord)) return(FALSE)
  idx <- base::seq.int(base::min(i, j), base::max(i, j))
  contigs <- base::as.character(genes_ord$contig[idx])
  if (base::length(base::unique(contigs[!base::is.na(contigs)])) > 1L) return(FALSE)
  directions <- base::as.character(genes_ord$direction[idx])
  directions <- directions[!base::is.na(directions) & base::nzchar(directions)]
  base::length(base::unique(directions)) <= 1L
}

.dnmb_rebasefinder_single_supported_role <- function(x) {
  x <- base::as.character(x)[1]
  if (base::is.na(x) || !base::nzchar(x)) return(NA_character_)
  roles <- base::unique(base::trimws(base::strsplit(x, "/", fixed = TRUE)[[1]]))
  roles <- roles[roles %in% c("M", "R", "S", "RM")]
  if (base::length(roles) == 1L) roles[[1]] else NA_character_
}

.dnmb_rebasefinder_role_supported <- function(role, supported_roles) {
  role <- base::as.character(role)[1]
  supported_roles <- base::as.character(supported_roles)[1]
  if (base::is.na(role) || !base::nzchar(role) ||
      base::is.na(supported_roles) || !base::nzchar(supported_roles)) return(FALSE)
  roles <- base::unique(base::trimws(base::strsplit(supported_roles, "/", fixed = TRUE)[[1]]))
  role %in% roles || (role %in% c("M", "R") && "RM" %in% roles)
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
                                                         min_alignment_coverage = 0.40,
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
  qtm_col <- pick(c("qtmscore", "qtm"))
  ttm_col <- pick(c("ttmscore", "ttm"))
  rmsd_col <- pick(c("rmsd"))
  lddt_col <- pick(c("lddt"))
  qcov_col <- pick(c("qcov", "query_coverage"))
  tcov_col <- pick(c("tcov", "target_coverage"))
  qstart_col <- pick(c("qstart"))
  qend_col <- pick(c("qend"))
  qlen_col <- pick(c("qlen"))
  tstart_col <- pick(c("tstart"))
  tend_col <- pick(c("tend"))
  tlen_col <- pick(c("tlen"))
  alnlen_col <- pick(c("alnlen", "alignment_length"))
  qaln_col <- pick(c("qaln", "qseq"))
  taln_col <- pick(c("taln", "tseq"))

  numeric_col <- function(column) {
    if (!base::is.na(column)) suppressWarnings(base::as.numeric(tbl[[column]])) else NA_real_
  }
  coverage_col <- function(column, start_col, end_col, length_col) {
    coverage <- numeric_col(column)
    if (base::length(coverage) == 1L && base::is.na(coverage)) {
      start <- numeric_col(start_col)
      end <- numeric_col(end_col)
      length <- numeric_col(length_col)
      coverage <- (base::abs(end - start) + 1) / length
    }
    coverage[base::is.finite(coverage) & coverage > 1] <- coverage[base::is.finite(coverage) & coverage > 1] / 100
    coverage
  }
  query_coverage <- coverage_col(qcov_col, qstart_col, qend_col, qlen_col)
  target_coverage <- coverage_col(tcov_col, tstart_col, tend_col, tlen_col)
  alignment_coverage <- base::ifelse(
    base::is.finite(query_coverage) & base::is.finite(target_coverage),
    base::pmin(query_coverage, target_coverage),
    base::ifelse(base::is.finite(query_coverage), query_coverage, target_coverage)
  )

  out <- data.frame(
    query = .dnmb_module_clean_annotation_key(tbl[[q_col]]),
    structure_hit = base::as.character(tbl[[t_col]]),
    structure_evalue = if (!base::is.na(e_col)) suppressWarnings(base::as.numeric(tbl[[e_col]])) else NA_real_,
    structure_bitscore = if (!base::is.na(b_col)) suppressWarnings(base::as.numeric(tbl[[b_col]])) else NA_real_,
    structure_probability = if (!base::is.na(p_col)) suppressWarnings(base::as.numeric(tbl[[p_col]])) else NA_real_,
    structure_tmscore = if (!base::is.na(tm_col)) suppressWarnings(base::as.numeric(tbl[[tm_col]])) else NA_real_,
    structure_query_tmscore = numeric_col(qtm_col),
    structure_target_tmscore = numeric_col(ttm_col),
    structure_rmsd = if (!base::is.na(rmsd_col)) suppressWarnings(base::as.numeric(tbl[[rmsd_col]])) else NA_real_,
    structure_lddt = numeric_col(lddt_col),
    structure_query_coverage = query_coverage,
    structure_target_coverage = target_coverage,
    structure_alignment_coverage = alignment_coverage,
    structure_query_start = numeric_col(qstart_col),
    structure_query_end = numeric_col(qend_col),
    structure_target_start = numeric_col(tstart_col),
    structure_target_end = numeric_col(tend_col),
    structure_alignment_length = numeric_col(alnlen_col),
    structure_query_alignment = if (!base::is.na(qaln_col)) base::as.character(tbl[[qaln_col]]) else NA_character_,
    structure_target_alignment = if (!base::is.na(taln_col)) base::as.character(tbl[[taln_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out <- out[!base::is.na(out$query) & base::nzchar(out$query), , drop = FALSE]
  if (!base::nrow(out)) return(NULL)
  out$structure_reference_id <- .dnmb_rebasefinder_structure_reference_id(out$structure_hit)
  out$structure_chain <- .dnmb_rebasefinder_structure_chain_id(out$structure_hit)
  refs <- .dnmb_rebasefinder_structure_reference_manifest(reference_manifest)
  ref_idx <- base::rep(NA_integer_, base::nrow(out))
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
  out$structure_family_raw <- out$structure_family
  out$structure_role_raw <- out$structure_role
  out$structure_reference_known <- !base::is.na(ref_idx)
  out$structure_family <- .dnmb_rebasefinder_canonical_structure_family(out$structure_family)
  out$structure_role <- .dnmb_rebasefinder_canonical_structure_role(out$structure_role)
  out$structure_chain_role <- .dnmb_rebasefinder_canonical_structure_role(out$structure_chain_role)
  excluded_chain <- base::tolower(out$structure_role) %in% c("exclude", "excluded", "non_rm", "ocr", "anti_restriction")
  out$structure_classified <- out$structure_reference_known |
    (!base::is.na(out$structure_family) & base::nzchar(out$structure_family) &
       !base::is.na(out$structure_role) & base::nzchar(out$structure_role))
  coverage_available <- base::is.finite(out$structure_query_coverage)
  coverage_pass <- !coverage_available | out$structure_query_coverage >= min_alignment_coverage
  out$structure_pass <- out$structure_classified & !excluded_chain & coverage_pass & (
    (!base::is.na(out$structure_evalue) & out$structure_evalue <= max_evalue) |
      (!base::is.na(out$structure_probability) & out$structure_probability >= min_probability) |
      (!base::is.na(out$structure_tmscore) & out$structure_tmscore >= min_tmscore)
  )
  out$structure_status <- base::ifelse(excluded_chain, "structure_excluded",
    base::ifelse(!coverage_pass, "structure_low_coverage",
    base::ifelse(!out$structure_classified, "structure_unclassified_review",
      base::ifelse(out$structure_pass, "structure_supported", "structure_weak"))))
  out <- out[base::order(
    out$query,
    !out$structure_pass,
    base::ifelse(base::is.na(out$structure_evalue), Inf, out$structure_evalue),
    -base::ifelse(base::is.na(out$structure_probability), -Inf, out$structure_probability),
    -base::ifelse(base::is.na(out$structure_tmscore), -Inf, out$structure_tmscore),
    -base::ifelse(base::is.na(out$structure_query_coverage), -Inf, out$structure_query_coverage),
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
    if ("structure_query_coverage" %in% base::names(row) &&
        !base::is.na(row$structure_query_coverage)) {
      base::paste0("qcov=", base::round(row$structure_query_coverage, 3))
    } else NULL,
    if (!base::is.na(row$structure_evalue)) base::paste0("fold_e=", format(row$structure_evalue, scientific = TRUE, digits = 3)) else NULL
  )
  if (!base::length(parts)) NA_character_ else base::paste(parts, collapse = "; ")
}

.dnmb_rebasefinder_canonical_structure_family <- function(x) {
  x <- base::trimws(base::as.character(x))
  out <- x
  out[base::grepl("^Type[ _-]?ISP\\b", x, ignore.case = TRUE)] <- "Type I"
	out[base::grepl("^Type[ _-]?II(?:S|G|L|C|P)?\\b", x, ignore.case = TRUE)] <- "Type II"
  out[base::grepl("^Type[ _-]?III\\b", x, ignore.case = TRUE)] <- "Type III"
  out[base::grepl("^Type[ _-]?IV\\b", x, ignore.case = TRUE)] <- "Type IV"
  out[base::grepl("^Type[ _-]?I\\b", x, ignore.case = TRUE)] <- "Type I"
  out
}

.dnmb_rebasefinder_canonical_structure_role <- function(x) {
  x <- base::toupper(base::trimws(base::as.character(x)))
  out <- x
  tokens <- base::strsplit(base::ifelse(base::is.na(x), "", x), "[/,+;[:space:]]+", perl = TRUE)
  out <- base::vapply(base::seq_along(tokens), function(i) {
    if (x[[i]] %in% c("EXCLUDE", "EXCLUDED", "NON_RM", "OCR", "ANTI_RESTRICTION")) {
      return(base::tolower(x[[i]]))
    }
    one <- tokens[[i]]
    one <- base::unique(one[one %in% c("M", "R", "S", "RM")])
    if (!base::length(one)) return(NA_character_)
    if ("RM" %in% one || base::all(c("R", "M") %in% one)) return("RM")
    one[[1]]
  }, character(1))
  out
}

.dnmb_rebasefinder_structure_family_compatible <- function(candidate_family, structure_family) {
  candidate_family <- .dnmb_rebasefinder_canonical_structure_family(candidate_family)
  structure_family <- .dnmb_rebasefinder_canonical_structure_family(structure_family)
  base::vapply(base::seq_along(candidate_family), function(i) {
    cf <- candidate_family[[i]]
    sf <- structure_family[[i]]
    if (base::is.na(cf) || !base::nzchar(cf) || base::is.na(sf) || !base::nzchar(sf)) return(TRUE)
    if (identical(cf, sf)) return(TRUE)
    FALSE
  }, logical(1))
}

.dnmb_rebasefinder_structure_role_compatible <- function(candidate_role, structure_role, structure_chain_role = NA_character_) {
  candidate_role <- .dnmb_rebasefinder_canonical_structure_role(candidate_role)
  structure_role <- .dnmb_rebasefinder_canonical_structure_role(structure_role)
  structure_chain_role <- .dnmb_rebasefinder_canonical_structure_role(structure_chain_role)
  base::vapply(base::seq_along(candidate_role), function(i) {
    cr <- candidate_role[[i]]
    if (base::is.na(cr) || !base::nzchar(cr)) return(TRUE)
	    roles <- stats::na.omit(c(structure_chain_role[[i]], structure_role[[i]]))
	    roles <- roles[base::nzchar(roles)]
	    if (!base::length(roles)) return(TRUE)
	    cr %in% roles
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
    hits$structure_supported <- pass
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
    add$structure_supported <- add$structure_pass %in% TRUE & add$structure_candidate_consistent %in% TRUE
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
  hits <- .dnmb_rebasefinder_reset_operon_context(hits, "typeiii")
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
  annotation_role <- base::mapply(
    .dnmb_rebasefinder_typeiii_context_role,
    text = anno_text,
    known_rm_type = NA_character_,
    known_role = NA_character_,
    USE.NAMES = FALSE
  )
  typeiii_specific_annotation <- !base::is.na(anno_text) & base::grepl(
    "type[ _-]?iii|\\bres[ _-]?(iii|3)\\b|\\bmod[ _-]?(iii|3)\\b|\\bres subunit\\b|\\bmod subunit\\b|pf04851|ipr006935|pf01555|ipr002295",
    anno_text,
    perl = TRUE
  )
  genes_ord$.hit_typing_eligible <- FALSE
  if ("typing_eligible" %in% base::names(hits)) {
    genes_ord$.hit_typing_eligible[have_hit] <- hits$typing_eligible[hit_idx[have_hit]] %in% TRUE
  }
  genes_ord$.hit_strict_blast <- FALSE
  strict_blast <- .dnmb_rebasefinder_strict_blast_support(hits)
  genes_ord$.hit_strict_blast[have_hit] <- strict_blast[hit_idx[have_hit]]
  genes_ord$.hit_structure_supported <- FALSE
  structure_supported <- .dnmb_rebasefinder_consistent_structure_support(hits)
  genes_ord$.hit_structure_supported[have_hit] <- structure_supported[hit_idx[have_hit]]
  trusted_typeiii_hit <- !base::is.na(genes_ord$.rm_type) & genes_ord$.rm_type == "Type III" &
    !base::is.na(genes_ord$.rm_role) & genes_ord$.rm_role %in% c("M", "R", "S") &
    (genes_ord$.hit_strict_blast | genes_ord$.hit_structure_supported |
       (typeiii_specific_annotation & !base::is.na(annotation_role) & annotation_role == genes_ord$.rm_role))
  context_role <- annotation_role
  context_role[trusted_typeiii_hit] <- genes_ord$.rm_role[trusted_typeiii_hit]
  genes_ord$.typeiii_blast_roles <- NA_character_
  if ("blast_supported_typeiii_roles" %in% base::names(hits)) {
    genes_ord$.typeiii_blast_roles[have_hit] <- hits$blast_supported_typeiii_roles[hit_idx[have_hit]]
  }
  alternate_role <- base::vapply(
    genes_ord$.typeiii_blast_roles,
    .dnmb_rebasefinder_single_supported_role,
    character(1)
  )
  use_alternate_role <- base::is.na(context_role) & !base::is.na(alternate_role)
  context_role[use_alternate_role] <- alternate_role[use_alternate_role]
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
  annotation_flags <- .dnmb_rebasefinder_annotation_flags(anno_text)
  canonical_existing <- !base::is.na(genes_ord$.rm_type) &
    genes_ord$.rm_type %in% c("Type I", "Type II", "Type III", "Type IV")
  family_conflict <- canonical_existing & genes_ord$.rm_type != "Type III"
  role_conflict <- canonical_existing & !base::is.na(genes_ord$.rm_role) &
    genes_ord$.rm_role %in% c("M", "R", "S", "RM") &
    !base::is.na(genes_ord$.typeiii_role) & genes_ord$.rm_role != genes_ord$.typeiii_role
  alternate_role_support <- base::mapply(
    .dnmb_rebasefinder_role_supported,
    genes_ord$.typeiii_role,
    genes_ord$.typeiii_blast_roles,
    USE.NAMES = FALSE
  )
  weak_cross_family_mod <- family_conflict &
    !genes_ord$.hit_strict_blast & !genes_ord$.hit_structure_supported &
    genes_ord$.typeiii_role %in% "M" & annotation_flags$dna_mtase &
    !annotation_flags$explicit_non_dna_mtase & !annotation_flags$other_defense
  incompatible_existing <- (family_conflict | role_conflict) &
    !alternate_role_support & !weak_cross_family_mod
  incompatible_existing[base::is.na(incompatible_existing)] <- FALSE
  genes_ord$.typeiii_role[incompatible_existing] <- NA_character_
  role_from_alternate <- alternate_role_support &
    (family_conflict | role_conflict | use_alternate_role) & !incompatible_existing
  invalid_typeiii_role <-
    (genes_ord$.typeiii_role == "M" &
       (annotation_flags$explicit_non_dna_mtase | annotation_flags$other_defense | annotation_flags$unrelated_metabolic)) |
    (genes_ord$.typeiii_role == "R" &
       (annotation_flags$other_defense | annotation_flags$repair_nuclease |
          annotation_flags$unrelated_helicase | annotation_flags$toxin_atpase |
          annotation_flags$unrelated_metabolic | annotation_flags$mobile_nuclease) &
       !annotation_flags$rm_rease)
  invalid_typeiii_role[base::is.na(invalid_typeiii_role)] <- FALSE
  genes_ord$.typeiii_role[invalid_typeiii_role] <- NA_character_
  evidence_mode <- if ("evidence_mode" %in% base::names(hits)) base::as.character(hits$evidence_mode) else base::rep(NA_character_, base::nrow(hits))
  genes_ord$.evidence_mode <- NA_character_
  genes_ord$.evidence_mode[have_hit] <- evidence_mode[hit_idx[have_hit]]
  operon_context_hit <- !base::is.na(genes_ord$.evidence_mode) & genes_ord$.evidence_mode == "operon_context"
  known_typeiii_role <- !operon_context_hit & trusted_typeiii_hit
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
    base::ifelse(role_from_alternate, "rebase_alternate",
    base::ifelse(!base::is.na(operon_source), operon_source,
                 base::ifelse(!base::is.na(context_role), "annotation",
                 base::ifelse(use_seq_role, "sequence_motif", NA_character_))
    ))
  )
  genes_ord$.typeiii_sequence_signals <- seq_signal_label

  is_typeiii_seed <- trusted_typeiii_hit & !base::is.na(genes_ord$.typeiii_role)
  if (!base::any(is_typeiii_seed)) {
    for (col in c("typeiii_context_status", "typeiii_context_roles", "typeiii_context_partners",
                  "typeiii_context_sources", "typeiii_context_anchor", "typeiii_operon_supported", "typeiii_context_window_bp")) {
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
    ctx$coherent_path <- base::vapply(ctx_idx, function(j) {
      .dnmb_rebasefinder_coherent_gene_path(genes_ord, seed_i, j)
    }, logical(1))
    ctx$operon_link <- ctx$same_direction & ctx$coherent_path &
      (base::is.na(ctx$distance_bp) | ctx$distance_bp <= max_operon_gap) &
      base::abs(ctx$relative_gene) <= max_intervening + 1L

    member_mask <- !base::is.na(ctx$.typeiii_role) & (ctx$operon_link | ctx$relative_gene == 0)
    source_rows <- ctx[member_mask, , drop = FALSE]
    roles_present <- base::sort(base::unique(stats::na.omit(source_rows$.typeiii_role)))
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
    if (base::nrow(source_rows)) {
      for (member_i in base::seq_len(base::nrow(source_rows))) {
        member <- source_rows[member_i, , drop = FALSE]
        partners <- source_rows[source_rows$locus_tag != member$locus_tag, , drop = FALSE]
        partner_label <- if (base::nrow(partners)) {
          partner_gap <- base::vapply(base::seq_len(base::nrow(partners)), function(j) {
            .dnmb_rebasefinder_intergenic_gap(member, partners[j, , drop = FALSE])
          }, numeric(1))
          base::paste(
            base::sprintf(
              "%s:%s:%s:%+dg:%sbp",
              partners$locus_tag,
              partners$.typeiii_role,
              partners$.typeiii_role_source,
              partners$.dnmb_order - member$.dnmb_order,
              base::ifelse(base::is.na(partner_gap), "NA", partner_gap)
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        }
        context_rows[[base::length(context_rows) + 1L]] <- data.frame(
          query = member$locus_tag,
          typeiii_context_status = status,
          typeiii_context_roles = base::paste(roles_present, collapse = "/"),
          typeiii_context_partners = partner_label,
          typeiii_context_sources = source_label,
          typeiii_context_anchor = seed$locus_tag,
          typeiii_operon_supported = has_m && has_r,
          typeiii_context_window_bp = base::max(ctx$end, na.rm = TRUE) - base::min(ctx$start, na.rm = TRUE) + 1,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  context_tbl <- do.call(rbind, context_rows)
  context_tbl <- context_tbl[base::order(
    context_tbl$query,
    !context_tbl$typeiii_operon_supported,
    context_tbl$typeiii_context_window_bp
  ), , drop = FALSE]
  context_tbl <- context_tbl[!base::duplicated(context_tbl$query), , drop = FALSE]

  for (col in base::setdiff(base::names(context_tbl), "query")) {
    if (!col %in% base::names(hits)) hits[[col]] <- .dnmb_na_vector_like(context_tbl[[col]], base::nrow(hits))
  }
  hidx <- base::match(hits$query, context_tbl$query)
  include <- !base::is.na(hidx)
  for (col in base::setdiff(base::names(context_tbl), "query")) hits[[col]][include] <- context_tbl[[col]][hidx[include]]
  if (!"operon_id" %in% base::names(hits)) hits$operon_id <- NA_character_
  supported_operon <- include & hits$typeiii_operon_supported %in% TRUE
  if (!"typeiii_context_previous_family" %in% base::names(hits)) {
    hits$typeiii_context_previous_family <- NA_character_
  }
  if (!"typeiii_context_reclassified" %in% base::names(hits)) {
    hits$typeiii_context_reclassified <- FALSE
  }
  strict_hit <- .dnmb_rebasefinder_strict_blast_support(hits)
  reclassify <- supported_operon & !strict_hit &
    !base::is.na(hits$family_id) & hits$family_id != "Type III" &
    hits$enzyme_role %in% c("M", "R")
  hits$typeiii_context_previous_family[reclassify] <- hits$family_id[reclassify]
  hits$typeiii_context_reclassified[reclassify] <- TRUE
  hits$family_id[reclassify] <- "Type III"
  hits$operon_id[supported_operon] <- base::paste0(
    "DNMB_TypeIII_",
    context_tbl$typeiii_context_anchor[hidx[supported_operon]]
  )
  missing_operon <- include & (base::is.na(hits$operon_id) | !base::nzchar(base::as.character(hits$operon_id)))
  hits$operon_id[missing_operon] <- base::paste0(
    "DNMB_TypeIII_",
    context_tbl$typeiii_context_anchor[hidx[missing_operon]]
  )

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
  candidate_idx <- base::which(
    genes_ord$locus_tag %in% context_tbl$query &
      !genes_ord$locus_tag %in% hits$query &
      !base::is.na(genes_ord$.typeiii_role) &
      genes_ord$.typeiii_role %in% c("M", "R")
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
    ctx_idx <- base::match(add$query, context_tbl$query)
    add$typeiii_context_status <- context_tbl$typeiii_context_status[ctx_idx]
    add$typeiii_context_roles <- context_tbl$typeiii_context_roles[ctx_idx]
    add$typeiii_context_partners <- context_tbl$typeiii_context_partners[ctx_idx]
    add$typeiii_context_sources <- context_tbl$typeiii_context_sources[ctx_idx]
    add$typeiii_context_anchor <- context_tbl$typeiii_context_anchor[ctx_idx]
    add$typeiii_operon_supported <- context_tbl$typeiii_operon_supported[ctx_idx]
    add$typeiii_context_window_bp <- context_tbl$typeiii_context_window_bp[ctx_idx]
    add$operon_id <- base::paste0("DNMB_TypeIII_", add$typeiii_context_anchor)
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
  hits <- .dnmb_rebasefinder_reset_operon_context(hits, "typei")
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
  annotation_role <- base::mapply(
    .dnmb_rebasefinder_typei_context_role,
    text = anno_text,
    known_rm_type = NA_character_,
    known_role = NA_character_,
    USE.NAMES = FALSE
  )
  typei_specific_annotation <- !base::is.na(anno_text) & base::grepl(
    "\\bhsd[mrs]\\b|type[ _-]?i restriction|class[ _-]?i[^|;]{0,35}(methyl|specificity)|tigr00348|tigr00497|rmtype1|ecor124i",
    anno_text,
    perl = TRUE
  )
  genes_ord$.hit_typing_eligible <- FALSE
  if ("typing_eligible" %in% base::names(hits)) {
    genes_ord$.hit_typing_eligible[have_hit] <- hits$typing_eligible[hit_idx[have_hit]] %in% TRUE
  }
  genes_ord$.hit_strict_blast <- FALSE
  strict_blast <- .dnmb_rebasefinder_strict_blast_support(hits)
  genes_ord$.hit_strict_blast[have_hit] <- strict_blast[hit_idx[have_hit]]
  genes_ord$.hit_structure_supported <- FALSE
  structure_supported <- .dnmb_rebasefinder_consistent_structure_support(hits)
  genes_ord$.hit_structure_supported[have_hit] <- structure_supported[hit_idx[have_hit]]
  trusted_typei_hit <- !base::is.na(genes_ord$.rm_type) & genes_ord$.rm_type == "Type I" &
    !base::is.na(genes_ord$.rm_role) & genes_ord$.rm_role %in% c("M", "R", "S") &
    (genes_ord$.hit_strict_blast | genes_ord$.hit_structure_supported |
       (typei_specific_annotation & !base::is.na(annotation_role) & annotation_role == genes_ord$.rm_role))
  context_role <- annotation_role
  context_role[trusted_typei_hit] <- genes_ord$.rm_role[trusted_typei_hit]
  genes_ord$.typei_blast_roles <- NA_character_
  if ("blast_supported_typei_roles" %in% base::names(hits)) {
    genes_ord$.typei_blast_roles[have_hit] <- hits$blast_supported_typei_roles[hit_idx[have_hit]]
  }
  alternate_role <- base::vapply(
    genes_ord$.typei_blast_roles,
    .dnmb_rebasefinder_single_supported_role,
    character(1)
  )
  use_alternate_role <- base::is.na(context_role) & !base::is.na(alternate_role)
  context_role[use_alternate_role] <- alternate_role[use_alternate_role]
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
  canonical_existing <- !base::is.na(genes_ord$.rm_type) &
    genes_ord$.rm_type %in% c("Type I", "Type II", "Type III", "Type IV")
  family_conflict <- canonical_existing & genes_ord$.rm_type != "Type I"
  role_conflict <- canonical_existing & !base::is.na(genes_ord$.rm_role) &
    genes_ord$.rm_role %in% c("M", "R", "S", "RM") &
    !base::is.na(genes_ord$.typei_role) & genes_ord$.rm_role != genes_ord$.typei_role
  alternate_role_support <- base::mapply(
    .dnmb_rebasefinder_role_supported,
    genes_ord$.typei_role,
    genes_ord$.typei_blast_roles,
    USE.NAMES = FALSE
  )
  incompatible_existing <- (family_conflict | role_conflict) & !alternate_role_support
  incompatible_existing[base::is.na(incompatible_existing)] <- FALSE
  genes_ord$.typei_role[incompatible_existing] <- NA_character_
  role_from_alternate <- alternate_role_support &
    (family_conflict | role_conflict | use_alternate_role) & !incompatible_existing
  annotation_flags <- .dnmb_rebasefinder_annotation_flags(anno_text)
  invalid_typei_role <-
    (genes_ord$.typei_role == "M" &
       (annotation_flags$explicit_non_dna_mtase | annotation_flags$other_defense | annotation_flags$unrelated_metabolic)) |
    (genes_ord$.typei_role == "R" &
       (annotation_flags$other_defense | annotation_flags$repair_nuclease |
          annotation_flags$unrelated_helicase | annotation_flags$toxin_atpase |
          annotation_flags$unrelated_metabolic | annotation_flags$mobile_nuclease) &
       !annotation_flags$rm_rease)
  invalid_typei_role[base::is.na(invalid_typei_role)] <- FALSE
  genes_ord$.typei_role[invalid_typei_role] <- NA_character_
  evidence_mode <- if ("evidence_mode" %in% base::names(hits)) base::as.character(hits$evidence_mode) else base::rep(NA_character_, base::nrow(hits))
  genes_ord$.evidence_mode <- NA_character_
  genes_ord$.evidence_mode[have_hit] <- evidence_mode[hit_idx[have_hit]]
  operon_context_hit <- !base::is.na(genes_ord$.evidence_mode) & genes_ord$.evidence_mode == "operon_context"
  known_typei_role <- !operon_context_hit & trusted_typei_hit
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
    base::ifelse(role_from_alternate, "rebase_alternate",
    base::ifelse(!base::is.na(operon_source), operon_source,
                 base::ifelse(!base::is.na(context_role), "annotation",
                              base::ifelse(use_seq_role, "sequence_motif", NA_character_))
    ))
  )
  genes_ord$.typei_sequence_signals <- seq_signal_label
  genes_ord$.trusted_typei_hit <- trusted_typei_hit
  genes_ord$.typei_specific_annotation <- typei_specific_annotation
  genes_ord$.extended_typei_r_candidate <- genes_ord$.typei_role == "R" &
    ((typei_specific_annotation & annotation_role == "R") |
       (!base::is.na(seq_role) & seq_role == "R"))
  genes_ord$.extended_typei_r_candidate[base::is.na(genes_ord$.extended_typei_r_candidate)] <- FALSE

  is_typei_seed <- trusted_typei_hit & !base::is.na(genes_ord$.typei_role)
  if (!base::any(is_typei_seed)) {
    for (col in c("typei_context_status", "typei_context_roles", "typei_context_partners",
                  "typei_context_sources", "typei_context_anchor", "typei_operon_supported", "typei_context_window_bp")) {
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
    ctx$coherent_path <- base::vapply(ctx_idx, function(j) {
      .dnmb_rebasefinder_coherent_gene_path(genes_ord, seed_i, j)
    }, logical(1))
    standard_neighbor_link <- base::abs(ctx$relative_gene) <=
      base::min(max_neighbors, max_intervening + 1L)
    extended_trusted_typei_link <- base::abs(ctx$relative_gene) <= max_neighbors &
      seed$.trusted_typei_hit %in% TRUE &
      ((ctx$.trusted_typei_hit %in% TRUE &
          (seed$.typei_specific_annotation %in% TRUE | ctx$.typei_specific_annotation %in% TRUE)) |
         ctx$.extended_typei_r_candidate %in% TRUE)
    ctx$operon_link <- ctx$same_direction & ctx$coherent_path &
      (base::is.na(ctx$distance_bp) | ctx$distance_bp <= max_operon_gap) &
      (standard_neighbor_link | extended_trusted_typei_link)

    member_mask <- !base::is.na(ctx$.typei_role) & (ctx$operon_link | ctx$relative_gene == 0)
    source_rows <- ctx[member_mask, , drop = FALSE]
    roles_present <- base::sort(base::unique(stats::na.omit(source_rows$.typei_role)))
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
    if (base::nrow(source_rows)) {
      for (member_i in base::seq_len(base::nrow(source_rows))) {
        member <- source_rows[member_i, , drop = FALSE]
        partners <- source_rows[source_rows$locus_tag != member$locus_tag, , drop = FALSE]
        partner_label <- if (base::nrow(partners)) {
          partner_gap <- base::vapply(base::seq_len(base::nrow(partners)), function(j) {
            .dnmb_rebasefinder_intergenic_gap(member, partners[j, , drop = FALSE])
          }, numeric(1))
          base::paste(
            base::sprintf(
              "%s:%s:%s:%+dg:%sbp",
              partners$locus_tag,
              partners$.typei_role,
              partners$.typei_role_source,
              partners$.dnmb_order - member$.dnmb_order,
              base::ifelse(base::is.na(partner_gap), "NA", partner_gap)
            ),
            collapse = " | "
          )
        } else {
          NA_character_
        }
        context_rows[[base::length(context_rows) + 1L]] <- data.frame(
          query = member$locus_tag,
          typei_context_status = status,
          typei_context_roles = base::paste(roles_present, collapse = "/"),
          typei_context_partners = partner_label,
          typei_context_sources = source_label,
          typei_context_anchor = seed$locus_tag,
          typei_operon_supported = has_m && has_r && has_s,
          typei_context_window_bp = base::max(ctx$end, na.rm = TRUE) - base::min(ctx$start, na.rm = TRUE) + 1,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  context_tbl <- do.call(rbind, context_rows)
  context_tbl <- context_tbl[base::order(
    context_tbl$query,
    !context_tbl$typei_operon_supported,
    context_tbl$typei_context_window_bp
  ), , drop = FALSE]
  context_tbl <- context_tbl[!base::duplicated(context_tbl$query), , drop = FALSE]

  for (col in base::setdiff(base::names(context_tbl), "query")) {
    if (!col %in% base::names(hits)) hits[[col]] <- .dnmb_na_vector_like(context_tbl[[col]], base::nrow(hits))
  }
  hidx <- base::match(hits$query, context_tbl$query)
  include <- !base::is.na(hidx)
  for (col in base::setdiff(base::names(context_tbl), "query")) hits[[col]][include] <- context_tbl[[col]][hidx[include]]
  if (!"operon_id" %in% base::names(hits)) hits$operon_id <- NA_character_
  supported_operon <- include & hits$typei_operon_supported %in% TRUE
  hits$operon_id[supported_operon] <- base::paste0(
    "DNMB_TypeI_",
    context_tbl$typei_context_anchor[hidx[supported_operon]]
  )
  missing_operon <- include & (base::is.na(hits$operon_id) | !base::nzchar(base::as.character(hits$operon_id)))
  hits$operon_id[missing_operon] <- base::paste0(
    "DNMB_TypeI_",
    context_tbl$typei_context_anchor[hidx[missing_operon]]
  )

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

  candidate_idx <- base::which(
    genes_ord$locus_tag %in% context_tbl$query &
      !genes_ord$locus_tag %in% hits$query &
      !base::is.na(genes_ord$.typei_role) &
      genes_ord$.typei_role %in% c("M", "R", "S")
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
    ctx_idx <- base::match(add$query, context_tbl$query)
    add$typei_context_status <- context_tbl$typei_context_status[ctx_idx]
    add$typei_context_roles <- context_tbl$typei_context_roles[ctx_idx]
    add$typei_context_partners <- context_tbl$typei_context_partners[ctx_idx]
    add$typei_context_sources <- context_tbl$typei_context_sources[ctx_idx]
    add$typei_context_anchor <- context_tbl$typei_context_anchor[ctx_idx]
    add$typei_operon_supported <- context_tbl$typei_operon_supported[ctx_idx]
    add$typei_context_window_bp <- context_tbl$typei_context_window_bp[ctx_idx]
    add$operon_id <- base::paste0("DNMB_TypeI_", add$typei_context_anchor)
    for (col in base::setdiff(base::names(hits), base::names(add))) {
      add[[col]] <- .dnmb_na_vector_like(hits[[col]], base::nrow(add))
    }
    hits <- base::rbind(hits, add[, base::names(hits), drop = FALSE])
  }

  hits
}

.dnmb_rebasefinder_add_typeii_context <- function(hits,
                                                  genes,
                                                  max_operon_gap = 5000,
                                                  max_intervening = 1L,
                                                  max_neighbors = 3L) {
  if (base::is.null(hits) || !base::is.data.frame(hits)) return(hits)
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  hits <- .dnmb_rebasefinder_reset_operon_context(hits, "typeii")
  genes_ord <- .dnmb_rebasefinder_gene_order(genes)
  if (!base::nrow(genes_ord)) return(hits)

  anno <- .dnmb_rebasefinder_gene_annotation_text(genes_ord)
  seqs <- if ("translation" %in% base::names(genes_ord)) {
    genes_ord$translation
  } else {
    base::rep(NA_character_, base::nrow(genes_ord))
  }
  mtase_signals <- base::mapply(
    .dnmb_rebasefinder_mtase_sequence_signals,
    translation = seqs,
    annotation = anno,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  rease_signals <- base::mapply(
    .dnmb_rebasefinder_rease_sequence_signals,
    translation = seqs,
    annotation = anno,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  mtase_verified <- base::vapply(mtase_signals, `[[`, logical(1), "verified")
  mtase_annotated <- base::vapply(mtase_signals, `[[`, logical(1), "annotation_support")
  rease_motor <- base::vapply(rease_signals, `[[`, logical(1), "motor")
  rease_nuclease <- base::vapply(rease_signals, `[[`, logical(1), "nuclease")
  rease_annotated <- base::vapply(rease_signals, `[[`, logical(1), "annotation_support")
  rease_component <- base::vapply(rease_signals, function(x) x$component %||% NA_character_, character(1))
  rease_signal_label <- base::vapply(rease_signals, function(x) {
    if (!base::length(x$signals)) NA_character_ else base::paste(x$signals, collapse = "+")
  }, character(1))
  rease_signal_count <- base::vapply(rease_signals, function(x) base::length(x$signals), integer(1))
  annotation_flags <- .dnmb_rebasefinder_annotation_flags(anno)

  hit_idx <- base::match(genes_ord$locus_tag, hits$query)
  have_hit <- !base::is.na(hit_idx)
  genes_ord$.rm_family <- NA_character_
  genes_ord$.rm_role <- NA_character_
  genes_ord$.typing_eligible <- FALSE
  genes_ord$.strict_blast_support <- FALSE
  genes_ord$.structure_supported <- FALSE
  if (base::nrow(hits)) {
    genes_ord$.rm_family[have_hit] <- base::as.character(hits$family_id[hit_idx[have_hit]])
    genes_ord$.rm_role[have_hit] <- base::as.character(hits$enzyme_role[hit_idx[have_hit]])
    if ("typing_eligible" %in% base::names(hits)) {
      genes_ord$.typing_eligible[have_hit] <- hits$typing_eligible[hit_idx[have_hit]] %in% TRUE
    }
    strict_blast <- .dnmb_rebasefinder_strict_blast_support(hits)
    genes_ord$.strict_blast_support[have_hit] <- strict_blast[hit_idx[have_hit]]
    structure_supported <- .dnmb_rebasefinder_consistent_structure_support(hits)
    genes_ord$.structure_supported[have_hit] <- structure_supported[hit_idx[have_hit]]
  }

  known_typeii_m <- !base::is.na(genes_ord$.rm_family) & genes_ord$.rm_family == "Type II" &
    !base::is.na(genes_ord$.rm_role) & genes_ord$.rm_role == "M" &
    (genes_ord$.strict_blast_support | genes_ord$.structure_supported) &
    !annotation_flags$explicit_non_dna_mtase & !annotation_flags$other_defense
  mtase_family_compatible <- base::is.na(genes_ord$.rm_family) |
    !base::nzchar(genes_ord$.rm_family) | genes_ord$.rm_family == "Type II"
  mtase_role_compatible <- base::is.na(genes_ord$.rm_role) |
    !base::nzchar(genes_ord$.rm_role) | genes_ord$.rm_role == "M"
  motif_m <- mtase_verified & mtase_annotated &
    mtase_family_compatible & mtase_role_compatible &
    !annotation_flags$explicit_non_dna_mtase & !annotation_flags$other_defense
  anchor_idx <- base::which(known_typeii_m | motif_m)
  if (!base::length(anchor_idx)) {
    return(.dnmb_rebasefinder_add_motif_evidence(hits, genes))
  }

  memberships <- list()
  for (anchor_i in anchor_idx) {
    anchor <- genes_ord[anchor_i, , drop = FALSE]
    same_contig <- genes_ord$contig == anchor$contig
    same_direction <- base::is.na(anchor$direction) | base::is.na(genes_ord$direction) |
      genes_ord$direction == anchor$direction
    relative_gene <- genes_ord$.dnmb_order - anchor$.dnmb_order
    nearby <- same_contig & same_direction &
      base::abs(relative_gene) <= base::min(max_neighbors, max_intervening + 1L)
    nearby[anchor_i] <- TRUE
    nearby_idx <- base::which(nearby)
    if (base::length(nearby_idx)) {
      gap <- base::vapply(nearby_idx, function(i) {
        .dnmb_rebasefinder_intergenic_gap(anchor, genes_ord[i, , drop = FALSE])
      }, numeric(1))
      coherent_path <- base::vapply(nearby_idx, function(i) {
        .dnmb_rebasefinder_coherent_gene_path(genes_ord, anchor_i, i)
      }, logical(1))
      nearby_idx <- nearby_idx[(base::is.na(gap) | gap <= max_operon_gap) & coherent_path]
    }
    rease_family_compatible <- base::is.na(genes_ord$.rm_family[nearby_idx]) |
      !base::nzchar(genes_ord$.rm_family[nearby_idx]) |
      genes_ord$.rm_family[nearby_idx] == "Type II"
    rease_role_compatible <- base::is.na(genes_ord$.rm_role[nearby_idx]) |
      !base::nzchar(genes_ord$.rm_role[nearby_idx]) |
      genes_ord$.rm_role[nearby_idx] == "R"
    r_pool <- nearby_idx[
      nearby_idx != anchor_i & rease_family_compatible & rease_role_compatible &
        !annotation_flags$repair_nuclease[nearby_idx] &
        !annotation_flags$unrelated_helicase[nearby_idx] &
        !annotation_flags$other_defense[nearby_idx] &
        !annotation_flags$toxin_atpase[nearby_idx] &
        !annotation_flags$unrelated_metabolic[nearby_idx] &
        (rease_annotated[nearby_idx] | rease_nuclease[nearby_idx] | rease_motor[nearby_idx] |
           genes_ord$.strict_blast_support[nearby_idx] | genes_ord$.structure_supported[nearby_idx])
    ]
    if (!base::length(r_pool)) next

    strong_r <- r_pool[
      rease_annotated[r_pool] | rease_nuclease[r_pool] |
        genes_ord$.strict_blast_support[r_pool] | genes_ord$.structure_supported[r_pool]
    ]
    if (!base::length(strong_r)) next
    validated_r <- r_pool[
      rease_annotated[r_pool] |
        genes_ord$.strict_blast_support[r_pool] |
        genes_ord$.structure_supported[r_pool] |
        (rease_nuclease[r_pool] & (rease_motor[r_pool] | rease_signal_count[r_pool] >= 2L))
    ]
    validated_operon <- base::length(validated_r) > 0L
    split_motor <- r_pool[rease_motor[r_pool]]
    keep_r <- base::unique(c(strong_r, split_motor))
    split_r <- base::length(split_motor) && base::any(!split_motor %in% strong_r)
    status <- if (!validated_operon) {
      "candidate_mr_unvalidated"
    } else if (split_r) {
      "complete_m_split_r"
    } else {
      "complete_mr"
    }
    anchor_hit <- hit_idx[[anchor_i]]
    existing_operon <- if (!base::is.na(anchor_hit) && "operon_id" %in% base::names(hits)) {
      base::as.character(hits$operon_id[[anchor_hit]])
    } else {
      NA_character_
    }
    if (base::is.na(existing_operon) || !base::nzchar(existing_operon)) {
      existing_operon <- base::paste0("DNMB_TypeII_", anchor$locus_tag)
    }
    member_idx <- c(anchor_i, keep_r)
    r_components <- rease_component[keep_r]
    missing_component <- base::is.na(r_components) | !base::nzchar(r_components)
    r_components[missing_component & genes_ord$.strict_blast_support[keep_r]] <- "R_homology"
    missing_component <- base::is.na(r_components) | !base::nzchar(r_components)
    r_components[missing_component & genes_ord$.structure_supported[keep_r]] <- "R_structure"
    missing_component <- base::is.na(r_components) | !base::nzchar(r_components)
    r_components[missing_component & rease_motor[keep_r]] <- "R_motor"
    r_components[base::is.na(r_components) | !base::nzchar(r_components)] <- "R"
    member_components <- c(
      "M",
      r_components
    )
    for (i in member_idx) {
      is_anchor <- i == anchor_i
      component <- member_components[[base::match(i, member_idx)]]
      source <- if (is_anchor) {
        if (known_typeii_m[[i]]) "rebase_hit" else "mtase_motif_pair"
      } else if (rease_annotated[[i]] && rease_motor[[i]]) {
        "annotation+motor_motif"
      } else if (rease_annotated[[i]]) {
        "endonuclease_annotation"
      } else if (rease_nuclease[[i]]) {
        "nuclease_motif"
      } else if (genes_ord$.strict_blast_support[[i]]) {
        "rebase_homology"
      } else if (genes_ord$.structure_supported[[i]]) {
        "structure_homology"
      } else {
        "split_motor_motif"
      }
      partner_idx <- member_idx[member_idx != i]
      partners <- base::paste(
        base::sprintf(
          "%s:%s:%+dg",
          genes_ord$locus_tag[partner_idx],
          member_components[base::match(partner_idx, member_idx)],
          genes_ord$.dnmb_order[partner_idx] - genes_ord$.dnmb_order[[i]]
        ),
        collapse = " | "
      )
      memberships[[base::length(memberships) + 1L]] <- base::data.frame(
        query = genes_ord$locus_tag[[i]],
        anchor = anchor$locus_tag,
        role = if (is_anchor) "M" else "R",
        component = component,
        source = source,
        status = status,
        validated = validated_operon,
        partners = partners,
        operon_id = existing_operon,
        window_bp = base::max(genes_ord$end[member_idx], na.rm = TRUE) -
          base::min(genes_ord$start[member_idx], na.rm = TRUE) + 1,
        signal = if (is_anchor) {
          base::paste(mtase_signals[[i]]$signals, collapse = "+")
        } else {
          rease_signal_label[[i]]
        },
        stringsAsFactors = FALSE
      )
    }
  }
  if (!base::length(memberships)) {
    return(.dnmb_rebasefinder_add_motif_evidence(hits, genes))
  }
  context <- base::do.call(base::rbind, memberships)
  context <- context[base::order(
    context$query,
    !context$validated,
    context$window_bp
  ), , drop = FALSE]
  context <- context[!base::duplicated(context$query), , drop = FALSE]

  required_hit_cols <- list(
    operon_id = NA_character_,
    typeii_context_status = NA_character_,
    typeii_context_roles = NA_character_,
    typeii_context_partners = NA_character_,
    typeii_context_sources = NA_character_,
    typeii_operon_candidate = NA,
    typeii_operon_supported = NA,
    typeii_context_window_bp = NA_real_,
    operon_component = NA_character_
  )
  for (col in base::names(required_hit_cols)) {
    if (!col %in% base::names(hits)) hits[[col]] <- base::rep(required_hit_cols[[col]], base::nrow(hits))
  }

  existing <- base::match(hits$query, context$query)
  present <- !base::is.na(existing)
  if (base::any(present)) {
    ctx <- context[existing[present], , drop = FALSE]
    hits$typeii_context_status[present] <- ctx$status
    hits$typeii_context_roles[present] <- "M/R"
    hits$typeii_context_partners[present] <- ctx$partners
    hits$typeii_context_sources[present] <- ctx$source
    hits$typeii_operon_candidate[present] <- TRUE
    hits$typeii_operon_supported[present] <- ctx$validated
    hits$typeii_context_window_bp[present] <- ctx$window_bp
    hits$operon_component[present] <- ctx$component
    missing_operon <- present & (base::is.na(hits$operon_id) | !base::nzchar(base::as.character(hits$operon_id)))
    hits$operon_id[missing_operon] <- context$operon_id[existing[missing_operon]]
    missing_family <- present & (base::is.na(hits$family_id) | !base::nzchar(base::as.character(hits$family_id)))
    missing_role <- present & (base::is.na(hits$enzyme_role) | !base::nzchar(base::as.character(hits$enzyme_role)))
    rescued <- missing_family | missing_role
    hits$family_id[missing_family] <- "Type II"
    hits$enzyme_role[missing_role] <- context$role[existing[missing_role]]
    missing_label <- rescued & (base::is.na(hits$hit_label) | !base::nzchar(base::as.character(hits$hit_label)))
    hits$hit_label[missing_label] <- base::paste0(
      "typeII_context:",
      context$role[existing[missing_label]],
      ":",
      context$source[existing[missing_label]]
    )
    low_evidence <- present & (
      base::is.na(hits$evidence_mode) |
        !base::nzchar(base::as.character(hits$evidence_mode)) |
        hits$evidence_mode %in% c("annotation_only", "low_confidence")
    )
    hits$evidence_mode[low_evidence] <- "operon_context_split_R"
    hits$typing_eligible[rescued] <- FALSE
    context_support <- base::paste0(
      "Type II operon_context=", context$status[existing[present]],
      "; anchor=", context$anchor[existing[present]],
      "; component=", context$component[existing[present]],
      "; motifs=", base::ifelse(
        base::is.na(context$signal[existing[present]]),
        "none",
        context$signal[existing[present]]
      )
    )
    present_idx <- base::which(present)
    hits$support[present_idx] <- .dnmb_rebasefinder_strip_typeii_support(hits$support[present_idx])
    hits$support[present_idx] <- base::ifelse(
      base::is.na(hits$support[present_idx]) | !base::nzchar(base::as.character(hits$support[present_idx])),
      context_support,
      base::paste(hits$support[present_idx], context_support, sep = "; ")
    )
  }

  new_context <- context[!context$query %in% hits$query, , drop = FALSE]
  if (base::nrow(new_context)) {
    add <- base::data.frame(
      query = new_context$query,
      source = "rebasefinder",
      family_system = "REBASEfinder",
      family_id = "Type II",
      hit_label = base::paste0("typeII_context:", new_context$role, ":", new_context$source),
      enzyme_role = new_context$role,
      evidence_mode = "operon_context_split_R",
      substrate_label = NA_character_,
      support = base::paste0(
        "Type II operon_context=", new_context$status,
        "; anchor=", new_context$anchor,
        "; component=", new_context$component,
        "; motifs=", base::ifelse(base::is.na(new_context$signal), "none", new_context$signal)
      ),
      typing_eligible = FALSE,
      stringsAsFactors = FALSE
    )
    add$operon_id <- new_context$operon_id
    add$typeii_context_status <- new_context$status
    add$typeii_context_roles <- "M/R"
    add$typeii_context_partners <- new_context$partners
    add$typeii_context_sources <- new_context$source
    add$typeii_operon_candidate <- TRUE
    add$typeii_operon_supported <- new_context$validated
    add$typeii_context_window_bp <- new_context$window_bp
    add$operon_component <- new_context$component
    for (col in base::setdiff(base::names(hits), base::names(add))) {
      add[[col]] <- .dnmb_na_vector_like(hits[[col]], base::nrow(add))
    }
    for (col in base::setdiff(base::names(add), base::names(hits))) {
      hits[[col]] <- .dnmb_na_vector_like(add[[col]], base::nrow(hits))
    }
    hits <- base::rbind(hits, add[, base::names(hits), drop = FALSE])
  }
  .dnmb_rebasefinder_add_motif_evidence(hits, genes)
}

.dnmb_rebasefinder_typeiv_candidate_mask <- function(genes) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(genes)) return(logical())
  text <- .dnmb_rebasefinder_gene_annotation_text(genes)
  flags <- .dnmb_rebasefinder_annotation_flags(text)
  explicit_typeiv <- base::grepl(
    "type[ _-]?iv restriction|\\bmrr\\b|\\bmcr[abc]?\\b|\\bgmrsd\\b|modified[ -]?dna.*restriction|methylated[ -]?dna.*restriction|restriction.*modified[ -]?dna|modification-dependent restriction|pvruts1i|pvurts1i|sauusi|mspji",
    text,
    perl = TRUE
  )
  not_other_system <- !flags$other_defense & !flags$repair_nuclease &
    !flags$mobile_nuclease & !flags$non_rm_restriction_fold
  (explicit_typeiv | flags$typeiv_specific) & not_other_system
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
  context_candidate <- evidence %in% c(
    "operon_context", "operon_context_split_R", "annotation_motif_candidate",
    "structure_supported", "structure_only"
  ) | base::grepl("^operon_context", evidence)
  family %in% c("Type I", "Type II", "Type III", "Type IV") &
    role %in% c("M", "R", "S") &
    blast_missing &
    context_candidate
}

.dnmb_rebasefinder_add_split_r_evidence <- function(hits,
                                                     genes,
                                                     max_pair_gap = 250,
                                                     min_combined_reference_coverage = 0.80,
                                                     max_single_reference_coverage = 0.75,
                                                     min_query_coverage = 0.75) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(hits) || !base::nrow(genes) ||
      !all(c("query", "family_id", "enzyme_role") %in% base::names(hits)) ||
      !"locus_tag" %in% base::names(genes)) return(hits)
  for (column in base::intersect(c("start", "end"), base::names(genes))) {
    genes[[column]] <- suppressWarnings(base::as.numeric(genes[[column]]))
  }

  ensure_col <- function(column, value) {
    if (!column %in% base::names(hits)) {
      hits[[column]] <<- base::rep(value, base::nrow(hits))
    }
  }
  ensure_col("split_r_pair_candidate", FALSE)
  ensure_col("split_r_pair_id", NA_character_)
  ensure_col("split_r_partner", NA_character_)
  ensure_col("split_r_reference_subject", NA_character_)
  ensure_col("split_r_combined_reference_coverage", NA_real_)
  ensure_col("split_r_genomic_gap_bp", NA_real_)
  ensure_col("split_r_status", NA_character_)
  ensure_col("operon_id", NA_character_)
  ensure_col("support", NA_character_)

  hits$query <- .dnmb_module_clean_annotation_key(hits$query)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  gene_idx <- base::match(hits$query, genes$locus_tag)
  subject <- base::rep(NA_character_, base::nrow(hits))
  for (column in c("curated_blast_subject_id", "curated_blast_resolved_subject_id",
                   "supplemental_blast_match", "blast_match", "hit_label")) {
    if (!column %in% base::names(hits)) next
    value <- base::as.character(hits[[column]])
    use <- (base::is.na(subject) | !base::nzchar(subject)) &
      !base::is.na(value) & base::nzchar(value)
    subject[use] <- value[use]
  }
  ref_cov <- if ("blast_reference_coverage" %in% base::names(hits)) {
    suppressWarnings(base::as.numeric(hits$blast_reference_coverage))
  } else if ("curated_blast_reference_coverage" %in% base::names(hits)) {
    suppressWarnings(base::as.numeric(hits$curated_blast_reference_coverage))
  } else {
    base::rep(NA_real_, base::nrow(hits))
  }
  query_cov <- if ("blast_query_coverage" %in% base::names(hits)) {
    suppressWarnings(base::as.numeric(hits$blast_query_coverage))
  } else if ("curated_blast_query_coverage" %in% base::names(hits)) {
    suppressWarnings(base::as.numeric(hits$curated_blast_query_coverage))
  } else {
    base::rep(NA_real_, base::nrow(hits))
  }

  eligible <- hits$family_id %in% "Type I" & hits$enzyme_role %in% "R" &
    !base::is.na(subject) & base::nzchar(subject) &
    !base::is.na(ref_cov) & ref_cov > 0 & ref_cov <= max_single_reference_coverage &
    !base::is.na(query_cov) & query_cov >= min_query_coverage &
    !base::is.na(gene_idx)
  groups <- base::split(base::which(eligible), subject[eligible])
  groups <- groups[base::lengths(groups) >= 2L]
  if (!base::length(groups)) return(hits)

  candidates <- list()
  for (indices in groups) {
    pairs <- utils::combn(indices, 2L)
    for (column in base::seq_len(base::ncol(pairs))) {
      pair <- pairs[, column]
      first <- genes[gene_idx[pair[[1]]], , drop = FALSE]
      second <- genes[gene_idx[pair[[2]]], , drop = FALSE]
      same_contig <- !"contig" %in% base::names(genes) ||
        base::is.na(first$contig) || base::is.na(second$contig) ||
        base::identical(base::as.character(first$contig), base::as.character(second$contig))
      same_direction <- !"direction" %in% base::names(genes) ||
        base::is.na(first$direction) || base::is.na(second$direction) ||
        base::identical(base::as.character(first$direction), base::as.character(second$direction))
      if (!same_contig || !same_direction) next
      gap <- .dnmb_rebasefinder_intergenic_gap(first, second)
      combined <- base::min(1, base::sum(ref_cov[pair]))
      if ((!base::is.na(gap) && gap > max_pair_gap) ||
          combined < min_combined_reference_coverage) next
      candidates[[base::length(candidates) + 1L]] <- base::data.frame(
        first = pair[[1]],
        second = pair[[2]],
        combined = combined,
        gap = gap,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!base::length(candidates)) return(hits)
  candidates <- base::do.call(base::rbind, candidates)
  candidates <- candidates[base::order(-candidates$combined, candidates$gap), , drop = FALSE]
  used <- integer()
  for (row in base::seq_len(base::nrow(candidates))) {
    pair <- base::as.integer(candidates[row, c("first", "second")])
    if (base::any(pair %in% used)) next
    used <- c(used, pair)
    loci <- hits$query[pair]
    pair_id <- base::paste0("DNMB_TypeI_splitR_", base::paste(base::sort(loci), collapse = "_"))
    hits$split_r_pair_candidate[pair] <- TRUE
    hits$split_r_pair_id[pair] <- pair_id
    hits$split_r_partner[pair] <- base::rev(loci)
    hits$split_r_reference_subject[pair] <- subject[pair]
    hits$split_r_combined_reference_coverage[pair] <- candidates$combined[[row]]
    hits$split_r_genomic_gap_bp[pair] <- candidates$gap[[row]]
    hits$split_r_status[pair] <- "complementary_reference_fragments_review"
    supported <- if ("typei_operon_supported" %in% base::names(hits)) {
      hits$typei_operon_supported[pair] %in% TRUE
    } else {
      base::rep(FALSE, base::length(pair))
    }
    hits$operon_id[pair[!supported]] <- pair_id
    for (position in base::seq_along(pair)) {
      i <- pair[[position]]
      existing <- base::as.character(hits$support[[i]])
      if (!base::is.na(existing) && base::grepl("split_R_pair=", existing, fixed = TRUE)) next
      detail <- base::sprintf(
        "split_R_pair=%s; partner=%s; shared_reference=%s; combined_reference_coverage=%.3f; genomic_gap_bp=%s",
        pair_id,
        hits$split_r_partner[[i]],
        subject[[i]],
        candidates$combined[[row]],
        if (base::is.na(candidates$gap[[row]])) "NA" else base::as.character(candidates$gap[[row]])
      )
      hits$support[[i]] <- if (base::is.na(existing) || !base::nzchar(existing)) {
        detail
      } else {
        base::paste(existing, detail, sep = "; ")
      }
    }
  }
  hits
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

.dnmb_rebasefinder_reference_url <- function() {
  "ftp://ftp.neb.com/pub/rebase/protein_seqs.txt"
}

.dnmb_rebasefinder_reference_state_path <- function(cache_root = NULL) {
  base::file.path(
    .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE),
    "rebase_reference_state.rds"
  )
}

.dnmb_rebasefinder_reference_needs_refresh <- function(raw_path,
                                                       remote_metadata,
                                                       state = NULL,
                                                       force = FALSE) {
  if (base::isTRUE(force) || !.dnmb_nonempty_file(raw_path)) return(TRUE)
  if (base::is.null(remote_metadata) || !base::isTRUE(remote_metadata$ok)) return(FALSE)

  local_info <- base::file.info(raw_path)
  remote_size <- base::suppressWarnings(base::as.numeric(remote_metadata$content_length %||% NA_real_))
  if (!base::is.na(remote_size) && !base::identical(base::as.numeric(local_info$size), remote_size)) {
    return(TRUE)
  }

  if (base::is.list(state)) {
    for (field in c("last_modified", "content_length", "etag")) {
      old <- base::as.character(state[[field]] %||% NA_character_)[1]
      new <- base::as.character(remote_metadata[[field]] %||% NA_character_)[1]
      if (!base::is.na(new) && base::nzchar(new) && !base::identical(old, new)) return(TRUE)
    }
    return(FALSE)
  }

  remote_time <- base::suppressWarnings(base::as.POSIXct(
    remote_metadata$last_modified %||% NA_character_,
    format = "%a, %d %b %Y %H:%M:%S",
    tz = "GMT"
  ))
  !base::is.na(remote_time) && remote_time > local_info$mtime
}

.dnmb_rebasefinder_refresh_reference <- function(cache_root = NULL,
                                                  force = FALSE,
                                                  check_max_age_days = 1,
                                                  verbose = FALSE) {
  cache_dir <- .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE)
  raw_path <- base::file.path(cache_dir, "REBASE_protein_seqs.txt")
  rds_path <- base::file.path(cache_dir, "rebase_data.rds")
  fasta_path <- base::file.path(cache_dir, "rebase_db.fasta")
  state_path <- .dnmb_rebasefinder_reference_state_path(cache_root)
  state <- if (base::file.exists(state_path)) {
    tryCatch(base::readRDS(state_path), error = function(e) NULL)
  } else {
    NULL
  }

  state_age <- if (base::is.list(state) && !base::is.null(state$checked_at)) {
    checked <- base::suppressWarnings(base::as.POSIXct(state$checked_at, tz = "UTC"))
    base::as.numeric(base::difftime(base::Sys.time(), checked, units = "days"))
  } else {
    Inf
  }
  if (!base::isTRUE(force) && .dnmb_nonempty_file(raw_path) && .dnmb_nonempty_file(rds_path) &&
      .dnmb_nonempty_file(fasta_path) &&
      !base::is.na(state_age) && state_age <= base::as.numeric(check_max_age_days)[1]) {
    return(list(ok = TRUE, refreshed = FALSE, detail = "REBASE reference check is current.", state = state))
  }

  remote <- .dnmb_remote_asset_metadata(.dnmb_rebasefinder_reference_url(), insecure = FALSE)
  raw_needs_download <- .dnmb_rebasefinder_reference_needs_refresh(
    raw_path = raw_path,
    remote_metadata = remote,
    state = state,
    force = force
  )
  data_needs_rebuild <- .dnmb_nonempty_file(raw_path) &&
    (!.dnmb_nonempty_file(rds_path) || !.dnmb_nonempty_file(fasta_path) ||
       base::file.info(raw_path)$mtime > base::file.info(rds_path)$mtime)
  needs_refresh <- base::isTRUE(raw_needs_download) || base::isTRUE(data_needs_rebuild)
  next_state <- list(
    url = .dnmb_rebasefinder_reference_url(),
    checked_at = base::format(base::Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    last_modified = remote$last_modified %||% NA_character_,
    content_length = remote$content_length %||% NA_real_,
    etag = remote$etag %||% NA_character_,
    local_size = if (.dnmb_nonempty_file(raw_path)) base::file.info(raw_path)$size else NA_real_
  )

  if (!base::isTRUE(needs_refresh)) {
    if (.dnmb_nonempty_file(raw_path) && .dnmb_nonempty_file(rds_path)) {
      saved_state <- next_state
      if (!base::isTRUE(remote$ok)) {
        saved_state <- state
        if (base::is.list(saved_state)) saved_state$checked_at <- next_state$checked_at
      }
      if (base::is.list(saved_state)) {
        tryCatch(base::saveRDS(saved_state, state_path), error = function(e) NULL)
      }
      detail <- if (base::isTRUE(remote$ok)) "REBASE reference is current." else "Offline; using the validated cached REBASE reference."
      return(list(ok = TRUE, refreshed = FALSE, detail = detail, state = saved_state))
    }
    return(list(ok = FALSE, refreshed = FALSE, detail = remote$error %||% "REBASE reference is unavailable.", state = state))
  }

  if (!requireNamespace("DefenseViz", quietly = TRUE)) {
    return(list(ok = FALSE, refreshed = FALSE, detail = "DefenseViz is required to parse the REBASE reference.", state = state))
  }
  stage_raw <- if (base::isTRUE(raw_needs_download)) {
    base::tempfile("rebase-reference-", tmpdir = cache_dir, fileext = ".txt")
  } else {
    raw_path
  }
  stage_rds <- base::tempfile("rebase-reference-", tmpdir = cache_dir, fileext = ".rds")
  stage_fasta <- base::tempfile("rebase-reference-", tmpdir = cache_dir, fileext = ".fasta")
  stage_state <- base::tempfile("rebase-reference-", tmpdir = cache_dir, fileext = ".state.rds")
  cleanup <- c(
    if (base::isTRUE(raw_needs_download)) stage_raw else character(),
    stage_rds,
    stage_fasta,
    stage_state
  )
  base::on.exit(base::unlink(cleanup, force = TRUE), add = TRUE)

  expected_size <- base::suppressWarnings(base::as.numeric(remote$content_length %||% NA_real_))
  if (base::isTRUE(raw_needs_download)) {
    download <- .dnmb_download_asset(.dnmb_rebasefinder_reference_url(), stage_raw, insecure = FALSE)
    downloaded_size <- if (.dnmb_nonempty_file(stage_raw)) base::file.info(stage_raw)$size else NA_real_
    size_ok <- .dnmb_nonempty_file(stage_raw) &&
      (base::is.na(expected_size) || base::identical(base::as.numeric(downloaded_size), expected_size))
    if (!base::isTRUE(download$ok) || !base::isTRUE(size_ok)) {
      return(list(ok = .dnmb_nonempty_file(rds_path), refreshed = FALSE,
                  detail = download$error %||% "REBASE download was incomplete; the previous cache was preserved.", state = state))
    }
  }

  parsed <- tryCatch(DefenseViz::parse_rebase_sequences(stage_raw), error = function(e) e)
  required_columns <- c("enzyme_name", "sequence")
  if (inherits(parsed, "error") || !base::is.data.frame(parsed) || !base::nrow(parsed) ||
      !base::all(required_columns %in% base::names(parsed))) {
    detail <- if (inherits(parsed, "error")) conditionMessage(parsed) else "Parsed REBASE data failed schema validation."
    return(list(ok = .dnmb_nonempty_file(rds_path), refreshed = FALSE, detail = detail, state = state))
  }
  rds_written <- tryCatch({
    base::saveRDS(parsed, stage_rds)
    .dnmb_nonempty_file(stage_rds)
  }, error = function(e) FALSE)
  if (!base::isTRUE(rds_written)) {
    return(list(ok = .dnmb_nonempty_file(rds_path), refreshed = FALSE,
                detail = "Could not stage the parsed REBASE data; the previous cache was preserved.", state = state))
  }

  fasta_data <- parsed
  fasta_data$locus_tag <- fasta_data$enzyme_name
  fasta_data$translation <- fasta_data$sequence
  write_fasta <- base::get("write_fasta_for_blast", envir = base::asNamespace("DefenseViz"))
  fasta_result <- tryCatch(
    write_fasta(fasta_data, stage_fasta, id_col = "locus_tag", seq_col = "translation"),
    error = function(e) NULL
  )
  if (base::is.null(fasta_result) || !.dnmb_nonempty_file(stage_fasta)) {
    return(list(ok = .dnmb_nonempty_file(rds_path), refreshed = FALSE,
                detail = "Could not create a FASTA file from the refreshed REBASE data.", state = state))
  }

  next_state$local_size <- base::file.info(stage_raw)$size
  next_state$sequence_count <- base::nrow(parsed)
  state_written <- tryCatch({
    base::saveRDS(next_state, stage_state)
    staged_state <- base::readRDS(stage_state)
    base::is.list(staged_state) && base::identical(staged_state$sequence_count, base::nrow(parsed))
  }, error = function(e) FALSE)
  if (!base::isTRUE(state_written)) {
    return(list(ok = .dnmb_nonempty_file(rds_path), refreshed = FALSE,
                detail = "Could not stage REBASE reference state; the previous cache was preserved.", state = state))
  }

  replacements <- c(
    if (base::isTRUE(raw_needs_download)) c(raw_path = stage_raw) else character(),
    rds_path = stage_rds,
    fasta_path = stage_fasta,
    state_path = stage_state
  )
  destinations <- c(
    if (base::isTRUE(raw_needs_download)) raw_path else character(),
    rds_path,
    fasta_path,
    state_path
  )
  db_sidecars <- base::list.files(cache_dir, pattern = "^rebase_db[.]fasta[.]", full.names = TRUE)
  commit <- .dnmb_transactional_replace(
    staged_paths = replacements,
    destination_paths = destinations,
    retire_paths = db_sidecars
  )
  if (!base::isTRUE(commit$ok)) {
    return(list(
      ok = .dnmb_nonempty_file(rds_path),
      refreshed = FALSE,
      detail = commit$detail,
      state = state
    ))
  }
  if (base::exists("rmscan_rebase_db", envir = .GlobalEnv, inherits = FALSE)) {
    base::rm("rmscan_rebase_db", envir = .GlobalEnv)
  }
  if (base::isTRUE(verbose)) {
    message("[REBASEfinder] Refreshed REBASE reference: ", base::nrow(parsed), " sequences")
  }
  list(ok = TRUE, refreshed = TRUE, detail = base::paste0("sequences=", base::nrow(parsed)), state = next_state)
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
    refresh <- .dnmb_rebasefinder_refresh_reference(cache_root = cache_root, verbose = verbose)
    rebase_data <- .dnmb_rebasefinder_cached_rebase_data(cache_root)
    if (!base::isTRUE(refresh$ok) && isTRUE(verbose)) {
      message("[REBASEfinder] Could not load/download REBASE data for supplemental BLAST: ", refresh$detail)
    }
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
  for (col in c("supplemental_blast_family", "supplemental_blast_role")) {
    if (!col %in% names(hits)) hits[[col]] <- NA_character_
  }
  canonical_blast_cols <- list(
    blast_identity = NA_real_,
    blast_evalue = NA_real_,
    blast_bitscore = NA_real_,
    blast_length = NA_real_,
    rec_seq = NA_character_
  )
  for (col in base::names(canonical_blast_cols)) {
    if (!col %in% base::names(hits)) {
      hits[[col]] <- base::rep(canonical_blast_cols[[col]], base::nrow(hits))
    }
  }
  if (!"supplemental_blast_context_conflict" %in% names(hits)) {
    hits$supplemental_blast_context_conflict <- FALSE
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
    hits$supplemental_blast_family[h_idx] <- if ("rm_type" %in% names(b)) b$rm_type else NA_character_
    hits$supplemental_blast_role[h_idx] <- if ("subunit" %in% names(b)) b$subunit else NA_character_

    replace_blast <- is.na(hits$blast_identity[h_idx])
    replace_idx <- h_idx[replace_blast]
    bb <- b[replace_blast, , drop = FALSE]
    hits$blast_identity[replace_idx] <- bb$best_match_identity
    hits$blast_evalue[replace_idx] <- if ("evalue" %in% names(bb)) bb$evalue else NA_real_
    hits$blast_bitscore[replace_idx] <- if ("bitscore" %in% names(bb)) bb$bitscore else NA_real_
    hits$blast_length[replace_idx] <- if ("length" %in% names(bb)) bb$length else NA_real_
    trusted_override <- !is.na(bb$best_match_identity) & bb$best_match_identity >= 0.50
    old_family <- base::gsub("_", " ", base::as.character(hits$family_id[replace_idx]), fixed = TRUE)
    new_family <- base::gsub("_", " ", base::as.character(bb$rm_type), fixed = TRUE)
    new_family[!new_family %in% c("Type I", "Type II", "Type III", "Type IV")] <- NA_character_
    old_role <- base::as.character(hits$enzyme_role[replace_idx])
    new_role <- base::as.character(bb$subunit)
    new_role[!new_role %in% c("M", "R", "S", "RM")] <- NA_character_
    family_compatible <- base::is.na(old_family) | !base::nzchar(old_family) |
      base::is.na(new_family) | !base::nzchar(new_family) | old_family == new_family
    role_compatible <- base::is.na(old_role) | !base::nzchar(old_role) |
      base::is.na(new_role) | old_role == new_role
    context_compatible <- family_compatible & role_compatible
    conflict <- trusted_override & !context_compatible
    hits$supplemental_blast_context_conflict[replace_idx] <- conflict
    add_support <- base::sprintf(
      "supplemental_REBASE_BLAST=%s; identity=%s%%",
      bb$best_rebase_match,
      round(bb$best_match_identity * 100, 1)
    )
    add_support[conflict] <- base::paste0(
      add_support[conflict],
      "; context_conflict=retained ", old_family[conflict], "/", old_role[conflict],
      " over ", new_family[conflict], "/", new_role[conflict]
    )
    hits$support[replace_idx] <- base::ifelse(
      is.na(hits$support[replace_idx]) | !nzchar(hits$support[replace_idx]),
      add_support,
      base::paste(hits$support[replace_idx], add_support, sep = "; ")
    )
  }

  list(hits = hits, status = .dnmb_rebasefinder_status_row(
    "supplemental_rebase_blast", "ok",
    base::sprintf("Supplemental REBASE BLAST matched %d/%d operon-context candidates", nrow(best), nrow(query_tbl))
  ))
}

.dnmb_rebasefinder_write_augmented_results <- function(output_dir, hits) {
  out_dir <- output_dir
  base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tsv <- base::file.path(out_dir, "DNMB_REBASEfinder_augmented_hits.tsv")
  utils::write.table(hits, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  xlsx <- base::file.path(out_dir, "DNMB_REBASEfinder_augmented_hits.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch({
      wb <- openxlsx::createWorkbook()
      tables <- if ("curation_tier" %in% base::names(hits)) {
        list(
          Curated = hits[hits$curation_keep %in% TRUE, , drop = FALSE],
          Review = hits[hits$curation_tier %in% "review", , drop = FALSE],
          Other_defense = hits[hits$curation_tier %in% "other_defense", , drop = FALSE],
          Excluded_noise = hits[hits$curation_tier %in% "excluded_noise", , drop = FALSE],
          All_evidence = hits
        )
      } else {
        list(REBASEfinder_hits = hits)
      }
      for (sheet in base::names(tables)) {
        tbl <- tables[[sheet]]
        openxlsx::addWorksheet(wb, sheet)
        openxlsx::writeData(wb, sheet, tbl, withFilter = base::nrow(tbl) > 0L)
        dnmb_write_hyperlink_columns(wb, sheet, tbl, startRow = 1L, startCol = 1L)
        .dnmb_rebasefinder_style_workbook_sheet(wb, sheet, tbl)
        openxlsx::freezePane(wb, sheet, firstRow = TRUE)
      }
      openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)
    }, error = function(e) NULL)
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
  operon_context <- !base::is.na(evidence) & base::grepl("^operon_context", evidence)
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
  roots <- c(
    output_dir,
    base::file.path(output_dir, "dnmb_module_rebasefinder")
  )
  if (!base::is.null(structure_validation_path) && !base::is.na(structure_validation_path) && base::nzchar(structure_validation_path)) {
    validation_root <- if (base::dir.exists(structure_validation_path)) {
      structure_validation_path
    } else {
      base::dirname(structure_validation_path)
    }
    roots <- c(
      roots,
      validation_root,
      base::file.path(validation_root, "dnmb_module_rebasefinder")
    )
  }

  roots <- base::unique(roots)
  independent_names <- c(
    "query_structures",
    "alphafold_query_structures",
    "esmfold_query_structures",
    "rebasefinder_query_structures"
  )
  independent_dirs <- base::unlist(base::lapply(
    roots,
    function(root) base::file.path(root, independent_names)
  ), use.names = FALSE)
  managed_dirs <- base::file.path(roots, "promod3_query_structures")
  # Independently predicted or user-supplied structures must win over managed
  # homology models, even when they live beside a validation result elsewhere.
  dirs <- c(independent_dirs, managed_dirs)
  base::unique(dirs[base::dir.exists(dirs)])
}

.dnmb_rebasefinder_homology_structure_path_map <- function(tbl) {
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!base::nrow(tbl)) return(stats::setNames(character(), character()))
  id_col <- base::intersect(c("query", "locus_tag"), base::names(tbl))
  path_col <- base::intersect(
    c("homology_model_path", "REBASEfinder_homology_model_path"),
    base::names(tbl)
  )
  if (!base::length(id_col) || !base::length(path_col)) {
    return(stats::setNames(character(), character()))
  }
  ids <- .dnmb_module_clean_annotation_key(tbl[[id_col[[1]]]])
  paths <- base::as.character(tbl[[path_col[[1]]]])
  keep <- !base::is.na(ids) & base::nzchar(ids) &
    !base::is.na(paths) & base::nzchar(paths) & base::file.exists(paths)
  if (!base::any(keep)) return(stats::setNames(character(), character()))
  paths <- paths[keep]
  ids <- ids[keep]
  stats::setNames(paths[!base::duplicated(ids)], ids[!base::duplicated(ids)])
}

.dnmb_rebasefinder_structure_file_for_query <- function(query, dirs, explicit_paths = NULL) {
  if (base::is.na(query) || !base::nzchar(query)) return(NA_character_)
  safe_query <- base::gsub("[^A-Za-z0-9_.-]+", "_", query)
  hit <- character()
  if (base::length(dirs)) {
    filenames <- base::paste0(
      base::rep(c(query, safe_query), each = 2L),
      c(".pdb", ".cif")
    )
    candidates <- base::unlist(base::lapply(
      dirs,
      function(dir) base::file.path(dir, filenames)
    ), use.names = FALSE)
    candidates <- base::unique(candidates)
    hit <- candidates[base::file.exists(candidates)]
  }
  if (base::length(hit)) {
    return(base::normalizePath(hit[[1]], winslash = "/", mustWork = FALSE))
  }
  explicit_names <- .dnmb_module_clean_annotation_key(base::names(explicit_paths))
  explicit_paths <- base::as.character(explicit_paths)
  explicit_idx <- base::match(.dnmb_module_clean_annotation_key(query), explicit_names)
  if (!base::is.na(explicit_idx)) {
    explicit <- explicit_paths[[explicit_idx]]
    if (!base::is.na(explicit) && base::nzchar(explicit) && base::file.exists(explicit)) {
      return(base::normalizePath(explicit, winslash = "/", mustWork = FALSE))
    }
  }
  NA_character_
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
  explicit_paths <- .dnmb_rebasefinder_homology_structure_path_map(hits)
  structure_file <- base::vapply(
    records$query,
    .dnmb_rebasefinder_structure_file_for_query,
    character(1),
    dirs = structure_dirs,
    explicit_paths = explicit_paths
  )
  structure_file_exists <- !base::is.na(structure_file) & base::nzchar(structure_file) & base::file.exists(structure_file)
  structure_model_coverage <- base::vapply(base::seq_len(base::nrow(records)), function(i) {
    if (!structure_file_exists[[i]] ||
        !base::grepl("[.]pdb$", structure_file[[i]], ignore.case = TRUE)) {
      return(NA_real_)
    }
    ca <- .dnmb_rebasefinder_read_pdb_ca(structure_file[[i]])
    if (!base::nrow(ca)) return(NA_real_)
    alignment <- .dnmb_rebasefinder_align_structure_chain(records$sequence[[i]], ca)
    if (base::is.null(alignment)) NA_real_ else alignment$sequence_coverage
  }, numeric(1))
  structure_model_scope <- base::ifelse(
    !structure_file_exists, "no_structure",
    base::ifelse(
      !base::is.finite(structure_model_coverage), "coverage_unknown",
      base::ifelse(
        structure_model_coverage >= 0.80, "near_full_model",
        base::ifelse(structure_model_coverage >= 0.30, "partial_model", "short_fragment_model")
      )
    )
  )
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
		  model_coverage_usable <- !base::is.finite(structure_model_coverage) |
		    structure_model_coverage >= 0.30
		  supported <- raw_supported & structure_consistency$structure_candidate_consistent &
		    model_coverage_usable
		  mismatch <- foldseek_hit_present & raw_supported & !supported
	  best_status[mismatch & !structure_consistency$structure_family_consistent] <- "structure_family_mismatch"
	  best_status[mismatch & structure_consistency$structure_family_consistent &
	                !structure_consistency$structure_role_consistent] <- "structure_role_mismatch"

		  coverage_status <- base::rep("structure_missing", base::nrow(records))
		  coverage_status[structure_file_exists] <- "structure_available_no_foldseek_hit"
		  coverage_status[foldseek_hit_present] <- "foldseek_checked_no_pass"
		  coverage_status[mismatch & !structure_consistency$structure_family_consistent] <-
		    "foldseek_family_mismatch"
		  coverage_status[mismatch & structure_consistency$structure_family_consistent &
		                    !structure_consistency$structure_role_consistent] <-
		    "foldseek_role_mismatch"
		  coverage_status[raw_supported & !model_coverage_usable] <- "foldseek_fragment_too_short"
		  coverage_status[supported] <- "foldseek_supported"
		  coverage_status[supported & structure_model_scope == "partial_model"] <-
		    "foldseek_supported_partial_model"

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
	    structure_model_coverage = structure_model_coverage,
	    structure_model_scope = structure_model_scope,
		    foldseek_hit_present = foldseek_hit_present,
	    structure_supported = supported,
	    structure_status = best_status,
	    structure_best_hit = best_hit,
	    structure_family = best_family,
	    structure_role = best_role,
	    structure_chain_role = best_chain_role,
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
  ensure_col("structure_model_coverage", NA_real_)
  ensure_col("structure_model_scope", NA_character_)
  ensure_col("foldseek_hit_present", FALSE)
  ensure_col("structure_supported", FALSE)
  ensure_col("structure_coverage_status", NA_character_)
  ensure_col("structure_family", NA_character_)
  ensure_col("structure_role", NA_character_)
  ensure_col("structure_chain_role", NA_character_)
  ensure_col("structure_family_consistent", NA)
  ensure_col("structure_role_consistent", NA)
  ensure_col("structure_candidate_consistent", NA)

  cov <- coverage_tbl[idx[ok], , drop = FALSE]
	  cov_status <- base::as.character(cov$coverage_status)
	  mapped_status <- base::rep("structure_missing", base::length(cov_status))
	  mapped_status[cov_status == "structure_available_no_foldseek_hit"] <-
	    "structure_available_no_foldseek_hit"
	  mapped_status[cov_status == "foldseek_checked_no_pass"] <- "structure_checked_no_pass"
	  mapped_status[cov_status == "foldseek_family_mismatch"] <- "structure_family_mismatch"
	  mapped_status[cov_status == "foldseek_role_mismatch"] <- "structure_role_mismatch"
	  mapped_status[cov_status == "foldseek_fragment_too_short"] <- "structure_fragment_too_short"
	  mapped_status[cov_status == "foldseek_supported"] <- "structure_supported"
	  mapped_status[cov_status == "foldseek_supported_partial_model"] <-
	    "structure_supported_partial_model"

	  status_missing <- base::is.na(hits$structure_status[ok]) | !base::nzchar(base::as.character(hits$structure_status[ok]))
	  repl_idx <- base::which(ok)[status_missing]
	  hits$structure_status[repl_idx] <- mapped_status[status_missing]

		  pass <- cov_status %in% c("foldseek_supported", "foldseek_supported_partial_model")
		  pass_idx <- base::which(ok)
		  if (base::any(pass)) hits$structure_status[pass_idx[pass]] <- "structure_supported"
		  partial_pass <- pass & cov_status == "foldseek_supported_partial_model"
		  if (base::any(partial_pass)) {
		    hits$structure_status[pass_idx[partial_pass]] <- "structure_supported_partial_model"
		  }
	  not_supported <- !pass
	  if (base::any(not_supported)) hits$structure_status[pass_idx[not_supported]] <- mapped_status[not_supported]
	  hits$structure_pass[pass_idx] <- pass
	  hits$structure_supported[pass_idx] <- pass
	  hits$structure_coverage_status[pass_idx] <- cov_status

		  if ("structure_file" %in% base::names(cov)) hits$structure_file[pass_idx] <- base::as.character(cov$structure_file)
		  if ("structure_file_exists" %in% base::names(cov)) hits$structure_file_exists[pass_idx] <- cov$structure_file_exists
		  if ("structure_model_coverage" %in% base::names(cov)) hits$structure_model_coverage[pass_idx] <- cov$structure_model_coverage
		  if ("structure_model_scope" %in% base::names(cov)) hits$structure_model_scope[pass_idx] <- base::as.character(cov$structure_model_scope)
		  if ("foldseek_hit_present" %in% base::names(cov)) hits$foldseek_hit_present[pass_idx] <- cov$foldseek_hit_present
		  if ("structure_family" %in% base::names(cov)) hits$structure_family[pass_idx] <- base::as.character(cov$structure_family)
		  if ("structure_role" %in% base::names(cov)) hits$structure_role[pass_idx] <- base::as.character(cov$structure_role)
		  if ("structure_chain_role" %in% base::names(cov)) hits$structure_chain_role[pass_idx] <- base::as.character(cov$structure_chain_role)
		  if ("structure_family_consistent" %in% base::names(cov)) hits$structure_family_consistent[pass_idx] <- cov$structure_family_consistent
	  if ("structure_role_consistent" %in% base::names(cov)) hits$structure_role_consistent[pass_idx] <- cov$structure_role_consistent
		  if ("structure_candidate_consistent" %in% base::names(cov)) hits$structure_candidate_consistent[pass_idx] <- cov$structure_candidate_consistent
		  if ("support" %in% base::names(hits) && "structure_model_coverage" %in% base::names(cov)) {
		    partial <- base::is.finite(cov$structure_model_coverage) & cov$structure_model_coverage < 0.80
		    for (i in base::which(partial)) {
		      row_idx <- pass_idx[[i]]
		      old <- base::as.character(hits$support[[row_idx]])
		      if (!base::is.na(old) && base::grepl("model_cov_fullseq=", old, fixed = TRUE)) next
		      note <- base::sprintf(
		        "model_cov_fullseq=%.3f; structure_scope=%s",
		        cov$structure_model_coverage[[i]], cov$structure_model_scope[[i]]
		      )
		      hits$support[[row_idx]] <- if (base::is.na(old) || !base::nzchar(old)) note else base::paste(old, note, sep = "; ")
		    }
		  }

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

.dnmb_rebasefinder_first_motif <- function(sequence, pattern) {
  sequence <- .dnmb_rebasefinder_normalize_protein(sequence)
  if (base::is.na(sequence)) {
    return(list(found = FALSE, match = NA_character_, position = NA_integer_))
  }
  hit <- base::regexpr(pattern, sequence, perl = TRUE)
  if (hit[[1]] < 1L) {
    return(list(found = FALSE, match = NA_character_, position = NA_integer_))
  }
  len <- base::attr(hit, "match.length")[[1]]
  list(
    found = TRUE,
    match = base::substr(sequence, hit[[1]], hit[[1]] + len - 1L),
    position = base::as.integer(hit[[1]])
  )
}

.dnmb_rebasefinder_mtase_annotation <- function(text) {
  text <- base::tolower(base::as.character(text)[1])
  if (base::is.na(text) || !base::nzchar(text)) return(FALSE)
  explicit_dna <- base::grepl(
    "dna[^|;]{0,35}(methylase|methyltransferase)|(?:methylase|methyltransferase)[^|;]{0,35}dna|site-specific dna-methyl|restriction[- ]modification[^|;]{0,35}methyl|\\bhsdm\\b|\\bmod subunit\\b|methyltransf[ _-]?d12|pf02086|n[46][ -]?(adenine|cytosine)[^|;]{0,25}methyl",
    text,
    perl = TRUE
  )
  excluded <- base::grepl(
    "rrna|ribosomal rna|trna|rna methyl|protein methyl|chemotaxis protein|ubiquinone|ribose methyl",
    text,
    perl = TRUE
  )
  explicit_dna && !excluded
}

.dnmb_rebasefinder_mtase_sequence_signals <- function(translation, annotation = NA_character_) {
  seq <- .dnmb_rebasefinder_normalize_protein(translation)
  empty <- base::is.na(seq)
  if (empty) {
    return(list(
      role = NA_character_, verified = FALSE, class = NA_character_, confidence = "none",
      signals = character(), sam_match = NA_character_, sam_position = NA_integer_,
      catalytic_match = NA_character_, catalytic_position = NA_integer_, annotation_support = FALSE,
      architecture = NA_character_, amino_pair = FALSE, c5_pair = FALSE,
      mmei_pair = FALSE
    ))
  }

  motif <- function(name) {
    .dnmb_rebasefinder_first_motif(seq, .dnmb_rebasefinder_motif_pattern(name))
  }
  positions <- function(name) {
    .dnmb_rebasefinder_motif_all_positions(seq, .dnmb_rebasefinder_motif_pattern(name))
  }
  sam <- motif("SAM")
  amino <- motif("Amino-IV")
  c5_pc <- motif("N5C-PC")
  c5_env <- motif("N5C-ENV")
  c5_qrr <- motif("N5C-QxRxR")
  mmei_x <- motif("MmeI-X")
  mmei_i <- motif("MmeI-I")
  plausible_length <- base::nchar(seq) >= 150L && base::nchar(seq) <= 2000L
  amino_pair <- isTRUE(sam$found) && isTRUE(amino$found)
  c5_pair <- .dnmb_rebasefinder_motif_ordered(
    list(positions("N5C-PC"), positions("N5C-ENV")),
    max_span = 250L
  )
  c5_extended <- c5_pair && .dnmb_rebasefinder_motif_ordered(
    list(positions("N5C-PC"), positions("N5C-ENV"), positions("N5C-QxRxR")),
    max_span = 450L
  )
  mmei_pair <- .dnmb_rebasefinder_motif_ordered(
    list(positions("MmeI-X"), positions("MmeI-I"), positions("Amino-IV")),
    max_span = 650L
  )
  catalytic <- if (isTRUE(amino$found)) amino else c5_pc
  verified <- plausible_length && (amino_pair || c5_pair || mmei_pair)
  annotation_support <- .dnmb_rebasefinder_mtase_annotation(annotation)
  cls <- if (mmei_pair) {
    "N6A_gamma_MmeI_like"
  } else if (isTRUE(amino$found)) {
    "N6A_or_N4C_amino_MTase"
  } else if (c5_pair) {
    "N5C_C5_MTase"
  } else {
    NA_character_
  }
  signal_flags <- c(
    SAM = sam$found,
    `Amino-IV` = amino$found,
    `N5C-PC` = c5_pc$found,
    `N5C-ENV` = c5_env$found,
    `N5C-QxRxR` = c5_qrr$found,
    `MmeI-X` = mmei_x$found,
    `MmeI-I` = mmei_i$found
  )
  architecture <- if (mmei_pair) {
    "MmeI_gamma_MTase_X-I-IV"
  } else if (c5_extended) {
    "C5_PCQ-ENV-QxRxR"
  } else if (c5_pair) {
    "C5_PCQ-ENV"
  } else if (amino_pair) {
    "amino_MTase_I-IV"
  } else {
    NA_character_
  }
  sam_call <- if (isTRUE(sam$found)) sam else if (isTRUE(mmei_i$found)) mmei_i else sam
  list(
    role = if (verified) "M" else NA_character_,
    verified = verified,
    class = cls,
    confidence = if (verified && annotation_support) "high" else if (verified) "medium" else "none",
    signals = base::names(signal_flags)[signal_flags],
    sam_match = sam_call$match,
    sam_position = sam_call$position,
    catalytic_match = catalytic$match,
    catalytic_position = catalytic$position,
    annotation_support = annotation_support,
    architecture = architecture,
    amino_pair = amino_pair,
    c5_pair = c5_pair,
    mmei_pair = mmei_pair
  )
}

.dnmb_rebasefinder_rease_annotation <- function(text) {
  text <- base::tolower(base::as.character(text)[1])
  if (base::is.na(text) || !base::nzchar(text)) return(FALSE)
  explicit <- base::grepl(
    "restriction[- ]?(modification )?(endonuclease|enzyme)|restriction subunit|atp-dependent endonuclease|\\bhsdr\\b|\\bres subunit\\b|modification-dependent restriction|methylated[- ]dna.*restriction",
    text,
    perl = TRUE
  )
  excluded <- base::grepl(
    "crispr|\\bcas[0-9a-z]*\\b|transposase|integrase|recombinase|rna-guided|homing endonuclease|endonuclease (iii|iv|v|viii)|exonuclease|ribonuclease|dnase i|uvr[abc]|\\bmut[ls]\\b|\\brec[abcdgq]\\b|add[ab] helicase",
    text,
    perl = TRUE
  )
  explicit && !excluded
}

.dnmb_rebasefinder_rease_sequence_signals <- function(translation, annotation = NA_character_) {
  seq <- .dnmb_rebasefinder_normalize_protein(translation)
  if (base::is.na(seq)) {
    return(list(role = NA_character_, signals = character(), motor = FALSE,
                raw_signals = character(), nuclease = FALSE,
                annotation_support = FALSE, component = NA_character_, confidence = "none",
                architecture = NA_character_, mmei_fusion = FALSE))
  }
  motif_positions <- function(name) {
    .dnmb_rebasefinder_motif_all_positions(seq, .dnmb_rebasefinder_motif_pattern(name))
  }
  positions <- list(
    `P-loop` = base::sort(base::unique(c(motif_positions("P-loop"), motif_positions("ResIII-WA")))),
    `Walker-B` = base::sort(base::unique(c(motif_positions("HsdR-WB"), motif_positions("ResIII-WB")))),
    `Motif-III` = base::sort(base::unique(c(motif_positions("HsdR-MIII"), motif_positions("ResIII-MIII")))),
    `PD-ExK` = base::sort(base::unique(c(motif_positions("HsdR-PD"), motif_positions("PD-ExK"), motif_positions("ResIII-PD")))),
    HNH = motif_positions("HNH"),
    `GIY-YIG` = motif_positions("GIY-YIG"),
    `PLD-HKD` = motif_positions("PLD-HKD"),
    `MmeI-PDExK` = motif_positions("MmeI-PDExK")
  )
  raw_flags <- base::vapply(positions, base::length, integer(1)) > 0L
  motor <- .dnmb_rebasefinder_motif_ordered(
    positions[c("P-loop", "Walker-B", "Motif-III")],
    max_span = 750L
  )
  mmei_fusion <- .dnmb_rebasefinder_motif_ordered(
    list(
      motif_positions("MmeI-X"),
      motif_positions("MmeI-I"),
      motif_positions("Amino-IV")
    ),
    max_span = 650L
  )
  annotation_support <- .dnmb_rebasefinder_rease_annotation(annotation)
  annotation_text <- base::tolower(base::as.character(annotation)[1])
  if (base::is.na(annotation_text)) annotation_text <- ""
  hnh_domain <- base::grepl("\\bhnh\\b|h-n-h|pf01844|pf13391", annotation_text, perl = TRUE)
  giy_domain <- base::grepl("giy[-_ ]?yig|pf01541|pf02820", annotation_text, perl = TRUE)
  pld_domain <- base::grepl("pf13091|pldc_2|pld[- ]like|phospholipase d/nuclease|hkd", annotation_text, perl = TRUE)
  mmei_domain <- mmei_fusion || base::grepl(
    "\\bmmei\\b|pf2046[5-7]|pf20473|ipr04681[689]",
    annotation_text,
    perl = TRUE
  )
  validated_nuclease <- c(
    `PD-ExK` = raw_flags[["PD-ExK"]] && annotation_support,
    HNH = raw_flags[["HNH"]] && hnh_domain,
    `GIY-YIG` = raw_flags[["GIY-YIG"]] && (giy_domain || annotation_support),
    `PLD-HKD-pair` = base::length(positions[["PLD-HKD"]]) >= 2L && pld_domain,
    `MmeI-PDExK` = raw_flags[["MmeI-PDExK"]] && mmei_domain
  )
  nuclease <- base::any(validated_nuclease)
  signals <- c(
    if (motor) c("P-loop", "Walker-B", "Motif-III") else character(),
    base::names(validated_nuclease)[validated_nuclease]
  )
  component <- if (motor && (nuclease || annotation_support)) {
    if (nuclease) "R_motor_nuclease" else "R_motor_annotation"
  } else if (motor) {
    "R_motor"
  } else if (nuclease) {
    "R_nuclease"
  } else if (annotation_support) {
    "R_annotation"
  } else {
    NA_character_
  }
  architecture <- if (mmei_fusion && validated_nuclease[["MmeI-PDExK"]]) {
    "MmeI_fused_REase-MTase"
  } else if (motor && nuclease) {
    "SF2_motor_nuclease"
  } else if (motor) {
    "SF2_motor"
  } else if (nuclease) {
    "nuclease_motif_domain"
  } else {
    NA_character_
  }
  list(
    role = if (!base::is.na(component)) "R" else NA_character_,
    signals = signals,
    raw_signals = base::names(raw_flags)[raw_flags],
    motor = motor,
    nuclease = nuclease,
    annotation_support = annotation_support,
    component = component,
    confidence = if (motor && nuclease && annotation_support) "high" else if (nuclease || (motor && annotation_support)) "medium" else if (annotation_support) "annotation_only" else "none",
    architecture = architecture,
    mmei_fusion = mmei_fusion
  )
}

.dnmb_rebasefinder_add_motif_evidence <- function(hits, genes) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  idx <- base::match(.dnmb_module_clean_annotation_key(hits$query), genes$locus_tag)
  anno <- .dnmb_rebasefinder_gene_annotation_text(genes)
  seqs <- if ("translation" %in% base::names(genes)) genes$translation else base::rep(NA_character_, base::nrow(genes))
  mt <- base::lapply(base::seq_len(base::nrow(hits)), function(i) {
    if (base::is.na(idx[[i]])) return(.dnmb_rebasefinder_mtase_sequence_signals(NA_character_))
    .dnmb_rebasefinder_mtase_sequence_signals(seqs[[idx[[i]]]], anno[[idx[[i]]]])
  })
  re <- base::lapply(base::seq_len(base::nrow(hits)), function(i) {
    if (base::is.na(idx[[i]])) return(.dnmb_rebasefinder_rease_sequence_signals(NA_character_))
    .dnmb_rebasefinder_rease_sequence_signals(seqs[[idx[[i]]]], anno[[idx[[i]]]])
  })
  hits$mtase_motif_verified <- base::vapply(mt, `[[`, logical(1), "verified")
  hits$mtase_motif_class <- base::vapply(mt, function(x) x$class %||% NA_character_, character(1))
  hits$mtase_motif_confidence <- base::vapply(mt, `[[`, character(1), "confidence")
  hits$mtase_motif_architecture <- base::vapply(mt, function(x) x$architecture %||% NA_character_, character(1))
  hits$mtase_sam_motif <- base::vapply(mt, function(x) x$sam_match %||% NA_character_, character(1))
  hits$mtase_sam_position <- base::vapply(mt, function(x) x$sam_position %||% NA_integer_, integer(1))
  hits$mtase_catalytic_motif <- base::vapply(mt, function(x) x$catalytic_match %||% NA_character_, character(1))
  hits$mtase_catalytic_position <- base::vapply(mt, function(x) x$catalytic_position %||% NA_integer_, integer(1))
  hits$rease_motif_hits <- base::vapply(re, function(x) {
    if (!base::length(x$signals)) NA_character_ else base::paste(x$signals, collapse = "+")
  }, character(1))
  hits$rease_motif_hits_raw <- base::vapply(re, function(x) {
    if (!base::length(x$raw_signals)) NA_character_ else base::paste(x$raw_signals, collapse = "+")
  }, character(1))
  hits$rease_motif_verified <- base::vapply(re, function(x) {
    base::isTRUE(x$nuclease) || (base::isTRUE(x$motor) && base::isTRUE(x$annotation_support))
  }, logical(1))
  hits$rease_motif_architecture <- base::vapply(
    re,
    function(x) x$architecture %||% NA_character_,
    character(1)
  )
  hits$rease_motor_verified <- base::vapply(re, `[[`, logical(1), "motor")
  hits$rease_nuclease_verified <- base::vapply(re, `[[`, logical(1), "nuclease")
  hits$rease_motif_confidence <- base::vapply(re, `[[`, character(1), "confidence")
  hits$rease_operon_component_raw <- base::vapply(re, function(x) x$component %||% NA_character_, character(1))
  role <- if ("enzyme_role" %in% base::names(hits)) base::as.character(hits$enzyme_role) else NA_character_
  hits$motif_role_conflict <- !base::is.na(hits$rease_operon_component_raw) & role %in% c("M", "S")
  hits$rease_operon_component <- hits$rease_operon_component_raw
  hits$rease_operon_component[hits$motif_role_conflict] <- NA_character_
  hits
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
  # DefenseViz reads this file with read.delim(). Embedded tabs/newlines in
  # annotation fields otherwise become fake genes and invalidate gene order.
  clean_field <- function(x) {
    if (base::is.list(x) && !base::is.data.frame(x)) {
      x <- base::vapply(x, function(value) {
        value <- base::as.character(value)
        value <- value[!base::is.na(value) & base::nzchar(value)]
        if (!base::length(value)) NA_character_ else base::paste(value, collapse = "; ")
      }, character(1))
    }
    if (base::is.factor(x)) x <- base::as.character(x)
    if (!base::is.character(x)) return(x)
    x <- base::gsub("[\r\n\t]+", " ", x, perl = TRUE)
    x <- base::gsub('"', "'", x, fixed = TRUE)
    x <- base::gsub("[[:space:]]{2,}", " ", x, perl = TRUE)
    base::trimws(x)
  }
  genes_df[] <- base::lapply(genes_df, clean_field)
  base::names(genes_df) <- base::gsub("[\r\n\t]+", " ", base::names(genes_df), perl = TRUE)

  tsv_path <- base::file.path(output_dir, "rebasefinder_input.tsv")
  utils::write.table(
    genes_df,
    file = tsv_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )

  roundtrip <- base::tryCatch(
    utils::read.delim(
      tsv_path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) e
  )
  valid <- base::is.data.frame(roundtrip) &&
    base::nrow(roundtrip) == base::nrow(genes_df) &&
    base::ncol(roundtrip) == base::ncol(genes_df) &&
    base::identical(base::as.character(roundtrip$locus_tag), base::as.character(genes_df$locus_tag))
  if (!isTRUE(valid)) {
    detail <- if (base::inherits(roundtrip, "error")) {
      base::conditionMessage(roundtrip)
    } else {
      base::sprintf(
        "expected %d rows x %d columns, read %d rows x %d columns",
        base::nrow(genes_df), base::ncol(genes_df),
        base::nrow(roundtrip), base::ncol(roundtrip)
      )
    }
    base::stop("REBASEfinder input round-trip validation failed: ", detail, call. = FALSE)
  }
  tsv_path
}

.dnmb_rebasefinder_empty_normalized_hits <- function() {
  out <- .dnmb_module_empty_optional_long_table()
  out$blast_identity <- numeric()
  out$blast_evalue <- numeric()
  out$blast_bitscore <- numeric()
  out$blast_length <- numeric()
  out$rec_seq <- character()
  out$operon_id <- character()
  out
}

.dnmb_rebasefinder_normalize_hits <- function(rm_comprehensive, id_col = "locus_tag") {
  if (base::is.null(rm_comprehensive) || !base::is.data.frame(rm_comprehensive) || !base::nrow(rm_comprehensive)) {
    return(.dnmb_rebasefinder_empty_normalized_hits())
  }

  tbl <- base::as.data.frame(rm_comprehensive, stringsAsFactors = FALSE)

  if (!id_col %in% base::names(tbl)) {
    return(.dnmb_rebasefinder_empty_normalized_hits())
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

  pipeline_reference_rec <- base::as.character(tbl[["rec_seq"]])
  has_pipeline_reference <- .dnmb_rebasefinder_valid_recognition(pipeline_reference_rec)
  tbl$reference_rec_seq <- base::ifelse(has_pipeline_reference, pipeline_reference_rec, NA_character_)
  tbl$recognition_match <- base::ifelse(
    has_pipeline_reference,
    .dnmb_rebasefinder_clean_rebase_subject(tbl$hit_label),
    NA_character_
  )
  tbl$recognition_match_subject_id <- base::ifelse(
    has_pipeline_reference,
    base::as.character(tbl$hit_label),
    NA_character_
  )
  tbl$recognition_donor <- tbl$recognition_match
  tbl$recognition_source <- base::ifelse(
    has_pipeline_reference,
    "defenseviz_rebase_match_reference",
    NA_character_
  )
  # DefenseViz recognition sites come from the matched REBASE reference, not
  # from direct evidence for the query strain.
  tbl$rec_seq <- base::rep(NA_character_, base::nrow(tbl))
  tbl$substrate_label <- base::rep(NA_character_, base::nrow(tbl))

  support_parts <- base::vapply(base::seq_len(base::nrow(tbl)), function(i) {
    parts <- c(
      if (!base::is.na(tbl[["rm_type"]][[i]]) && base::nzchar(tbl[["rm_type"]][[i]])) base::paste0("rm_type=", tbl[["rm_type"]][[i]]) else NULL,
      if (!base::is.na(tbl[["subunit"]][[i]]) && base::nzchar(tbl[["subunit"]][[i]])) base::paste0("subunit=", tbl[["subunit"]][[i]]) else NULL,
      if (!base::is.na(tbl[["blast_match"]][[i]]) && base::nzchar(tbl[["blast_match"]][[i]])) base::paste0("rebase_match=", tbl[["blast_match"]][[i]]) else NULL,
      if (!base::is.na(tbl[["blast_identity"]][[i]])) base::paste0("identity=", base::round(tbl[["blast_identity"]][[i]] * 100, 1), "%") else NULL,
      if (!base::is.na(tbl[["reference_rec_seq"]][[i]]) && base::nzchar(tbl[["reference_rec_seq"]][[i]])) base::paste0("reference_rec_seq=", tbl[["reference_rec_seq"]][[i]]) else NULL,
      if (!base::is.na(tbl[["operon_id"]][[i]]) && base::nzchar(tbl[["operon_id"]][[i]])) base::paste0("operon=", tbl[["operon_id"]][[i]]) else NULL
    )
    if (!base::length(parts)) return(NA_character_)
    base::paste(parts, collapse = "; ")
  }, character(1))
  tbl$support <- support_parts

  tbl$typing_eligible <- !base::is.na(tbl[["blast_identity"]]) & tbl[["blast_identity"]] >= 0.5

  extra_cols <- c("blast_identity", "blast_evalue", "blast_bitscore", "blast_length",
                  "rec_seq", "reference_rec_seq", "recognition_match",
                  "recognition_match_subject_id", "recognition_donor",
                  "recognition_source", "operon_id", "partner_locus_tag")
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
                                         homology_modeling = TRUE,
                                         homology_max_candidates = 24L,
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

  if (base::isTRUE(compare_rebase) && base::isTRUE(install)) {
    reference_refresh <- .dnmb_rebasefinder_refresh_reference(cache_root = cache_root, verbose = verbose)
    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
      "rebase_reference",
      if (base::isTRUE(reference_refresh$refreshed)) "updated" else if (base::isTRUE(reference_refresh$ok)) "cached" else "failed",
      reference_refresh$detail
    ))
  }

  input_path <- .dnmb_rebasefinder_prepare_input(genes, output_dir)

  # Patch detect_methyltransferases: deduplicate after union to prevent
  # row explosion (keyword fallback can produce duplicates via left_join)
  detect_patch_original <- NULL
  base::tryCatch({
    current_detect <- DefenseViz:::detect_methyltransferases
    orig_detect <- base::attr(current_detect, "dnmb_original_detect", exact = TRUE)
    if (base::is.null(orig_detect)) orig_detect <- current_detect
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
    base::attr(patched_detect, "dnmb_original_detect") <- orig_detect
    utils::assignInNamespace("detect_methyltransferases", patched_detect, ns = "DefenseViz")
    detect_patch_original <- orig_detect
    if (isTRUE(verbose)) message("[REBASEfinder] Patched detect_methyltransferases (dedup + cap)")
  }, error = function(e) {
    if (isTRUE(verbose)) message("[REBASEfinder] Could not patch detect_methyltransferases: ", conditionMessage(e))
  })
  if (!base::is.null(detect_patch_original)) {
    base::on.exit(
      tryCatch(
        utils::assignInNamespace("detect_methyltransferases", detect_patch_original, ns = "DefenseViz"),
        error = function(e) NULL
      ),
      add = TRUE
    )
  }

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
	  metadata_blast <- if (base::is.data.frame(pipeline_result$blast_filtered)) {
	    base::as.data.frame(pipeline_result$blast_filtered, stringsAsFactors = FALSE)
	  } else {
	    base::data.frame()
	  }
	  rebase_data <- .dnmb_rebasefinder_cached_rebase_data(cache_root)
	  # Context expansion must be seeded by coverage-aware homology, a matching
	  # RM annotation, or structure support. Identity-only DefenseViz calls are
	  # retained as raw evidence but cannot seed neighboring R/M/S candidates.
	  hits <- .dnmb_rebasefinder_select_primary_blast(
	    hits,
	    genes,
	    blast_tbl = metadata_blast,
	    rebase_data = rebase_data,
	    min_identity = blast_min_identity,
	    min_length = blast_min_length,
	    max_evalue = 0.001
	  )
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
  hits <- .dnmb_rebasefinder_add_typeii_context(
    hits,
    genes,
    max_operon_gap = max_operon_gap,
    max_intervening = max_intervening,
    max_neighbors = max(3L, as.integer(typeiii_context_neighbors))
  )
  hits <- .dnmb_rebasefinder_add_typeiv_candidates(hits, genes)
  hits <- .dnmb_rebasefinder_add_motif_evidence(hits, genes)
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

	  supplemental_raw <- base::file.path(
	    output_dir,
	    "supplemental_rebase_blast",
	    "supplemental_blast_results.txt"
	  )
	  if (base::file.exists(supplemental_raw) && base::file.info(supplemental_raw)$size > 0) {
	    parse_blast_results <- base::get("parse_blast_results", envir = base::asNamespace("DefenseViz"))
	    supplemental_tbl <- tryCatch(
	      parse_blast_results(supplemental_raw, verbose = FALSE),
	      error = function(e) base::data.frame()
	    )
	    if (base::is.data.frame(supplemental_tbl) && base::nrow(supplemental_tbl)) {
	      metadata_blast <- dplyr::bind_rows(metadata_blast, supplemental_tbl)
	    }
	  }
	  hits <- .dnmb_rebasefinder_select_primary_blast(
	    hits,
	    genes,
	    blast_tbl = metadata_blast,
	    rebase_data = rebase_data,
	    min_identity = blast_min_identity,
	    min_length = blast_min_length,
	    max_evalue = 0.001
	  )
	  # Supplemental full-length hits can validate a context-only partner. Refresh
	  # operon membership after the final coverage-aware BLAST selection.
	  hits <- .dnmb_rebasefinder_add_typeiii_context(
	    hits, genes,
	    max_operon_gap = max_operon_gap,
	    max_intervening = max_intervening,
	    max_neighbors = typeiii_context_neighbors
	  )
	  hits <- .dnmb_rebasefinder_add_typei_context(
	    hits, genes,
	    max_operon_gap = max_operon_gap,
	    max_intervening = max_intervening,
	    max_neighbors = max(4L, as.integer(typeiii_context_neighbors))
	  )
		  hits <- .dnmb_rebasefinder_add_typeii_context(
		    hits, genes,
		    max_operon_gap = max_operon_gap,
		    max_intervening = max_intervening,
		    max_neighbors = max(3L, as.integer(typeiii_context_neighbors))
		  )
		  hits <- .dnmb_rebasefinder_add_split_r_evidence(hits, genes)
		  hits <- .dnmb_rebasefinder_add_motif_evidence(hits, genes)
	  hits <- .dnmb_rebasefinder_enrich_blast_metadata(
	    hits,
	    blast_tbl = metadata_blast,
	    rebase_data = rebase_data,
	    cache_root = cache_root,
	    query_genbank = genbank,
	    min_identity = blast_min_identity,
	    min_length = blast_min_length,
	    max_evalue = 0.001,
	    verbose = verbose
	  )
	  n_reference_rec <- if ("reference_rec_seq" %in% base::names(hits)) {
	    base::sum(.dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq))
	  } else 0L
	  n_pacbio <- if ("recognition_match_pacbio_available" %in% base::names(hits)) {
	    base::sum(hits$recognition_match_pacbio_available %in% TRUE)
	  } else 0L
	  query_pacbio <- if ("query_strain_pacbio_available" %in% base::names(hits)) {
	    base::any(hits$query_strain_pacbio_available %in% TRUE)
	  } else FALSE
	  status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
	    "rebase_metadata", "ok",
	    base::sprintf(
	      "%d candidates have a REBASE reference recognition sequence; %d reference matches and query strain=%s link to PacBio methylome data",
	      n_reference_rec, n_pacbio, if (query_pacbio) "yes" else "no"
	    )
	  ))
		  hits <- .dnmb_rebasefinder_add_sequence_partial_status(hits, genes)
		  hits <- .dnmb_rebasefinder_curate_hits(hits, genes)
		  augmented_files <- list()
		  if (base::isTRUE(homology_modeling)) {
		    homology <- .dnmb_rebasefinder_run_homology_models(
		      hits = hits,
		      genes = genes,
		      output_dir = output_dir,
		      cache_root = cache_root,
		      install = install,
		      cpu = cpu,
		      max_candidates = homology_max_candidates,
		      verbose = verbose
		    )
		    hits <- homology$hits
		    hits <- .dnmb_rebasefinder_curate_hits(hits, genes)
		    if (base::length(homology$files)) {
		      augmented_files$homology_templates_tsv <- homology$files$tsv %||% NA_character_
		      augmented_files$homology_templates_xlsx <- homology$files$xlsx %||% NA_character_
		      augmented_files$promod3_model_manifest <- homology$files$model_manifest %||% NA_character_
		    }
		    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
		      "promod3_homology_modeling", homology$status, homology$detail
		    ))
		  } else {
		    .dnmb_rebasefinder_clear_promod3_outputs(output_dir)
		    status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
		      "promod3_homology_modeling", "disabled", "Local homology modeling disabled"
		    ))
		  }
		  structure_hits <- hits[!hits$curation_tier %in% c("excluded_noise", "other_defense"), , drop = FALSE]
	  augmented_files$structure_queries_faa <- .dnmb_rebasefinder_write_structure_queries(output_dir, genes, structure_hits)
	  structure_coverage <- .dnmb_rebasefinder_write_structure_coverage(
	    output_dir,
	    genes,
	    structure_hits,
    fasta_path = augmented_files$structure_queries_faa,
	    structure_validation_path = structure_path,
	    structure_tbl = structure_tbl
	  )
	  if (!base::is.null(structure_coverage)) {
	    hits <- .dnmb_rebasefinder_apply_structure_coverage(hits, structure_coverage$table)
	    hits <- .dnmb_rebasefinder_curate_hits(hits, genes)
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
	  curation_workbook <- .dnmb_rebasefinder_write_curation_workbook(
	    output_dir,
	    hits,
	    genes,
	    raw_defenseviz = rm_comprehensive
	  )
	  if (!base::is.null(curation_workbook)) {
	    augmented_files$curated_analysis_xlsx <- curation_workbook$xlsx
	    augmented_files$raw_defenseviz_xlsx <- curation_workbook$raw_xlsx
	  }
	  tier_counts <- base::table(base::factor(
	    hits$curation_tier,
	    levels = c("high", "medium", "review", "other_defense", "excluded_noise")
	  ))
	  status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
	    "evidence_curation", "ok",
	    base::sprintf(
	      "Ranked candidates: %d high, %d medium, %d review, %d other-defense, %d excluded noise",
	      tier_counts[["high"]], tier_counts[["medium"]],
	      tier_counts[["review"]], tier_counts[["other_defense"]], tier_counts[["excluded_noise"]]
	    )
	  ))
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
  all_hits <- hits
  hits <- all_hits[all_hits$curation_keep %in% TRUE, , drop = FALSE]
  output_table <- .dnmb_rebasefinder_output_table(genes = genes, hits = hits)

  n_hits <- base::nrow(hits)
  n_high <- base::sum(hits$curation_tier == "high", na.rm = TRUE)
  n_medium <- base::sum(hits$curation_tier == "medium", na.rm = TRUE)
  n_typeiii <- base::sum(hits$family_id == "Type III", na.rm = TRUE)
  n_typeiii_context <- if ("typeiii_operon_supported" %in% base::names(hits)) {
    base::sum(hits$typeiii_operon_supported %in% TRUE, na.rm = TRUE)
  } else 0L
  n_typei <- base::sum(hits$family_id == "Type I", na.rm = TRUE)
  n_typei_context <- if ("typei_operon_supported" %in% base::names(hits)) {
    base::sum(hits$typei_operon_supported %in% TRUE, na.rm = TRUE)
  } else 0L
  n_typeii <- base::sum(hits$family_id == "Type II", na.rm = TRUE)
  n_typeii_context <- if ("typeii_operon_supported" %in% base::names(hits)) {
    base::sum(hits$typeii_operon_supported %in% TRUE, na.rm = TRUE)
  } else 0L
  n_typeiv <- base::sum(hits$family_id == "Type IV", na.rm = TRUE)

  status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
    "rebasefinder_run", "ok",
    base::sprintf(
      "%d curated R-M candidates (%d high, %d medium; %d Type I/%d operon-supported; %d Type II/%d operon-supported; %d Type III/%d operon-supported; %d Type IV)",
      n_hits, n_high, n_medium, n_typei, n_typei_context, n_typeii, n_typeii_context,
      n_typeiii, n_typeiii_context, n_typeiv
    )
  ))

  list(
    ok = TRUE,
    status = status,
    files = augmented_files,
    hits = hits,
    all_hits = all_hits,
    output_table = output_table,
    pipeline_result = pipeline_result
  )
}
