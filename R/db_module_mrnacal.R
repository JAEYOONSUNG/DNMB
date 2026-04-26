.dnmb_mrnacal_module_name <- function() {
  "mrnacal"
}

.dnmb_mrnacal_default_version <- function() {
  "current"
}

.dnmb_mrnacal_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_mrnacal_empty_status <- function() {
  tibble::tibble(
    component = character(),
    status = character(),
    detail = character()
  )
}

.dnmb_mrnacal_normalize_dna <- function(x) {
  x <- base::toupper(base::as.character(x)[1])
  x <- base::gsub("[^ACGTUN]", "", x)
  x <- base::gsub("U", "T", x, fixed = TRUE)
  x
}

.dnmb_mrnacal_dna_to_rna <- function(x) {
  base::gsub("T", "U", .dnmb_mrnacal_normalize_dna(x), fixed = TRUE)
}

.dnmb_mrnacal_revcomp <- function(x) {
  x <- .dnmb_mrnacal_normalize_dna(x)
  if (!base::nzchar(x)) {
    return("")
  }
  base::as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
}

.dnmb_mrnacal_segment <- function(seq, start, end, circular = FALSE) {
  seq <- .dnmb_mrnacal_normalize_dna(seq)
  n <- base::nchar(seq)
  start <- base::as.integer(start)
  end <- base::as.integer(end)
  if (!base::nzchar(seq) || base::is.na(start) || base::is.na(end) || end < start) {
    return("")
  }
  if (base::isTRUE(circular) && n > 0L) {
    pos <- ((base::seq.int(start, end) - 1L) %% n) + 1L
    return(base::paste0(base::substring(seq, pos, pos), collapse = ""))
  }
  left <- base::max(1L, start)
  right <- base::min(n, end)
  if (right < left) {
    return("")
  }
  base::substr(seq, left, right)
}

.dnmb_mrnacal_circular_records <- function(genbank = NULL) {
  if (base::is.null(genbank) || !base::file.exists(genbank)) {
    return(logical())
  }
  lines <- base::readLines(genbank, warn = FALSE)
  loci <- base::grep("^LOCUS", lines, value = TRUE)
  if (!base::length(loci)) {
    return(logical())
  }
  base::grepl("\\bcircular\\b", loci, ignore.case = TRUE)
}

.dnmb_mrnacal_contig_table <- function(genes, genbank = NULL) {
  if (!base::is.data.frame(genes)) {
    genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  }
  rows <- list()

  if (!base::is.null(genbank) && base::file.exists(genbank)) {
    records <- .dnmb_prophage_parse_genbank_records(genbank)
    circular <- .dnmb_mrnacal_circular_records(genbank)
    if (base::nrow(records)) {
      for (i in base::seq_len(base::nrow(records))) {
        rows[[base::length(rows) + 1L]] <- data.frame(
          contig_number = i,
          contig = records$definition[[i]] %||% records$accession[[i]] %||% records$locus[[i]],
          accession = records$accession[[i]] %||% NA_character_,
          locus = records$locus[[i]] %||% NA_character_,
          sequence = .dnmb_mrnacal_normalize_dna(records$sequence[[i]]),
          circular = if (base::length(circular) >= i) base::isTRUE(circular[[i]]) else FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!base::length(rows)) {
    contig_names <- base::ls(envir = .GlobalEnv, pattern = "^contig_[0-9]+_seq$")
    if (base::length(contig_names)) {
      contig_numbers <- base::as.integer(base::sub("^contig_([0-9]+)_seq$", "\\1", contig_names))
      contig_names <- contig_names[base::order(contig_numbers)]
      contig_numbers <- base::sort(contig_numbers)
      unique_contigs <- if ("contig" %in% base::names(genes)) base::unique(base::as.character(genes$contig)) else character()
      for (i in base::seq_along(contig_names)) {
        rows[[base::length(rows) + 1L]] <- data.frame(
          contig_number = contig_numbers[[i]],
          contig = if (base::length(unique_contigs) >= i) unique_contigs[[i]] else contig_names[[i]],
          accession = NA_character_,
          locus = contig_names[[i]],
          sequence = .dnmb_mrnacal_normalize_dna(base::get(contig_names[[i]], envir = .GlobalEnv)),
          circular = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!base::length(rows)) {
    return(data.frame())
  }
  out <- dplyr::bind_rows(rows)
  out <- out[!base::is.na(out$sequence) & base::nzchar(out$sequence), , drop = FALSE]
  out$length_bp <- base::nchar(out$sequence)
  rownames(out) <- NULL
  out
}

.dnmb_mrnacal_match_contig <- function(row, contigs) {
  if (!base::nrow(contigs)) {
    return(NULL)
  }
  if ("contig_number" %in% base::names(row)) {
    idx <- base::match(base::as.integer(row$contig_number[[1]]), base::as.integer(contigs$contig_number))
    if (!base::is.na(idx)) {
      return(contigs[idx, , drop = FALSE])
    }
  }
  contig_value <- if ("contig" %in% base::names(row)) base::as.character(row$contig[[1]]) else NA_character_
  if (!base::is.na(contig_value) && base::nzchar(contig_value)) {
    keys <- c("contig", "accession", "locus")
    for (key in keys[keys %in% base::names(contigs)]) {
      idx <- base::match(contig_value, base::as.character(contigs[[key]]))
      if (!base::is.na(idx)) {
        return(contigs[idx, , drop = FALSE])
      }
    }
    idx <- base::which(base::startsWith(base::as.character(contigs$contig), contig_value) |
                         base::startsWith(contig_value, base::as.character(contigs$contig)))
    if (base::length(idx)) {
      return(contigs[idx[[1]], , drop = FALSE])
    }
  }
  NULL
}

.dnmb_mrnacal_transcript_window <- function(row, contigs, upstream = 60L, downstream = 60L) {
  ctg <- .dnmb_mrnacal_match_contig(row, contigs)
  if (base::is.null(ctg)) {
    return(NULL)
  }
  seq <- ctg$sequence[[1]]
  circular <- base::isTRUE(ctg$circular[[1]])
  start <- base::as.integer(row$start[[1]])
  end <- base::as.integer(row$end[[1]])
  direction <- base::as.character(row$direction[[1]])
  if (base::is.na(start) || base::is.na(end) || end < start || !direction %in% c("+", "-")) {
    return(NULL)
  }

  upstream <- base::as.integer(upstream)[1]
  downstream <- base::as.integer(downstream)[1]
  if (base::is.na(upstream) || upstream < 0L) upstream <- 60L
  if (base::is.na(downstream) || downstream < 3L) downstream <- 60L

  if (identical(direction, "+")) {
    up_seq <- .dnmb_mrnacal_segment(seq, start - upstream, start - 1L, circular = circular)
    cds_seq <- .dnmb_mrnacal_segment(seq, start, start + downstream - 1L, circular = circular)
    anchor <- start
  } else {
    up_genomic <- .dnmb_mrnacal_segment(seq, end + 1L, end + upstream, circular = circular)
    cds_genomic <- .dnmb_mrnacal_segment(seq, end - downstream + 1L, end, circular = circular)
    up_seq <- .dnmb_mrnacal_revcomp(up_genomic)
    cds_seq <- .dnmb_mrnacal_revcomp(cds_genomic)
    anchor <- end
  }

  combined_dna <- base::paste0(up_seq, cds_seq)
  list(
    sequence_dna = combined_dna,
    sequence_rna = .dnmb_mrnacal_dna_to_rna(combined_dna),
    upstream_len_observed = base::nchar(up_seq),
    downstream_len_observed = base::nchar(cds_seq),
    start_codon = base::substr(cds_seq, 1L, 3L),
    anchor = anchor,
    contig_number = ctg$contig_number[[1]],
    contig = ctg$contig[[1]],
    circular = circular
  )
}

.dnmb_mrnacal_all_codons <- function() {
  bases <- c("A", "C", "G", "T")
  codons <- character()
  for (b1 in bases) for (b2 in bases) for (b3 in bases) {
    codons <- c(codons, base::paste0(b1, b2, b3))
  }
  codons
}

.dnmb_mrnacal_genetic_code <- function() {
  c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )
}

.dnmb_mrnacal_codon_counts <- function(cds_dna) {
  cds_dna <- .dnmb_mrnacal_normalize_dna(cds_dna)
  counts <- stats::setNames(base::integer(64L), .dnmb_mrnacal_all_codons())
  if (!base::nzchar(cds_dna) || base::nchar(cds_dna) < 3L) {
    return(counts)
  }
  n_codons <- base::floor(base::nchar(cds_dna) / 3L)
  if (n_codons < 1L) {
    return(counts)
  }
  starts <- base::seq.int(1L, n_codons * 3L, by = 3L)
  ends <- starts + 2L
  codons <- base::substring(cds_dna, starts, ends)
  codons <- codons[base::nchar(codons) == 3L & !base::grepl("[^ACGT]", codons)]
  if (!base::length(codons)) {
    return(counts)
  }
  tab <- base::table(codons)
  for (cdn in base::names(tab)) {
    if (cdn %in% base::names(counts)) {
      counts[[cdn]] <- counts[[cdn]] + base::as.integer(tab[[cdn]])
    }
  }
  counts
}

.dnmb_mrnacal_codon_counts_total <- function(cds_list) {
  total <- stats::setNames(base::integer(64L), .dnmb_mrnacal_all_codons())
  for (cds in cds_list) {
    if (base::is.na(cds) || !base::nzchar(cds)) next
    total <- total + .dnmb_mrnacal_codon_counts(cds)
  }
  total
}

.dnmb_mrnacal_cai_weights <- function(counts) {
  code <- .dnmb_mrnacal_genetic_code()
  codons <- .dnmb_mrnacal_all_codons()
  weights <- stats::setNames(base::rep(NA_real_, 64L), codons)
  excluded <- c("M", "W", "*")
  for (aa in base::unique(code)) {
    family_codons <- codons[code == aa]
    if (aa %in% excluded || base::length(family_codons) < 2L) {
      next
    }
    family_counts <- counts[family_codons]
    family_max <- base::max(family_counts, na.rm = TRUE)
    if (!base::is.finite(family_max) || family_max <= 0L) {
      next
    }
    w <- base::pmax(family_counts, 0.5) / family_max
    weights[family_codons] <- base::pmin(1, w)
  }
  weights
}

.dnmb_mrnacal_cai <- function(cds_dna, weights) {
  if (!base::length(weights) || base::all(base::is.na(weights))) {
    return(NA_real_)
  }
  counts <- .dnmb_mrnacal_codon_counts(cds_dna)
  inform <- !base::is.na(weights) & counts > 0L
  if (!base::any(inform)) {
    return(NA_real_)
  }
  w <- weights[inform]
  n <- counts[inform]
  base::exp(base::sum(n * base::log(w)) / base::sum(n))
}

.dnmb_mrnacal_parse_anticodon_dna <- function(s) {
  if (base::is.na(s) || !base::nzchar(s)) {
    return(NA_character_)
  }
  m <- stringr::str_match(base::as.character(s)[1], "seq:([acgtuACGTU]{3})")[, 2]
  if (base::is.na(m)) {
    return(NA_character_)
  }
  base::toupper(base::gsub("U", "T", m, fixed = TRUE))
}

.dnmb_mrnacal_aa3_to_aa1 <- function(aa3) {
  table <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Glu = "E", Gln = "Q", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Sec = "U", Pyl = "O"
  )
  table[base::as.character(aa3)]
}

.dnmb_mrnacal_infer_anticodon_from_seq <- function(row) {
  if (!"rearranged_nt_seq" %in% base::names(row) || !"product" %in% base::names(row)) {
    return(NA_character_)
  }
  seq_dna <- .dnmb_mrnacal_normalize_dna(row$rearranged_nt_seq[[1]])
  product <- base::as.character(row$product[[1]])
  if (!base::nzchar(seq_dna) || base::nchar(seq_dna) < 70L || base::nchar(seq_dna) > 110L) {
    return(NA_character_)
  }
  match_aa <- stringr::str_match(product, "^tRNA-([A-Z][a-z]{2})")[, 2]
  if (base::is.na(match_aa)) {
    return(NA_character_)
  }
  expected_aa <- base::as.character(.dnmb_mrnacal_aa3_to_aa1(match_aa))
  if (base::is.na(expected_aa)) {
    return(NA_character_)
  }
  code <- .dnmb_mrnacal_genetic_code()
  L <- base::nchar(seq_dna)
  offsets <- base::unique(c(L - 42L, L - 41L, L - 43L, L - 40L, 34L, 35L, 33L))
  offsets <- offsets[offsets >= 1L & (offsets + 2L) <= L]
  for (offset in offsets) {
    ac <- base::toupper(base::substr(seq_dna, offset, offset + 2L))
    if (base::nchar(ac) != 3L || base::grepl("[^ACGT]", ac)) next
    codon <- .dnmb_mrnacal_revcomp(ac)
    aa <- code[codon]
    if (!base::is.na(aa) && base::identical(base::as.character(aa), expected_aa)) {
      return(ac)
    }
  }
  NA_character_
}

.dnmb_mrnacal_resolve_anticodon <- function(row, contigs = NULL) {
  s <- if ("anticodon" %in% base::names(row)) base::as.character(row$anticodon[[1]]) else NA_character_
  if (!base::is.na(s) && base::nzchar(s)) {
    ac <- .dnmb_mrnacal_parse_anticodon_dna(s)
    if (!base::is.na(ac)) {
      return(ac)
    }
    if (!base::is.null(contigs) && base::is.data.frame(contigs) && base::nrow(contigs)) {
      is_complement <- base::grepl("complement\\(", s)
      pos <- stringr::str_match(s, "(\\d+)\\.\\.(\\d+)")
      ac_start <- suppressWarnings(base::as.integer(pos[, 2]))
      ac_end <- suppressWarnings(base::as.integer(pos[, 3]))
      if (!base::is.na(ac_start) && !base::is.na(ac_end) && ac_end >= ac_start) {
        ctg <- .dnmb_mrnacal_match_contig(row, contigs)
        if (!base::is.null(ctg)) {
          seq_str <- ctg$sequence[[1]]
          ac_seg <- .dnmb_mrnacal_segment(seq_str, ac_start, ac_end, circular = base::isTRUE(ctg$circular[[1]]))
          if (base::nzchar(ac_seg) && base::nchar(ac_seg) == 3L) {
            if (is_complement) {
              ac_seg <- .dnmb_mrnacal_revcomp(ac_seg)
            }
            return(base::toupper(ac_seg))
          }
        }
      }
    }
  }
  .dnmb_mrnacal_infer_anticodon_from_seq(row)
}

.dnmb_mrnacal_trna_gcn <- function(genes, contigs = NULL) {
  empty <- stats::setNames(base::integer(64L), .dnmb_mrnacal_all_codons())
  if (!base::is.data.frame(genes) || !base::nrow(genes)) {
    return(empty)
  }
  has_anticodon <- "anticodon" %in% base::names(genes)
  has_product <- "product" %in% base::names(genes)
  has_seq <- "rearranged_nt_seq" %in% base::names(genes)
  if (!has_anticodon && !(has_product && has_seq)) {
    return(empty)
  }
  rows <- genes
  if (has_product) {
    is_trna <- !base::is.na(rows$product) & base::grepl("^tRNA-[A-Z][a-z]{2}", rows$product)
    rows <- rows[is_trna, , drop = FALSE]
  } else if (has_anticodon) {
    rows <- rows[!base::is.na(rows$anticodon) & base::nzchar(rows$anticodon), , drop = FALSE]
  }
  if (!base::nrow(rows)) {
    return(empty)
  }
  acs <- base::vapply(base::seq_len(base::nrow(rows)), function(i) {
    .dnmb_mrnacal_resolve_anticodon(rows[i, , drop = FALSE], contigs = contigs)
  }, character(1))
  acs <- acs[!base::is.na(acs)]
  if (!base::length(acs)) {
    return(empty)
  }
  codon_for_ac <- base::vapply(acs, .dnmb_mrnacal_revcomp, character(1))
  gcn_by_codon <- empty
  tab <- base::table(codon_for_ac)
  for (cdn in base::names(tab)) {
    if (cdn %in% base::names(gcn_by_codon)) {
      gcn_by_codon[[cdn]] <- gcn_by_codon[[cdn]] + base::as.integer(tab[[cdn]])
    }
  }
  gcn_by_codon
}

.dnmb_mrnacal_tai_weights <- function(gcn_by_codon) {
  codons <- .dnmb_mrnacal_all_codons()
  W <- stats::setNames(base::rep(0, 64L), codons)
  s_GU <- 0.41
  s_UG <- 0.68
  for (canonical in base::names(gcn_by_codon)) {
    g <- gcn_by_codon[[canonical]]
    if (base::is.na(g) || g <= 0L) next
    W[[canonical]] <- W[[canonical]] + g
    nt3 <- base::substr(canonical, 3L, 3L)
    prefix <- base::substr(canonical, 1L, 2L)
    anticodon_5p <- switch(nt3, A = "T", C = "G", G = "C", T = "A", NA_character_)
    if (base::is.na(anticodon_5p)) next
    if (anticodon_5p == "G") {
      partner <- base::paste0(prefix, "T")
      if (partner %in% codons) {
        W[[partner]] <- W[[partner]] + (1 - s_GU) * g
      }
    } else if (anticodon_5p == "T") {
      partner <- base::paste0(prefix, "G")
      if (partner %in% codons) {
        W[[partner]] <- W[[partner]] + (1 - s_UG) * g
      }
    }
  }
  code <- .dnmb_mrnacal_genetic_code()
  excluded <- c("M", "W", "*")
  inform <- !(code %in% excluded)
  W_inform <- W[inform]
  W_max <- base::max(W_inform, na.rm = TRUE)
  if (!base::is.finite(W_max) || W_max <= 0) {
    return(stats::setNames(base::rep(NA_real_, 64L), codons))
  }
  w <- W / W_max
  zero_w <- inform & (W <= 0)
  if (base::any(zero_w)) {
    pos <- inform & (W > 0)
    if (base::any(pos)) {
      gm <- base::exp(base::mean(base::log(w[pos])))
      w[zero_w] <- gm
    } else {
      w[zero_w] <- NA_real_
    }
  }
  w[!inform] <- NA_real_
  w
}

.dnmb_mrnacal_tai <- function(cds_dna, weights) {
  if (!base::length(weights) || base::all(base::is.na(weights))) {
    return(NA_real_)
  }
  counts <- .dnmb_mrnacal_codon_counts(cds_dna)
  inform <- !base::is.na(weights) & counts > 0L
  if (!base::any(inform)) {
    return(NA_real_)
  }
  w <- weights[inform]
  n <- counts[inform]
  w_safe <- base::pmax(w, 1e-6)
  base::exp(base::sum(n * base::log(w_safe)) / base::sum(n))
}

.dnmb_mrnacal_full_cds <- function(row, contigs) {
  if ("rearranged_nt_seq" %in% base::names(row)) {
    s <- .dnmb_mrnacal_normalize_dna(row$rearranged_nt_seq[[1]])
    if (base::nzchar(s) && base::nchar(s) >= 3L) {
      return(s)
    }
  }
  ctg <- .dnmb_mrnacal_match_contig(row, contigs)
  if (base::is.null(ctg)) {
    return("")
  }
  seq <- ctg$sequence[[1]]
  circular <- base::isTRUE(ctg$circular[[1]])
  start <- base::as.integer(row$start[[1]])
  end <- base::as.integer(row$end[[1]])
  direction <- base::as.character(row$direction[[1]])
  if (base::is.na(start) || base::is.na(end) || end < start || !direction %in% c("+", "-")) {
    return("")
  }
  cds <- .dnmb_mrnacal_segment(seq, start, end, circular = circular)
  if (direction == "-") {
    cds <- .dnmb_mrnacal_revcomp(cds)
  }
  cds
}

.dnmb_mrnacal_select_cai_reference <- function(genes) {
  if (!base::is.data.frame(genes) || !base::nrow(genes)) {
    return(list(mask = base::logical(), label = "all_cds"))
  }
  if (!"product" %in% base::names(genes)) {
    return(list(mask = base::rep(TRUE, base::nrow(genes)), label = "all_cds"))
  }
  product <- base::as.character(genes$product)
  is_ribo <- !base::is.na(product) &
    base::grepl("^([35]0S\\s+)?ribosomal protein\\s+[SL][0-9]+", product, ignore.case = TRUE)
  if (base::sum(is_ribo) >= 8L) {
    return(list(mask = is_ribo, label = "ribosomal_proteins"))
  }
  list(mask = !base::is.na(product) & base::nzchar(product), label = "all_cds")
}

.dnmb_mrnacal_sd_seed <- function(genes, translation_domain = "bacteria", sd_seed = NULL) {
  if (!base::is.null(sd_seed) && base::nzchar(base::as.character(sd_seed)[1])) {
    return(.dnmb_mrnacal_normalize_dna(sd_seed))
  }
  if (all(c("product", "rearranged_nt_seq") %in% base::names(genes))) {
    rrna_cols <- c("product", "rearranged_nt_seq")
    rrna_src <- genes[, rrna_cols, drop = FALSE]
    rrna <- rrna_src[.dnmb_rrna_product_mask(rrna_src$product), , drop = FALSE]
    rrna <- rrna[!base::is.na(rrna$rearranged_nt_seq) & base::nzchar(base::as.character(rrna$rearranged_nt_seq)), , drop = FALSE]
    if (base::nrow(rrna)) {
      seeds <- base::vapply(rrna$rearranged_nt_seq, function(x) {
        x <- .dnmb_mrnacal_normalize_dna(x)
        tail <- base::substr(x, base::max(1L, base::nchar(x) - 11L), base::nchar(x))
        .dnmb_mrnacal_revcomp(tail)
      }, character(1))
      seeds <- seeds[!base::is.na(seeds) & base::nzchar(seeds)]
      if (base::length(seeds)) {
        return(names(sort(table(seeds), decreasing = TRUE))[[1]])
      }
    }
  }
  if (identical(translation_domain, "archaea")) {
    return("GGAGG")
  }
  "AGGAGG"
}

.dnmb_mrnacal_seed_kmers <- function(seed, min_len = 4L, max_len = 9L) {
  seed <- .dnmb_mrnacal_normalize_dna(seed)
  n <- base::nchar(seed)
  if (!base::nzchar(seed) || n < min_len) {
    return(c("AGGAGG", "GGAGG", "AGGA"))
  }
  kmers <- character()
  for (k in base::seq.int(min_len, base::min(max_len, n))) {
    for (i in base::seq_len(n - k + 1L)) {
      kmers <- c(kmers, base::substr(seed, i, i + k - 1L))
    }
  }
  unique(c(seed, kmers, "AGGAGG", "GGAGG", "AGGA"))
}

.dnmb_mrnacal_hamming <- function(a, b) {
  aa <- base::strsplit(a, "", fixed = TRUE)[[1]]
  bb <- base::strsplit(b, "", fixed = TRUE)[[1]]
  base::sum(aa != bb)
}

.dnmb_mrnacal_best_rbs <- function(sequence_dna, upstream_len, sd_seed = "AGGAGG") {
  sequence_dna <- .dnmb_mrnacal_normalize_dna(sequence_dna)
  upstream_len <- base::as.integer(upstream_len)[1]
  if (!base::nzchar(sequence_dna) || base::is.na(upstream_len) || upstream_len < 4L) {
    return(list(
      motif = NA_character_, seed = sd_seed, mismatches = NA_integer_,
      spacer = NA_integer_, rbs_start = NA_integer_, rbs_end = NA_integer_,
      score = 0
    ))
  }

  upstream <- base::substr(sequence_dna, 1L, upstream_len)
  kmers <- .dnmb_mrnacal_seed_kmers(sd_seed)
  best <- NULL
  scan_left <- base::max(1L, upstream_len - 35L)
  scan_right <- base::max(1L, upstream_len - 3L)
  for (motif_seed in kmers) {
    k <- base::nchar(motif_seed)
    if (k > upstream_len || scan_right - k + 1L < scan_left) {
      next
    }
    for (pos in base::seq.int(scan_left, scan_right - k + 1L)) {
      hit <- base::substr(upstream, pos, pos + k - 1L)
      mismatches <- .dnmb_mrnacal_hamming(hit, motif_seed)
      max_mismatch <- base::max(1L, base::floor(k / 3L))
      if (mismatches > max_mismatch) {
        next
      }
      spacer <- upstream_len - (pos + k - 1L)
      if (spacer < 3L || spacer > 15L) {
        next
      }
      motif_quality <- (k - mismatches) / k
      spacing_quality <- base::exp(-((spacer - 7)^2) / (2 * 2.6^2))
      score <- 100 * (0.72 * motif_quality + 0.28 * spacing_quality)
      candidate <- list(
        motif = hit,
        seed = motif_seed,
        mismatches = mismatches,
        spacer = spacer,
        rbs_start = pos,
        rbs_end = pos + k - 1L,
        score = score
      )
      if (base::is.null(best) ||
          candidate$score > best$score ||
          (base::identical(candidate$score, best$score) && base::nchar(candidate$motif) > base::nchar(best$motif))) {
        best <- candidate
      }
    }
  }

  if (base::is.null(best)) {
    return(list(
      motif = NA_character_, seed = sd_seed, mismatches = NA_integer_,
      spacer = NA_integer_, rbs_start = NA_integer_, rbs_end = NA_integer_,
      score = 0
    ))
  }
  best$score <- base::round(best$score, 2)
  best
}

.dnmb_mrnacal_start_score <- function(start_codon) {
  start_codon <- .dnmb_mrnacal_normalize_dna(start_codon)
  score <- switch(
    start_codon,
    ATG = 100,
    GTG = 72,
    TTG = 55,
    CTG = 35,
    10
  )
  base::as.numeric(score)
}

.dnmb_mrnacal_internal_sd <- function(cds_dna, sd_seed = "AGGAGG", scan_len = 60L) {
  cds_dna <- .dnmb_mrnacal_normalize_dna(cds_dna)
  empty <- list(
    internal_sd_count = 0L,
    internal_sd_motifs = NA_character_,
    internal_sd_positions = NA_character_,
    internal_sd_min_position = NA_integer_,
    internal_sd_penalty = 0
  )
  if (!base::nzchar(cds_dna)) {
    return(empty)
  }
  scan_len <- base::as.integer(scan_len)[1]
  if (base::is.na(scan_len) || scan_len < 6L) {
    scan_len <- 60L
  }
  region <- base::substr(cds_dna, 1L, base::min(scan_len, base::nchar(cds_dna)))
  if (base::nchar(region) < 6L) {
    return(empty)
  }
  seed_kmers <- base::unique(.dnmb_mrnacal_seed_kmers(sd_seed, min_len = 6L, max_len = 6L))
  seed_kmers <- seed_kmers[base::nchar(seed_kmers) == 6L]
  if (!base::length(seed_kmers)) {
    seed_kmers <- "AGGAGG"
  }
  n <- base::nchar(region)
  hits <- list()
  used <- base::rep(FALSE, n)
  for (pos in base::seq_len(n - 5L)) {
    if (base::any(used[pos:(pos + 5L)])) {
      next
    }
    window <- base::substr(region, pos, pos + 5L)
    best_seed <- NA_character_
    best_mm <- 7L
    for (seed in seed_kmers) {
      mm <- .dnmb_mrnacal_hamming(window, seed)
      if (mm < best_mm) {
        best_mm <- mm
        best_seed <- seed
      }
    }
    if (best_mm <= 1L) {
      hits[[base::length(hits) + 1L]] <- list(
        position = pos,
        motif = window,
        seed = best_seed,
        mismatches = best_mm
      )
      used[pos:(pos + 5L)] <- TRUE
    }
  }
  if (!base::length(hits)) {
    return(empty)
  }
  count <- base::length(hits)
  motifs <- vapply(hits, function(h) h$motif, character(1))
  positions <- vapply(hits, function(h) h$position, integer(1))
  penalty <- base::min(20, 5 * count + 3 * base::sum(vapply(hits, function(h) h$mismatches == 0L, logical(1))))
  list(
    internal_sd_count = base::as.integer(count),
    internal_sd_motifs = base::paste(motifs, collapse = ","),
    internal_sd_positions = base::paste(positions, collapse = ","),
    internal_sd_min_position = base::as.integer(base::min(positions)),
    internal_sd_penalty = base::round(penalty, 2)
  )
}

.dnmb_mrnacal_early_context <- function(sequence_dna, translation = NA_character_) {
  sequence_dna <- .dnmb_mrnacal_normalize_dna(sequence_dna)
  cds <- if (base::nchar(sequence_dna) >= 3L) sequence_dna else ""
  codon_count <- base::floor(base::nchar(cds) / 3L)
  codons <- if (codon_count > 0L) {
    base::substring(cds, base::seq(1L, codon_count * 3L, by = 3L), base::seq(3L, codon_count * 3L, by = 3L))
  } else {
    character()
  }
  first_codons <- codons[base::seq_len(base::min(10L, base::length(codons)))]
  early <- if (base::length(first_codons) >= 2L) {
    first_codons[base::seq.int(2L, base::min(8L, base::length(first_codons)))]
  } else {
    character()
  }
  ncs45_codons <- if (base::length(codons) >= 2L) {
    codons[base::seq.int(2L, base::min(16L, base::length(codons)))]
  } else {
    character()
  }
  ncs45_sequence <- base::paste(ncs45_codons, collapse = "")
  ncs45_chars <- if (base::nzchar(ncs45_sequence)) base::strsplit(ncs45_sequence, "", fixed = TRUE)[[1]] else character()
  ncs45_at_fraction <- if (base::length(ncs45_chars)) base::mean(ncs45_chars %in% c("A", "T")) else NA_real_
  ncs45_lysine_count <- base::sum(ncs45_codons %in% c("AAA", "AAG"), na.rm = TRUE)
  ncs45_aaa_count <- base::sum(ncs45_codons %in% "AAA", na.rm = TRUE)
  ncs45_aag_count <- base::sum(ncs45_codons %in% "AAG", na.rm = TRUE)
  lys_codons <- early %in% c("AAA", "AAG")
  aaa_codons <- early %in% "AAA"
  aag_codons <- early %in% "AAG"
  second_codon <- if (base::length(first_codons) >= 2L) first_codons[[2]] else NA_character_
  early_nt <- base::paste(early, collapse = "")
  early_nt_chars <- if (base::nzchar(early_nt)) base::strsplit(early_nt, "", fixed = TRUE)[[1]] else character()
  early_at_fraction <- if (base::length(early_nt_chars)) base::mean(early_nt_chars %in% c("A", "T")) else NA_real_
  aaa_run <- if (base::any(aaa_codons)) {
    r <- base::rle(aaa_codons)
    base::max(r$lengths[r$values], na.rm = TRUE)
  } else {
    0L
  }
  early_poly_a_run <- .dnmb_mrnacal_max_at_run(early_nt)

  protein <- base::toupper(base::gsub("[^A-Z*]", "", base::as.character(translation)[1]))
  aa_2_8 <- if (!base::is.na(protein) && base::nchar(protein) >= 2L) {
    base::substr(protein, 2L, base::min(8L, base::nchar(protein)))
  } else {
    ""
  }
  aa_chars <- if (base::nzchar(aa_2_8)) base::strsplit(aa_2_8, "", fixed = TRUE)[[1]] else character()
  basic_count <- base::sum(aa_chars %in% c("K", "R"), na.rm = TRUE)
  lys_count <- base::sum(lys_codons, na.rm = TRUE)
  aaa_count <- base::sum(aaa_codons, na.rm = TRUE)
  aag_count <- base::sum(aag_codons, na.rm = TRUE)
  second_lys <- !base::is.na(second_codon) && second_codon %in% c("AAA", "AAG")
  skik_like <- base::grepl("^S?KI[KRN]", aa_2_8)
  poly_a_penalty <- 12 * base::max(0L, aaa_run - 1L) + 8 * base::max(0L, early_poly_a_run - 8L)
  early_at_for_score <- if (base::is.na(early_at_fraction)) 0.5 else early_at_fraction

  score <- 18 +
    14 * base::min(3L, lys_count) +
    6 * aag_count +
    12 * base::as.integer(second_lys) +
    5 * basic_count +
    12 * base::as.integer(skik_like) +
    18 * early_at_for_score -
    poly_a_penalty
  score <- base::min(100, base::max(0, score))
  list(
    score = base::round(score, 2),
    second_codon = second_codon,
    early_codons = base::paste(first_codons, collapse = " "),
    lysine_codon_count_2_8 = lys_count,
    aaa_count_2_8 = aaa_count,
    aag_count_2_8 = aag_count,
    early_coding_at_fraction = base::round(early_at_fraction, 4),
    ncs45_sequence = if (base::nzchar(ncs45_sequence)) ncs45_sequence else NA_character_,
    ncs45_at_fraction = base::round(ncs45_at_fraction, 4),
    ncs45_lysine_codon_count = ncs45_lysine_count,
    ncs45_aaa_count = ncs45_aaa_count,
    ncs45_aag_count = ncs45_aag_count,
    early_aaa_run = aaa_run,
    early_poly_a_run = early_poly_a_run,
    early_poly_a_penalty = base::round(poly_a_penalty, 2),
    basic_aa_count_2_8 = basic_count,
    skik_like = skik_like,
    aa_2_8 = aa_2_8
  )
}

.dnmb_mrnacal_weighted_mean <- function(values, weights, fallback = NA_real_) {
  values <- suppressWarnings(base::as.numeric(values))
  weights <- suppressWarnings(base::as.numeric(weights))
  keep <- !base::is.na(values) & !base::is.na(weights) & weights > 0
  if (!base::any(keep)) {
    return(fallback)
  }
  base::sum(values[keep] * weights[keep]) / base::sum(weights[keep])
}

.dnmb_mrnacal_max_at_run <- function(sequence_dna) {
  sequence_dna <- .dnmb_mrnacal_normalize_dna(sequence_dna)
  if (!base::nzchar(sequence_dna)) {
    return(0L)
  }
  hits <- base::gregexpr("[AT]+", sequence_dna, perl = TRUE)[[1]]
  if (base::identical(hits[[1]], -1L)) {
    return(0L)
  }
  base::max(base::attr(hits, "match.length"), na.rm = TRUE)
}

.dnmb_mrnacal_upstream_context <- function(sequence_dna, upstream_len, rbs_start = NA_integer_) {
  sequence_dna <- .dnmb_mrnacal_normalize_dna(sequence_dna)
  upstream_len <- suppressWarnings(base::as.integer(upstream_len)[1])
  rbs_start <- suppressWarnings(base::as.integer(rbs_start)[1])
  if (!base::nzchar(sequence_dna) || base::is.na(upstream_len) || upstream_len < 1L) {
    return(list(
      upstream20_sequence = NA_character_, upstream20_at_fraction = NA_real_,
      upstream20_score = NA_real_, enhancer_sequence = NA_character_,
      enhancer_at_fraction = NA_real_, enhancer_at_run = NA_integer_,
      upstream_au_score = NA_real_
    ))
  }

  up_end <- upstream_len
  up_start <- base::max(1L, upstream_len - 19L)
  upstream20 <- base::substr(sequence_dna, up_start, up_end)
  upstream20_chars <- base::strsplit(upstream20, "", fixed = TRUE)[[1]]
  upstream20_at <- base::mean(upstream20_chars %in% c("A", "T"))
  upstream20_score <- 100 * upstream20_at

  if (!base::is.na(rbs_start) && rbs_start > 1L) {
    enhancer_start <- base::max(1L, rbs_start - 12L)
    enhancer_end <- rbs_start - 1L
  } else {
    enhancer_start <- base::max(1L, upstream_len - 35L)
    enhancer_end <- base::max(1L, upstream_len - 19L)
  }
  enhancer <- if (enhancer_end >= enhancer_start) {
    base::substr(sequence_dna, enhancer_start, enhancer_end)
  } else {
    ""
  }
  enhancer_chars <- if (base::nzchar(enhancer)) base::strsplit(enhancer, "", fixed = TRUE)[[1]] else character()
  enhancer_at <- if (base::length(enhancer_chars)) base::mean(enhancer_chars %in% c("A", "T")) else NA_real_
  enhancer_run <- .dnmb_mrnacal_max_at_run(enhancer)
  upstream_au_score <- if (!base::is.na(enhancer_at)) {
    100 * (
      0.58 * enhancer_at +
        0.27 * base::min(1, enhancer_run / 6) +
        0.15 * upstream20_at
    )
  } else {
    upstream20_score
  }
  list(
    upstream20_sequence = upstream20,
    upstream20_at_fraction = base::round(upstream20_at, 4),
    upstream20_score = base::round(upstream20_score, 2),
    enhancer_sequence = if (base::nzchar(enhancer)) enhancer else NA_character_,
    enhancer_at_fraction = base::round(enhancer_at, 4),
    enhancer_at_run = enhancer_run,
    upstream_au_score = base::round(base::max(0, base::min(100, upstream_au_score)), 2)
  )
}

.dnmb_mrnacal_local_tir_window <- function(sequence_dna, upstream_len, upstream_flank = 35L, downstream_flank = 30L) {
  sequence_dna <- .dnmb_mrnacal_normalize_dna(sequence_dna)
  upstream_len <- suppressWarnings(base::as.integer(upstream_len)[1])
  if (!base::nzchar(sequence_dna) || base::is.na(upstream_len)) {
    return(list(sequence_dna = "", sequence_rna = "", offset = 0L, upstream_len = 0L, downstream_len = 0L))
  }
  start_coord <- upstream_len + 1L
  left <- base::max(1L, start_coord - base::as.integer(upstream_flank)[1])
  right <- base::min(base::nchar(sequence_dna), upstream_len + base::as.integer(downstream_flank)[1])
  local_dna <- if (right >= left) base::substr(sequence_dna, left, right) else ""
  offset <- left - 1L
  list(
    sequence_dna = local_dna,
    sequence_rna = .dnmb_mrnacal_dna_to_rna(local_dna),
    offset = offset,
    upstream_len = base::max(0L, upstream_len - offset),
    downstream_len = base::max(0L, right - upstream_len)
  )
}

.dnmb_mrnacal_write_fasta <- function(ids, seqs, path) {
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- base::file(path, "w")
  on.exit(base::close(con), add = TRUE)
  for (i in base::seq_along(ids)) {
    base::writeLines(base::paste0(">", ids[[i]]), con = con)
    base::writeLines(seqs[[i]], con = con)
  }
  invisible(path)
}

.dnmb_mrnacal_find_tool <- function(path = NULL, tool = "RNAfold") {
  if (!base::is.null(path) && base::nzchar(base::as.character(path)[1])) {
    path <- base::path.expand(base::as.character(path)[1])
    if (base::dir.exists(path)) {
      candidate <- base::file.path(path, tool)
      if (base::file.exists(candidate)) {
        return(candidate)
      }
    }
    if (base::file.exists(path)) {
      if (base::identical(base::basename(path), tool)) {
        return(path)
      }
      sibling <- base::file.path(base::dirname(path), tool)
      if (base::file.exists(sibling)) {
        return(sibling)
      }
    }
    found_path <- Sys.which(path)
    if (base::nzchar(found_path)) {
      if (base::identical(base::basename(found_path), tool)) {
        return(found_path)
      }
      sibling <- base::file.path(base::dirname(found_path), tool)
      if (base::file.exists(sibling)) {
        return(sibling)
      }
    }
    if (base::identical(base::basename(path), tool)) {
      return(path)
    }
  }
  direct_candidates <- base::file.path(c("/usr/local/bin", "/usr/bin", "/opt/conda/bin"), tool)
  direct_candidates <- direct_candidates[base::file.exists(direct_candidates)]
  if (base::length(direct_candidates)) {
    return(direct_candidates[[1]])
  }
  found <- Sys.which(tool)
  if (base::nzchar(found)) {
    return(found)
  }
  ""
}

.dnmb_mrnacal_parse_rnafold <- function(lines, id_map) {
  rows <- list()
  i <- 1L
  while (i <= base::length(lines)) {
    line <- base::trimws(lines[[i]])
    if (!base::startsWith(line, ">")) {
      i <- i + 1L
      next
    }
    id <- base::sub("^>", "", line)
    seq_line <- if (i + 1L <= base::length(lines)) base::trimws(lines[[i + 1L]]) else NA_character_
    struct_line <- if (i + 2L <= base::length(lines)) base::trimws(lines[[i + 2L]]) else NA_character_
    structure <- NA_character_
    mfe <- NA_real_
    if (!base::is.na(struct_line) && base::nzchar(struct_line)) {
      structure <- base::strsplit(struct_line, "\\s+")[[1]][[1]]
      energy_txt <- stringr::str_match(struct_line, "\\((\\s*[-+]?[0-9]+\\.?[0-9]*)\\s*\\)")[, 2]
      mfe <- suppressWarnings(base::as.numeric(energy_txt))
    }
    rows[[base::length(rows) + 1L]] <- data.frame(
      fold_id = id,
      locus_tag = id_map[[id]] %||% id,
      fold_sequence = seq_line,
      fold_structure = structure,
      fold_mfe = mfe,
      stringsAsFactors = FALSE
    )
    i <- i + 3L
  }
  if (!base::length(rows)) {
    return(data.frame(
      fold_id = character(),
      locus_tag = character(),
      fold_sequence = character(),
      fold_structure = character(),
      fold_mfe = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

.dnmb_mrnacal_run_rnafold_batch <- function(tbl, output_dir, rnafold_path = NULL, cpu = 1L) {
  tool <- .dnmb_mrnacal_find_tool(rnafold_path, "RNAfold")
  if (!base::nzchar(tool)) {
    return(list(ok = FALSE, error = "RNAfold not found in PATH.", table = data.frame()))
  }
  ids <- base::sprintf("MRNACAL_%06d", base::seq_len(base::nrow(tbl)))
  names(ids) <- tbl$locus_tag
  id_map <- stats::setNames(tbl$locus_tag, ids)
  fasta <- base::file.path(output_dir, "mrnacal_fold_sequences.fa")
  .dnmb_mrnacal_write_fasta(ids, tbl$fold_sequence, fasta)

  cpu <- suppressWarnings(base::as.integer(cpu)[1])
  if (base::is.na(cpu) || cpu < 1L) {
    cpu <- 1L
  }
  args <- c("--noPS", if (cpu > 1L) c("-j", base::as.character(cpu)) else character(), fasta)
  run <- dnmb_run_external(tool, args = args, required = FALSE)
  if (!base::isTRUE(run$ok) && cpu > 1L) {
    run <- dnmb_run_external(tool, args = c("--noPS", fasta), required = FALSE)
  }
  if (!base::isTRUE(run$ok)) {
    return(list(ok = FALSE, error = run$error %||% "RNAfold failed.", table = data.frame(), tool = tool))
  }
  parsed <- .dnmb_mrnacal_parse_rnafold(run$stdout, id_map = id_map)
  out_path <- base::file.path(output_dir, "mrnacal_fold_results.tsv")
  readr::write_tsv(parsed, out_path)
  list(ok = TRUE, table = parsed, file = out_path, tool = tool)
}

.dnmb_mrnacal_fasta_lines <- function(ids, seqs) {
  lines <- character()
  for (i in base::seq_along(ids)) {
    lines <- c(lines, base::paste0(">", ids[[i]]), base::as.character(seqs[[i]]))
  }
  lines
}

.dnmb_mrnacal_parse_lunp_file <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  lines <- base::readLines(path, warn = FALSE)
  lines <- lines[base::nzchar(base::trimws(lines))]
  data_lines <- lines[!base::startsWith(base::trimws(lines), "#")]
  if (!base::length(data_lines)) {
    return(data.frame())
  }
  tbl <- tryCatch(
    utils::read.table(
      text = base::paste(data_lines, collapse = "\n"),
      header = FALSE,
      fill = TRUE,
      stringsAsFactors = FALSE
    ),
    error = function(e) data.frame()
  )
  if (!base::nrow(tbl) || base::ncol(tbl) < 2L) {
    return(data.frame())
  }
  names(tbl) <- c("position", base::paste0("l", base::seq_len(base::ncol(tbl) - 1L)))
  tbl$position <- suppressWarnings(base::as.integer(tbl$position))
  for (col in base::setdiff(base::names(tbl), "position")) {
    tbl[[col]] <- suppressWarnings(base::as.numeric(tbl[[col]]))
  }
  tbl
}

.dnmb_mrnacal_lunp_region <- function(lunp, start, end) {
  if (!base::is.data.frame(lunp) || !base::nrow(lunp)) {
    return(NA_real_)
  }
  start <- suppressWarnings(base::as.integer(start)[1])
  end <- suppressWarnings(base::as.integer(end)[1])
  if (base::is.na(start) || base::is.na(end) || end < start) {
    return(NA_real_)
  }
  start <- base::max(1L, start)
  region_len <- end - start + 1L
  col <- base::paste0("l", region_len)
  if (col %in% base::names(lunp)) {
    idx <- base::match(end, lunp$position)
    if (!base::is.na(idx)) {
      val <- suppressWarnings(base::as.numeric(lunp[[col]][[idx]]))
      if (!base::is.na(val)) {
        return(val)
      }
    }
  }

  max_l <- suppressWarnings(base::as.integer(base::sub("^l", "", base::grep("^l[0-9]+$", base::names(lunp), value = TRUE))))
  max_l <- max_l[!base::is.na(max_l)]
  if (!base::length(max_l)) {
    return(NA_real_)
  }
  k <- base::max(max_l[max_l <= region_len], na.rm = TRUE)
  if (base::is.infinite(k) || base::is.na(k) || k < 1L) {
    return(NA_real_)
  }
  col <- base::paste0("l", k)
  ends <- base::seq.int(start + k - 1L, end)
  idx <- base::match(ends, lunp$position)
  vals <- suppressWarnings(base::as.numeric(lunp[[col]][idx[!base::is.na(idx)]]))
  vals <- vals[!base::is.na(vals)]
  if (!base::length(vals)) {
    return(NA_real_)
  }
  base::mean(vals)
}

.dnmb_mrnacal_run_rnaplfold_batch <- function(tbl,
                                              output_dir,
                                              rnaplfold_path = NULL,
                                              rnafold_path = NULL,
                                              window = 70L,
                                              span = 40L,
                                              ulength = 20L) {
  tool <- .dnmb_mrnacal_find_tool(rnaplfold_path %||% rnafold_path, "RNAplfold")
  if (!base::nzchar(tool)) {
    return(list(ok = FALSE, error = "RNAplfold not found in PATH.", table = data.frame()))
  }
  if (!base::nrow(tbl)) {
    return(list(ok = TRUE, table = data.frame(), tool = tool))
  }

  ids <- base::sprintf("MRNACAL_%06d", base::seq_len(base::nrow(tbl)))
  names(ids) <- tbl$locus_tag
  id_map <- stats::setNames(tbl$locus_tag, ids)
  work_dir <- base::file.path(output_dir, "mrnacal_rnaplfold")
  base::dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  seq_col <- if ("plfold_sequence" %in% base::names(tbl)) "plfold_sequence" else "fold_sequence"
  seqs <- base::as.character(tbl[[seq_col]])
  seq_len <- base::nchar(seqs)
  window <- suppressWarnings(base::as.integer(window)[1])
  span <- suppressWarnings(base::as.integer(span)[1])
  ulength <- suppressWarnings(base::as.integer(ulength)[1])
  if (base::is.na(window) || window < 10L) window <- 45L
  if (base::is.na(span) || span < 10L) span <- 35L
  if (base::is.na(ulength) || ulength < 1L) ulength <- 20L
  min_len <- base::max(10L, base::min(seq_len[seq_len >= 10L], na.rm = TRUE))
  window <- base::min(window, min_len)
  span <- base::min(span, window)
  ulength <- base::min(ulength, window)

  old_wd <- base::getwd()
  on.exit(base::setwd(old_wd), add = TRUE)
  base::setwd(work_dir)
  run <- tryCatch(
    base::system2(
      tool,
      args = c("--cutoff", "1.1", "-W", base::as.character(window), "-L", base::as.character(span), "-u", base::as.character(ulength)),
      input = .dnmb_mrnacal_fasta_lines(ids, seqs),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) structure(character(), status = 1L, error = conditionMessage(e))
  )
  status_code <- attr(run, "status")
  if (base::is.null(status_code)) {
    status_code <- 0L
  }
  if (!base::identical(base::as.integer(status_code), 0L)) {
    detail <- base::paste(run[base::nzchar(run)], collapse = "\n")
    if (!base::nzchar(detail)) {
      detail <- attr(run, "error") %||% "RNAplfold failed."
    }
    return(list(ok = FALSE, error = detail, table = data.frame(), tool = tool))
  }

  rows <- vector("list", base::nrow(tbl))
  for (i in base::seq_len(base::nrow(tbl))) {
    id <- ids[[i]]
    lunp <- .dnmb_mrnacal_parse_lunp_file(base::file.path(work_dir, base::paste0(id, "_lunp")))
    rbs_start <- if ("rbs_plfold_start" %in% base::names(tbl)) tbl$rbs_plfold_start[[i]] else tbl$rbs_start[[i]]
    rbs_end <- if ("rbs_plfold_end" %in% base::names(tbl)) tbl$rbs_plfold_end[[i]] else tbl$rbs_end[[i]]
    start_access_start <- if ("start_plfold_access_start" %in% base::names(tbl)) tbl$start_plfold_access_start[[i]] else tbl$start_access_start[[i]]
    start_access_end <- if ("start_plfold_access_end" %in% base::names(tbl)) tbl$start_plfold_access_end[[i]] else tbl$start_access_end[[i]]
    tir_access_start <- if ("tir_plfold_access_start" %in% base::names(tbl)) tbl$tir_plfold_access_start[[i]] else NA_integer_
    tir_access_end <- if ("tir_plfold_access_end" %in% base::names(tbl)) tbl$tir_plfold_access_end[[i]] else NA_integer_
    standby_start <- if ("standby_plfold_start" %in% base::names(tbl)) tbl$standby_plfold_start[[i]] else NA_integer_
    standby_end <- if ("standby_plfold_end" %in% base::names(tbl)) tbl$standby_plfold_end[[i]] else NA_integer_
    downstream_start <- if ("downstream_plfold_access_start" %in% base::names(tbl)) tbl$downstream_plfold_access_start[[i]] else NA_integer_
    downstream_end <- if ("downstream_plfold_access_end" %in% base::names(tbl)) tbl$downstream_plfold_access_end[[i]] else NA_integer_
    rbs_prob <- .dnmb_mrnacal_lunp_region(lunp, rbs_start, rbs_end)
    start_prob <- .dnmb_mrnacal_lunp_region(lunp, start_access_start, start_access_end)
    tir_prob <- .dnmb_mrnacal_lunp_region(lunp, tir_access_start, tir_access_end)
    standby_prob <- .dnmb_mrnacal_lunp_region(lunp, standby_start, standby_end)
    downstream_prob <- .dnmb_mrnacal_lunp_region(lunp, downstream_start, downstream_end)
    rows[[i]] <- data.frame(
      fold_id = id,
      locus_tag = id_map[[id]],
      rbs_plfold_unpaired_probability = rbs_prob,
      start_plfold_unpaired_probability = start_prob,
      tir_plfold_unpaired_probability = tir_prob,
      standby_plfold_unpaired_probability = standby_prob,
      downstream_plfold_unpaired_probability = downstream_prob,
      plfold_window = window,
      plfold_span = span,
      plfold_ulength = ulength,
      stringsAsFactors = FALSE
    )
  }
  out <- dplyr::bind_rows(rows)
  out$plfold_accessibility_score <- base::round(100 * base::mapply(
    function(rbs, start, tir, standby, downstream) {
      .dnmb_mrnacal_weighted_mean(
        c(rbs, start, tir, standby, downstream),
        c(0.26, 0.20, 0.30, 0.08, 0.16),
        fallback = NA_real_
      )
    },
    out$rbs_plfold_unpaired_probability,
    out$start_plfold_unpaired_probability,
    out$tir_plfold_unpaired_probability,
    out$standby_plfold_unpaired_probability,
    out$downstream_plfold_unpaired_probability
  ), 2)
  out_path <- base::file.path(output_dir, "mrnacal_plfold_results.tsv")
  readr::write_tsv(out, out_path)
  base::unlink(Sys.glob(base::file.path(work_dir, "*_dp.ps")), force = TRUE)
  list(ok = TRUE, table = out, file = out_path, tool = tool)
}

.dnmb_mrnacal_parse_rnaduplex <- function(lines, tbl, ids, anti_sd_rna) {
  rows <- vector("list", base::nrow(tbl))
  for (i in base::seq_len(base::nrow(tbl))) {
    line <- if (base::length(lines) >= i) base::trimws(lines[[i]]) else NA_character_
    energy <- NA_real_
    anti_from <- anti_to <- target_from <- target_to <- NA_integer_
    structure <- NA_character_
    if (!base::is.na(line) && base::nzchar(line)) {
      energy_txt <- stringr::str_match(line, "\\(\\s*([-+]?[0-9]+\\.?[0-9]*)\\s*\\)\\s*$")[, 2]
      energy <- suppressWarnings(base::as.numeric(energy_txt))
      range <- stringr::str_match(line, "([0-9]+)\\s*,\\s*([0-9]+)\\s*:\\s*([0-9]+)\\s*,\\s*([0-9]+)")
      anti_from <- suppressWarnings(base::as.integer(range[, 2]))
      anti_to <- suppressWarnings(base::as.integer(range[, 3]))
      target_from <- suppressWarnings(base::as.integer(range[, 4]))
      target_to <- suppressWarnings(base::as.integer(range[, 5]))
      structure <- base::strsplit(line, "\\s+")[[1]][[1]]
    }
    offset <- suppressWarnings(base::as.integer(tbl$duplex_target_offset[[i]]))
    duplex_start <- if (!base::is.na(target_from) && !base::is.na(offset)) offset + target_from - 1L else NA_integer_
    duplex_end <- if (!base::is.na(target_to) && !base::is.na(offset)) offset + target_to - 1L else NA_integer_
    motif <- if (!base::is.na(duplex_start) && !base::is.na(duplex_end) && duplex_end >= duplex_start) {
      base::substr(tbl$sequence_dna[[i]], duplex_start, duplex_end)
    } else {
      NA_character_
    }
    rows[[i]] <- data.frame(
      fold_id = ids[[i]],
      locus_tag = tbl$locus_tag[[i]],
      anti_sd_sequence = anti_sd_rna,
      duplex_structure = structure,
      duplex_energy = energy,
      duplex_anti_sd_start = anti_from,
      duplex_anti_sd_end = anti_to,
      duplex_target_start = duplex_start,
      duplex_target_end = duplex_end,
      duplex_motif = motif,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

.dnmb_mrnacal_duplex_score <- function(energy) {
  energy <- suppressWarnings(base::as.numeric(energy))
  if (base::is.na(energy)) {
    return(NA_real_)
  }
  base::round(base::pmax(0, base::pmin(100, 100 * ((-energy) - 1.5) / 7.5)), 2)
}

.dnmb_mrnacal_run_rnaduplex_batch <- function(tbl,
                                              output_dir,
                                              sd_seed = "AGGAGG",
                                              rnaduplex_path = NULL,
                                              rnafold_path = NULL) {
  tool <- .dnmb_mrnacal_find_tool(rnaduplex_path %||% rnafold_path, "RNAduplex")
  if (!base::nzchar(tool)) {
    return(list(ok = FALSE, error = "RNAduplex not found in PATH.", table = data.frame()))
  }
  if (!base::nrow(tbl)) {
    return(list(ok = TRUE, table = data.frame(), tool = tool))
  }

  anti_sd_rna <- .dnmb_mrnacal_dna_to_rna(.dnmb_mrnacal_revcomp(sd_seed))
  if (!base::nzchar(anti_sd_rna)) {
    anti_sd_rna <- "CCUCCU"
  }
  ids <- base::sprintf("MRNACAL_%06d", base::seq_len(base::nrow(tbl)))
  input <- character()
  work <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  work$duplex_target_offset <- NA_integer_
  for (i in base::seq_len(base::nrow(work))) {
    upstream_len <- suppressWarnings(base::as.integer(work$window_upstream[[i]]))
    sequence <- .dnmb_mrnacal_normalize_dna(work$sequence_dna[[i]])
    scan_left <- base::max(1L, upstream_len - 35L)
    scan_right <- base::max(1L, upstream_len - 3L)
    if (base::is.na(upstream_len) || upstream_len < 4L || scan_right < scan_left) {
      target <- ""
      scan_left <- 1L
    } else {
      target <- .dnmb_mrnacal_dna_to_rna(base::substr(sequence, scan_left, scan_right))
    }
    if (!base::nzchar(target)) {
      target <- "A"
    }
    work$duplex_target_offset[[i]] <- scan_left
    input <- c(input, anti_sd_rna, target)
  }
  run <- tryCatch(
    base::system2(tool, input = input, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(character(), status = 1L, error = conditionMessage(e))
  )
  status_code <- attr(run, "status")
  if (base::is.null(status_code)) {
    status_code <- 0L
  }
  if (!base::identical(base::as.integer(status_code), 0L)) {
    detail <- base::paste(run[base::nzchar(run)], collapse = "\n")
    if (!base::nzchar(detail)) {
      detail <- attr(run, "error") %||% "RNAduplex failed."
    }
    return(list(ok = FALSE, error = detail, table = data.frame(), tool = tool))
  }

  parsed <- .dnmb_mrnacal_parse_rnaduplex(run, work, ids = ids, anti_sd_rna = anti_sd_rna)
  parsed$duplex_score <- vapply(parsed$duplex_energy, .dnmb_mrnacal_duplex_score, numeric(1))
  out_path <- base::file.path(output_dir, "mrnacal_rnaduplex_results.tsv")
  readr::write_tsv(parsed, out_path)
  list(ok = TRUE, table = parsed, file = out_path, tool = tool)
}

.dnmb_mrnacal_pair_table <- function(structure) {
  structure <- base::as.character(structure)[1]
  if (base::is.na(structure) || !base::nzchar(structure)) {
    return(data.frame(i = integer(), j = integer()))
  }
  chars <- base::strsplit(structure, "", fixed = TRUE)[[1]]
  stack <- integer()
  pairs <- list()
  for (i in base::seq_along(chars)) {
    if (identical(chars[[i]], "(")) {
      stack <- c(stack, i)
    } else if (identical(chars[[i]], ")") && base::length(stack)) {
      left <- stack[[base::length(stack)]]
      stack <- stack[-base::length(stack)]
      pairs[[base::length(pairs) + 1L]] <- data.frame(i = left, j = i)
    }
  }
  if (!base::length(pairs)) {
    return(data.frame(i = integer(), j = integer()))
  }
  dplyr::bind_rows(pairs)
}

.dnmb_mrnacal_unpaired_fraction <- function(structure, start, end) {
  structure <- base::as.character(structure)[1]
  start <- base::as.integer(start)[1]
  end <- base::as.integer(end)[1]
  if (base::is.na(structure) || !base::nzchar(structure) || base::is.na(start) || base::is.na(end)) {
    return(NA_real_)
  }
  n <- base::nchar(structure)
  start <- base::max(1L, start)
  end <- base::min(n, end)
  if (end < start) {
    return(NA_real_)
  }
  chars <- base::strsplit(base::substr(structure, start, end), "", fixed = TRUE)[[1]]
  base::mean(chars == ".")
}

.dnmb_mrnacal_fold_score <- function(mfe, seq_len) {
  mfe <- suppressWarnings(base::as.numeric(mfe))
  seq_len <- suppressWarnings(base::as.numeric(seq_len))
  if (base::is.na(mfe) || base::is.na(seq_len) || seq_len <= 0) {
    return(NA_real_)
  }
  mfe_per_nt <- mfe / seq_len
  base::round(base::pmax(0, base::pmin(100, 100 * (mfe_per_nt + 0.6) / 0.6)), 2)
}

.dnmb_mrnacal_band <- function(score) {
  score <- suppressWarnings(base::as.numeric(score))
  ifelse(
    base::is.na(score), NA_character_,
    ifelse(score >= 70, "high",
           ifelse(score >= 50, "medium",
                  ifelse(score >= 35, "low", "poor")))
  )
}

.dnmb_mrnacal_normalize_hits <- function(results) {
  if (base::is.null(results) || !base::is.data.frame(results) || !base::nrow(results)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  tbl <- base::as.data.frame(results, stringsAsFactors = FALSE)
  data.frame(
    query = .dnmb_module_clean_annotation_key(tbl$locus_tag),
    source = "mRNAcal",
    family_system = "translation_initiation",
    family_id = base::as.character(tbl$tir_score_band),
    hit_label = base::paste0("Translation efficiency: ", base::as.character(tbl$tir_score_band)),
    enzyme_role = "translation_efficiency",
    evidence_mode = "sequence_model",
    substrate_label = NA_character_,
    support = base::as.character(tbl$support),
    typing_eligible = !base::is.na(tbl$tir_score),
    tir_score = suppressWarnings(base::as.numeric(tbl$tir_score)),
    rbs_score = suppressWarnings(base::as.numeric(tbl$rbs_score)),
    duplex_score = suppressWarnings(base::as.numeric(tbl$duplex_score)),
    accessibility_score = suppressWarnings(base::as.numeric(tbl$accessibility_score)),
    plfold_accessibility_score = suppressWarnings(base::as.numeric(tbl$plfold_accessibility_score)),
    upstream_au_score = suppressWarnings(base::as.numeric(tbl$upstream_au_score)),
    fold_mfe = suppressWarnings(base::as.numeric(tbl$fold_mfe)),
    stringsAsFactors = FALSE
  )
}

.dnmb_mrnacal_output_table <- function(genes, results) {
  if (!base::is.data.frame(genes)) {
    genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  }
  base_cols <- base::intersect(dnmb_backbone_columns(), base::names(genes))
  out <- base::as.data.frame(genes[, base_cols, drop = FALSE], stringsAsFactors = FALSE)
  out$locus_tag <- .dnmb_module_clean_annotation_key(out$locus_tag)

  cols <- c(
    "family_id", "hit_label", "tir_score", "tir_score_band",
    "rbs_motif", "rbs_seed", "rbs_mismatches", "rbs_spacer", "rbs_start", "rbs_end", "rbs_score",
    "anti_sd_sequence", "duplex_structure", "duplex_energy", "duplex_score",
    "duplex_target_start", "duplex_target_end", "duplex_motif",
    "start_codon", "start_codon_score", "early_k_score",
    "second_codon", "early_codons", "lysine_codon_count_2_8",
    "aaa_count_2_8", "aag_count_2_8", "early_coding_at_fraction",
    "ncs45_sequence", "ncs45_at_fraction", "ncs45_lysine_codon_count",
    "ncs45_aaa_count", "ncs45_aag_count",
    "early_aaa_run", "early_poly_a_run", "early_poly_a_penalty",
    "basic_aa_count_2_8", "skik_like",
    "upstream20_sequence", "upstream20_at_fraction", "upstream20_score",
    "upstream_enhancer_sequence", "upstream_enhancer_at_fraction",
    "upstream_enhancer_at_run", "upstream_au_score",
    "fold_mfe", "fold_mfe_per_nt", "fold_score",
    "rbs_unpaired_fraction", "start_unpaired_fraction",
    "rbs_plfold_unpaired_probability", "start_plfold_unpaired_probability",
    "tir_plfold_unpaired_probability", "standby_plfold_unpaired_probability",
    "downstream_plfold_unpaired_probability",
    "plfold_accessibility_score", "accessibility_score", "accessibility_method",
    "plfold_window", "plfold_span", "plfold_ulength",
    "plfold_local_upstream", "plfold_local_downstream",
    "fold_structure", "fold_sequence",
    "internal_sd_count", "internal_sd_motifs", "internal_sd_min_position", "internal_sd_penalty",
    "tir_score_percentile",
    "cai", "cai_score", "tai", "tai_score",
    "codon_efficiency_score", "cai_reference_set", "cai_reference_size",
    "trna_total_gcn", "informative_codon_count",
    "window_upstream", "window_downstream", "support"
  )
  for (col in cols) {
    out[[col]] <- NA
  }
  if (base::is.null(results) || !base::is.data.frame(results) || !base::nrow(results)) {
    return(out)
  }
  results$locus_tag <- .dnmb_module_clean_annotation_key(results$locus_tag)
  idx <- base::match(out$locus_tag, results$locus_tag)
  keep <- !base::is.na(idx)
  for (col in base::intersect(cols, base::names(results))) {
    out[[col]][keep] <- results[[col]][idx[keep]]
  }
  out
}

dnmb_run_mrnacal_module <- function(genes,
                                    output_dir,
                                    genbank = NULL,
                                    upstream = 60L,
                                    downstream = 60L,
                                    rnafold_path = NULL,
                                    require_rnafold = TRUE,
                                    sd_seed = NULL,
                                    translation_domain = NULL,
                                    cpu = 1L,
                                    top_folds = 12L,
                                    verbose = FALSE) {
  if (!base::is.data.frame(genes)) {
    genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  }
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- list()

  required <- c("locus_tag", "start", "end", "direction")
  missing <- base::setdiff(required, base::names(genes))
  if (base::length(missing)) {
    return(list(
      ok = FALSE,
      status = .dnmb_mrnacal_status_row("mrnacal_input", "failed", base::paste("Missing columns:", base::paste(missing, collapse = ", "))),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame(),
      files = list()
    ))
  }

  rnafold_tool <- .dnmb_mrnacal_find_tool(rnafold_path, "RNAfold")
  if (!base::nzchar(rnafold_tool) && base::isTRUE(require_rnafold)) {
    return(list(
      ok = FALSE,
      status = .dnmb_mrnacal_status_row("RNAfold", "failed", "RNAfold not found in PATH."),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame(),
      files = list()
    ))
  }

  genbank <- genbank %||% .dnmb_module_detect_genbank(base::getwd())
  contig_cols <- base::intersect(c("contig_number", "contig"), base::names(genes))
  contig_hint <- if (base::length(contig_cols)) genes[, contig_cols, drop = FALSE] else genes[0, , drop = FALSE]
  contigs <- .dnmb_mrnacal_contig_table(contig_hint, genbank = genbank)
  if (!base::nrow(contigs)) {
    return(list(
      ok = FALSE,
      status = .dnmb_mrnacal_status_row("mrnacal_genome_sequence", "missing", "No contig sequences were available from GenBank or global contig objects."),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame(),
      files = list()
    ))
  }

  translation_domain <- translation_domain %||% .dnmb_detect_translation_domain(target = genes, gb_path = genbank)
  seed_cols <- base::intersect(c("product", "rearranged_nt_seq"), base::names(genes))
  seed_genes <- if (base::length(seed_cols)) genes[, seed_cols, drop = FALSE] else genes[0, , drop = FALSE]
  seed <- .dnmb_mrnacal_sd_seed(seed_genes, translation_domain = translation_domain, sd_seed = sd_seed)
  status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("mrnacal_sd_seed", "ok", seed)

  trna_gcn <- .dnmb_mrnacal_trna_gcn(genes, contigs = contigs)
  trna_total <- base::as.integer(base::sum(trna_gcn, na.rm = TRUE))
  tai_weights <- .dnmb_mrnacal_tai_weights(trna_gcn)
  if (trna_total > 0L) {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("mrnacal_trna", "ok", base::paste0(trna_total, " tRNA genes annotated"))
  } else {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("mrnacal_trna", "missing", "No tRNA anticodon annotations found; tAI will be NA.")
  }

  calc_cols <- base::intersect(
    c("locus_tag", "start", "end", "direction", "translation", "contig_number", "contig", "product", "rearranged_nt_seq"),
    base::names(genes)
  )
  gene_tbl <- base::as.data.frame(genes[, calc_cols, drop = FALSE], stringsAsFactors = FALSE)
  gene_tbl$locus_tag <- .dnmb_module_clean_annotation_key(gene_tbl$locus_tag)
  gene_tbl <- gene_tbl[!base::duplicated(gene_tbl$locus_tag), , drop = FALSE]
  if ("translation" %in% base::names(gene_tbl)) {
    tr <- .dnmb_normalize_translation(gene_tbl$translation)
    gene_tbl <- gene_tbl[!base::is.na(tr) & base::nzchar(tr), , drop = FALSE]
    gene_tbl$translation <- tr[!base::is.na(tr) & base::nzchar(tr)]
  }
  gene_tbl <- gene_tbl[!base::is.na(gene_tbl$locus_tag) & base::nzchar(gene_tbl$locus_tag), , drop = FALSE]
  if (!base::nrow(gene_tbl)) {
    empty <- .dnmb_mrnacal_output_table(genes, data.frame())
    return(list(
      ok = TRUE,
      status = dplyr::bind_rows(status, .dnmb_mrnacal_status_row("mrnacal_cds", "empty", "No protein-coding genes with translations were available.")),
      results = data.frame(),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = empty,
      files = list()
    ))
  }

  windows <- lapply(base::seq_len(base::nrow(gene_tbl)), function(i) {
    .dnmb_mrnacal_transcript_window(gene_tbl[i, , drop = FALSE], contigs, upstream = upstream, downstream = downstream)
  })
  keep <- !vapply(windows, base::is.null, logical(1))
  gene_tbl <- gene_tbl[keep, , drop = FALSE]
  windows <- windows[keep]
  if (!base::nrow(gene_tbl)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_mrnacal_status_row("mrnacal_windows", "failed", "No transcript windows could be extracted.")),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame(),
      files = list()
    ))
  }

  full_cds_list <- vapply(base::seq_len(base::nrow(gene_tbl)), function(i) {
    .dnmb_mrnacal_full_cds(gene_tbl[i, , drop = FALSE], contigs)
  }, character(1))
  cai_ref_info <- .dnmb_mrnacal_select_cai_reference(gene_tbl)
  cai_ref_mask <- cai_ref_info$mask
  cai_ref_label <- cai_ref_info$label
  if (base::length(cai_ref_mask) != base::length(full_cds_list)) {
    cai_ref_mask <- base::rep(TRUE, base::length(full_cds_list))
    cai_ref_label <- "all_cds"
  }
  cai_ref_seqs <- full_cds_list[cai_ref_mask]
  cai_ref_seqs <- cai_ref_seqs[!base::is.na(cai_ref_seqs) & base::nzchar(cai_ref_seqs)]
  cai_ref_size <- base::length(cai_ref_seqs)
  cai_weights <- if (cai_ref_size > 0L) {
    .dnmb_mrnacal_cai_weights(.dnmb_mrnacal_codon_counts_total(cai_ref_seqs))
  } else {
    stats::setNames(base::rep(NA_real_, 64L), .dnmb_mrnacal_all_codons())
  }
  if (cai_ref_size > 0L) {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("mrnacal_cai", "ok", base::paste0("reference=", cai_ref_label, " (", cai_ref_size, " CDS)"))
  } else {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("mrnacal_cai", "missing", "No reference CDS available; CAI will be NA.")
  }

  base_rows <- vector("list", base::nrow(gene_tbl))
  for (i in base::seq_len(base::nrow(gene_tbl))) {
    win <- windows[[i]]
    row <- gene_tbl[i, , drop = FALSE]
    rbs <- .dnmb_mrnacal_best_rbs(win$sequence_dna, win$upstream_len_observed, sd_seed = seed)
    start_score <- .dnmb_mrnacal_start_score(win$start_codon)
    cds_dna <- base::substr(win$sequence_dna, win$upstream_len_observed + 1L, base::nchar(win$sequence_dna))
    full_cds_dna <- full_cds_list[[i]]
    cai_value <- .dnmb_mrnacal_cai(full_cds_dna, cai_weights)
    tai_value <- .dnmb_mrnacal_tai(full_cds_dna, tai_weights)
    cds_codon_counts <- .dnmb_mrnacal_codon_counts(full_cds_dna)
    informative_codons <- base::sum(cds_codon_counts[!base::is.na(cai_weights)])
    early <- .dnmb_mrnacal_early_context(cds_dna, translation = if ("translation" %in% base::names(row)) row$translation[[1]] else NA_character_)
    internal_sd <- .dnmb_mrnacal_internal_sd(cds_dna, sd_seed = seed, scan_len = 60L)
    upstream_context <- .dnmb_mrnacal_upstream_context(win$sequence_dna, win$upstream_len_observed, rbs_start = rbs$rbs_start)
    local <- .dnmb_mrnacal_local_tir_window(win$sequence_dna, win$upstream_len_observed, upstream_flank = 35L, downstream_flank = 30L)
    start_coord <- win$upstream_len_observed + 1L
    seq_n <- base::nchar(win$sequence_dna)
    start_access_start <- base::max(1L, start_coord - 5L)
    start_access_end <- base::min(seq_n, start_coord + 14L)
    tir_access_start <- base::max(1L, start_coord - 18L)
    tir_access_end <- base::min(seq_n, start_coord + 9L)
    downstream_access_start <- base::min(seq_n, start_coord + 1L)
    downstream_access_end <- base::min(seq_n, start_coord + 30L)
    if (downstream_access_end < downstream_access_start) {
      downstream_access_start <- NA_integer_
      downstream_access_end <- NA_integer_
    }
    if (!base::is.na(rbs$rbs_start) && rbs$rbs_start > 1L) {
      standby_start <- base::max(1L, rbs$rbs_start - 8L)
      standby_end <- rbs$rbs_start - 1L
    } else {
      standby_start <- base::max(1L, start_coord - 35L)
      standby_end <- base::max(1L, start_coord - 19L)
    }
    to_local <- function(x) {
      x <- suppressWarnings(base::as.integer(x)[1])
      if (base::is.na(x)) NA_integer_ else x - local$offset
    }
    base_rows[[i]] <- data.frame(
      locus_tag = row$locus_tag[[1]],
      fold_sequence = win$sequence_rna,
      plfold_sequence = local$sequence_rna,
      plfold_offset = local$offset,
      plfold_local_upstream = local$upstream_len,
      plfold_local_downstream = local$downstream_len,
      sequence_dna = win$sequence_dna,
      window_upstream = win$upstream_len_observed,
      window_downstream = win$downstream_len_observed,
      start_access_start = start_access_start,
      start_access_end = start_access_end,
      tir_access_start = tir_access_start,
      tir_access_end = tir_access_end,
      downstream_access_start = downstream_access_start,
      downstream_access_end = downstream_access_end,
      standby_start = standby_start,
      standby_end = standby_end,
      rbs_plfold_start = to_local(rbs$rbs_start),
      rbs_plfold_end = to_local(rbs$rbs_end),
      start_plfold_access_start = to_local(start_access_start),
      start_plfold_access_end = to_local(start_access_end),
      tir_plfold_access_start = to_local(tir_access_start),
      tir_plfold_access_end = to_local(tir_access_end),
      downstream_plfold_access_start = to_local(downstream_access_start),
      downstream_plfold_access_end = to_local(downstream_access_end),
      standby_plfold_start = to_local(standby_start),
      standby_plfold_end = to_local(standby_end),
      start_codon = win$start_codon,
      start_codon_score = start_score,
      rbs_motif = rbs$motif,
      rbs_seed = rbs$seed,
      rbs_mismatches = rbs$mismatches,
      rbs_spacer = rbs$spacer,
      rbs_start = rbs$rbs_start,
      rbs_end = rbs$rbs_end,
      rbs_score = rbs$score,
      early_k_score = early$score,
      second_codon = early$second_codon,
      early_codons = early$early_codons,
      lysine_codon_count_2_8 = early$lysine_codon_count_2_8,
      aaa_count_2_8 = early$aaa_count_2_8,
      aag_count_2_8 = early$aag_count_2_8,
      early_coding_at_fraction = early$early_coding_at_fraction,
      ncs45_sequence = early$ncs45_sequence,
      ncs45_at_fraction = early$ncs45_at_fraction,
      ncs45_lysine_codon_count = early$ncs45_lysine_codon_count,
      ncs45_aaa_count = early$ncs45_aaa_count,
      ncs45_aag_count = early$ncs45_aag_count,
      early_aaa_run = early$early_aaa_run,
      early_poly_a_run = early$early_poly_a_run,
      early_poly_a_penalty = early$early_poly_a_penalty,
      basic_aa_count_2_8 = early$basic_aa_count_2_8,
      skik_like = early$skik_like,
      internal_sd_count = internal_sd$internal_sd_count,
      internal_sd_motifs = internal_sd$internal_sd_motifs,
      internal_sd_positions = internal_sd$internal_sd_positions,
      internal_sd_min_position = internal_sd$internal_sd_min_position,
      internal_sd_penalty = internal_sd$internal_sd_penalty,
      cai = cai_value,
      cai_score = base::round(100 * cai_value, 2),
      tai = tai_value,
      tai_score = base::round(100 * tai_value, 2),
      cai_reference_set = cai_ref_label,
      cai_reference_size = cai_ref_size,
      trna_total_gcn = trna_total,
      informative_codon_count = base::as.integer(informative_codons),
      upstream20_sequence = upstream_context$upstream20_sequence,
      upstream20_at_fraction = upstream_context$upstream20_at_fraction,
      upstream20_score = upstream_context$upstream20_score,
      upstream_enhancer_sequence = upstream_context$enhancer_sequence,
      upstream_enhancer_at_fraction = upstream_context$enhancer_at_fraction,
      upstream_enhancer_at_run = upstream_context$enhancer_at_run,
      upstream_au_score = upstream_context$upstream_au_score,
      stringsAsFactors = FALSE
    )
  }
  results <- dplyr::bind_rows(base_rows)
  results <- results[base::nchar(results$fold_sequence) >= 10L, , drop = FALSE]
  if (!base::nrow(results)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_mrnacal_status_row("mrnacal_windows", "failed", "Extracted transcript windows were too short.")),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame(),
      files = list()
    ))
  }

  invisible(base::gc(verbose = FALSE))
  fold <- .dnmb_mrnacal_run_rnafold_batch(results[, c("locus_tag", "fold_sequence"), drop = FALSE], output_dir, rnafold_path = rnafold_tool, cpu = cpu)
  if (!base::isTRUE(fold$ok)) {
    if (base::isTRUE(require_rnafold)) {
      return(list(
        ok = FALSE,
        status = dplyr::bind_rows(status, .dnmb_mrnacal_status_row("RNAfold", "failed", fold$error)),
        hits = .dnmb_module_empty_optional_long_table(),
        output_table = data.frame(),
        files = list()
      ))
    }
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("RNAfold", "missing", fold$error)
    fold_tbl <- data.frame(
      locus_tag = results$locus_tag,
      fold_sequence = results$fold_sequence,
      fold_structure = NA_character_,
      fold_mfe = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("RNAfold", "ok", fold$tool %||% "RNAfold")
    fold_tbl <- fold$table
  }
  results <- dplyr::left_join(results, fold_tbl[, c("locus_tag", "fold_structure", "fold_mfe"), drop = FALSE], by = "locus_tag")
  results$fold_mfe_per_nt <- suppressWarnings(results$fold_mfe / base::nchar(results$fold_sequence))
  results$fold_score <- mapply(.dnmb_mrnacal_fold_score, results$fold_mfe, base::nchar(results$fold_sequence))
  results$rbs_unpaired_fraction <- mapply(.dnmb_mrnacal_unpaired_fraction, results$fold_structure, results$rbs_start, results$rbs_end)
  results$start_unpaired_fraction <- mapply(.dnmb_mrnacal_unpaired_fraction, results$fold_structure, results$start_access_start, results$start_access_end)
  structure_accessibility_score <- base::round(100 * rowMeans(
    cbind(results$rbs_unpaired_fraction, results$start_unpaired_fraction),
    na.rm = TRUE
  ), 2)
  structure_accessibility_score[base::is.nan(structure_accessibility_score)] <- NA_real_

  plfold <- .dnmb_mrnacal_run_rnaplfold_batch(
    results[, c(
      "locus_tag", "plfold_sequence",
      "rbs_plfold_start", "rbs_plfold_end",
      "start_plfold_access_start", "start_plfold_access_end",
      "tir_plfold_access_start", "tir_plfold_access_end",
      "downstream_plfold_access_start", "downstream_plfold_access_end",
      "standby_plfold_start", "standby_plfold_end"
    ), drop = FALSE],
    output_dir = output_dir,
    rnafold_path = rnafold_tool
  )
  if (base::isTRUE(plfold$ok)) {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("RNAplfold", "ok", plfold$tool %||% "RNAplfold")
    results <- dplyr::left_join(
      results,
      plfold$table[, c(
        "locus_tag",
        "rbs_plfold_unpaired_probability",
        "start_plfold_unpaired_probability",
        "tir_plfold_unpaired_probability",
        "standby_plfold_unpaired_probability",
        "downstream_plfold_unpaired_probability",
        "plfold_accessibility_score",
        "plfold_window", "plfold_span", "plfold_ulength"
      ), drop = FALSE],
      by = "locus_tag"
    )
  } else {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("RNAplfold", "missing", plfold$error %||% "RNAplfold failed.")
    results$rbs_plfold_unpaired_probability <- NA_real_
    results$start_plfold_unpaired_probability <- NA_real_
    results$tir_plfold_unpaired_probability <- NA_real_
    results$standby_plfold_unpaired_probability <- NA_real_
    results$downstream_plfold_unpaired_probability <- NA_real_
    results$plfold_accessibility_score <- NA_real_
    results$plfold_window <- NA_integer_
    results$plfold_span <- NA_integer_
    results$plfold_ulength <- NA_integer_
  }
  results$accessibility_method <- ifelse(
    !base::is.na(results$plfold_accessibility_score),
    "RNAplfold_unpaired_probability",
    ifelse(!base::is.na(structure_accessibility_score), "RNAfold_dotbracket_fraction", NA_character_)
  )
  results$accessibility_score <- ifelse(
    !base::is.na(results$plfold_accessibility_score),
    results$plfold_accessibility_score,
    structure_accessibility_score
  )

  duplex <- .dnmb_mrnacal_run_rnaduplex_batch(
    results[, c("locus_tag", "sequence_dna", "window_upstream"), drop = FALSE],
    output_dir = output_dir,
    sd_seed = seed,
    rnafold_path = rnafold_tool
  )
  duplex_cols <- c(
    "locus_tag", "anti_sd_sequence", "duplex_structure", "duplex_energy", "duplex_score",
    "duplex_target_start", "duplex_target_end", "duplex_motif"
  )
  if (base::isTRUE(duplex$ok)) {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("RNAduplex", "ok", duplex$tool %||% "RNAduplex")
    results <- dplyr::left_join(results, duplex$table[, duplex_cols, drop = FALSE], by = "locus_tag")
  } else {
    status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("RNAduplex", "missing", duplex$error %||% "RNAduplex failed.")
    for (col in base::setdiff(duplex_cols, "locus_tag")) {
      results[[col]] <- if (col %in% c("duplex_energy", "duplex_score", "duplex_target_start", "duplex_target_end")) NA_real_ else NA_character_
    }
  }

  fold_score_for_total <- ifelse(base::is.na(results$fold_score), 50, results$fold_score)
  accessibility_for_total <- ifelse(base::is.na(results$accessibility_score), 50, results$accessibility_score)
  duplex_for_total <- ifelse(base::is.na(results$duplex_score), results$rbs_score, results$duplex_score)
  upstream_au_for_total <- ifelse(base::is.na(results$upstream_au_score), 50, results$upstream_au_score)
  internal_sd_penalty_for_total <- ifelse(base::is.na(results$internal_sd_penalty), 0, results$internal_sd_penalty)
  raw_tir_score <- 0.20 * results$rbs_score +
    0.18 * duplex_for_total +
    0.32 * accessibility_for_total +
    0.10 * upstream_au_for_total +
    0.10 * results$start_codon_score +
    0.07 * results$early_k_score +
    0.03 * fold_score_for_total -
    internal_sd_penalty_for_total
  results$tir_score <- base::round(base::pmin(100, base::pmax(0, raw_tir_score)), 2)
  results$tir_score_percentile <- if (base::sum(!base::is.na(results$tir_score)) >= 2L) {
    base::round(100 * stats::ecdf(results$tir_score)(results$tir_score), 2)
  } else {
    NA_real_
  }
  results$tir_score_band <- .dnmb_mrnacal_band(results$tir_score)
  results$family_id <- results$tir_score_band
  results$hit_label <- base::paste0("Translation efficiency: ", results$tir_score_band)
  cai_for_ce <- ifelse(base::is.na(results$cai_score), NA_real_, results$cai_score)
  tai_for_ce <- ifelse(base::is.na(results$tai_score), NA_real_, results$tai_score)
  results$codon_efficiency_score <- base::round(base::rowMeans(
    base::cbind(cai_for_ce, tai_for_ce),
    na.rm = TRUE
  ), 2)
  results$codon_efficiency_score[base::is.nan(results$codon_efficiency_score)] <- NA_real_
  results$support <- base::paste0(
    "window=-", results$window_upstream, "/+", results$window_downstream,
    "; sd_seed=", results$rbs_seed,
    "; spacer=", results$rbs_spacer,
    "; anti_sd=", results$anti_sd_sequence,
    "; accessibility=", results$accessibility_method,
    "; plfold_window=-", results$plfold_local_upstream, "/+", results$plfold_local_downstream,
    "; downstream_access=", base::round(100 * results$downstream_plfold_unpaired_probability, 1),
    "; upstream_AU=", results$upstream_au_score,
    "; earlyK=", results$lysine_codon_count_2_8,
    "; NCS45_K=", results$ncs45_lysine_codon_count,
    "; internal_SD=", results$internal_sd_count,
    "; CAI=", base::round(results$cai, 3),
    "; tAI=", base::round(results$tai, 3),
    "; codon_eff=", results$codon_efficiency_score,
    "; cai_ref=", results$cai_reference_set,
    "; tRNA_GCN=", results$trna_total_gcn,
    "; SKIK_like=", results$skik_like,
    "; pct=", results$tir_score_percentile,
    "; fold=RNAfold_MFE",
    "; score=0.20*RBS+0.18*antiSD_duplex+0.32*local_RNAplfold_access+0.10*upstreamAU+0.10*start+0.07*earlyK+0.03*MFE-internal_SD_penalty (CAI/tAI reported separately)"
  )

  result_path <- base::file.path(output_dir, "mrnacal_translation_efficiency.tsv")
  readr::write_tsv(results, result_path)
  summary_path <- base::file.path(output_dir, "mrnacal_score_summary.tsv")
  score_summary <- results |>
    dplyr::count(.data$tir_score_band, name = "gene_count") |>
    dplyr::arrange(dplyr::desc(.data$gene_count))
  readr::write_tsv(score_summary, summary_path)

  hits <- .dnmb_mrnacal_normalize_hits(results)
  output_table <- .dnmb_mrnacal_output_table(genes, results)
  files <- list(results = result_path, summary = summary_path)
  if (base::isTRUE(fold$ok) && !base::is.null(fold$file)) {
    files$fold_results <- fold$file
    files$fold_sequences <- base::file.path(output_dir, "mrnacal_fold_sequences.fa")
  }
  if (base::isTRUE(plfold$ok) && !base::is.null(plfold$file)) {
    files$plfold_results <- plfold$file
  }
  if (base::isTRUE(duplex$ok) && !base::is.null(duplex$file)) {
    files$duplex_results <- duplex$file
  }

  status[[base::length(status) + 1L]] <- .dnmb_mrnacal_status_row("mrnacal_results", "ok", base::paste0(base::nrow(results), " CDS windows scored"))
  list(
    ok = TRUE,
    status = dplyr::bind_rows(status),
    results = results,
    hits = hits,
    output_table = output_table,
    files = files
  )
}

.dnmb_mrnacal_arc_plot <- function(row) {
  sequence <- base::as.character(row$fold_sequence[[1]])
  structure <- base::as.character(row$fold_structure[[1]])
  if (base::is.na(sequence) || !base::nzchar(sequence) || base::is.na(structure) || !base::nzchar(structure)) {
    return(NULL)
  }
  n <- base::nchar(sequence)
  chars <- base::strsplit(sequence, "", fixed = TRUE)[[1]]
  pos <- data.frame(
    x = base::seq_len(n),
    y = 0,
    nt = chars,
    region = "UTR",
    stringsAsFactors = FALSE
  )
  up <- base::as.integer(row$window_upstream[[1]])
  if (!base::is.na(up) && up < n) {
    pos$region[pos$x > up] <- "CDS"
    pos$region[pos$x >= up + 1L & pos$x <= base::min(n, up + 3L)] <- "start"
  }
  rbs_start <- suppressWarnings(base::as.integer(row$rbs_start[[1]]))
  rbs_end <- suppressWarnings(base::as.integer(row$rbs_end[[1]]))
  if (!base::is.na(rbs_start) && !base::is.na(rbs_end)) {
    pos$region[pos$x >= rbs_start & pos$x <= rbs_end] <- "RBS"
  }
  pairs <- .dnmb_mrnacal_pair_table(structure)
  title <- base::paste0(
    row$locus_tag[[1]],
    " | score=", base::round(base::as.numeric(row$tir_score[[1]]), 1),
    " | MFE=", base::round(base::as.numeric(row$fold_mfe[[1]]), 1)
  )
  ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = pairs,
      ggplot2::aes(x = .data$i, xend = .data$j, y = 0, yend = 0),
      curvature = -0.45,
      linewidth = 0.25,
      alpha = 0.38,
      color = "#4B5563",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = pos,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$region),
      shape = 21,
      size = 1.65,
      color = "white",
      stroke = 0.15
    ) +
    ggplot2::scale_fill_manual(values = c(UTR = "#94A3B8", CDS = "#64748B", RBS = "#16A34A", start = "#DC2626")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::coord_cartesian(ylim = c(-0.12, 1.0), clip = "off") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 8, hjust = 0),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}

.dnmb_plot_mrnacal_module <- function(genbank_table, output_dir = getwd(), top_n = 12L) {
  if (!is.data.frame(genbank_table) || !"mRNAcal_tir_score" %in% names(genbank_table)) {
    return(NULL)
  }
  tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
  tbl$tir_score <- suppressWarnings(as.numeric(tbl$mRNAcal_tir_score))
  keep <- !is.na(tbl$tir_score)
  if (!any(keep)) {
    return(NULL)
  }
  tbl <- tbl[keep, , drop = FALSE]
  tbl$tir_score_band <- factor(
    as.character(tbl$mRNAcal_tir_score_band),
    levels = c("high", "medium", "low", "poor")
  )
  band_cols <- c(high = "#15803D", medium = "#2563EB", low = "#D97706", poor = "#B91C1C")
  p_hist <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$tir_score, fill = .data$tir_score_band)) +
    ggplot2::geom_histogram(binwidth = 5, color = "white", linewidth = 0.2, boundary = 0) +
    ggplot2::scale_fill_manual(values = band_cols, drop = FALSE, na.value = "#9CA3AF") +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    ggplot2::labs(title = "Translation Initiation Score", x = "Composite score", y = "Genes", fill = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))

  scatter_tbl <- tbl
  scatter_tbl$fold_mfe <- suppressWarnings(as.numeric(scatter_tbl$mRNAcal_fold_mfe))
  scatter_tbl$duplex_energy <- if ("mRNAcal_duplex_energy" %in% names(scatter_tbl)) suppressWarnings(as.numeric(scatter_tbl$mRNAcal_duplex_energy)) else NA_real_
  scatter_tbl$local_accessibility <- if ("mRNAcal_accessibility_score" %in% names(scatter_tbl)) suppressWarnings(as.numeric(scatter_tbl$mRNAcal_accessibility_score)) else NA_real_
  p_scatter <- ggplot2::ggplot(scatter_tbl, ggplot2::aes(x = .data$local_accessibility, y = .data$tir_score, color = .data$tir_score_band)) +
    ggplot2::geom_point(alpha = 0.72, size = 1.4, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = band_cols, drop = FALSE, na.value = "#9CA3AF") +
    ggplot2::labs(title = "Local Accessibility vs Score", x = "RNAplfold/local accessibility score", y = "Composite score", color = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))

  component_cols <- c(
    RBS = "mRNAcal_rbs_score",
    AntiSD = "mRNAcal_duplex_score",
    Accessibility = "mRNAcal_accessibility_score",
    UpstreamAU = "mRNAcal_upstream_au_score",
    Start = "mRNAcal_start_codon_score",
    EarlyK = "mRNAcal_early_k_score",
    MFE = "mRNAcal_fold_score"
  )
  top_component <- tbl[order(-tbl$tir_score), , drop = FALSE]
  top_component <- utils::head(top_component, min(20L, nrow(top_component)))
  comp_rows <- list()
  for (nm in names(component_cols)) {
    col <- component_cols[[nm]]
    if (col %in% names(top_component)) {
      comp_rows[[length(comp_rows) + 1L]] <- data.frame(
        locus_tag = top_component$locus_tag,
        component = nm,
        score = suppressWarnings(as.numeric(top_component[[col]])),
        stringsAsFactors = FALSE
      )
    }
  }
  comp_tbl <- if (length(comp_rows)) dplyr::bind_rows(comp_rows) else data.frame()
  comp_tbl$locus_tag <- factor(comp_tbl$locus_tag, levels = rev(unique(top_component$locus_tag)))
  p_comp <- ggplot2::ggplot(comp_tbl, ggplot2::aes(x = .data$score, y = .data$locus_tag, fill = .data$component)) +
    ggplot2::geom_col(position = "dodge", width = 0.78, na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
    ggplot2::scale_fill_manual(values = c(RBS = "#16A34A", AntiSD = "#0F766E", Accessibility = "#0891B2", UpstreamAU = "#CA8A04", Start = "#DC2626", EarlyK = "#9333EA", MFE = "#64748B")) +
    ggplot2::labs(title = "Top Gene Components", x = "Component score", y = NULL, fill = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  summary_pdf <- file.path(plot_dir, "mRNAcal_translation_efficiency.pdf")
  summary_plot <- cowplot::plot_grid(p_hist, p_scatter, p_comp, ncol = 1, rel_heights = c(0.86, 0.86, 1.2))
  ggplot2::ggsave(summary_pdf, summary_plot, width = 10, height = 13, bg = "white")

  fold_pdf <- NULL
  if (all(c("mRNAcal_fold_structure", "mRNAcal_fold_sequence") %in% names(tbl))) {
    fold_tbl <- tbl[!is.na(tbl$mRNAcal_fold_structure) & nzchar(as.character(tbl$mRNAcal_fold_structure)), , drop = FALSE]
    if (nrow(fold_tbl)) {
      fold_tbl <- fold_tbl[order(-fold_tbl$tir_score), , drop = FALSE]
      fold_tbl <- utils::head(fold_tbl, min(as.integer(top_n), nrow(fold_tbl)))
      local_rows <- data.frame(
        locus_tag = fold_tbl$locus_tag,
        fold_sequence = fold_tbl$mRNAcal_fold_sequence,
        fold_structure = fold_tbl$mRNAcal_fold_structure,
        tir_score = fold_tbl$mRNAcal_tir_score,
        fold_mfe = fold_tbl$mRNAcal_fold_mfe,
        window_upstream = fold_tbl$mRNAcal_window_upstream,
        rbs_start = fold_tbl$mRNAcal_rbs_start %||% NA,
        rbs_end = fold_tbl$mRNAcal_rbs_end %||% NA,
        stringsAsFactors = FALSE
      )
      fold_plots <- lapply(seq_len(nrow(local_rows)), function(i) .dnmb_mrnacal_arc_plot(local_rows[i, , drop = FALSE]))
      fold_plots <- Filter(Negate(is.null), fold_plots)
      if (length(fold_plots)) {
        fold_pdf <- file.path(plot_dir, "mRNAcal_top_folds.pdf")
        fold_grid <- cowplot::plot_grid(plotlist = fold_plots, ncol = 1)
        ggplot2::ggsave(fold_pdf, fold_grid, width = 12, height = max(6, min(32, length(fold_plots) * 1.15)), bg = "white")
      }
    }
  }

  out <- list(pdf = summary_pdf)
  if (!is.null(fold_pdf)) {
    out$fold_pdf <- fold_pdf
  }
  out
}
