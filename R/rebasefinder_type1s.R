.dnmb_type1s_gold_url <- function() {
  "https://rebase.neb.com/rebase/Type_I_S_subunit_Gold_Standards_Protein.txt"
}

.dnmb_type1s_gold_filename <- function() {
  "Type_I_S_subunit_Gold_Standards_Protein.txt"
}

.dnmb_type1s_parser_version <- function() {
  "type1s-v1.2.0"
}

.dnmb_type1s_database_version <- function() {
  "separate-trd-v1.2.0"
}

.dnmb_type1s_prediction_version <- function() {
  "type1s-predict-v1.4.0"
}

.dnmb_type1s_spacer_model_version <- function() {
  "scaffold-ensemble-tael-v1.1.0"
}

.dnmb_type1s_clean_id <- function(x) {
  x <- base::trimws(base::as.character(x))
  x[base::is.na(x) | !base::nzchar(x)] <- NA_character_
  x
}

.dnmb_type1s_clean_protein <- function(x) {
  x <- base::toupper(base::as.character(x))
  x <- base::gsub("[[:space:]<>*.-]+", "", x)
  x <- base::gsub("[^A-Z]", "", x)
  x <- base::chartr("BJOUZ", "XXXXX", x)
  x[!base::nzchar(x)] <- NA_character_
  x
}

.dnmb_type1s_reverse_complement <- function(x) {
  x <- base::toupper(base::as.character(x))
  base::vapply(x, function(site) {
    if (base::is.na(site) || !base::nzchar(site)) return(NA_character_)
    complemented <- base::chartr(
      "ACGTRYWSMKBDHVN",
      "TGCAYRWSKMVHDBN",
      site
    )
    base::paste0(base::rev(base::strsplit(complemented, "", fixed = TRUE)[[1]]), collapse = "")
  }, character(1))
}

.dnmb_type1s_parse_motif <- function(recognition) {
  recognition <- base::toupper(base::gsub("[[:space:]^]", "", base::as.character(recognition)[1]))
  if (base::is.na(recognition) || !base::nzchar(recognition) || recognition == "?" ||
      base::grepl("[/,;]", recognition)) {
    return(NULL)
  }

  compact_match <- base::regexec(
    "^([ACGTRYWSMKBDHVN]+)\\(([0-9]{1,2})\\)([ACGTRYWSMKBDHVN]+)$",
    recognition,
    perl = TRUE
  )
  compact_parts <- base::regmatches(recognition, compact_match)[[1]]
  if (base::length(compact_parts) == 4L) {
    spacer <- base::as.integer(compact_parts[[3]])
    if (!base::is.na(spacer) && spacer >= 3L && spacer <= 15L) {
      return(list(
        left = compact_parts[[2]],
        spacer = spacer,
        right = compact_parts[[4]],
        recognition = base::paste0(compact_parts[[2]], base::strrep("N", spacer), compact_parts[[4]])
      ))
    }
  }

  if (!base::grepl("^[ACGTRYWSMKBDHVN]+$", recognition)) return(NULL)
  runs <- base::gregexpr("N+", recognition, perl = TRUE)[[1]]
  if (runs[[1]] < 0L) return(NULL)
  run_lengths <- base::attr(runs, "match.length")
  eligible <- base::which(
    run_lengths >= 3L & run_lengths <= 15L & runs > 1L &
      (runs + run_lengths) <= base::nchar(recognition)
  )
  if (!base::length(eligible)) return(NULL)

  longest <- base::max(run_lengths[eligible])
  candidates <- eligible[run_lengths[eligible] == longest]
  if (base::length(candidates) > 1L) {
    centers <- runs[candidates] + (run_lengths[candidates] - 1) / 2
    motif_center <- (base::nchar(recognition) + 1) / 2
    candidates <- candidates[base::which.min(base::abs(centers - motif_center))]
  }
  selected <- candidates[[1]]
  left <- base::substr(recognition, 1L, runs[selected] - 1L)
  right <- base::substr(
    recognition,
    runs[selected] + run_lengths[selected],
    base::nchar(recognition)
  )
  if (!base::nzchar(left) || !base::nzchar(right)) return(NULL)

  list(
    left = left,
    spacer = base::as.integer(run_lengths[selected]),
    right = right,
    recognition = base::paste0(left, base::strrep("N", run_lengths[selected]), right)
  )
}

.dnmb_type1s_valid_gold_file <- function(path) {
  if (base::is.null(path) || !base::file.exists(path)) return(FALSE)
  info <- base::file.info(path)
  if (base::is.na(info$size) || info$size < 100000) return(FALSE)
  headers <- base::grep("^>", base::readLines(path, warn = FALSE))
  base::length(headers) >= 500L
}

.dnmb_type1s_nonempty_file <- function(path) {
  if (base::is.null(path) || !base::file.exists(path)) return(FALSE)
  size <- base::file.info(path)$size
  !base::is.na(size) && size > 0L
}

.dnmb_type1s_longest_tetrapeptide_repeat <- function(sequence,
                                                      motifs = NULL) {
  sequence <- .dnmb_type1s_clean_protein(sequence)
  if (base::is.na(sequence) || base::nchar(sequence) < 8L) return(0L)
  protein_length <- base::nchar(sequence)
  longest <- 0L
  for (start in base::seq_len(protein_length - 7L)) {
    motif <- base::substr(sequence, start, start + 3L)
    if (!base::is.null(motifs) && !motif %in% motifs) next
    repeat_count <- 1L
    cursor <- start + 4L
    while (cursor + 3L <= protein_length &&
           base::identical(base::substr(sequence, cursor, cursor + 3L), motif)) {
      repeat_count <- repeat_count + 1L
      cursor <- cursor + 4L
    }
    if (repeat_count >= 2L) longest <- base::max(longest, repeat_count)
  }
  base::as.integer(longest)
}

.dnmb_type1s_scaffold_features <- function(sequence) {
  sequence <- .dnmb_type1s_clean_protein(sequence)
  if (base::is.na(sequence)) return(NULL)
  protein_length <- base::nchar(sequence)
  if (protein_length < 300L || protein_length > 700L) return(NULL)

  first_end <- base::round(protein_length * 0.11)
  central_start <- base::round(protein_length * 0.30) + 1L
  central_end <- base::round(protein_length * 0.66)
  distal_start <- base::round(protein_length * 0.82) + 1L
  if (first_end < 1L || central_end < central_start || distal_start > protein_length) {
    return(NULL)
  }
  first <- base::substr(sequence, 1L, first_end)
  central <- base::substr(sequence, central_start, central_end)
  distal <- base::substr(sequence, distal_start, protein_length)
  scaffold <- base::paste0(first, central, distal)
  residues <- base::strsplit(scaffold, "", fixed = TRUE)[[1]]
  helix_friendly <- base::mean(residues %in% c("A", "E", "K", "L", "Q", "R"))
  helix_breaker <- base::mean(residues %in% c("G", "P"))

  list(
    sequence = scaffold,
    length = base::nchar(scaffold),
    central_sequence = central,
    central_length = base::nchar(central),
    helix_friendly_fraction = helix_friendly,
    helix_breaker_fraction = helix_breaker,
    tetrapeptide_repeat_count = .dnmb_type1s_longest_tetrapeptide_repeat(central),
    tael_like_repeat_count = .dnmb_type1s_longest_tetrapeptide_repeat(
      central,
      motifs = c("TAEL", "LEAT", "SEAL", "TSEL")
    ),
    first_start = 1L,
    first_end = first_end,
    central_start = central_start,
    central_end = central_end,
    distal_start = distal_start,
    distal_end = protein_length,
    method = "conserved_scaffold_fixed_windows_v1"
  )
}

.dnmb_type1s_ensure_scaffold_columns <- function(table) {
  table <- base::as.data.frame(table, stringsAsFactors = FALSE)
  required <- c(
    "spacer_scaffold_sequence", "spacer_scaffold_length",
    "spacer_central_length", "spacer_scaffold_helix_friendly_fraction",
    "spacer_scaffold_helix_breaker_fraction", "spacer_tetrapeptide_repeat_count",
    "spacer_tael_like_repeat_count"
  )
  if (base::all(required %in% base::names(table))) return(table)
  if (!"sequence" %in% base::names(table)) return(table)
  features <- base::lapply(table$sequence, .dnmb_type1s_scaffold_features)
  value <- function(name, mode = c("character", "integer", "numeric")) {
    mode <- base::match.arg(mode)
    default <- switch(mode, character = NA_character_, integer = NA_integer_, numeric = NA_real_)
    base::vapply(features, function(x) {
      if (base::is.null(x)) return(default)
      x[[name]]
    }, switch(mode, character = character(1), integer = integer(1), numeric = numeric(1)))
  }
  table$spacer_scaffold_sequence <- value("sequence", "character")
  table$spacer_scaffold_length <- value("length", "integer")
  table$spacer_central_length <- value("central_length", "integer")
  table$spacer_scaffold_helix_friendly_fraction <- value(
    "helix_friendly_fraction",
    "numeric"
  )
  table$spacer_scaffold_helix_breaker_fraction <- value(
    "helix_breaker_fraction",
    "numeric"
  )
  table$spacer_tetrapeptide_repeat_count <- value("tetrapeptide_repeat_count", "integer")
  table$spacer_tael_like_repeat_count <- value("tael_like_repeat_count", "integer")
  table
}

.dnmb_type1s_download_gold <- function(cache_dir, force = FALSE, verbose = TRUE) {
  base::dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- base::file.path(cache_dir, .dnmb_type1s_gold_filename())
  if (!base::isTRUE(force) && .dnmb_type1s_valid_gold_file(destination)) return(destination)

  temporary <- base::tempfile("type1s-gold-", tmpdir = cache_dir, fileext = ".txt")
  on.exit(base::unlink(temporary, force = TRUE), add = TRUE)
  url <- .dnmb_type1s_gold_url()
  if (base::isTRUE(verbose)) message("[Type I-S] Downloading official REBASE Gold Standard...")

  downloaded <- base::tryCatch({
    utils::download.file(url, temporary, mode = "wb", quiet = !base::isTRUE(verbose), method = "libcurl")
    .dnmb_type1s_valid_gold_file(temporary)
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!base::isTRUE(downloaded)) {
    curl <- dnmb_run_external(
      "curl",
      args = c("-L", "--fail", "--silent", "--show-error", url, "-o", temporary),
      required = FALSE
    )
    downloaded <- base::isTRUE(curl$ok) && .dnmb_type1s_valid_gold_file(temporary)
  }
  if (!base::isTRUE(downloaded)) {
    if (.dnmb_type1s_valid_gold_file(destination)) return(destination)
    base::stop("Could not download the official REBASE Type I S Gold Standard file.", call. = FALSE)
  }

  if (!base::file.rename(temporary, destination)) {
    copied <- base::file.copy(temporary, destination, overwrite = TRUE)
    if (!base::isTRUE(copied)) base::stop("Could not cache the REBASE Type I S Gold file.", call. = FALSE)
  }
  destination
}

.dnmb_type1s_extract_windows <- function(sequence) {
  sequence <- .dnmb_type1s_clean_protein(sequence)
  if (base::is.na(sequence)) return(NULL)
  protein_length <- base::nchar(sequence)
  if (protein_length < 300L || protein_length > 700L) return(NULL)

  trd1_start <- base::max(1L, base::round(protein_length * 0.11) + 1L)
  trd1_end <- base::min(protein_length, base::round(protein_length * 0.30))
  trd2_start <- base::max(1L, base::round(protein_length * 0.66) + 1L)
  trd2_end <- base::min(protein_length, base::round(protein_length * 0.82))
  if (trd1_end <= trd1_start || trd2_end <= trd2_start) return(NULL)
  scaffold <- .dnmb_type1s_scaffold_features(sequence)
  if (base::is.null(scaffold)) return(NULL)

  list(
    sequence = sequence,
    length = protein_length,
    trd1 = base::substr(sequence, trd1_start, trd1_end),
    trd2 = base::substr(sequence, trd2_start, trd2_end),
    trd1_start = trd1_start,
    trd1_end = trd1_end,
    trd2_start = trd2_start,
    trd2_end = trd2_end,
    spacer_scaffold = scaffold$sequence,
    spacer_scaffold_length = scaffold$length,
    spacer_central_length = scaffold$central_length,
    spacer_scaffold_helix_friendly_fraction = scaffold$helix_friendly_fraction,
    spacer_scaffold_helix_breaker_fraction = scaffold$helix_breaker_fraction,
    spacer_tetrapeptide_repeat_count = scaffold$tetrapeptide_repeat_count,
    spacer_tael_like_repeat_count = scaffold$tael_like_repeat_count,
    spacer_scaffold_method = scaffold$method,
    boundary_method = "gold_cv_specificity_core_window_v1"
  )
}

.dnmb_type1s_parse_gold_file <- function(path) {
  if (!.dnmb_type1s_valid_gold_file(path)) {
    base::stop("Invalid REBASE Type I S Gold Standard file: ", path, call. = FALSE)
  }
  lines <- base::readLines(path, warn = FALSE)
  header_index <- base::grep("^>", lines)
  record_end <- c(header_index[-1L] - 1L, base::length(lines))

  records <- base::lapply(base::seq_along(header_index), function(i) {
    header <- base::trimws(base::sub("^>", "", lines[[header_index[[i]]]]))
    fields <- base::strsplit(header, "[[:space:]]+", perl = TRUE)[[1]]
    if (base::length(fields) < 2L) return(NULL)
    motif <- .dnmb_type1s_parse_motif(fields[[2]])
    if (base::is.null(motif)) return(NULL)
    sequence_lines <- if (header_index[[i]] < record_end[[i]]) {
      lines[(header_index[[i]] + 1L):record_end[[i]]]
    } else {
      character()
    }
    sequence <- .dnmb_type1s_clean_protein(base::paste(sequence_lines, collapse = ""))
    windows <- .dnmb_type1s_extract_windows(sequence)
    if (base::is.null(windows)) return(NULL)
    base::data.frame(
      source_enzyme = fields[[1]],
      sequence = windows$sequence,
      recognition = motif$recognition,
      left_half = motif$left,
      right_half = motif$right,
      spacer_length = motif$spacer,
      trd1_sequence = windows$trd1,
      trd2_sequence = windows$trd2,
      spacer_scaffold_sequence = windows$spacer_scaffold,
      spacer_scaffold_length = windows$spacer_scaffold_length,
      spacer_central_length = windows$spacer_central_length,
      spacer_scaffold_helix_friendly_fraction = windows$spacer_scaffold_helix_friendly_fraction,
      spacer_scaffold_helix_breaker_fraction = windows$spacer_scaffold_helix_breaker_fraction,
      spacer_tetrapeptide_repeat_count = windows$spacer_tetrapeptide_repeat_count,
      spacer_tael_like_repeat_count = windows$spacer_tael_like_repeat_count,
      stringsAsFactors = FALSE
    )
  })
  reference <- base::do.call(base::rbind, records[!base::vapply(records, base::is.null, logical(1))])
  if (base::is.null(reference) || !base::nrow(reference)) {
    base::stop("No usable canonical Type I S references were parsed.", call. = FALSE)
  }

  sequence_label_count <- stats::ave(
    reference$recognition,
    reference$sequence,
    FUN = function(x) base::length(base::unique(x))
  )
  reference <- reference[sequence_label_count == 1L, , drop = FALSE]
  reference <- reference[!base::duplicated(reference$sequence), , drop = FALSE]
  spacer_annotation <- function(trd_sequence) {
    groups <- base::split(reference$spacer_length, trd_sequence)
    spacer_sets <- base::vapply(
      groups,
      function(values) base::paste(base::sort(base::unique(values)), collapse = ","),
      character(1)
    )
    counts <- base::vapply(groups, function(values) base::length(base::unique(values)), integer(1))
    list(
      set = base::unname(spacer_sets[trd_sequence]),
      ambiguous = base::unname(counts[trd_sequence]) > 1L
    )
  }
  trd1_spacer <- spacer_annotation(reference$trd1_sequence)
  trd2_spacer <- spacer_annotation(reference$trd2_sequence)
  reference$trd1_spacer_set <- trd1_spacer$set
  reference$trd1_spacer_ambiguous <- trd1_spacer$ambiguous
  reference$trd2_spacer_set <- trd2_spacer$set
  reference$trd2_spacer_ambiguous <- trd2_spacer$ambiguous
  reference$reference_id <- base::sprintf("R%04d", base::seq_len(base::nrow(reference)))
  reference <- reference[, c(
    "reference_id", "source_enzyme", "sequence", "recognition", "left_half",
    "right_half", "spacer_length", "trd1_sequence", "trd2_sequence",
    "spacer_scaffold_sequence", "spacer_scaffold_length",
    "spacer_central_length", "spacer_scaffold_helix_friendly_fraction",
    "spacer_scaffold_helix_breaker_fraction", "spacer_tetrapeptide_repeat_count",
    "spacer_tael_like_repeat_count",
    "trd1_spacer_set", "trd1_spacer_ambiguous",
    "trd2_spacer_set", "trd2_spacer_ambiguous"
  ), drop = FALSE]
  base::rownames(reference) <- NULL
  reference
}

.dnmb_type1s_reference <- function(cache_dir, force = FALSE, download = TRUE, verbose = TRUE) {
  gold_path <- base::file.path(cache_dir, .dnmb_type1s_gold_filename())
  if ((base::isTRUE(force) || !.dnmb_type1s_valid_gold_file(gold_path)) && base::isTRUE(download)) {
    gold_path <- .dnmb_type1s_download_gold(cache_dir, force = force, verbose = verbose)
  }
  if (!.dnmb_type1s_valid_gold_file(gold_path)) {
    base::stop(
      "Official REBASE Type I S Gold file is unavailable. Expected: ", gold_path,
      call. = FALSE
    )
  }

  source_md5 <- base::unname(tools::md5sum(gold_path))
  cache_key <- base::paste(
    .dnmb_type1s_parser_version(),
    .dnmb_type1s_database_version(),
    .dnmb_type1s_spacer_model_version(),
    base::substr(source_md5, 1L, 12L),
    sep = "-"
  )
  reference_dir <- base::file.path(cache_dir, "type1s", cache_key)
  reference_rds <- base::file.path(reference_dir, "reference.rds")
  if (!base::isTRUE(force) && base::file.exists(reference_rds)) {
    cached <- base::tryCatch(base::readRDS(reference_rds), error = function(e) NULL)
    required_columns <- c(
      "reference_id", "source_enzyme", "sequence", "recognition",
      "left_half", "right_half", "spacer_length", "trd1_sequence",
      "trd2_sequence", "trd1_spacer_set", "trd1_spacer_ambiguous",
      "trd2_spacer_set", "trd2_spacer_ambiguous"
    )
    if (base::is.list(cached) && base::is.data.frame(cached$table) &&
        base::nrow(cached$table) &&
        base::all(required_columns %in% base::names(cached$table))) {
      # Cached payloads store absolute paths for provenance. Re-anchor those
      # paths so a copied cache remains usable in its new location.
      cached$source_path <- base::normalizePath(gold_path, winslash = "/", mustWork = FALSE)
      cached$source_url <- .dnmb_type1s_gold_url()
      cached$reference_dir <- reference_dir
      cached$version <- cache_key
      cached$parser_version <- .dnmb_type1s_parser_version()
      cached$database_version <- .dnmb_type1s_database_version()
      cached$spacer_model_version <- .dnmb_type1s_spacer_model_version()
      cached$table <- .dnmb_type1s_ensure_scaffold_columns(cached$table)
      return(cached)
    }
  }

  base::dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  if (base::isTRUE(force) && base::dir.exists(reference_dir)) {
    base::unlink(reference_dir, recursive = TRUE, force = TRUE)
    base::dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  }
  reference <- .dnmb_type1s_parse_gold_file(gold_path)
  payload <- list(
    table = reference,
    source_path = base::normalizePath(gold_path, winslash = "/", mustWork = FALSE),
    source_url = .dnmb_type1s_gold_url(),
    source_md5 = source_md5,
    version = cache_key,
    parser_version = .dnmb_type1s_parser_version(),
    database_version = .dnmb_type1s_database_version(),
    spacer_model_version = .dnmb_type1s_spacer_model_version(),
    built_at = base::format(base::Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    reference_dir = reference_dir
  )
  base::saveRDS(payload, reference_rds)
  payload
}

.dnmb_type1s_write_fasta <- function(ids, sequences, path) {
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)
  valid <- !base::is.na(sequences) & base::nzchar(sequences)
  lines <- base::as.vector(base::rbind(base::paste0(">", ids[valid]), sequences[valid]))
  temporary <- base::tempfile("type1s-fasta-", tmpdir = base::dirname(path))
  on.exit(base::unlink(temporary, force = TRUE), add = TRUE)
  base::writeLines(lines, temporary)
  if (!.dnmb_type1s_nonempty_file(temporary)) {
    base::stop("Could not write a non-empty Type I-S FASTA file.", call. = FALSE)
  }
  if (base::file.exists(path)) base::unlink(path, force = TRUE)
  if (!base::file.rename(temporary, path)) {
    copied <- base::file.copy(temporary, path, overwrite = TRUE)
    if (!base::isTRUE(copied)) base::stop("Could not install Type I-S FASTA file.", call. = FALSE)
  }
  invisible(path)
}

.dnmb_type1s_prepare_reference_files <- function(reference, reference_dir) {
  table <- .dnmb_type1s_ensure_scaffold_columns(reference$table)
  paths <- list(
    trd1_fasta = base::file.path(reference_dir, "type1s_trd1.faa"),
    trd2_fasta = base::file.path(reference_dir, "type1s_trd2.faa"),
    full_fasta = base::file.path(reference_dir, "type1s_full.faa"),
    scaffold_fasta = base::file.path(reference_dir, "type1s_spacer_scaffold.faa"),
    trd1_metadata = base::file.path(reference_dir, "type1s_trd1_metadata.tsv"),
    trd2_metadata = base::file.path(reference_dir, "type1s_trd2_metadata.tsv"),
    scaffold_metadata = base::file.path(reference_dir, "type1s_spacer_scaffold_metadata.tsv")
  )
  if (!.dnmb_type1s_nonempty_file(paths$trd1_fasta)) {
    .dnmb_type1s_write_fasta(
      base::paste0(table$reference_id, "_T1"),
      table$trd1_sequence,
      paths$trd1_fasta
    )
  }
  if (!.dnmb_type1s_nonempty_file(paths$trd2_fasta)) {
    .dnmb_type1s_write_fasta(
      base::paste0(table$reference_id, "_T2"),
      table$trd2_sequence,
      paths$trd2_fasta
    )
  }
  if (!.dnmb_type1s_nonempty_file(paths$full_fasta)) {
    .dnmb_type1s_write_fasta(table$reference_id, table$sequence, paths$full_fasta)
  }
  if (!.dnmb_type1s_nonempty_file(paths$scaffold_fasta)) {
    .dnmb_type1s_write_fasta(
      table$reference_id,
      table$spacer_scaffold_sequence,
      paths$scaffold_fasta
    )
  }

  write_metadata <- function(position, path) {
    metadata_columns <- c(
      "subject_id", "reference_id", "source_enzyme", "source_position",
      "position_half_site", "canonical_half_site", "cross_position_half_site",
      "spacer_length", "spacer_set", "spacer_ambiguous", "recognition",
      "protein_length", "window_start", "window_end", "boundary_method",
      "evidence", "trd_sequence"
    )
    existing_valid <- FALSE
    if (.dnmb_type1s_nonempty_file(path)) {
      existing <- base::tryCatch(
        utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
      )
      existing_valid <- base::is.data.frame(existing) &&
        base::nrow(existing) == base::nrow(table) &&
        base::all(metadata_columns %in% base::names(existing))
    }
    if (existing_valid) return(invisible(path))
    is_trd1 <- position == "T1"
    half_site <- if (is_trd1) table$left_half else table$right_half
    ambiguity_column <- if (is_trd1) "trd1_spacer_ambiguous" else "trd2_spacer_ambiguous"
    spacer_set_column <- if (is_trd1) "trd1_spacer_set" else "trd2_spacer_set"
    spacer_ambiguous <- if (ambiguity_column %in% base::names(table)) {
      table[[ambiguity_column]]
    } else {
      base::rep(FALSE, base::nrow(table))
    }
    spacer_set <- if (spacer_set_column %in% base::names(table)) {
      table[[spacer_set_column]]
    } else {
      base::as.character(table$spacer_length)
    }
    protein_length <- base::nchar(table$sequence)
    window_start <- if (is_trd1) {
      base::round(protein_length * 0.11) + 1L
    } else {
      base::round(protein_length * 0.66) + 1L
    }
    window_end <- if (is_trd1) {
      base::round(protein_length * 0.30)
    } else {
      base::round(protein_length * 0.82)
    }
    metadata <- base::data.frame(
      subject_id = base::paste0(table$reference_id, "_", position),
      reference_id = table$reference_id,
      source_enzyme = table$source_enzyme,
      source_position = position,
      position_half_site = half_site,
      canonical_half_site = if (is_trd1) half_site else .dnmb_type1s_reverse_complement(half_site),
      cross_position_half_site = .dnmb_type1s_reverse_complement(half_site),
      spacer_length = table$spacer_length,
      spacer_set = spacer_set,
      spacer_ambiguous = spacer_ambiguous,
      recognition = table$recognition,
      protein_length = protein_length,
      window_start = window_start,
      window_end = window_end,
      boundary_method = "gold_cv_specificity_core_window_v1",
      evidence = "REBASE_Type_I_S_Gold_Standard",
      trd_sequence = if (is_trd1) table$trd1_sequence else table$trd2_sequence,
      stringsAsFactors = FALSE
    )
    temporary <- base::tempfile("type1s-metadata-", tmpdir = base::dirname(path))
    on.exit(base::unlink(temporary, force = TRUE), add = TRUE)
    utils::write.table(
      metadata,
      file = temporary,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      na = ""
    )
    if (!.dnmb_type1s_nonempty_file(temporary)) {
      base::stop("Could not write Type I-S metadata.", call. = FALSE)
    }
    if (base::file.exists(path)) base::unlink(path, force = TRUE)
    if (!base::file.rename(temporary, path)) {
      copied <- base::file.copy(temporary, path, overwrite = TRUE)
      if (!base::isTRUE(copied)) base::stop("Could not install Type I-S metadata.", call. = FALSE)
    }
    invisible(path)
  }
  write_metadata("T1", paths$trd1_metadata)
  write_metadata("T2", paths$trd2_metadata)

  scaffold_columns <- c(
    "subject_id", "reference_id", "source_enzyme", "spacer_length",
    "recognition", "protein_length", "scaffold_length", "central_length",
    "helix_friendly_fraction", "helix_breaker_fraction",
    "tetrapeptide_repeat_count", "tael_like_repeat_count", "segment_1",
    "segment_2", "segment_3", "scaffold_method", "spacer_model_version",
    "evidence", "scaffold_sequence"
  )
  scaffold_valid <- FALSE
  if (.dnmb_type1s_nonempty_file(paths$scaffold_metadata)) {
    existing <- base::tryCatch(
      utils::read.delim(paths$scaffold_metadata, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    scaffold_valid <- base::is.data.frame(existing) &&
      base::nrow(existing) == base::nrow(table) &&
      base::all(scaffold_columns %in% base::names(existing)) &&
      base::isTRUE(base::all(
        existing$spacer_model_version == .dnmb_type1s_spacer_model_version()
      ))
  }
  if (!scaffold_valid) {
    protein_length <- base::nchar(table$sequence)
    first_end <- base::round(protein_length * 0.11)
    central_start <- base::round(protein_length * 0.30) + 1L
    central_end <- base::round(protein_length * 0.66)
    distal_start <- base::round(protein_length * 0.82) + 1L
    scaffold_metadata <- base::data.frame(
      subject_id = table$reference_id,
      reference_id = table$reference_id,
      source_enzyme = table$source_enzyme,
      spacer_length = table$spacer_length,
      recognition = table$recognition,
      protein_length = protein_length,
      scaffold_length = table$spacer_scaffold_length,
      central_length = table$spacer_central_length,
      helix_friendly_fraction = table$spacer_scaffold_helix_friendly_fraction,
      helix_breaker_fraction = table$spacer_scaffold_helix_breaker_fraction,
      tetrapeptide_repeat_count = table$spacer_tetrapeptide_repeat_count,
      tael_like_repeat_count = table$spacer_tael_like_repeat_count,
      segment_1 = base::paste0("1-", first_end),
      segment_2 = base::paste0(central_start, "-", central_end),
      segment_3 = base::paste0(distal_start, "-", protein_length),
      scaffold_method = "conserved_scaffold_fixed_windows_v1",
      spacer_model_version = .dnmb_type1s_spacer_model_version(),
      evidence = "REBASE_Type_I_S_Gold_Standard",
      scaffold_sequence = table$spacer_scaffold_sequence,
      stringsAsFactors = FALSE
    )
    temporary <- base::tempfile("type1s-scaffold-metadata-", tmpdir = reference_dir)
    on.exit(base::unlink(temporary, force = TRUE), add = TRUE)
    utils::write.table(
      scaffold_metadata,
      file = temporary,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      na = ""
    )
    if (!.dnmb_type1s_nonempty_file(temporary)) {
      base::stop("Could not write Type I-S scaffold metadata.", call. = FALSE)
    }
    if (base::file.exists(paths$scaffold_metadata)) {
      base::unlink(paths$scaffold_metadata, force = TRUE)
    }
    if (!base::file.rename(temporary, paths$scaffold_metadata)) {
      copied <- base::file.copy(temporary, paths$scaffold_metadata, overwrite = TRUE)
      if (!base::isTRUE(copied)) {
        base::stop("Could not install Type I-S scaffold metadata.", call. = FALSE)
      }
    }
  }
  paths
}

.dnmb_type1s_prepare_blast_db <- function(reference, reference_dir) {
  files <- .dnmb_type1s_prepare_reference_files(reference, reference_dir)
  trd1_prefix <- base::file.path(reference_dir, "type1s_trd1")
  trd2_prefix <- base::file.path(reference_dir, "type1s_trd2")
  full_prefix <- base::file.path(reference_dir, "type1s_full")
  scaffold_prefix <- base::file.path(reference_dir, "type1s_spacer_scaffold")

  build_one <- function(fasta, prefix) {
    required <- base::paste0(prefix, c(".phr", ".pin", ".psq"))
    if (base::all(base::vapply(required, .dnmb_type1s_nonempty_file, logical(1)))) {
      return(invisible(prefix))
    }
    stage <- base::tempfile("type1s-db-")
    base::dir.create(stage, recursive = TRUE, showWarnings = FALSE)
    on.exit(base::unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
    staged_fasta <- base::file.path(stage, "reference.faa")
    staged_prefix <- base::file.path(stage, "reference")
    base::file.copy(fasta, staged_fasta, overwrite = TRUE)
    run <- dnmb_run_external(
      "makeblastdb",
      args = c(
        "-in", staged_fasta,
        "-dbtype", "prot",
        "-blastdb_version", "4",
        "-out", staged_prefix
      ),
      required = FALSE
    )
    staged_required <- base::paste0(staged_prefix, c(".phr", ".pin", ".psq"))
    if (!base::isTRUE(run$ok) ||
        !base::all(base::vapply(staged_required, .dnmb_type1s_nonempty_file, logical(1)))) {
      base::stop(run$error %||% "makeblastdb failed for Type I-S reference.", call. = FALSE)
    }
    staged_files <- base::list.files(stage, pattern = "^reference\\.", full.names = TRUE)
    for (source in staged_files) {
      suffix <- base::substring(base::basename(source), base::nchar("reference") + 1L)
      copied <- base::file.copy(source, base::paste0(prefix, suffix), overwrite = TRUE)
      if (!base::isTRUE(copied)) {
        base::stop("Could not install Type I-S BLAST database artifact.", call. = FALSE)
      }
    }
    if (!base::all(base::vapply(required, .dnmb_type1s_nonempty_file, logical(1)))) {
      base::stop("Type I-S BLAST database was not created completely.", call. = FALSE)
    }
    invisible(prefix)
  }

  build_one(files$trd1_fasta, trd1_prefix)
  build_one(files$trd2_fasta, trd2_prefix)
  build_one(files$full_fasta, full_prefix)
  build_one(files$scaffold_fasta, scaffold_prefix)
  list(
    trd1 = trd1_prefix,
    trd2 = trd2_prefix,
    full = full_prefix,
    scaffold = scaffold_prefix
  )
}

.dnmb_type1s_prepare_diamond_db <- function(reference, reference_dir) {
  files <- .dnmb_type1s_prepare_reference_files(reference, reference_dir)
  trd1_prefix <- base::file.path(reference_dir, "type1s_trd1_diamond")
  trd2_prefix <- base::file.path(reference_dir, "type1s_trd2_diamond")
  full_prefix <- base::file.path(reference_dir, "type1s_full_diamond")
  scaffold_prefix <- base::file.path(reference_dir, "type1s_spacer_scaffold_diamond")
  build_one <- function(fasta, prefix) {
    database <- base::paste0(prefix, ".dmnd")
    if (.dnmb_type1s_nonempty_file(database)) return(invisible(prefix))
    staged_prefix <- base::tempfile("type1s-diamond-", tmpdir = reference_dir)
    staged_database <- base::paste0(staged_prefix, ".dmnd")
    on.exit(base::unlink(c(staged_prefix, staged_database), force = TRUE), add = TRUE)
    run <- dnmb_run_external(
      "diamond",
      args = c("makedb", "--in", fasta, "--db", staged_prefix),
      required = FALSE
    )
    if (!base::isTRUE(run$ok) || !.dnmb_type1s_nonempty_file(staged_database)) {
      base::stop(run$error %||% "diamond makedb failed for Type I-S reference.", call. = FALSE)
    }
    if (base::file.exists(database)) base::unlink(database, force = TRUE)
    if (!base::file.rename(staged_database, database)) {
      copied <- base::file.copy(staged_database, database, overwrite = TRUE)
      if (!base::isTRUE(copied)) base::stop("Could not install Type I-S DIAMOND database.", call. = FALSE)
    }
    invisible(prefix)
  }
  build_one(files$trd1_fasta, trd1_prefix)
  build_one(files$trd2_fasta, trd2_prefix)
  build_one(files$full_fasta, full_prefix)
  build_one(files$scaffold_fasta, scaffold_prefix)
  list(
    trd1 = trd1_prefix,
    trd2 = trd2_prefix,
    full = full_prefix,
    scaffold = scaffold_prefix
  )
}

.dnmb_type1s_write_manifest <- function(reference, files, databases) {
  artifacts <- base::data.frame(
    artifact = base::names(files),
    backend = "source",
    path = base::unlist(files, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  if (!base::is.null(databases$diamond)) {
    diamond <- base::paste0(base::unlist(databases$diamond, use.names = FALSE), ".dmnd")
    artifacts <- base::rbind(
      artifacts,
      base::data.frame(
        artifact = base::paste0(base::names(databases$diamond), "_database"),
        backend = "diamond",
        path = diamond,
        stringsAsFactors = FALSE
      )
    )
  }
  if (!base::is.null(databases$blastp)) {
    for (position in base::names(databases$blastp)) {
      prefix <- databases$blastp[[position]]
      candidates <- base::paste0(prefix, c(".phr", ".pin", ".psq"))
      artifacts <- base::rbind(
        artifacts,
        base::data.frame(
          artifact = base::paste0(position, "_database", c("_phr", "_pin", "_psq")),
          backend = "blastp",
          path = candidates,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  artifacts$path <- base::normalizePath(artifacts$path, winslash = "/", mustWork = FALSE)
  artifacts$exists <- base::file.exists(artifacts$path)
  artifacts$size <- NA_real_
  artifacts$md5 <- NA_character_
  if (base::any(artifacts$exists)) {
    info <- base::file.info(artifacts$path[artifacts$exists])
    artifacts$size[artifacts$exists] <- info$size
    artifacts$md5[artifacts$exists] <- base::unname(tools::md5sum(artifacts$path[artifacts$exists]))
  }
  artifacts$reference_version <- reference$version
  artifacts$source_url <- reference$source_url %||% .dnmb_type1s_gold_url()
  artifacts$database_version <- reference$database_version
  artifacts$spacer_model_version <- .dnmb_type1s_spacer_model_version()
  artifacts$source_md5 <- reference$source_md5
  # Keep the signed manifest content-derived. Wall-clock timestamps make an
  # otherwise identical database rebuild look like an artifact change to the
  # REBASEfinder stage cache. Build time remains available in reference.rds;
  # the manifest records the official source, versions, and checksums needed
  # to reproduce and audit the installed files.
  artifacts <- artifacts[, c(
    "reference_version", "database_version", "spacer_model_version",
    "source_url", "source_md5",
    "artifact", "backend", "path", "exists", "size", "md5"
  ), drop = FALSE]

  manifest <- base::file.path(reference$reference_dir, "type1s_database_manifest.tsv")
  temporary <- base::tempfile("type1s-manifest-", tmpdir = reference$reference_dir)
  on.exit(base::unlink(temporary, force = TRUE), add = TRUE)
  utils::write.table(
    artifacts,
    file = temporary,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  if (!.dnmb_type1s_nonempty_file(temporary)) {
    base::stop("Could not write the Type I-S database manifest.", call. = FALSE)
  }
  if (base::file.exists(manifest) &&
      base::identical(
        base::unname(tools::md5sum(temporary)),
        base::unname(tools::md5sum(manifest))
      )) {
    return(list(path = manifest, table = artifacts))
  }
  if (base::file.exists(manifest)) base::unlink(manifest, force = TRUE)
  if (!base::file.rename(temporary, manifest)) {
    copied <- base::file.copy(temporary, manifest, overwrite = TRUE)
    if (!base::isTRUE(copied)) base::stop("Could not install the Type I-S database manifest.", call. = FALSE)
  }
  list(path = manifest, table = artifacts)
}

#' Build separate Type I HsdS TRD specificity-core databases
#'
#' Builds position-specific TRD1 and TRD2 specificity-core protein-window
#' databases from the official REBASE Type I S Gold Standard. The windows are
#' validation-selected coordinates, not experimentally curated domain
#' boundaries. Each database has a companion metadata table that maps the
#' protein window and coordinates to its displayed/canonical half-site,
#' reverse-complemented cross-position half-site, spacer set/ambiguity, full
#' motif, evidence, and source enzyme. Full-HsdS and conserved structural-
#' scaffold databases are retained for spacer prediction. The scaffold joins
#' the three conserved regions outside the two specificity cores and records
#' central-region length, helix-friendly/breaker composition, and tandem
#' four-residue repeat diagnostics.
#'
#' @param cache_dir Optional REBASEfinder cache directory.
#' @param download Download the official Gold reference when absent.
#' @param force Re-download the reference and rebuild all database files.
#' @param backend Protein-search database format: `"auto"`, `"diamond"`,
#'   `"blastp"`, or `"both"`.
#' @param verbose Print progress messages.
#'
#' @return A list with database paths, reference provenance, and per-position
#'   sequence/label statistics, plus a checksummed artifact manifest.
#' @export
dnmb_build_type1s_trd_databases <- function(cache_dir = NULL,
                                             download = TRUE,
                                             force = FALSE,
                                             backend = c("auto", "diamond", "blastp", "both"),
                                             verbose = TRUE) {
  backend <- base::match.arg(backend)
  if (base::is.null(cache_dir) || !base::nzchar(base::as.character(cache_dir)[1])) {
    cache_dir <- base::file.path(.dnmb_db_cache_root(create = TRUE), "rebasefinder", "cache")
  }
  cache_dir <- base::path.expand(base::as.character(cache_dir)[1])
  base::dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  reference <- .dnmb_type1s_reference(
    cache_dir,
    force = base::isTRUE(force),
    download = base::isTRUE(download),
    verbose = verbose
  )
  files <- .dnmb_type1s_prepare_reference_files(reference, reference$reference_dir)

  diamond_found <- base::isTRUE(dnmb_detect_binary("diamond", required = FALSE)$found)
  blast_found <- base::isTRUE(dnmb_detect_binary("blastp", required = FALSE)$found) &&
    base::isTRUE(dnmb_detect_binary("makeblastdb", required = FALSE)$found)
  selected <- backend
  if (backend == "auto") {
    selected <- if (diamond_found) "diamond" else if (blast_found) "blastp" else "none"
  }
  if (selected %in% c("diamond", "both") && !diamond_found) {
    base::stop("DIAMOND is required for backend='", backend, "'.", call. = FALSE)
  }
  if (selected %in% c("blastp", "both") && !blast_found) {
    base::stop("blastp and makeblastdb are required for backend='", backend, "'.", call. = FALSE)
  }
  if (selected == "none") {
    base::stop("Neither DIAMOND nor BLASTP/makeblastdb is available.", call. = FALSE)
  }

  databases <- list()
  if (selected %in% c("diamond", "both")) {
    diamond_result <- base::tryCatch(
      .dnmb_type1s_prepare_diamond_db(reference, reference$reference_dir),
      error = function(e) e
    )
    if (base::inherits(diamond_result, "error")) {
      if (backend == "auto" && blast_found) {
        if (base::isTRUE(verbose)) {
          message(
            "[Type I-S] DIAMOND database build failed; falling back to BLASTP: ",
            base::conditionMessage(diamond_result)
          )
        }
        selected <- "blastp"
      } else {
        base::stop(base::conditionMessage(diamond_result), call. = FALSE)
      }
    } else {
      databases$diamond <- diamond_result
    }
  }
  if (selected %in% c("blastp", "both")) {
    databases$blastp <- .dnmb_type1s_prepare_blast_db(reference, reference$reference_dir)
  }
  manifest <- .dnmb_type1s_write_manifest(reference, files, databases)

  summarize_position <- function(sequence, half_site, spacer_length) {
    half_groups <- base::split(half_site, sequence)
    spacer_groups <- base::split(spacer_length, sequence)
    half_conflicts <- base::vapply(
      half_groups,
      function(labels) base::length(base::unique(labels)) > 1L,
      logical(1)
    )
    spacer_conflicts <- base::vapply(
      spacer_groups,
      function(spacers) base::length(base::unique(spacers)) > 1L,
      logical(1)
    )
    base::list(
      records = base::length(sequence),
      unique_sequences = base::length(half_groups),
      conflicting_half_site_sequences = base::sum(half_conflicts),
      conflicting_spacer_sequences = base::sum(spacer_conflicts)
    )
  }
  table <- .dnmb_type1s_ensure_scaffold_columns(reference$table)
  result <- list(
    reference_version = reference$version,
    database_version = reference$database_version,
    spacer_model_version = .dnmb_type1s_spacer_model_version(),
    source_path = reference$source_path,
    source_url = reference$source_url %||% .dnmb_type1s_gold_url(),
    source_md5 = reference$source_md5,
    reference_records = base::nrow(table),
    statistics = list(
      trd1 = summarize_position(table$trd1_sequence, table$left_half, table$spacer_length),
      trd2 = summarize_position(table$trd2_sequence, table$right_half, table$spacer_length),
      spacer_scaffold = base::list(
        records = base::nrow(table),
        unique_sequences = base::length(base::unique(table$spacer_scaffold_sequence)),
        median_length = stats::median(table$spacer_scaffold_length),
        spacer_model_version = .dnmb_type1s_spacer_model_version()
      )
    ),
    files = base::lapply(files, base::normalizePath, winslash = "/", mustWork = FALSE),
    databases = databases,
    manifest_path = base::normalizePath(manifest$path, winslash = "/", mustWork = FALSE),
    manifest = manifest$table
  )
  base::class(result) <- c("dnmb_type1s_databases", "list")
  result
}

.dnmb_type1s_run_blast <- function(query_fasta, db_prefix, output_path, cpu = 1L) {
  run <- dnmb_run_external(
    "blastp",
    args = c(
      "-query", query_fasta,
      "-db", db_prefix,
      "-out", output_path,
      "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
      "-evalue", "1e-3",
      "-max_target_seqs", "30",
      "-num_threads", base::as.character(base::max(1L, base::as.integer(cpu)[1]))
    ),
    required = FALSE
  )
  if (!base::isTRUE(run$ok)) {
    base::stop(run$error %||% "blastp failed for Type I-S prediction.", call. = FALSE)
  }
  if (!base::file.exists(output_path) || base::file.info(output_path)$size == 0L) {
    return(base::data.frame())
  }
  hits <- utils::read.table(output_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  base::names(hits) <- c(
    "query_id", "subject_id", "identity", "alignment_length", "query_length",
    "subject_length", "query_start", "query_end", "subject_start", "subject_end",
    "evalue", "bitscore"
  )
  hits$query_coverage <- 100 * hits$alignment_length / hits$query_length
  hits$subject_coverage <- 100 * hits$alignment_length / hits$subject_length
  hits <- hits[
    hits$identity >= 25 & hits$alignment_length >= 35 &
      hits$query_coverage >= 30 & hits$subject_coverage >= 30 & hits$evalue <= 1e-3,
    , drop = FALSE
  ]
  hits[base::order(hits$query_id, -hits$bitscore, -hits$query_coverage, -hits$identity), , drop = FALSE]
}

.dnmb_type1s_run_diamond <- function(query_fasta, db_prefix, output_path, cpu = 1L) {
  run <- dnmb_run_external(
    "diamond",
    args = c(
      "blastp",
      "--query", query_fasta,
      "--db", db_prefix,
      "--out", output_path,
      "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
      "--evalue", "1e-3",
      "--max-target-seqs", "30",
      "--threads", base::as.character(base::max(1L, base::as.integer(cpu)[1])),
      "--sensitive"
    ),
    required = FALSE
  )
  if (!base::isTRUE(run$ok)) {
    base::stop(run$error %||% "DIAMOND failed for Type I-S prediction.", call. = FALSE)
  }
  if (!base::file.exists(output_path) || base::file.info(output_path)$size == 0L) {
    return(base::data.frame())
  }
  hits <- utils::read.table(output_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  base::names(hits) <- c(
    "query_id", "subject_id", "identity", "alignment_length", "query_length",
    "subject_length", "query_start", "query_end", "subject_start", "subject_end",
    "evalue", "bitscore"
  )
  hits$query_coverage <- 100 * hits$alignment_length / hits$query_length
  hits$subject_coverage <- 100 * hits$alignment_length / hits$subject_length
  hits <- hits[
    hits$identity >= 25 & hits$alignment_length >= 35 &
      hits$query_coverage >= 30 & hits$subject_coverage >= 30 & hits$evalue <= 1e-3,
    , drop = FALSE
  ]
  hits[base::order(hits$query_id, -hits$bitscore, -hits$query_coverage, -hits$identity), , drop = FALSE]
}

.dnmb_type1s_empty_prediction <- function(ids) {
  ids <- base::as.character(ids)
  n <- base::length(ids)
  base::data.frame(
    locus_tag = ids,
    type1s_prediction_status = base::rep("not_run", n),
    type1s_left_half = base::rep(NA_character_, n),
    type1s_right_half = base::rep(NA_character_, n),
    type1s_spacer_length = base::rep(NA_integer_, n),
    type1s_predicted_recognition = base::rep(NA_character_, n),
    type1s_halfsite_confidence = base::rep(NA_character_, n),
    type1s_spacer_confidence = base::rep(NA_character_, n),
    type1s_overall_confidence = base::rep(NA_character_, n),
    type1s_prediction_eligible = base::rep(FALSE, n),
    type1s_trd1_identity = base::rep(NA_real_, n),
    type1s_trd1_coverage = base::rep(NA_real_, n),
    type1s_trd1_subject_coverage = base::rep(NA_real_, n),
    type1s_trd1_source = base::rep(NA_character_, n),
    type1s_trd1_source_position = base::rep(NA_character_, n),
    type1s_trd1_spacer_ambiguous = base::rep(NA, n),
    type1s_trd2_identity = base::rep(NA_real_, n),
    type1s_trd2_coverage = base::rep(NA_real_, n),
    type1s_trd2_subject_coverage = base::rep(NA_real_, n),
    type1s_trd2_source = base::rep(NA_character_, n),
    type1s_trd2_source_position = base::rep(NA_character_, n),
    type1s_trd2_spacer_ambiguous = base::rep(NA, n),
    type1s_spacer_source = base::rep(NA_character_, n),
    type1s_spacer_source_identity = base::rep(NA_real_, n),
    type1s_spacer_source_coverage = base::rep(NA_real_, n),
    type1s_spacer_source_subject_coverage = base::rep(NA_real_, n),
    type1s_spacer_method = base::rep(NA_character_, n),
    type1s_spacer_neighbor_count = base::rep(NA_integer_, n),
    type1s_spacer_supporter_count = base::rep(NA_integer_, n),
    type1s_spacer_voter_ids = base::rep(NA_character_, n),
    type1s_spacer_supporter_ids = base::rep(NA_character_, n),
    type1s_spacer_vote_support = base::rep(NA_real_, n),
    type1s_spacer_vote_margin = base::rep(NA_real_, n),
    type1s_spacer_full_vote = base::rep(NA_integer_, n),
    type1s_spacer_scaffold_vote = base::rep(NA_integer_, n),
    type1s_spacer_scaffold_vote_support = base::rep(NA_real_, n),
    type1s_spacer_scaffold_vote_margin = base::rep(NA_real_, n),
    type1s_spacer_model_agreement = base::rep(NA, n),
    type1s_spacer_model_version = base::rep(.dnmb_type1s_spacer_model_version(), n),
    type1s_spacer_tael_ruler_applied = base::rep(FALSE, n),
    type1s_spacer_tael_base_spacer = base::rep(NA_integer_, n),
    type1s_spacer_tael_query_repeat_count = base::rep(NA_integer_, n),
    type1s_spacer_tael_reference_repeat_count = base::rep(NA_integer_, n),
    type1s_spacer_tael_reference_id = base::rep(NA_character_, n),
    type1s_spacer_tael_reference_identity = base::rep(NA_real_, n),
    type1s_spacer_tael_reference_reciprocal_coverage = base::rep(NA_real_, n),
    type1s_spacer_tael_base_method = base::rep(NA_character_, n),
    type1s_spacer_scaffold_length = base::rep(NA_integer_, n),
    type1s_spacer_central_length = base::rep(NA_integer_, n),
    type1s_spacer_scaffold_helix_friendly_fraction = base::rep(NA_real_, n),
    type1s_spacer_scaffold_helix_breaker_fraction = base::rep(NA_real_, n),
    type1s_spacer_tetrapeptide_repeat_count = base::rep(NA_integer_, n),
    type1s_spacer_tael_like_repeat_count = base::rep(NA_integer_, n),
    type1s_search_backend = base::rep(NA_character_, n),
    type1s_boundary_method = base::rep(NA_character_, n),
    type1s_prediction_version = base::rep(NA_character_, n),
    type1s_reference_version = base::rep(NA_character_, n),
    type1s_note = base::rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

.dnmb_type1s_select_trd_hit <- function(hits, query_id, query_position, reference_table) {
  query_hits <- hits[hits$query_id == query_id, , drop = FALSE]
  if (!base::nrow(query_hits)) return(NULL)
  query_hits$subject_position <- base::sub("^.*_(T[12])$", "\\1", query_hits$subject_id)
  query_hits$reference_id <- base::sub("_(T[12])$", "", query_hits$subject_id)
  same <- query_hits[query_hits$subject_position == query_position, , drop = FALSE]
  cross <- query_hits[query_hits$subject_position != query_position, , drop = FALSE]
  same <- if (base::nrow(same)) same[1L, , drop = FALSE] else NULL
  cross <- if (base::nrow(cross)) cross[1L, , drop = FALSE] else NULL

  # Scores from the two position-specific reference windows are not on exactly
  # the same scale. These priors were selected on the frozen Gold validation
  # split, separately for N- and C-terminal TRDs.
  cross_factor <- if (query_position == "T1") 0.85 else 0.975
  chosen <- same
  if (base::is.null(chosen) || (!base::is.null(cross) && cross$bitscore[[1]] > chosen$bitscore[[1]] * cross_factor)) {
    chosen <- cross
  }
  if (base::is.null(chosen)) return(NULL)
  ref_index <- base::match(chosen$reference_id[[1]], reference_table$reference_id)
  if (base::is.na(ref_index)) return(NULL)
  cross_position <- chosen$subject_position[[1]] != query_position
  spacer_ambiguity_column <- if (chosen$subject_position[[1]] == "T1") {
    "trd1_spacer_ambiguous"
  } else {
    "trd2_spacer_ambiguous"
  }
  spacer_set_column <- if (chosen$subject_position[[1]] == "T1") {
    "trd1_spacer_set"
  } else {
    "trd2_spacer_set"
  }
  spacer_ambiguous <- if (spacer_ambiguity_column %in% base::names(reference_table)) {
    base::isTRUE(reference_table[[spacer_ambiguity_column]][[ref_index]])
  } else {
    FALSE
  }
  spacer_set <- if (spacer_set_column %in% base::names(reference_table)) {
    base::as.character(reference_table[[spacer_set_column]][[ref_index]])
  } else {
    base::as.character(reference_table$spacer_length[[ref_index]])
  }
  half_site <- if (query_position == "T1") {
    if (cross_position) .dnmb_type1s_reverse_complement(reference_table$right_half[[ref_index]]) else reference_table$left_half[[ref_index]]
  } else {
    if (cross_position) .dnmb_type1s_reverse_complement(reference_table$left_half[[ref_index]]) else reference_table$right_half[[ref_index]]
  }
  list(
    half_site = half_site,
    identity = chosen$identity[[1]],
    coverage = chosen$query_coverage[[1]],
    subject_coverage = chosen$subject_coverage[[1]],
    source = reference_table$source_enzyme[[ref_index]],
    source_position = chosen$subject_position[[1]],
    spacer = reference_table$spacer_length[[ref_index]],
    spacer_set = spacer_set,
    spacer_ambiguous = spacer_ambiguous,
    reference_id = reference_table$reference_id[[ref_index]],
    bitscore = chosen$bitscore[[1]]
  )
}

.dnmb_type1s_hit_is_high <- function(hit) {
  if (base::is.null(hit)) return(FALSE)
  subject <- hit$subject_coverage
  if (base::is.null(subject) || base::is.na(subject)) subject <- hit$coverage
  base::isTRUE(hit$identity >= 90) &&
    base::isTRUE(base::min(hit$coverage, subject) >= 80)
}

.dnmb_type1s_choose_spacer <- function(left_hit,
                                       right_hit,
                                       full_spacer_hit = NULL,
                                       weighted_spacer_hit = NULL) {
  consistent_trd_sources <- !base::is.null(left_hit) && !base::is.null(right_hit) &&
    !base::isTRUE(left_hit$spacer_ambiguous) &&
    !base::isTRUE(right_hit$spacer_ambiguous) &&
    base::identical(base::as.integer(left_hit$spacer), base::as.integer(right_hit$spacer))
  close_trd_sources <- consistent_trd_sources &&
    .dnmb_type1s_hit_is_high(left_hit) &&
    .dnmb_type1s_hit_is_high(right_hit)
  consensus <- function() {
    reference_ids <- base::paste(
      base::unique(c(left_hit$reference_id, right_hit$reference_id)),
      collapse = ","
    )
    list(
      spacer = base::as.integer(left_hit$spacer),
      identity = base::min(left_hit$identity, right_hit$identity),
      coverage = base::min(left_hit$coverage, right_hit$coverage),
      subject_coverage = base::min(left_hit$subject_coverage, right_hit$subject_coverage),
      source = base::paste(base::unique(c(left_hit$source, right_hit$source)), collapse = "+"),
      method = "trd_source_consensus",
      neighbor_count = 2L,
      supporter_count = 2L,
      voter_ids = reference_ids,
      supporter_ids = reference_ids,
      vote_support = 1,
      vote_margin = 1,
      can_be_high = TRUE
    )
  }

  # Spacer is strongly influenced by conserved HsdS context, so remote TRD
  # labels must not override whole-protein evidence. Preserve either calibrated
  # close-homolog transfer, then use the development-selected full/scaffold
  # ensemble. Remote TRD consensus is only a last-resort exploratory fallback.
  if (close_trd_sources) return(consensus())
  if (.dnmb_type1s_hit_is_high(left_hit) &&
      .dnmb_type1s_hit_is_high(right_hit) &&
      .dnmb_type1s_hit_is_high(full_spacer_hit)) {
    return(full_spacer_hit)
  }
  if (!base::is.null(weighted_spacer_hit)) return(weighted_spacer_hit)
  if (!base::is.null(full_spacer_hit)) return(full_spacer_hit)
  if (consistent_trd_sources) return(consensus())
  NULL
}

.dnmb_type1s_select_spacer_1nn <- function(full_hits, reference_table) {
  if (base::is.null(full_hits) || !base::is.data.frame(full_hits) || !base::nrow(full_hits)) {
    return(NULL)
  }
  ref_index <- base::match(full_hits$subject_id, reference_table$reference_id)
  valid <- base::which(
    !base::is.na(ref_index) &
      !base::is.na(full_hits$identity) &
      !base::is.na(full_hits$query_coverage)
  )
  if (!base::length(valid)) return(NULL)
  row <- valid[[1]]
  r <- ref_index[[row]]
  subject_coverage <- full_hits$subject_coverage[[row]]
  if (base::is.null(subject_coverage) || base::is.na(subject_coverage)) {
    subject_coverage <- full_hits$query_coverage[[row]]
  }
  list(
    spacer = base::as.integer(reference_table$spacer_length[[r]]),
    identity = full_hits$identity[[row]],
    coverage = full_hits$query_coverage[[row]],
    subject_coverage = subject_coverage,
    source = reference_table$source_enzyme[[r]],
    reference_id = reference_table$reference_id[[r]],
    tael_like_repeat_count = if (
      "spacer_tael_like_repeat_count" %in% base::names(reference_table)
    ) {
      base::as.integer(reference_table$spacer_tael_like_repeat_count[[r]])
    } else {
      NA_integer_
    },
    method = "whole_hsds_1nn",
    neighbor_count = 1L,
    supporter_count = 1L,
    voter_ids = reference_table$reference_id[[r]],
    supporter_ids = reference_table$reference_id[[r]],
    vote_support = 1,
    vote_margin = 1,
    can_be_high = TRUE
  )
}

.dnmb_type1s_canonical_tael_spacer <- function(repeat_count) {
  repeat_count <- base::as.integer(repeat_count)[1]
  if (base::is.na(repeat_count)) return(NA_integer_)
  if (repeat_count == 2L) return(6L)
  if (repeat_count == 3L) return(7L)
  NA_integer_
}

.dnmb_type1s_apply_tael_ruler <- function(hit,
                                          nearest_full_hit,
                                          query_repeat_count,
                                          mode = c("close", "remote"),
                                          remote_margin_max = 0.40) {
  mode <- base::match.arg(mode)
  if (base::is.null(hit) || base::is.null(nearest_full_hit)) return(hit)

  query_repeat_count <- base::as.integer(query_repeat_count)[1]
  reference_repeat_count <- base::as.integer(
    nearest_full_hit$tael_like_repeat_count %||% NA_integer_
  )[1]
  reference_spacer <- base::as.integer(nearest_full_hit$spacer %||% NA_integer_)[1]
  query_spacer <- .dnmb_type1s_canonical_tael_spacer(query_repeat_count)
  canonical_reference_spacer <- .dnmb_type1s_canonical_tael_spacer(
    reference_repeat_count
  )
  canonical_pair <- !base::is.na(query_spacer) &&
    !base::is.na(canonical_reference_spacer) &&
    base::abs(query_repeat_count - reference_repeat_count) == 1L &&
    base::identical(reference_spacer, canonical_reference_spacer)
  if (!canonical_pair) return(hit)

  adjusted_spacer <- reference_spacer + query_repeat_count - reference_repeat_count
  if (!base::identical(base::as.integer(adjusted_spacer), query_spacer)) return(hit)
  if (base::identical(base::as.integer(hit$spacer), query_spacer)) return(hit)
  reciprocal_coverage <- base::min(
    nearest_full_hit$coverage %||% NA_real_,
    nearest_full_hit$subject_coverage %||% nearest_full_hit$coverage %||% NA_real_
  )
  eligible <- if (mode == "close") {
    base::identical(hit$method %||% "", "whole_hsds_1nn") &&
      .dnmb_type1s_hit_is_high(nearest_full_hit)
  } else {
    margin <- base::as.numeric(hit$vote_margin %||% NA_real_)[1]
    base::identical(hit$method %||% "", "whole_hsds_scaffold_ensemble") &&
      base::isTRUE(nearest_full_hit$identity >= 70) &&
      base::isTRUE(reciprocal_coverage >= 80) &&
      !base::is.na(margin) && margin <= remote_margin_max
  }
  if (!base::isTRUE(eligible)) return(hit)

  hit$tael_ruler_applied <- TRUE
  hit$tael_base_spacer <- base::as.integer(hit$spacer)
  hit$tael_query_repeat_count <- query_repeat_count
  hit$tael_reference_repeat_count <- reference_repeat_count
  hit$tael_reference_id <- nearest_full_hit$reference_id %||% NA_character_
  hit$tael_reference_identity <- nearest_full_hit$identity %||% NA_real_
  hit$tael_reference_reciprocal_coverage <- reciprocal_coverage
  hit$tael_base_method <- hit$method %||% NA_character_
  hit$spacer <- query_spacer
  hit$method <- if (mode == "close") {
    "whole_hsds_1nn_tael_ruler"
  } else {
    "whole_hsds_scaffold_tael_ruler"
  }
  hit$can_be_high <- base::identical(mode, "close")
  hit
}

.dnmb_type1s_spacer_vote <- function(hits, reference_table, k = 5L) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(NULL)
  }
  k <- base::max(1L, base::as.integer(k)[1])
  ref_index <- base::match(hits$subject_id, reference_table$reference_id)
  valid <- !base::is.na(ref_index) &
    !base::is.na(hits$identity) &
    !base::is.na(hits$query_coverage)
  if (!base::any(valid)) return(NULL)
  hits <- hits[valid, , drop = FALSE]
  ref_index <- ref_index[valid]
  unique_rows <- !base::duplicated(hits$subject_id)
  hits <- hits[unique_rows, , drop = FALSE]
  ref_index <- ref_index[unique_rows]
  take <- base::seq_len(base::min(k, base::nrow(hits)))
  hits <- hits[take, , drop = FALSE]
  ref_index <- ref_index[take]

  subject_coverage <- hits$subject_coverage
  if (base::is.null(subject_coverage)) subject_coverage <- hits$query_coverage
  subject_coverage[base::is.na(subject_coverage)] <- hits$query_coverage[base::is.na(subject_coverage)]
  reciprocal_coverage <- base::pmin(hits$query_coverage, subject_coverage)
  weights <- hits$identity * reciprocal_coverage
  spacers <- base::as.integer(reference_table$spacer_length[ref_index])
  valid_vote <- !base::is.na(spacers) & !base::is.na(weights) & weights > 0
  if (!base::any(valid_vote)) return(NULL)
  hits <- hits[valid_vote, , drop = FALSE]
  ref_index <- ref_index[valid_vote]
  subject_coverage <- subject_coverage[valid_vote]
  reciprocal_coverage <- reciprocal_coverage[valid_vote]
  weights <- weights[valid_vote]
  spacers <- spacers[valid_vote]

  vote_totals <- base::tapply(weights, spacers, base::sum)
  winning_total <- base::max(vote_totals)
  tied_spacers <- base::as.integer(base::names(vote_totals)[vote_totals == winning_total])
  winner <- spacers[base::which(spacers %in% tied_spacers)[[1]]]
  supporting <- base::which(spacers == winner)
  nearest_supporter <- supporting[[1]]
  ordered_totals <- base::sort(base::as.numeric(vote_totals), decreasing = TRUE)
  total_vote <- base::sum(ordered_totals)
  runner_up <- if (base::length(ordered_totals) >= 2L) ordered_totals[[2]] else 0

  list(
    hits = hits,
    ref_index = ref_index,
    subject_coverage = subject_coverage,
    reciprocal_coverage = reciprocal_coverage,
    weights = weights,
    spacers = spacers,
    vote_totals = vote_totals,
    normalized_vote = vote_totals / total_vote,
    winner = winner,
    supporting = supporting,
    nearest_supporter = nearest_supporter,
    vote_support = winning_total / total_vote,
    vote_margin = (winning_total - runner_up) / total_vote
  )
}

.dnmb_type1s_select_spacer_knn <- function(full_hits, reference_table, k = 5L) {
  vote <- .dnmb_type1s_spacer_vote(full_hits, reference_table, k = k)
  if (base::is.null(vote)) return(NULL)
  nearest_supporter <- vote$nearest_supporter
  supporting <- vote$supporting

  list(
    spacer = vote$winner,
    identity = vote$hits$identity[[nearest_supporter]],
    coverage = vote$hits$query_coverage[[nearest_supporter]],
    subject_coverage = vote$subject_coverage[[nearest_supporter]],
    source = reference_table$source_enzyme[[vote$ref_index[[nearest_supporter]]]],
    method = base::paste0("whole_hsds_weighted_knn", k),
    neighbor_count = base::length(vote$spacers),
    supporter_count = base::length(supporting),
    voter_ids = base::paste(reference_table$reference_id[vote$ref_index], collapse = ","),
    supporter_ids = base::paste(
      reference_table$reference_id[vote$ref_index[supporting]],
      collapse = ","
    ),
    vote_support = vote$vote_support,
    vote_margin = vote$vote_margin,
    full_vote = vote$winner,
    can_be_high = FALSE
  )
}

.dnmb_type1s_select_spacer_ensemble <- function(full_hits,
                                                scaffold_hits,
                                                reference_table,
                                                full_k = 5L,
                                                scaffold_k = 20L,
                                                full_weight = 0.70) {
  full_vote <- .dnmb_type1s_spacer_vote(full_hits, reference_table, k = full_k)
  scaffold_vote <- .dnmb_type1s_spacer_vote(
    scaffold_hits,
    reference_table,
    k = scaffold_k
  )
  if (base::is.null(scaffold_vote)) {
    return(.dnmb_type1s_select_spacer_knn(full_hits, reference_table, k = full_k))
  }
  full_weight <- base::min(1, base::max(0, base::as.numeric(full_weight)[1]))
  if (base::is.null(full_vote)) full_weight <- 0

  spacer_classes <- base::sort(base::unique(base::c(
    if (!base::is.null(full_vote)) base::as.integer(base::names(full_vote$normalized_vote)) else integer(),
    base::as.integer(base::names(scaffold_vote$normalized_vote))
  )))
  full_scores <- stats::setNames(base::rep(0, base::length(spacer_classes)), spacer_classes)
  scaffold_scores <- full_scores
  if (!base::is.null(full_vote)) {
    full_scores[base::names(full_vote$normalized_vote)] <- full_vote$normalized_vote
  }
  scaffold_scores[base::names(scaffold_vote$normalized_vote)] <- scaffold_vote$normalized_vote
  ensemble_scores <- full_weight * full_scores + (1 - full_weight) * scaffold_scores
  winning_total <- base::max(ensemble_scores)
  tied <- base::as.integer(base::names(ensemble_scores)[ensemble_scores == winning_total])
  winner <- if (!base::is.null(full_vote) && full_vote$winner %in% tied) {
    full_vote$winner
  } else {
    tied[[1]]
  }
  ordered_totals <- base::sort(base::as.numeric(ensemble_scores), decreasing = TRUE)
  runner_up <- if (base::length(ordered_totals) >= 2L) ordered_totals[[2]] else 0

  full_supporting <- if (base::is.null(full_vote)) integer() else {
    base::which(full_vote$spacers == winner)
  }
  scaffold_supporting <- base::which(scaffold_vote$spacers == winner)
  use_full_provenance <- base::length(full_supporting) > 0L
  provenance_vote <- if (use_full_provenance) full_vote else scaffold_vote
  provenance_row <- if (use_full_provenance) full_supporting[[1]] else scaffold_supporting[[1]]
  full_voters <- if (base::is.null(full_vote)) character() else {
    base::paste0("F:", reference_table$reference_id[full_vote$ref_index])
  }
  scaffold_voters <- base::paste0(
    "S:",
    reference_table$reference_id[scaffold_vote$ref_index]
  )
  supporter_ids <- base::unique(base::c(
    if (base::length(full_supporting)) {
      base::paste0("F:", reference_table$reference_id[full_vote$ref_index[full_supporting]])
    } else {
      character()
    },
    if (base::length(scaffold_supporting)) {
      base::paste0(
        "S:",
        reference_table$reference_id[scaffold_vote$ref_index[scaffold_supporting]]
      )
    } else {
      character()
    }
  ))

  list(
    spacer = base::as.integer(winner),
    identity = provenance_vote$hits$identity[[provenance_row]],
    coverage = provenance_vote$hits$query_coverage[[provenance_row]],
    subject_coverage = provenance_vote$subject_coverage[[provenance_row]],
    source = reference_table$source_enzyme[[provenance_vote$ref_index[[provenance_row]]]],
    method = "whole_hsds_scaffold_ensemble",
    neighbor_count = base::length(full_voters) + base::length(scaffold_voters),
    supporter_count = base::length(supporter_ids),
    voter_ids = base::paste(base::c(full_voters, scaffold_voters), collapse = ","),
    supporter_ids = base::paste(supporter_ids, collapse = ","),
    vote_support = winning_total,
    vote_margin = winning_total - runner_up,
    full_vote = if (base::is.null(full_vote)) NA_integer_ else full_vote$winner,
    scaffold_vote = scaffold_vote$winner,
    scaffold_vote_support = scaffold_vote$vote_support,
    scaffold_vote_margin = scaffold_vote$vote_margin,
    model_agreement = !base::is.null(full_vote) &&
      base::identical(base::as.integer(full_vote$winner), base::as.integer(scaffold_vote$winner)),
    spacer_model_version = .dnmb_type1s_spacer_model_version(),
    provenance_branch = if (use_full_provenance) "full" else "scaffold",
    can_be_high = FALSE
  )
}

.dnmb_type1s_confidence <- function(left_hit, right_hit, spacer_hit, exact = FALSE) {
  if (base::isTRUE(exact)) {
    return(list(half = "known", spacer = "known", overall = "known", eligible = TRUE))
  }
  reciprocal_coverage <- function(hit) {
    if (base::is.null(hit)) return(NA_real_)
    subject_coverage <- hit$subject_coverage
    if (base::is.null(subject_coverage) || base::is.na(subject_coverage)) {
      subject_coverage <- hit$coverage
    }
    base::min(hit$coverage, subject_coverage)
  }
  both_halves <- !base::is.null(left_hit) && !base::is.null(right_hit)
  half_identity <- if (both_halves) base::min(left_hit$identity, right_hit$identity) else NA_real_
  half_coverage <- if (both_halves) {
    base::min(reciprocal_coverage(left_hit), reciprocal_coverage(right_hit))
  } else {
    NA_real_
  }
  half <- if (!both_halves) {
    "partial"
  } else if (half_coverage >= 80 && half_identity >= 90) {
    "high"
  } else if (half_coverage >= 70 && half_identity >= 70) {
    "medium"
  } else if (half_coverage >= 50 && half_identity >= 50) {
    "low"
  } else {
    "exploratory"
  }

  spacer_identity <- if (!base::is.null(spacer_hit)) spacer_hit$identity else NA_real_
  spacer_coverage <- reciprocal_coverage(spacer_hit)
  spacer_method <- if (!base::is.null(spacer_hit)) spacer_hit$method %||% "" else ""
  spacer_consensus <- base::identical(spacer_method, "trd_source_consensus")
  can_be_high <- !base::is.null(spacer_hit) &&
    spacer_method %in% c(
      "trd_source_consensus", "whole_hsds_1nn",
      "whole_hsds_1nn_tael_ruler"
    ) &&
    base::identical(spacer_hit$can_be_high, TRUE)
  spacer <- if (base::is.na(spacer_identity) || base::is.na(spacer_coverage)) {
    "unsupported"
  } else if (can_be_high && spacer_coverage >= 80 && spacer_identity >= 90) {
    "high"
  } else if (spacer_consensus && spacer_coverage >= 70 && spacer_identity >= 70) {
    "medium"
  } else if (spacer_consensus) {
    "low"
  } else if (spacer_coverage >= 70 && spacer_identity >= 70) {
    "medium"
  } else if (spacer_coverage >= 70 && spacer_identity >= 50) {
    "low"
  } else {
    "exploratory"
  }

  overall <- if (!both_halves || base::is.null(spacer_hit)) {
    "partial"
  } else if (can_be_high && half_coverage >= 80 && half_identity >= 90 &&
             spacer_coverage >= 80 && spacer_identity >= 90) {
    "high"
  } else if (half_coverage >= 70 && half_identity >= 70 && spacer_coverage >= 70 && spacer_identity >= 70) {
    "medium"
  } else {
    "low"
  }
  list(half = half, spacer = spacer, overall = overall, eligible = overall %in% c("known", "high"))
}

.dnmb_type1s_predict_core <- function(candidates,
                                      reference,
                                      cache_dir,
                                      cpu = 1L,
                                      verbose = TRUE,
                                      backend = c("auto", "diamond", "blastp")) {
  backend <- base::match.arg(backend)
  candidates <- base::as.data.frame(candidates, stringsAsFactors = FALSE)
  output <- .dnmb_type1s_empty_prediction(candidates$locus_tag)
  output$type1s_prediction_version <- .dnmb_type1s_prediction_version()
  output$type1s_reference_version <- reference$version
  reference_table <- reference$table
  sequences <- .dnmb_type1s_clean_protein(candidates$translation)
  windows <- base::lapply(sequences, .dnmb_type1s_extract_windows)

  for (i in base::seq_len(base::nrow(output))) {
    if (base::is.null(windows[[i]])) {
      output$type1s_prediction_status[[i]] <- "abstained_architecture"
      output$type1s_note[[i]] <- "Canonical full-length Type I HsdS requires 300-700 aa; partial/circular variants are not predicted."
    } else {
      output$type1s_boundary_method[[i]] <- windows[[i]]$boundary_method
      output$type1s_spacer_scaffold_length[[i]] <- windows[[i]]$spacer_scaffold_length
      output$type1s_spacer_central_length[[i]] <- windows[[i]]$spacer_central_length
      output$type1s_spacer_scaffold_helix_friendly_fraction[[i]] <-
        windows[[i]]$spacer_scaffold_helix_friendly_fraction
      output$type1s_spacer_scaffold_helix_breaker_fraction[[i]] <-
        windows[[i]]$spacer_scaffold_helix_breaker_fraction
      output$type1s_spacer_tetrapeptide_repeat_count[[i]] <-
        windows[[i]]$spacer_tetrapeptide_repeat_count
      output$type1s_spacer_tael_like_repeat_count[[i]] <-
        windows[[i]]$spacer_tael_like_repeat_count
    }
  }

  exact_index <- base::match(sequences, reference_table$sequence)
  exact_rows <- base::which(!base::is.na(exact_index) & !base::vapply(windows, base::is.null, logical(1)))
  for (i in exact_rows) {
    r <- exact_index[[i]]
    confidence <- .dnmb_type1s_confidence(NULL, NULL, NULL, exact = TRUE)
    output$type1s_prediction_status[[i]] <- "known_gold_match"
    output$type1s_left_half[[i]] <- reference_table$left_half[[r]]
    output$type1s_right_half[[i]] <- reference_table$right_half[[r]]
    output$type1s_spacer_length[[i]] <- reference_table$spacer_length[[r]]
    output$type1s_predicted_recognition[[i]] <- reference_table$recognition[[r]]
    output$type1s_halfsite_confidence[[i]] <- confidence$half
    output$type1s_spacer_confidence[[i]] <- confidence$spacer
    output$type1s_overall_confidence[[i]] <- confidence$overall
    output$type1s_prediction_eligible[[i]] <- confidence$eligible
    output$type1s_trd1_identity[[i]] <- 100
    output$type1s_trd1_coverage[[i]] <- 100
    output$type1s_trd1_subject_coverage[[i]] <- 100
    output$type1s_trd1_source[[i]] <- reference_table$source_enzyme[[r]]
    output$type1s_trd1_source_position[[i]] <- "T1"
    output$type1s_trd1_spacer_ambiguous[[i]] <- if ("trd1_spacer_ambiguous" %in% base::names(reference_table)) {
      base::isTRUE(reference_table$trd1_spacer_ambiguous[[r]])
    } else {
      FALSE
    }
    output$type1s_trd2_identity[[i]] <- 100
    output$type1s_trd2_coverage[[i]] <- 100
    output$type1s_trd2_subject_coverage[[i]] <- 100
    output$type1s_trd2_source[[i]] <- reference_table$source_enzyme[[r]]
    output$type1s_trd2_source_position[[i]] <- "T2"
    output$type1s_trd2_spacer_ambiguous[[i]] <- if ("trd2_spacer_ambiguous" %in% base::names(reference_table)) {
      base::isTRUE(reference_table$trd2_spacer_ambiguous[[r]])
    } else {
      FALSE
    }
    output$type1s_spacer_source[[i]] <- reference_table$source_enzyme[[r]]
    output$type1s_spacer_source_identity[[i]] <- 100
    output$type1s_spacer_source_coverage[[i]] <- 100
    output$type1s_spacer_source_subject_coverage[[i]] <- 100
    output$type1s_spacer_method[[i]] <- "exact_gold"
    output$type1s_spacer_neighbor_count[[i]] <- 1L
    output$type1s_spacer_supporter_count[[i]] <- 1L
    output$type1s_spacer_voter_ids[[i]] <- reference_table$reference_id[[r]]
    output$type1s_spacer_supporter_ids[[i]] <- reference_table$reference_id[[r]]
    output$type1s_spacer_vote_support[[i]] <- 1
    output$type1s_spacer_vote_margin[[i]] <- 1
    output$type1s_spacer_full_vote[[i]] <- reference_table$spacer_length[[r]]
    output$type1s_spacer_scaffold_vote[[i]] <- reference_table$spacer_length[[r]]
    output$type1s_spacer_scaffold_vote_support[[i]] <- 1
    output$type1s_spacer_scaffold_vote_margin[[i]] <- 1
    output$type1s_spacer_model_agreement[[i]] <- TRUE
    output$type1s_note[[i]] <- "Exact protein match to official REBASE Type I S Gold Standard."
    output$type1s_search_backend[[i]] <- "exact_match"
  }

  search_rows <- base::which(
    base::is.na(exact_index) & !base::vapply(windows, base::is.null, logical(1))
  )
  if (!base::length(search_rows)) return(output)
  if (base::isTRUE(verbose)) message("[Type I-S] Predicting ", base::length(search_rows), " non-identical HsdS sequence(s)...")
  search_failed <- function(message, backend) {
    output$type1s_prediction_status[search_rows] <- "search_failed"
    output$type1s_search_backend[search_rows] <- backend
    output$type1s_note[search_rows] <- base::paste0("Protein search failed: ", message)
    output
  }

  diamond_available <- base::isTRUE(dnmb_detect_binary("diamond", required = FALSE)$found)
  blast_available <- base::isTRUE(dnmb_detect_binary("blastp", required = FALSE)$found) &&
    base::isTRUE(dnmb_detect_binary("makeblastdb", required = FALSE)$found)
  search_backend <- if (backend == "auto") {
    if (diamond_available) "diamond" else if (blast_available) "blastp" else "none"
  } else {
    backend
  }
  if (search_backend == "none" ||
      (search_backend == "diamond" && !diamond_available) ||
      (search_backend == "blastp" && !blast_available)) {
    detail <- if (backend == "diamond") {
      "Requested backend 'diamond' is unavailable."
    } else if (backend == "blastp") {
      "Requested backend 'blastp' requires both blastp and makeblastdb."
    } else {
      "Neither DIAMOND nor BLASTP/makeblastdb is available."
    }
    return(search_failed(detail, search_backend))
  }
  prepare_backend <- function(backend) {
    if (backend == "diamond") {
      return(list(
        db = .dnmb_type1s_prepare_diamond_db(reference, reference$reference_dir),
        run = .dnmb_type1s_run_diamond
      ))
    }
    list(
      db = .dnmb_type1s_prepare_blast_db(reference, reference$reference_dir),
      run = .dnmb_type1s_run_blast
    )
  }
  allow_blast_fallback <- backend == "auto" && blast_available
  setup <- base::tryCatch(prepare_backend(search_backend), error = function(e) e)
  if (base::inherits(setup, "error")) {
    if (search_backend == "diamond" && allow_blast_fallback) {
      if (base::isTRUE(verbose)) {
        message("[Type I-S] DIAMOND setup failed; falling back to BLASTP: ", base::conditionMessage(setup))
      }
      search_backend <- "blastp"
      setup <- base::tryCatch(prepare_backend(search_backend), error = function(e) e)
      if (base::inherits(setup, "error")) {
        if (base::isTRUE(verbose)) message("[Type I-S] BLASTP setup failed: ", base::conditionMessage(setup))
        return(search_failed(base::conditionMessage(setup), search_backend))
      }
    } else {
      if (base::isTRUE(verbose)) message("[Type I-S] Protein search setup failed: ", base::conditionMessage(setup))
      return(search_failed(base::conditionMessage(setup), search_backend))
    }
  }
  work_dir <- base::tempfile("dnmb-type1s-")
  base::dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(base::unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  query_ids <- base::sprintf("Q%04d", base::seq_along(search_rows))
  trd_query <- base::file.path(work_dir, "query_trd.faa")
  full_query <- base::file.path(work_dir, "query_full.faa")
  scaffold_query <- base::file.path(work_dir, "query_spacer_scaffold.faa")
  .dnmb_type1s_write_fasta(
    base::c(base::paste0(query_ids, "_T1"), base::paste0(query_ids, "_T2")),
    base::c(
      base::vapply(windows[search_rows], `[[`, character(1), "trd1"),
      base::vapply(windows[search_rows], `[[`, character(1), "trd2")
    ),
    trd_query
  )
  .dnmb_type1s_write_fasta(query_ids, sequences[search_rows], full_query)
  .dnmb_type1s_write_fasta(
    query_ids,
    base::vapply(windows[search_rows], `[[`, character(1), "spacer_scaffold"),
    scaffold_query
  )
  run_pair <- function(search_setup) {
    trd1_hits <- search_setup$run(
      trd_query,
      search_setup$db$trd1,
      base::file.path(work_dir, "trd1_database_hits.tsv"),
      cpu = cpu
    )
    trd2_hits <- search_setup$run(
      trd_query,
      search_setup$db$trd2,
      base::file.path(work_dir, "trd2_database_hits.tsv"),
      cpu = cpu
    )
    trd_hits <- dplyr::bind_rows(trd1_hits, trd2_hits)
    if (base::nrow(trd_hits)) {
      trd_hits <- trd_hits[
        base::order(
          trd_hits$query_id,
          -trd_hits$bitscore,
          -trd_hits$query_coverage,
          -trd_hits$identity
        ),
        , drop = FALSE
      ]
    }
    scaffold_hits <- if (!base::is.null(search_setup$db$scaffold)) {
      search_setup$run(
        scaffold_query,
        search_setup$db$scaffold,
        base::file.path(work_dir, "spacer_scaffold_hits.tsv"),
        cpu = cpu
      )
    } else {
      base::data.frame()
    }
    list(
      trd = trd_hits,
      full = search_setup$run(
        full_query,
        search_setup$db$full,
        base::file.path(work_dir, "full_hits.tsv"),
        cpu = cpu
      ),
      scaffold = scaffold_hits
    )
  }
  search_result <- base::tryCatch(run_pair(setup), error = function(e) e)
  if (base::inherits(search_result, "error") && search_backend == "diamond" && allow_blast_fallback) {
    if (base::isTRUE(verbose)) {
      message("[Type I-S] DIAMOND search failed; falling back to BLASTP: ", base::conditionMessage(search_result))
    }
    search_backend <- "blastp"
    setup <- base::tryCatch(prepare_backend(search_backend), error = function(e) e)
    search_result <- if (base::inherits(setup, "error")) {
      setup
    } else {
      base::tryCatch(run_pair(setup), error = function(e) e)
    }
  }
  if (base::inherits(search_result, "error")) {
    if (base::isTRUE(verbose)) message("[Type I-S] Protein search failed: ", base::conditionMessage(search_result))
    return(search_failed(base::conditionMessage(search_result), search_backend))
  }
  trd_hits <- search_result$trd
  full_hits <- search_result$full
  scaffold_hits <- search_result$scaffold

  for (j in base::seq_along(search_rows)) {
    row <- search_rows[[j]]
    query_id <- query_ids[[j]]
    left_hit <- .dnmb_type1s_select_trd_hit(trd_hits, base::paste0(query_id, "_T1"), "T1", reference_table)
    right_hit <- .dnmb_type1s_select_trd_hit(trd_hits, base::paste0(query_id, "_T2"), "T2", reference_table)

    candidate_full <- full_hits[full_hits$query_id == query_id, , drop = FALSE]
    candidate_scaffold <- scaffold_hits[
      scaffold_hits$query_id == query_id,
      , drop = FALSE
    ]
    full_spacer_hit <- .dnmb_type1s_select_spacer_1nn(
      candidate_full,
      reference_table
    )
    nearest_full_hit <- full_spacer_hit
    weighted_spacer_hit <- .dnmb_type1s_select_spacer_ensemble(
      candidate_full,
      candidate_scaffold,
      reference_table,
      full_k = 5L,
      scaffold_k = 20L,
      full_weight = 0.70
    )
    full_spacer_hit <- .dnmb_type1s_apply_tael_ruler(
      full_spacer_hit,
      full_spacer_hit,
      windows[[row]]$spacer_tael_like_repeat_count,
      mode = "close"
    )
    weighted_spacer_hit <- .dnmb_type1s_apply_tael_ruler(
      weighted_spacer_hit,
      nearest_full_hit,
      windows[[row]]$spacer_tael_like_repeat_count,
      mode = "remote"
    )
    spacer_hit <- .dnmb_type1s_choose_spacer(
      left_hit,
      right_hit,
      full_spacer_hit,
      weighted_spacer_hit
    )
    confidence <- .dnmb_type1s_confidence(left_hit, right_hit, spacer_hit)
    output$type1s_left_half[[row]] <- if (!base::is.null(left_hit)) left_hit$half_site else NA_character_
    output$type1s_right_half[[row]] <- if (!base::is.null(right_hit)) right_hit$half_site else NA_character_
    output$type1s_spacer_length[[row]] <- if (!base::is.null(spacer_hit)) spacer_hit$spacer else NA_integer_
    complete <- !base::is.null(left_hit) && !base::is.null(right_hit) && !base::is.null(spacer_hit)
    output$type1s_predicted_recognition[[row]] <- if (complete) {
      base::paste0(left_hit$half_site, base::strrep("N", spacer_hit$spacer), right_hit$half_site)
    } else {
      NA_character_
    }
    output$type1s_prediction_status[[row]] <- if (complete) "predicted_complete" else if (!base::is.null(left_hit) || !base::is.null(right_hit)) "predicted_partial" else "no_supported_hit"
    output$type1s_halfsite_confidence[[row]] <- confidence$half
    output$type1s_spacer_confidence[[row]] <- confidence$spacer
    output$type1s_overall_confidence[[row]] <- confidence$overall
    output$type1s_prediction_eligible[[row]] <- confidence$eligible
    if (!base::is.null(left_hit)) {
      output$type1s_trd1_identity[[row]] <- left_hit$identity
      output$type1s_trd1_coverage[[row]] <- left_hit$coverage
      output$type1s_trd1_subject_coverage[[row]] <- left_hit$subject_coverage
      output$type1s_trd1_source[[row]] <- left_hit$source
      output$type1s_trd1_source_position[[row]] <- left_hit$source_position
      output$type1s_trd1_spacer_ambiguous[[row]] <- base::isTRUE(left_hit$spacer_ambiguous)
    }
    if (!base::is.null(right_hit)) {
      output$type1s_trd2_identity[[row]] <- right_hit$identity
      output$type1s_trd2_coverage[[row]] <- right_hit$coverage
      output$type1s_trd2_subject_coverage[[row]] <- right_hit$subject_coverage
      output$type1s_trd2_source[[row]] <- right_hit$source
      output$type1s_trd2_source_position[[row]] <- right_hit$source_position
      output$type1s_trd2_spacer_ambiguous[[row]] <- base::isTRUE(right_hit$spacer_ambiguous)
    }
    if (!base::is.null(spacer_hit)) {
      output$type1s_spacer_source[[row]] <- spacer_hit$source
      output$type1s_spacer_source_identity[[row]] <- spacer_hit$identity
      output$type1s_spacer_source_coverage[[row]] <- spacer_hit$coverage
      output$type1s_spacer_source_subject_coverage[[row]] <- spacer_hit$subject_coverage
      output$type1s_spacer_method[[row]] <- spacer_hit$method
      output$type1s_spacer_neighbor_count[[row]] <- spacer_hit$neighbor_count %||% 2L
      output$type1s_spacer_supporter_count[[row]] <- spacer_hit$supporter_count %||% 2L
      output$type1s_spacer_voter_ids[[row]] <- spacer_hit$voter_ids %||% NA_character_
      output$type1s_spacer_supporter_ids[[row]] <- spacer_hit$supporter_ids %||% NA_character_
      output$type1s_spacer_vote_support[[row]] <- spacer_hit$vote_support %||% NA_real_
      output$type1s_spacer_vote_margin[[row]] <- spacer_hit$vote_margin %||% NA_real_
      output$type1s_spacer_full_vote[[row]] <- spacer_hit$full_vote %||% NA_integer_
      output$type1s_spacer_scaffold_vote[[row]] <- spacer_hit$scaffold_vote %||% NA_integer_
      output$type1s_spacer_scaffold_vote_support[[row]] <-
        spacer_hit$scaffold_vote_support %||% NA_real_
      output$type1s_spacer_scaffold_vote_margin[[row]] <-
        spacer_hit$scaffold_vote_margin %||% NA_real_
      output$type1s_spacer_model_agreement[[row]] <-
        spacer_hit$model_agreement %||% NA
      output$type1s_spacer_tael_ruler_applied[[row]] <-
        base::isTRUE(spacer_hit$tael_ruler_applied)
      output$type1s_spacer_tael_base_spacer[[row]] <-
        spacer_hit$tael_base_spacer %||% NA_integer_
      output$type1s_spacer_tael_query_repeat_count[[row]] <-
        spacer_hit$tael_query_repeat_count %||% NA_integer_
      output$type1s_spacer_tael_reference_repeat_count[[row]] <-
        spacer_hit$tael_reference_repeat_count %||% NA_integer_
      output$type1s_spacer_tael_reference_id[[row]] <-
        spacer_hit$tael_reference_id %||% NA_character_
      output$type1s_spacer_tael_reference_identity[[row]] <-
        spacer_hit$tael_reference_identity %||% NA_real_
      output$type1s_spacer_tael_reference_reciprocal_coverage[[row]] <-
        spacer_hit$tael_reference_reciprocal_coverage %||% NA_real_
      output$type1s_spacer_tael_base_method[[row]] <-
        spacer_hit$tael_base_method %||% NA_character_
    }
    output$type1s_note[[row]] <- if (
      !base::is.null(spacer_hit) && base::isTRUE(spacer_hit$tael_ruler_applied)
    ) {
      base::paste0(
        if (confidence$eligible) "High-confidence close-homolog transfer" else "Exploratory structural spacer prediction",
        " with conservative Type-IC TAEL molecular-ruler correction; inspect base ensemble votes and ruler provenance."
      )
    } else if (confidence$eligible) {
      "High-confidence close-homolog transfer (>=90% identity and >=80% reciprocal coverage) under official REBASE Gold calibration."
    } else if (confidence$half == "high") {
      "Half-sites are high confidence; spacer or complete motif remains below the high-confidence threshold."
    } else {
      "Exploratory prediction; inspect per-TRD identity/coverage and validate experimentally."
    }
    output$type1s_search_backend[[row]] <- search_backend
  }
  output
}

#' Predict Type I HsdS recognition sequences
#'
#' Predicts the two half-sites and intervening spacer for canonical full-length
#' Type I specificity (HsdS) proteins. The reference is the official REBASE
#' Type I S Gold Standard. Validation-selected TRD1 and TRD2 specificity-core
#' windows are stored in physically separate position-specific databases; they
#' are not claimed to be curated domain boundaries. Each query window searches
#' both databases so reverse-complement-aware cross-position reuse remains
#' possible. Spacer length preserves close-homolog transfer first. Remote
#' predictions combine a reciprocal-coverage-weighted five-neighbor full-HsdS
#' vote (70 percent) with a twenty-neighbor conserved-scaffold vote (30
#' percent). This sequence-derived structural proxy captures coiled-coil
#' length, helix-breaking residues, repeat phase, and conserved-family context.
#' A narrowly gated Type-IC TAEL-repeat molecular-ruler correction distinguishes
#' N6 from N7 only for canonical two-versus-three-repeat reference pairs. Remote
#' ensemble and ruler calls can never enter the high-confidence tier;
#' low-evidence results remain explicitly exploratory.
#'
#' @param candidates Data frame containing Type I HsdS candidates.
#' @param seq_col Protein-sequence column.
#' @param id_col Stable candidate identifier column.
#' @param cache_dir Optional REBASEfinder cache directory.
#' @param download Download the official Gold reference when absent.
#' @param force_reference Re-download and rebuild the reference.
#' @param cpu Number of protein-search threads.
#' @param backend Protein-search backend. `"auto"` prefers DIAMOND and falls
#'   back to BLASTP; explicit values make runs backend-reproducible.
#' @param verbose Print progress messages.
#'
#' @return Data frame with half-site, spacer, provenance, and calibrated
#'   confidence columns.
#' @export
dnmb_predict_type1s_recognition <- function(candidates,
                                            seq_col = "translation",
                                            id_col = "locus_tag",
                                            cache_dir = NULL,
                                            download = TRUE,
                                            force_reference = FALSE,
                                            cpu = 1L,
                                            backend = c("auto", "diamond", "blastp"),
                                            verbose = TRUE) {
  backend <- base::match.arg(backend)
  candidates <- base::as.data.frame(candidates, stringsAsFactors = FALSE)
  missing <- base::setdiff(c(id_col, seq_col), base::names(candidates))
  if (base::length(missing)) {
    base::stop("Type I-S candidates are missing columns: ", base::paste(missing, collapse = ", "), call. = FALSE)
  }
  input <- base::data.frame(
    locus_tag = .dnmb_type1s_clean_id(candidates[[id_col]]),
    translation = .dnmb_type1s_clean_protein(candidates[[seq_col]]),
    stringsAsFactors = FALSE
  )
  input <- input[!base::is.na(input$locus_tag) & !base::duplicated(input$locus_tag), , drop = FALSE]
  if (!base::nrow(input)) return(.dnmb_type1s_empty_prediction(character()))
  if (base::is.null(cache_dir) || !base::nzchar(base::as.character(cache_dir)[1])) {
    cache_dir <- base::file.path(.dnmb_db_cache_root(create = TRUE), "rebasefinder", "cache")
  }
  base::dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  reference <- .dnmb_type1s_reference(
    cache_dir,
    force = base::isTRUE(force_reference),
    download = base::isTRUE(download),
    verbose = verbose
  )
  .dnmb_type1s_predict_core(
    input,
    reference,
    cache_dir = cache_dir,
    cpu = cpu,
    verbose = verbose,
    backend = backend
  )
}

.dnmb_type1s_candidate_ids <- function(genes, rm_comprehensive = NULL) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  ids <- character()
  if (!base::is.null(rm_comprehensive) && base::is.data.frame(rm_comprehensive) && base::nrow(rm_comprehensive)) {
    rm <- base::as.data.frame(rm_comprehensive, stringsAsFactors = FALSE)
    id_col <- base::intersect(c("locus_tag", "protein_id"), base::names(rm))[1]
    if (!base::is.na(id_col) && base::all(c("rm_type", "subunit") %in% base::names(rm))) {
      type_norm <- base::gsub(
        "[^a-z0-9]",
        "",
        base::tolower(base::as.character(rm$rm_type))
      )
      subunit_norm <- base::tolower(base::trimws(base::as.character(rm$subunit)))
      keep <- type_norm %in% c("typei", "i") & subunit_norm %in% c("s", "specificity", "specificitysubunit")
      ids <- base::c(ids, base::as.character(rm[[id_col]][keep]))
    }
  }
  if (base::all(c("locus_tag", "product") %in% base::names(genes))) {
    product <- base::tolower(base::as.character(genes$product))
    product[base::is.na(product)] <- ""
    strong_hsds <- base::grepl("(^|[^a-z0-9])hsds([^a-z0-9]|$)", product, perl = TRUE) |
      base::grepl(
        "type[ -]?i([^a-z0-9]|$).*specificity|specificity.*type[ -]?i([^a-z0-9]|$)",
        product,
        perl = TRUE
      ) |
      base::grepl(
        "restriction(?:-modification)?(?: system| enzyme| endonuclease)?.*subunit[ -]?s([^a-z0-9]|$)",
        product,
        perl = TRUE
      )
    ids <- base::c(ids, base::as.character(genes$locus_tag[strong_hsds]))
  }
  ids <- .dnmb_type1s_clean_id(ids)
  base::unique(ids[!base::is.na(ids)])
}

.dnmb_type1s_candidate_rows <- function(gene_ids, candidate_ids) {
  gene_ids <- .dnmb_type1s_clean_id(gene_ids)
  candidate_ids <- base::unique(.dnmb_type1s_clean_id(candidate_ids))
  candidate_ids <- candidate_ids[!base::is.na(candidate_ids)]
  exact <- !base::is.na(gene_ids) & gene_ids %in% candidate_ids
  if (!base::length(candidate_ids) || base::all(exact)) return(exact)

  candidate_has_exact_gene <- candidate_ids %in% gene_ids
  fallback_ids <- candidate_ids[!candidate_has_exact_gene]
  if (!base::length(fallback_ids)) return(exact)

  exact_keys <- .dnmb_module_clean_annotation_key(candidate_ids[candidate_has_exact_gene])
  fallback_keys <- .dnmb_module_clean_annotation_key(fallback_ids)
  keep_fallback <- !base::is.na(fallback_keys) & !fallback_keys %in% exact_keys
  fallback_keys <- fallback_keys[keep_fallback]
  if (!base::length(fallback_keys)) return(exact)

  gene_keys <- .dnmb_module_clean_annotation_key(gene_ids)
  unique_gene_key <- !base::duplicated(gene_keys) & !base::duplicated(gene_keys, fromLast = TRUE)
  unique_candidate_key <- !base::duplicated(fallback_keys) &
    !base::duplicated(fallback_keys, fromLast = TRUE)
  fallback_match <- base::match(gene_keys, fallback_keys)
  use_fallback <- !base::is.na(fallback_match) & unique_gene_key &
    unique_candidate_key[fallback_match]
  exact | use_fallback
}
