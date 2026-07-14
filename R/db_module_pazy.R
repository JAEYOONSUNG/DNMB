.dnmb_pazy_module_name <- function() {
  "pazy"
}

.dnmb_pazy_default_version <- function() {
  "current"
}

.dnmb_pazy_default_metadata_url <- function() {
  "https://api.pazy.eu/api/proteins/?page_size=500"
}

.dnmb_pazy_default_fasta_url <- function() {
  "https://api.pazy.eu/api/proteins/fasta/"
}

.dnmb_pazy_fasta_contract_version <- function() {
  2L
}

.dnmb_pazy_makedb_timeout <- function() {
  # PAZy currently contains only a few hundred proteins. A database build that
  # runs for more than five minutes is wedged rather than merely slow.
  .dnmb_external_timeout_seconds(base::getOption("dnmb.pazy_makedb_timeout", 300))
}

.dnmb_pazy_search_timeout <- function() {
  .dnmb_external_timeout_seconds(base::getOption("dnmb.pazy_search_timeout", 1800))
}

.dnmb_pazy_network_timeout <- function() {
  .dnmb_external_timeout_seconds(base::getOption("dnmb.pazy_network_timeout", 300))
}

.dnmb_pazy_remote_asset_state <- function(metadata_url, fasta_url, enabled = TRUE) {
  if (!base::isTRUE(enabled)) return(NULL)
  urls <- c(metadata = metadata_url, fasta = fasta_url)
  checked_at <- base::format(base::Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  rows <- base::lapply(base::names(urls), function(asset_name) {
    metadata <- .dnmb_remote_asset_metadata(
      urls[[asset_name]],
      insecure = FALSE,
      timeout = .dnmb_pazy_network_timeout()
    )
    tibble::tibble(
      asset_name = asset_name,
      url = metadata$url %||% NA_character_,
      ok = base::isTRUE(metadata$ok),
      last_modified = metadata$last_modified %||% NA_character_,
      content_length = metadata$content_length %||% NA_real_,
      etag = metadata$etag %||% NA_character_,
      checked_at = checked_at
    )
  })
  dplyr::bind_rows(rows)
}

.dnmb_pazy_remote_update_needed <- function(manifest, remote_state, max_age_days = 30) {
  if (base::is.null(manifest) || base::is.null(remote_state) ||
      !base::is.data.frame(remote_state) || !base::nrow(remote_state) ||
      !base::all(remote_state$ok %in% TRUE)) {
    return(FALSE)
  }
  previous <- manifest$remote_asset_state
  if (base::is.null(previous) || !base::is.data.frame(previous) || !base::nrow(previous)) {
    return(TRUE)
  }
  fields <- c("asset_name", "last_modified", "content_length", "etag")
  previous <- previous[, base::intersect(fields, base::names(previous)), drop = FALSE]
  current <- remote_state[, fields, drop = FALSE]
  if (!base::all(fields %in% base::names(previous))) return(TRUE)
  previous <- previous[base::order(previous$asset_name), , drop = FALSE]
  current <- current[base::order(current$asset_name), , drop = FALSE]
  base::rownames(previous) <- NULL
  base::rownames(current) <- NULL
  if (!base::identical(previous, current)) return(TRUE)

  previous_state <- manifest$remote_asset_state
  checked_at <- if ("checked_at" %in% base::names(previous_state)) {
    previous_state$checked_at[[1]]
  } else {
    manifest$written_at %||% NA_character_
  }
  checked_time <- base::suppressWarnings(base::as.POSIXct(checked_at, tz = "UTC"))
  age_days <- base::as.numeric(base::difftime(base::Sys.time(), checked_time, units = "days"))
  base::is.na(age_days) || age_days > base::as.numeric(max_age_days)[1]
}

.dnmb_pazy_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = as.character(component)[1],
    status = as.character(status)[1],
    detail = as.character(detail)[1]
  )
}

.dnmb_pazy_empty_status <- function() {
  tibble::tibble(
    component = character(),
    status = character(),
    detail = character()
  )
}

.dnmb_pazy_trace <- function(path, text) {
  base::cat(text, "\n", file = path, append = TRUE)
}

.dnmb_pazy_asset_layout <- function(module_dir) {
  reference_fasta <- base::file.path(module_dir, "pazy_reference.faa")
  metadata_json <- base::file.path(module_dir, "pazy_metadata.json")
  metadata_tsv <- base::file.path(module_dir, "pazy_metadata.tsv")
  blast_db_prefix <- base::file.path(module_dir, "pazy_reference")
  list(
    module_dir = module_dir,
    reference_fasta = reference_fasta,
    metadata_json = metadata_json,
    metadata_tsv = metadata_tsv,
    blast_db_prefix = blast_db_prefix,
    blast_db_files = .dnmb_pazy_find_blast_db_files(blast_db_prefix)
  )
}

.dnmb_pazy_find_blast_db_files <- function(blast_db_prefix) {
  # V4: .phr, .pin, .psq
  # V5: adds .pdb, .pot, .ptf, .pto
  all_exts <- c(".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto")
  candidates <- base::paste0(blast_db_prefix, all_exts)
  found <- candidates[base::file.exists(candidates)]
  if (!base::length(found)) {
    # Fallback to V4 set even if not yet created (install step)
    return(base::paste0(blast_db_prefix, c(".phr", ".pin", ".psq")))
  }
  found
}

.dnmb_pazy_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for PAZy must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_pazy_fetch_json_text <- function(url) {
  run <- dnmb_run_external(
    "curl",
    args = c("-L", "-sS", url),
    required = FALSE,
    timeout = .dnmb_pazy_network_timeout()
  )
  if (!base::isTRUE(run$ok) || !base::length(run$stdout)) {
    run <- dnmb_run_external(
      "wget",
      args = c("-qO-", url),
      required = FALSE,
      timeout = .dnmb_pazy_network_timeout()
    )
  }
  if (!base::isTRUE(run$ok) || !base::length(run$stdout)) {
    base::stop(run$error %||% base::paste("Failed to fetch PAZy JSON from", url), call. = FALSE)
  }
  base::paste(run$stdout, collapse = "\n")
}

.dnmb_pazy_fetch_metadata_pages <- function(url = .dnmb_pazy_default_metadata_url()) {
  pages <- list()
  current_url <- as.character(url)[1]
  while (!base::is.na(current_url) && base::nzchar(current_url)) {
    payload <- jsonlite::fromJSON(.dnmb_pazy_fetch_json_text(current_url), simplifyVector = TRUE)
    pages[[base::length(pages) + 1L]] <- payload$results
    next_page_url <- payload[["next"]]
    current_url <- if (base::is.null(next_page_url) || !base::nzchar(next_page_url)) "" else base::as.character(next_page_url)[1]
  }
  dplyr::bind_rows(pages)
}

.dnmb_pazy_metadata_scalar <- function(row, field, default) {
  if (!field %in% base::names(row)) return(default)
  value <- row[[field]]
  if (base::is.list(value)) value <- value[[1]]
  if (base::is.null(value) || !base::length(value)) return(default)
  value[[1]]
}

.dnmb_pazy_flatten_metadata <- function(metadata_tbl) {
  if (base::is.null(metadata_tbl) || !base::is.data.frame(metadata_tbl) || !base::nrow(metadata_tbl)) {
    return(data.frame())
  }

  rows <- lapply(base::seq_len(base::nrow(metadata_tbl)), function(i) {
    row <- metadata_tbl[i, , drop = FALSE]
    substrates <- row$substrates[[1]]
    accessions <- row$accessions[[1]]
    organism <- row$organism[[1]]
    literature <- row$literature[[1]]
    amino_acid_sequence <- base::as.character(
      .dnmb_pazy_metadata_scalar(row, "amino_acid_sequence", NA_character_)
    )
    sequence_length <- base::suppressWarnings(base::as.integer(
      .dnmb_pazy_metadata_scalar(row, "sequence_length", NA_integer_)
    ))
    sequence_md5 <- base::as.character(
      .dnmb_pazy_metadata_scalar(row, "sequence_md5", NA_character_)
    )
    tibble::tibble(
      pazy_id = base::paste0("PAZY_", base::as.character(row$id[[1]])),
      pazy_name = base::as.character(row$name[[1]]),
      verified = base::isTRUE(row$verified[[1]]),
      substrate_names = if (base::is.data.frame(substrates) && base::nrow(substrates)) base::paste(substrates$name, collapse = "; ") else NA_character_,
      substrate_abbreviations = if (base::is.data.frame(substrates) && base::nrow(substrates)) base::paste(substrates$abbreviation, collapse = "; ") else NA_character_,
      organism_name = if (base::is.list(organism)) base::as.character(organism$scientific_name %||% NA_character_) else NA_character_,
      domain = if (base::is.list(organism)) base::as.character(organism$domain %||% NA_character_) else NA_character_,
      phylum = if (base::is.list(organism)) base::as.character(organism$phylum %||% NA_character_) else NA_character_,
      accession_list = if (base::is.data.frame(accessions) && base::nrow(accessions)) base::paste(accessions$accession, collapse = "; ") else NA_character_,
      accession_db_list = if (base::is.data.frame(accessions) && base::nrow(accessions)) base::paste(accessions$database, collapse = "; ") else NA_character_,
      literature_doi = if (base::is.data.frame(literature) && base::nrow(literature)) base::paste(literature$doi, collapse = "; ") else NA_character_,
      amino_acid_sequence = amino_acid_sequence,
      sequence_length = sequence_length,
      sequence_md5 = sequence_md5
    )
  })

  dplyr::bind_rows(rows)
}

.dnmb_pazy_normalize_record_id <- function(source_id) {
  source_id <- base::trimws(base::as.character(source_id)[1])
  if (base::is.na(source_id) || !base::nzchar(source_id)) {
    return(NA_character_)
  }
  normalized <- if (base::startsWith(source_id, "PAZY_")) {
    source_id
  } else {
    base::paste0("PAZY_", source_id)
  }
  normalized <- base::gsub("[^A-Za-z0-9_.:-]", "_", normalized)
  if (base::identical(normalized, "PAZY_") || !base::nzchar(normalized)) NA_character_ else normalized
}

.dnmb_pazy_clean_residues <- function(sequence) {
  sequence <- base::as.character(sequence)
  sequence[base::is.na(sequence)] <- ""
  raw <- base::paste(sequence, collapse = "")
  without_space <- base::gsub("[[:space:]]", "", raw)
  upper <- base::toupper(without_space)
  terminal_stop <- base::grepl("[*]+$", upper)
  cleaned <- base::sub("[*]+$", "", upper)
  valid <- base::nzchar(cleaned) &&
    base::grepl("^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$", cleaned)
  list(
    sequence = cleaned,
    valid = base::isTRUE(valid),
    whitespace_fixed = !base::identical(raw, without_space),
    whitespace_removed = base::nchar(raw, type = "chars") - base::nchar(without_space, type = "chars"),
    case_fixed = !base::identical(without_space, upper),
    terminal_stop_fixed = base::isTRUE(terminal_stop),
    invalid_symbols = if (base::isTRUE(valid) || !base::nzchar(cleaned)) {
      character()
    } else {
      base::unique(base::strsplit(base::gsub("[ACDEFGHIKLMNPQRSTVWYBXZJUO]", "", cleaned), "", fixed = TRUE)[[1]])
    }
  )
}

.dnmb_pazy_sequence_md5 <- function(sequence) {
  sequence <- base::as.character(sequence)[1]
  if (base::is.na(sequence)) return(NA_character_)
  path <- base::tempfile("dnmb-pazy-sequence-")
  base::on.exit(base::unlink(path, force = TRUE), add = TRUE)
  base::writeBin(base::charToRaw(sequence), path)
  base::unname(tools::md5sum(path)[[1]])
}

.dnmb_pazy_metadata_sequences <- function(metadata_flat) {
  empty <- base::data.frame(
    pazy_id = character(),
    pazy_name = character(),
    sequence = character(),
    stringsAsFactors = FALSE
  )
  if (base::is.null(metadata_flat) || !base::is.data.frame(metadata_flat) ||
      !base::all(c("pazy_id", "amino_acid_sequence") %in% base::names(metadata_flat))) {
    return(empty)
  }
  ids <- base::as.character(metadata_flat$pazy_id)
  names <- if ("pazy_name" %in% base::names(metadata_flat)) {
    base::as.character(metadata_flat$pazy_name)
  } else {
    ids
  }
  raw_sequences <- base::as.character(metadata_flat$amino_acid_sequence)
  normalized <- base::lapply(raw_sequences, .dnmb_pazy_clean_residues)
  normalized_sequences <- base::vapply(normalized, `[[`, character(1), "sequence")
  valid <- base::vapply(normalized, function(x) base::isTRUE(x$valid), logical(1))
  if ("sequence_length" %in% base::names(metadata_flat)) {
    expected_length <- base::suppressWarnings(base::as.integer(metadata_flat$sequence_length))
    valid <- valid & (base::is.na(expected_length) | base::nchar(normalized_sequences) == expected_length)
  }
  if ("sequence_md5" %in% base::names(metadata_flat)) {
    expected_md5 <- base::tolower(base::trimws(base::as.character(metadata_flat$sequence_md5)))
    verify <- base::which(valid & !base::is.na(expected_md5) & base::nzchar(expected_md5))
    if (base::length(verify)) {
      observed_md5 <- base::vapply(
        normalized_sequences[verify],
        .dnmb_pazy_sequence_md5,
        character(1)
      )
      valid[verify] <- observed_md5 == expected_md5[verify]
    }
  }
  if (!base::any(valid)) return(empty)
  out <- base::data.frame(
    pazy_id = ids[valid],
    pazy_name = names[valid],
    sequence = normalized_sequences[valid],
    stringsAsFactors = FALSE
  )
  out <- out[!base::is.na(out$pazy_id) & base::nzchar(out$pazy_id), , drop = FALSE]
  out <- out[!base::duplicated(out$pazy_id), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_pazy_fasta_audit_summary <- function(audit) {
  fields <- c(
    "ok", "clean", "n_input_records", "n_output_records",
    "n_whitespace_fixed", "n_whitespace_removed", "n_case_fixed", "n_terminal_stop_fixed",
    "n_empty_dropped", "n_invalid_dropped", "n_metadata_recovered",
    "n_metadata_added", "n_exact_duplicates_dropped",
    "n_conflicting_duplicates", "n_empty_headers", "n_orphan_lines",
    "n_ids_normalized", "corrected_ids", "dropped_ids", "detail"
  )
  audit[base::intersect(fields, base::names(audit))]
}

.dnmb_pazy_fasta_audit_detail <- function(audit) {
  base::sprintf(
    paste0(
      "records=%d/%d; whitespace_fixed=%d(%d chars); empty_dropped=%d; ",
      "invalid_dropped=%d; metadata_recovered=%d; duplicates_dropped=%d"
    ),
    audit$n_output_records,
    audit$n_input_records,
    audit$n_whitespace_fixed,
    audit$n_whitespace_removed,
    audit$n_empty_dropped,
    audit$n_invalid_dropped,
    audit$n_metadata_recovered + audit$n_metadata_added,
    audit$n_exact_duplicates_dropped
  )
}

.dnmb_pazy_sanitize_fasta <- function(source_fasta, dest_fasta = NULL, metadata_flat = NULL) {
  source_fasta <- base::path.expand(base::as.character(source_fasta)[1])
  if (base::is.na(source_fasta) || !base::nzchar(source_fasta) || !base::file.exists(source_fasta)) {
    return(list(ok = FALSE, clean = FALSE, detail = base::paste0("PAZy FASTA not found: ", source_fasta)))
  }

  lines <- base::readLines(source_fasta, warn = FALSE)
  header_index <- base::which(base::startsWith(lines, ">"))
  orphan_index <- if (base::length(header_index)) {
    base::seq_len(header_index[[1]] - 1L)
  } else {
    base::seq_along(lines)
  }
  orphan_lines <- if (base::length(orphan_index)) {
    lines[orphan_index][base::nzchar(base::trimws(lines[orphan_index]))]
  } else {
    character()
  }

  metadata_sequences <- .dnmb_pazy_metadata_sequences(metadata_flat)
  metadata_names <- if (!base::is.null(metadata_flat) && base::is.data.frame(metadata_flat) &&
                        base::all(c("pazy_id", "pazy_name") %in% base::names(metadata_flat))) {
    stats::setNames(base::as.character(metadata_flat$pazy_name), base::as.character(metadata_flat$pazy_id))
  } else {
    character()
  }
  canonical_sequences <- stats::setNames(metadata_sequences$sequence, metadata_sequences$pazy_id)

  records <- list()
  record_ids <- character()
  corrected_ids <- character()
  dropped_ids <- character()
  whitespace_fixed <- 0L
  whitespace_removed <- 0L
  case_fixed <- 0L
  terminal_stop_fixed <- 0L
  empty_dropped <- 0L
  invalid_dropped <- 0L
  metadata_recovered <- 0L
  metadata_added <- 0L
  exact_duplicates <- 0L
  conflicting_duplicates <- 0L
  conflicting_ids <- character()
  empty_headers <- 0L
  ids_normalized <- 0L

  if (base::length(header_index)) {
    record_end <- base::c(header_index[-1L] - 1L, base::length(lines))
    for (i in base::seq_along(header_index)) {
      header <- base::trimws(base::sub("^>", "", lines[header_index[[i]]]))
      if (!base::nzchar(header)) {
        empty_headers <- empty_headers + 1L
        next
      }
      source_id <- base::strsplit(header, "[[:space:]]+", perl = TRUE)[[1]][[1]]
      current_id <- .dnmb_pazy_normalize_record_id(source_id)
      if (base::is.na(current_id) || !base::nzchar(current_id)) {
        empty_headers <- empty_headers + 1L
        next
      }
      if (!base::identical(source_id, current_id)) ids_normalized <- ids_normalized + 1L

      sequence_lines <- if (record_end[[i]] > header_index[[i]]) {
        lines[(header_index[[i]] + 1L):record_end[[i]]]
      } else {
        character()
      }
      cleaned <- .dnmb_pazy_clean_residues(sequence_lines)
      whitespace_fixed <- whitespace_fixed + base::as.integer(cleaned$whitespace_fixed)
      whitespace_removed <- whitespace_removed + base::as.integer(cleaned$whitespace_removed)
      case_fixed <- case_fixed + base::as.integer(cleaned$case_fixed)
      terminal_stop_fixed <- terminal_stop_fixed + base::as.integer(cleaned$terminal_stop_fixed)
      if (cleaned$whitespace_fixed || cleaned$case_fixed || cleaned$terminal_stop_fixed) {
        corrected_ids <- base::c(corrected_ids, current_id)
      }

      canonical <- base::unname(canonical_sequences[current_id])
      has_canonical <- !base::is.na(canonical) && base::nzchar(canonical)
      if (!base::isTRUE(cleaned$valid)) {
        if (has_canonical) {
          cleaned$sequence <- canonical
          cleaned$valid <- TRUE
          metadata_recovered <- metadata_recovered + 1L
          corrected_ids <- base::c(corrected_ids, current_id)
        } else {
          if (!base::nzchar(cleaned$sequence)) {
            empty_dropped <- empty_dropped + 1L
          } else {
            invalid_dropped <- invalid_dropped + 1L
          }
          dropped_ids <- base::c(dropped_ids, current_id)
          next
        }
      } else if (has_canonical && !base::identical(cleaned$sequence, canonical)) {
        cleaned$sequence <- canonical
        metadata_recovered <- metadata_recovered + 1L
        corrected_ids <- base::c(corrected_ids, current_id)
      }

      duplicate_index <- base::match(current_id, record_ids)
      if (!base::is.na(duplicate_index)) {
        if (base::identical(records[[duplicate_index]]$sequence, cleaned$sequence)) {
          exact_duplicates <- exact_duplicates + 1L
          dropped_ids <- base::c(dropped_ids, current_id)
        } else {
          conflicting_duplicates <- conflicting_duplicates + 1L
          conflicting_ids <- base::c(conflicting_ids, current_id)
        }
        next
      }

      original_label <- base::trimws(base::sub("^[^[:space:]]+[[:space:]]*", "", header))
      label <- base::unname(metadata_names[current_id])
      if (!base::length(label) || base::is.na(label) || !base::nzchar(base::trimws(label))) {
        label <- original_label
      }
      if (base::is.na(label) || !base::nzchar(base::trimws(label))) label <- current_id
      label <- base::trimws(base::gsub("[[:space:]>]+", " ", base::as.character(label)[1]))
      records[[base::length(records) + 1L]] <- list(
        id = current_id,
        label = label,
        sequence = cleaned$sequence
      )
      record_ids <- base::c(record_ids, current_id)
    }
  }

  if (base::nrow(metadata_sequences)) {
    missing_metadata <- base::setdiff(metadata_sequences$pazy_id, record_ids)
    for (current_id in missing_metadata) {
      label <- metadata_sequences$pazy_name[base::match(current_id, metadata_sequences$pazy_id)]
      if (base::is.na(label) || !base::nzchar(base::trimws(label))) label <- current_id
      records[[base::length(records) + 1L]] <- list(
        id = current_id,
        label = base::trimws(base::gsub("[[:space:]>]+", " ", label)),
        sequence = base::unname(canonical_sequences[current_id])
      )
      record_ids <- base::c(record_ids, current_id)
      metadata_added <- metadata_added + 1L
      corrected_ids <- base::c(corrected_ids, current_id)
    }
  }

  fatal_detail <- character()
  if (!base::length(header_index) && !base::nrow(metadata_sequences)) {
    fatal_detail <- base::c(fatal_detail, "no FASTA headers")
  }
  if (base::length(orphan_lines)) {
    fatal_detail <- base::c(fatal_detail, "residues before the first FASTA header")
  }
  if (empty_headers > 0L) {
    fatal_detail <- base::c(fatal_detail, base::paste0(empty_headers, " empty FASTA header(s)"))
  }
  if (conflicting_duplicates > 0L) {
    fatal_detail <- base::c(
      fatal_detail,
      base::paste0("conflicting duplicate ID(s): ", base::paste(base::unique(conflicting_ids), collapse = ", "))
    )
  }
  if (!base::length(records)) {
    fatal_detail <- base::c(fatal_detail, "no valid PAZy protein records")
  }
  ok <- !base::length(fatal_detail)
  clean <- base::isTRUE(ok) &&
    whitespace_fixed == 0L && case_fixed == 0L && terminal_stop_fixed == 0L &&
    empty_dropped == 0L && invalid_dropped == 0L && metadata_recovered == 0L &&
    metadata_added == 0L && exact_duplicates == 0L && empty_headers == 0L &&
    !base::length(orphan_lines) && ids_normalized == 0L

  audit <- list(
    ok = ok,
    clean = clean,
    n_input_records = base::length(header_index),
    n_output_records = base::length(records),
    n_whitespace_fixed = whitespace_fixed,
    n_whitespace_removed = whitespace_removed,
    n_case_fixed = case_fixed,
    n_terminal_stop_fixed = terminal_stop_fixed,
    n_empty_dropped = empty_dropped,
    n_invalid_dropped = invalid_dropped,
    n_metadata_recovered = metadata_recovered,
    n_metadata_added = metadata_added,
    n_exact_duplicates_dropped = exact_duplicates,
    n_conflicting_duplicates = conflicting_duplicates,
    n_empty_headers = empty_headers,
    n_orphan_lines = base::length(orphan_lines),
    n_ids_normalized = ids_normalized,
    corrected_ids = base::unique(corrected_ids),
    dropped_ids = base::unique(dropped_ids),
    detail = if (ok) "" else base::paste(fatal_detail, collapse = "; ")
  )
  if (ok) audit$detail <- .dnmb_pazy_fasta_audit_detail(audit)

  if (base::isTRUE(ok) && !base::is.null(dest_fasta)) {
    dest_fasta <- base::path.expand(base::as.character(dest_fasta)[1])
    base::dir.create(base::dirname(dest_fasta), recursive = TRUE, showWarnings = FALSE)
    output <- base::unlist(base::lapply(records, function(record) {
      sequence_start <- base::seq.int(1L, base::nchar(record$sequence), by = 80L)
      sequence_lines <- base::substring(record$sequence, sequence_start, sequence_start + 79L)
      base::c(base::paste0(">", record$id, " ", record$label), sequence_lines)
    }), use.names = FALSE)
    base::writeLines(output, dest_fasta, useBytes = TRUE)
  }
  audit
}

.dnmb_pazy_rewrite_fasta_with_ids <- function(source_fasta, dest_fasta, metadata_flat) {
  audit <- .dnmb_pazy_sanitize_fasta(source_fasta, dest_fasta, metadata_flat = metadata_flat)
  if (!base::isTRUE(audit$ok)) {
    base::stop("Unsafe PAZy FASTA: ", audit$detail, call. = FALSE)
  }
  base::invisible(dest_fasta)
}

.dnmb_pazy_blast_db_ready <- function(db_prefix) {
  base::all(base::vapply(
    base::paste0(db_prefix, c(".phr", ".pin", ".psq")),
    .dnmb_nonempty_file,
    logical(1)
  ))
}

.dnmb_pazy_prepare_blast_db <- function(fasta_path, db_prefix, trace_log = NULL) {
  fasta_audit <- .dnmb_pazy_sanitize_fasta(fasta_path)
  if (!base::isTRUE(fasta_audit$ok) || !base::isTRUE(fasta_audit$clean)) {
    return(list(
      ok = FALSE,
      status = "failed",
      detail = base::paste0("Unsafe PAZy FASTA: ", fasta_audit$detail %||% fasta_path)
    ))
  }
  output_files <- base::paste0(db_prefix, c(".phr", ".pin", ".psq"))
  if (.dnmb_pazy_blast_db_ready(db_prefix)) {
    return(list(ok = TRUE, status = "cached", detail = db_prefix))
  }
  stage_dir <- tempfile("dnmb-pazy-makeblastdb-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged_faa <- base::file.path(stage_dir, "pazy_reference.faa")
  staged_db <- base::file.path(stage_dir, "pazy_reference")
  copied_fasta <- base::file.copy(fasta_path, staged_faa, overwrite = TRUE)
  if (!base::isTRUE(copied_fasta) || !.dnmb_nonempty_file(staged_faa)) {
    return(list(ok = FALSE, status = "failed", detail = "Could not stage the PAZy FASTA for makeblastdb."))
  }
  args <- c("-in", staged_faa, "-dbtype", "prot", "-out", staged_db)
  if (!base::is.null(trace_log) && base::nzchar(trace_log)) {
    .dnmb_pazy_trace(trace_log, base::sprintf("[%s] makeblastdb %s", base::Sys.time(), .dnmb_format_command("makeblastdb", args)))
  }
  run <- dnmb_run_external(
    "makeblastdb",
    args = args,
    required = FALSE,
    timeout = .dnmb_pazy_makedb_timeout()
  )
  staged_files <- base::paste0(staged_db, c(".phr", ".pin", ".psq"))
  ok <- base::isTRUE(run$ok) && base::all(base::vapply(staged_files, .dnmb_nonempty_file, logical(1)))
  if (ok) {
    # Copy all BLAST DB files including v5 LMDB files (.pdb, .pot, .ptf, .pto)
    all_staged <- base::list.files(stage_dir, pattern = base::paste0("^", base::basename(staged_db), "\\."),
                                   full.names = TRUE)
    dest_dir <- base::dirname(db_prefix)
    copied <- base::vapply(all_staged, function(sf) {
      dest_name <- base::file.path(dest_dir, sub(base::basename(staged_db), base::basename(db_prefix), base::basename(sf), fixed = TRUE))
      base::isTRUE(base::file.copy(sf, dest_name, overwrite = TRUE))
    }, logical(1))
    ok <- base::length(copied) > 0L && base::all(copied) &&
      .dnmb_pazy_blast_db_ready(db_prefix)
  }
  list(
    ok = ok,
    status = if (ok) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    detail = if (ok) db_prefix else (run$error %||% db_prefix),
    command = run
  )
}

dnmb_pazy_install_module <- function(version = .dnmb_pazy_default_version(),
                                     cache_root = NULL,
                                     install = TRUE,
                                     metadata_url = .dnmb_pazy_default_metadata_url(),
                                     fasta_url = .dnmb_pazy_default_fasta_url(),
                                     asset_urls = NULL,
                                     force = FALSE) {
  module <- .dnmb_pazy_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_pazy_asset_layout(module_dir)
  trace_log <- base::file.path(module_dir, "pazy_install_trace.log")
  asset_urls <- .dnmb_pazy_normalize_asset_urls(asset_urls)
  metadata_source <- asset_urls$metadata_tsv %||% asset_urls$metadata_json %||% metadata_url
  fasta_source <- asset_urls$fasta_faa %||% asset_urls$reference_fasta %||% fasta_url
  status <- .dnmb_pazy_empty_status()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  remote_state <- .dnmb_pazy_remote_asset_state(
    metadata_url = metadata_source,
    fasta_url = fasta_source,
    enabled = base::isTRUE(install) &&
      !base::file.exists(metadata_source) && !base::dir.exists(metadata_source) &&
      !base::file.exists(fasta_source) && !base::dir.exists(fasta_source)
  )
  remote_refresh_needed <- .dnmb_pazy_remote_update_needed(manifest, remote_state)

  local_reference_available <- base::all(base::vapply(
    c(layout$reference_fasta, layout$metadata_tsv),
    .dnmb_nonempty_file,
    logical(1)
  ))
  cache_assets_ready <- local_reference_available &&
    .dnmb_nonempty_file(layout$metadata_json) &&
    .dnmb_pazy_blast_db_ready(layout$blast_db_prefix)
  cache_audit <- if (.dnmb_nonempty_file(layout$reference_fasta)) {
    .dnmb_pazy_sanitize_fasta(layout$reference_fasta)
  } else {
    list(ok = FALSE, clean = FALSE, detail = "PAZy reference FASTA is missing or empty.")
  }
  contract_current <- !base::is.null(manifest) &&
    base::identical(
      base::suppressWarnings(base::as.integer(manifest$fasta_contract_version)[1]),
      .dnmb_pazy_fasta_contract_version()
    )
  current_reference_md5 <- if (.dnmb_nonempty_file(layout$reference_fasta)) {
    base::unname(tools::md5sum(layout$reference_fasta)[[1]])
  } else {
    NA_character_
  }
  manifest_reference_md5 <- if (!base::is.null(manifest$reference_md5)) {
    base::as.character(manifest$reference_md5)[1]
  } else {
    NA_character_
  }
  reference_signature_current <- !base::is.na(current_reference_md5) &&
    !base::is.na(manifest_reference_md5) &&
    base::identical(current_reference_md5, manifest_reference_md5)
  cache_safe <- cache_assets_ready && !base::is.null(manifest) &&
    base::isTRUE(manifest$install_ok) && contract_current &&
    base::isTRUE(cache_audit$ok) && base::isTRUE(cache_audit$clean) &&
    reference_signature_current

  if (cache_safe && !base::isTRUE(force) && !base::isTRUE(remote_refresh_needed)) {
    return(list(
      ok = TRUE,
      status = .dnmb_pazy_status_row("pazy_install", "cached", module_dir),
      files = c(
        list(reference_fasta = layout$reference_fasta, metadata_tsv = layout$metadata_tsv),
        stats::setNames(as.list(layout$blast_db_files), base::basename(layout$blast_db_files))
      ),
      manifest = manifest
    ))
  }

  explicit_assets <- base::length(asset_urls) > 0L
  local_cache_repair <- local_reference_available && !base::isTRUE(force) &&
    !base::isTRUE(remote_refresh_needed) && !base::isTRUE(explicit_assets) &&
    !base::isTRUE(cache_safe)
  manifest_metadata_source <- metadata_source
  manifest_fasta_source <- fasta_source
  if (local_cache_repair) {
    # Migrate old installations without network access. The original generation
    # remains untouched unless sanitation and makeblastdb both pass.
    metadata_source <- layout$metadata_tsv
    fasta_source <- layout$reference_fasta
    manifest_metadata_source <- manifest$metadata_url %||% manifest_metadata_source
    manifest_fasta_source <- manifest$fasta_url %||% manifest_fasta_source
    status <- dplyr::bind_rows(
      status,
      .dnmb_pazy_status_row("pazy_cache_audit", "repair", cache_audit$detail %||% module_dir)
    )
  }

  if (!base::isTRUE(install) && !base::all(base::file.exists(c(layout$reference_fasta, layout$metadata_tsv)))) {
    return(list(
      ok = FALSE,
      status = .dnmb_pazy_status_row("pazy_install", "missing", "PAZy reference is missing and module_install is FALSE."),
      files = list(),
      manifest = NULL
    ))
  }

  if (base::interactive() && base::isTRUE(install) && !base::all(base::file.exists(c(layout$reference_fasta, layout$metadata_tsv)))) {
    answer <- utils::askYesNo("PAZy reference data are missing. Download the current PAZy protein set from the official API?")
    if (!base::isTRUE(answer)) {
      return(list(
        ok = FALSE,
        status = .dnmb_pazy_status_row("pazy_install", "declined", "User declined PAZy download."),
        files = list(),
        manifest = NULL
      ))
    }
  }

  stage_dir <- base::tempfile(".dnmb-pazy-generation-", tmpdir = module_dir)
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  base::on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  stage_layout <- .dnmb_pazy_asset_layout(stage_dir)

  metadata_flat <- if (base::dir.exists(metadata_source)) {
    utils::read.delim(base::file.path(metadata_source, "pazy_metadata.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  } else if (base::file.exists(metadata_source) && grepl("\\.tsv$", metadata_source, ignore.case = TRUE)) {
    utils::read.delim(metadata_source, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  } else if (base::file.exists(metadata_source) && grepl("\\.json$", metadata_source, ignore.case = TRUE)) {
    jsonlite::fromJSON(metadata_source, simplifyVector = TRUE)
  } else {
    .dnmb_pazy_flatten_metadata(.dnmb_pazy_fetch_metadata_pages(metadata_source))
  }
  utils::write.table(metadata_flat, file = stage_layout$metadata_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  jsonlite::write_json(metadata_flat, path = stage_layout$metadata_json, auto_unbox = TRUE, pretty = TRUE)
  status <- dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_metadata", "ok", layout$metadata_tsv))

  raw_fasta <- base::tempfile("dnmb-pazy-raw-", tmpdir = stage_dir, fileext = ".faa")
  copied_fasta <- TRUE
  if (base::dir.exists(fasta_source)) {
    copied_fasta <- base::file.copy(base::file.path(fasta_source, "pazy_reference.faa"), raw_fasta, overwrite = TRUE)
  } else if (base::file.exists(fasta_source)) {
    copied_fasta <- base::file.copy(fasta_source, raw_fasta, overwrite = TRUE)
  } else {
    download <- .dnmb_download_asset(
      fasta_source,
      raw_fasta,
      insecure = FALSE,
      timeout = .dnmb_pazy_network_timeout()
    )
    copied_fasta <- base::isTRUE(download$ok)
    if (!base::isTRUE(copied_fasta) || !.dnmb_nonempty_file(raw_fasta)) {
      return(list(
        ok = FALSE,
        status = dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_fasta", "failed", download$error %||% fasta_source)),
        files = list(),
        manifest = NULL
      ))
    }
  }
  if (!base::isTRUE(copied_fasta) || !.dnmb_nonempty_file(raw_fasta)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_fasta", "failed", fasta_source)),
      files = list(),
      manifest = manifest
    ))
  }
  fasta_audit <- .dnmb_pazy_sanitize_fasta(
    raw_fasta,
    stage_layout$reference_fasta,
    metadata_flat = metadata_flat
  )
  if (!base::isTRUE(fasta_audit$ok) || !.dnmb_nonempty_file(stage_layout$reference_fasta)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(
        status,
        .dnmb_pazy_status_row(
          "pazy_fasta",
          "failed",
          base::paste0("Unsafe PAZy FASTA: ", fasta_audit$detail %||% fasta_source)
        )
      ),
      files = list(),
      manifest = manifest
    ))
  }
  status <- dplyr::bind_rows(
    status,
    .dnmb_pazy_status_row(
      "pazy_fasta",
      if (local_cache_repair || !base::isTRUE(fasta_audit$clean)) "repaired" else "ok",
      .dnmb_pazy_fasta_audit_detail(fasta_audit)
    )
  )

  blast_prepare <- .dnmb_pazy_prepare_blast_db(
    stage_layout$reference_fasta,
    stage_layout$blast_db_prefix,
    trace_log = trace_log
  )
  status <- dplyr::bind_rows(
    status,
    .dnmb_pazy_status_row(
      "pazy_prepare",
      blast_prepare$status,
      blast_prepare$detail %||% layout$blast_db_prefix
    )
  )
  if (!base::isTRUE(blast_prepare$ok)) {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  staged_db_files <- .dnmb_pazy_find_blast_db_files(stage_layout$blast_db_prefix)
  staged_db_files <- staged_db_files[base::file.exists(staged_db_files)]
  staged_assets <- base::c(
    stage_layout$reference_fasta,
    stage_layout$metadata_json,
    stage_layout$metadata_tsv,
    staged_db_files
  )
  destination_assets <- base::file.path(module_dir, base::basename(staged_assets))
  old_db_files <- base::list.files(
    module_dir,
    pattern = "^pazy_reference[.](phr|pin|psq|pdb|pot|ptf|pto)$",
    full.names = TRUE
  )
  commit <- .dnmb_transactional_replace(
    staged_paths = staged_assets,
    destination_paths = destination_assets,
    retire_paths = base::setdiff(old_db_files, destination_assets)
  )
  if (!base::isTRUE(commit$ok)) {
    status <- dplyr::bind_rows(
      status,
      .dnmb_pazy_status_row("pazy_commit", "failed", commit$detail)
    )
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }
  layout <- .dnmb_pazy_asset_layout(module_dir)

  manifest <- list(
    install_ok = TRUE,
    module = module,
    version = version,
    module_dir = module_dir,
    metadata_url = manifest_metadata_source,
    fasta_url = manifest_fasta_source,
    reference_fasta = layout$reference_fasta,
    metadata_tsv = layout$metadata_tsv,
    blast_db_prefix = layout$blast_db_prefix,
    metadata_count = base::nrow(metadata_flat),
    reference_count = fasta_audit$n_output_records,
    reference_md5 = base::unname(tools::md5sum(layout$reference_fasta)[[1]]),
    fasta_contract_version = .dnmb_pazy_fasta_contract_version(),
    fasta_audit = .dnmb_pazy_fasta_audit_summary(fasta_audit),
    remote_asset_state = remote_state %||% manifest$remote_asset_state
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_pazy_default_version(),
    cache_root = cache_root
  )

  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_install", "ok", module_dir)),
    files = c(
      list(reference_fasta = layout$reference_fasta, metadata_tsv = layout$metadata_tsv, trace_log = trace_log),
      stats::setNames(as.list(layout$blast_db_files), base::basename(layout$blast_db_files))
    ),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_pazy_get_module <- function(version = .dnmb_pazy_default_version(),
                                 cache_root = NULL,
                                 required = FALSE) {
  manifest <- dnmb_db_read_manifest(.dnmb_pazy_module_name(), version, cache_root = cache_root, required = FALSE)
  layout <- .dnmb_pazy_asset_layout(.dnmb_db_module_dir(.dnmb_pazy_module_name(), version, cache_root = cache_root, create = FALSE))
  fasta_audit <- if (.dnmb_nonempty_file(layout$reference_fasta)) {
    .dnmb_pazy_sanitize_fasta(layout$reference_fasta)
  } else {
    list(ok = FALSE, clean = FALSE)
  }
  reference_md5 <- if (.dnmb_nonempty_file(layout$reference_fasta)) {
    base::unname(tools::md5sum(layout$reference_fasta)[[1]])
  } else {
    NA_character_
  }
  manifest_md5 <- if (!base::is.null(manifest$reference_md5)) {
    base::as.character(manifest$reference_md5)[1]
  } else {
    NA_character_
  }
  ok <- !base::is.null(manifest) && base::isTRUE(manifest$install_ok) &&
    base::identical(
      base::suppressWarnings(base::as.integer(manifest$fasta_contract_version)[1]),
      .dnmb_pazy_fasta_contract_version()
    ) &&
    base::isTRUE(fasta_audit$ok) && base::isTRUE(fasta_audit$clean) &&
    !base::is.na(reference_md5) && base::identical(reference_md5, manifest_md5) &&
    .dnmb_nonempty_file(layout$metadata_tsv) &&
    .dnmb_pazy_blast_db_ready(layout$blast_db_prefix)
  if (required && !ok) {
    base::stop("PAZy module is not installed for version `", version, "`.", call. = FALSE)
  }
  list(
    ok = ok,
    manifest = manifest,
    files = list(
      reference_fasta = layout$reference_fasta,
      metadata_tsv = layout$metadata_tsv,
      blast_db_prefix = layout$blast_db_prefix,
      blast_db_files = layout$blast_db_files
    )
  )
}

dnmb_pazy_parse_metadata <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

dnmb_pazy_parse_blast_tabular <- function(path, metadata) {
  if (!base::file.exists(path) || !base::isTRUE(base::file.info(path)$size > 0)) {
    return(data.frame())
  }
  tbl <- utils::read.delim(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (base::ncol(tbl) < 12L) {
    return(data.frame())
  }
  tbl <- tbl[, base::seq_len(base::min(13L, base::ncol(tbl))), drop = FALSE]
  base::colnames(tbl) <- c(
    "query", "subject_id", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    if (base::ncol(tbl) >= 13L) "subject_title" else NULL
  )
  tbl$pident <- base::as.numeric(tbl$pident)
  tbl$length <- base::as.integer(tbl$length)
  tbl$qlen <- base::as.integer(tbl$qlen)
  tbl$slen <- base::as.integer(tbl$slen)
  tbl$evalue <- base::as.numeric(tbl$evalue)
  tbl$bitscore <- base::as.numeric(tbl$bitscore)
  tbl$qcov <- tbl$length / tbl$qlen
  tbl$scov <- tbl$length / tbl$slen
  tbl <- dplyr::left_join(tbl, metadata, by = c("subject_id" = "pazy_id"))
  tbl
}

dnmb_pazy_normalize_hits <- function(hits) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  support <- base::paste0(
    "pazy_id=", hits$subject_id,
    "; evalue=", base::signif(hits$evalue, 4),
    "; bitscore=", base::signif(hits$bitscore, 4),
    "; pident=", base::sprintf("%.1f", hits$pident),
    "; qcov=", base::sprintf("%.3f", hits$qcov),
    "; scov=", base::sprintf("%.3f", hits$scov),
    "; verified=", ifelse(hits$verified %in% TRUE, "TRUE", "FALSE"),
    "; substrate=", hits$substrate_abbreviations
  )

  out <- base::data.frame(
    query = .dnmb_module_clean_annotation_key(hits$query),
    source = "pazy",
    family_system = "PAZy",
    family_id = base::as.character(hits$pazy_name),
    hit_label = base::as.character(hits$pazy_name),
    enzyme_role = "plastic-active enzyme",
    evidence_mode = "direct",
    substrate_label = base::as.character(hits$substrate_abbreviations),
    support = support,
    pazy_id = base::as.character(hits$subject_id),
    verified = base::as.logical(hits$verified),
    accession_list = base::as.character(hits$accession_list),
    organism_name = base::as.character(hits$organism_name),
    pident = base::as.numeric(hits$pident),
    alignment_length = base::as.integer(hits$length),
    query_length = base::as.integer(hits$qlen),
    subject_length = base::as.integer(hits$slen),
    evalue = base::as.numeric(hits$evalue),
    bitscore = base::as.numeric(hits$bitscore),
    qcov = base::as.numeric(hits$qcov),
    scov = base::as.numeric(hits$scov),
    typing_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  out <- out[, c(.dnmb_module_optional_long_columns(), base::setdiff(base::names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_pazy_append_evalue_threshold <- function() {
  1e-10
}

.dnmb_pazy_append_qcov_threshold <- function() {
  0.5
}

.dnmb_pazy_append_scov_threshold <- function() {
  0.5
}

.dnmb_pazy_hits_for_output <- function(hits,
                                       evalue_threshold = .dnmb_pazy_append_evalue_threshold(),
                                       qcov_threshold = .dnmb_pazy_append_qcov_threshold(),
                                       scov_threshold = .dnmb_pazy_append_scov_threshold()) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  out <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::all(c("evalue", "qcov", "scov") %in% base::names(out))) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out$evalue <- suppressWarnings(base::as.numeric(out$evalue))
  out$qcov <- suppressWarnings(base::as.numeric(out$qcov))
  out$scov <- suppressWarnings(base::as.numeric(out$scov))
  out <- out[
    !base::is.na(out$evalue) &
      !base::is.na(out$qcov) &
      !base::is.na(out$scov) &
      out$evalue <= base::as.numeric(evalue_threshold)[1] &
      out$qcov >= base::as.numeric(qcov_threshold)[1] &
      out$scov >= base::as.numeric(scov_threshold)[1],
    ,
    drop = FALSE
  ]
  base::rownames(out) <- NULL
  if (!base::nrow(out)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out
}

.dnmb_pazy_output_table <- function(genes,
                                    hits,
                                    evalue_threshold = .dnmb_pazy_append_evalue_threshold(),
                                    qcov_threshold = .dnmb_pazy_append_qcov_threshold(),
                                    scov_threshold = .dnmb_pazy_append_scov_threshold()) {
  filtered_hits <- .dnmb_pazy_hits_for_output(
    hits = hits,
    evalue_threshold = evalue_threshold,
    qcov_threshold = qcov_threshold,
    scov_threshold = scov_threshold
  )
  out <- .dnmb_module_output_table(genes = genes, hits = filtered_hits)
  drop_cols <- base::intersect(
    c("enzyme_role", "evidence_mode", "support", "typing_eligible", "verified"),
    base::names(out)
  )
  if (base::length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  base_cols <- base::intersect(dnmb_backbone_columns(), base::names(out))
  pazy_cols <- c("family_id", "hit_label", "substrate_label", base::setdiff(base::names(out), c(base_cols, "family_id", "hit_label", "substrate_label")))
  out[, c(base_cols, pazy_cols), drop = FALSE]
}

.dnmb_pazy_reference_state <- function(reference_fasta) {
  if (!.dnmb_nonempty_file(reference_fasta)) return(NULL)
  audit <- .dnmb_pazy_sanitize_fasta(reference_fasta)
  if (!base::isTRUE(audit$ok) || !base::isTRUE(audit$clean)) return(NULL)
  list(
    fasta_contract_version = .dnmb_pazy_fasta_contract_version(),
    reference_md5 = base::unname(tools::md5sum(reference_fasta)[[1]]),
    reference_size = base::unname(base::as.numeric(base::file.info(reference_fasta)$size))
  )
}

.dnmb_pazy_dmnd_state_path <- function(dmnd_db) {
  base::paste0(dmnd_db, ".state.rds")
}

.dnmb_pazy_dmnd_is_current <- function(dmnd_db, reference_fasta) {
  state_path <- .dnmb_pazy_dmnd_state_path(dmnd_db)
  expected <- .dnmb_pazy_reference_state(reference_fasta)
  if (!.dnmb_nonempty_file(dmnd_db) || base::is.null(expected) ||
      !.dnmb_nonempty_file(state_path)) {
    return(FALSE)
  }
  stored <- base::tryCatch(base::readRDS(state_path), error = function(e) NULL)
  base::is.list(stored) &&
    base::identical(stored$fasta_contract_version, expected$fasta_contract_version) &&
    base::identical(stored$reference_md5, expected$reference_md5) &&
    base::identical(stored$reference_size, expected$reference_size)
}

.dnmb_pazy_prepare_diamond_db <- function(reference_fasta, dmnd_db, trace_log = NULL) {
  expected_state <- .dnmb_pazy_reference_state(reference_fasta)
  if (base::is.null(expected_state)) {
    return(list(ok = FALSE, status = "failed", detail = "PAZy reference FASTA is missing or empty."))
  }
  if (.dnmb_pazy_dmnd_is_current(dmnd_db, reference_fasta)) {
    return(list(ok = TRUE, status = "cached", detail = dmnd_db, command = NULL))
  }

  stage_prefix <- base::tempfile(".pazy-reference-", tmpdir = base::dirname(dmnd_db))
  stage_candidates <- c(stage_prefix, base::paste0(stage_prefix, ".dmnd"))
  stage_state <- base::paste0(stage_prefix, ".state.rds")
  base::on.exit(base::unlink(c(stage_candidates, stage_state), force = TRUE), add = TRUE)
  args <- c("makedb", "--in", reference_fasta, "-d", stage_prefix)
  if (!base::is.null(trace_log) && base::nzchar(trace_log)) {
    .dnmb_pazy_trace(
      trace_log,
      base::sprintf("[%s] diamond %s", base::Sys.time(), .dnmb_format_command("diamond", args))
    )
  }
  run <- dnmb_run_external(
    "diamond",
    args = args,
    required = FALSE,
    timeout = .dnmb_pazy_makedb_timeout()
  )
  generated <- stage_candidates[base::vapply(stage_candidates, .dnmb_nonempty_file, logical(1))]
  if (!base::isTRUE(run$ok) || base::length(generated) != 1L) {
    return(list(
      ok = FALSE,
      status = if (!base::nzchar(run$resolved_command)) "missing" else "failed",
      detail = run$error %||% "DIAMOND did not create a complete PAZy database.",
      command = run
    ))
  }

  base::saveRDS(expected_state, stage_state)
  commit <- .dnmb_transactional_replace(
    staged_paths = c(generated[[1]], stage_state),
    destination_paths = c(dmnd_db, .dnmb_pazy_dmnd_state_path(dmnd_db))
  )
  ok <- base::isTRUE(commit$ok) && .dnmb_pazy_dmnd_is_current(dmnd_db, reference_fasta)
  list(
    ok = ok,
    status = if (ok) "ok" else "failed",
    detail = if (ok) dmnd_db else commit$detail,
    command = run
  )
}

dnmb_run_pazy_module <- function(genes,
                                 output_dir,
                                 version = .dnmb_pazy_default_version(),
                                 cache_root = NULL,
                                 install = TRUE,
                                 metadata_url = .dnmb_pazy_default_metadata_url(),
                                 fasta_url = .dnmb_pazy_default_fasta_url(),
                                 asset_urls = NULL,
                                 cpu = 1L,
                                 genbank = NULL,
                                 search_backend = c("diamond", "blastp")) {
  search_backend <- base::match.arg(search_backend)
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- base::file.path(output_dir, "pazy_module_trace.log")
  fasta_path <- base::file.path(output_dir, "pazy_query_proteins.faa")
  existing_faa <- dnmb_resolve_query_faa(genbank = genbank, output_dir = output_dir, fallback_filename = base::basename(fasta_path))
  if (!base::is.null(existing_faa) && .dnmb_can_reuse_query_fasta(existing_faa, genes)) {
    proteins <- .dnmb_prepare_query_proteins(genes)
    fasta <- list(path = existing_faa, n = base::nrow(proteins), proteins = proteins)
    fasta_path <- existing_faa
  } else {
    fasta <- .dnmb_write_query_fasta(genes, fasta_path)
  }

  status <- .dnmb_pazy_status_row("pazy_query_fasta", if (fasta$n) "ok" else "empty", base::paste0("proteins=", fasta$n))
  if (!fasta$n) {
    return(list(ok = TRUE, status = status, files = list(query_fasta = fasta_path), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  install_result <- dnmb_pazy_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    metadata_url = metadata_url,
    fasta_url = fasta_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = install_result$files, hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  module <- dnmb_pazy_get_module(version = version, cache_root = cache_root, required = TRUE)
  blast_out <- base::file.path(output_dir, "pazy_blastp.tsv")
  command <- NULL

  # --- diamond backend ---
  if (base::identical(search_backend, "diamond")) {
    diamond_check <- dnmb_run_external("diamond", args = "version", required = FALSE)
    if (base::nzchar(diamond_check$resolved_command)) {
      dmnd_db <- base::file.path(output_dir, "pazy_reference.dmnd")
      diamond_db <- .dnmb_pazy_prepare_diamond_db(
        reference_fasta = module$files$reference_fasta,
        dmnd_db = dmnd_db,
        trace_log = trace_log
      )
      if (!base::isTRUE(diamond_db$ok)) {
        .dnmb_pazy_trace(
          trace_log,
          base::sprintf(
            "[%s] diamond makedb %s, falling back to blastp",
            base::Sys.time(),
            if (base::isTRUE(diamond_db$command$timed_out)) "timed out" else "failed"
          )
        )
        search_backend <- "blastp"
      }
      if (base::identical(search_backend, "diamond")) {
        stage_out <- base::file.path(output_dir, "pazy_blastp.tsv")
        d_args <- c(
          "blastp", "-q", fasta_path, "-d", dmnd_db,
          "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
          "-k", "5", "--threads", base::as.character(base::as.integer(cpu)[1]),
          "--sensitive", "-o", stage_out
        )
        .dnmb_pazy_trace(trace_log, base::sprintf("[%s] diamond blastp --sensitive threads=%s", base::Sys.time(), cpu))
        command <- dnmb_run_external(
          "diamond",
          args = d_args,
          required = FALSE,
          timeout = .dnmb_pazy_search_timeout()
        )
      }
    } else {
      .dnmb_pazy_trace(trace_log, base::sprintf("[%s] diamond not found, falling back to blastp", base::Sys.time()))
      search_backend <- "blastp"
    }
  }

  # --- blastp backend (fallback) ---
  if (base::identical(search_backend, "blastp") || base::is.null(command)) {
    stage_dir <- tempfile("dnmb-pazy-blastp-")
    base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
    stage_query <- base::file.path(stage_dir, "query.faa")
    stage_db <- base::file.path(stage_dir, "pazy_reference")
    stage_out <- base::file.path(stage_dir, "pazy_blastp.tsv")
    base::file.copy(fasta_path, stage_query, overwrite = TRUE)
    base::file.copy(module$files$reference_fasta, base::file.path(stage_dir, "pazy_reference.faa"), overwrite = TRUE)
    db_src <- module$files$blast_db_files
    if (base::is.character(db_src)) {
      db_src <- db_src[base::file.exists(db_src)]
      db_dest <- base::paste0(stage_db, ".", tools::file_ext(db_src))
      base::file.copy(db_src, db_dest, overwrite = TRUE)
    }
    args <- c(
      "-query", stage_query, "-db", stage_db,
      "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle",
      "-max_target_seqs", "5",
      "-num_threads", base::as.character(base::as.integer(cpu)[1]),
      "-out", stage_out
    )
    .dnmb_pazy_trace(trace_log, base::sprintf("[%s] blastp %s", base::Sys.time(), .dnmb_format_command("blastp", args)))
    command <- dnmb_run_external(
      "blastp",
      args = args,
      required = FALSE,
      timeout = .dnmb_pazy_search_timeout()
    )
    if (base::isTRUE(command$ok) && base::file.exists(stage_out)) {
      base::file.copy(stage_out, blast_out, overwrite = TRUE)
    }
  }

  if (base::isTRUE(command$ok) && !base::file.exists(blast_out)) {
    base::writeLines(character(0), blast_out)
  }
  metadata <- dnmb_pazy_parse_metadata(module$files$metadata_tsv)
  base::names(metadata)[base::match("pazy_id", base::names(metadata))] <- "pazy_id"
  parsed <- if (base::isTRUE(command$ok) && base::file.exists(blast_out)) dnmb_pazy_parse_blast_tabular(blast_out, metadata = metadata) else data.frame()
  hits <- dnmb_pazy_normalize_hits(parsed)
  status <- dplyr::bind_rows(
    status,
    .dnmb_pazy_status_row(
      "pazy_blastp",
      if (base::isTRUE(command$ok)) if (base::nrow(hits)) "ok" else "empty" else if (!base::nzchar(command$resolved_command)) "missing" else "failed",
      if (base::isTRUE(command$ok)) blast_out else (command$error %||% "PAZy blastp failed")
    )
  )

  out_table <- .dnmb_pazy_output_table(genes = genes, hits = hits)

  # Write pazy_merged.tsv for downstream visualization (Panel A/B/C pub figure)
  merged_path <- base::file.path(output_dir, "pazy_merged.tsv")
  base::tryCatch({
    if (base::is.data.frame(parsed) && base::nrow(parsed) > 0L) {
      backbone <- c("locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction")
      base_genes <- genes[, base::intersect(backbone, base::names(genes)), drop = FALSE]
      # Rename parsed columns to match _pub expectations
      p <- parsed
      p$locus_tag <- .dnmb_module_clean_annotation_key(p$query)
      col_map <- c(pazy_name = "family_id", subject_id = "pazy_id",
                   substrate_abbreviations = "substrate_label",
                   length = "alignment_length", qlen = "query_length", slen = "subject_length")
      for (old_nm in base::names(col_map)) {
        if (old_nm %in% base::names(p)) base::names(p)[base::names(p) == old_nm] <- col_map[[old_nm]]
      }
      # Best hit per locus_tag
      p <- p[base::order(p$evalue, -p$bitscore), , drop = FALSE]
      p <- p[!base::duplicated(p$locus_tag), , drop = FALSE]
      merge_cols <- base::setdiff(base::names(p), c("query", backbone))
      merged <- base::merge(base_genes, p[, c("locus_tag", merge_cols), drop = FALSE],
                            by = "locus_tag", all.x = TRUE)
      rename_idx <- !base::names(merged) %in% backbone
      base::names(merged)[rename_idx] <- base::paste0("PAZy_", base::names(merged)[rename_idx])
      utils::write.table(merged, merged_path, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }, error = function(e) NULL)

  list(
    ok = base::isTRUE(command$ok),
    status = status,
    files = c(install_result$files, list(query_fasta = fasta_path, blast_tsv = blast_out, trace_log = trace_log)),
    hits = hits,
    output_table = out_table,
    module_result = parsed
  )
}
