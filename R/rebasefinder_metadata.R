.dnmb_rebasefinder_header_field <- function(header, field) {
  match <- base::regmatches(
    header,
    base::regexec(base::paste0("(?:^|\\t)", field, ":([^\\t]*)"), header, perl = TRUE)
  )[[1]]
  if (base::length(match) < 2L) return(NA_character_)
  value <- base::trimws(match[[2]])
  if (!base::nzchar(value)) NA_character_ else value
}

.dnmb_rebasefinder_reference_header_path <- function(cache_root = NULL) {
  path <- base::file.path(
    .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = FALSE),
    "REBASE_protein_seqs.txt"
  )
  if (base::file.exists(path)) path else NA_character_
}

.dnmb_rebasefinder_reference_header_metadata <- function(subject_ids,
                                                         cache_root = NULL,
                                                         raw_path = NULL,
                                                         rebase_data = NULL) {
  subject_ids <- base::unique(base::as.character(subject_ids))
  subject_ids <- subject_ids[!base::is.na(subject_ids) & base::nzchar(subject_ids)]
  n_subjects <- base::length(subject_ids)
  empty <- base::data.frame(
    subject_id = subject_ids,
    enzyme_name = .dnmb_rebasefinder_clean_rebase_subject(subject_ids),
    enz_type = base::rep(NA_character_, n_subjects),
    rec_seq = base::rep(NA_character_, n_subjects),
    genbank_accession = base::rep(NA_character_, n_subjects),
    rebase_locus = base::rep(NA_character_, n_subjects),
    protein_id = base::rep(NA_character_, n_subjects),
    uniprot_id = base::rep(NA_character_, n_subjects),
    stringsAsFactors = FALSE
  )
  if (!base::length(subject_ids)) return(empty)

  if (base::is.null(raw_path)) raw_path <- .dnmb_rebasefinder_reference_header_path(cache_root)
  if (base::is.na(raw_path) || !base::file.exists(raw_path)) return(empty)
  info <- base::file.info(raw_path)
  signature <- base::paste(info$size, base::as.numeric(info$mtime), sep = ":")
  cache_path <- base::file.path(
    .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE),
    "rebase_reference_header_metadata.rds"
  )
  cached <- if (base::file.exists(cache_path)) {
    tryCatch(base::readRDS(cache_path), error = function(e) NULL)
  } else {
    NULL
  }
  cache_data <- base::data.frame()
  if (base::is.list(cached) && identical(cached$signature, signature) &&
      base::is.data.frame(cached$data) &&
      base::all(c(".raw_sequence", ".raw_occurrence", ".raw_org_name") %in%
                  base::names(cached$data))) {
    cache_data <- cached$data
  }
  base_names <- .dnmb_rebasefinder_clean_rebase_subject(subject_ids)
  known_names <- if (base::nrow(cache_data)) {
    base::unique(base::as.character(cache_data$enzyme_name))
  } else {
    character()
  }
  target_names <- base::setdiff(base::unique(base_names), known_names)

  if (base::length(target_names)) {
    found <- list()
    pending <- NULL
    append_pending <- function() {
      if (base::is.null(pending)) return(invisible(NULL))
      found[[base::length(found) + 1L]] <<- base::data.frame(
        enzyme_name = pending$enzyme_name,
        enz_type = .dnmb_rebasefinder_header_field(pending$header, "EnzType"),
        rec_seq = .dnmb_rebasefinder_header_field(pending$header, "RecSeq"),
        genbank_accession = .dnmb_rebasefinder_header_field(pending$header, "GenBank"),
        rebase_locus = .dnmb_rebasefinder_header_field(pending$header, "Locus"),
        protein_id = .dnmb_rebasefinder_header_field(pending$header, "ProteinId"),
        uniprot_id = .dnmb_rebasefinder_header_field(pending$header, "UniProt"),
        .raw_sequence = base::toupper(base::gsub("\\s+", "", pending$sequence, perl = TRUE)),
        .raw_occurrence = NA_integer_,
        .raw_org_name = .dnmb_rebasefinder_header_field(pending$header, "OrgName"),
        stringsAsFactors = FALSE
      )
      pending <<- NULL
      invisible(NULL)
    }
    con <- base::file(raw_path, open = "r")
    on.exit(base::close(con), add = TRUE)
    repeat {
      lines <- base::readLines(con, n = 100000L, warn = FALSE)
      if (!base::length(lines)) break
      header_idx <- base::which(base::startsWith(lines, ">"))
      first_header <- if (base::length(header_idx)) header_idx[[1]] else base::length(lines) + 1L
      if (!base::is.null(pending)) {
        before_header <- if (first_header > 1L) lines[base::seq_len(first_header - 1L)] else character()
        pending$sequence <- base::paste0(pending$sequence, base::paste(before_header, collapse = ""))
        if (base::length(header_idx)) append_pending()
      }
      if (!base::length(header_idx)) next
      for (j in base::seq_along(header_idx)) {
        idx <- header_idx[[j]]
        header <- lines[[idx]]
        enzyme <- base::sub("^>REBASE:([^\\t ]+).*$", "\\1", header, perl = TRUE)
        if (!enzyme %in% target_names) next
        next_idx <- if (j < base::length(header_idx)) header_idx[[j + 1L]] else base::length(lines) + 1L
        seq_lines <- if (idx + 1L < next_idx) lines[base::seq.int(idx + 1L, next_idx - 1L)] else character()
        pending <- list(
          enzyme_name = enzyme,
          header = header,
          sequence = base::paste(seq_lines, collapse = "")
        )
        if (j < base::length(header_idx)) append_pending()
      }
    }
    append_pending()
    found_data <- if (base::length(found)) {
      base::do.call(base::rbind, found)
    } else {
      base::data.frame()
    }
    found_names <- if (base::nrow(found_data)) {
      base::unique(base::as.character(found_data$enzyme_name))
    } else {
      character()
    }
    missing_names <- base::setdiff(target_names, found_names)
    missing_data <- if (base::length(missing_names)) {
      base::data.frame(
        enzyme_name = missing_names,
        enz_type = NA_character_,
        rec_seq = NA_character_,
        genbank_accession = NA_character_,
        rebase_locus = NA_character_,
        protein_id = NA_character_,
        uniprot_id = NA_character_,
        .raw_sequence = NA_character_,
        .raw_occurrence = 0L,
        .raw_org_name = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      base::data.frame()
    }
    if (base::nrow(found_data)) {
      new_data <- found_data
      new_data$.raw_occurrence <- stats::ave(
        base::seq_len(base::nrow(new_data)),
        new_data$enzyme_name,
        FUN = base::seq_along
      )
    } else {
      new_data <- base::data.frame()
    }
    new_data <- base::rbind(new_data, missing_data)
    if (base::nrow(new_data)) {
      cache_data <- base::rbind(cache_data, new_data)
      tryCatch(
        base::saveRDS(list(signature = signature, data = cache_data), cache_path),
        error = function(e) NULL
      )
    }
  }

  if (base::is.null(rebase_data)) {
    rebase_data <- .dnmb_rebasefinder_cached_rebase_data(cache_root)
  }
  normalize_value <- function(x) {
    x <- base::tolower(base::trimws(base::as.character(x)))
    x[base::is.na(x) | !base::nzchar(x)] <- NA_character_
    x
  }
  normalize_sequence <- function(x) {
    x <- base::toupper(base::gsub("\\s+", "", base::as.character(x), perl = TRUE))
    x[base::is.na(x) | !base::nzchar(x)] <- NA_character_
    x
  }
  value_equal <- function(x, target, normalizer = normalize_value) {
    x <- normalizer(x)
    target <- normalizer(target)[[1]]
    !base::is.na(target) & !base::is.na(x) & x == target
  }
  rebase_names <- if (base::is.data.frame(rebase_data) &&
                      "enzyme_name" %in% base::names(rebase_data)) {
    base::as.character(rebase_data$enzyme_name)
  } else {
    character()
  }
  resolved <- base::lapply(base::seq_along(subject_ids), function(i) {
    subject <- subject_ids[[i]]
    enzyme <- base_names[[i]]
    choices <- cache_data[cache_data$enzyme_name == enzyme, , drop = FALSE]
    suffix <- if (base::grepl("_[0-9]+$", subject, perl = TRUE)) {
      suppressWarnings(base::as.integer(base::sub("^.*_([0-9]+)$", "\\1", subject, perl = TRUE)))
    } else {
      NA_integer_
    }
    reference_idx <- if (!base::is.na(suffix) && suffix <= base::length(rebase_names) &&
                         rebase_names[[suffix]] == enzyme) suffix else NA_integer_
    chosen <- 1L
    if (base::nrow(choices) && !base::is.na(reference_idx)) {
      target <- rebase_data[reference_idx, , drop = FALSE]
      score <- base::rep(0, base::nrow(choices))
      if ("sequence" %in% base::names(target)) {
        score <- score + 1000 * value_equal(
          choices$.raw_sequence, target$sequence, normalize_sequence
        )
      }
      if ("rec_seq" %in% base::names(target)) {
        score <- score + 100 * value_equal(choices$rec_seq, target$rec_seq)
      }
      if ("enz_type" %in% base::names(target)) {
        score <- score + 10 * value_equal(choices$enz_type, target$enz_type)
      }
      if ("org_name" %in% base::names(target)) {
        score <- score + 5 * value_equal(choices$.raw_org_name, target$org_name)
      }
      best <- base::which(score == base::max(score))

      # Exact duplicate fingerprints retain their raw order under DefenseViz's
      # stable sort. Use that rank only after sequence and metadata matching.
      rebase_candidates <- base::which(rebase_names == enzyme)
      equivalent <- rebase_candidates
      raw_equivalent <- best
      fields <- list(
        sequence = list(raw = ".raw_sequence", normalizer = normalize_sequence),
        rec_seq = list(raw = "rec_seq", normalizer = normalize_value),
        enz_type = list(raw = "enz_type", normalizer = normalize_value),
        org_name = list(raw = ".raw_org_name", normalizer = normalize_value)
      )
      for (field in base::names(fields)) {
        if (!field %in% base::names(rebase_data)) next
        target_value <- rebase_data[[field]][[reference_idx]]
        normalizer <- fields[[field]]$normalizer
        normalized_target <- normalizer(target_value)[[1]]
        if (base::is.na(normalized_target)) next
        equivalent <- equivalent[
          value_equal(rebase_data[[field]][equivalent], target_value, normalizer)
        ]
        raw_field <- fields[[field]]$raw
        raw_equivalent <- raw_equivalent[
          value_equal(choices[[raw_field]][raw_equivalent], target_value, normalizer)
        ]
      }
      if (base::length(raw_equivalent)) {
        target_rank <- base::match(reference_idx, equivalent)
        chosen <- raw_equivalent[[base::min(
          if (base::is.na(target_rank)) 1L else target_rank,
          base::length(raw_equivalent)
        )]]
      } else if (base::length(best)) {
        chosen <- best[[1]]
      }
    }
    if (base::nrow(choices)) {
      row <- choices[chosen, base::setdiff(base::names(empty), "subject_id"), drop = FALSE]
    } else {
      row <- empty[0, base::setdiff(base::names(empty), "subject_id"), drop = FALSE]
      row[1, ] <- NA
      row$enzyme_name <- enzyme
    }
    base::data.frame(subject_id = subject, row, stringsAsFactors = FALSE)
  })
  base::do.call(base::rbind, resolved)
}

.dnmb_rebasefinder_reference_donor_metadata <- function(donor_names,
                                                        preferred_accessions = NA_character_,
                                                        cache_root = NULL,
                                                        rebase_data = NULL) {
  donor_names <- base::as.character(donor_names)
  preferred_accessions <- base::rep_len(base::as.character(preferred_accessions), base::length(donor_names))
  out <- base::data.frame(
    subject_id = base::rep(NA_character_, base::length(donor_names)),
    enzyme_name = .dnmb_rebasefinder_clean_rebase_subject(donor_names),
    enz_type = base::rep(NA_character_, base::length(donor_names)),
    rec_seq = base::rep(NA_character_, base::length(donor_names)),
    genbank_accession = base::rep(NA_character_, base::length(donor_names)),
    rebase_locus = base::rep(NA_character_, base::length(donor_names)),
    protein_id = base::rep(NA_character_, base::length(donor_names)),
    uniprot_id = base::rep(NA_character_, base::length(donor_names)),
    stringsAsFactors = FALSE
  )
  valid <- !base::is.na(donor_names) & base::nzchar(donor_names)
  if (!base::any(valid)) return(out)

  if (base::is.null(rebase_data)) {
    rebase_data <- .dnmb_rebasefinder_cached_rebase_data(cache_root)
  }
  request_ids <- base::lapply(base::seq_along(donor_names), function(i) {
    if (!valid[[i]]) return(character())
    enzyme <- .dnmb_rebasefinder_clean_rebase_subject(donor_names[[i]])
    positions <- if (base::is.data.frame(rebase_data) && "enzyme_name" %in% base::names(rebase_data)) {
      base::which(base::as.character(rebase_data$enzyme_name) == enzyme)
    } else {
      integer()
    }
    if (base::length(positions) > 1L) base::paste0(enzyme, "_", positions) else donor_names[[i]]
  })
  all_ids <- base::unique(base::unlist(request_ids, use.names = FALSE))
  metadata <- .dnmb_rebasefinder_reference_header_metadata(
    all_ids,
    cache_root = cache_root,
    rebase_data = rebase_data
  )

  for (i in base::which(valid)) {
    choices <- metadata[base::match(request_ids[[i]], metadata$subject_id), , drop = FALSE]
    choices <- choices[!base::is.na(choices$subject_id), , drop = FALSE]
    if (!base::nrow(choices)) next
    preferred <- preferred_accessions[[i]]
    preferred_match <- if (!base::is.na(preferred) && base::nzchar(preferred)) {
      base::which(.dnmb_rebasefinder_accession_key(choices$genbank_accession) ==
                    .dnmb_rebasefinder_accession_key(preferred))
    } else {
      integer()
    }
    if (base::length(preferred_match)) {
      out[i, ] <- choices[preferred_match[[1]], , drop = FALSE]
      next
    }
    if (base::nrow(choices) == 1L) {
      out[i, ] <- choices[1, , drop = FALSE]
      next
    }
    accessions <- base::unique(choices$genbank_accession[
      !base::is.na(choices$genbank_accession) & base::nzchar(choices$genbank_accession)
    ])
    if (base::length(accessions) == 1L) {
      out$genbank_accession[[i]] <- accessions[[1]]
    }
  }
  out
}

.dnmb_rebasefinder_parse_rebase_pacbio_index <- function(path) {
  if (!base::file.exists(path)) return(base::data.frame())
  html <- base::paste(base::readLines(path, warn = FALSE), collapse = "\n")
  blocks <- base::strsplit(html, "(?i)<tr[^>]*>", perl = TRUE)[[1]]
  rows <- base::lapply(blocks, function(block) {
    id_match <- base::regmatches(block, base::regexec("pacbioget\\?([0-9]+)", block, perl = TRUE))[[1]]
    if (base::length(id_match) < 2L) return(NULL)
    org_match <- base::regmatches(
      block,
      base::regexec("pacbioget\\?[0-9]+[^>]*>([^<]+)</[Aa]>", block, perl = TRUE)
    )[[1]]
    if (base::length(org_match) < 2L) return(NULL)
    pacbio_link <- base::regexpr("pacbioget\\?[0-9]+[^>]*>[^<]+</[Aa]>", block, perl = TRUE)
    after <- if (pacbio_link[[1]] > 0L) {
      base::substr(
        block,
        pacbio_link[[1]] + base::attr(pacbio_link, "match.length")[[1]],
        base::nchar(block)
      )
    } else {
      block
    }
    plain <- .dnmb_rebasefinder_html_text(after)
    accessions <- base::regmatches(
      plain,
      base::gregexpr("(?:NZ_)?[A-Z]{1,6}_?[A-Z0-9]*[0-9]{5,}(?:\\.[0-9]+)?", plain, perl = TRUE)
    )[[1]]
    accessions <- base::unique(accessions[base::nzchar(accessions)])
    if (!base::length(accessions)) return(NULL)
    base::data.frame(
      nucleotide_accession = accessions,
      organism = .dnmb_rebasefinder_html_text(org_match[[2]]),
      rebase_org_id = id_match[[2]],
      pacbio_url = base::paste0("https://rebase.neb.com/cgi-bin/pacbioget?", id_match[[2]]),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(base::is.null), rows)
  if (!base::length(rows)) return(base::data.frame())
  out <- base::do.call(base::rbind, rows)
  out$.accession_key <- .dnmb_rebasefinder_accession_key(out$nucleotide_accession)
  out <- out[!base::duplicated(out$.accession_key), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_rebasefinder_download_rebase_pacbio_index <- function(cache_root = NULL,
                                                            force = FALSE,
                                                            max_age_days = 31L) {
  cache_path <- base::file.path(
    .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE),
    "rebase_pacbio_index.rds"
  )
  if (!isTRUE(force) && base::file.exists(cache_path)) {
    age <- base::as.numeric(base::difftime(base::Sys.time(), base::file.info(cache_path)$mtime, units = "days"))
    cached <- tryCatch(base::readRDS(cache_path), error = function(e) NULL)
    if (base::is.data.frame(cached) && base::nrow(cached) && age <= max_age_days) return(cached)
  }
  tmp <- base::tempfile(fileext = ".html")
  on.exit(base::unlink(tmp, force = TRUE), add = TRUE)
  old_timeout <- base::getOption("timeout")
  on.exit(base::options(timeout = old_timeout), add = TRUE)
  base::options(timeout = base::max(60, old_timeout))
  ok <- tryCatch({
    utils::download.file(
      "https://rebase.neb.com/cgi-bin/pblist",
      tmp,
      quiet = TRUE,
      mode = "wb",
      method = "libcurl"
    )
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !base::file.exists(tmp) || base::file.info(tmp)$size < 1000) {
    if (base::file.exists(cache_path)) return(tryCatch(base::readRDS(cache_path), error = function(e) base::data.frame()))
    return(base::data.frame())
  }
  index <- .dnmb_rebasefinder_parse_rebase_pacbio_index(tmp)
  if (base::nrow(index)) tryCatch(base::saveRDS(index, cache_path), error = function(e) NULL)
  index
}

.dnmb_rebasefinder_match_rebase_pacbio <- function(accessions, index) {
  accessions <- base::as.character(accessions)
  n <- base::length(accessions)
  out <- base::data.frame(
    query_accession = accessions,
    pacbio_available = base::rep(FALSE, n),
    pacbio_organism = base::rep(NA_character_, n),
    pacbio_url = base::rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
  if (!base::is.data.frame(index) || !base::nrow(index)) return(out)
  if (!".accession_key" %in% base::names(index)) {
    index$.accession_key <- .dnmb_rebasefinder_accession_key(index$nucleotide_accession)
  }
  idx <- base::match(.dnmb_rebasefinder_accession_key(accessions), index$.accession_key)
  found <- !base::is.na(idx)
  out$pacbio_available[found] <- TRUE
  out$pacbio_organism[found] <- base::as.character(index$organism[idx[found]])
  out$pacbio_url[found] <- base::as.character(index$pacbio_url[idx[found]])
  out
}

.dnmb_rebasefinder_rebase_link <- function(enzyme) {
  enzyme <- base::as.character(enzyme)
  out <- base::rep(NA_character_, base::length(enzyme))
  keep <- !base::is.na(enzyme) & base::nzchar(enzyme) &
    !base::grepl("_context:", enzyme, fixed = TRUE)
  out[keep] <- base::paste0(
    "https://rebase.neb.com/rebase/enz/",
    utils::URLencode(.dnmb_rebasefinder_clean_rebase_subject(enzyme[keep]), reserved = TRUE),
    ".html"
  )
  out
}

.dnmb_rebasefinder_ncbi_link <- function(accession) {
  accession <- base::as.character(accession)
  out <- base::rep(NA_character_, base::length(accession))
  keep <- !base::is.na(accession) & base::nzchar(accession)
  out[keep] <- base::paste0(
    "https://www.ncbi.nlm.nih.gov/nuccore/",
    utils::URLencode(accession[keep], reserved = TRUE)
  )
  out
}

.dnmb_rebasefinder_interval_union_length <- function(start, end) {
  start <- suppressWarnings(base::as.numeric(start))
  end <- suppressWarnings(base::as.numeric(end))
  valid <- !base::is.na(start) & !base::is.na(end)
  if (!base::any(valid)) return(NA_real_)
  intervals <- base::cbind(base::pmin(start[valid], end[valid]), base::pmax(start[valid], end[valid]))
  intervals <- intervals[base::order(intervals[, 1], intervals[, 2]), , drop = FALSE]
  total <- 0
  current_start <- intervals[1, 1]
  current_end <- intervals[1, 2]
  if (base::nrow(intervals) > 1L) {
    for (i in 2:base::nrow(intervals)) {
      if (intervals[i, 1] <= current_end + 1) {
        current_end <- base::max(current_end, intervals[i, 2])
      } else {
        total <- total + current_end - current_start + 1
        current_start <- intervals[i, 1]
        current_end <- intervals[i, 2]
      }
    }
  }
  total + current_end - current_start + 1
}

.dnmb_rebasefinder_collapse_blast_hsps <- function(blast_tbl) {
  if (!base::is.data.frame(blast_tbl) || !base::nrow(blast_tbl)) return(blast_tbl)
  key <- base::paste(blast_tbl$query_id, blast_tbl$rebase_enzyme, sep = "\r")
  groups <- base::split(base::seq_len(base::nrow(blast_tbl)), key)
  rows <- base::lapply(groups, function(idx) {
    if (base::length(idx) == 1L) {
      row <- blast_tbl[idx, , drop = FALSE]
      row$dnmb_aligned_query_length <- if (!base::is.na(row$qstart) && !base::is.na(row$qend)) {
        base::abs(row$qend - row$qstart) + 1
      } else row$length
      row$dnmb_aligned_subject_length <- if (!base::is.na(row$sstart) && !base::is.na(row$send)) {
        base::abs(row$send - row$sstart) + 1
      } else row$length
      row$dnmb_hsp_count <- 1L
      return(row)
    }
    complete_intervals <- base::all(!base::is.na(blast_tbl$qstart[idx])) &
      base::all(!base::is.na(blast_tbl$qend[idx])) &
      base::all(!base::is.na(blast_tbl$sstart[idx])) &
      base::all(!base::is.na(blast_tbl$send[idx]))
    best <- idx[base::order(
      -base::ifelse(base::is.na(blast_tbl$bitscore[idx]), -Inf, blast_tbl$bitscore[idx]),
      base::ifelse(base::is.na(blast_tbl$evalue[idx]), Inf, blast_tbl$evalue[idx])
    )][1]
    row <- blast_tbl[best, , drop = FALSE]
    if (!complete_intervals) {
      row$dnmb_aligned_query_length <- row$length
      row$dnmb_aligned_subject_length <- row$length
      row$dnmb_hsp_count <- 1L
      return(row)
    }
    q_union <- .dnmb_rebasefinder_interval_union_length(blast_tbl$qstart[idx], blast_tbl$qend[idx])
    s_union <- .dnmb_rebasefinder_interval_union_length(blast_tbl$sstart[idx], blast_tbl$send[idx])
    weights <- suppressWarnings(base::as.numeric(blast_tbl$length[idx]))
    weights[base::is.na(weights) | weights <= 0] <- 1
    row$pct_identity <- if (base::all(base::is.na(blast_tbl$pct_identity[idx]))) {
      NA_real_
    } else {
      stats::weighted.mean(blast_tbl$pct_identity[idx], weights, na.rm = TRUE)
    }
    row$bitscore <- if (base::all(base::is.na(blast_tbl$bitscore[idx]))) NA_real_ else base::sum(blast_tbl$bitscore[idx], na.rm = TRUE)
    row$evalue <- if (base::all(base::is.na(blast_tbl$evalue[idx]))) NA_real_ else base::min(blast_tbl$evalue[idx], na.rm = TRUE)
    row$length <- base::min(q_union, s_union)
    row$qstart <- base::min(base::pmin(blast_tbl$qstart[idx], blast_tbl$qend[idx]))
    row$qend <- base::max(base::pmax(blast_tbl$qstart[idx], blast_tbl$qend[idx]))
    row$sstart <- base::min(base::pmin(blast_tbl$sstart[idx], blast_tbl$send[idx]))
    row$send <- base::max(base::pmax(blast_tbl$sstart[idx], blast_tbl$send[idx]))
    row$dnmb_aligned_query_length <- q_union
    row$dnmb_aligned_subject_length <- s_union
    row$dnmb_hsp_count <- base::length(idx)
    row
  })
  out <- base::do.call(base::rbind, rows)
  base::rownames(out) <- NULL
  out
}

.dnmb_rebasefinder_genbank_accessions <- function(genbank) {
  if (base::is.null(genbank)) return(character())
  paths <- base::as.character(genbank)
  paths <- paths[!base::is.na(paths) & base::nzchar(paths) & base::file.exists(paths)]
  if (!base::length(paths)) return(character())
  accessions <- base::unlist(base::lapply(paths, function(path) {
    lines <- tryCatch(base::readLines(path, warn = FALSE), error = function(e) character())
    extract_field <- function(field) {
      pattern <- base::paste0("^", field, "[[:space:]]+([^[:space:]]+)")
      matches <- base::regmatches(lines, base::regexec(pattern, lines, perl = TRUE))
      values <- base::vapply(matches, function(value) {
        if (base::length(value) >= 2L) value[[2L]] else NA_character_
      }, character(1))
      values <- values[
        !base::is.na(values) &
          base::grepl("^[A-Za-z]{1,6}_?[A-Za-z0-9]+(?:\\.[0-9]+)?$", values, perl = TRUE)
      ]
      base::unique(values)
    }
    version <- extract_field("VERSION")
    accession <- extract_field("ACCESSION")
    if (base::length(version)) version else accession
  }), use.names = FALSE)
  base::unique(accessions[!base::is.na(accessions) & base::nzchar(accessions)])
}

.dnmb_rebasefinder_standardize_blast <- function(blast_tbl) {
  if (!base::is.data.frame(blast_tbl) || !base::nrow(blast_tbl)) return(base::data.frame())
  blast_tbl <- base::as.data.frame(blast_tbl, stringsAsFactors = FALSE)
  if (!"query_id" %in% base::names(blast_tbl) && "qseqid" %in% base::names(blast_tbl)) blast_tbl$query_id <- blast_tbl$qseqid
  if (!"subject_id" %in% base::names(blast_tbl) && "rebase_enzyme" %in% base::names(blast_tbl)) blast_tbl$subject_id <- blast_tbl$rebase_enzyme
  if (!"rebase_enzyme" %in% base::names(blast_tbl) && "subject_id" %in% base::names(blast_tbl)) blast_tbl$rebase_enzyme <- blast_tbl$subject_id
  if (!"pct_identity" %in% base::names(blast_tbl) && "pident" %in% base::names(blast_tbl)) blast_tbl$pct_identity <- blast_tbl$pident / 100
  n <- base::nrow(blast_tbl)
  if (!"query_id" %in% base::names(blast_tbl)) blast_tbl$query_id <- base::rep(NA_character_, n)
  if (!"rebase_enzyme" %in% base::names(blast_tbl)) blast_tbl$rebase_enzyme <- base::rep(NA_character_, n)
  for (col in c("pct_identity", "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send")) {
    if (!col %in% base::names(blast_tbl)) blast_tbl[[col]] <- base::rep(NA_real_, n)
  }
  .dnmb_rebasefinder_collapse_blast_hsps(blast_tbl)
}

.dnmb_rebasefinder_enrich_blast_metadata <- function(hits,
                                                      blast_tbl = NULL,
                                                      rebase_data = NULL,
                                                      cache_root = NULL,
                                                      query_genbank = NULL,
                                                      min_identity = 0.10,
                                                      min_length = 50,
                                                      max_evalue = 0.001,
                                                      verbose = FALSE) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  blast_tbl <- .dnmb_rebasefinder_standardize_blast(blast_tbl)
  if (base::nrow(blast_tbl)) {
    keep_blast <- !base::is.na(blast_tbl$pct_identity) &
      blast_tbl$pct_identity >= min_identity
    if ("length" %in% base::names(blast_tbl)) {
      keep_blast <- keep_blast & !base::is.na(blast_tbl$length) &
        blast_tbl$length >= min_length
    }
    if ("evalue" %in% base::names(blast_tbl)) {
      keep_blast <- keep_blast & !base::is.na(blast_tbl$evalue) &
        blast_tbl$evalue <= max_evalue
    }
    blast_tbl <- blast_tbl[keep_blast, , drop = FALSE]
  }
  bairoch <- tryCatch(
    .dnmb_rebasefinder_download_bairoch(cache_root = cache_root),
    error = function(e) NULL
  )

  supplemental_conflict <- if ("supplemental_blast_context_conflict" %in% base::names(hits)) {
    !base::is.na(hits$supplemental_blast_context_conflict) &
      hits$supplemental_blast_context_conflict
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  fallback_subject <- if ("supplemental_blast_match" %in% base::names(hits)) {
    base::ifelse(
      !supplemental_conflict & !base::is.na(hits$supplemental_blast_match) &
        base::nzchar(hits$supplemental_blast_match),
      base::as.character(hits$supplemental_blast_match),
      base::as.character(hits$hit_label)
    )
  } else {
    base::as.character(hits$hit_label)
  }
  selected_subject <- if ("curated_blast_resolved_subject_id" %in% base::names(hits)) {
    base::as.character(hits$curated_blast_resolved_subject_id)
  } else if ("curated_blast_subject_id" %in% base::names(hits)) {
    base::as.character(hits$curated_blast_subject_id)
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  selected_promoted <- if ("curated_blast_promoted" %in% base::names(hits)) {
    hits$curated_blast_promoted %in% TRUE
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  has_selected <- selected_promoted & !base::is.na(selected_subject) & base::nzchar(selected_subject)
  primary_subject <- base::ifelse(has_selected, selected_subject, fallback_subject)
  subjects <- base::unique(c(primary_subject, if (base::nrow(blast_tbl)) blast_tbl$rebase_enzyme else character()))
  headers <- .dnmb_rebasefinder_reference_header_metadata(
    subjects,
    cache_root = cache_root,
    rebase_data = rebase_data
  )
  header_rec <- headers[.dnmb_rebasefinder_valid_recognition(headers$rec_seq), c("enzyme_name", "rec_seq"), drop = FALSE]
  if (base::nrow(header_rec)) {
    header_rec$is_gold_standard <- TRUE
    if (base::is.null(rebase_data) || !base::is.data.frame(rebase_data)) {
      recognition_data <- header_rec
    } else {
      recognition_data <- base::rbind(
        base::data.frame(enzyme_name = rebase_data$enzyme_name, rec_seq = rebase_data$rec_seq,
                         is_gold_standard = if ("is_gold_standard" %in% base::names(rebase_data)) rebase_data$is_gold_standard else FALSE),
        base::data.frame(enzyme_name = header_rec$enzyme_name, rec_seq = header_rec$rec_seq, is_gold_standard = TRUE)
      )
    }
  } else {
    recognition_data <- rebase_data
  }
  hits <- .dnmb_rebasefinder_enrich_recognition(hits, recognition_data, bairoch)

  ensure <- function(name, value) {
    if (!name %in% base::names(hits)) hits[[name]] <<- base::rep(value, base::nrow(hits))
  }
  ensure("reference_rec_seq", NA_character_)
  ensure("recognition_match_identity", NA_real_)
  ensure("recognition_match_evalue", NA_real_)
  ensure("recognition_match_bitscore", NA_real_)
  ensure("recognition_match_length", NA_real_)
  ensure("recognition_match_query_coverage", NA_real_)
  ensure("recognition_match_reference_coverage", NA_real_)
  ensure("recognition_match_alignment_quality", NA_character_)
  ensure("recognition_donor", NA_character_)
  ensure("recognition_match_subject_id", NA_character_)
  ensure("recognition_transfer_eligible", FALSE)
  ensure("recognition_transferred", FALSE)

  lookup <- .dnmb_rebasefinder_recognition_lookup(recognition_data, bairoch)
  if (base::nrow(blast_tbl) && base::nrow(lookup)) {
    candidates <- list()
    for (i in base::seq_len(base::nrow(blast_tbl))) {
      subject_id <- base::as.character(blast_tbl$rebase_enzyme[[i]])
      subject <- .dnmb_rebasefinder_clean_rebase_subject(subject_id)
      query <- .dnmb_module_clean_annotation_key(blast_tbl$query_id[[i]])
      hit_i <- base::match(query, .dnmb_module_clean_annotation_key(hits$query))
      selected <- if (!base::is.na(hit_i)) primary_subject[[hit_i]] else NA_character_
      same_primary_record_name <- !base::is.na(selected) & base::nzchar(selected) &
        subject == .dnmb_rebasefinder_clean_rebase_subject(selected)
      primary_related <- !base::is.na(selected) & base::nzchar(selected) & (
        subject_id == selected |
          same_primary_record_name |
          .dnmb_rebasefinder_system_key(subject) == .dnmb_rebasefinder_system_key(selected)
      )
      resolved_subject_id <- if (same_primary_record_name) selected else subject_id
      duplicate_subject <- base::is.data.frame(rebase_data) &&
        "enzyme_name" %in% base::names(rebase_data) &&
        base::sum(base::as.character(rebase_data$enzyme_name) == subject, na.rm = TRUE) > 1L
      resolved_duplicate <- !duplicate_subject || base::grepl("_[0-9]+$", resolved_subject_id, perl = TRUE)
      if (duplicate_subject && !resolved_duplicate) next
      header_idx <- base::match(resolved_subject_id, headers$subject_id)
      header_has_rec <- !base::is.na(header_idx) &&
        .dnmb_rebasefinder_valid_recognition(headers$rec_seq[[header_idx]])
      if (header_has_rec) {
        donor_name <- subject
        recognition_site <- headers$rec_seq[[header_idx]]
        recognition_source <- "rebase_protein_header_blast_exact_record"
      } else {
        exact <- base::which(lookup$enzyme_name == subject)
        relation <- "exact"
        donor <- exact
        if (!base::length(donor)) {
          donor <- base::which(lookup$system_key == .dnmb_rebasefinder_system_key(subject))
          relation <- "cognate_system"
        }
        if (!base::length(donor)) next
        donor <- donor[[1]]
        donor_name <- lookup$enzyme_name[[donor]]
        recognition_site <- lookup$rec_seq[[donor]]
        recognition_source <- base::paste0(lookup$recognition_source[[donor]], "_blast_", relation)
      }
      candidates[[base::length(candidates) + 1L]] <- base::data.frame(
        query = query,
        subject_id = resolved_subject_id,
        subject = subject,
        primary_related = primary_related,
        donor = donor_name,
        rec_seq = recognition_site,
        source = recognition_source,
        exact_record = header_has_rec && resolved_duplicate,
        identity = blast_tbl$pct_identity[[i]],
        evalue = if ("evalue" %in% base::names(blast_tbl)) blast_tbl$evalue[[i]] else NA_real_,
        bitscore = if ("bitscore" %in% base::names(blast_tbl)) blast_tbl$bitscore[[i]] else NA_real_,
        length = if ("length" %in% base::names(blast_tbl)) blast_tbl$length[[i]] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    if (base::length(candidates)) {
      candidates <- base::do.call(base::rbind, candidates)
      hit_idx <- base::match(candidates$query, .dnmb_module_clean_annotation_key(hits$query))
      query_len <- if ("query_aa_len" %in% base::names(hits)) {
        suppressWarnings(base::as.numeric(hits$query_aa_len[hit_idx]))
      } else {
        base::rep(NA_real_, base::nrow(candidates))
      }
      refs <- .dnmb_rebasefinder_reference_rows(
        candidates$subject_id,
        rebase_data,
        alignment_length = candidates$length
      )
      candidates$query_coverage <- base::pmin(1, candidates$length / query_len)
      candidates$reference_coverage <- base::pmin(1, candidates$length / refs$reference_aa_len)
      candidates$alignment_quality <- .dnmb_rebasefinder_alignment_quality(
        candidates$identity,
        candidates$query_coverage,
        candidates$reference_coverage,
        candidates$evalue
      )
      quality_rank <- base::match(candidates$alignment_quality, c("weak", "moderate", "strong"))
      min_coverage <- base::pmin(candidates$query_coverage, candidates$reference_coverage)
      candidates <- candidates[base::order(
        candidates$query,
        !candidates$primary_related,
        -quality_rank,
        -base::ifelse(base::is.na(min_coverage), -Inf, min_coverage),
        -base::ifelse(base::is.na(candidates$bitscore), -Inf, candidates$bitscore),
        -candidates$identity
      ), , drop = FALSE]
      candidates <- candidates[!base::duplicated(candidates$query), , drop = FALSE]
      idx <- base::match(hits$query, candidates$query)
      found <- !base::is.na(idx)
      if (base::any(found)) {
        cand <- candidates[idx[found], , drop = FALSE]
        found_idx <- base::which(found)
        hits$recognition_match[found] <- cand$subject
        hits$recognition_match_subject_id[found] <- cand$subject_id
        replace_recognition <- cand$exact_record %in% TRUE |
          !.dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq[found_idx])
        if (base::any(replace_recognition)) {
          replace_idx <- found_idx[replace_recognition]
          hits$recognition_donor[replace_idx] <- cand$donor[replace_recognition]
          hits$reference_rec_seq[replace_idx] <- cand$rec_seq[replace_recognition]
          hits$recognition_source[replace_idx] <- cand$source[replace_recognition]
        }
        hits$recognition_match_identity[found] <- cand$identity
        hits$recognition_match_evalue[found] <- cand$evalue
        hits$recognition_match_bitscore[found] <- cand$bitscore
        hits$recognition_match_length[found] <- cand$length
        hits$recognition_match_query_coverage[found] <- cand$query_coverage
        hits$recognition_match_reference_coverage[found] <- cand$reference_coverage
        hits$recognition_match_alignment_quality[found] <- cand$alignment_quality
        primary_clean <- .dnmb_rebasefinder_clean_rebase_subject(primary_subject[found])
        related_to_primary <- !supplemental_conflict[found] & cand$primary_related &
          cand$alignment_quality %in% c("strong", "moderate")
        # Recognition sites from homologs remain reference-only. Even a close
        # match can change specificity through a small target-recognition region.
        hits$recognition_transfer_eligible[found] <- FALSE
        hits$recognition_transferred[found] <- FALSE
        if (base::any(related_to_primary)) {
          hits$recognition_source[base::which(found)[related_to_primary]] <- base::paste0(
            hits$recognition_source[base::which(found)[related_to_primary]],
            "_primary_or_cognate_reference"
          )
        }
      }
    }
  }
  direct_rec <- .dnmb_rebasefinder_valid_recognition(hits$rec_seq) &
    !.dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq)
  hits$reference_rec_seq[direct_rec] <- hits$rec_seq[direct_rec]

  primary_meta <- headers[base::match(primary_subject, headers$subject_id), , drop = FALSE]
  recognition_meta <- headers[base::match(hits$recognition_match_subject_id, headers$subject_id), , drop = FALSE]
  # Header metadata are indexed by the original BLAST subject; clean-name fallback
  # covers unsuffixed exact records and cognate enrichment.
  missing_rec_meta <- base::is.na(recognition_meta$subject_id) & !base::is.na(hits$recognition_match)
  if (base::any(missing_rec_meta)) {
    clean_idx <- base::match(
      .dnmb_rebasefinder_clean_rebase_subject(hits$recognition_match[missing_rec_meta]),
      headers$enzyme_name
    )
    recognition_meta[missing_rec_meta, ] <- headers[clean_idx, , drop = FALSE]
  }
  donor_meta <- .dnmb_rebasefinder_reference_donor_metadata(
    base::rep(NA_character_, base::nrow(hits)),
    cache_root = cache_root,
    rebase_data = rebase_data
  )
  same_record_donor <- !base::is.na(hits$recognition_donor) &
    !base::is.na(hits$recognition_match) &
    .dnmb_rebasefinder_clean_rebase_subject(hits$recognition_donor) ==
      .dnmb_rebasefinder_clean_rebase_subject(hits$recognition_match)
  if (base::any(same_record_donor)) {
    donor_meta[same_record_donor, ] <- recognition_meta[same_record_donor, , drop = FALSE]
  }
  cognate_donor <- !same_record_donor & !base::is.na(hits$recognition_donor) &
    base::nzchar(hits$recognition_donor)
  if (base::any(cognate_donor)) {
    donor_meta[cognate_donor, ] <- .dnmb_rebasefinder_reference_donor_metadata(
      hits$recognition_donor[cognate_donor],
      preferred_accessions = recognition_meta$genbank_accession[cognate_donor],
      cache_root = cache_root,
      rebase_data = rebase_data
    )
  }
  hits$rebase_reference_accession <- primary_meta$genbank_accession
  hits$rebase_reference_locus <- primary_meta$rebase_locus
  hits$rebase_reference_protein_id <- primary_meta$protein_id
  hits$rebase_match_link <- .dnmb_rebasefinder_rebase_link(primary_subject)
  hits$rebase_reference_ncbi_link <- .dnmb_rebasefinder_ncbi_link(primary_meta$genbank_accession)
  hits$recognition_match_reference_accession <- recognition_meta$genbank_accession
  hits$recognition_match_reference_locus <- recognition_meta$rebase_locus
  hits$recognition_match_reference_protein_id <- recognition_meta$protein_id
  hits$recognition_match_rebase_link <- .dnmb_rebasefinder_rebase_link(hits$recognition_match)
  hits$recognition_match_reference_ncbi_link <- .dnmb_rebasefinder_ncbi_link(recognition_meta$genbank_accession)
  hits$recognition_donor_accession <- donor_meta$genbank_accession
  hits$recognition_donor_locus <- donor_meta$rebase_locus
  hits$recognition_donor_protein_id <- donor_meta$protein_id
  hits$recognition_donor_rebase_link <- .dnmb_rebasefinder_rebase_link(hits$recognition_donor)
  hits$recognition_donor_ncbi_link <- .dnmb_rebasefinder_ncbi_link(donor_meta$genbank_accession)
  # Backward-compatible reference fields now point to the record that supplied
  # the recognition sequence, not merely the matched homolog.
  hits$recognition_reference_accession <- donor_meta$genbank_accession
  hits$recognition_reference_ncbi_link <- hits$recognition_donor_ncbi_link

  query_accessions <- .dnmb_rebasefinder_genbank_accessions(query_genbank)
  accessions <- base::unique(c(
    primary_meta$genbank_accession,
    recognition_meta$genbank_accession,
    donor_meta$genbank_accession,
    query_accessions
  ))
  accessions <- accessions[!base::is.na(accessions) & base::nzchar(accessions)]
  pacbio_index <- if (base::length(accessions)) {
    tryCatch(.dnmb_rebasefinder_download_rebase_pacbio_index(cache_root), error = function(e) base::data.frame())
  } else {
    base::data.frame()
  }
  primary_pb <- .dnmb_rebasefinder_match_rebase_pacbio(primary_meta$genbank_accession, pacbio_index)
  rec_pb <- .dnmb_rebasefinder_match_rebase_pacbio(recognition_meta$genbank_accession, pacbio_index)
  donor_pb <- .dnmb_rebasefinder_match_rebase_pacbio(donor_meta$genbank_accession, pacbio_index)
  query_pb <- .dnmb_rebasefinder_match_rebase_pacbio(query_accessions, pacbio_index)
  hits$rebase_match_pacbio_available <- primary_pb$pacbio_available
  hits$rebase_match_pacbio_organism <- primary_pb$pacbio_organism
  hits$rebase_match_pacbio_link <- primary_pb$pacbio_url
  hits$recognition_match_pacbio_available <- rec_pb$pacbio_available
  hits$recognition_match_pacbio_organism <- rec_pb$pacbio_organism
  hits$recognition_match_pacbio_link <- rec_pb$pacbio_url
  hits$recognition_donor_pacbio_available <- donor_pb$pacbio_available
  hits$recognition_donor_pacbio_organism <- donor_pb$pacbio_organism
  hits$recognition_donor_pacbio_link <- donor_pb$pacbio_url
  query_found <- if (base::nrow(query_pb)) base::which(query_pb$pacbio_available %in% TRUE) else integer()
  query_pick <- if (base::length(query_found)) query_found[[1]] else NA_integer_
  hits$query_strain_accessions <- base::paste(query_accessions, collapse = ";")
  hits$query_strain_pacbio_available <- !base::is.na(query_pick)
  hits$query_strain_pacbio_organism <- if (!base::is.na(query_pick)) query_pb$pacbio_organism[[query_pick]] else NA_character_
  hits$query_strain_pacbio_link <- if (!base::is.na(query_pick)) query_pb$pacbio_url[[query_pick]] else NA_character_
  hits$query_strain_pacbio_index_link <- base::rep(
    "https://rebase.neb.com/cgi-bin/pblist",
    base::nrow(hits)
  )
  if (isTRUE(verbose)) {
    message(
      "[REBASEfinder] Recognition metadata: ",
      base::sum(.dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq)), " reference sites; ",
      base::sum(hits$recognition_match_pacbio_available %in% TRUE), " PacBio-linked references"
    )
  }
  hits
}
