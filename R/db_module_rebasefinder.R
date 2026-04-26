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
    inferred <- base::ifelse(base::grepl("^M[0-9]?\\.", hit_name), "M",
                base::ifelse(base::grepl("^R[0-9]?\\.", hit_name), "R",
                base::ifelse(base::grepl("^S[0-9]?\\.", hit_name), "S", NA_character_)))
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

  # Override DefenseViz cache dir → DNMB's db_modules/rebasefinder/cache
  # This works regardless of DefenseViz version (no env var dependency)
  rebase_cache <- .dnmb_rebasefinder_cache_dir(cache_root = cache_root, create = TRUE)
  base::dir.create(rebase_cache, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    cache_path <- rebase_cache  # capture value
    override_fn <- local({
      path <- cache_path
      function() {
        if (!base::dir.exists(path)) base::dir.create(path, recursive = TRUE, showWarnings = FALSE)
        path
      }
    })
    utils::assignInNamespace("get_rebase_cache_dir", override_fn, ns = "DefenseViz")
    if (isTRUE(verbose)) message("[REBASEfinder] REBASE cache: ", rebase_cache)
  }, error = function(e) {
    if (isTRUE(verbose)) message("[REBASEfinder] Could not override cache dir: ", conditionMessage(e))
  })

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
  output_table <- .dnmb_rebasefinder_output_table(genes = genes, hits = hits)

  n_hits <- base::nrow(hits)
  n_high <- base::sum(hits$typing_eligible, na.rm = TRUE)

  status <- dplyr::bind_rows(status, .dnmb_rebasefinder_status_row(
    "rebasefinder_run", "ok",
    base::sprintf("%d R-M candidates (%d high-confidence)", n_hits, n_high)
  ))

  list(
    ok = n_hits > 0L || base::nrow(rm_comprehensive) == 0L,
    status = status,
    hits = hits,
    output_table = output_table,
    pipeline_result = pipeline_result
  )
}
