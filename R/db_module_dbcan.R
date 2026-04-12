.dnmb_dbcan_module_name <- function() {
  "dbcan"
}

.dnmb_dbcan_default_version <- function() {
  "current"
}

.dnmb_dbcan_default_base_url <- function() {
  "https://pro.unl.edu/dbCAN2/download/Databases"
}

.dnmb_dbcan_release_page_url <- function() {
  "https://pro.unl.edu/dbCAN2/"
}

.dnmb_dbcan_status_row <- function(component, status, detail) {
  tibble::tibble(
    component = as.character(component),
    status = as.character(status),
    detail = as.character(detail)
  )
}

.dnmb_dbcan_empty_status <- function() {
  .dnmb_dbcan_status_row(character(), character(), character())
}

.dnmb_dbcan_trace <- function(path, text) {
  cat(text, "\n", file = path, append = TRUE)
}

.dnmb_dbcan_asset_layout <- function(module_dir) {
  hmm_path <- file.path(module_dir, "dbCAN.txt")
  mapping_path <- file.path(module_dir, "fam-substrate-mapping.tsv")
  list(
    module_dir = module_dir,
    hmm_path = hmm_path,
    mapping_path = mapping_path,
    hmmpress_files = paste0(hmm_path, c(".h3f", ".h3i", ".h3m", ".h3p"))
  )
}

.dnmb_dbcan_fallback_release_info <- function() {
  list(
    version = "V14",
    release_note = "Fallback release info",
    source = "fallback",
    dbcan_release_date = NA_character_
  )
}

.dnmb_dbcan_release_info <- function(page_url = .dnmb_dbcan_release_page_url()) {
  html <- character()
  curl_detection <- dnmb_detect_binary("curl", required = FALSE)
  if (isTRUE(curl_detection$found)) {
    run <- dnmb_run_external("curl", c("-k", "-L", "-sS", page_url), required = FALSE)
    if (isTRUE(run$ok) && length(run$stdout)) {
      html <- run$stdout
    }
  }
  if (!length(html)) {
    wget_detection <- dnmb_detect_binary("wget", required = FALSE)
    if (isTRUE(wget_detection$found)) {
      run <- dnmb_run_external("wget", c("--no-check-certificate", "-qO-", page_url), required = FALSE)
      if (isTRUE(run$ok) && length(run$stdout)) {
        html <- run$stdout
      }
    }
  }
  if (!length(html)) {
    return(.dnmb_dbcan_fallback_release_info())
  }

  pattern <- "([0-9]{1,2}/[0-9]{1,2}/[0-9]{4}).*dbCAN HMMdb v([0-9]+) is released \\(based on CAZyDB ([0-9/]+)\\)"
  match_line <- html[grepl(pattern, html, perl = TRUE)]
  if (!length(match_line)) {
    return(.dnmb_dbcan_fallback_release_info())
  }

  groups <- regmatches(match_line[[1]], regexec(pattern, match_line[[1]], perl = TRUE))[[1]]
  if (length(groups) < 4L) {
    return(.dnmb_dbcan_fallback_release_info())
  }

  list(
    version = paste0("V", groups[[3]]),
    release_note = trimws(gsub("<[^>]+>", "", match_line[[1]])),
    source = page_url,
    dbcan_release_date = groups[[2]],
    cazydb_release_date = groups[[4]]
  )
}

.dnmb_dbcan_aws_s3_base <- function() {
  # dbCAN's canonical AWS S3 distribution (see
  # dbcan/constants/databases_constants.py::AWS_S3_URL in the
  # `dbcan` python package). The legacy HTTP mirror under
  # `https://bcb.unl.edu/dbCAN2/download/...` is dead — every path
  # 302-redirects to a landing page and downloads save HTML instead
  # of the actual HMM/DIAMOND/TSV files.
  "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2_9-13-2025"
}

.dnmb_dbcan_asset_urls <- function(base_url = .dnmb_dbcan_default_base_url(),
                                   version = .dnmb_dbcan_default_version(),
                                   asset_urls = NULL,
                                   release_info = NULL) {
  if (!is.null(asset_urls)) {
    urls <- unlist(asset_urls, use.names = TRUE)
    urls <- stats::setNames(as.character(urls), names(urls))
    if (!length(urls) || is.null(names(urls)) || any(!nzchar(names(urls)))) {
      stop("`asset_urls` must be a named character vector or list.", call. = FALSE)
    }
    return(urls)
  }

  # Ignore the legacy versioned HTTP layout entirely — the S3 mirror
  # ships a single `dbCAN.hmm` (no per-version HMMdb-Vxx.txt split) and
  # the same `fam-substrate-mapping.tsv`. DNMB saves the HMM locally
  # as `dbCAN.txt` which hmmsearch still happily accepts.
  s3 <- .dnmb_dbcan_aws_s3_base()
  stats::setNames(
    c(
      paste0(s3, "/dbCAN.hmm"),
      paste0(s3, "/fam-substrate-mapping.tsv")
    ),
    c("dbCAN.txt", "fam-substrate-mapping.tsv")
  )
}

.dnmb_dbcan_is_hmm <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(FALSE)
  }
  first_lines <- tryCatch(readLines(path, n = 5L, warn = FALSE), error = function(...) character())
  any(grepl("^HMMER", first_lines))
}

.dnmb_dbcan_is_mapping <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(FALSE)
  }
  first_lines <- tryCatch(readLines(path, n = 5L, warn = FALSE), error = function(...) character())
  !any(grepl("^<!DOCTYPE html|^<html", first_lines, ignore.case = TRUE)) && any(grepl("\t", first_lines))
}

.dnmb_dbcan_hmmpress_summary <- function(hmmpress, args = character()) {
  paste(
    c(
      paste0("path=", if (isTRUE(hmmpress$found) && nzchar(hmmpress$path)) hmmpress$path else "missing"),
      if (!is.null(hmmpress$version)) paste0("version=", hmmpress$version),
      if (length(args)) paste0("args=", .dnmb_format_command("hmmpress", args))
    ),
    collapse = "; "
  )
}

.dnmb_dbcan_detect_hmmer <- function(binary) {
  detection <- dnmb_detect_binary(binary, required = FALSE)
  version <- NA_character_
  if (isTRUE(detection$found)) {
    run <- dnmb_run_external(binary, args = "-h", required = FALSE)
    version <- .dnmb_parse_tool_version(c(run$stdout, run$stderr))
  }
  list(
    binary = binary,
    path = detection$path,
    found = isTRUE(detection$found),
    version = version
  )
}

.dnmb_dbcan_prepare_hmm_db <- function(hmm_path, force = FALSE, trace_log = NULL) {
  pressed_files <- paste0(hmm_path, c(".h3f", ".h3i", ".h3m", ".h3p"))
  hmmpress <- .dnmb_dbcan_detect_hmmer("hmmpress")
  args <- c("-f", hmm_path)
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] hmmpress %s", Sys.time(), .dnmb_dbcan_hmmpress_summary(hmmpress, args)))
  }

  if (!isTRUE(force) && all(file.exists(pressed_files))) {
    return(list(
      ok = TRUE,
      status = "cached",
      detail = paste("dbCAN HMMER database files already exist", .dnmb_dbcan_hmmpress_summary(hmmpress, args), sep = "; "),
      hmm_path = hmm_path,
      hmmpress_files = pressed_files,
      command = NULL,
      hmmpress = hmmpress,
      args = args
    ))
  }

  run <- dnmb_run_external("hmmpress", args = args, required = FALSE)
  ok <- isTRUE(run$ok) && all(file.exists(pressed_files))
  detail <- if (ok) {
    paste(hmm_path, .dnmb_dbcan_hmmpress_summary(hmmpress, args), sep = "; ")
  } else {
    stderr_text <- .dnmb_compact_output(c(run$stderr, run$stdout), max_lines = 6L)
    paste(
      c(
        paste0("path=", if (isTRUE(hmmpress$found) && nzchar(hmmpress$path)) hmmpress$path else "missing"),
        paste0("version=", hmmpress$version %||% "unknown"),
        paste0("args=", .dnmb_format_command("hmmpress", args)),
        if (nzchar(stderr_text)) paste0("stderr=", stderr_text)
      ),
      collapse = "; "
    )
  }
  status <- if (ok) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed"
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] hmmpress_result status=%s detail=%s", Sys.time(), status, detail))
  }

  list(
    ok = ok,
    status = status,
    detail = detail,
    hmm_path = hmm_path,
    hmmpress_files = pressed_files[file.exists(pressed_files)],
    command = run,
    hmmpress = hmmpress,
    args = args
  )
}

.dnmb_dbcan_prepare_pul_blast_db <- function(pul_faa, trace_log = NULL) {
  output_files <- paste0(pul_faa, c(".phr", ".pin", ".psq"))
  if (all(file.exists(output_files))) {
    return(list(
      ok = TRUE,
      status = "cached",
      detail = pul_faa,
      files = output_files
    ))
  }

  stage_dir <- tempfile("dnmb-dbcan-pul-blast-")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged_faa <- file.path(stage_dir, "PUL.faa")
  file.copy(pul_faa, staged_faa, overwrite = TRUE)
  args <- c("-in", staged_faa, "-dbtype", "prot", "-out", staged_faa)
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] pul_makeblastdb args=%s", Sys.time(), .dnmb_format_command("makeblastdb", args)))
  }
  run <- dnmb_run_external("makeblastdb", args = args, required = FALSE)
  staged_files <- paste0(staged_faa, c(".phr", ".pin", ".psq"))
  ok <- isTRUE(run$ok) && all(file.exists(staged_files))
  if (ok) {
    file.copy(staged_files, output_files, overwrite = TRUE)
  }
  list(
    ok = ok,
    status = if (ok) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed",
    detail = if (ok) pul_faa else (run$error %||% pul_faa),
    files = output_files[file.exists(output_files)],
    command = run,
    args = args
  )
}

.dnmb_dbcan_manifest_diagnostics <- function(layout, prepare_result = NULL, release_info = NULL) {
  hmmpress <- prepare_result$hmmpress %||% .dnmb_dbcan_detect_hmmer("hmmpress")
  args <- prepare_result$args %||% c("-f", layout$hmm_path)
  info <- release_info %||% .dnmb_dbcan_fallback_release_info()

  list(
    hmmpress_path = if (isTRUE(hmmpress$found) && nzchar(hmmpress$path)) hmmpress$path else "",
    hmmpress_version = hmmpress$version %||% NA_character_,
    hmmpress_args = as.character(args),
    resolved_release_version = info$version %||% NA_character_,
    dbcan_release_date = info$dbcan_release_date %||% NA_character_,
    cazydb_release_date = info$cazydb_release_date %||% NA_character_,
    release_note = info$release_note %||% NA_character_
  )
}

.dnmb_dbcan_remote_asset_state <- function(urls, enabled = TRUE) {
  if (!isTRUE(enabled) || !length(urls)) {
    return(NULL)
  }

  rows <- lapply(names(urls), function(asset_name) {
    metadata <- .dnmb_remote_asset_metadata(urls[[asset_name]], insecure = TRUE)
    tibble::tibble(
      asset_name = asset_name,
      url = metadata$url %||% NA_character_,
      ok = isTRUE(metadata$ok),
      last_modified = metadata$last_modified %||% NA_character_,
      content_length = metadata$content_length %||% NA_real_,
      etag = metadata$etag %||% NA_character_
    )
  })

  dplyr::bind_rows(rows)
}

.dnmb_dbcan_remote_update_needed <- function(manifest, remote_state) {
  if (is.null(manifest) || is.null(remote_state) || !is.data.frame(remote_state) || !nrow(remote_state)) {
    return(FALSE)
  }
  if (!all(remote_state$ok %in% TRUE)) {
    return(FALSE)
  }

  previous <- manifest$remote_asset_state
  if (is.null(previous) || !is.data.frame(previous) || !nrow(previous)) {
    return(TRUE)
  }

  previous <- previous[, intersect(c("asset_name", "last_modified", "content_length", "etag"), names(previous)), drop = FALSE]
  current <- remote_state[, c("asset_name", "last_modified", "content_length", "etag"), drop = FALSE]
  previous <- previous[order(previous$asset_name), , drop = FALSE]
  current <- current[order(current$asset_name), , drop = FALSE]
  rownames(previous) <- NULL
  rownames(current) <- NULL

  !identical(previous, current)
}

.dnmb_dbcan_profile_family <- function(profile_id) {
  clean <- trimws(as.character(profile_id)[1])
  clean <- sub("\\.hmm$", "", clean, ignore.case = TRUE)
  clean <- sub("\\.txt$", "", clean, ignore.case = TRUE)
  if (grepl("^GT2_", clean)) {
    return("GT2")
  }
  clean
}

.dnmb_dbcan_empty_hits <- function() {
  tibble::tibble(
    query = character(),
    profile_id = character(),
    profile_length = integer(),
    gene_length = integer(),
    evalue = numeric(),
    profile_start = integer(),
    profile_end = integer(),
    gene_start = integer(),
    gene_end = integer(),
    coverage = numeric(),
    family_id = character(),
    hit_label = character(),
    substrate_label = character()
  )
}

.dnmb_dbcan_filter_overlaps <- function(hits) {
  if (is.null(hits) || !nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }

  split_hits <- split(as.data.frame(hits, stringsAsFactors = FALSE), hits$query)
  filtered <- lapply(split_hits, function(tbl) {
    tbl <- tbl[order(tbl$gene_start, tbl$gene_end, tbl$evalue), , drop = FALSE]
    if (nrow(tbl) <= 1L) {
      return(tbl)
    }

    i <- 1L
    while (i < nrow(tbl)) {
      len1 <- tbl$gene_end[[i]] - tbl$gene_start[[i]]
      len2 <- tbl$gene_end[[i + 1L]] - tbl$gene_start[[i + 1L]]
      overlap <- tbl$gene_end[[i]] - tbl$gene_start[[i + 1L]]
      if (!is.na(overlap) && overlap > 0 && (
        (!is.na(len1) && len1 > 0 && overlap / len1 > 0.5) ||
          (!is.na(len2) && len2 > 0 && overlap / len2 > 0.5)
      )) {
        drop_index <- if (!is.na(tbl$evalue[[i]]) && !is.na(tbl$evalue[[i + 1L]]) && tbl$evalue[[i]] < tbl$evalue[[i + 1L]]) {
          i + 1L
        } else {
          i
        }
        tbl <- tbl[-drop_index, , drop = FALSE]
        if (i > 1L) {
          i <- i - 1L
        }
      } else {
        i <- i + 1L
      }
    }
    tbl
  })

  tibble::as_tibble(dplyr::bind_rows(filtered))
}

dnmb_dbcan_parse_domtblout <- function(path,
                                       evalue_threshold = 1e-15,
                                       coverage_threshold = 0.35) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(.dnmb_dbcan_empty_hits())
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (!length(lines)) {
    return(.dnmb_dbcan_empty_hits())
  }

  parsed <- lapply(lines, function(line) {
    fields <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(fields) < 19L) {
      return(NULL)
    }

    profile_id <- fields[[4]]
    profile_length <- suppressWarnings(as.integer(fields[[6]]))
    gene_id <- fields[[1]]
    gene_length <- suppressWarnings(as.integer(fields[[3]]))
    evalue <- suppressWarnings(as.numeric(fields[[13]]))
    profile_start <- suppressWarnings(as.integer(fields[[16]]))
    profile_end <- suppressWarnings(as.integer(fields[[17]]))
    gene_start <- suppressWarnings(as.integer(fields[[18]]))
    gene_end <- suppressWarnings(as.integer(fields[[19]]))
    coverage <- if (is.na(profile_length) || profile_length <= 0L || is.na(profile_start) || is.na(profile_end)) {
      NA_real_
    } else {
      (profile_end - profile_start) / profile_length
    }
    family_id <- .dnmb_dbcan_profile_family(profile_id)

    tibble::tibble(
      query = .dnmb_module_clean_annotation_key(gene_id),
      profile_id = profile_id,
      profile_length = profile_length,
      gene_length = gene_length,
      evalue = evalue,
      profile_start = profile_start,
      profile_end = profile_end,
      gene_start = gene_start,
      gene_end = gene_end,
      coverage = coverage,
      family_id = family_id,
      hit_label = family_id,
      substrate_label = NA_character_
    )
  })

  hits <- dplyr::bind_rows(parsed)
  if (!nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }

  hits <- .dnmb_dbcan_filter_overlaps(hits)
  hits <- hits[!is.na(hits$evalue) & hits$evalue <= evalue_threshold & !is.na(hits$coverage) & hits$coverage >= coverage_threshold, , drop = FALSE]
  rownames(hits) <- NULL
  if (!nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }

  hits <- hits[, names(.dnmb_dbcan_empty_hits()), drop = FALSE]
  tibble::as_tibble(hits)
}

dnmb_dbcan_normalize_hits <- function(hits) {
  if (is.null(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  required <- c("query", "family_id", "hit_label", "evalue", "coverage", "profile_length", "gene_length", "profile_start", "profile_end", "gene_start", "gene_end")
  missing <- setdiff(required, names(hits))
  if (length(missing)) {
    stop("`hits` is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  support <- paste0(
    "evalue=", signif(hits$evalue, 4),
    "; coverage=", sprintf("%.3f", hits$coverage),
    "; profile_length=", hits$profile_length,
    "; gene_length=", hits$gene_length,
    "; profile_start=", hits$profile_start,
    "; profile_end=", hits$profile_end,
    "; gene_start=", hits$gene_start,
    "; gene_end=", hits$gene_end
  )

  out <- data.frame(
    query = .dnmb_module_clean_annotation_key(hits$query),
    source = rep("dbcan", nrow(hits)),
    family_system = rep("dbCAN", nrow(hits)),
    family_id = as.character(hits$family_id),
    hit_label = as.character(hits$hit_label),
    enzyme_role = rep("CAZyme", nrow(hits)),
    evidence_mode = rep("direct", nrow(hits)),
    substrate_label = as.character(hits$substrate_label),
    support = support,
    profile_id = as.character(hits$profile_id),
    evalue = as.numeric(hits$evalue),
    coverage = as.numeric(hits$coverage),
    profile_length = as.integer(hits$profile_length),
    gene_length = as.integer(hits$gene_length),
    profile_start = as.integer(hits$profile_start),
    profile_end = as.integer(hits$profile_end),
    gene_start = as.integer(hits$gene_start),
    gene_end = as.integer(hits$gene_end),
    typing_eligible = rep(TRUE, nrow(hits)),
    stringsAsFactors = FALSE
  )
  out <- out[, c(.dnmb_module_optional_long_columns(), setdiff(names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  hit_order <- order(out$query, hits$evalue, -hits$coverage)
  out <- out[hit_order, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_dbcan_output_table <- function(genes, hits, cgc_genes = NULL, substrate = NULL) {
  out <- .dnmb_module_output_table(genes = genes, hits = hits)
  out$dbcan_hit <- out$family_id
  if (!is.null(cgc_genes) && is.data.frame(cgc_genes) && nrow(cgc_genes)) {
    cgc_genes <- as.data.frame(cgc_genes, stringsAsFactors = FALSE)
    cgc_genes$query <- .dnmb_module_clean_annotation_key(cgc_genes$query)
    join_cols <- setdiff(names(cgc_genes), "query")
    match_rows <- match(out$locus_tag, cgc_genes$query)
    include <- !is.na(match_rows)
    for (column_name in join_cols) {
      if (!column_name %in% names(out)) {
        out[[column_name]] <- .dnmb_na_vector_like(cgc_genes[[column_name]], nrow(out))
      }
      out[[column_name]][include] <- cgc_genes[[column_name]][match_rows[include]]
    }
  }
  if (!is.null(substrate) && is.data.frame(substrate) && nrow(substrate) && "dbcan_cgc_id" %in% names(out)) {
    substrate <- as.data.frame(substrate, stringsAsFactors = FALSE)
    match_rows <- match(out$dbcan_cgc_id, substrate$dbcan_cgc_id)
    include <- !is.na(match_rows)
    join_cols <- setdiff(names(substrate), "dbcan_cgc_id")
    for (column_name in join_cols) {
      if (!column_name %in% names(out)) {
        out[[column_name]] <- .dnmb_na_vector_like(substrate[[column_name]], nrow(out))
      }
      out[[column_name]][include] <- substrate[[column_name]][match_rows[include]]
    }
  }
  drop_cols <- intersect(
    c(
      "hit_label",
      "enzyme_role",
      "evidence_mode",
      "substrate_label",
      "typing_eligible"
    ),
    names(out)
  )
  if (length(drop_cols)) {
    out[drop_cols] <- NULL
  }

  base_cols <- intersect(dnmb_backbone_columns(), names(out))
  dbcan_cols <- c("dbcan_hit", setdiff(names(out), c(base_cols, "dbcan_hit")))
  out[, c(base_cols, dbcan_cols), drop = FALSE]
}

.dnmb_dbcan_standalone_layout <- function(module_dir) {
  file.path(module_dir, "standalone_bundle")
}

.dnmb_dbcan_supporting_urls <- function(base_url = .dnmb_dbcan_default_base_url()) {
  c(
    "fam-substrate-mapping.tsv" = "https://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv",
    "PUL.faa" = "https://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa",
    "dbCAN-PUL.xlsx" = "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx",
    "dbCAN-PUL.tar.gz" = "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz",
    "dbCAN_sub.hmm" = "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm",
    "CAZyDB.fa" = "https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa",
    "tcdb.fa" = "https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa",
    "tf-1.hmm" = "https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm",
    "tf-2.hmm" = "https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm",
    "stp.hmm" = "https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm"
  )
}

.dnmb_dbcan_supporting_paths <- function(bundle_dir) {
  list(
    bundle_dir = bundle_dir,
    hmm_txt = file.path(bundle_dir, "dbCAN.txt"),
    hmm_alias = file.path(bundle_dir, "dbCAN-HMMdb-V14.txt"),
    fam_substrate_mapping = file.path(bundle_dir, "fam-substrate-mapping.tsv"),
    pul_faa = file.path(bundle_dir, "PUL.faa"),
    pul_excel = file.path(bundle_dir, "dbCAN-PUL.xlsx"),
    pul_tar = file.path(bundle_dir, "dbCAN-PUL.tar.gz"),
    pul_dir = file.path(bundle_dir, "dbCAN-PUL"),
    tcdb_fa = file.path(bundle_dir, "tcdb.fa"),
    tcdb_dmnd = file.path(bundle_dir, "tcdb.dmnd"),
    tf1_hmm = file.path(bundle_dir, "tf-1.hmm"),
    tf2_hmm = file.path(bundle_dir, "tf-2.hmm"),
    stp_hmm = file.path(bundle_dir, "stp.hmm"),
    dbcan_sub_hmm = file.path(bundle_dir, "dbCAN_sub.hmm"),
    cazy_fa = file.path(bundle_dir, "CAZyDB.fa"),
    cazy_dmnd = file.path(bundle_dir, "CAZy.dmnd")
  )
}

.dnmb_dbcan_standalone_ready <- function(paths) {
  required <- c(
    paths$hmm_txt,
    paste0(paths$hmm_txt, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$hmm_alias,
    paste0(paths$hmm_alias, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$fam_substrate_mapping,
    paths$pul_faa,
    paths$pul_excel,
    paths$pul_dir,
    paths$tcdb_dmnd,
    paths$tf1_hmm,
    paths$tf2_hmm,
    paths$stp_hmm,
    paste0(paths$tf1_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paste0(paths$tf2_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paste0(paths$stp_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$cazy_dmnd,
    paths$dbcan_sub_hmm,
    paste0(paths$dbcan_sub_hmm, c(".h3f", ".h3i", ".h3m", ".h3p"))
  )
  all(file.exists(required))
}

.dnmb_dbcan_stage_hmm_for_standalone <- function(module_dir, bundle_dir, version_label = NULL) {
  layout <- .dnmb_dbcan_asset_layout(module_dir)
  release_info <- .dnmb_dbcan_release_info()
  version_label <- version_label %||% release_info$version %||% "V14"
  txt_base <- file.path(bundle_dir, "dbCAN.txt")
  alias_base <- file.path(bundle_dir, paste0("dbCAN-HMMdb-", version_label, ".txt"))
  source_files <- c(layout$hmm_path, layout$hmmpress_files)
  target_sets <- list(
    c(txt_base, paste0(txt_base, c(".h3f", ".h3i", ".h3m", ".h3p"))),
    c(alias_base, paste0(alias_base, c(".h3f", ".h3i", ".h3m", ".h3p")))
  )
  target_files <- unlist(target_sets, use.names = FALSE)
  if (all(file.exists(target_files))) {
    return(list(ok = TRUE, alias = alias_base, files = target_files))
  }
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- logical()
  for (target_set in target_sets) {
    copied <- c(copied, file.copy(source_files, target_set, overwrite = TRUE))
  }
  list(ok = all(copied), alias = alias_base, files = target_files[file.exists(target_files)])
}

.dnmb_dbcan_prepare_reference_indexes <- function(paths) {
  if (!file.exists(paths$cazy_dmnd)) {
    dnmb_run_external("diamond", c("makedb", "--in", paths$cazy_fa, "-d", sub("\\.dmnd$", "", paths$cazy_dmnd)), required = TRUE)
  }
  if (!file.exists(paste0(paths$dbcan_sub_hmm, ".h3f"))) {
    dnmb_run_external("hmmpress", c("-f", paths$dbcan_sub_hmm), required = TRUE)
  }
}

.dnmb_dbcan_install_supporting_bundle <- function(module_dir,
                                                  trace_log,
                                                  base_url = .dnmb_dbcan_default_base_url(),
                                                  force = FALSE) {
  bundle_dir <- .dnmb_dbcan_standalone_layout(module_dir)
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  module_layout <- .dnmb_dbcan_asset_layout(module_dir)
  paths <- .dnmb_dbcan_supporting_paths(bundle_dir)
  urls <- .dnmb_dbcan_supporting_urls(base_url = base_url)

  status <- .dnmb_dbcan_empty_status()
  if (!isTRUE(force) && .dnmb_dbcan_standalone_ready(paths)) {
    return(list(ok = TRUE, status = .dnmb_dbcan_status_row("dbcan_support_bundle", "cached", bundle_dir), paths = paths))
  }

  for (asset_name in names(urls)) {
    dest <- switch(
      asset_name,
      "PUL.faa" = paths$pul_faa,
      "fam-substrate-mapping.tsv" = paths$fam_substrate_mapping,
      "dbCAN-PUL.xlsx" = paths$pul_excel,
      "dbCAN-PUL.tar.gz" = paths$pul_tar,
      "dbCAN_sub.hmm" = paths$dbcan_sub_hmm,
      "CAZyDB.fa" = paths$cazy_fa,
      "tcdb.fa" = paths$tcdb_fa,
      "tf-1.hmm" = paths$tf1_hmm,
      "tf-2.hmm" = paths$tf2_hmm,
      "stp.hmm" = paths$stp_hmm,
      file.path(bundle_dir, asset_name)
    )
    result <- .dnmb_download_asset(urls[[asset_name]], dest, insecure = TRUE)
    ok <- isTRUE(result$ok) && file.exists(dest)
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row(
        paste0("dbcan_support_download:", asset_name),
        if (ok) "ok" else "failed",
        if (ok) dest else (result$error %||% asset_name)
      )
    )
    if (!ok) {
      return(list(ok = FALSE, status = status, paths = paths))
    }
  }

  if (!file.exists(paths$pul_dir) && file.exists(paths$pul_tar)) {
    untar(paths$pul_tar, exdir = bundle_dir)
  }
  if (!file.exists(paste0(paths$pul_faa, ".phr"))) {
    pul_prepare <- .dnmb_dbcan_prepare_pul_blast_db(paths$pul_faa, trace_log = trace_log)
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row("dbcan_support_bundle:makeblastdb", pul_prepare$status, pul_prepare$detail)
    )
    if (!isTRUE(pul_prepare$ok)) {
      return(list(ok = FALSE, status = status, paths = paths))
    }
  }
  if (!file.exists(paths$tcdb_dmnd)) {
    dnmb_run_external("diamond", c("makedb", "--in", paths$tcdb_fa, "-d", file.path(bundle_dir, "tcdb")), required = TRUE)
  }
  if (!file.exists(paste0(paths$tf1_hmm, ".h3f"))) {
    dnmb_run_external("hmmpress", c("-f", paths$tf1_hmm), required = TRUE)
  }
  if (!file.exists(paste0(paths$tf2_hmm, ".h3f"))) {
    dnmb_run_external("hmmpress", c("-f", paths$tf2_hmm), required = TRUE)
  }
  if (!file.exists(paste0(paths$stp_hmm, ".h3f"))) {
    dnmb_run_external("hmmpress", c("-f", paths$stp_hmm), required = TRUE)
  }
  hmm_alias <- .dnmb_dbcan_stage_hmm_for_standalone(module_dir = module_dir, bundle_dir = bundle_dir)
  status <- dplyr::bind_rows(
    status,
    .dnmb_dbcan_status_row(
      "dbcan_support_bundle:hmm_alias",
      if (isTRUE(hmm_alias$ok)) "ok" else "failed",
      if (isTRUE(hmm_alias$ok)) hmm_alias$alias else module_layout$hmm_path
    )
  )
  if (!isTRUE(hmm_alias$ok)) {
    return(list(ok = FALSE, status = status, paths = paths))
  }
  .dnmb_dbcan_prepare_reference_indexes(paths)
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] dbcan_support_bundle=%s", Sys.time(), bundle_dir))

  list(
    ok = .dnmb_dbcan_standalone_ready(paths),
    status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_support_bundle", "ok", bundle_dir)),
    paths = paths
  )
}

.dnmb_dbcan_can_run_standalone <- function(genes) {
  # The DNMB "standalone" path was written against dbCAN 3.x CLI
  # (`run_dbcan <query> protein --tools all --cgc_substrate ...`).
  # Bioconda's current `dbcan` package is 5.2+, which replaced that
  # interface with subcommands like `run_dbcan easy_substrate
  # --input_raw_data <fasta> --mode protein --db_dir <db> --output_dir
  # <out>`, AND the legacy HTTP mirrors for the supporting bundle are
  # dead (everything 302-redirects to a landing page). Until the
  # module is ported to the new CLI + S3 bundle layout, force the
  # hmmsearch fallback path so the module at least produces HMM-level
  # hits from `dbCAN.hmm` (downloaded from the S3 mirror). Users can
  # re-enable the old path via `options(dnmb.dbcan.standalone = TRUE)`
  # once the rewrite lands.
  if (!isTRUE(getOption("dnmb.dbcan.standalone", FALSE))) {
    return(FALSE)
  }
  required <- c("locus_tag", "contig", "start", "end", "direction")
  all(required %in% names(genes))
}

.dnmb_dbcan_safe_contig_map <- function(genes) {
  contigs <- unique(as.character(genes$contig))
  contigs <- contigs[!is.na(contigs) & nzchar(contigs)]
  tibble::tibble(
    contig = contigs,
    safe_contig = paste0("ctg", seq_along(contigs))
  )
}

.dnmb_dbcan_write_cluster_tsv <- function(genes, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- as.data.frame(genes, stringsAsFactors = FALSE)
  out$locus_tag <- .dnmb_module_clean_annotation_key(out$locus_tag)
  out <- out[!is.na(out$locus_tag) & nzchar(out$locus_tag), c("contig", "locus_tag", "start", "end", "direction"), drop = FALSE]
  mapping <- .dnmb_dbcan_safe_contig_map(out)
  out$contig <- mapping$safe_contig[match(out$contig, mapping$contig)]
  utils::write.table(out, file = path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  list(path = path, contig_map = mapping)
}

.dnmb_dbcan_restore_cgc_ids <- function(tbl, contig_map) {
  if (is.null(tbl) || !is.data.frame(tbl) || !nrow(tbl) || is.null(contig_map) || !nrow(contig_map) || !"dbcan_cgc_id" %in% names(tbl)) {
    return(tbl)
  }
  out <- as.data.frame(tbl, stringsAsFactors = FALSE)
  out$dbcan_cgc_id <- as.character(out$dbcan_cgc_id)
  for (i in seq_len(nrow(contig_map))) {
    safe <- paste0("^", gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", contig_map$safe_contig[[i]]), "\\|")
    repl <- paste0(contig_map$contig[[i]], "|")
    out$dbcan_cgc_id <- sub(safe, repl, out$dbcan_cgc_id, perl = TRUE)
  }
  out
}

.dnmb_dbcan_stage_dir_for_synteny <- function(source_dir) {
  stage_parent <- tempfile("dnmb-dbcan-synteny-db-")
  dir.create(stage_parent, recursive = TRUE, showWarnings = FALSE)
  stage_dir <- file.path(stage_parent, "db")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  source_pul_dir <- normalizePath(file.path(source_dir, "dbCAN-PUL"), winslash = "/", mustWork = TRUE)
  target_pul_dir <- file.path(stage_dir, "dbCAN-PUL")
  linked <- suppressWarnings(file.symlink(source_pul_dir, target_pul_dir))
  if (!isTRUE(linked)) {
    dir.create(target_pul_dir, recursive = TRUE, showWarnings = FALSE)
    copy_ok <- file.copy(list.files(source_pul_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE), target_pul_dir, recursive = TRUE)
    if (!all(copy_ok)) {
      stop("Failed to stage dbCAN-PUL database directory for synteny plotting.", call. = FALSE)
    }
  }
  list(
    db_dir = stage_dir,
    cleanup = function() unlink(stage_parent, recursive = TRUE, force = TRUE)
  )
}

.dnmb_dbcan_copy_dir_contents <- function(src_dir, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(src_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- file.copy(entries, dest_dir, recursive = TRUE, overwrite = TRUE)
  if (!all(ok)) {
    stop("Failed to copy staged run_dbcan outputs.", call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_dbcan_stage_for_run <- function(query_fasta, cluster_path, bundle_dir) {
  stage_parent <- tempfile("dnmb-dbcan-run-")
  dir.create(stage_parent, recursive = TRUE, showWarnings = FALSE)
  stage_query <- file.path(stage_parent, "query.faa")
  stage_cluster <- file.path(stage_parent, "cluster.tsv")
  stage_db <- file.path(stage_parent, "db")
  stage_out <- file.path(stage_parent, "run_dbcan")
  dir.create(stage_out, recursive = TRUE, showWarnings = FALSE)
  file.copy(query_fasta, stage_query, overwrite = TRUE)
  file.copy(cluster_path, stage_cluster, overwrite = TRUE)
  linked <- suppressWarnings(file.symlink(normalizePath(bundle_dir, winslash = "/", mustWork = TRUE), stage_db))
  if (!isTRUE(linked)) {
    .dnmb_dbcan_copy_dir_contents(bundle_dir, stage_db)
  }
  list(
    query_fasta = stage_query,
    cluster_tsv = stage_cluster,
    db_dir = stage_db,
    out_dir = stage_out,
    cleanup = function() unlink(stage_parent, recursive = TRUE, force = TRUE)
  )
}

.dnmb_dbcan_run_syntenic_plot <- function(run_dbcan_dir, bundle_dir, trace_log = NULL) {
  required_files <- c(
    file.path(run_dbcan_dir, "PUL_blast.out"),
    file.path(run_dbcan_dir, "cgc_standard.out"),
    file.path(run_dbcan_dir, "substrate.out")
  )
  if (!all(file.exists(required_files))) {
    return(list(
      ok = FALSE,
      status = .dnmb_dbcan_status_row("dbcan_synteny", "missing_inputs", run_dbcan_dir),
      command = NULL,
      files = list()
    ))
  }

  stage <- .dnmb_dbcan_stage_dir_for_synteny(bundle_dir)
  on.exit(stage$cleanup(), add = TRUE)
  args <- c(
    "syntenic_plot",
    "-b", "PUL_blast.out",
    "--cgc", "cgc_standard.out",
    "-i", "substrate.out",
    "--db", stage$db_dir
  )
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] syntenic_plot args=%s", Sys.time(), .dnmb_format_command("syntenic_plot", args)))
  }
  command <- dnmb_run_external("syntenic_plot", args = args, wd = run_dbcan_dir, required = FALSE)
  synteny_dir <- file.path(run_dbcan_dir, "synteny.pdf")
  ok <- isTRUE(command$ok) && dir.exists(synteny_dir)
  status <- .dnmb_dbcan_status_row(
    "dbcan_synteny",
    if (ok) "ok" else if (!nzchar(command$resolved_command)) "missing" else "failed",
    if (ok) synteny_dir else (command$error %||% synteny_dir)
  )

  list(
    ok = ok,
    status = status,
    command = command,
    files = if (dir.exists(synteny_dir)) list(synteny_dir = synteny_dir) else list()
  )
}

dnmb_dbcan_parse_hmmer_table <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(.dnmb_dbcan_empty_hits())
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) <= 1L) {
    return(.dnmb_dbcan_empty_hits())
  }
  lines <- lines[-1L]
  parsed <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 10L) {
      return(NULL)
    }
    tibble::tibble(
      query = .dnmb_module_clean_annotation_key(fields[[3]]),
      profile_id = fields[[1]],
      profile_length = suppressWarnings(as.integer(fields[[2]])),
      gene_length = suppressWarnings(as.integer(fields[[4]])),
      evalue = suppressWarnings(as.numeric(fields[[5]])),
      profile_start = suppressWarnings(as.integer(fields[[6]])),
      profile_end = suppressWarnings(as.integer(fields[[7]])),
      gene_start = suppressWarnings(as.integer(fields[[8]])),
      gene_end = suppressWarnings(as.integer(fields[[9]])),
      coverage = suppressWarnings(as.numeric(fields[[10]])),
      family_id = .dnmb_dbcan_profile_family(fields[[1]]),
      hit_label = .dnmb_dbcan_profile_family(fields[[1]]),
      substrate_label = NA_character_
    )
  })
  hits <- dplyr::bind_rows(parsed)
  if (!nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }
  hits <- hits[, names(.dnmb_dbcan_empty_hits()), drop = FALSE]
  tibble::as_tibble(hits)
}

dnmb_dbcan_parse_cgc_standard <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(tibble::tibble(
      query = character(),
      dbcan_cgc_id = character(),
      dbcan_cgc_gene_type = character(),
      dbcan_cgc_protein_family = character()
    ))
  }

  lines <- readLines(path, warn = FALSE)
  if (length(lines) <= 1L) {
    return(tibble::tibble(
      query = character(),
      dbcan_cgc_id = character(),
      dbcan_cgc_gene_type = character(),
      dbcan_cgc_protein_family = character()
    ))
  }

  parsed <- lapply(lines[-1L], function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 8L) {
      return(NULL)
    }
    tibble::tibble(
      query = .dnmb_module_clean_annotation_key(fields[[4]]),
      dbcan_cgc_id = paste(fields[[3]], fields[[1]], sep = "|"),
      dbcan_cgc_gene_type = fields[[2]],
      dbcan_cgc_protein_family = fields[[8]]
    )
  })

  dplyr::bind_rows(parsed)
}

dnmb_dbcan_parse_substrate_table <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(tibble::tibble(
      dbcan_cgc_id = character(),
      dbcan_pul_id = character(),
      dbcan_pul_substrate = character(),
      dbcan_pul_bitscore = numeric(),
      dbcan_pul_signature_pairs = character(),
      dbcan_sub_substrate = character(),
      dbcan_sub_score = character()
    ))
  }

  tbl <- utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    check.names = FALSE
  )
  if (!nrow(tbl)) {
    return(tibble::tibble(
      dbcan_cgc_id = character(),
      dbcan_pul_id = character(),
      dbcan_pul_substrate = character(),
      dbcan_pul_bitscore = numeric(),
      dbcan_pul_signature_pairs = character(),
      dbcan_sub_substrate = character(),
      dbcan_sub_score = character()
    ))
  }
  names(tbl)[seq_len(min(7L, ncol(tbl)))] <- c(
    "dbcan_cgc_id",
    "dbcan_pul_id",
    "dbcan_pul_substrate",
    "dbcan_pul_bitscore",
    "dbcan_pul_signature_pairs",
    "dbcan_sub_substrate",
    "dbcan_sub_score"
  )[seq_len(min(7L, ncol(tbl)))]
  out <- tibble::as_tibble(tbl)
  out$dbcan_pul_bitscore <- suppressWarnings(as.numeric(out$dbcan_pul_bitscore))
  out[, c("dbcan_cgc_id", "dbcan_pul_id", "dbcan_pul_substrate", "dbcan_pul_bitscore", "dbcan_pul_signature_pairs", "dbcan_sub_substrate", "dbcan_sub_score")]
}

.dnmb_dbcan_augment_hits_with_standalone <- function(hits, cgc_genes = NULL, substrate = NULL) {
  hits <- as.data.frame(hits, stringsAsFactors = FALSE)
  if (!is.null(cgc_genes) && is.data.frame(cgc_genes) && nrow(cgc_genes)) {
    hits <- dplyr::left_join(hits, cgc_genes, by = c("query" = "query"))
  }
  if (!is.null(substrate) && is.data.frame(substrate) && nrow(substrate) && "dbcan_cgc_id" %in% names(hits)) {
    hits <- dplyr::left_join(hits, substrate, by = "dbcan_cgc_id")
  }
  hits
}

dnmb_dbcan_run_standalone <- function(query_fasta,
                                      genes,
                                      output_dir,
                                      version = .dnmb_dbcan_default_version(),
                                      cache_root = NULL,
                                      cpu = 1L,
                                      install = FALSE,
                                      base_url = .dnmb_dbcan_default_base_url()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_dbcan_empty_status()
  trace_log <- file.path(output_dir, "dbcan_trace.log")
  module_dir <- .dnmb_db_module_dir(.dnmb_dbcan_module_name(), version, cache_root = cache_root, create = TRUE)
  bundle_dir <- .dnmb_dbcan_standalone_layout(module_dir)

  module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  if (!isTRUE(module$ok) && isTRUE(install)) {
    install_result <- dnmb_dbcan_install_module(version = version, cache_root = cache_root, base_url = base_url, force = FALSE, prepare = TRUE)
    status <- dplyr::bind_rows(status, install_result$status)
    module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  }
  if (!isTRUE(module$ok)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_module", "missing", paste0("dbCAN module not installed for version ", version))),
      files = list(query_fasta = query_fasta),
      command = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      cgc_genes = tibble::tibble(),
      substrate = tibble::tibble(),
      manifest = module$manifest
    ))
  }

  support_bundle <- .dnmb_dbcan_install_supporting_bundle(
    module_dir = module_dir,
    trace_log = trace_log,
    base_url = base_url,
    force = FALSE
  )
  status <- dplyr::bind_rows(status, support_bundle$status)
  if (!isTRUE(support_bundle$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = list(query_fasta = query_fasta, trace_log = trace_log),
      command = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      cgc_genes = tibble::tibble(),
      substrate = tibble::tibble(),
      manifest = module$manifest
    ))
  }

  cluster_path <- file.path(output_dir, "dbcan_cluster.tsv")
  cluster_info <- .dnmb_dbcan_write_cluster_tsv(genes, cluster_path)
  run_dbcan_output <- file.path(output_dir, "run_dbcan")
  dir.create(run_dbcan_output, recursive = TRUE, showWarnings = FALSE)
  stage <- .dnmb_dbcan_stage_for_run(query_fasta = query_fasta, cluster_path = cluster_info$path, bundle_dir = bundle_dir)
  on.exit(stage$cleanup(), add = TRUE)
  args <- c(
    stage$query_fasta,
    "protein",
    "--cluster", stage$cluster_tsv,
    "--tools", "all",
    "--cgc_substrate",
    "--cgc_sig_genes", "all",
    "--db_dir", stage$db_dir,
    "--out_dir", stage$out_dir
  )
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] run_dbcan args=%s", Sys.time(), .dnmb_format_command("run_dbcan", args)))
  command <- dnmb_run_external("run_dbcan", args = args, required = FALSE)
  if (dir.exists(stage$out_dir)) {
    .dnmb_dbcan_copy_dir_contents(stage$out_dir, run_dbcan_output)
  }

  hmmer_path <- file.path(run_dbcan_output, "hmmer.out")
  overview_path <- file.path(run_dbcan_output, "overview.txt")
  cgc_standard_path <- file.path(run_dbcan_output, "cgc_standard.out")
  substrate_path <- file.path(run_dbcan_output, "substrate.out")

  raw_hits <- dnmb_dbcan_parse_hmmer_table(hmmer_path)
  cgc_genes <- .dnmb_dbcan_restore_cgc_ids(
    dnmb_dbcan_parse_cgc_standard(cgc_standard_path),
    contig_map = cluster_info$contig_map
  )
  substrate <- .dnmb_dbcan_restore_cgc_ids(
    dnmb_dbcan_parse_substrate_table(substrate_path),
    contig_map = cluster_info$contig_map
  )
  hits <- dnmb_dbcan_normalize_hits(.dnmb_dbcan_augment_hits_with_standalone(raw_hits, cgc_genes = cgc_genes, substrate = substrate))
  synteny <- .dnmb_dbcan_run_syntenic_plot(
    run_dbcan_dir = run_dbcan_output,
    bundle_dir = bundle_dir,
    trace_log = trace_log
  )

  ok <- file.exists(overview_path) && file.exists(cgc_standard_path)
  status <- dplyr::bind_rows(
    status,
    .dnmb_dbcan_status_row("run_dbcan", if (ok) "ok" else if (isTRUE(command$ok)) "partial" else "failed", run_dbcan_output),
    .dnmb_dbcan_status_row("dbcan_overview", if (file.exists(overview_path)) "ok" else "missing", overview_path),
    .dnmb_dbcan_status_row("dbcan_cgc", if (file.exists(cgc_standard_path)) "ok" else "missing", cgc_standard_path),
    .dnmb_dbcan_status_row("dbcan_substrate", if (file.exists(substrate_path)) "ok" else "missing", substrate_path),
    synteny$status
  )

  list(
    ok = ok || (isTRUE(command$ok) && file.exists(overview_path)),
    status = status,
    files = c(
      list(
        query_fasta = query_fasta,
        cluster_tsv = cluster_path,
        run_dbcan_dir = run_dbcan_output,
        hmmer_out = hmmer_path,
        overview = overview_path,
        cgc_standard = cgc_standard_path,
        substrate = substrate_path,
        trace_log = trace_log
      ),
      synteny$files,
      module$files
    ),
    command = command,
    raw_hits = raw_hits,
    hits = hits,
    cgc_genes = cgc_genes,
    substrate = substrate,
    manifest = module$manifest
  )
}

dnmb_dbcan_install_module <- function(version = .dnmb_dbcan_default_version(),
                                      cache_root = NULL,
                                      base_url = .dnmb_dbcan_default_base_url(),
                                      asset_urls = NULL,
                                      force = FALSE,
                                      prepare = TRUE) {
  module <- .dnmb_dbcan_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_dbcan_asset_layout(module_dir)
  install_trace_log <- file.path(module_dir, "dbcan_install_trace.log")
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  release_info <- .dnmb_dbcan_release_info()
  urls <- .dnmb_dbcan_asset_urls(base_url = base_url, version = version, asset_urls = asset_urls, release_info = release_info)
  status <- .dnmb_dbcan_empty_status()
  remote_asset_state <- .dnmb_dbcan_remote_asset_state(
    urls = urls,
    enabled = identical(trimws(as.character(version)[1]), "current") && is.null(asset_urls)
  )
  remote_refresh_needed <- .dnmb_dbcan_remote_update_needed(manifest, remote_asset_state)
  if (!is.null(remote_asset_state) && nrow(remote_asset_state)) {
    remote_detail <- paste(
      sprintf(
        "%s(last_modified=%s; content_length=%s)",
        remote_asset_state$asset_name,
        ifelse(is.na(remote_asset_state$last_modified), "NA", remote_asset_state$last_modified),
        ifelse(is.na(remote_asset_state$content_length), "NA", as.character(remote_asset_state$content_length))
      ),
      collapse = "; "
    )
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row(
        "dbcan_remote_state",
        if (all(remote_asset_state$ok %in% TRUE)) {
          if (remote_refresh_needed) "update_detected" else "ok"
        } else {
          "unavailable"
        },
        remote_detail
      )
    )
    .dnmb_dbcan_trace(install_trace_log, sprintf("[%s] remote_asset_state %s", Sys.time(), remote_detail))
  }
  if (isTRUE(remote_refresh_needed) && !isTRUE(force)) {
    force <- TRUE
    .dnmb_dbcan_trace(install_trace_log, sprintf("[%s] forcing refresh because current dbCAN metadata changed", Sys.time()))
  }

  required_assets <- c("dbCAN.txt")
  required_ready <- file.exists(layout$hmm_path) && .dnmb_dbcan_is_hmm(layout$hmm_path)
  manifest_ready <- !is.null(manifest)
  if (manifest_ready && required_ready && !isTRUE(force)) {
    prepare_result <- if (isTRUE(prepare)) {
      .dnmb_dbcan_prepare_hmm_db(layout$hmm_path, force = FALSE, trace_log = install_trace_log)
    } else {
      NULL
    }
    if (!is.null(prepare_result)) {
      status <- dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_prepare", prepare_result$status, prepare_result$detail))
    }
    cached_ok <- if (isTRUE(prepare) && !is.null(prepare_result)) {
      isTRUE(prepare_result$ok)
    } else {
      isTRUE(!is.null(manifest$install_ok) && manifest$install_ok) || .dnmb_dbcan_is_hmm(layout$hmm_path)
    }
    manifest_fields <- unclass(manifest)
    manifest_fields[c("module", "version", "module_dir", "manifest_path", "written_at")] <- NULL
    manifest_fields <- utils::modifyList(
      manifest_fields,
      c(
        .dnmb_dbcan_manifest_diagnostics(layout, prepare_result = prepare_result, release_info = release_info),
        list(
          hmmpress_files = layout$hmmpress_files[file.exists(layout$hmmpress_files)],
          prepared_with_hmmpress = isTRUE(!is.null(prepare_result) && prepare_result$status %in% c("ok", "cached") && length(prepare_result$hmmpress_files) == 4L),
          remote_asset_state = remote_asset_state,
          install_ok = cached_ok
        )
      )
    )
    dnmb_db_write_manifest(
      module = module,
      version = version,
      cache_root = cache_root,
      manifest = manifest_fields,
      overwrite = TRUE
    )
    .dnmb_db_autoprune_default_versions(
      module = module,
      version = version,
      default_version = .dnmb_dbcan_default_version(),
      cache_root = cache_root
    )
    manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
    return(list(
      ok = cached_ok,
      status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_install", if (cached_ok) "cached" else "failed", module_dir)),
      manifest = manifest,
      files = c(
        list(hmm_path = layout$hmm_path, trace_log = install_trace_log),
        if (file.exists(layout$mapping_path) && .dnmb_dbcan_is_mapping(layout$mapping_path)) list(mapping_path = layout$mapping_path),
        stats::setNames(as.list(layout$hmmpress_files[file.exists(layout$hmmpress_files)]), basename(layout$hmmpress_files[file.exists(layout$hmmpress_files)]))
      )
    ))
  }

  download_results <- list()
  for (asset_name in names(urls)) {
    asset_path <- switch(
      asset_name,
      "dbCAN.txt" = layout$hmm_path,
      "fam-substrate-mapping.tsv" = layout$mapping_path,
      file.path(module_dir, asset_name)
    )
    result <- .dnmb_download_asset(urls[[asset_name]], asset_path, insecure = TRUE)
    download_results[[asset_name]] <- result
    required_asset <- asset_name %in% required_assets
    asset_ok <- if (identical(asset_name, "dbCAN.txt")) {
      .dnmb_dbcan_is_hmm(asset_path)
    } else if (identical(asset_name, "fam-substrate-mapping.tsv")) {
      .dnmb_dbcan_is_mapping(asset_path)
    } else {
      isTRUE(result$ok) || file.exists(asset_path)
    }
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row(
        paste0("dbcan_download:", asset_name),
        if (asset_ok) "ok" else if (required_asset) "failed" else "optional_missing",
        if (asset_ok) asset_path else (result$error %||% asset_name)
      )
    )
  }

  ok <- .dnmb_dbcan_is_hmm(layout$hmm_path)
  prepare_result <- NULL
  if (ok && isTRUE(prepare)) {
    prepare_result <- .dnmb_dbcan_prepare_hmm_db(layout$hmm_path, force = isTRUE(force), trace_log = install_trace_log)
    status <- dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_prepare", prepare_result$status, prepare_result$detail))
    ok <- isTRUE(prepare_result$ok)
  }

  mapping_available <- file.exists(layout$mapping_path) && .dnmb_dbcan_is_mapping(layout$mapping_path)
  manifest_payload <- c(
    list(
      source = "dbcan",
      base_url = base_url,
      asset_urls = urls,
      asset_files = names(urls),
      required_assets = required_assets,
      hmm_path = layout$hmm_path,
      mapping_path = if (mapping_available) layout$mapping_path else NA_character_,
      hmmpress_files = layout$hmmpress_files[file.exists(layout$hmmpress_files)],
      prepared_with_hmmpress = isTRUE(!is.null(prepare_result) && prepare_result$status %in% c("ok", "cached") && length(prepare_result$hmmpress_files) == 4L),
      mapping_available = mapping_available,
      remote_asset_state = remote_asset_state
    ),
    .dnmb_dbcan_manifest_diagnostics(layout, prepare_result = prepare_result, release_info = release_info),
    list(
      download_methods = vapply(download_results, function(x) x$method %||% NA_character_, character(1)),
      install_ok = ok
    )
  )
  manifest_path <- dnmb_db_write_manifest(
    module = module,
    version = version,
    cache_root = cache_root,
    manifest = manifest_payload,
    overwrite = TRUE
  )
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_dbcan_default_version(),
    cache_root = cache_root
  )
  manifest <- readRDS(manifest_path)
  class(manifest) <- c("dnmb_db_manifest", class(manifest))

  list(
    ok = ok,
    status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_install", if (ok) "ok" else "failed", module_dir)),
    manifest = manifest,
    files = c(
      list(hmm_path = layout$hmm_path, trace_log = install_trace_log),
      if (mapping_available) list(mapping_path = layout$mapping_path),
      stats::setNames(as.list(layout$hmmpress_files[file.exists(layout$hmmpress_files)]), basename(layout$hmmpress_files[file.exists(layout$hmmpress_files)]))
    )
  )
}

dnmb_dbcan_get_module <- function(version = .dnmb_dbcan_default_version(),
                                  cache_root = NULL,
                                  required = FALSE) {
  module <- .dnmb_dbcan_module_name()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = FALSE)
  layout <- .dnmb_dbcan_asset_layout(module_dir)
  ok <- !is.null(manifest) && .dnmb_dbcan_is_hmm(layout$hmm_path)

  if (!ok && isTRUE(required)) {
    stop("dbCAN module is not installed for version `", version, "`.", call. = FALSE)
  }

  list(
    ok = ok,
    module = module,
    version = version,
    module_dir = module_dir,
    manifest = manifest,
    files = c(
      list(hmm_path = layout$hmm_path),
      if (file.exists(layout$mapping_path) && .dnmb_dbcan_is_mapping(layout$mapping_path)) list(mapping_path = layout$mapping_path),
      stats::setNames(as.list(layout$hmmpress_files[file.exists(layout$hmmpress_files)]), basename(layout$hmmpress_files[file.exists(layout$hmmpress_files)]))
    )
  )
}

dnmb_dbcan_run_hmmsearch <- function(query_fasta,
                                     output_dir,
                                     version = .dnmb_dbcan_default_version(),
                                     cache_root = NULL,
                                     cpu = 1L,
                                     install = FALSE,
                                     base_url = .dnmb_dbcan_default_base_url(),
                                     asset_urls = NULL,
                                     evalue_threshold = 1e-15,
                                     coverage_threshold = 0.35) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_dbcan_empty_status()
  trace_log <- file.path(output_dir, "dbcan_trace.log")
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] enter dnmb_dbcan_run_hmmsearch query=%s", Sys.time(), query_fasta))

  module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  if (!isTRUE(module$ok) && isTRUE(install)) {
    install_result <- dnmb_dbcan_install_module(
      version = version,
      cache_root = cache_root,
      base_url = base_url,
      asset_urls = asset_urls,
      force = FALSE,
      prepare = TRUE
    )
    status <- dplyr::bind_rows(status, install_result$status)
    module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  }

  if (!isTRUE(module$ok)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_module", "missing", paste0("dbCAN module not installed for version ", version))),
      files = list(query_fasta = query_fasta),
      command = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      manifest = module$manifest
    ))
  }

  prepare_result <- .dnmb_dbcan_prepare_hmm_db(module$files$hmm_path, force = FALSE, trace_log = trace_log)
  status <- dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_prepare", prepare_result$status, prepare_result$detail))
  if (!isTRUE(prepare_result$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = c(list(query_fasta = query_fasta), module$files),
      command = prepare_result$command,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      manifest = module$manifest
    ))
  }

  domtblout <- file.path(output_dir, "dbcan_hmmsearch.domtblout")
  stdout_log <- file.path(output_dir, "dbcan_hmmsearch.stdout.log")
  stderr_log <- file.path(output_dir, "dbcan_hmmsearch.stderr.log")
  hmm_stdout <- file.path(output_dir, "dbcan_hmmsearch.out")
  args <- c(
    "--domtblout", domtblout,
    "--cpu", as.character(as.integer(cpu)),
    "-o", hmm_stdout,
    module$files$hmm_path,
    query_fasta
  )
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] hmmsearch args=%s", Sys.time(), .dnmb_format_command("hmmsearch", args)))
  command <- dnmb_run_external("hmmsearch", args = args, required = FALSE)
  writeLines(command$stdout, con = stdout_log)
  writeLines(command$stderr, con = stderr_log)

  raw_hits <- dnmb_dbcan_parse_domtblout(
    domtblout,
    evalue_threshold = evalue_threshold,
    coverage_threshold = coverage_threshold
  )
  hits <- dnmb_dbcan_normalize_hits(raw_hits)
  status <- dplyr::bind_rows(
    status,
    .dnmb_dbcan_status_row(
      "hmmsearch",
      if (isTRUE(command$ok)) if (nrow(raw_hits)) "ok" else "empty" else if (!nzchar(command$resolved_command)) "missing" else "failed",
      if (isTRUE(command$ok)) domtblout else (command$error %||% "hmmsearch failed")
    )
  )

  list(
    ok = isTRUE(command$ok),
    status = status,
    files = c(
      list(
        query_fasta = query_fasta,
        domtblout = domtblout,
        hmm_stdout = hmm_stdout,
        stdout_log = stdout_log,
        stderr_log = stderr_log,
        trace_log = trace_log
      ),
      module$files
    ),
    command = command,
    raw_hits = raw_hits,
    hits = hits,
    manifest = module$manifest
  )
}

dnmb_run_dbcan_module <- function(genes,
                                  output_dir,
                                  version = .dnmb_dbcan_default_version(),
                                  cache_root = NULL,
                                  install = TRUE,
                                  base_url = .dnmb_dbcan_default_base_url(),
                                  asset_urls = NULL,
                                  cpu = 1L,
                                  genbank = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- file.path(output_dir, "dbcan_module_trace.log")
  fasta_path <- file.path(output_dir, "dbcan_query_proteins.faa")
  existing_faa <- dnmb_resolve_query_faa(genbank = genbank, output_dir = output_dir, fallback_filename = basename(fasta_path))
  if (!is.null(existing_faa) && .dnmb_can_reuse_query_fasta(existing_faa, genes)) {
    proteins <- .dnmb_prepare_query_proteins(genes)
    fasta <- list(path = existing_faa, n = nrow(proteins), proteins = proteins)
    fasta_path <- existing_faa
  } else {
    fasta <- .dnmb_write_query_fasta(genes, fasta_path)
  }
  status <- .dnmb_dbcan_status_row(
    "dbcan_query_fasta",
    if (fasta$n) "ok" else "empty",
    paste0("proteins=", fasta$n)
  )
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] query_faa=%s existing=%s proteins=%s", Sys.time(), fasta_path, !is.null(existing_faa), fasta$n))

  if (!fasta$n) {
    return(list(
      ok = TRUE,
      status = status,
      files = list(query_fasta = fasta_path, trace_log = trace_log),
      manifest = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      query_proteins = fasta$proteins
    ))
  }

  standalone_possible <- .dnmb_dbcan_can_run_standalone(genes) &&
    isTRUE(dnmb_detect_binary("run_dbcan", required = FALSE)$found)
  if (standalone_possible) {
    search <- dnmb_dbcan_run_standalone(
      query_fasta = fasta_path,
      genes = genes,
      output_dir = output_dir,
      version = version,
      cache_root = cache_root,
      cpu = cpu,
      install = install,
      base_url = base_url
    )
  } else {
    search <- dnmb_dbcan_run_hmmsearch(
      query_fasta = fasta_path,
      output_dir = output_dir,
      version = version,
      cache_root = cache_root,
      cpu = cpu,
      install = install,
      base_url = base_url,
      asset_urls = asset_urls
    )
  }

  list(
    ok = isTRUE(search$ok),
    status = dplyr::bind_rows(status, search$status),
    files = c(search$files, list(module_trace_log = trace_log)),
    manifest = search$manifest,
    raw_hits = search$raw_hits,
    hits = search$hits,
    cgc_genes = search$cgc_genes %||% tibble::tibble(),
    substrate = search$substrate %||% tibble::tibble(),
    output_table = .dnmb_dbcan_output_table(
      genes = genes,
      hits = search$hits,
      cgc_genes = search$cgc_genes %||% NULL,
      substrate = search$substrate %||% NULL
    ),
    query_proteins = fasta$proteins
  )
}
