.dnmb_merops_module_name <- function() {
  "merops"
}

.dnmb_merops_default_version <- function() {
  "current"
}

.dnmb_merops_default_base_url <- function() {
  "https://ftp.ebi.ac.uk/pub/databases/merops/current_release"
}

.dnmb_merops_status_row <- function(component, status, detail) {
  tibble::tibble(
    component = as.character(component),
    status = as.character(status),
    detail = as.character(detail)
  )
}

.dnmb_merops_empty_status <- function() {
  .dnmb_merops_status_row(character(), character(), character())
}

.dnmb_merops_trace <- function(path, text) {
  cat(text, "\n", file = path, append = TRUE)
}

.dnmb_merops_asset_layout <- function(module_dir) {
  fasta_path <- file.path(module_dir, "pepunit.lib")
  mapping_path <- file.path(module_dir, "dnld_list.txt")
  blast_db_prefix <- file.path(module_dir, "merops_reference")
  list(
    module_dir = module_dir,
    fasta_path = fasta_path,
    mapping_path = mapping_path,
    blast_db_prefix = blast_db_prefix,
    blast_db_files = .dnmb_merops_find_blast_db_files(blast_db_prefix)
  )
}

.dnmb_merops_find_blast_db_files <- function(blast_db_prefix) {
  all_exts <- c(".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto")
  candidates <- paste0(blast_db_prefix, all_exts)
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    return(paste0(blast_db_prefix, c(".phr", ".pin", ".psq")))
  }
  found
}

.dnmb_merops_trim_url_slash <- function(x) {
  sub("/+$", "", as.character(x)[1])
}

.dnmb_merops_asset_urls <- function(base_url = .dnmb_merops_default_base_url(),
                                     version = .dnmb_merops_default_version(),
                                     asset_urls = NULL) {
  if (!is.null(asset_urls)) {
    urls <- unlist(asset_urls, use.names = TRUE)
    urls <- stats::setNames(as.character(urls), names(urls))
    if (!length(urls) || is.null(names(urls)) || any(!nzchar(names(urls)))) {
      stop("`asset_urls` must be a named character vector or list.", call. = FALSE)
    }
    return(urls)
  }

  base_url <- .dnmb_merops_trim_url_slash(base_url)
  version <- trimws(as.character(version)[1])
  base_prefix <- if (is.na(version) || !nzchar(version) || identical(version, "current")) {
    base_url
  } else {
    paste(base_url, version, sep = "/")
  }

  stats::setNames(
    paste(base_prefix, c("pepunit.lib", "dnld_list.txt"), sep = "/"),
    c("pepunit.lib", "dnld_list.txt")
  )
}

.dnmb_merops_is_current_version <- function(version) {
  version <- trimws(as.character(version)[1])
  is.na(version) || !nzchar(version) || identical(version, "current")
}

.dnmb_merops_remote_asset_state <- function(urls, enabled = TRUE) {
  if (!isTRUE(enabled) || !length(urls)) {
    return(NULL)
  }

  rows <- lapply(names(urls), function(asset_name) {
    metadata <- .dnmb_remote_asset_metadata(urls[[asset_name]])
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

.dnmb_merops_remote_update_needed <- function(manifest, remote_state) {
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

.dnmb_merops_is_fasta <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(FALSE)
  }
  first_lines <- tryCatch(readLines(path, n = 3L, warn = FALSE), error = function(error) character())
  any(grepl("^>", first_lines))
}

.dnmb_merops_has_space <- function(path) {
  grepl("[[:space:]]", as.character(path)[1])
}

.dnmb_merops_stage_file <- function(source, staging_dir, target_name = basename(source)) {
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(staging_dir, target_name)
  if (file.exists(target)) {
    unlink(target, force = TRUE)
  }
  linked <- suppressWarnings(file.symlink(source, target))
  if (!isTRUE(linked)) {
    copied <- file.copy(source, target, overwrite = TRUE)
    if (!isTRUE(copied)) {
      stop("Failed to stage file for BLAST execution: ", source, call. = FALSE)
    }
  }
  file.path(
    normalizePath(dirname(target), winslash = "/", mustWork = TRUE),
    basename(target)
  )
}

.dnmb_merops_prefix_files <- function(db_prefix) {
  db_prefix <- normalizePath(db_prefix, winslash = "/", mustWork = FALSE)
  prefix_dir <- dirname(db_prefix)
  prefix_name <- basename(db_prefix)
  if (!dir.exists(prefix_dir)) {
    return(character())
  }
  list.files(
    prefix_dir,
    pattern = paste0("^", prefix_name, "\\."),
    full.names = TRUE
  )
}

.dnmb_merops_prepare_makeblastdb_paths <- function(fasta_path, db_prefix) {
  paths_with_spaces <- c(fasta_path, db_prefix)
  if (!any(vapply(paths_with_spaces, .dnmb_merops_has_space, logical(1)))) {
    return(list(
      staged = FALSE,
      cleanup = function() NULL,
      fasta_path = fasta_path,
      db_prefix = db_prefix,
      blast_db_files = paste0(db_prefix, c(".phr", ".pin", ".psq"))
    ))
  }

  staging_dir <- tempfile("dnmb-merops-makeblastdb-")
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  staged_fasta <- .dnmb_merops_stage_file(fasta_path, staging_dir, target_name = basename(fasta_path))
  staged_prefix <- file.path(staging_dir, basename(db_prefix))
  list(
    staged = TRUE,
    cleanup = function() unlink(staging_dir, recursive = TRUE, force = TRUE),
    fasta_path = staged_fasta,
    db_prefix = staged_prefix,
    blast_db_files = paste0(staged_prefix, c(".phr", ".pin", ".psq"))
  )
}

.dnmb_merops_finalize_makeblastdb_paths <- function(run_paths, final_db_prefix) {
  if (!isTRUE(run_paths$staged)) {
    return(paste0(final_db_prefix, c(".phr", ".pin", ".psq")))
  }

  staged_db_files <- .dnmb_merops_prefix_files(run_paths$db_prefix)
  final_db_files <- file.path(dirname(final_db_prefix), basename(staged_db_files))
  dir.create(dirname(final_db_prefix), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(staged_db_files, final_db_files, overwrite = TRUE)
  if (!all(copied)) {
    stop("Failed to copy staged BLAST database files into the cache directory.", call. = FALSE)
  }
  paste0(final_db_prefix, c(".phr", ".pin", ".psq"))
}

.dnmb_merops_prepare_blastp_paths <- function(query_fasta, db_prefix, blast_out) {
  paths_with_spaces <- c(query_fasta, db_prefix, blast_out)
  if (!any(vapply(paths_with_spaces, .dnmb_merops_has_space, logical(1)))) {
    return(list(
      staged = FALSE,
      cleanup = function() NULL,
      query_fasta = query_fasta,
      db_prefix = db_prefix,
      blast_out = blast_out
    ))
  }

  staging_dir <- tempfile("dnmb-merops-blastp-")
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  staged_query <- .dnmb_merops_stage_file(query_fasta, staging_dir, target_name = basename(query_fasta))
  prefix_files <- .dnmb_merops_prefix_files(db_prefix)
  for (path in prefix_files) {
    .dnmb_merops_stage_file(path, staging_dir, target_name = basename(path))
  }
  staged_prefix <- file.path(staging_dir, basename(db_prefix))
  staged_out <- file.path(staging_dir, basename(blast_out))
  list(
    staged = TRUE,
    cleanup = function() unlink(staging_dir, recursive = TRUE, force = TRUE),
    query_fasta = staged_query,
    db_prefix = staged_prefix,
    blast_out = staged_out
  )
}

.dnmb_merops_makeblastdb_args <- function(fasta_path, db_prefix, capabilities = .dnmb_makeblastdb_capabilities()) {
  args <- c("-in", fasta_path, "-dbtype", "prot")
  if (isTRUE(capabilities$supports_blastdb_version)) {
    args <- c(args, "-blastdb_version", "4")
  }
  if (isTRUE(capabilities$supports_parse_seqids)) {
    args <- c(args, "-parse_seqids")
  }
  c(args, "-out", db_prefix)
}

.dnmb_merops_makeblastdb_summary <- function(capabilities, args = character()) {
  paste(
    c(
      paste0("path=", if (isTRUE(capabilities$found) && nzchar(capabilities$path)) capabilities$path else "missing"),
      paste0("version=", capabilities$version %||% "unknown"),
      paste0("supports_parse_seqids=", isTRUE(capabilities$supports_parse_seqids)),
      paste0("supports_blastdb_version=", isTRUE(capabilities$supports_blastdb_version)),
      if (length(args)) paste0("args=", .dnmb_format_command("makeblastdb", args))
    ),
    collapse = "; "
  )
}

.dnmb_merops_makeblastdb_failure_detail <- function(run, capabilities, args) {
  stderr_text <- .dnmb_compact_output(c(run$stderr, run$stdout), max_lines = 6L)
  if (!nzchar(stderr_text) && nzchar(run$error %||% "")) {
    stderr_text <- .dnmb_compact_output(run$error, max_lines = 6L)
  }

  paste(
    c(
      paste0("path=", if (isTRUE(capabilities$found) && nzchar(capabilities$path)) capabilities$path else "missing"),
      paste0("version=", capabilities$version %||% "unknown"),
      paste0("args=", .dnmb_format_command("makeblastdb", args)),
      if (nzchar(stderr_text)) paste0("stderr=", stderr_text)
    ),
    collapse = "; "
  )
}

.dnmb_merops_install_detail <- function(layout, ok, prepare_result = NULL) {
  if (isTRUE(ok)) {
    return(layout$module_dir)
  }
  if (!file.exists(layout$fasta_path)) {
    return(paste0("Required asset missing: ", layout$fasta_path))
  }
  if (!.dnmb_merops_is_fasta(layout$fasta_path)) {
    return(paste0("MEROPS FASTA missing or malformed: ", layout$fasta_path))
  }
  if (!is.null(prepare_result) && !isTRUE(prepare_result$ok) && nzchar(prepare_result$detail %||% "")) {
    return(prepare_result$detail)
  }
  paste0("MEROPS install incomplete: ", layout$module_dir)
}

.dnmb_merops_manifest_diagnostics <- function(layout, prepare_result = NULL) {
  capabilities <- prepare_result$capabilities %||% .dnmb_makeblastdb_capabilities()
  args <- prepare_result$args %||% .dnmb_merops_makeblastdb_args(
    fasta_path = layout$fasta_path,
    db_prefix = layout$blast_db_prefix,
    capabilities = capabilities
  )

  list(
    makeblastdb_path = if (isTRUE(capabilities$found) && nzchar(capabilities$path)) capabilities$path else "",
    makeblastdb_version = capabilities$version %||% NA_character_,
    makeblastdb_args = as.character(args),
    blast_capability_supports_blastdb_version = isTRUE(capabilities$supports_blastdb_version),
    blast_capability_supports_parse_seqids = isTRUE(capabilities$supports_parse_seqids)
  )
}

.dnmb_merops_prepare_blast_db <- function(fasta_path, db_prefix, force = FALSE, trace_log = NULL) {
  db_files <- paste0(db_prefix, c(".phr", ".pin", ".psq"))
  capabilities <- .dnmb_makeblastdb_capabilities()
  run_paths <- .dnmb_merops_prepare_makeblastdb_paths(fasta_path = fasta_path, db_prefix = db_prefix)
  on.exit(run_paths$cleanup(), add = TRUE)
  args <- .dnmb_merops_makeblastdb_args(
    fasta_path = run_paths$fasta_path,
    db_prefix = run_paths$db_prefix,
    capabilities = capabilities
  )
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_merops_trace(trace_log, sprintf("[%s] makeblastdb %s", Sys.time(), .dnmb_merops_makeblastdb_summary(capabilities, args)))
  }

  if (!isTRUE(force) && all(file.exists(db_files))) {
    return(list(
      ok = TRUE,
      status = "cached",
      detail = paste("BLAST database files already exist", .dnmb_merops_makeblastdb_summary(capabilities, args), sep = "; "),
      blast_db_prefix = db_prefix,
      blast_db_files = db_files,
      command = NULL,
      capabilities = capabilities,
      args = args
    ))
  }

  dir.create(dirname(db_prefix), recursive = TRUE, showWarnings = FALSE)
  run <- dnmb_run_external(
    command = "makeblastdb",
    args = args,
    required = FALSE
  )
  staged_ok <- isTRUE(run$ok) && all(file.exists(run_paths$blast_db_files))
  if (staged_ok) {
    db_files <- .dnmb_merops_finalize_makeblastdb_paths(run_paths, final_db_prefix = db_prefix)
  }
  ok <- isTRUE(run$ok) && all(file.exists(db_files))
  status <- if (ok) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed"
  detail <- if (ok) {
    paste(db_prefix, .dnmb_merops_makeblastdb_summary(capabilities, args), sep = "; ")
  } else {
    .dnmb_merops_makeblastdb_failure_detail(run, capabilities, args)
  }
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_merops_trace(trace_log, sprintf("[%s] makeblastdb_result status=%s detail=%s", Sys.time(), status, detail))
  }

  list(
    ok = ok,
    status = status,
    detail = detail,
    blast_db_prefix = db_prefix,
    blast_db_files = db_files[file.exists(db_files)],
    command = run,
    capabilities = capabilities,
    args = args
  )
}

.dnmb_merops_subject_accession <- function(x) {
  value <- trimws(as.character(x))
  value[nchar(value) == 0L] <- NA_character_
  value <- sub("^gnl\\|", "", value)
  value <- sub("^lcl\\|", "", value)
  value <- sub("[[:space:]].*$", "", value)
  value <- sub("\\|.*$", "", value)
  value[nchar(value) == 0L] <- NA_character_
  value
}

.dnmb_merops_family_token <- function(values) {
  values <- trimws(as.character(values))
  values <- values[nzchar(values)]
  if (!length(values)) {
    return(NA_character_)
  }

  upper <- toupper(values)
  matched <- upper[grepl("^[ACGIMPNSTUZIQJLSR][0-9]{1,3}(?:\\.[0-9]+)?[A-Z]?$", upper)]
  if (!length(matched)) {
    return(NA_character_)
  }
  matched[[1]]
}

.dnmb_merops_family_from_title <- function(title) {
  title <- trimws(as.character(title)[1])
  if (is.na(title) || !nzchar(title)) {
    return(NA_character_)
  }

  anchored <- regmatches(
    title,
    gregexpr("#([ACGIMPNSTUZIQJLSR][0-9]{1,3}[A-Z]?)#", title, perl = TRUE)
  )[[1]]
  if (length(anchored)) {
    anchored <- gsub("^#|#$", "", anchored)
    anchored <- anchored[nzchar(anchored)]
    if (length(anchored)) {
      return(anchored[[1]])
    }
  }

  tokens <- unlist(strsplit(gsub("[^[:alnum:].#]+", " ", title), "\\s+", perl = TRUE), use.names = FALSE)
  .dnmb_merops_family_token(tokens)
}

dnmb_merops_parse_family_mapping <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(tibble::tibble(
      subject_accession = character(),
      family_id = character(),
      raw_line = character()
    ))
  }

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  if (!length(lines)) {
    return(tibble::tibble(
      subject_accession = character(),
      family_id = character(),
      raw_line = character()
    ))
  }

  parsed <- lapply(lines, function(line) {
    fields <- trimws(strsplit(line, "\t", fixed = FALSE)[[1]])
    fields <- fields[nzchar(fields)]
    if (!length(fields)) {
      return(NULL)
    }

    accession <- .dnmb_merops_subject_accession(fields[[1]])
    family_id <- .dnmb_merops_family_token(fields)
    if (is.na(accession) || is.na(family_id)) {
      return(NULL)
    }

    tibble::tibble(
      subject_accession = accession,
      family_id = family_id,
      raw_line = line
    )
  })

  mapping <- dplyr::bind_rows(parsed)
  if (!nrow(mapping)) {
    return(tibble::tibble(
      subject_accession = character(),
      family_id = character(),
      raw_line = character()
    ))
  }

  mapping <- mapping[!duplicated(mapping$subject_accession), , drop = FALSE]
  rownames(mapping) <- NULL
  mapping
}

dnmb_merops_install_module <- function(version = .dnmb_merops_default_version(),
                                       cache_root = NULL,
                                       base_url = .dnmb_merops_default_base_url(),
                                       asset_urls = NULL,
                                       force = FALSE,
                                       prepare = TRUE) {
  module <- .dnmb_merops_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_merops_asset_layout(module_dir)
  install_trace_log <- file.path(module_dir, "merops_install_trace.log")
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  urls <- .dnmb_merops_asset_urls(base_url = base_url, version = version, asset_urls = asset_urls)
  status <- .dnmb_merops_empty_status()
  remote_asset_state <- .dnmb_merops_remote_asset_state(
    urls = urls,
    enabled = .dnmb_merops_is_current_version(version) && is.null(asset_urls)
  )
  remote_refresh_needed <- .dnmb_merops_remote_update_needed(manifest, remote_asset_state)
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
      .dnmb_merops_status_row(
        "merops_remote_state",
        if (all(remote_asset_state$ok %in% TRUE)) {
          if (remote_refresh_needed) "update_detected" else "ok"
        } else {
          "unavailable"
        },
        remote_detail
      )
    )
    .dnmb_merops_trace(install_trace_log, sprintf("[%s] remote_asset_state %s", Sys.time(), remote_detail))
  }
  if (isTRUE(remote_refresh_needed) && !isTRUE(force)) {
    force <- TRUE
    .dnmb_merops_trace(install_trace_log, sprintf("[%s] forcing refresh because current_release metadata changed", Sys.time()))
  }

  required_assets <- c("pepunit.lib")
  required_ready <- all(file.exists(file.path(module_dir, required_assets)))
  manifest_ready <- !is.null(manifest)
  if (manifest_ready && required_ready && !isTRUE(force)) {
    prepare_result <- if (isTRUE(prepare)) {
      .dnmb_merops_prepare_blast_db(layout$fasta_path, layout$blast_db_prefix, force = FALSE, trace_log = install_trace_log)
    } else {
      NULL
    }
    if (!is.null(prepare_result)) {
      status <- dplyr::bind_rows(status, .dnmb_merops_status_row("merops_prepare", prepare_result$status, prepare_result$detail))
    }
    cached_ok <- if (isTRUE(prepare) && !is.null(prepare_result)) {
      isTRUE(prepare_result$ok)
    } else {
      isTRUE(!is.null(manifest$install_ok) && manifest$install_ok) || .dnmb_merops_is_fasta(layout$fasta_path)
    }
    manifest_fields <- unclass(manifest)
    manifest_fields[c("module", "version", "module_dir", "manifest_path", "written_at")] <- NULL
    manifest_fields <- utils::modifyList(
      manifest_fields,
      c(
        .dnmb_merops_manifest_diagnostics(layout, prepare_result = prepare_result),
        list(
          blast_db_files = layout$blast_db_files[file.exists(layout$blast_db_files)],
          prepared_with_makeblastdb = isTRUE(!is.null(prepare_result) && prepare_result$status %in% c("ok", "cached") && length(prepare_result$blast_db_files) == 3L),
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
      default_version = .dnmb_merops_default_version(),
      cache_root = cache_root
    )
    manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
    return(list(
      ok = cached_ok,
      status = dplyr::bind_rows(status, .dnmb_merops_status_row("merops_install", if (cached_ok) "cached" else "failed", .dnmb_merops_install_detail(layout, cached_ok, prepare_result = prepare_result))),
      manifest = manifest,
      files = c(
        list(
          fasta_path = layout$fasta_path,
          blast_db_prefix = layout$blast_db_prefix,
          trace_log = install_trace_log
        ),
        if (file.exists(layout$mapping_path)) list(mapping_path = layout$mapping_path),
        stats::setNames(as.list(layout$blast_db_files[file.exists(layout$blast_db_files)]), basename(layout$blast_db_files[file.exists(layout$blast_db_files)]))
      )
    ))
  }

  download_results <- list()
  for (asset_name in names(urls)) {
    asset_path <- file.path(module_dir, asset_name)
    result <- .dnmb_download_asset(urls[[asset_name]], asset_path)
    download_results[[asset_name]] <- result
    required_asset <- asset_name %in% required_assets
    asset_ok <- isTRUE(result$ok) || (!required_asset && file.exists(asset_path))
    status <- dplyr::bind_rows(
      status,
      .dnmb_merops_status_row(
        paste0("merops_download:", asset_name),
        if (asset_ok) "ok" else if (required_asset) "failed" else "optional_missing",
        if (asset_ok) asset_path else (result$error %||% asset_name)
      )
    )
  }

  ok <- .dnmb_merops_is_fasta(layout$fasta_path)
  prepare_result <- NULL
  if (ok && isTRUE(prepare)) {
    prepare_result <- .dnmb_merops_prepare_blast_db(layout$fasta_path, layout$blast_db_prefix, force = isTRUE(force), trace_log = install_trace_log)
    status <- dplyr::bind_rows(status, .dnmb_merops_status_row("merops_prepare", prepare_result$status, prepare_result$detail))
    ok <- isTRUE(prepare_result$ok)
  }

  mapping_available <- file.exists(layout$mapping_path)
  mapping <- if (mapping_available) dnmb_merops_parse_family_mapping(layout$mapping_path) else NULL
  manifest_payload <- c(
    list(
      source = "merops",
      base_url = base_url,
      asset_urls = urls,
      asset_files = names(urls),
      required_assets = required_assets,
      primary_fasta = layout$fasta_path,
      family_mapping_path = if (mapping_available) layout$mapping_path else NA_character_,
      family_mapping_rows = if (mapping_available) nrow(mapping) else 0L,
      blast_db_prefix = layout$blast_db_prefix,
      blast_db_files = layout$blast_db_files[file.exists(layout$blast_db_files)],
      prepared_with_makeblastdb = isTRUE(!is.null(prepare_result) && prepare_result$status %in% c("ok", "cached") && length(prepare_result$blast_db_files) == 3L),
      mapping_available = mapping_available,
      remote_asset_state = remote_asset_state
    ),
    .dnmb_merops_manifest_diagnostics(layout, prepare_result = prepare_result),
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
    default_version = .dnmb_merops_default_version(),
    cache_root = cache_root
  )
  manifest <- readRDS(manifest_path)
  class(manifest) <- c("dnmb_db_manifest", class(manifest))

  list(
    ok = ok,
    status = dplyr::bind_rows(status, .dnmb_merops_status_row("merops_install", if (ok) "ok" else "failed", .dnmb_merops_install_detail(layout, ok, prepare_result = prepare_result))),
    manifest = manifest,
    files = c(
      list(
        fasta_path = layout$fasta_path,
        blast_db_prefix = layout$blast_db_prefix,
        trace_log = install_trace_log
      ),
      if (mapping_available) list(mapping_path = layout$mapping_path),
      stats::setNames(as.list(layout$blast_db_files[file.exists(layout$blast_db_files)]), basename(layout$blast_db_files[file.exists(layout$blast_db_files)]))
    )
  )
}

dnmb_merops_get_module <- function(version = .dnmb_merops_default_version(),
                                   cache_root = NULL,
                                   required = FALSE) {
  module <- .dnmb_merops_module_name()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = FALSE)
  layout <- .dnmb_merops_asset_layout(module_dir)
  ok <- !is.null(manifest) && .dnmb_merops_is_fasta(layout$fasta_path)

  if (!ok && isTRUE(required)) {
    stop("MEROPS module is not installed for version `", version, "`.", call. = FALSE)
  }

  list(
    ok = ok,
    module = module,
    version = version,
    module_dir = module_dir,
    manifest = manifest,
    files = c(
      list(
        fasta_path = layout$fasta_path,
        blast_db_prefix = layout$blast_db_prefix
      ),
      if (file.exists(layout$mapping_path)) list(mapping_path = layout$mapping_path),
      stats::setNames(as.list(layout$blast_db_files[file.exists(layout$blast_db_files)]), basename(layout$blast_db_files[file.exists(layout$blast_db_files)]))
    )
  )
}

.dnmb_merops_empty_blast_hits <- function() {
  tibble::tibble(
    query = character(),
    subject_id = character(),
    subject_accession = character(),
    hit_label = character(),
    pident = numeric(),
    length = integer(),
    qlen = integer(),
    slen = integer(),
    qstart = integer(),
    qend = integer(),
    sstart = integer(),
    send = integer(),
    evalue = numeric(),
    bitscore = numeric(),
    qcov = numeric(),
    scov = numeric(),
    family_id = character(),
    family_confirmed = logical(),
    family_mapping_source = character()
  )
}

dnmb_merops_parse_blast_tabular <- function(path, family_mapping = NULL) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(.dnmb_merops_empty_blast_hits())
  }

  hits <- tryCatch(
    utils::read.delim(
      path,
      header = FALSE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      fill = TRUE,
      comment.char = ""
    ),
    error = function(...) NULL
  )
  if (is.null(hits) || !nrow(hits)) {
    return(.dnmb_merops_empty_blast_hits())
  }

  while (ncol(hits) < 13L) {
    hits[[ncol(hits) + 1L]] <- NA_character_
  }
  hits <- hits[, seq_len(13L), drop = FALSE]
  names(hits) <- c(
    "query",
    "subject_id",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "subject_title"
  )

  hits$query <- .dnmb_module_clean_annotation_key(hits$query)
  hits$subject_id <- as.character(hits$subject_id)
  hits$subject_title <- trimws(as.character(hits$subject_title))
  hits$subject_title[!nzchar(hits$subject_title)] <- NA_character_
  hits$subject_accession <- .dnmb_merops_subject_accession(hits$subject_id)
  hits$hit_label <- ifelse(is.na(hits$subject_title), hits$subject_id, hits$subject_title)
  numeric_columns <- c("pident", "evalue", "bitscore")
  length_columns <- c("length", "qlen", "slen", "qstart", "qend", "sstart", "send")
  hits[numeric_columns] <- lapply(hits[numeric_columns], as.numeric)
  hits[length_columns] <- lapply(hits[length_columns], as.numeric)
  hits$qcov <- ifelse(is.na(hits$qlen) | hits$qlen <= 0L, NA_real_, hits$length / hits$qlen)
  hits$scov <- ifelse(is.na(hits$slen) | hits$slen <= 0L, NA_real_, hits$length / hits$slen)

  if (is.character(family_mapping) && length(family_mapping) == 1L && file.exists(family_mapping)) {
    family_mapping <- dnmb_merops_parse_family_mapping(family_mapping)
  }
  title_family_id <- vapply(hits$hit_label, .dnmb_merops_family_from_title, character(1))
  if (is.null(family_mapping) || !is.data.frame(family_mapping) || !nrow(family_mapping)) {
    hits$family_id <- title_family_id
    hits$family_confirmed <- !is.na(hits$family_id) & nzchar(hits$family_id)
    hits$family_mapping_source <- ifelse(hits$family_confirmed, "subject_title", "absent")
  } else {
    mapping <- tibble::as_tibble(family_mapping)
    mapping <- mapping[!duplicated(mapping$subject_accession), c("subject_accession", "family_id"), drop = FALSE]
    mapped <- mapping$family_id[match(hits$subject_accession, mapping$subject_accession)]
    hits$family_id <- ifelse(!is.na(mapped) & nzchar(mapped), as.character(mapped), title_family_id)
    hits$family_confirmed <- !is.na(hits$family_id) & nzchar(hits$family_id)
    hits$family_mapping_source <- ifelse(
      !is.na(mapped) & nzchar(mapped),
      "dnld_list",
      ifelse(hits$family_confirmed, "subject_title", "unmapped")
    )
  }

  hits <- hits[, names(.dnmb_merops_empty_blast_hits()), drop = FALSE]
  tibble::as_tibble(hits)
}

dnmb_merops_normalize_hits <- function(hits) {
  if (is.null(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  required_columns <- c("query", "family_id", "hit_label", "subject_accession", "evalue", "bitscore", "pident", "qcov", "scov", "family_confirmed", "family_mapping_source")
  missing_columns <- setdiff(required_columns, names(hits))
  if (length(missing_columns)) {
    stop("`hits` is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  support <- paste0(
    "subject_accession=", hits$subject_accession,
    "; evalue=", signif(hits$evalue, 4),
    "; bitscore=", signif(hits$bitscore, 4),
    "; pident=", sprintf("%.1f", hits$pident),
    "; qcov=", sprintf("%.3f", hits$qcov),
    "; scov=", sprintf("%.3f", hits$scov),
    "; family_confirmed=", ifelse(hits$family_confirmed, "TRUE", "FALSE"),
    "; family_mapping_source=", hits$family_mapping_source
  )

  out <- data.frame(
    query = .dnmb_module_clean_annotation_key(hits$query),
    source = rep("merops", nrow(hits)),
    family_system = rep("MEROPS", nrow(hits)),
    family_id = ifelse(hits$family_confirmed, as.character(hits$family_id), NA_character_),
    hit_label = as.character(hits$hit_label),
    enzyme_role = rep("protease", nrow(hits)),
    evidence_mode = rep("direct", nrow(hits)),
    substrate_label = rep(NA_character_, nrow(hits)),
    support = support,
    subject_accession = as.character(hits$subject_accession),
    pident = as.numeric(hits$pident),
    alignment_length = as.numeric(hits$length),
    query_length = as.numeric(hits$qlen),
    subject_length = as.numeric(hits$slen),
    query_start = as.numeric(hits$qstart),
    query_end = as.numeric(hits$qend),
    subject_start = as.numeric(hits$sstart),
    subject_end = as.numeric(hits$send),
    evalue = as.numeric(hits$evalue),
    bitscore = as.numeric(hits$bitscore),
    qcov = as.numeric(hits$qcov),
    scov = as.numeric(hits$scov),
    family_confirmed = as.logical(hits$family_confirmed),
    family_mapping_source = as.character(hits$family_mapping_source),
    typing_eligible = rep(TRUE, nrow(hits)),
    stringsAsFactors = FALSE
  )
  out <- out[, c(.dnmb_module_optional_long_columns(), setdiff(names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_merops_append_evalue_threshold <- function() {
  1e-2
}

.dnmb_merops_hits_for_output <- function(hits,
                                         evalue_threshold = .dnmb_merops_append_evalue_threshold()) {
  if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  out <- as.data.frame(hits, stringsAsFactors = FALSE)
  if (!"evalue" %in% names(out)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out$evalue <- suppressWarnings(as.numeric(out$evalue))
  out <- out[!is.na(out$evalue) & out$evalue < as.numeric(evalue_threshold)[1], , drop = FALSE]
  rownames(out) <- NULL
  if (!nrow(out)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out
}

.dnmb_merops_output_table <- function(genes,
                                      hits,
                                      evalue_threshold = .dnmb_merops_append_evalue_threshold()) {
  filtered_hits <- .dnmb_merops_hits_for_output(hits, evalue_threshold = evalue_threshold)
  out <- .dnmb_module_output_table(genes = genes, hits = filtered_hits)
  drop_cols <- intersect(
    c(
      "enzyme_role",
      "evidence_mode",
      "substrate_label",
      "support",
      "typing_eligible",
      "family_confirmed",
      "family_mapping_source"
    ),
    names(out)
  )
  if (length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  out
}

.dnmb_merops_run_diamond <- function(query_fasta, reference_fasta, output_dir, cpu = 1L, trace_log = NULL) {
  blast_out <- file.path(output_dir, "merops_blastp.tsv")
  # diamond -d takes a prefix and appends .dmnd automatically
  dmnd_prefix <- file.path(output_dir, "merops_reference")
  dmnd_file <- paste0(dmnd_prefix, ".dmnd")

  # Check diamond availability
  diamond_check <- dnmb_run_external("diamond", args = "version", required = FALSE)
  if (!nzchar(diamond_check$resolved_command)) {
    return(list(ok = FALSE, resolved_command = "", error = "diamond not found", blast_out = blast_out))
  }

  # diamond makedb (uses staging dir to avoid spaces in paths)
  if (!file.exists(dmnd_file)) {
    stage <- tempfile("dnmb-merops-diamond-")
    dir.create(stage, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
    stage_ref <- file.path(stage, "ref.faa")
    stage_db <- file.path(stage, "ref")
    # Clean reference FASTA: remove spaces from sequence lines (MEROPS has them)
    ref_lines <- readLines(reference_fasta, warn = FALSE)
    ref_lines <- ifelse(startsWith(ref_lines, ">"), ref_lines, gsub(" ", "", ref_lines, useBytes = TRUE))
    writeLines(ref_lines, stage_ref)
    makedb_args <- c("makedb", "--in", stage_ref, "-d", stage_db)
    if (!is.null(trace_log)) .dnmb_merops_trace(trace_log, sprintf("[%s] diamond makedb", Sys.time()))
    makedb_run <- dnmb_run_external("diamond", args = makedb_args, required = FALSE)
    if (!isTRUE(makedb_run$ok)) {
      return(list(ok = FALSE, resolved_command = diamond_check$resolved_command, error = "diamond makedb failed", blast_out = blast_out))
    }
    file.copy(paste0(stage_db, ".dmnd"), dmnd_file, overwrite = TRUE)
  }

  # diamond blastp — same outfmt 6 columns as NCBI blastp
  # Use staging dir for query to avoid path issues
  stage2 <- tempfile("dnmb-merops-dsearch-")
  dir.create(stage2, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage2, recursive = TRUE, force = TRUE), add = TRUE)
  stage_query <- file.path(stage2, "query.faa")
  stage_out <- file.path(stage2, "hits.tsv")
  file.copy(query_fasta, stage_query, overwrite = TRUE)
  file.copy(dmnd_file, file.path(stage2, "ref.dmnd"), overwrite = TRUE)
  search_args <- c(
    "blastp",
    "-q", stage_query,
    "-d", file.path(stage2, "ref"),
    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
    "-k", "5",
    "--threads", as.character(as.integer(cpu)),
    "--sensitive",
    "-o", stage_out
  )
  if (!is.null(trace_log)) .dnmb_merops_trace(trace_log, sprintf("[%s] diamond blastp --sensitive threads=%s", Sys.time(), cpu))
  search_run <- dnmb_run_external("diamond", args = search_args, required = FALSE)
  if (isTRUE(search_run$ok) && file.exists(stage_out)) {
    file.copy(stage_out, blast_out, overwrite = TRUE)
  }

  list(
    ok = isTRUE(search_run$ok),
    resolved_command = diamond_check$resolved_command,
    error = if (!isTRUE(search_run$ok)) paste(search_run$stderr, collapse = " ") else NULL,
    blast_out = blast_out
  )
}

dnmb_merops_run_blastp <- function(query_fasta,
                                   output_dir,
                                   version = .dnmb_merops_default_version(),
                                   cache_root = NULL,
                                   cpu = 1L,
                                   install = FALSE,
                                   base_url = .dnmb_merops_default_base_url(),
                                   asset_urls = NULL,
                                   search_backend = c("diamond", "blastp")) {
  search_backend <- match.arg(search_backend)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_merops_empty_status()
  trace_log <- file.path(output_dir, "merops_trace.log")
  .dnmb_merops_trace(trace_log, sprintf("[%s] enter dnmb_merops_run_blastp query=%s backend=%s", Sys.time(), query_fasta, search_backend))

  module <- dnmb_merops_get_module(version = version, cache_root = cache_root, required = FALSE)
  if (!isTRUE(module$ok) && isTRUE(install)) {
    install_result <- dnmb_merops_install_module(
      version = version,
      cache_root = cache_root,
      base_url = base_url,
      asset_urls = asset_urls,
      force = FALSE,
      prepare = TRUE
    )
    status <- dplyr::bind_rows(status, install_result$status)
    module <- dnmb_merops_get_module(version = version, cache_root = cache_root, required = FALSE)
  }

  if (!isTRUE(module$ok)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_merops_status_row("merops_module", "missing", paste0("MEROPS module not installed for version ", version))),
      files = list(query_fasta = query_fasta),
      command = NULL,
      raw_hits = .dnmb_merops_empty_blast_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      manifest = module$manifest
    ))
  }

  # --- diamond backend ---
  if (identical(search_backend, "diamond")) {
    diamond_result <- .dnmb_merops_run_diamond(
      query_fasta = query_fasta,
      reference_fasta = module$files$fasta_path,
      output_dir = output_dir,
      cpu = cpu,
      trace_log = trace_log
    )
    if (isTRUE(diamond_result$ok) || !nzchar(diamond_result$resolved_command %||% "")) {
      # diamond succeeded or diamond not found — return result
      # If diamond not found, fall through to blastp
      if (isTRUE(diamond_result$ok) || nzchar(diamond_result$resolved_command %||% "")) {
        blast_out <- diamond_result$blast_out
        status <- dplyr::bind_rows(status, .dnmb_merops_status_row("merops_prepare", "diamond", "diamond makedb"))
        family_mapping <- if (!is.null(module$files$mapping_path) && file.exists(module$files$mapping_path)) {
          dnmb_merops_parse_family_mapping(module$files$mapping_path)
        } else {
          NULL
        }
        raw_hits <- if (file.exists(blast_out)) {
          dnmb_merops_parse_blast_tabular(blast_out, family_mapping = family_mapping)
        } else {
          .dnmb_merops_empty_blast_hits()
        }
        hits <- dnmb_merops_normalize_hits(raw_hits)
        status <- dplyr::bind_rows(
          status,
          .dnmb_merops_status_row(
            "diamond_blastp",
            if (isTRUE(diamond_result$ok)) if (nrow(raw_hits)) "ok" else "empty" else "failed",
            if (isTRUE(diamond_result$ok)) blast_out else (diamond_result$error %||% "diamond failed")
          )
        )
        return(list(
          ok = isTRUE(diamond_result$ok),
          status = status,
          files = c(list(query_fasta = query_fasta, blast_tsv = blast_out, trace_log = trace_log), module$files),
          command = diamond_result,
          raw_hits = raw_hits,
          hits = hits,
          manifest = module$manifest
        ))
      }
    }
    # diamond binary not found — fall through to blastp
    .dnmb_merops_trace(trace_log, sprintf("[%s] diamond not found, falling back to blastp", Sys.time()))
  }

  # --- blastp backend ---
  prepare_result <- .dnmb_merops_prepare_blast_db(
    fasta_path = module$files$fasta_path,
    db_prefix = module$files$blast_db_prefix,
    force = FALSE,
    trace_log = trace_log
  )
  status <- dplyr::bind_rows(status, .dnmb_merops_status_row("merops_prepare", prepare_result$status, prepare_result$detail))
  if (!isTRUE(prepare_result$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = c(list(query_fasta = query_fasta), module$files),
      command = prepare_result$command,
      raw_hits = .dnmb_merops_empty_blast_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      manifest = module$manifest
    ))
  }

  blast_out <- file.path(output_dir, "merops_blastp.tsv")
  stdout_log <- file.path(output_dir, "merops_blastp.stdout.log")
  stderr_log <- file.path(output_dir, "merops_blastp.stderr.log")
  run_paths <- .dnmb_merops_prepare_blastp_paths(
    query_fasta = query_fasta,
    db_prefix = module$files$blast_db_prefix,
    blast_out = blast_out
  )
  on.exit(run_paths$cleanup(), add = TRUE)
  args <- c(
    "-query", run_paths$query_fasta,
    "-db", run_paths$db_prefix,
    "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle",
    "-max_target_seqs", "5",
    "-num_threads", as.character(as.integer(cpu)),
    "-out", run_paths$blast_out
  )
  .dnmb_merops_trace(trace_log, sprintf("[%s] blastp args=%s", Sys.time(), .dnmb_format_command("blastp", args)))
  command <- dnmb_run_external("blastp", args = args, required = FALSE)
  writeLines(command$stdout, con = stdout_log)
  writeLines(command$stderr, con = stderr_log)
  if (isTRUE(command$ok) && isTRUE(run_paths$staged) && file.exists(run_paths$blast_out)) {
    file.copy(run_paths$blast_out, blast_out, overwrite = TRUE)
  }

  family_mapping <- if (!is.null(module$files$mapping_path) && file.exists(module$files$mapping_path)) {
    dnmb_merops_parse_family_mapping(module$files$mapping_path)
  } else {
    NULL
  }
  raw_hits <- if (file.exists(blast_out)) {
    dnmb_merops_parse_blast_tabular(blast_out, family_mapping = family_mapping)
  } else {
    .dnmb_merops_empty_blast_hits()
  }
  hits <- dnmb_merops_normalize_hits(raw_hits)
  status <- dplyr::bind_rows(
    status,
    .dnmb_merops_status_row(
      "blastp",
      if (isTRUE(command$ok)) if (nrow(raw_hits)) "ok" else "empty" else if (!nzchar(command$resolved_command)) "missing" else "failed",
      if (isTRUE(command$ok)) blast_out else (command$error %||% "blastp failed")
    )
  )

  list(
    ok = isTRUE(command$ok),
    status = status,
    files = c(
      list(
        query_fasta = query_fasta,
        blast_tsv = blast_out,
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

dnmb_run_merops_module <- function(genes,
                                   output_dir,
                                   version = .dnmb_merops_default_version(),
                                   cache_root = NULL,
                                   install = TRUE,
                                   base_url = .dnmb_merops_default_base_url(),
                                   asset_urls = NULL,
                                   cpu = 1L,
                                   genbank = NULL,
                                   search_backend = c("diamond", "blastp")) {
  search_backend <- match.arg(search_backend)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- file.path(output_dir, "merops_module_trace.log")
  fasta_path <- file.path(output_dir, "merops_query_proteins.faa")
  existing_faa <- dnmb_resolve_query_faa(genbank = genbank, output_dir = output_dir, fallback_filename = basename(fasta_path))
  if (!is.null(existing_faa) && .dnmb_can_reuse_query_fasta(existing_faa, genes)) {
    proteins <- .dnmb_prepare_query_proteins(genes)
    fasta <- list(path = existing_faa, n = nrow(proteins), proteins = proteins)
    fasta_path <- existing_faa
  } else {
    fasta <- .dnmb_write_query_fasta(genes, fasta_path)
  }
  status <- .dnmb_merops_status_row(
    "merops_query_fasta",
    if (fasta$n) "ok" else "empty",
    paste0("proteins=", fasta$n)
  )
  .dnmb_merops_trace(trace_log, sprintf("[%s] query_faa=%s existing=%s proteins=%s", Sys.time(), fasta_path, !is.null(existing_faa), fasta$n))

  if (!fasta$n) {
    return(list(
      ok = TRUE,
      status = status,
      files = list(query_fasta = fasta_path, trace_log = trace_log),
      manifest = NULL,
      raw_hits = .dnmb_merops_empty_blast_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      query_proteins = fasta$proteins
    ))
  }

  search <- dnmb_merops_run_blastp(
    query_fasta = fasta_path,
    output_dir = output_dir,
    version = version,
    cache_root = cache_root,
    cpu = cpu,
    install = install,
    base_url = base_url,
    asset_urls = asset_urls,
    search_backend = search_backend
  )

  list(
    ok = isTRUE(search$ok),
    status = dplyr::bind_rows(status, search$status),
    files = c(search$files, list(module_trace_log = trace_log)),
    manifest = search$manifest,
    raw_hits = search$raw_hits,
    hits = search$hits,
    query_proteins = fasta$proteins
  )
}
