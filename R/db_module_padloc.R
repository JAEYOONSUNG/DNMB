.dnmb_padloc_module_name <- function() {
  "padloc"
}

.dnmb_padloc_default_version <- function() {
  "current"
}

.dnmb_padloc_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = as.character(component)[1],
    status = as.character(status)[1],
    detail = as.character(detail)[1]
  )
}

.dnmb_padloc_empty_status <- function() {
  tibble::tibble(
    component = character(),
    status = character(),
    detail = character()
  )
}

.dnmb_padloc_trace <- function(path, text) {
  cat(paste0(text, "\n"), file = path, append = TRUE)
}

.dnmb_padloc_asset_layout <- function(module_dir) {
  list(
    module_dir = module_dir,
    data_dir = file.path(module_dir, "data"),
    cli_path = file.path(module_dir, "bin", "padloc"),
    install_trace_log = file.path(module_dir, "padloc_install_trace.log")
  )
}

.dnmb_padloc_default_db_repo <- function() {
  "https://github.com/padlocbio/padloc-db"
}

.dnmb_padloc_copy_dir_contents <- function(src_dir, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(src_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  if (!length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- file.copy(entries, dest_dir, recursive = TRUE, copy.mode = TRUE)
  if (!all(ok)) {
    stop("Failed to copy PADLOC assets from ", src_dir, " to ", dest_dir, call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_padloc_compile_db <- function(data_dir) {
  hmm_dir <- file.path(data_dir, "hmm")
  if (dir.exists(hmm_dir) && !file.exists(file.path(hmm_dir, "padlocdb.hmm"))) {
    hmm_files <- list.files(hmm_dir, pattern = "\\.hmm$", full.names = TRUE)
    if (length(hmm_files)) {
      compiled_hmm <- file.path(hmm_dir, "padlocdb.hmm")
      con <- file(compiled_hmm, open = "w")
      on.exit(close(con), add = TRUE)
      for (hmm_file in hmm_files) {
        writeLines(readLines(hmm_file, warn = FALSE), con = con)
      }
    }
  }
  cm_dir <- file.path(data_dir, "cm")
  if (dir.exists(cm_dir) && !file.exists(file.path(cm_dir, "padlocdb.cm"))) {
    cm_files <- list.files(cm_dir, pattern = "\\.cm$", full.names = TRUE)
    if (length(cm_files)) {
      compiled_cm <- file.path(cm_dir, "padlocdb.cm")
      con_cm <- file(compiled_cm, open = "w")
      on.exit(close(con_cm), add = TRUE)
      for (cm_file in cm_files) {
        writeLines(readLines(cm_file, warn = FALSE), con = con_cm)
      }
    }
  }
  invisible(data_dir)
}

.dnmb_padloc_release_info <- function(db_repo = .dnmb_padloc_default_db_repo()) {
  api_url <- "https://api.github.com/repos/padlocbio/padloc-db/releases/latest"
  tag_name <- tryCatch({
    payload <- jsonlite::fromJSON(api_url)
    as.character(payload$tag_name %||% "")
  }, error = function(e) "")
  if (!nzchar(tag_name)) {
    tag_name <- "v1.4.0"
  }
  list(
    tag_name = tag_name,
    archive_url = paste0(db_repo, "/archive/refs/tags/", tag_name, ".tar.gz")
  )
}

.dnmb_padloc_locate_release_root <- function(extract_dir) {
  candidates <- list.dirs(extract_dir, recursive = TRUE, full.names = TRUE)
  matches <- candidates[
    dir.exists(file.path(candidates, "hmm")) &
      dir.exists(file.path(candidates, "sys"))
  ]
  if (!length(matches)) {
    return("")
  }
  matches[[1]]
}

.dnmb_padloc_install_db_from_release <- function(layout) {
  info <- .dnmb_padloc_release_info()
  stage_dir <- tempfile("dnmb-padloc-db-")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  archive_path <- file.path(stage_dir, paste0(info$tag_name, ".tar.gz"))
  download <- .dnmb_download_asset(info$archive_url, archive_path, insecure = FALSE)
  if (!isTRUE(download$ok) || !file.exists(archive_path)) {
    return(list(ok = FALSE, detail = download$error %||% info$archive_url, tag_name = info$tag_name))
  }
  untar_ok <- tryCatch({
    utils::untar(archive_path, exdir = stage_dir)
    TRUE
  }, error = function(e) FALSE)
  if (!untar_ok) {
    return(list(ok = FALSE, detail = archive_path, tag_name = info$tag_name))
  }
  root_dir <- .dnmb_padloc_locate_release_root(stage_dir)
  if (!nzchar(root_dir)) {
    return(list(ok = FALSE, detail = stage_dir, tag_name = info$tag_name))
  }
  .dnmb_padloc_copy_dir_contents(root_dir, layout$data_dir)
  .dnmb_padloc_compile_db(layout$data_dir)
  list(ok = .dnmb_padloc_has_db(layout), detail = layout$data_dir, tag_name = info$tag_name)
}

.dnmb_padloc_default_data_dir <- function(runner) {
  command <- path.expand(as.character(runner$command %||% "")[1])
  if (!nzchar(command) || !file.exists(command)) {
    return("")
  }
  candidate <- normalizePath(file.path(dirname(command), "..", "data"), winslash = "/", mustWork = FALSE)
  if (dir.exists(candidate)) candidate else ""
}

.dnmb_padloc_candidate_roots <- function() {
  conda_detect <- tryCatch(dnmb_detect_binary("conda", required = FALSE)$path, error = function(e) "")
  conda_exe <- Sys.getenv("CONDA_EXE", unset = "")
  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  roots <- c(
    if (nzchar(conda_detect)) dirname(dirname(conda_detect)) else NULL,
    if (nzchar(conda_exe)) dirname(dirname(conda_exe)) else NULL,
    if (nzchar(conda_prefix)) conda_prefix else NULL,
    "/opt/homebrew/Caskroom/miniforge/base",
    "/opt/biotools",
    path.expand("~/miniforge3"),
    path.expand("~/miniconda3"),
    path.expand("~/mambaforge"),
    path.expand("~/Library/r-miniconda-arm64")
  )
  roots <- unique(path.expand(as.character(roots)))
  roots[nzchar(roots)]
}

.dnmb_padloc_detect_runner <- function() {
  direct <- tryCatch(dnmb_detect_binary("padloc", required = FALSE)$path, error = function(e) "")
  if (nzchar(direct)) {
    return(list(command = direct, args = character(), detail = direct))
  }

  roots <- .dnmb_padloc_candidate_roots()
  path_candidates <- unique(c(
    file.path(roots, "envs", "padloc", "bin", "padloc"),
    file.path(roots, "bin", "padloc")
  ))
  path_candidates <- path_candidates[file.exists(path_candidates)]
  if (length(path_candidates)) {
    return(list(command = path_candidates[[1]], args = character(), detail = path_candidates[[1]]))
  }

  conda_candidates <- unique(c(
    tryCatch(dnmb_detect_binary("conda", required = FALSE)$path, error = function(e) ""),
    Sys.getenv("CONDA_EXE", unset = ""),
    file.path(roots, "bin", "conda"),
    file.path(roots, "condabin", "conda")
  ))
  conda_candidates <- conda_candidates[nzchar(conda_candidates) & file.exists(conda_candidates)]
  if (length(conda_candidates)) {
    return(list(
      command = conda_candidates[[1]],
      args = c("run", "-n", "padloc", "padloc"),
      detail = paste(c(conda_candidates[[1]], "run", "-n", "padloc", "padloc"), collapse = " ")
    ))
  }

  list(command = "", args = character(), detail = "")
}

.dnmb_padloc_has_db <- function(layout) {
  file.exists(file.path(layout$data_dir, "hmm", "padlocdb.hmm")) &&
    dir.exists(file.path(layout$data_dir, "sys"))
}

dnmb_padloc_install_module <- function(version = .dnmb_padloc_default_version(),
                                       cache_root = NULL,
                                       install = TRUE) {
  module <- .dnmb_padloc_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_padloc_asset_layout(module_dir)
  status <- .dnmb_padloc_empty_status()
  runner <- .dnmb_padloc_detect_runner()

  if (!nzchar(runner$command)) {
    return(list(
      ok = FALSE,
      status = .dnmb_padloc_status_row(
        "padloc_cli",
        "missing",
        "PADLOC executable not found. Install PADLOC or use the DNMB Docker image."
      ),
      files = list(),
      manifest = NULL
    ))
  }

  .dnmb_write_exec_wrapper(layout$cli_path, runner$command, runner$args)
  status <- dplyr::bind_rows(
    status,
    .dnmb_padloc_status_row("padloc_cli", "ok", runner$detail)
  )

  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  if (!is.null(manifest) &&
      isTRUE(manifest$install_ok) &&
      file.exists(layout$cli_path) &&
      .dnmb_padloc_has_db(layout)) {
    return(list(
      ok = TRUE,
      status = dplyr::bind_rows(status, .dnmb_padloc_status_row("padloc_install", "cached", module_dir)),
      files = list(cli = layout$cli_path, data_dir = layout$data_dir, trace_log = layout$install_trace_log),
      manifest = manifest
    ))
  }

  if (!isTRUE(install) && !.dnmb_padloc_has_db(layout)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_padloc_status_row("padloc_install", "missing", "PADLOC DB is missing and module_install is FALSE.")),
      files = list(cli = layout$cli_path, data_dir = layout$data_dir, trace_log = layout$install_trace_log),
      manifest = manifest
    ))
  }

  dir.create(layout$data_dir, recursive = TRUE, showWarnings = FALSE)
  source_data_dir <- .dnmb_padloc_default_data_dir(runner)
  copy_ok <- FALSE
  db_version <- NA_character_
  if (nzchar(source_data_dir) &&
      file.exists(file.path(source_data_dir, "data.tar.gz")) &&
      !.dnmb_padloc_has_db(layout)) {
    .dnmb_padloc_trace(
      layout$install_trace_log,
      sprintf("[%s] PADLOC DB source detected but uncompiled: %s", Sys.time(), source_data_dir)
    )
  }
  if (nzchar(source_data_dir) &&
      dir.exists(source_data_dir) &&
      file.exists(file.path(source_data_dir, "hmm", "padlocdb.hmm")) &&
      dir.exists(file.path(source_data_dir, "sys"))) {
    .dnmb_padloc_trace(
      layout$install_trace_log,
      sprintf("[%s] Copying PADLOC DB from %s", Sys.time(), source_data_dir)
    )
    .dnmb_padloc_copy_dir_contents(source_data_dir, layout$data_dir)
    copy_ok <- .dnmb_padloc_has_db(layout)
  }
  if (!copy_ok) {
    .dnmb_padloc_trace(
      layout$install_trace_log,
      sprintf("[%s] Downloading PADLOC DB release from GitHub", Sys.time())
    )
    release_install <- .dnmb_padloc_install_db_from_release(layout)
    copy_ok <- isTRUE(release_install$ok)
    db_version <- release_install$tag_name %||% NA_character_
    if (!copy_ok) {
      .dnmb_padloc_trace(
        layout$install_trace_log,
        sprintf("[%s] PADLOC DB download failed: %s", Sys.time(), release_install$detail %||% "unknown error")
      )
    }
  }

  ok <- copy_ok && .dnmb_padloc_has_db(layout)
  status <- dplyr::bind_rows(
    status,
    .dnmb_padloc_status_row(
      "padloc_db",
      if (ok) "ok" else "failed",
      if (ok) layout$data_dir else layout$data_dir
    )
  )
  if (!ok) {
    return(list(
      ok = FALSE,
      status = status,
      files = list(cli = layout$cli_path, data_dir = layout$data_dir, trace_log = layout$install_trace_log),
      manifest = manifest
    ))
  }

  if (is.na(db_version) || !nzchar(db_version)) {
    version_run <- dnmb_run_external(layout$cli_path, args = c("--db-version", "--data", layout$data_dir), required = FALSE)
    db_version <- .dnmb_parse_tool_version(c(version_run$stdout, version_run$stderr))
  }
  manifest <- list(
    install_ok = TRUE,
    cli_path = layout$cli_path,
    data_dir = layout$data_dir,
    runner_command = runner$command,
    runner_args = runner$args,
    db_version = db_version
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_padloc_default_version(),
    cache_root = cache_root
  )

  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_padloc_status_row("padloc_install", "ok", module_dir)),
    files = list(cli = layout$cli_path, data_dir = layout$data_dir, trace_log = layout$install_trace_log),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_padloc_get_module <- function(version = .dnmb_padloc_default_version(),
                                   cache_root = NULL,
                                   required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_padloc_module_name(), version, cache_root = cache_root, required = FALSE)
  layout <- .dnmb_padloc_asset_layout(.dnmb_db_module_dir(.dnmb_padloc_module_name(), version, cache_root = cache_root, create = FALSE))
  ok <- !is.null(manifest) &&
    isTRUE(manifest$install_ok) &&
    file.exists(layout$cli_path) &&
    .dnmb_padloc_has_db(layout)
  if (required && !ok) {
    stop("PADLOC module is not installed for version `", version, "`.", call. = FALSE)
  }
  list(
    ok = ok,
    manifest = manifest,
    files = list(
      cli = layout$cli_path,
      data_dir = layout$data_dir,
      trace_log = layout$install_trace_log
    )
  )
}

dnmb_padloc_parse_hits <- function(path, id_map) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  tbl <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(tbl) || !"target.name" %in% names(tbl)) {
    return(data.frame())
  }

  map_tbl <- as.data.frame(id_map, stringsAsFactors = FALSE)
  map_tbl$query_id <- as.character(map_tbl$query_id)
  map_tbl$locus_tag <- as.character(map_tbl$locus_tag)
  tbl$query <- map_tbl$locus_tag[match(as.character(tbl[["target.name"]]), map_tbl$query_id)]
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$query)
  tbl <- tbl[!is.na(tbl$query) & nzchar(tbl$query), , drop = FALSE]
  rownames(tbl) <- NULL
  tbl
}

dnmb_padloc_normalize_hits <- function(hits_tbl) {
  if (is.null(hits_tbl) || !is.data.frame(hits_tbl) || !nrow(hits_tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  tbl <- as.data.frame(hits_tbl, stringsAsFactors = FALSE)
  support <- vapply(seq_len(nrow(tbl)), function(i) {
    parts <- c(
      if (!is.na(tbl[["system.number"]][[i]]) && nzchar(as.character(tbl[["system.number"]][[i]]))) paste0("system=", tbl[["system.number"]][[i]]) else NA_character_,
      if (!is.na(tbl[["hmm.accession"]][[i]]) && nzchar(as.character(tbl[["hmm.accession"]][[i]]))) paste0("hmm=", tbl[["hmm.accession"]][[i]]) else NA_character_,
      if (!is.na(tbl[["protein.name"]][[i]]) && nzchar(as.character(tbl[["protein.name"]][[i]]))) paste0("protein=", tbl[["protein.name"]][[i]]) else NA_character_,
      if (!is.na(tbl[["domain.iE.value"]][[i]])) paste0("domain_iE=", signif(as.numeric(tbl[["domain.iE.value"]][[i]]), 4)) else NA_character_,
      if (!is.na(tbl[["target.coverage"]][[i]])) paste0("target_cov=", sprintf("%.3f", as.numeric(tbl[["target.coverage"]][[i]]))) else NA_character_,
      if (!is.na(tbl[["hmm.coverage"]][[i]])) paste0("hmm_cov=", sprintf("%.3f", as.numeric(tbl[["hmm.coverage"]][[i]]))) else NA_character_
    )
    parts <- parts[!is.na(parts)]
    if (!length(parts)) NA_character_ else paste(parts, collapse = "; ")
  }, character(1))

  out <- data.frame(
    query = .dnmb_module_clean_annotation_key(tbl$query),
    source = "padloc",
    family_system = "PADLOC",
    family_id = as.character(tbl[["system"]]),
    hit_label = as.character(tbl[["protein.name"]]),
    enzyme_role = as.character(tbl[["hmm.name"]]),
    evidence_mode = "direct",
    substrate_label = NA_character_,
    support = support,
    typing_eligible = TRUE,
    system_number = suppressWarnings(as.integer(tbl[["system.number"]])),
    system = as.character(tbl[["system"]]),
    seqid = as.character(tbl[["seqid"]]),
    target_name = as.character(tbl[["target.name"]]),
    hmm_accession = as.character(tbl[["hmm.accession"]]),
    hmm_name = as.character(tbl[["hmm.name"]]),
    protein_name = as.character(tbl[["protein.name"]]),
    full_seq_evalue = suppressWarnings(as.numeric(tbl[["full.seq.E.value"]])),
    domain_ievalue = suppressWarnings(as.numeric(tbl[["domain.iE.value"]])),
    target_coverage = suppressWarnings(as.numeric(tbl[["target.coverage"]])),
    hmm_coverage = suppressWarnings(as.numeric(tbl[["hmm.coverage"]])),
    start = suppressWarnings(as.integer(tbl[["start"]])),
    end = suppressWarnings(as.integer(tbl[["end"]])),
    strand = as.character(tbl[["strand"]]),
    target_description = as.character(tbl[["target.description"]]),
    relative_position = suppressWarnings(as.integer(tbl[["relative.position"]])),
    contig_end = suppressWarnings(as.integer(tbl[["contig.end"]])),
    all_domains = as.character(tbl[["all.domains"]]),
    best_hits = as.character(tbl[["best.hits"]]),
    stringsAsFactors = FALSE
  )
  out <- out[, c(.dnmb_module_optional_long_columns(), setdiff(names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_padloc_hits_for_output <- function(hits) {
  if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out <- as.data.frame(hits, stringsAsFactors = FALSE)
  out$domain_ievalue <- suppressWarnings(as.numeric(out$domain_ievalue))
  out$target_coverage <- suppressWarnings(as.numeric(out$target_coverage))
  out$hmm_coverage <- suppressWarnings(as.numeric(out$hmm_coverage))
  out$system_number <- suppressWarnings(as.integer(out$system_number))
  out <- out[order(
    out$query,
    ifelse(is.na(out$domain_ievalue), Inf, out$domain_ievalue),
    -ifelse(is.na(out$target_coverage), -Inf, out$target_coverage),
    -ifelse(is.na(out$hmm_coverage), -Inf, out$hmm_coverage),
    out$system_number
  ), , drop = FALSE]
  out <- out[!duplicated(out$query), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_padloc_output_table <- function(genes, hits) {
  selected_hits <- .dnmb_padloc_hits_for_output(hits)
  out <- .dnmb_module_output_table(genes = genes, hits = selected_hits)
  drop_cols <- intersect(c("family_id", "hit_label", "enzyme_role", "evidence_mode", "substrate_label", "typing_eligible"), names(out))
  if (length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  ordered <- c(
    intersect(dnmb_backbone_columns(), names(out)),
    intersect(
      c("system", "system_number", "protein_name", "target_description", "hmm_accession", "hmm_name",
        "full_seq_evalue", "domain_ievalue", "target_coverage", "hmm_coverage", "relative_position",
        "seqid", "support"),
      names(out)
    ),
    setdiff(
      names(out),
      c(intersect(dnmb_backbone_columns(), names(out)),
        "system", "system_number", "protein_name", "target_description", "hmm_accession", "hmm_name",
        "full_seq_evalue", "domain_ievalue", "target_coverage", "hmm_coverage", "relative_position",
        "seqid", "support")
    )
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_padloc_module <- function(genes,
                                   output_dir,
                                   version = .dnmb_padloc_default_version(),
                                   cache_root = NULL,
                                   install = TRUE,
                                   cpu = 1L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- file.path(output_dir, "padloc_module_trace.log")
  status <- .dnmb_padloc_empty_status()

  install_result <- dnmb_padloc_install_module(
    version = version,
    cache_root = cache_root,
    install = install
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!isTRUE(install_result$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = c(install_result$files, list(trace_log = trace_log)),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame()
    ))
  }

  module <- dnmb_padloc_get_module(version = version, cache_root = cache_root, required = TRUE)
  input <- .dnmb_write_padloc_input(genes = genes, output_dir = output_dir, prefix = "padloc_query_proteins")
  status <- dplyr::bind_rows(
    status,
    .dnmb_padloc_status_row("padloc_input", if (nrow(input$map)) "ok" else "empty", paste0("proteins=", nrow(input$map)))
  )
  if (!nrow(input$map)) {
    return(list(
      ok = TRUE,
      status = status,
      files = c(module$files, list(faa = input$faa, gff = input$gff, id_map = input$map_path, trace_log = trace_log)),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame()
    ))
  }

  args <- c(
    "--faa", input$faa,
    "--gff", input$gff,
    "--outdir", output_dir,
    "--data", module$files$data_dir,
    "--cpu", as.character(as.integer(cpu)[1])
  )
  .dnmb_padloc_trace(trace_log, sprintf("[%s] %s", Sys.time(), .dnmb_format_command(module$files$cli, args)))
  run <- dnmb_run_external(module$files$cli, args = args, required = FALSE)
  if (length(run$stdout)) {
    .dnmb_padloc_trace(trace_log, paste(run$stdout, collapse = "\n"))
  }
  if (length(run$stderr)) {
    .dnmb_padloc_trace(trace_log, paste(run$stderr, collapse = "\n"))
  }

  csv_candidates <- list.files(output_dir, pattern = "_padloc\\.csv$", full.names = TRUE)
  csv_candidates <- csv_candidates[!grepl("_system_summary\\.csv$", csv_candidates)]
  csv_path <- if (length(csv_candidates)) csv_candidates[[1]] else NA_character_
  gff_candidates <- list.files(output_dir, pattern = "_padloc\\.gff$", full.names = TRUE)
  domtbl_candidates <- list.files(output_dir, pattern = "\\.domtblout$", full.names = TRUE)
  summary_candidates <- list.files(output_dir, pattern = "_padloc_system_summary\\.csv$", full.names = TRUE)

  raw_hits <- if (!is.na(csv_path) && file.exists(csv_path)) dnmb_padloc_parse_hits(csv_path, id_map = input$map) else data.frame()
  hits <- dnmb_padloc_normalize_hits(raw_hits)
  output_table <- .dnmb_padloc_output_table(genes = genes, hits = hits)

  ok <- !is.na(csv_path) && file.exists(csv_path) && isTRUE(run$ok)
  status <- dplyr::bind_rows(
    status,
    .dnmb_padloc_status_row(
      "padloc_run",
      if (ok) "ok" else "failed",
      if (ok) csv_path else (run$error %||% output_dir)
    )
  )
  status <- status[!(is.na(status$component) & is.na(status$status) & is.na(status$detail)), , drop = FALSE]

  list(
    ok = ok,
    status = status,
    files = c(
      module$files,
      list(
        faa = input$faa,
        gff = input$gff,
        id_map = input$map_path,
        csv = csv_path,
        padloc_gff = if (length(gff_candidates)) gff_candidates[[1]] else NA_character_,
        domtblout = if (length(domtbl_candidates)) domtbl_candidates[[1]] else NA_character_,
        system_summary = if (length(summary_candidates)) summary_candidates[[1]] else NA_character_,
        trace_log = trace_log
      )
    ),
    raw_hits = raw_hits,
    hits = hits,
    output_table = output_table
  )
}
