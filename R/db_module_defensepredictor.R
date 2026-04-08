.dnmb_defensepredictor_module_name <- function() {
  "defensepredictor"
}

.dnmb_defensepredictor_default_version <- function() {
  "current"
}

.dnmb_defensepredictor_default_threshold <- function() {
  4.0
}

.dnmb_defensepredictor_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = as.character(component)[1],
    status = as.character(status)[1],
    detail = as.character(detail)[1]
  )
}

.dnmb_defensepredictor_empty_status <- function() {
  tibble::tibble(
    component = character(),
    status = character(),
    detail = character()
  )
}

.dnmb_defensepredictor_trace <- function(path, text) {
  cat(paste0(text, "\n"), file = path, append = TRUE)
}

.dnmb_defensepredictor_asset_layout <- function(module_dir) {
  list(
    module_dir = module_dir,
    env_dir = file.path(module_dir, "venv"),
    env_python = file.path(module_dir, "venv", "bin", "python"),
    env_pip = file.path(module_dir, "venv", "bin", "pip"),
    cli_path = file.path(module_dir, "bin", "defense_predictor"),
    download_cli = file.path(module_dir, "bin", "defense_predictor_download"),
    install_trace_log = file.path(module_dir, "defensepredictor_install_trace.log")
  )
}

.dnmb_defensepredictor_common_roots <- function() {
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
    path.expand("~/mambaforge")
  )
  roots <- unique(path.expand(as.character(roots)))
  roots[nzchar(roots)]
}

.dnmb_defensepredictor_detect_external <- function() {
  cli <- tryCatch(dnmb_detect_binary("defense_predictor", required = FALSE)$path, error = function(e) "")
  download_cli <- tryCatch(dnmb_detect_binary("defense_predictor_download", required = FALSE)$path, error = function(e) "")
  if (nzchar(cli) && nzchar(download_cli)) {
    return(list(cli = cli, download_cli = download_cli, detail = cli, mode = "external"))
  }

  roots <- .dnmb_defensepredictor_common_roots()
  cli_candidates <- unique(c(
    file.path(roots, "bin", "defense_predictor"),
    file.path(roots, "envs", "base", "bin", "defense_predictor")
  ))
  download_candidates <- unique(c(
    file.path(roots, "bin", "defense_predictor_download"),
    file.path(roots, "envs", "base", "bin", "defense_predictor_download")
  ))
  cli_candidates <- cli_candidates[file.exists(cli_candidates)]
  download_candidates <- download_candidates[file.exists(download_candidates)]
  if (length(cli_candidates) && length(download_candidates)) {
    return(list(
      cli = cli_candidates[[1]],
      download_cli = download_candidates[[1]],
      detail = cli_candidates[[1]],
      mode = "external"
    ))
  }

  list(cli = "", download_cli = "", detail = "", mode = "missing")
}

.dnmb_defensepredictor_prepare_env <- function(layout, force = FALSE) {
  py <- .dnmb_defensefinder_candidate_python()
  if (!nzchar(py)) {
    return(.dnmb_defensepredictor_status_row("defensepredictor_python", "missing", "No python3 executable with pip found in PATH."))
  }
  if (dir.exists(layout$env_dir) && isTRUE(force)) {
    unlink(layout$env_dir, recursive = TRUE, force = TRUE)
  }
  if (dir.exists(layout$env_dir)) {
    if (!file.exists(layout$env_python)) {
      unlink(layout$env_dir, recursive = TRUE, force = TRUE)
    } else {
      test <- dnmb_run_external(layout$env_python, args = c("-c", "print('ok')"), required = FALSE)
      if (!isTRUE(test$ok)) {
        unlink(layout$env_dir, recursive = TRUE, force = TRUE)
      }
    }
  }
  if (!file.exists(layout$env_python)) {
    run <- dnmb_run_external(py, args = c("-m", "venv", layout$env_dir), required = FALSE)
    if (!isTRUE(run$ok)) {
      return(.dnmb_defensepredictor_status_row("defensepredictor_python", "failed", run$error %||% py))
    }
  }
  .dnmb_defensepredictor_status_row("defensepredictor_python", "ok", layout$env_python)
}

.dnmb_defensepredictor_install_cli <- function(layout) {
  dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", "--upgrade", "pip"), required = FALSE)
  run <- dnmb_run_external(layout$env_pip, args = c("install", "--no-cache-dir", "defense_predictor"), required = FALSE)
  if (!isTRUE(run$ok)) {
    return(.dnmb_defensepredictor_status_row("defensepredictor_cli", "failed", run$error %||% layout$env_pip))
  }
  direct_cli <- file.path(layout$env_dir, "bin", "defense_predictor")
  direct_download <- file.path(layout$env_dir, "bin", "defense_predictor_download")
  if (!file.exists(direct_cli) || !file.exists(direct_download)) {
    return(.dnmb_defensepredictor_status_row("defensepredictor_cli", "failed", layout$env_dir))
  }
  .dnmb_write_exec_wrapper(layout$cli_path, direct_cli)
  .dnmb_write_exec_wrapper(layout$download_cli, direct_download)
  .dnmb_defensepredictor_status_row("defensepredictor_cli", "ok", direct_cli)
}

.dnmb_defensepredictor_download_weights <- function(layout) {
  run <- dnmb_run_external(layout$download_cli, required = FALSE)
  .dnmb_defensepredictor_status_row(
    "defensepredictor_weights",
    if (isTRUE(run$ok)) "ok" else "failed",
    if (isTRUE(run$ok)) layout$download_cli else (run$error %||% layout$download_cli)
  )
}

dnmb_defensepredictor_install_module <- function(version = .dnmb_defensepredictor_default_version(),
                                                 cache_root = NULL,
                                                 install = TRUE,
                                                 force = FALSE) {
  module <- .dnmb_defensepredictor_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_defensepredictor_asset_layout(module_dir)
  status <- .dnmb_defensepredictor_empty_status()

  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  if (!isTRUE(force) &&
      !is.null(manifest) &&
      isTRUE(manifest$install_ok) &&
      .dnmb_exec_wrapper_ok(layout$cli_path) &&
      .dnmb_exec_wrapper_ok(layout$download_cli)) {
    return(list(
      ok = TRUE,
      status = .dnmb_defensepredictor_status_row("defensepredictor_install", "cached", module_dir),
      files = list(cli = layout$cli_path, download_cli = layout$download_cli, trace_log = layout$install_trace_log),
      manifest = manifest
    ))
  }

  external <- .dnmb_defensepredictor_detect_external()
  if (nzchar(external$cli) && nzchar(external$download_cli)) {
    .dnmb_write_exec_wrapper(layout$cli_path, external$cli)
    .dnmb_write_exec_wrapper(layout$download_cli, external$download_cli)
    status <- dplyr::bind_rows(
      status,
      .dnmb_defensepredictor_status_row("defensepredictor_cli", "ok", external$detail)
    )
    weights_status <- .dnmb_defensepredictor_download_weights(layout)
    status <- dplyr::bind_rows(status, weights_status)
    if (weights_status$status != "ok") {
      return(list(
        ok = FALSE,
        status = status,
        files = list(cli = layout$cli_path, download_cli = layout$download_cli, trace_log = layout$install_trace_log),
        manifest = manifest
      ))
    }
    manifest <- list(
      install_ok = TRUE,
      mode = "external",
      cli_path = layout$cli_path,
      download_cli = layout$download_cli,
      external_cli = external$cli,
      external_download_cli = external$download_cli
    )
    dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
    .dnmb_db_autoprune_default_versions(
      module = module,
      version = version,
      default_version = .dnmb_defensepredictor_default_version(),
      cache_root = cache_root
    )
    return(list(
      ok = TRUE,
      status = dplyr::bind_rows(status, .dnmb_defensepredictor_status_row("defensepredictor_install", "ok", module_dir)),
      files = list(cli = layout$cli_path, download_cli = layout$download_cli, trace_log = layout$install_trace_log),
      manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
    ))
  }

  if (!isTRUE(install)) {
    return(list(
      ok = FALSE,
      status = .dnmb_defensepredictor_status_row("defensepredictor_install", "missing", "DefensePredictor is missing and module_install is FALSE."),
      files = list(),
      manifest = manifest
    ))
  }

  env_status <- .dnmb_defensepredictor_prepare_env(layout, force = force)
  status <- dplyr::bind_rows(status, env_status)
  if (env_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  cli_status <- .dnmb_defensepredictor_install_cli(layout)
  status <- dplyr::bind_rows(status, cli_status)
  if (cli_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  weights_status <- .dnmb_defensepredictor_download_weights(layout)
  status <- dplyr::bind_rows(status, weights_status)
  if (weights_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  manifest <- list(
    install_ok = TRUE,
    mode = "venv",
    env_python = layout$env_python,
    cli_path = layout$cli_path,
    download_cli = layout$download_cli
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_defensepredictor_default_version(),
    cache_root = cache_root
  )
  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_defensepredictor_status_row("defensepredictor_install", "ok", module_dir)),
    files = list(cli = layout$cli_path, download_cli = layout$download_cli, trace_log = layout$install_trace_log),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_defensepredictor_get_module <- function(version = .dnmb_defensepredictor_default_version(),
                                             cache_root = NULL,
                                             required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_defensepredictor_module_name(), version, cache_root = cache_root, required = FALSE)
  layout <- .dnmb_defensepredictor_asset_layout(.dnmb_db_module_dir(.dnmb_defensepredictor_module_name(), version, cache_root = cache_root, create = FALSE))
  ok <- !is.null(manifest) &&
    isTRUE(manifest$install_ok) &&
    file.exists(layout$cli_path) &&
    file.exists(layout$download_cli)
  if (required && !ok) {
    stop("DefensePredictor module is not installed for version `", version, "`.", call. = FALSE)
  }
  list(
    ok = ok,
    manifest = manifest,
    files = list(
      cli = layout$cli_path,
      download_cli = layout$download_cli,
      trace_log = layout$install_trace_log
    )
  )
}

dnmb_defensepredictor_parse_output <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.dnmb_defensepredictor_score_band <- function(score) {
  score <- suppressWarnings(as.numeric(score))
  out <- ifelse(
    is.na(score),
    NA_character_,
    ifelse(score >= 10, "DP_10plus",
      ifelse(score >= 8, "DP_8to10",
        ifelse(score >= 6, "DP_6to8",
          ifelse(score >= 4, "DP_4to6", "DP_lt4")
        )
      )
    )
  )
  out
}

dnmb_defensepredictor_normalize_hits <- function(hits_tbl) {
  if (is.null(hits_tbl) || !is.data.frame(hits_tbl) || !nrow(hits_tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  tbl <- as.data.frame(hits_tbl, stringsAsFactors = FALSE)
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$locus_tag %||% tbl$product_accession)
  tbl$mean_log_odds <- suppressWarnings(as.numeric(tbl$mean_log_odds))
  tbl$sd_log_odds <- suppressWarnings(as.numeric(tbl$sd_log_odds))
  tbl$min_log_odds <- suppressWarnings(as.numeric(tbl$min_log_odds))
  tbl$max_log_odds <- suppressWarnings(as.numeric(tbl$max_log_odds))
  tbl$score_band <- .dnmb_defensepredictor_score_band(tbl$mean_log_odds)

  support <- vapply(seq_len(nrow(tbl)), function(i) {
    parts <- c(
      if (!is.na(tbl$mean_log_odds[[i]])) paste0("mean=", sprintf("%.3f", tbl$mean_log_odds[[i]])) else NA_character_,
      if (!is.na(tbl$sd_log_odds[[i]])) paste0("sd=", sprintf("%.3f", tbl$sd_log_odds[[i]])) else NA_character_,
      if (!is.na(tbl$min_log_odds[[i]])) paste0("min=", sprintf("%.3f", tbl$min_log_odds[[i]])) else NA_character_,
      if (!is.na(tbl$max_log_odds[[i]])) paste0("max=", sprintf("%.3f", tbl$max_log_odds[[i]])) else NA_character_,
      if (!is.na(tbl$product_accession[[i]]) && nzchar(as.character(tbl$product_accession[[i]]))) paste0("product_accession=", tbl$product_accession[[i]]) else NA_character_,
      if (!is.na(tbl$genomic_accession[[i]]) && nzchar(as.character(tbl$genomic_accession[[i]]))) paste0("contig=", tbl$genomic_accession[[i]]) else NA_character_
    )
    parts <- parts[!is.na(parts)]
    if (!length(parts)) NA_character_ else paste(parts, collapse = "; ")
  }, character(1))

  out <- data.frame(
    query = tbl$query,
    source = "defensepredictor",
    family_system = "DefensePredictor",
    family_id = tbl$score_band,
    hit_label = tbl$name %||% tbl$product_accession,
    enzyme_role = NA_character_,
    evidence_mode = ifelse(tbl$mean_log_odds >= 8, "high_score", ifelse(tbl$mean_log_odds >= 4, "candidate", "background")),
    substrate_label = NA_character_,
    support = support,
    typing_eligible = tbl$mean_log_odds >= .dnmb_defensepredictor_default_threshold(),
    protein_context_id = as.character(tbl$protein_context_id %||% NA_character_),
    score_band = tbl$score_band,
    mean_log_odds = tbl$mean_log_odds,
    sd_log_odds = tbl$sd_log_odds,
    min_log_odds = tbl$min_log_odds,
    max_log_odds = tbl$max_log_odds,
    feature_class = as.character(tbl$class %||% NA_character_),
    assembly = as.character(tbl$assembly %||% NA_character_),
    assembly_unit = as.character(tbl$assembly_unit %||% NA_character_),
    seq_type = as.character(tbl$seq_type %||% NA_character_),
    chromosome = as.character(tbl$chromosome %||% NA_character_),
    genomic_accession = as.character(tbl$genomic_accession %||% NA_character_),
    product_accession = as.character(tbl$product_accession %||% NA_character_),
    nonredundant_refseq = as.character(tbl[["non-redundant_refseq"]] %||% NA_character_),
    related_accession = as.character(tbl$related_accession %||% NA_character_),
    name = as.character(tbl$name %||% NA_character_),
    symbol = as.character(tbl$symbol %||% NA_character_),
    GeneID = as.character(tbl$GeneID %||% NA_character_),
    feature_interval_length = suppressWarnings(as.numeric(tbl$feature_interval_length %||% NA)),
    product_length = suppressWarnings(as.numeric(tbl$product_length %||% NA)),
    attributes = as.character(tbl$attributes %||% NA_character_),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$query) & nzchar(out$query), , drop = FALSE]
  out <- out[, c(.dnmb_module_optional_long_columns(), setdiff(names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_defensepredictor_hits_for_output <- function(hits,
                                                   threshold = .dnmb_defensepredictor_default_threshold()) {
  if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out <- as.data.frame(hits, stringsAsFactors = FALSE)
  out$mean_log_odds <- suppressWarnings(as.numeric(out$mean_log_odds))
  out <- out[!is.na(out$mean_log_odds) & out$mean_log_odds >= as.numeric(threshold)[1], , drop = FALSE]
  out <- out[order(out$query, -out$mean_log_odds), , drop = FALSE]
  out <- out[!duplicated(out$query), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_defensepredictor_output_table <- function(genes,
                                                hits,
                                                threshold = .dnmb_defensepredictor_default_threshold()) {
  selected_hits <- .dnmb_defensepredictor_hits_for_output(hits, threshold = threshold)
  out <- .dnmb_module_output_table(genes = genes, hits = selected_hits)
  drop_cols <- intersect(c("family_id", "hit_label", "enzyme_role", "evidence_mode", "substrate_label", "typing_eligible"), names(out))
  if (length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  ordered <- c(
    intersect(dnmb_backbone_columns(), names(out)),
    intersect(
      c("score_band", "mean_log_odds", "sd_log_odds", "min_log_odds", "max_log_odds",
        "protein_context_id", "product_accession", "genomic_accession", "name", "symbol",
        "feature_class", "feature_interval_length", "product_length", "support"),
      names(out)
    ),
    setdiff(
      names(out),
      c(intersect(dnmb_backbone_columns(), names(out)),
        "score_band", "mean_log_odds", "sd_log_odds", "min_log_odds", "max_log_odds",
        "protein_context_id", "product_accession", "genomic_accession", "name", "symbol",
        "feature_class", "feature_interval_length", "product_length", "support")
    )
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_defensepredictor_module <- function(genes,
                                             output_dir,
                                             version = .dnmb_defensepredictor_default_version(),
                                             cache_root = NULL,
                                             install = TRUE,
                                             threshold = .dnmb_defensepredictor_default_threshold()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- file.path(output_dir, "defensepredictor_module_trace.log")
  status <- .dnmb_defensepredictor_empty_status()

  install_result <- dnmb_defensepredictor_install_module(
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

  module <- dnmb_defensepredictor_get_module(version = version, cache_root = cache_root, required = TRUE)
  input <- .dnmb_write_defensepredictor_input(genes = genes, output_dir = output_dir, prefix = "defense_predictor_query")
  status <- dplyr::bind_rows(
    status,
    .dnmb_defensepredictor_status_row("defensepredictor_input", if (nrow(input$map)) "ok" else "empty", paste0("proteins=", nrow(input$map)))
  )
  if (!nrow(input$map)) {
    return(list(
      ok = TRUE,
      status = status,
      files = c(module$files, list(feature_table = input$feature_table, cds_fna = input$cds_fna, protein_faa = input$protein_faa, id_map = input$map_path, trace_log = trace_log)),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame()
    ))
  }

  output_csv <- file.path(output_dir, "defense_predictor_output.csv")
  args <- c(
    "--ncbi_feature_table", input$feature_table,
    "--ncbi_cds_from_genomic", input$cds_fna,
    "--ncbi_protein_fasta", input$protein_faa,
    "--output", output_csv
  )
  .dnmb_defensepredictor_trace(trace_log, sprintf("[%s] %s", Sys.time(), .dnmb_format_command(module$files$cli, args)))
  run <- dnmb_run_external(module$files$cli, args = args, required = FALSE)
  if (length(run$stdout)) {
    .dnmb_defensepredictor_trace(trace_log, paste(run$stdout, collapse = "\n"))
  }
  if (length(run$stderr)) {
    .dnmb_defensepredictor_trace(trace_log, paste(run$stderr, collapse = "\n"))
  }

  raw_hits <- dnmb_defensepredictor_parse_output(output_csv)
  hits <- dnmb_defensepredictor_normalize_hits(raw_hits)
  output_table <- .dnmb_defensepredictor_output_table(genes = genes, hits = hits, threshold = threshold)
  ok <- file.exists(output_csv) && isTRUE(run$ok)
  status <- dplyr::bind_rows(
    status,
    .dnmb_defensepredictor_status_row(
      "defensepredictor_run",
      if (ok) "ok" else "failed",
      if (ok) output_csv else (run$error %||% output_dir)
    )
  )
  status <- status[!(is.na(status$component) & is.na(status$status) & is.na(status$detail)), , drop = FALSE]

  list(
    ok = ok,
    status = status,
    files = c(
      module$files,
      list(
        feature_table = input$feature_table,
        cds_fna = input$cds_fna,
        protein_faa = input$protein_faa,
        id_map = input$map_path,
        output_csv = output_csv,
        trace_log = trace_log
      )
    ),
    raw_hits = raw_hits,
    hits = hits,
    output_table = output_table
  )
}
