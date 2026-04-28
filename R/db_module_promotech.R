.dnmb_promotech_module_name <- function() {
  "promotech"
}

.dnmb_promotech_default_version <- function() {
  "current"
}

.dnmb_promotech_default_repo_url <- function() {
  "https://github.com/BioinformaticsLabAtMUN/Promotech.git"
}

.dnmb_promotech_default_model_base_url <- function() {
  "https://www.cs.mun.ca/~lourdes/public/PromoTech_models"
}

.dnmb_promotech_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_promotech_empty_status <- function() {
  tibble::tibble(
    component = character(),
    status = character(),
    detail = character()
  )
}

.dnmb_promotech_asset_layout <- function(module_dir) {
  list(
    module_dir = module_dir,
    repo_dir = base::file.path(module_dir, "Promotech"),
    script_path = base::file.path(module_dir, "Promotech", "promotech.py"),
    models_dir = base::file.path(module_dir, "Promotech", "models"),
    runner_path = base::file.path(module_dir, "dnmb_promotech_runner.py")
  )
}

.dnmb_promotech_required_repo_entries <- function() {
  c("promotech.py", "genome", "core", "models", "sequences")
}

.dnmb_promotech_copy_dir_contents <- function(src_dir, dest_dir) {
  base::dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entry_names <- .dnmb_promotech_required_repo_entries()
  entries <- base::file.path(src_dir, entry_names)
  entries <- entries[base::file.exists(entries)]
  if (!base::length(entries)) {
    base::stop("Promotech source does not contain required runtime files: ", src_dir, call. = FALSE)
  }
  ok <- base::file.copy(entries, dest_dir, recursive = TRUE, copy.mode = TRUE)
  if (!base::all(ok)) {
    base::stop("Failed to copy Promotech assets from ", src_dir, " to ", dest_dir, call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_promotech_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for Promotech must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_promotech_model_asset_name <- function(model = "RF-HOT") {
  model <- base::as.character(model)[1]
  if (model %in% c("RF-HOT", "RF-TETRA")) {
    return(base::paste0(model, ".model"))
  }
  if (identical(model, "GRU")) {
    return("GRU-0.h5")
  }
  if (identical(model, "LSTM")) {
    return("LSTM-3.h5")
  }
  base::paste0(model, ".model")
}

.dnmb_promotech_default_model_source <- function(model = "RF-HOT",
                                                model_base_url = .dnmb_promotech_default_model_base_url()) {
  model <- base::as.character(model)[1]
  model_base_url <- sub("/+$", "", base::as.character(model_base_url)[1])
  if (model %in% c("RF-HOT", "RF-TETRA")) {
    return(base::paste0(model_base_url, "/", model, ".zip"))
  }
  base::paste0(model_base_url, "/", .dnmb_promotech_model_asset_name(model))
}

.dnmb_promotech_extract_or_copy_model <- function(source, target, force = FALSE) {
  if (base::file.exists(target) && !base::isTRUE(force)) {
    return(list(ok = TRUE, detail = target, status = "cached"))
  }

  source <- base::as.character(source)[1]
  if (base::is.na(source) || !base::nzchar(source)) {
    return(list(ok = FALSE, detail = "Promotech model source is empty.", status = "failed"))
  }

  base::dir.create(base::dirname(target), recursive = TRUE, showWarnings = FALSE)
  is_zip <- grepl("\\.zip($|[?#])", source, ignore.case = TRUE)
  download_dest <- if (is_zip) {
    base::file.path(base::dirname(target), base::paste0(tools::file_path_sans_ext(base::basename(target)), ".zip"))
  } else {
    target
  }

  dl <- .dnmb_download_asset(source, download_dest, insecure = FALSE)
  if (!base::isTRUE(dl$ok) || !base::file.exists(download_dest)) {
    return(list(ok = FALSE, detail = dl$error %||% source, status = "failed"))
  }

  if (!is_zip) {
    return(list(ok = TRUE, detail = target, status = "ok"))
  }

  stage_dir <- base::tempfile("dnmb-promotech-model-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  unzip_ok <- tryCatch({
    utils::unzip(download_dest, exdir = stage_dir, overwrite = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!base::isTRUE(unzip_ok)) {
    return(list(ok = FALSE, detail = base::paste0("Failed to unzip Promotech model archive: ", download_dest), status = "failed"))
  }

  target_name <- base::basename(target)
  extracted <- base::list.files(stage_dir, pattern = base::paste0("^", gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", target_name), "$"), recursive = TRUE, full.names = TRUE)
  if (!base::length(extracted)) {
    return(list(ok = FALSE, detail = base::paste0("Model archive did not contain ", target_name), status = "failed"))
  }
  ok <- base::file.copy(extracted[[1]], target, overwrite = TRUE, copy.mode = TRUE)
  list(
    ok = base::isTRUE(ok) && base::file.exists(target),
    detail = if (base::isTRUE(ok) && base::file.exists(target)) target else extracted[[1]],
    status = if (base::isTRUE(ok) && base::file.exists(target)) "ok" else "failed"
  )
}

.dnmb_promotech_prepare_model_assets <- function(layout,
                                                asset_urls = list(),
                                                force = FALSE,
                                                model = NULL,
                                                download_model = FALSE,
                                                model_base_url = .dnmb_promotech_default_model_base_url()) {
  rows <- list()
  base::dir.create(layout$models_dir, recursive = TRUE, showWarnings = FALSE)

  if (!base::is.null(asset_urls$models_dir) && base::dir.exists(base::path.expand(asset_urls$models_dir))) {
    source_dir <- base::normalizePath(base::path.expand(asset_urls$models_dir), winslash = "/", mustWork = TRUE)
    target_dir <- base::normalizePath(layout$models_dir, winslash = "/", mustWork = FALSE)
    if (!identical(source_dir, target_dir)) {
      if (base::isTRUE(force) && base::dir.exists(target_dir)) {
        base::unlink(target_dir, recursive = TRUE, force = TRUE)
        base::dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      }
      entries <- base::list.files(source_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
      ok <- if (base::length(entries)) base::file.copy(entries, target_dir, recursive = TRUE, copy.mode = TRUE) else logical()
      if (base::length(entries) && !base::all(ok)) {
        rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", "failed", base::paste0("Failed to copy models from ", source_dir))
      } else {
        rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", "ok", target_dir)
      }
    } else {
      rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", "cached", target_dir)
    }
  }

  single_models <- c(
    rf_hot_model = "RF-HOT.model",
    rf_tetra_model = "RF-TETRA.model",
    rf_hot_zip = "RF-HOT.model",
    rf_tetra_zip = "RF-TETRA.model",
    gru_model = "GRU-0.h5",
    lstm_model = "LSTM-3.h5"
  )
  for (asset_name in names(single_models)) {
    source <- asset_urls[[asset_name]]
    if (base::is.null(source) || !base::nzchar(base::as.character(source)[1])) {
      next
    }
    target <- base::file.path(layout$models_dir, single_models[[asset_name]])
    result <- .dnmb_promotech_extract_or_copy_model(source, target = target, force = force)
    rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", result$status, result$detail)
  }

  generic_model <- asset_urls$model_file
  if (!base::is.null(generic_model) && base::nzchar(base::as.character(generic_model)[1])) {
    source <- base::as.character(generic_model)[1]
    target <- base::file.path(layout$models_dir, base::basename(source))
    result <- .dnmb_promotech_extract_or_copy_model(source, target = target, force = force)
    rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", result$status, result$detail)
  }

  model_url <- asset_urls$model_url
  if (!base::is.null(model_url) && base::nzchar(base::as.character(model_url)[1]) && !base::is.null(model)) {
    target <- .dnmb_promotech_model_file(layout$models_dir, model = model)
    result <- .dnmb_promotech_extract_or_copy_model(model_url, target = target, force = force)
    rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", result$status, result$detail)
  }

  download_override <- asset_urls$download_model %||% asset_urls$download_models
  if (!base::is.null(download_override)) {
    download_model <- isTRUE(download_override) || identical(tolower(base::as.character(download_override)[1]), "true")
  }
  model_base_url <- asset_urls$model_base_url %||% asset_urls$models_base_url %||% model_base_url
  if (base::isTRUE(download_model) && !base::is.null(model)) {
    target <- .dnmb_promotech_model_file(layout$models_dir, model = model)
    source <- .dnmb_promotech_default_model_source(model = model, model_base_url = model_base_url)
    result <- .dnmb_promotech_extract_or_copy_model(source, target = target, force = force)
    rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", result$status, result$detail)
    if (model %in% c("GRU", "LSTM")) {
      tokenizer_target <- base::file.path(layout$models_dir, "tokenizer.data")
      tokenizer_source <- base::paste0(sub("/+$", "", base::as.character(model_base_url)[1]), "/tokenizer.data")
      tokenizer_result <- .dnmb_promotech_extract_or_copy_model(tokenizer_source, target = tokenizer_target, force = force)
      rows[[base::length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_models", tokenizer_result$status, tokenizer_result$detail)
    }
  }

  if (!base::length(rows)) {
    return(.dnmb_promotech_empty_status())
  }
  dplyr::bind_rows(rows)
}

.dnmb_promotech_has_model_assets <- function(asset_urls = list()) {
  if (!base::is.list(asset_urls) || !base::length(asset_urls)) {
    return(FALSE)
  }
  base::any(base::names(asset_urls) %in% c("models_dir", "model_file", "model_url", "model_base_url", "models_base_url", "download_model", "download_models", "rf_hot_model", "rf_tetra_model", "rf_hot_zip", "rf_tetra_zip", "gru_model", "lstm_model"))
}

.dnmb_promotech_model_file <- function(models_dir, model = "RF-HOT") {
  base::file.path(models_dir, .dnmb_promotech_model_asset_name(model))
}

.dnmb_promotech_prepare_repo <- function(layout,
                                         repo_source = .dnmb_promotech_default_repo_url(),
                                         force = FALSE) {
  if (base::dir.exists(layout$repo_dir) &&
      base::file.exists(layout$script_path) &&
      !base::isTRUE(force)) {
    return(.dnmb_promotech_status_row("promotech_repo", "cached", layout$repo_dir))
  }

  if (base::dir.exists(layout$repo_dir)) {
    base::unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }
  base::dir.create(base::dirname(layout$repo_dir), recursive = TRUE, showWarnings = FALSE)

  repo_source <- base::path.expand(base::as.character(repo_source)[1])
  if (base::dir.exists(repo_source)) {
    .dnmb_promotech_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_promotech_status_row("promotech_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", "--filter=blob:none", "--sparse", repo_source, layout$repo_dir),
    required = FALSE
  )
  if (base::isTRUE(run$ok)) {
    sparse_run <- dnmb_run_external(
      "git",
      args = c("-C", layout$repo_dir, "sparse-checkout", "set", .dnmb_promotech_required_repo_entries()),
      required = FALSE
    )
    if (!base::isTRUE(sparse_run$ok)) {
      base::unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
      run <- sparse_run
    }
  }
  if (!base::isTRUE(run$ok)) {
    run <- dnmb_run_external(
      "git",
      args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
      required = FALSE
    )
  }
  .dnmb_promotech_status_row(
    "promotech_repo",
    if (base::isTRUE(run$ok) && base::file.exists(layout$script_path)) "ok" else "failed",
    if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% repo_source)
  )
}

.dnmb_promotech_cached_install_state <- function(layout, manifest = NULL) {
  if (base::is.null(manifest) || !base::isTRUE(manifest$install_ok)) {
    return(list(ok = FALSE, detail = "Promotech manifest is missing or marked incomplete."))
  }
  if (!base::dir.exists(layout$repo_dir)) {
    return(list(ok = FALSE, detail = base::paste0("Promotech repo is missing: ", layout$repo_dir)))
  }
  if (!base::file.exists(layout$script_path)) {
    return(list(ok = FALSE, detail = base::paste0("Promotech script is missing: ", layout$script_path)))
  }
  list(ok = TRUE, detail = layout$repo_dir)
}

dnmb_promotech_install_module <- function(version = .dnmb_promotech_default_version(),
                                          cache_root = NULL,
                                          install = TRUE,
                                          repo_url = .dnmb_promotech_default_repo_url(),
                                          asset_urls = NULL,
                                          model = "RF-HOT",
                                          download_model = FALSE,
                                          model_base_url = .dnmb_promotech_default_model_base_url(),
                                          force = FALSE) {
  module <- .dnmb_promotech_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_promotech_asset_layout(module_dir)
  asset_urls <- .dnmb_promotech_normalize_asset_urls(asset_urls)
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  cache_state <- .dnmb_promotech_cached_install_state(layout, manifest = manifest)

  if (base::isTRUE(cache_state$ok) && !base::isTRUE(force) && !base::isTRUE(download_model) && !base::isTRUE(.dnmb_promotech_has_model_assets(asset_urls))) {
    return(list(
      ok = TRUE,
      status = .dnmb_promotech_status_row("promotech_install", "cached", module_dir),
      files = list(repo_dir = layout$repo_dir, script_path = layout$script_path, models_dir = layout$models_dir),
      manifest = manifest
    ))
  }
  if (base::isTRUE(cache_state$ok) && !base::isTRUE(force)) {
    model_status <- .dnmb_promotech_prepare_model_assets(
      layout,
      asset_urls = asset_urls,
      force = FALSE,
      model = model,
      download_model = download_model,
      model_base_url = model_base_url
    )
    if (base::any(model_status$status %in% "failed", na.rm = TRUE)) {
      return(list(
        ok = FALSE,
        status = dplyr::bind_rows(.dnmb_promotech_status_row("promotech_install", "cached", module_dir), model_status),
        files = list(repo_dir = layout$repo_dir, script_path = layout$script_path, models_dir = layout$models_dir),
        manifest = manifest
      ))
    }
    return(list(
      ok = TRUE,
      status = dplyr::bind_rows(.dnmb_promotech_status_row("promotech_install", "cached", module_dir), model_status),
      files = list(repo_dir = layout$repo_dir, script_path = layout$script_path, models_dir = layout$models_dir),
      manifest = manifest
    ))
  }

  if (!base::isTRUE(install)) {
    return(list(
      ok = FALSE,
      status = .dnmb_promotech_status_row("promotech_install", "missing", cache_state$detail),
      files = list(),
      manifest = NULL
    ))
  }

  repo_source <- asset_urls$repo_dir %||% asset_urls$repo_url %||% repo_url
  repo_status <- .dnmb_promotech_prepare_repo(layout, repo_source = repo_source, force = force)
  if (!repo_status$status %in% c("ok", "cached")) {
    return(list(ok = FALSE, status = repo_status, files = list(), manifest = NULL))
  }
  model_status <- .dnmb_promotech_prepare_model_assets(
    layout,
    asset_urls = asset_urls,
    force = force,
    model = model,
    download_model = download_model,
    model_base_url = model_base_url
  )
  if (base::any(model_status$status %in% "failed", na.rm = TRUE)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(repo_status, model_status),
      files = list(repo_dir = layout$repo_dir, script_path = layout$script_path, models_dir = layout$models_dir),
      manifest = NULL
    ))
  }

  manifest <- list(
    install_ok = TRUE,
    repo_dir = layout$repo_dir,
    script_path = layout$script_path,
    models_dir = layout$models_dir,
    repo_url = repo_url,
    model = model,
    model_base_url = model_base_url
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_promotech_default_version(),
    cache_root = cache_root
  )
  list(
    ok = TRUE,
    status = dplyr::bind_rows(repo_status, model_status, .dnmb_promotech_status_row("promotech_install", "ok", module_dir)),
    files = list(repo_dir = layout$repo_dir, script_path = layout$script_path, models_dir = layout$models_dir),
    manifest = manifest
  )
}

dnmb_promotech_get_module <- function(version = .dnmb_promotech_default_version(),
                                      cache_root = NULL,
                                      required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_promotech_module_name(), version, cache_root = cache_root, required = FALSE)
  if (base::is.null(manifest)) {
    if (base::isTRUE(required)) {
      base::stop("Promotech module is not installed.", call. = FALSE)
    }
    return(NULL)
  }
  files <- list(
    repo_dir = manifest$repo_dir,
    script_path = manifest$script_path,
    models_dir = manifest$models_dir
  )
  if (base::isTRUE(required) && (!base::dir.exists(files$repo_dir) || !base::file.exists(files$script_path))) {
    base::stop("Promotech module cache is incomplete.", call. = FALSE)
  }
  list(manifest = manifest, files = files)
}

.dnmb_promotech_prediction_candidates <- function(output_dir, asset_urls = list(), predictions = NULL) {
  candidates <- c(
    predictions,
    asset_urls$predictions,
    asset_urls$predictions_path,
    asset_urls$prediction_csv,
    base::file.path(output_dir, "genome_predictions.csv"),
    base::file.path(output_dir, "promotech_predictions.tsv"),
    base::file.path(getwd(), "genome_predictions.csv"),
    base::file.path(getwd(), "promotech_predictions.tsv")
  )
  candidates <- base::path.expand(base::as.character(candidates))
  candidates <- candidates[!base::is.na(candidates) & base::nzchar(candidates)]
  candidates
}

.dnmb_promotech_first_existing_prediction <- function(output_dir, asset_urls = list(), predictions = NULL) {
  candidates <- .dnmb_promotech_prediction_candidates(output_dir, asset_urls = asset_urls, predictions = predictions)
  existing <- candidates[base::file.exists(candidates)]
  if (base::length(existing)) {
    return(existing[[1]])
  }
  ""
}

.dnmb_promotech_read_predictions <- function(path, threshold = NULL) {
  if (base::is.null(path) || !base::nzchar(path) || !base::file.exists(path)) {
    return(data.frame())
  }
  tbl <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  if (!base::nrow(tbl)) {
    return(data.frame())
  }
  names(tbl) <- base::tolower(base::trimws(names(tbl)))
  required <- c("chrom", "start", "end", "score", "strand", "sequence")
  missing <- base::setdiff(required, names(tbl))
  if (base::length(missing)) {
    base::stop("Promotech predictions are missing columns: ", base::paste(missing, collapse = ", "), call. = FALSE)
  }
  tbl$start <- suppressWarnings(base::as.integer(tbl$start))
  tbl$end <- suppressWarnings(base::as.integer(tbl$end))
  tbl$score <- suppressWarnings(base::as.numeric(tbl$score))
  tbl$strand <- base::as.character(tbl$strand)
  tbl$chrom <- base::as.character(tbl$chrom)
  tbl$sequence <- base::as.character(tbl$sequence)
  tbl <- tbl[!base::is.na(tbl$start) & !base::is.na(tbl$end) & !base::is.na(tbl$score), , drop = FALSE]
  if (!base::is.null(threshold) && base::length(threshold) && !base::is.na(threshold[[1]])) {
    tbl <- tbl[tbl$score >= base::as.numeric(threshold)[1], , drop = FALSE]
  }
  tbl$promoter_start <- tbl$start + 1L
  tbl$promoter_end <- tbl$end + 1L
  tbl$promoter_id <- base::sprintf("Promotech_%06d", base::seq_len(base::nrow(tbl)))
  tbl
}

.dnmb_promotech_contig_keys <- function(x) {
  x <- base::trimws(base::as.character(x))
  x <- x[!base::is.na(x) & base::nzchar(x)]
  if (!base::length(x)) {
    return(character())
  }
  base::unique(c(x, base::sub("\\.\\d+$", "", x)))
}

.dnmb_promotech_escape_qualifier <- function(x) {
  x <- base::as.character(x)[1] %||% ""
  x <- base::gsub("[\r\n\t]+", " ", x)
  x <- base::gsub("\\\\", "\\\\\\\\", x)
  base::gsub("\"", "'", x)
}

.dnmb_promotech_feature_location <- function(start, end, strand) {
  start <- suppressWarnings(base::as.integer(start)[1])
  end <- suppressWarnings(base::as.integer(end)[1])
  if (base::is.na(start) || base::is.na(end)) {
    return(NA_character_)
  }
  if (start > end) {
    tmp <- start
    start <- end
    end <- tmp
  }
  interval <- base::paste0(start, "..", end)
  if (identical(base::as.character(strand)[1], "-")) {
    return(base::paste0("complement(", interval, ")"))
  }
  interval
}

.dnmb_promotech_feature_inputs <- function(predictions, hits = NULL) {
  if (!base::is.data.frame(predictions) || !base::nrow(predictions)) {
    return(data.frame())
  }
  pred <- base::as.data.frame(predictions, stringsAsFactors = FALSE)
  if (!"promoter_start" %in% base::names(pred) && "start" %in% base::names(pred)) {
    pred$promoter_start <- suppressWarnings(base::as.integer(pred$start)) + 1L
  }
  if (!"promoter_end" %in% base::names(pred) && "end" %in% base::names(pred)) {
    pred$promoter_end <- suppressWarnings(base::as.integer(pred$end)) + 1L
  }
  if (!"promoter_id" %in% base::names(pred)) {
    pred$promoter_id <- base::sprintf("Promotech_%06d", base::seq_len(base::nrow(pred)))
  }
  for (column in c("chrom", "strand", "sequence")) {
    if (!column %in% base::names(pred)) {
      pred[[column]] <- NA_character_
    }
  }
  if (!"score" %in% base::names(pred)) {
    pred$score <- NA_real_
  }
  pred$promoter_start <- suppressWarnings(base::as.integer(pred$promoter_start))
  pred$promoter_end <- suppressWarnings(base::as.integer(pred$promoter_end))
  pred$score <- suppressWarnings(base::as.numeric(pred$score))

  pred$target_gene <- NA_character_
  pred$distance_to_gene <- NA_integer_
  if (base::is.data.frame(hits) && base::nrow(hits) && "promoter_id" %in% base::names(hits)) {
    hit_tbl <- base::as.data.frame(hits, stringsAsFactors = FALSE)
    if (!"distance_to_gene" %in% base::names(hit_tbl)) {
      hit_tbl$distance_to_gene <- NA_integer_
    }
    hit_tbl$distance_to_gene <- suppressWarnings(base::as.integer(hit_tbl$distance_to_gene))
    hit_tbl <- hit_tbl[!base::is.na(hit_tbl$promoter_id), , drop = FALSE]
    if (base::nrow(hit_tbl)) {
      hit_tbl <- hit_tbl[base::order(hit_tbl$promoter_id, hit_tbl$distance_to_gene), , drop = FALSE]
      hit_tbl <- hit_tbl[!base::duplicated(hit_tbl$promoter_id), , drop = FALSE]
      matched <- base::match(pred$promoter_id, hit_tbl$promoter_id)
      query_col <- if ("query" %in% base::names(hit_tbl)) "query" else if ("locus_tag" %in% base::names(hit_tbl)) "locus_tag" else NA_character_
      if (!base::is.na(query_col)) {
        pred$target_gene <- base::as.character(hit_tbl[[query_col]][matched])
      }
      pred$distance_to_gene <- hit_tbl$distance_to_gene[matched]
    }
  }

  pred <- pred[!base::is.na(pred$promoter_start) & !base::is.na(pred$promoter_end), , drop = FALSE]
  base::rownames(pred) <- NULL
  pred
}

.dnmb_promotech_feature_lines <- function(predictions,
                                          hits = NULL,
                                          feature_key = "promoter") {
  pred <- .dnmb_promotech_feature_inputs(predictions, hits = hits)
  if (!base::nrow(pred)) {
    return(character())
  }
  loc <- mapply(
    .dnmb_promotech_feature_location,
    pred$promoter_start,
    pred$promoter_end,
    pred$strand,
    USE.NAMES = FALSE
  )
  keep <- !base::is.na(loc)
  if (!base::any(keep)) {
    return(character())
  }
  pred <- pred[keep, , drop = FALSE]
  loc <- loc[keep]
  score <- suppressWarnings(base::as.numeric(pred$score))
  score_text <- ifelse(
    base::is.na(score),
    "NA",
    base::format(base::round(score, 5), scientific = FALSE, trim = TRUE)
  )
  direction <- ifelse(base::as.character(pred$strand) == "-", "LEFT", "RIGHT")
  notes <- vapply(base::seq_len(base::nrow(pred)), function(idx) {
    note_parts <- c(
      "Promotech promoter",
      base::paste0("score=", score_text[[idx]]),
      if (!base::is.na(pred$chrom[[idx]]) && base::nzchar(base::as.character(pred$chrom[[idx]]))) base::paste0("contig=", pred$chrom[[idx]]) else NA_character_,
      if (!base::is.na(pred$target_gene[[idx]]) && base::nzchar(base::as.character(pred$target_gene[[idx]]))) base::paste0("target=", pred$target_gene[[idx]]) else NA_character_,
      if (!base::is.na(pred$distance_to_gene[[idx]])) base::paste0("distance_to_gene=", pred$distance_to_gene[[idx]]) else NA_character_
    )
    note_parts <- note_parts[!base::is.na(note_parts) & base::nzchar(note_parts)]
    .dnmb_promotech_escape_qualifier(base::paste(note_parts, collapse = "; "))
  }, character(1))
  lines <- base::matrix("", nrow = base::nrow(pred), ncol = 4L)
  lines[, 1L] <- base::sprintf("     %-15s%s", feature_key, loc)
  label_text <- ifelse(
    base::is.na(score),
    base::as.character(pred$promoter_id),
    base::paste0(base::as.character(pred$promoter_id), " (score=", score_text, ")")
  )
  lines[, 2L] <- base::paste0("                     /label=\"", vapply(label_text, .dnmb_promotech_escape_qualifier, character(1)), "\"")
  lines[, 3L] <- base::paste0("                     /note=\"", notes, "\"")
  lines[, 4L] <- base::paste0("                     /note=\"color: #ffcc99; direction: ", direction, "\"")
  base::as.vector(base::t(lines))
}

.dnmb_promotech_genbank_records <- function(lines) {
  loci <- base::grep("^LOCUS\\s+", lines)
  if (!base::length(loci)) {
    return(data.frame())
  }
  ends <- c(loci[-1L] - 1L, base::length(lines))
  records <- lapply(base::seq_along(loci), function(i) {
    rec <- lines[loci[[i]]:ends[[i]]]
    locus <- base::sub("^LOCUS\\s+([^ ]+).*$", "\\1", rec[[1]])
    acc_idx <- base::grep("^ACCESSION\\s+", rec)
    accession <- if (base::length(acc_idx)) {
      base::strsplit(base::trimws(base::sub("^ACCESSION\\s+", "", rec[[acc_idx[[1]]]])), "\\s+")[[1]][[1]]
    } else {
      locus
    }
    data.frame(
      record_index = i,
      start_line = loci[[i]],
      end_line = ends[[i]],
      locus = locus,
      accession = accession,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(records)
}

.dnmb_promotech_record_prediction_mask <- function(predictions, record, single_record = FALSE) {
  pred_chrom <- predictions$chrom %||% rep(NA_character_, base::nrow(predictions))
  record_keys <- .dnmb_promotech_contig_keys(c(record$locus[[1]], record$accession[[1]]))
  mask <- vapply(pred_chrom, function(chrom) {
    keys <- .dnmb_promotech_contig_keys(chrom)
    base::length(keys) && base::length(base::intersect(keys, record_keys))
  }, logical(1))
  if (!base::any(mask) && base::isTRUE(single_record)) {
    mask <- rep(TRUE, base::nrow(predictions))
  }
  mask
}

.dnmb_promotech_insert_features_into_genbank <- function(genbank,
                                                        predictions,
                                                        hits = NULL,
                                                        output_path) {
  if (base::is.null(genbank) || !base::nzchar(base::as.character(genbank)[1]) || !base::file.exists(base::as.character(genbank)[1])) {
    return(list(ok = FALSE, path = NA_character_, count = 0L, unmatched = base::nrow(predictions), detail = "GenBank file is missing."))
  }
  pred <- .dnmb_promotech_feature_inputs(predictions, hits = hits)
  if (!base::nrow(pred)) {
    return(list(ok = FALSE, path = NA_character_, count = 0L, unmatched = 0L, detail = "No promoter predictions to annotate."))
  }
  lines <- base::readLines(genbank, warn = FALSE)
  records <- .dnmb_promotech_genbank_records(lines)
  if (!base::nrow(records)) {
    return(list(ok = FALSE, path = NA_character_, count = 0L, unmatched = base::nrow(pred), detail = "No LOCUS records found in GenBank."))
  }

  out <- character()
  cursor <- 1L
  annotated <- rep(FALSE, base::nrow(pred))
  inserted_count <- 0L
  single_record <- base::nrow(records) == 1L
  for (idx in base::seq_len(base::nrow(records))) {
    rec_meta <- records[idx, , drop = FALSE]
    if (cursor < rec_meta$start_line[[1]]) {
      out <- c(out, lines[cursor:(rec_meta$start_line[[1]] - 1L)])
    }
    rec_lines <- lines[rec_meta$start_line[[1]]:rec_meta$end_line[[1]]]
    mask <- .dnmb_promotech_record_prediction_mask(pred, rec_meta, single_record = single_record)
    feature_lines <- .dnmb_promotech_feature_lines(pred[mask, , drop = FALSE], hits = hits)
    if (base::length(feature_lines)) {
      annotated[mask] <- TRUE
      inserted_count <- inserted_count + base::sum(mask)
    }
    insert_at <- base::grep("^(ORIGIN|BASE COUNT|CONTIG)", rec_lines)[1]
    if (base::is.na(insert_at)) {
      insert_at <- base::grep("^//", rec_lines)[1]
    }
    if (base::is.na(insert_at)) {
      insert_at <- base::length(rec_lines) + 1L
    }
    before <- if (insert_at > 1L) rec_lines[base::seq_len(insert_at - 1L)] else character()
    after <- if (insert_at <= base::length(rec_lines)) rec_lines[insert_at:base::length(rec_lines)] else character()
    out <- c(out, before, feature_lines, after)
    cursor <- rec_meta$end_line[[1]] + 1L
  }
  if (cursor <= base::length(lines)) {
    out <- c(out, lines[cursor:base::length(lines)])
  }

  base::dir.create(base::dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  base::writeLines(out, con = output_path)
  list(
    ok = base::file.exists(output_path),
    path = output_path,
    count = inserted_count,
    unmatched = base::sum(!annotated),
    detail = output_path
  )
}

.dnmb_promotech_write_genbank_artifacts <- function(predictions,
                                                   hits = NULL,
                                                   genbank = NULL,
                                                   output_dir,
                                                   snippet_path = NULL,
                                                   annotated_genbank_path = NULL) {
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pred <- .dnmb_promotech_feature_inputs(predictions, hits = hits)
  snippet_path <- snippet_path %||% base::file.path(output_dir, "promotech_promoter_feature_for_gb")
  annotated_genbank_path <- annotated_genbank_path %||% base::file.path(output_dir, "promotech_promoters_annotated.gbk")
  if (!base::nrow(pred)) {
    return(list(
      status = .dnmb_promotech_status_row("promotech_genbank_features", "empty", "No Promotech predictions to annotate."),
      files = list(feature_snippet = snippet_path, annotated_genbank = NA_character_)
    ))
  }

  feature_lines <- .dnmb_promotech_feature_lines(pred, hits = hits)
  base::writeLines(feature_lines, con = snippet_path)
  gb_result <- .dnmb_promotech_insert_features_into_genbank(
    genbank = genbank,
    predictions = pred,
    hits = hits,
    output_path = annotated_genbank_path
  )
  detail <- if (base::isTRUE(gb_result$ok)) {
    base::paste0(
      "snippet=", snippet_path,
      "; annotated_genbank=", gb_result$path,
      "; inserted=", gb_result$count,
      "; unmatched=", gb_result$unmatched
    )
  } else {
    base::paste0("snippet=", snippet_path, "; annotated_genbank=skipped; ", gb_result$detail)
  }
  list(
    status = .dnmb_promotech_status_row(
      "promotech_genbank_features",
      if (base::file.exists(snippet_path)) "ok" else "failed",
      detail
    ),
    files = list(
      feature_snippet = snippet_path,
      annotated_genbank = if (base::isTRUE(gb_result$ok)) gb_result$path else NA_character_
    )
  )
}

.dnmb_promotech_prepare_input_fasta <- function(genbank, output_dir) {
  if (base::is.null(genbank) || !base::nzchar(base::as.character(genbank)[1]) || !base::file.exists(base::as.character(genbank)[1])) {
    base::stop("Promotech live prediction requires a GenBank file path.", call. = FALSE)
  }
  records <- .dnmb_prophage_parse_genbank_records(genbank)
  if (!base::length(records)) {
    base::stop("No sequences could be parsed from GenBank for Promotech.", call. = FALSE)
  }
  input_dir <- base::file.path(output_dir, "subjects")
  base::dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  fasta_paths <- character()
  record_count <- if (base::is.data.frame(records)) base::nrow(records) else base::length(records)
  for (idx in base::seq_len(record_count)) {
    if (base::is.data.frame(records)) {
      seq_value <- records$sequence[[idx]]
      contig_id <- records$accession[[idx]] %||% records$locus[[idx]] %||% base::sprintf("contig_%03d", idx)
    } else {
      rec <- records[[idx]]
      seq_value <- rec$sequence
      contig_id <- rec$id %||% rec$accession %||% rec$locus %||% base::sprintf("contig_%03d", idx)
    }
    seq <- base::toupper(base::gsub("[^ACGTN]", "N", seq_value))
    if (base::nchar(seq) < 41L) {
      next
    }
    safe_id <- base::gsub("[^[:alnum:]_.-]+", "_", contig_id)
    path <- base::file.path(input_dir, base::paste0(safe_id, ".fna"))
    base::writeLines(c(base::paste0(">", contig_id), seq), con = path)
    fasta_paths <- c(fasta_paths, path)
  }
  if (!base::length(fasta_paths)) {
    base::stop("No contigs >= 41 bp are available for Promotech.", call. = FALSE)
  }
  fasta_paths
}

.dnmb_promotech_write_runner <- function(path) {
  lines <- c(
    "import argparse, os, sys",
    "parser = argparse.ArgumentParser()",
    "parser.add_argument('--repo', required=True)",
    "parser.add_argument('--fasta', required=True)",
    "parser.add_argument('--out-dir', required=True)",
    "parser.add_argument('--model', default='RF-HOT')",
    "parser.add_argument('--threshold', type=float, default=0.5)",
    "parser.add_argument('--test-samples', type=int, default=None)",
    "args = parser.parse_args()",
    "os.chdir(args.repo)",
    "sys.path.insert(0, args.repo)",
    "import pandas as pd",
    "if not hasattr(pd.DataFrame, 'append'):",
    "    def _dnmb_append(self, other=None, ignore_index=False, **kwargs):",
    "        return pd.concat([self, pd.DataFrame([other])], ignore_index=ignore_index)",
    "    pd.DataFrame.append = _dnmb_append",
    "if args.model == 'RF-HOT':",
    "    import types",
    "    skbio_stub = types.ModuleType('skbio')",
    "    class _DNMBSequence(str):",
    "        def kmer_frequencies(self, *a, **k):",
    "            raise RuntimeError('scikit-bio is required for RF-TETRA features, not RF-HOT')",
    "    skbio_stub.Sequence = _DNMBSequence",
    "    sys.modules.setdefault('skbio', skbio_stub)",
    "import genome.process_genome as promotech_genome",
    "if args.model in ('GRU', 'LSTM'):",
    "    import tensorflow as tf",
    "    promotech_genome.tf = tf",
    "os.makedirs(args.out_dir, exist_ok=True)",
    "promotech_genome.parseGenome40NTSequences(fasta_file_path=args.fasta, out_dir=args.out_dir, test_sample_size=args.test_samples, data_type=args.model)",
    "promotech_genome.predictGenomeSequences(input_dir=args.out_dir, out_dir=args.out_dir, model_type=args.model, threshold=args.threshold)"
  )
  base::writeLines(lines, con = path)
  path
}

.dnmb_promotech_python_status <- function(python = "python3",
                                          model = "RF-HOT",
                                          model_file = NULL) {
  rows <- list()
  py_path <- Sys.which(python)
  if (!base::nzchar(py_path)) {
    rows[[length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_dependency_python", "missing", base::paste0(python, " not found in PATH."))
    return(dplyr::bind_rows(rows))
  }
  rows[[length(rows) + 1L]] <- .dnmb_promotech_status_row("promotech_dependency_python", "ok", py_path)

  modules <- c("numpy", "pandas", "joblib", "Bio", "progressbar")
  if (model %in% c("RF-HOT", "RF-TETRA")) {
    modules <- c(modules, "sklearn")
  } else {
    modules <- c(modules, "tensorflow")
  }
  for (module_name in modules) {
    check <- dnmb_run_external(py_path, args = c("-c", base::sprintf("import %s", module_name)), required = FALSE)
    rows[[length(rows) + 1L]] <- .dnmb_promotech_status_row(
      base::paste0("promotech_dependency_", module_name),
      if (base::isTRUE(check$ok)) "ok" else "missing",
      if (base::isTRUE(check$ok)) module_name else base::paste0("Python module ", module_name, " is missing.")
    )
  }
  if (!base::is.null(model_file) && base::nzchar(model_file)) {
    rows[[length(rows) + 1L]] <- .dnmb_promotech_status_row(
      "promotech_model",
      if (base::file.exists(model_file)) "ok" else "missing",
      if (base::file.exists(model_file)) model_file else base::paste0("Promotech model is missing: ", model_file)
    )
  }
  dplyr::bind_rows(rows)
}

.dnmb_promotech_run_live <- function(module,
                                     genbank,
                                     output_dir,
                                     model = "RF-HOT",
                                     threshold = 0.5,
                                     test_samples = NULL,
                                     python = "python3") {
  input_fastas <- .dnmb_promotech_prepare_input_fasta(genbank, output_dir)
  runner <- .dnmb_promotech_write_runner(base::file.path(output_dir, "dnmb_promotech_runner.py"))
  py_path <- Sys.which(python)
  combined <- list()
  run_results <- list()
  for (fasta in input_fastas) {
    contig_out <- base::file.path(output_dir, "promotech_runs", tools::file_path_sans_ext(base::basename(fasta)))
    base::dir.create(contig_out, recursive = TRUE, showWarnings = FALSE)
    pred_path <- base::file.path(contig_out, "genome_predictions.csv")
    # Per-contig cache reuse: if a previous run already produced
    # genome_predictions.csv for this same FASTA (same content, same
    # model) and the CSV is newer than the FASTA, skip the python
    # invocation and reload the cached predictions. This makes
    # repeated docker runs on the same genome essentially free for
    # PromoTech and avoids regenerating the multi-GB intermediate
    # RF-HOT / SEQS .data files.
    fasta_mtime <- tryCatch(base::file.info(fasta)$mtime, error = function(e) NA)
    pred_mtime <- if (base::file.exists(pred_path)) {
      tryCatch(base::file.info(pred_path)$mtime, error = function(e) NA)
    } else {
      NA
    }
    cache_hit <- !is.na(pred_mtime) && !is.na(fasta_mtime) && pred_mtime >= fasta_mtime
    if (isTRUE(cache_hit)) {
      run_results[[length(run_results) + 1L]] <- list(ok = TRUE, cached = TRUE,
                                                       contig = base::basename(contig_out))
      pred <- .dnmb_promotech_read_predictions(pred_path, threshold = threshold)
      if (base::nrow(pred)) {
        combined[[length(combined) + 1L]] <- pred
      }
      next
    }
    args <- c(
      runner,
      "--repo", module$files$repo_dir,
      "--fasta", fasta,
      "--out-dir", contig_out,
      "--model", model,
      "--threshold", base::as.character(threshold)
    )
    if (!base::is.null(test_samples) && !base::is.na(test_samples[[1]])) {
      args <- c(args, "--test-samples", base::as.character(base::as.integer(test_samples)[1]))
    }
    run <- dnmb_run_external(py_path, args = args, required = FALSE, wd = module$files$repo_dir)
    run_results[[length(run_results) + 1L]] <- run
    if (base::file.exists(pred_path)) {
      pred <- .dnmb_promotech_read_predictions(pred_path, threshold = threshold)
      if (base::nrow(pred)) {
        combined[[length(combined) + 1L]] <- pred
      }
      # Drop the multi-GB intermediate matrices (forward/reverse
      # one-hot encoded sequences) once the predictions CSV is
      # written. They aren't needed downstream and bloat the
      # per-genome module dir by ~3.5 GB.
      stale_files <- c("RF-HOT.data", "RF-HOT-INV.data",
                       "SEQS.data", "SEQS-INV.data")
      for (sf in stale_files) {
        sp <- base::file.path(contig_out, sf)
        if (base::file.exists(sp)) {
          tryCatch(base::file.remove(sp), error = function(e) FALSE)
        }
      }
    }
  }
  predictions <- if (base::length(combined)) dplyr::bind_rows(combined) else data.frame()
  output_path <- base::file.path(output_dir, "promotech_predictions.tsv")
  utils::write.table(predictions, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  list(ok = base::all(vapply(run_results, function(x) base::isTRUE(x$ok), logical(1))),
       predictions = predictions,
       prediction_path = output_path,
       runs = run_results)
}

.dnmb_promotech_map_predictions_to_genes <- function(predictions,
                                                     genes,
                                                     max_distance = 300L) {
  if (!base::is.data.frame(predictions) || !base::nrow(predictions)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  required <- c("locus_tag", "contig", "start", "end", "direction")
  if (base::length(base::setdiff(required, names(genes)))) {
    return(.dnmb_module_empty_optional_long_table())
  }
  genes$start <- suppressWarnings(base::as.integer(genes$start))
  genes$end <- suppressWarnings(base::as.integer(genes$end))
  genes$direction <- base::as.character(genes$direction)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  genes <- genes[!base::is.na(genes$locus_tag) & !base::is.na(genes$start) & !base::is.na(genes$end), , drop = FALSE]
  if (!base::nrow(genes)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  predictions <- base::as.data.frame(predictions, stringsAsFactors = FALSE)
  predictions$promoter_start <- suppressWarnings(base::as.integer(predictions$promoter_start))
  predictions$promoter_end <- suppressWarnings(base::as.integer(predictions$promoter_end))
  predictions$strand <- base::as.character(predictions$strand)
  predictions$chrom <- base::as.character(predictions$chrom)
  predictions$promoter_id <- base::as.character(predictions$promoter_id)
  predictions$sequence <- base::as.character(predictions$sequence)
  predictions$score <- suppressWarnings(base::as.numeric(predictions$score))
  predictions <- predictions[
    !base::is.na(predictions$promoter_start) &
      !base::is.na(predictions$promoter_end) &
      !base::is.na(predictions$promoter_id),
    ,
    drop = FALSE
  ]
  if (!base::nrow(predictions)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  max_distance <- base::as.integer(max_distance)[1]
  if (base::is.na(max_distance) || max_distance < 0L) {
    max_distance <- 300L
  }

  unique_contigs <- unique(base::as.character(genes$contig))
  contig_lookup <- stats::setNames(unique_contigs, unique_contigs)
  contig_keys <- lapply(unique_contigs, .dnmb_promotech_contig_keys)
  mapped_contig <- vapply(predictions$chrom, function(chrom) {
    chrom <- base::as.character(chrom)[1]
    if (!base::is.na(chrom) && chrom %in% unique_contigs) {
      return(chrom)
    }
    keys <- .dnmb_promotech_contig_keys(chrom)
    hits <- vapply(contig_keys, function(x) base::length(base::intersect(keys, x)) > 0L, logical(1))
    if (base::sum(hits) == 1L) {
      return(unique_contigs[which(hits)[[1]]])
    }
    if (base::length(unique_contigs) == 1L) {
      return(unique_contigs[[1]])
    }
    NA_character_
  }, character(1))

  .make_rows <- function(pred_idx, gene_idx, distance) {
    pred <- predictions[pred_idx, , drop = FALSE]
    gene <- genes[gene_idx, , drop = FALSE]
    score <- base::round(base::as.numeric(pred$score), 5)
    support <- base::paste0(
      "score=", score,
      "; interval=", pred$promoter_start, "-", pred$promoter_end,
      "; strand=", pred$strand,
      "; distance_to_gene=", distance
    )
    base::data.frame(
      query = gene$locus_tag,
      source = "promotech",
      family_system = "Promotech",
      family_id = "promoter",
      hit_label = "Promotech promoter",
      enzyme_role = "transcription",
      evidence_mode = "ml_prediction",
      substrate_label = "promoter",
      support = support,
      typing_eligible = TRUE,
      promoter_id = pred$promoter_id,
      promoter_chrom = pred$chrom,
      promoter_start = pred$promoter_start,
      promoter_end = pred$promoter_end,
      promoter_strand = pred$strand,
      promoter_score = base::as.numeric(pred$score),
      distance_to_gene = base::as.integer(distance),
      promoter_sequence = pred$sequence,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  for (contig in unique_contigs) {
    contig_genes <- genes[base::as.character(genes$contig) == contig, , drop = FALSE]
    if (!base::nrow(contig_genes)) {
      next
    }
    contig_pred_idx <- which(!base::is.na(mapped_contig) & mapped_contig == contig)
    if (!base::length(contig_pred_idx)) {
      next
    }

    for (strand in c("+", "-")) {
      pred_idx <- contig_pred_idx[predictions$strand[contig_pred_idx] == strand]
      if (!base::length(pred_idx)) {
        next
      }
      candidate_genes <- contig_genes
      same_strand <- candidate_genes$direction == strand
      if (base::any(same_strand, na.rm = TRUE)) {
        candidate_genes <- candidate_genes[same_strand, , drop = FALSE]
      }
      if (!base::nrow(candidate_genes)) {
        next
      }

      if (identical(strand, "-")) {
        ord <- base::order(candidate_genes$end, candidate_genes$start)
        gene_ordered <- candidate_genes[ord, , drop = FALSE]
        idx <- base::findInterval(predictions$promoter_start[pred_idx], gene_ordered$end)
        valid <- idx >= 1L
        if (base::any(valid)) {
          pred_valid <- pred_idx[valid]
          gene_valid <- idx[valid]
          dist <- predictions$promoter_start[pred_valid] - gene_ordered$end[gene_valid]
          keep <- !base::is.na(dist) & dist >= 0L & dist <= max_distance
          if (base::any(keep)) {
            rows[[base::length(rows) + 1L]] <- .make_rows(
              pred_valid[keep],
              as.integer(base::rownames(gene_ordered)[gene_valid[keep]]),
              dist[keep]
            )
          }
        }
      } else {
        ord <- base::order(candidate_genes$start, candidate_genes$end)
        gene_ordered <- candidate_genes[ord, , drop = FALSE]
        idx <- base::findInterval(predictions$promoter_end[pred_idx] - 1L, gene_ordered$start) + 1L
        valid <- idx <= base::nrow(gene_ordered)
        if (base::any(valid)) {
          pred_valid <- pred_idx[valid]
          gene_valid <- idx[valid]
          dist <- gene_ordered$start[gene_valid] - predictions$promoter_end[pred_valid]
          keep <- !base::is.na(dist) & dist >= 0L & dist <= max_distance
          if (base::any(keep)) {
            rows[[base::length(rows) + 1L]] <- .make_rows(
              pred_valid[keep],
              as.integer(base::rownames(gene_ordered)[gene_valid[keep]]),
              dist[keep]
            )
          }
        }
      }
    }
  }
  if (!base::length(rows)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  hits <- dplyr::bind_rows(rows)
  hits <- hits[base::order(hits$query, -hits$promoter_score, hits$distance_to_gene), , drop = FALSE]
  base::rownames(hits) <- NULL
  hits
}

.dnmb_promotech_output_table <- function(genes, hits) {
  out <- .dnmb_module_output_table(genes = genes, hits = hits)
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(c("family_id", "hit_label", "promoter_id", "promoter_score", "promoter_start", "promoter_end", "promoter_strand", "distance_to_gene", "promoter_sequence", "support"), base::names(out)),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "family_id", "hit_label", "promoter_id", "promoter_score", "promoter_start", "promoter_end", "promoter_strand", "distance_to_gene", "promoter_sequence", "support"))
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_promotech_module <- function(genes,
                                      output_dir,
                                      version = .dnmb_promotech_default_version(),
                                      cache_root = NULL,
                                      install = TRUE,
                                      repo_url = .dnmb_promotech_default_repo_url(),
                                      asset_urls = NULL,
                                      genbank = NULL,
                                      predictions = NULL,
                                      model = "RF-HOT",
                                      threshold = 0.5,
                                      max_distance = 300L,
                                      test_samples = NULL,
                                      python = "python3",
                                      download_model = TRUE,
                                      model_base_url = .dnmb_promotech_default_model_base_url()) {
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_promotech_empty_status()
  trace_log <- base::file.path(output_dir, "promotech_module_trace.log")
  asset_urls <- .dnmb_promotech_normalize_asset_urls(asset_urls)

  prediction_path <- .dnmb_promotech_first_existing_prediction(output_dir, asset_urls = asset_urls, predictions = predictions)
  predictions_tbl <- data.frame()
  live_result <- NULL

  if (base::nzchar(prediction_path)) {
    predictions_tbl <- .dnmb_promotech_read_predictions(prediction_path, threshold = threshold)
    status <- dplyr::bind_rows(status, .dnmb_promotech_status_row("promotech_predictions", "imported", prediction_path))
  } else {
    install_result <- dnmb_promotech_install_module(
      version = version,
      cache_root = cache_root,
      install = install,
      repo_url = repo_url,
      asset_urls = asset_urls,
      model = model,
      download_model = download_model,
      model_base_url = model_base_url
    )
    status <- dplyr::bind_rows(status, install_result$status)
    if (!base::isTRUE(install_result$ok)) {
      return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), predictions = data.frame(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
    }
    module <- dnmb_promotech_get_module(version = version, cache_root = cache_root, required = TRUE)
    model_file <- .dnmb_promotech_model_file(module$files$models_dir, model = model)
    model_status <- .dnmb_promotech_status_row(
      "promotech_model",
      if (base::file.exists(model_file)) "ok" else "missing",
      if (base::file.exists(model_file)) model_file else base::paste0("Promotech model is missing: ", model_file)
    )
    status <- dplyr::bind_rows(status, model_status)
    if (!base::file.exists(model_file)) {
      status <- dplyr::bind_rows(status, .dnmb_promotech_status_row("promotech_run", "failed", model_status$detail[[1]]))
      return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), predictions = data.frame(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
    }
    dependency_status <- .dnmb_promotech_python_status(python = python, model = model, model_file = NULL)
    status <- dplyr::bind_rows(status, dependency_status)
    blocking <- dependency_status$status %in% c("missing", "failed")
    if (base::any(blocking)) {
      detail <- base::paste(dependency_status$detail[blocking], collapse = "; ")
      status <- dplyr::bind_rows(status, .dnmb_promotech_status_row("promotech_run", "failed", detail))
      return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), predictions = data.frame(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
    }
    live_result <- .dnmb_promotech_run_live(
      module = module,
      genbank = genbank,
      output_dir = output_dir,
      model = model,
      threshold = threshold,
      test_samples = test_samples,
      python = python
    )
    prediction_path <- live_result$prediction_path
    predictions_tbl <- live_result$predictions
    status <- dplyr::bind_rows(
      status,
      .dnmb_promotech_status_row(
        "promotech_run",
        if (base::isTRUE(live_result$ok)) "ok" else "partial",
        prediction_path
      )
    )
  }

  utils::write.table(predictions_tbl, file = base::file.path(output_dir, "promotech_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  hits <- .dnmb_promotech_map_predictions_to_genes(predictions_tbl, genes = genes, max_distance = max_distance)
  utils::write.table(hits, file = base::file.path(output_dir, "promotech_mapped_hits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  output_table <- .dnmb_promotech_output_table(genes = genes, hits = hits)
  genbank_artifacts <- .dnmb_promotech_write_genbank_artifacts(
    predictions = predictions_tbl,
    hits = hits,
    genbank = genbank,
    output_dir = output_dir
  )

  status <- dplyr::bind_rows(
    status,
    .dnmb_promotech_status_row(
      "promotech_mapping",
      if (base::nrow(hits)) "ok" else "empty",
      base::paste0(base::nrow(hits), " mapped promoters from ", base::nrow(predictions_tbl), " predictions")
    ),
    genbank_artifacts$status
  )

  list(
    ok = base::nrow(predictions_tbl) > 0L || (!base::is.null(live_result) && base::isTRUE(live_result$ok)),
    status = status,
    files = list(
      trace_log = trace_log,
      predictions = base::file.path(output_dir, "promotech_predictions.tsv"),
      mapped_hits = base::file.path(output_dir, "promotech_mapped_hits.tsv"),
      feature_snippet = genbank_artifacts$files$feature_snippet,
      annotated_genbank = genbank_artifacts$files$annotated_genbank
    ),
    predictions = predictions_tbl,
    hits = hits,
    output_table = output_table
  )
}
