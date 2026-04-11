.dnmb_clean_module_name <- function() {
  "clean"
}

.dnmb_clean_default_version <- function() {
  "split100"
}

.dnmb_clean_default_base_url <- function() {
  "https://github.com/tttianhao/CLEAN.git"
}

.dnmb_clean_default_esm_url <- function() {
  "https://github.com/facebookresearch/esm.git"
}

.dnmb_clean_default_pretrained_bundle_url <- function() {
  "https://drive.google.com/file/d/1kwYd4VtzYuMvJMWXy6Vks91DSUAOcKpZ/view?usp=sharing"
}

.dnmb_clean_default_pvalue_threshold <- function() {
  1e-5
}

.dnmb_clean_default_distance_threshold <- function() {
  0.5
}

.dnmb_clean_default_nk_random <- function() {
  20L
}

.dnmb_clean_expasy_url <- function(ec_number) {
  ec_number <- trimws(as.character(ec_number))
  out <- ifelse(is.na(ec_number) | !nzchar(ec_number), NA_character_, paste0("https://purl.expasy.org/enzyme/EC/", ec_number))
  out
}

.dnmb_clean_brenda_url <- function(ec_number) {
  ec_number <- trimws(as.character(ec_number))
  out <- ifelse(is.na(ec_number) | !nzchar(ec_number), NA_character_, paste0("https://www.brenda-enzymes.org/enzyme.php?ecno=", ec_number))
  out
}

.dnmb_clean_required_pretrained_files <- function(version = .dnmb_clean_default_version()) {
  version <- as.character(version)[1]
  if (identical(version, "split100")) {
    return(c("split100.pth", "100.pt", "gmm_ensumble.pkl"))
  }
  c(paste0(version, ".pth"), "gmm_ensumble.pkl")
}

.dnmb_clean_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = as.character(component)[1],
    status = as.character(status)[1],
    detail = as.character(detail)[1]
  )
}

.dnmb_clean_asset_layout <- function(module_dir) {
  repo_dir <- file.path(module_dir, "CLEAN")
  app_dir <- file.path(repo_dir, "app")
  env_dir <- file.path(module_dir, "conda_env")
  data_dir <- file.path(app_dir, "data")
  list(
    module_dir = module_dir,
    repo_dir = repo_dir,
    app_dir = app_dir,
    env_dir = env_dir,
    env_python = file.path(env_dir, "bin", "python"),
    data_dir = data_dir,
    inputs_dir = file.path(data_dir, "inputs"),
    pretrained_dir = file.path(data_dir, "pretrained"),
    results_dir = file.path(app_dir, "results", "inputs"),
    esm_repo_dir = file.path(app_dir, "esm"),
    infer_script = file.path(app_dir, "CLEAN_infer_fasta.py"),
    build_script = file.path(app_dir, "build.py"),
    manifest_path = file.path(module_dir, "manifest.rds")
  )
}

.dnmb_clean_copy_dir_contents <- function(src_dir, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(src_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  if (!length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- file.copy(entries, dest_dir, recursive = TRUE, copy.mode = TRUE)
  if (!all(ok)) {
    stop("Failed to copy CLEAN assets from ", src_dir, " to ", dest_dir, call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_clean_normalize_asset_urls <- function(asset_urls = NULL) {
  if (is.null(asset_urls)) {
    return(list())
  }
  if (!is.list(asset_urls)) {
    stop("`asset_urls` for CLEAN must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_clean_trace <- function(path, text) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(NULL))
  }
  cat(paste0(text, "\n"), file = path, append = TRUE)
  invisible(NULL)
}

.dnmb_clean_emit <- function(text, trace_log = NULL) {
  message(text)
  .dnmb_clean_trace(trace_log, sprintf("[%s] %s", Sys.time(), text))
  invisible(NULL)
}

.dnmb_clean_prepare_repo <- function(layout,
                                     repo_source = .dnmb_clean_default_base_url(),
                                     force = FALSE,
                                     trace_log = NULL) {
  repo_source <- as.character(repo_source)[1]
  if (dir.exists(layout$repo_dir) && !isTRUE(force)) {
    return(.dnmb_clean_status_row("clean_repo", "cached", layout$repo_dir))
  }

  if (dir.exists(layout$repo_dir) && isTRUE(force)) {
    unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }

  if (dir.exists(repo_source)) {
    .dnmb_clean_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_clean_status_row("clean_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
    required = FALSE
  )
  status <- if (isTRUE(run$ok)) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed"
  .dnmb_clean_status_row("clean_repo", status, if (isTRUE(run$ok)) layout$repo_dir else run$error %||% repo_source)
}

.dnmb_clean_prepare_esm_repo <- function(layout,
                                         esm_source = .dnmb_clean_default_esm_url(),
                                         force = FALSE) {
  esm_source <- as.character(esm_source)[1]
  if (dir.exists(layout$esm_repo_dir) && !isTRUE(force)) {
    return(.dnmb_clean_status_row("clean_esm_repo", "cached", layout$esm_repo_dir))
  }

  if (dir.exists(layout$esm_repo_dir) && isTRUE(force)) {
    unlink(layout$esm_repo_dir, recursive = TRUE, force = TRUE)
  }

  if (dir.exists(esm_source)) {
    .dnmb_clean_copy_dir_contents(esm_source, layout$esm_repo_dir)
    return(.dnmb_clean_status_row("clean_esm_repo", "ok", layout$esm_repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", esm_source, layout$esm_repo_dir),
    required = FALSE
  )
  status <- if (isTRUE(run$ok)) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed"
  .dnmb_clean_status_row("clean_esm_repo", status, if (isTRUE(run$ok)) layout$esm_repo_dir else run$error %||% esm_source)
}

.dnmb_clean_verify_python <- function(python_path) {
  if (is.null(python_path) || !nzchar(python_path) || !file.exists(python_path)) {
    return(FALSE)
  }
  run <- dnmb_run_external(python_path, args = c("-c", "import torch, pandas, sklearn, esm, gdown, CLEAN"), required = FALSE)
  isTRUE(run$ok)
}

.dnmb_clean_install_into_python <- function(layout, python_path, install = TRUE) {
  python_path <- path.expand(as.character(python_path)[1])
  if (!nzchar(python_path) || !file.exists(python_path)) {
    return(.dnmb_clean_status_row("clean_python_setup", "missing", "Python executable not found."))
  }
  if (!grepl("python", basename(python_path), ignore.case = TRUE)) {
    return(.dnmb_clean_status_row("clean_python_setup", "cached", python_path))
  }
  if (.dnmb_clean_verify_python(python_path)) {
    return(.dnmb_clean_status_row("clean_python_setup", "cached", python_path))
  }
  if (!isTRUE(install)) {
    return(.dnmb_clean_status_row("clean_python_setup", "missing", "Python environment is present but CLEAN dependencies are missing."))
  }

  install_cmds <- list(
    list(command = python_path, args = c("-m", "pip", "install", "--upgrade", "pip"), wd = NULL),
    list(command = python_path, args = c("-m", "pip", "install", "torch"), wd = NULL),
    list(command = python_path, args = c("-m", "pip", "install", "-r", "requirements.txt"), wd = layout$app_dir),
    list(command = python_path, args = c("-m", "pip", "install", "gdown"), wd = NULL),
    list(command = python_path, args = c("build.py", "install"), wd = layout$app_dir)
  )

  for (step in install_cmds) {
    run <- dnmb_run_external(step$command, args = step$args, wd = step$wd, required = FALSE)
    if (!isTRUE(run$ok)) {
      return(.dnmb_clean_status_row("clean_python_setup", "failed", run$error %||% "python environment setup failed"))
    }
  }

  .dnmb_clean_status_row(
    "clean_python_setup",
    if (.dnmb_clean_verify_python(python_path)) "ok" else "failed",
    python_path
  )
}

.dnmb_clean_candidate_python_paths <- function() {
  home <- path.expand("~")
  candidates <- c(
    file.path(home, "miniforge3", "bin", "python3.10"),
    file.path(home, "Library", "r-miniconda-arm64", "bin", "python3.10"),
    file.path(home, "conda", "bin", "python3.10"),
    file.path(home, "miniconda3", "bin", "python3.10"),
    "/opt/homebrew/bin/python3.10",
    "/usr/local/bin/python3.10",
    "/opt/homebrew/bin/python3.11",
    "/usr/local/bin/python3.11",
    "/opt/homebrew/bin/python3.12",
    "/usr/local/bin/python3.12",
    "/opt/homebrew/bin/python3.13",
    "/usr/local/bin/python3.13",
    "/opt/homebrew/bin/python3",
    "/usr/local/bin/python3",
    dnmb_detect_binary("python3.10", required = FALSE)$path,
    dnmb_detect_binary("python3.11", required = FALSE)$path,
    dnmb_detect_binary("python3.12", required = FALSE)$path,
    dnmb_detect_binary("python3.13", required = FALSE)$path,
    dnmb_detect_binary("python3", required = FALSE)$path,
    dnmb_detect_binary("python", required = FALSE)$path
  )
  candidates <- unique(candidates[nzchar(candidates)])
  candidates[file.exists(candidates)]
}

.dnmb_clean_prepare_env <- function(layout,
                                    conda_prefix = NULL,
                                    python_path = NULL,
                                    install = TRUE) {
  if (!is.null(python_path) && nzchar(as.character(python_path)[1]) && file.exists(as.character(python_path)[1])) {
    resolved <- path.expand(as.character(python_path)[1])
    return(list(
      status = .dnmb_clean_status_row("clean_python", "ok", resolved),
      python = resolved,
      managed_env = FALSE
    ))
  }

  if (!is.null(conda_prefix) && nzchar(as.character(conda_prefix)[1])) {
    env_python <- file.path(path.expand(as.character(conda_prefix)[1]), "bin", "python")
    if (file.exists(env_python)) {
      return(list(
        status = .dnmb_clean_status_row("clean_python", "ok", path.expand(env_python)),
        python = path.expand(env_python),
        managed_env = FALSE
      ))
    }
  }

  if (.dnmb_clean_verify_python(layout$env_python)) {
    return(list(
      status = .dnmb_clean_status_row("clean_python", "cached", layout$env_python),
      python = normalizePath(layout$env_python, winslash = "/", mustWork = FALSE),
      managed_env = TRUE
    ))
  }

  if (!isTRUE(install)) {
    return(list(
      status = .dnmb_clean_status_row("clean_python", "missing", "CLEAN Python environment is missing."),
      python = "",
      managed_env = TRUE
    ))
  }

  # If the managed venv already exists (e.g. seeded by the Dockerfile
  # with torch/fair-esm/pandas/... pre-installed), the only thing missing
  # is the CLEAN package itself. Try running build.py install inside the
  # existing venv before any fallback that would recreate the venv —
  # recreating wipes out the seeded dependencies and then the pip
  # install -r requirements.txt step inside a plain biotools python
  # fails at wheel build, which is what we saw in the full-pipeline run.
  if (file.exists(layout$env_python) && file.exists(layout$build_script)) {
    build_run <- dnmb_run_external(
      layout$env_python,
      args = c("build.py", "install"),
      wd = layout$app_dir,
      required = FALSE
    )
    if (isTRUE(build_run$ok) && .dnmb_clean_verify_python(layout$env_python)) {
      return(list(
        status = .dnmb_clean_status_row("clean_python", "ok", layout$env_python),
        python = normalizePath(layout$env_python, winslash = "/", mustWork = FALSE),
        managed_env = TRUE
      ))
    }
  }

  install_with_env_python <- function(env_python) {
    install_cmds <- list(
      list(command = env_python, args = c("-m", "pip", "install", "--upgrade", "pip"), wd = NULL),
      list(command = env_python, args = c("-m", "pip", "install", "torch"), wd = NULL),
      list(command = env_python, args = c("-m", "pip", "install", "-r", "requirements.txt"), wd = layout$app_dir),
      list(command = env_python, args = c("-m", "pip", "install", "gdown"), wd = NULL),
      list(command = env_python, args = c("build.py", "install"), wd = layout$app_dir)
    )

    for (step in install_cmds) {
      run <- dnmb_run_external(step$command, args = step$args, wd = step$wd, required = FALSE)
      if (!isTRUE(run$ok)) {
        return(run$error %||% "environment setup failed")
      }
    }
    NULL
  }

  dir.create(dirname(layout$env_dir), recursive = TRUE, showWarnings = FALSE)
  conda <- dnmb_detect_binary("conda", required = FALSE)
  conda_error <- NULL
  if (isTRUE(conda$found)) {
    create_run <- dnmb_run_external(
      "conda",
      args = c("create", "-p", layout$env_dir, "python=3.10", "-y"),
      required = FALSE
    )
    if (isTRUE(create_run$ok)) {
      install_error <- install_with_env_python(layout$env_python)
      if (is.null(install_error)) {
        ok <- .dnmb_clean_verify_python(layout$env_python)
        return(list(
          status = .dnmb_clean_status_row("clean_python", if (ok) "ok" else "failed", if (ok) layout$env_python else "Python env created but imports failed."),
          python = if (ok) path.expand(layout$env_python) else "",
          managed_env = TRUE
        ))
      }
      conda_error <- install_error
    } else {
      conda_error <- create_run$error %||% "conda create failed"
    }
  }

  py_candidates <- .dnmb_clean_candidate_python_paths()
  if (!length(py_candidates)) {
    return(list(
      status = .dnmb_clean_status_row("clean_python", "failed", conda_error %||% "No usable Python or conda installation found."),
      python = "",
      managed_env = TRUE
    ))
  }

  venv_error <- NULL
  for (py_fallback in py_candidates) {
    venv_run <- dnmb_run_external(py_fallback, args = c("-m", "venv", layout$env_dir), required = FALSE)
    if (isTRUE(venv_run$ok)) {
      install_error <- install_with_env_python(layout$env_python)
      if (!is.null(install_error)) {
        venv_error <- install_error
        unlink(layout$env_dir, recursive = TRUE, force = TRUE)
        next
      }

      ok <- .dnmb_clean_verify_python(layout$env_python)
      return(list(
        status = .dnmb_clean_status_row("clean_python", if (ok) "ok" else "failed", if (ok) layout$env_python else "Python env created but imports failed."),
        python = if (ok) path.expand(layout$env_python) else "",
        managed_env = TRUE
      ))
    }
    venv_error <- venv_run$error
  }
  list(
    status = .dnmb_clean_status_row("clean_python", "failed", paste(c(conda_error, venv_error), collapse = " | ")),
    python = "",
    managed_env = TRUE
  )
}

.dnmb_clean_copy_pretrained_dir <- function(src_dir, dest_dir) {
  src_dir <- normalizePath(src_dir, winslash = "/", mustWork = TRUE)
  .dnmb_clean_copy_dir_contents(src_dir, dest_dir)
}

.dnmb_clean_download_pretrained_bundle <- function(layout, python_path, bundle_url) {
  zip_path <- file.path(layout$module_dir, "clean_pretrained_bundle.zip")
  stage_dir <- file.path(layout$module_dir, "clean_pretrained_bundle")
  unlink(stage_dir, recursive = TRUE, force = TRUE)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  run <- dnmb_run_external(
    python_path,
    args = c("-m", "gdown", "--fuzzy", bundle_url, "-O", zip_path),
    required = FALSE
  )
  if (!isTRUE(run$ok) || !file.exists(zip_path)) {
    return(.dnmb_clean_status_row("clean_pretrained", "failed", run$error %||% bundle_url))
  }

  utils::unzip(zip_path, exdir = stage_dir)
  top_dirs <- list.dirs(stage_dir, recursive = FALSE, full.names = TRUE)
  source_dir <- if (length(top_dirs) == 1L) top_dirs[[1L]] else stage_dir
  .dnmb_clean_copy_dir_contents(source_dir, layout$pretrained_dir)
  .dnmb_clean_status_row("clean_pretrained", "ok", layout$pretrained_dir)
}

.dnmb_clean_prepare_pretrained <- function(layout,
                                           version = .dnmb_clean_default_version(),
                                           python_path = NULL,
                                           install = TRUE,
                                           pretrained_dir = NULL,
                                           pretrained_bundle_url = .dnmb_clean_default_pretrained_bundle_url()) {
  dir.create(layout$pretrained_dir, recursive = TRUE, showWarnings = FALSE)
  required_files <- .dnmb_clean_required_pretrained_files(version = version)
  existing <- file.exists(file.path(layout$pretrained_dir, required_files))
  if (all(existing)) {
    return(.dnmb_clean_status_row("clean_pretrained", "cached", layout$pretrained_dir))
  }

  if (!is.null(pretrained_dir) && nzchar(as.character(pretrained_dir)[1]) && dir.exists(as.character(pretrained_dir)[1])) {
    .dnmb_clean_copy_pretrained_dir(as.character(pretrained_dir)[1], layout$pretrained_dir)
    existing <- file.exists(file.path(layout$pretrained_dir, required_files))
    return(.dnmb_clean_status_row("clean_pretrained", if (all(existing)) "ok" else "failed", layout$pretrained_dir))
  }

  if (!isTRUE(install) || is.null(python_path) || !nzchar(python_path)) {
    return(.dnmb_clean_status_row(
      "clean_pretrained",
      "missing",
      paste("Missing pretrained assets:", paste(required_files[!existing], collapse = ", "))
    ))
  }

  .dnmb_clean_download_pretrained_bundle(layout, python_path = python_path, bundle_url = pretrained_bundle_url)
}

dnmb_clean_install_module <- function(version = .dnmb_clean_default_version(),
                                      cache_root = NULL,
                                      install = TRUE,
                                      base_url = .dnmb_clean_default_base_url(),
                                      asset_urls = NULL) {
  module <- .dnmb_clean_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_clean_asset_layout(module_dir)
  asset_urls <- .dnmb_clean_normalize_asset_urls(asset_urls)

  repo_source <- asset_urls$repo %||% base_url
  esm_source <- asset_urls$esm_repo %||% .dnmb_clean_default_esm_url()
  python_override <- asset_urls$python %||% NULL
  conda_prefix <- asset_urls$conda_prefix %||% NULL
  pretrained_dir <- asset_urls$pretrained_dir %||% NULL
  pretrained_bundle_url <- asset_urls$pretrained_bundle_url %||% asset_urls$pretrained_bundle %||% .dnmb_clean_default_pretrained_bundle_url()

  status <- dplyr::bind_rows(
    .dnmb_clean_prepare_repo(layout, repo_source = repo_source, force = FALSE),
    .dnmb_clean_prepare_esm_repo(layout, esm_source = esm_source, force = FALSE)
  )

  env_info <- .dnmb_clean_prepare_env(
    layout = layout,
    conda_prefix = conda_prefix,
    python_path = python_override,
    install = install
  )
  status <- dplyr::bind_rows(status, env_info$status)
  python_setup_status <- .dnmb_clean_install_into_python(layout, env_info$python, install = install)
  status <- dplyr::bind_rows(status, python_setup_status)

  pretrained_status <- .dnmb_clean_prepare_pretrained(
    layout = layout,
    version = version,
    python_path = env_info$python,
    install = install,
    pretrained_dir = pretrained_dir,
    pretrained_bundle_url = pretrained_bundle_url
  )
  status <- dplyr::bind_rows(status, pretrained_status)

  ready <- all(c(
    file.exists(layout$infer_script),
    file.exists(layout$build_script),
    file.exists(layout$esm_repo_dir),
    file.exists(file.path(layout$pretrained_dir, .dnmb_clean_required_pretrained_files(version)))
  )) && nzchar(env_info$python)

  manifest <- list(
    install_ok = ready,
    module = module,
    version = version,
    module_dir = module_dir,
    repo_dir = if (dir.exists(layout$repo_dir)) normalizePath(layout$repo_dir, winslash = "/", mustWork = FALSE) else layout$repo_dir,
    app_dir = if (dir.exists(layout$app_dir)) normalizePath(layout$app_dir, winslash = "/", mustWork = FALSE) else layout$app_dir,
    env_python = env_info$python,
    pretrained_dir = if (dir.exists(layout$pretrained_dir)) normalizePath(layout$pretrained_dir, winslash = "/", mustWork = FALSE) else layout$pretrained_dir,
    required_pretrained = .dnmb_clean_required_pretrained_files(version),
    repo_source = repo_source,
    esm_source = esm_source,
    pretrained_bundle_url = pretrained_bundle_url
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_clean_default_version(),
    cache_root = cache_root
  )

  list(
    ok = ready,
    status = status,
    files = list(
      module_dir = module_dir,
      repo_dir = layout$repo_dir,
      app_dir = layout$app_dir,
      env_python = env_info$python,
      pretrained_dir = layout$pretrained_dir
    ),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root)
  )
}

dnmb_clean_resolve_module <- function(version = .dnmb_clean_default_version(),
                                      cache_root = NULL,
                                      required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_clean_module_name(), version, cache_root = cache_root, required = required)
  if (is.null(manifest)) {
    return(NULL)
  }
  if (isTRUE(required) && !isTRUE(manifest$install_ok)) {
    stop("CLEAN module is not installed or incomplete for version `", version, "`.", call. = FALSE)
  }

  layout <- .dnmb_clean_asset_layout(manifest$module_dir)
  list(
    module_dir = manifest$module_dir,
    repo_dir = manifest$repo_dir %||% layout$repo_dir,
    app_dir = manifest$app_dir %||% layout$app_dir,
    env_python = manifest$env_python %||% layout$env_python,
    pretrained_dir = manifest$pretrained_dir %||% layout$pretrained_dir,
    results_dir = layout$results_dir,
    inputs_dir = layout$inputs_dir,
    infer_script = layout$infer_script,
    manifest = manifest
  )
}

.dnmb_clean_parse_prediction_token <- function(token) {
  token <- trimws(as.character(token)[1])
  if (is.na(token) || !nzchar(token)) {
    return(list(ec = NA_character_, distance = NA_real_, raw = NA_character_))
  }

  matched <- regexec("^EC:([^/]+)/([0-9.]+)$", token, perl = TRUE)
  pieces <- regmatches(token, matched)[[1]]
  if (length(pieces) == 3L) {
    return(list(
      ec = pieces[[2]],
      distance = suppressWarnings(as.numeric(pieces[[3]])),
      raw = token
    ))
  }

  list(ec = sub("^EC:", "", token), distance = NA_real_, raw = token)
}

dnmb_clean_parse_maxsep <- function(path) {
  if (!file.exists(path)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  rows <- lapply(lines, function(line) {
    fields <- strsplit(line, ",", fixed = TRUE)[[1]]
    query <- .dnmb_module_clean_annotation_key(fields[[1]])
    preds <- trimws(fields[-1])
    preds <- preds[nzchar(preds)]
    parsed <- lapply(preds, .dnmb_clean_parse_prediction_token)
    best <- if (length(parsed)) parsed[[1]] else list(ec = NA_character_, distance = NA_real_, raw = NA_character_)
    tibble::tibble(
      query = query,
      source = "CLEAN",
      family_system = "EC",
      family_id = best$ec,
      hit_label = best$ec,
      enzyme_role = "enzyme",
      evidence_mode = "contrastive",
      substrate_label = NA_character_,
      support = if (length(preds)) paste(preds, collapse = "; ") else NA_character_,
      typing_eligible = !is.na(best$ec) & nzchar(best$ec),
      best_hit_label = best$ec,
      best_distance = best$distance,
      prediction_count = length(preds),
      all_predictions = if (length(preds)) paste(preds, collapse = "; ") else NA_character_
    )
  })

  out <- dplyr::bind_rows(rows)
  out <- out[, c(.dnmb_module_optional_long_columns(), setdiff(names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_clean_parse_pvalue_detail_token <- function(token) {
  token <- trimws(as.character(token)[1])
  if (is.na(token) || !nzchar(token)) {
    return(list(ec = NA_character_, pvalue = NA_real_, distance = NA_real_, raw = NA_character_))
  }
  matched <- regexec("^EC:([^/]+)/p=([^/]+)/d=([^/]+)$", token, perl = TRUE)
  pieces <- regmatches(token, matched)[[1]]
  if (length(pieces) == 4L) {
    return(list(
      ec = pieces[[2]],
      pvalue = suppressWarnings(as.numeric(pieces[[3]])),
      distance = suppressWarnings(as.numeric(pieces[[4]])),
      raw = token
    ))
  }
  list(ec = NA_character_, pvalue = NA_real_, distance = NA_real_, raw = token)
}

dnmb_clean_parse_pvalue_details <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  tbl <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(tbl)) {
    return(data.frame())
  }
  required <- c("query", "best_hit_label", "best_pvalue", "best_distance", "prediction_count", "all_predictions")
  missing <- setdiff(required, names(tbl))
  if (length(missing)) {
    stop("CLEAN pvalue detail file is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$query)
  names(tbl)[match("best_hit_label", names(tbl))] <- "pvalue_best_hit_label"
  names(tbl)[match("best_pvalue", names(tbl))] <- "pvalue_best_pvalue"
  names(tbl)[match("best_distance", names(tbl))] <- "pvalue_best_distance"
  names(tbl)[match("prediction_count", names(tbl))] <- "pvalue_prediction_count"
  names(tbl)[match("all_predictions", names(tbl))] <- "pvalue_all_predictions"
  tbl
}

.dnmb_clean_write_pvalue_helper <- function(path) {
  writeLines(
    c(
      "import argparse",
      "import csv",
      "import os",
      "import numpy as np",
      "import pandas as pd",
      "import torch",
      "from CLEAN.utils import get_ec_id_dict, retrive_esm1b_embedding, model_embedding_test, seed_everything",
      "from CLEAN.model import LayerNormNet",
      "from CLEAN.distance_map import get_dist_map_test, get_random_nk_dist_map",
      "from CLEAN.evaluate import random_nk_model",
      "",
      "parser = argparse.ArgumentParser()",
      "parser.add_argument('--fasta_data', required=True)",
      "parser.add_argument('--train_data', default='split100')",
      "parser.add_argument('--p_value', type=float, default=1e-5)",
      "parser.add_argument('--nk_random', type=int, default=20)",
      "args = parser.parse_args()",
      "",
      "use_cuda = torch.cuda.is_available()",
      "device = torch.device('cuda:0' if use_cuda else 'cpu')",
      "dtype = torch.float32",
      "test_data = 'inputs/' + args.fasta_data",
      "fasta_path = './data/' + test_data + '.fasta'",
      "csv_path = './data/' + test_data + '.csv'",
      "with open(csv_path, 'w', newline='') as csvfile:",
      "    writer = csv.writer(csvfile, delimiter='\\t')",
      "    writer.writerow(['Entry', 'EC number', 'Sequence'])",
      "    with open(fasta_path, 'r') as fastafile:",
      "        for line in fastafile:",
      "            if line.startswith('>'):",
      "                writer.writerow([line.strip()[1:], ' ', ' '])",
      "ids = []",
      "with open(fasta_path, 'r') as fastafile:",
      "    for line in fastafile:",
      "        if line.startswith('>'):",
      "            ids.append(line.strip()[1:])",
      "missing = [seq_id for seq_id in ids if not os.path.exists('./data/esm_data/' + seq_id + '.pt')]",
      "if missing:",
      "    retrive_esm1b_embedding(test_data)",
      "id_ec_train, ec_id_dict_train = get_ec_id_dict('./data/' + args.train_data + '.csv')",
      "id_ec_test, _ = get_ec_id_dict('./data/' + test_data + '.csv')",
      "model = LayerNormNet(512, 128, device, dtype)",
      "checkpoint = torch.load('./data/pretrained/' + args.train_data + '.pth', map_location=device)",
      "model.load_state_dict(checkpoint)",
      "model.eval()",
      "if args.train_data == 'split70':",
      "    emb_train = torch.load('./data/pretrained/70.pt', map_location=device)",
      "elif args.train_data == 'split100':",
      "    emb_train = torch.load('./data/pretrained/100.pt', map_location=device)",
      "else:",
      "    raise RuntimeError('Unsupported CLEAN train_data for pvalue helper: ' + args.train_data)",
      "emb_test = model_embedding_test(id_ec_test, model, device, dtype)",
      "eval_dist = get_dist_map_test(emb_train, emb_test, ec_id_dict_train, id_ec_test, device, dtype)",
      "seed_everything()",
      "eval_df = pd.DataFrame.from_dict(eval_dist)",
      "rand_nk_ids, rand_nk_emb_train = random_nk_model(id_ec_train, ec_id_dict_train, emb_train, n=args.nk_random, weighted=True)",
      "random_nk_dist_map = get_random_nk_dist_map(emb_train, rand_nk_emb_train, ec_id_dict_train, rand_nk_ids, device, dtype)",
      "nk = len(random_nk_dist_map.keys())",
      "threshold = args.p_value * nk",
      "out_path = 'results/' + test_data + '_pvalue_details.tsv'",
      "with open(out_path, 'w', newline='') as handle:",
      "    writer = csv.writer(handle, delimiter='\\t')",
      "    writer.writerow(['query', 'best_hit_label', 'best_pvalue', 'best_distance', 'prediction_count', 'all_predictions'])",
      "    for query in eval_df.columns:",
      "        preds = []",
      "        smallest = eval_df[query].nsmallest(10)",
      "        limit = min(10, len(smallest))",
      "        for i in range(limit):",
      "            ec_id = smallest.index[i]",
      "            dist_i = float(smallest.iloc[i])",
      "            rand_dists = np.sort([random_nk_dist_map[rand_id][ec_id] for rand_id in random_nk_dist_map.keys()])",
      "            rank = int(np.searchsorted(rand_dists, dist_i, side='left'))",
      "            empirical_p = float(rank) / float(nk) if nk else 1.0",
      "            if rank <= threshold or i == 0:",
      "                preds.append((ec_id, empirical_p, dist_i))",
      "            else:",
      "                break",
      "        if preds:",
      "            best_ec, best_p, best_d = preds[0]",
      "            all_preds = '; '.join(['EC:%s/p=%s/d=%.4f' % (ec, ('%.6g' % p), d) for ec, p, d in preds])",
      "            writer.writerow([query, best_ec, ('%.6g' % best_p), ('%.4f' % best_d), len(preds), all_preds])",
      "        else:",
      "            writer.writerow([query, '', '', '', 0, ''])",
      "try:",
      "    import os",
      "    os.remove('./data/' + test_data + '.csv')",
      "except OSError:",
      "    pass"
    ),
    con = path
  )
  invisible(path)
}

.dnmb_clean_run_pvalue_details <- function(module,
                                           sample_name,
                                           p_value = .dnmb_clean_default_pvalue_threshold(),
                                           nk_random = .dnmb_clean_default_nk_random()) {
  result_path <- file.path(module$results_dir, paste0(sample_name, "_pvalue_details.tsv"))
  if (file.exists(result_path)) {
    return(list(
      ok = TRUE,
      status = .dnmb_clean_status_row("clean_pvalue", "cached", result_path),
      path = result_path
    ))
  }

  if (!grepl("python", basename(module$env_python), ignore.case = TRUE)) {
    return(list(
      ok = FALSE,
      status = .dnmb_clean_status_row("clean_pvalue", "missing", "CLEAN pvalue helper requires a Python interpreter."),
      path = result_path
    ))
  }

  script_path <- file.path(module$app_dir, "dnmb_clean_pvalue_helper.py")
  .dnmb_clean_write_pvalue_helper(script_path)
  run <- dnmb_run_external(
    module$env_python,
    args = c(
      basename(script_path),
      "--fasta_data", sample_name,
      "--train_data", .dnmb_clean_default_version(),
      "--p_value", format(as.numeric(p_value)[1], scientific = TRUE),
      "--nk_random", as.character(as.integer(nk_random)[1])
    ),
    wd = module$app_dir,
    env = c(
      PATH = paste(
        c(dirname(module$env_python), strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]),
        collapse = .Platform$path.sep
      )
    ),
    required = FALSE
  )
  list(
    ok = isTRUE(run$ok) && file.exists(result_path),
    status = .dnmb_clean_status_row(
      "clean_pvalue",
      if (isTRUE(run$ok) && file.exists(result_path)) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed",
      if (isTRUE(run$ok) && file.exists(result_path)) result_path else (run$error %||% result_path)
    ),
    path = result_path,
    command = run
  )
}

.dnmb_clean_merge_prediction_strings <- function(maxsep_predictions, pvalue_predictions) {
  if (is.null(maxsep_predictions)) maxsep_predictions <- character(0)
  if (is.null(pvalue_predictions)) pvalue_predictions <- character(0)
  maxsep_predictions <- as.character(maxsep_predictions)
  pvalue_predictions <- as.character(pvalue_predictions)
  n <- max(length(maxsep_predictions), length(pvalue_predictions))
  if (n == 0L) return(character(0))
  maxsep_predictions <- rep_len(maxsep_predictions, n)
  pvalue_predictions <- rep_len(pvalue_predictions, n)

  merge_one <- function(maxsep_str, pvalue_str) {
    parse_maxsep <- function(x) {
      x <- trimws(as.character(x))
      if (is.na(x) || !nzchar(x)) {
        return(data.frame())
      }
      toks <- strsplit(x, ";\\s*", perl = TRUE)[[1]]
      toks <- toks[nzchar(toks)]
      rows <- lapply(toks, function(tok) {
        m <- regexec("^EC:([^/]+)/([0-9.]+)$", tok, perl = TRUE)
        p <- regmatches(tok, m)[[1]]
        if (length(p) == 3L) {
          data.frame(ec = p[[2]], score = suppressWarnings(as.numeric(p[[3]])), stringsAsFactors = FALSE)
        } else {
          NULL
        }
      })
      rows <- Filter(Negate(is.null), rows)
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }
    parse_pvalue <- function(x) {
      x <- trimws(as.character(x))
      if (is.na(x) || !nzchar(x)) {
        return(data.frame())
      }
      toks <- strsplit(x, ";\\s*", perl = TRUE)[[1]]
      toks <- toks[nzchar(toks)]
      rows <- lapply(toks, function(tok) {
        m <- regexec("^EC:([^/]+)/p=([^/]+)/d=([^/]+)$", tok, perl = TRUE)
        p <- regmatches(tok, m)[[1]]
        if (length(p) == 4L) {
          data.frame(
            ec = p[[2]],
            pvalue = suppressWarnings(as.numeric(p[[3]])),
            distance = suppressWarnings(as.numeric(p[[4]])),
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      })
      rows <- Filter(Negate(is.null), rows)
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }

    mx <- parse_maxsep(maxsep_str)
    pv <- parse_pvalue(pvalue_str)
    if (!nrow(mx) && !nrow(pv)) {
      return(NA_character_)
    }
    ordered_ec <- unique(c(if (nrow(mx)) mx$ec else character(), if (nrow(pv)) pv$ec else character()))
    parts <- vapply(ordered_ec, function(ec) {
      score <- if (nrow(mx) && ec %in% mx$ec) mx$score[match(ec, mx$ec)] else NA_real_
      pval <- if (nrow(pv) && ec %in% pv$ec) pv$pvalue[match(ec, pv$ec)] else NA_real_
      dist <- if (nrow(pv) && ec %in% pv$ec) pv$distance[match(ec, pv$ec)] else NA_real_
      out <- paste0("EC:", ec)
      if (!is.na(score)) out <- paste0(out, "/score=", format(signif(score, 4), trim = TRUE, scientific = FALSE))
      if (!is.na(pval)) out <- paste0(out, "/p=", format(signif(pval, 4), trim = TRUE, scientific = TRUE))
      if (!is.na(dist)) out <- paste0(out, "/dist=", format(signif(dist, 4), trim = TRUE, scientific = FALSE))
      out
    }, character(1))
    paste(parts, collapse = "; ")
  }

  vapply(seq_len(n), function(i) merge_one(maxsep_predictions[[i]], pvalue_predictions[[i]]), character(1))
}

.dnmb_clean_output_table <- function(genes,
                                     hits,
                                     pvalue_threshold = .dnmb_clean_default_pvalue_threshold(),
                                     distance_threshold = .dnmb_clean_default_distance_threshold()) {
  out <- .dnmb_module_output_table(genes = genes, hits = hits)
  merged_preds <- .dnmb_clean_merge_prediction_strings(out$all_predictions, out$pvalue_all_predictions)
  out$predictions <- if (length(merged_preds) == nrow(out)) merged_preds else rep(NA_character_, nrow(out))

  keep_cols <- c(
    "best_hit_label",
    "best_distance",
    "pvalue_best_pvalue",
    "predictions"
  )
  for (col in keep_cols) {
    if (!col %in% names(out)) {
      out[[col]] <- NA
    }
  }

  pvals <- suppressWarnings(as.numeric(out$pvalue_best_pvalue))
  dists <- suppressWarnings(as.numeric(out$best_distance))
  drop_rows <- is.na(pvals) |
    !(pvals <= as.numeric(pvalue_threshold)[1]) |
    is.na(dists) |
    !(dists >= as.numeric(distance_threshold)[1])
  if (any(drop_rows)) {
    for (col in keep_cols) {
      if (is.numeric(out[[col]])) {
        out[[col]][drop_rows] <- NA_real_
      } else if (is.integer(out[[col]])) {
        out[[col]][drop_rows] <- NA_integer_
      } else if (is.logical(out[[col]])) {
        out[[col]][drop_rows] <- NA
      } else {
        out[[col]][drop_rows] <- NA_character_
      }
    }
  }

  out$expasy_link <- .dnmb_clean_expasy_url(out$best_hit_label)
  out$brenda_link <- .dnmb_clean_brenda_url(out$best_hit_label)

  base_cols <- intersect(dnmb_backbone_columns(), names(out))
  keep_cols <- c(
    keep_cols,
    "expasy_link",
    "brenda_link"
  )
  out[, c(base_cols, keep_cols), drop = FALSE]
}

dnmb_run_clean_module <- function(genes,
                                  output_dir,
                                  version = .dnmb_clean_default_version(),
                                  cache_root = NULL,
                                  install = TRUE,
                                  base_url = .dnmb_clean_default_base_url(),
                                  asset_urls = NULL,
                                  cpu = 1L,
                                  genbank = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- file.path(output_dir, "clean_module_trace.log")
  .dnmb_clean_trace(trace_log, sprintf("[%s] CLEAN run started", Sys.time()))
  fasta_path <- file.path(output_dir, "clean_query_proteins.faa")
  existing_faa <- dnmb_resolve_query_faa(genbank = genbank, output_dir = output_dir, fallback_filename = basename(fasta_path))
  if (!is.null(existing_faa) && .dnmb_can_reuse_query_fasta(existing_faa, genes)) {
    proteins <- .dnmb_prepare_query_proteins(genes)
    fasta <- list(path = existing_faa, n = nrow(proteins), proteins = proteins)
    fasta_path <- existing_faa
  } else {
    fasta <- .dnmb_write_query_fasta(genes, fasta_path)
  }

  status <- .dnmb_clean_status_row("clean_query_fasta", if (fasta$n) "ok" else "empty", paste0("proteins=", fasta$n))
  if (!fasta$n) {
    return(list(
      ok = TRUE,
      status = status,
      files = list(query_fasta = fasta_path, trace_log = trace_log),
      manifest = NULL,
      raw_hits = data.frame(),
      hits = .dnmb_module_empty_optional_long_table(),
      query_proteins = fasta$proteins
    ))
  }

  install_result <- dnmb_clean_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    base_url = base_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!isTRUE(install_result$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = c(install_result$files, list(query_fasta = fasta_path, trace_log = trace_log)),
      manifest = install_result$manifest,
      raw_hits = data.frame(),
      hits = .dnmb_module_empty_optional_long_table(),
      query_proteins = fasta$proteins
    ))
  }

  module <- dnmb_clean_resolve_module(version = version, cache_root = cache_root, required = TRUE)
  sample_name <- paste0("dnmb_clean_", format(Sys.time(), "%Y%m%d%H%M%S"))
  target_fasta <- file.path(module$inputs_dir, paste0(sample_name, ".fasta"))
  dir.create(module$inputs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(module$results_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(fasta_path, target_fasta, overwrite = TRUE)

  # --- Optimized ESM embedding: skip already-computed .pt files ---
  esm_data_dir <- file.path(module$app_dir, "data", "esm_data")
  dir.create(esm_data_dir, recursive = TRUE, showWarnings = FALSE)

  # Read query FASTA headers to find which proteins need embedding
  fasta_lines <- readLines(target_fasta, warn = FALSE)
  all_ids <- sub("^>", "", fasta_lines[grepl("^>", fasta_lines)])
  all_ids <- trimws(sub(" .*", "", all_ids))
  existing_pts <- tools::file_path_sans_ext(list.files(esm_data_dir, pattern = "\\.pt$"))
  new_ids <- setdiff(all_ids, existing_pts)

  n_threads <- max(1L, as.integer(cpu))
  # Scale toks_per_batch by available memory (default 8192 for ~8GB, scale up for more)
  sys_mem_gb <- tryCatch({
    mem_info <- system("free -g 2>/dev/null || sysctl -n hw.memsize 2>/dev/null", intern = TRUE)
    mem_val <- as.numeric(gsub("[^0-9]", "", mem_info[1]))
    if (mem_val > 1e9) mem_val / 1e9 else mem_val  # bytes vs GB
  }, error = function(e) 8)
  toks_batch <- as.character(max(4096L, min(32768L, as.integer(sys_mem_gb / 8 * 8192))))
  clean_env <- c(
    OMP_NUM_THREADS = as.character(n_threads),
    MKL_NUM_THREADS = as.character(n_threads),
    TORCH_NUM_THREADS = as.character(n_threads),
    TOKENIZERS_PARALLELISM = "false",
    PATH = paste(
      c(dirname(module$env_python), strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]),
      collapse = .Platform$path.sep
    )
  )

  if (length(new_ids) > 0L && length(new_ids) < length(all_ids)) {
    # Write a subset FASTA with only new proteins
    subset_fasta <- file.path(module$app_dir, "data", "inputs", paste0(sample_name, "_new.fasta"))
    fasta_idx <- which(grepl("^>", fasta_lines))
    keep <- sub("^>", "", sub(" .*", "", fasta_lines[fasta_idx])) %in% new_ids
    subset_lines <- character()
    for (i in seq_along(fasta_idx)) {
      if (keep[i]) {
        start <- fasta_idx[i]
        end <- if (i < length(fasta_idx)) fasta_idx[i + 1] - 1L else length(fasta_lines)
        subset_lines <- c(subset_lines, fasta_lines[start:end])
      }
    }
    writeLines(subset_lines, subset_fasta)
    .dnmb_clean_emit(sprintf("[CLEAN] ESM embedding: %d/%d proteins cached, computing %d new (toks_per_batch=%s, threads=%d)",
                    length(all_ids) - length(new_ids), length(all_ids), length(new_ids), toks_batch, n_threads), trace_log)
    if (length(new_ids) >= 500L) {
      .dnmb_clean_emit(sprintf("[CLEAN] Large uncached set detected (%d proteins). CPU-only ESM extraction can take hours.", length(new_ids)), trace_log)
    }

    esm_started <- Sys.time()
    esm_command <- dnmb_run_external(
      module$env_python,
      args = c(
        file.path("esm", "scripts", "extract.py"),
        "esm1b_t33_650M_UR50S",
        subset_fasta,
        esm_data_dir,
        "--include", "mean",
        "--toks_per_batch", toks_batch,
        "--nogpu"
      ),
      wd = module$app_dir,
      env = clean_env,
      required = FALSE,
      stream_stderr = TRUE
    )
    .dnmb_clean_trace(
      trace_log,
      sprintf("[%s] CLEAN ESM extract finished in %.1f min (ok=%s)", Sys.time(), as.numeric(difftime(Sys.time(), esm_started, units = "mins")), isTRUE(esm_command$ok))
    )
    unlink(subset_fasta, force = TRUE)
  } else if (length(new_ids) == 0L) {
    .dnmb_clean_emit(sprintf("[CLEAN] ESM embedding: all %d proteins cached, skipping extract.py", length(all_ids)), trace_log)
  } else {
    # All new — run on full FASTA with optimized params
    .dnmb_clean_emit(sprintf("[CLEAN] ESM embedding: computing %d proteins (toks_per_batch=%s, threads=%d)", length(all_ids), toks_batch, n_threads), trace_log)
    if (length(all_ids) >= 500L) {
      .dnmb_clean_emit(sprintf("[CLEAN] No reusable ESM cache detected for %d proteins. CPU-only ESM extraction can take hours.", length(all_ids)), trace_log)
    }
    esm_started <- Sys.time()
    esm_command <- dnmb_run_external(
      module$env_python,
      args = c(
        file.path("esm", "scripts", "extract.py"),
        "esm1b_t33_650M_UR50S",
        target_fasta,
        esm_data_dir,
        "--include", "mean",
        "--toks_per_batch", toks_batch,
        "--nogpu"
      ),
      wd = module$app_dir,
      env = clean_env,
      required = FALSE,
      stream_stderr = TRUE
    )
    .dnmb_clean_trace(
      trace_log,
      sprintf("[%s] CLEAN ESM extract finished in %.1f min (ok=%s)", Sys.time(), as.numeric(difftime(Sys.time(), esm_started, units = "mins")), isTRUE(esm_command$ok))
    )
  }

  # --- Run CLEAN inference (embeddings already computed above) ---
  # Write a wrapper that skips ESM re-extraction and goes straight to inference
  wrapper_path <- file.path(module$app_dir, paste0("_dnmb_infer_", sample_name, ".py"))
  writeLines(c(
    "import os, sys, csv",
    "from CLEAN.infer import infer_maxsep",
    "",
    paste0("test_data = 'inputs/", sample_name, "'"),
    "train_data = 'split100'",
    "",
    "# Write dummy CSV (ESM embeddings already exist)",
    "csvfile = open('./data/' + test_data + '.csv', 'w', newline='')",
    "csvwriter = csv.writer(csvfile, delimiter='\\t')",
    "csvwriter.writerow(['Entry', 'EC number', 'Sequence'])",
    "fastafile = open('./data/' + test_data + '.fasta', 'r')",
    "seq = ''",
    "name = ''",
    "for line in fastafile.readlines():",
    "    if line[0] == '>':",
    "        if name:",
    "            csvwriter.writerow([name, 'EC:0.0.0.0', seq])",
    "        name = line[1:].strip().split()[0]",
    "        seq = ''",
    "    else:",
    "        seq += line.strip()",
    "if name:",
    "    csvwriter.writerow([name, 'EC:0.0.0.0', seq])",
    "csvfile.close()",
    "fastafile.close()",
    "",
    "infer_maxsep(train_data, test_data, report_metrics=False, pretrained=True, gmm='./data/pretrained/gmm_ensumble.pkl')",
    "os.remove('./data/' + test_data + '.csv')"
  ), wrapper_path)
  on.exit(unlink(wrapper_path, force = TRUE), add = TRUE)

  .dnmb_clean_emit("[CLEAN] Starting CLEAN inference from cached embeddings", trace_log)
  infer_started <- Sys.time()
  command <- dnmb_run_external(
    module$env_python,
    args = c(basename(wrapper_path)),
    wd = module$app_dir,
    env = clean_env,
    required = FALSE
  )
  .dnmb_clean_trace(
    trace_log,
    sprintf("[%s] CLEAN inference finished in %.1f min (ok=%s)", Sys.time(), as.numeric(difftime(Sys.time(), infer_started, units = "mins")), isTRUE(command$ok))
  )

  result_csv <- file.path(module$results_dir, paste0(sample_name, "_maxsep.csv"))
  hits <- if (isTRUE(command$ok) && file.exists(result_csv)) dnmb_clean_parse_maxsep(result_csv) else .dnmb_module_empty_optional_long_table()
  status <- dplyr::bind_rows(
    status,
    .dnmb_clean_status_row(
      "clean_infer",
      if (isTRUE(command$ok)) if (nrow(hits)) "ok" else "empty" else if (!nzchar(command$resolved_command)) "missing" else "failed",
      if (isTRUE(command$ok)) result_csv else (command$error %||% "CLEAN inference failed")
    )
  )

  pvalue_result <- if (isTRUE(command$ok)) .dnmb_clean_run_pvalue_details(module, sample_name = sample_name) else list(ok = FALSE, status = .dnmb_clean_status_row("clean_pvalue", "missing", "maxsep inference failed"), path = NULL)
  status <- dplyr::bind_rows(status, pvalue_result$status)
  pvalue_hits <- if (isTRUE(pvalue_result$ok) && file.exists(pvalue_result$path)) dnmb_clean_parse_pvalue_details(pvalue_result$path) else data.frame()
  if (nrow(pvalue_hits)) {
    hits <- dplyr::left_join(hits, pvalue_hits, by = "query")
  }

  output_table <- .dnmb_clean_output_table(genes = genes, hits = hits)
  list(
    ok = isTRUE(command$ok),
    status = status,
    files = c(
      list(
        query_fasta = fasta_path,
        clean_input_fasta = target_fasta,
        clean_result_csv = result_csv,
        clean_pvalue_tsv = pvalue_result$path,
        trace_log = trace_log
      ),
      install_result$files
    ),
    manifest = install_result$manifest,
    raw_hits = hits,
    hits = hits,
    output_table = output_table,
    query_proteins = fasta$proteins
  )
}
