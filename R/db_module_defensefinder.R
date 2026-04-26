.dnmb_defensefinder_module_name <- function() {
  "defensefinder"
}

.dnmb_defensefinder_default_version <- function() {
  "current"
}

.dnmb_defensefinder_default_repo_url <- function() {
  "https://github.com/mdmparis/defense-finder.git"
}

.dnmb_defensefinder_default_models_version <- function() {
  "2.0.2"
}

.dnmb_defensefinder_default_casfinder_version <- function() {
  "3.1.0"
}

.dnmb_defensefinder_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_defensefinder_empty_status <- function() {
  .dnmb_defensefinder_status_row(character(), character(), character())
}

.dnmb_defensefinder_trace <- function(path, text) {
  base::cat(base::paste0(text, "\n"), file = path, append = TRUE)
}

.dnmb_defensefinder_trace_result <- function(path, result, label = "defensefinder_run") {
  if (base::is.null(result) || !base::is.list(result)) {
    return(invisible(FALSE))
  }
  lines <- c(
    base::sprintf("[%s] %s status=%s ok=%s", base::Sys.time(), label, result$status %||% NA_integer_, base::isTRUE(result$ok)),
    if (!base::is.null(result$error) && base::nzchar(base::as.character(result$error))) {
      base::paste0("error: ", base::as.character(result$error))
    } else {
      character()
    },
    if (base::length(result$stderr)) {
      c("stderr:", result$stderr)
    } else {
      character()
    },
    if (base::length(result$stdout)) {
      c("stdout:", result$stdout)
    } else {
      character()
    }
  )
  .dnmb_defensefinder_trace(path, base::paste(lines, collapse = "\n"))
  invisible(TRUE)
}

.dnmb_defensefinder_asset_layout <- function(module_dir) {
  list(
    module_dir = module_dir,
    repo_dir = base::file.path(module_dir, "defense-finder"),
    env_dir = base::file.path(module_dir, "venv"),
    env_python = base::file.path(module_dir, "venv", "bin", "python"),
    env_pip = base::file.path(module_dir, "venv", "bin", "pip"),
    macsydata_path = base::file.path(module_dir, "venv", "bin", "macsydata"),
    cli_path = base::file.path(module_dir, "venv", "bin", "defense-finder"),
    models_dir = base::file.path(module_dir, "models"),
    install_trace_log = base::file.path(module_dir, "defensefinder_install_trace.log")
  )
}

.dnmb_defensefinder_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for DefenseFinder must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_defensefinder_casfinder_envvar <- function() {
  "DNMB_DEFENSEFINDER_CASFINDER_DIR"
}

.dnmb_defensefinder_copy_dir_contents <- function(src_dir, dest_dir) {
  base::dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- base::list.files(src_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  if (!base::length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- base::file.copy(entries, dest_dir, recursive = TRUE, copy.mode = TRUE)
  if (!base::all(ok)) {
    base::stop("Failed to copy DefenseFinder assets from ", src_dir, " to ", dest_dir, call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_defensefinder_path_first <- function(...) {
  candidates <- unlist(list(...), use.names = FALSE)
  candidates <- base::trimws(base::as.character(candidates))
  candidates <- candidates[!base::is.na(candidates) & base::nzchar(candidates)]
  if (!base::length(candidates)) {
    return("")
  }
  base::path.expand(candidates[[1]])
}

.dnmb_defensefinder_na_if_empty <- function(x) {
  x <- base::as.character(x)
  x[!base::nzchar(base::trimws(x))] <- NA_character_
  x
}

.dnmb_defensefinder_model_parts <- function(model_fqn) {
  path_parts <- base::strsplit(base::as.character(model_fqn), "/", fixed = TRUE)
  model_root <- vapply(
    path_parts,
    function(parts) if (base::length(parts) >= 2L) parts[[2L]] else NA_character_,
    character(1)
  )
  subtype <- vapply(
    path_parts,
    function(parts) if (base::length(parts) >= 1L) parts[[base::length(parts)]] else NA_character_,
    character(1)
  )
  type <- vapply(
    path_parts,
    function(parts) if (base::length(parts) >= 4L) parts[[base::length(parts) - 1L]] else NA_character_,
    character(1)
  )
  activity <- ifelse(model_root %in% c("AntiDefenseFinder", "ADF"), "Anti-defense", "Defense")
  data.frame(
    model_root = .dnmb_defensefinder_na_if_empty(model_root),
    type = .dnmb_defensefinder_na_if_empty(type),
    subtype = .dnmb_defensefinder_na_if_empty(subtype),
    activity = .dnmb_defensefinder_na_if_empty(activity),
    stringsAsFactors = FALSE
  )
}

.dnmb_defensefinder_annotate_models <- function(tbl) {
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !base::nrow(tbl) || !"model_fqn" %in% base::names(tbl)) {
    return(tbl)
  }
  annotated <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  model_info <- .dnmb_defensefinder_model_parts(annotated$model_fqn)
  for (column_name in c("type", "subtype", "activity")) {
    current <- annotated[[column_name]]
    if (base::is.null(current)) {
      current <- base::rep(NA_character_, base::nrow(annotated))
    }
    annotated[[column_name]] <- dplyr::coalesce(
      .dnmb_defensefinder_na_if_empty(current),
      .dnmb_defensefinder_na_if_empty(model_info[[column_name]])
    )
  }
  annotated
}

.dnmb_defensefinder_model_definition_version <- function(model_xml) {
  if (!base::file.exists(model_xml)) {
    return(NA_character_)
  }
  lines <- tryCatch(base::readLines(model_xml, warn = FALSE), error = function(e) character())
  model_line <- lines[grep("<model ", lines, fixed = TRUE)][1]
  if (base::is.na(model_line) || !base::nzchar(model_line)) {
    return(NA_character_)
  }
  version_value <- sub('.*vers="([^"]+)".*', "\\1", model_line)
  if (identical(version_value, model_line)) {
    return(NA_character_)
  }
  version_value
}

.dnmb_defensefinder_casfinder_model_version <- function(casfinder_dir) {
  casfinder_dir <- .dnmb_defensefinder_path_first(casfinder_dir)
  if (!base::dir.exists(casfinder_dir)) {
    return(NA_character_)
  }
  preferred_xml <- base::file.path(casfinder_dir, "definitions", "CAS_Class2-Subtype-VI-X.xml")
  xml_path <- if (base::file.exists(preferred_xml)) {
    preferred_xml
  } else {
    xml_files <- base::list.files(base::file.path(casfinder_dir, "definitions"), pattern = "\\.xml$", full.names = TRUE)
    xml_files <- xml_files[base::length(xml_files) > 0L]
    if (!base::length(xml_files)) {
      return(NA_character_)
    }
    xml_files[[1]]
  }
  .dnmb_defensefinder_model_definition_version(xml_path)
}

.dnmb_defensefinder_casfinder_release_version <- function(casfinder_dir) {
  casfinder_dir <- .dnmb_defensefinder_path_first(casfinder_dir)
  metadata_path <- base::file.path(casfinder_dir, "metadata.yml")
  if (!base::file.exists(metadata_path)) {
    return(NA_character_)
  }
  lines <- tryCatch(base::readLines(metadata_path, warn = FALSE), error = function(e) character())
  vers_line <- lines[grep("^vers:", lines)][1]
  if (base::is.na(vers_line) || !base::nzchar(vers_line)) {
    return(NA_character_)
  }
  sub("^vers:\\s*", "", vers_line)
}

.dnmb_defensefinder_casfinder_is_compatible <- function(casfinder_dir) {
  model_version <- .dnmb_defensefinder_casfinder_model_version(casfinder_dir)
  base::identical(base::as.character(model_version), "2.0")
}

.dnmb_defensefinder_find_compatible_casfinder <- function(asset_urls = list()) {
  env_source <- base::Sys.getenv(.dnmb_defensefinder_casfinder_envvar(), unset = "")
  candidates <- c(
    asset_urls$casfinder_dir %||% "",
    env_source,
    "~/.macsyfinder/models/CasFinder"
  )
  candidates <- unique(base::path.expand(base::trimws(base::as.character(candidates))))
  candidates <- candidates[base::nzchar(candidates)]
  for (candidate in candidates) {
    if (.dnmb_defensefinder_casfinder_is_compatible(candidate)) {
      return(candidate)
    }
  }
  ""
}

.dnmb_defensefinder_ensure_compatible_casfinder <- function(layout, asset_urls = list()) {
  target_dir <- base::file.path(layout$models_dir, "CasFinder")
  if (.dnmb_defensefinder_casfinder_is_compatible(target_dir)) {
    return(.dnmb_defensefinder_status_row(
      "defensefinder_casfinder",
      "ok",
      base::paste0(
        target_dir,
        " (release ",
        .dnmb_defensefinder_casfinder_release_version(target_dir) %||% "unknown",
        ", model vers ",
        .dnmb_defensefinder_casfinder_model_version(target_dir) %||% "unknown",
        ")"
      )
    ))
  }

  source_dir <- .dnmb_defensefinder_find_compatible_casfinder(asset_urls = asset_urls)
  if (!base::nzchar(source_dir)) {
    return(.dnmb_defensefinder_status_row(
      "defensefinder_casfinder",
      "failed",
      base::paste0(
        "No compatible CasFinder source found. Set ",
        .dnmb_defensefinder_casfinder_envvar(),
        " or provide asset_urls$casfinder_dir with a CasFinder bundle using model vers 2.0."
      )
    ))
  }

  if (base::dir.exists(target_dir)) {
    unlink(target_dir, recursive = TRUE, force = TRUE)
  }
  .dnmb_defensefinder_copy_dir_contents(source_dir, target_dir)

  .dnmb_defensefinder_status_row(
    "defensefinder_casfinder",
    if (.dnmb_defensefinder_casfinder_is_compatible(target_dir)) "ok" else "failed",
    base::paste0(
      target_dir,
      " <= ",
      source_dir,
      " (release ",
      .dnmb_defensefinder_casfinder_release_version(target_dir) %||% "unknown",
      ", model vers ",
      .dnmb_defensefinder_casfinder_model_version(target_dir) %||% "unknown",
      ")"
    )
  )
}

.dnmb_defensefinder_cached_install_state <- function(layout, manifest = NULL) {
  if (base::is.null(manifest) || !base::isTRUE(manifest$install_ok)) {
    return(list(ok = FALSE, detail = "DefenseFinder manifest is missing or marked incomplete."))
  }

  if (!base::dir.exists(layout$repo_dir)) {
    return(list(ok = FALSE, detail = base::paste0("DefenseFinder repo is missing: ", layout$repo_dir)))
  }
  if (!base::dir.exists(layout$models_dir)) {
    return(list(ok = FALSE, detail = base::paste0("DefenseFinder models directory is missing: ", layout$models_dir)))
  }

  required_paths <- c(layout$env_python, layout$env_pip, layout$macsydata_path, layout$cli_path)
  missing <- required_paths[!base::file.exists(required_paths)]
  if (base::length(missing)) {
    return(list(ok = FALSE, detail = base::paste0("DefenseFinder runtime is missing: ", missing[[1]])))
  }

  py_run <- dnmb_run_external(
    layout$env_python,
    args = c("-c", "import sys; print(sys.executable)"),
    required = FALSE
  )
  if (!base::isTRUE(py_run$ok)) {
    return(list(ok = FALSE, detail = py_run$error %||% layout$env_python))
  }

  cli_run <- dnmb_run_external(layout$cli_path, args = "-h", required = FALSE)
  if (!base::isTRUE(cli_run$ok)) {
    return(list(ok = FALSE, detail = cli_run$error %||% layout$cli_path))
  }

  list(ok = TRUE, detail = layout$cli_path)
}

.dnmb_defensefinder_install_ready <- function(layout, manifest = NULL) {
  base::isTRUE(.dnmb_defensefinder_cached_install_state(layout, manifest = manifest)$ok)
}

.dnmb_defensefinder_candidate_python <- function() {
  candidates <- c(
    dnmb_detect_binary("python3.12", required = FALSE)$path,
    dnmb_detect_binary("python3.11", required = FALSE)$path,
    dnmb_detect_binary("python3.13", required = FALSE)$path,
    dnmb_detect_binary("python3.10", required = FALSE)$path,
    dnmb_detect_binary("python3", required = FALSE)$path
  )
  candidates <- unique(candidates[nzchar(candidates)])
  # Only use python that has pip (system python often lacks it)
  for (py in candidates) {
    test <- dnmb_run_external(py, args = c("-m", "pip", "--version"), required = FALSE)
    if (base::isTRUE(test$ok)) return(py)
  }
  if (base::length(candidates)) candidates[[1]] else ""
}

.dnmb_defensefinder_prepare_repo <- function(layout,
                                             repo_source = .dnmb_defensefinder_default_repo_url(),
                                             force = FALSE) {
  repo_source <- base::as.character(repo_source)[1]
  if (base::dir.exists(layout$repo_dir) && !base::isTRUE(force)) {
    run <- dnmb_run_external("git", args = c("-C", layout$repo_dir, "pull", "--ff-only"), required = FALSE)
    return(.dnmb_defensefinder_status_row(
      "defensefinder_repo",
      if (base::isTRUE(run$ok)) "updated" else "cached",
      if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% layout$repo_dir)
    ))
  }

  if (base::dir.exists(layout$repo_dir) && base::isTRUE(force)) {
    unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }

  if (base::dir.exists(repo_source)) {
    .dnmb_defensefinder_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_defensefinder_status_row("defensefinder_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
    required = FALSE
  )
  .dnmb_defensefinder_status_row(
    "defensefinder_repo",
    if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% repo_source)
  )
}

.dnmb_defensefinder_prepare_env <- function(layout, force = FALSE) {
  py <- .dnmb_defensefinder_candidate_python()
  if (!base::nzchar(py)) {
    return(.dnmb_defensefinder_status_row("defensefinder_python", "missing", "No python3 executable found in PATH."))
  }
  if (base::dir.exists(layout$env_dir) && base::isTRUE(force)) {
    unlink(layout$env_dir, recursive = TRUE, force = TRUE)
  }
  # Check if existing venv works (may be from different platform/python version)
  if (base::dir.exists(layout$env_dir)) {
    if (!base::file.exists(layout$env_python)) {
      # Broken symlink or missing python — nuke and recreate
      unlink(layout$env_dir, recursive = TRUE, force = TRUE)
    } else {
      test <- dnmb_run_external(layout$env_python, args = c("-c", "print('ok')"), required = FALSE)
      if (!base::isTRUE(test$ok)) {
        unlink(layout$env_dir, recursive = TRUE, force = TRUE)
      }
    }
  }
  if (!base::file.exists(layout$env_python)) {
    run <- dnmb_run_external(py, args = c("-m", "venv", layout$env_dir), required = FALSE)
    if (!base::isTRUE(run$ok)) {
      return(.dnmb_defensefinder_status_row("defensefinder_python", "failed", run$error %||% py))
    }
  }
  .dnmb_defensefinder_status_row("defensefinder_python", "ok", layout$env_python)
}

.dnmb_defensefinder_install_cli <- function(layout) {
  # pip upgrade is optional — don't fail if it doesn't work
  dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", "--upgrade", "pip"), required = FALSE)
  # Install defense-finder from cloned repo
  run <- dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", "-U", layout$repo_dir), required = FALSE)
  if (!base::isTRUE(run$ok)) {
    return(.dnmb_defensefinder_status_row("defensefinder_cli", "failed", run$error %||% layout$env_python))
  }
  .dnmb_defensefinder_status_row(
    "defensefinder_cli",
    if (base::file.exists(layout$cli_path)) "ok" else "failed",
    if (base::file.exists(layout$cli_path)) layout$cli_path else layout$env_dir
  )
}

.dnmb_defensefinder_update_models <- function(layout) {
  stage_models <- base::tempfile("dnmb-defensefinder-models-")
  base::dir.create(stage_models, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_models, recursive = TRUE, force = TRUE), add = TRUE)
  macsydata_bin <- if (base::file.exists(layout$macsydata_path)) layout$macsydata_path else ""
  if (!base::nzchar(macsydata_bin)) {
    return(.dnmb_defensefinder_status_row(
      "defensefinder_models",
      "failed",
      "macsydata is missing from the DefenseFinder virtualenv."
    ))
  }

  env_path <- base::paste(
    c(base::dirname(dnmb_detect_binary("hmmsearch", required = TRUE)$path), strsplit(base::Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]),
    collapse = .Platform$path.sep
  )
  defense_models_pkg <- base::paste0("defense-finder-models==", .dnmb_defensefinder_default_models_version())
  casfinder_pkg <- base::paste0("CasFinder==", .dnmb_defensefinder_default_casfinder_version())
  commands <- list(
    dnmb_run_external(
      macsydata_bin,
      args = c("install", "--force", "--target", stage_models, "--org", "mdmparis", defense_models_pkg),
      env = c(PATH = env_path),
      required = FALSE
    ),
    dnmb_run_external(
      macsydata_bin,
      args = c("install", "--force", "--target", stage_models, casfinder_pkg),
      env = c(PATH = env_path),
      required = FALSE
    )
  )
  run_ok <- base::all(vapply(commands, function(x) base::isTRUE(x$ok), logical(1)))
  if (run_ok) {
    if (base::dir.exists(layout$models_dir)) {
      unlink(layout$models_dir, recursive = TRUE, force = TRUE)
    }
    base::dir.create(layout$models_dir, recursive = TRUE, showWarnings = FALSE)
    .dnmb_defensefinder_copy_dir_contents(stage_models, layout$models_dir)
  }
  detail <- if (run_ok) {
    base::paste0(
      layout$models_dir,
      " (defense-finder-models ",
      .dnmb_defensefinder_default_models_version(),
      "; CasFinder ",
      .dnmb_defensefinder_default_casfinder_version(),
      ")"
    )
  } else {
    failures <- vapply(commands, function(x) x$error %||% "macsydata install failed", character(1))
    base::paste(failures, collapse = " || ")
  }
  .dnmb_defensefinder_status_row(
    "defensefinder_models",
    if (run_ok && base::length(base::list.files(layout$models_dir))) "ok" else "failed",
    detail
  )
}

dnmb_defensefinder_install_module <- function(version = .dnmb_defensefinder_default_version(),
                                              cache_root = NULL,
                                              install = TRUE,
                                              repo_url = .dnmb_defensefinder_default_repo_url(),
                                              asset_urls = NULL,
                                              force = FALSE) {
  module <- .dnmb_defensefinder_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_defensefinder_asset_layout(module_dir)
  asset_urls <- .dnmb_defensefinder_normalize_asset_urls(asset_urls)
  status <- .dnmb_defensefinder_empty_status()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  cache_state <- .dnmb_defensefinder_cached_install_state(layout, manifest = manifest)
  ready <- base::isTRUE(cache_state$ok)

  if (!ready && !base::isTRUE(force) && !base::is.null(manifest) && base::isTRUE(manifest$install_ok)) {
    status <- dplyr::bind_rows(
      status,
      .dnmb_defensefinder_status_row(
        "defensefinder_install",
        "stale",
        cache_state$detail %||% "Cached DefenseFinder runtime is not executable in the current environment."
      )
    )
  }

  if (ready && !base::isTRUE(force)) {
    casfinder_status <- .dnmb_defensefinder_ensure_compatible_casfinder(layout, asset_urls = asset_urls)
    status <- dplyr::bind_rows(status, casfinder_status)
    if (casfinder_status$status == "ok") {
      return(list(
        ok = TRUE,
        status = dplyr::bind_rows(status, .dnmb_defensefinder_status_row("defensefinder_install", "cached", module_dir)),
        files = list(repo_dir = layout$repo_dir, cli = layout$cli_path, models_dir = layout$models_dir),
        manifest = manifest
      ))
    }
  }

  if (!base::isTRUE(install)) {
    if (ready) {
      return(list(
        ok = FALSE,
        status = dplyr::bind_rows(status, .dnmb_defensefinder_status_row(
          "defensefinder_install",
          "invalid",
          "DefenseFinder assets exist but CasFinder is incompatible and module_install is FALSE."
        )),
        files = list(repo_dir = layout$repo_dir, cli = layout$cli_path, models_dir = layout$models_dir),
        manifest = manifest
      ))
    }
    return(list(ok = FALSE, status = .dnmb_defensefinder_status_row("defensefinder_install", "missing", "DefenseFinder is missing and module_install is FALSE."), files = list(), manifest = NULL))
  }

  if (base::interactive()) {
    answer <- utils::askYesNo("DefenseFinder resources are missing or may be outdated. Install or update the official DefenseFinder module now?")
    if (!base::isTRUE(answer)) {
      return(list(ok = FALSE, status = .dnmb_defensefinder_status_row("defensefinder_install", "declined", "User declined DefenseFinder install/update."), files = list(), manifest = NULL))
    }
  }

  repo_status <- .dnmb_defensefinder_prepare_repo(layout, repo_source = asset_urls$repo_dir %||% asset_urls$repo_url %||% repo_url, force = force)
  status <- dplyr::bind_rows(status, repo_status)
  if (!repo_status$status %in% c("ok", "cached", "updated")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  env_status <- .dnmb_defensefinder_prepare_env(layout, force = force)
  status <- dplyr::bind_rows(status, env_status)
  if (env_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  cli_status <- .dnmb_defensefinder_install_cli(layout)
  status <- dplyr::bind_rows(status, cli_status)
  if (cli_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  models_status <- .dnmb_defensefinder_update_models(layout)
  status <- dplyr::bind_rows(status, models_status)
  if (models_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  casfinder_status <- .dnmb_defensefinder_ensure_compatible_casfinder(layout, asset_urls = asset_urls)
  status <- dplyr::bind_rows(status, casfinder_status)
  if (casfinder_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  manifest <- list(
    install_ok = TRUE,
    repo_dir = layout$repo_dir,
    env_python = layout$env_python,
    cli_path = layout$cli_path,
    models_dir = layout$models_dir,
    models_version = .dnmb_defensefinder_default_models_version(),
    casfinder_version = .dnmb_defensefinder_casfinder_release_version(base::file.path(layout$models_dir, "CasFinder")),
    casfinder_model_version = .dnmb_defensefinder_casfinder_model_version(base::file.path(layout$models_dir, "CasFinder"))
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_defensefinder_default_version(),
    cache_root = cache_root
  )
  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_defensefinder_status_row("defensefinder_install", "ok", module_dir)),
    files = list(repo_dir = layout$repo_dir, cli = layout$cli_path, models_dir = layout$models_dir),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_defensefinder_get_module <- function(version = .dnmb_defensefinder_default_version(),
                                          cache_root = NULL,
                                          required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_defensefinder_module_name(), version, cache_root = cache_root, required = FALSE)
  if (base::is.null(manifest)) {
    if (base::isTRUE(required)) {
      base::stop("DefenseFinder module is not installed.", call. = FALSE)
    }
    return(list(ok = FALSE, manifest = NULL))
  }
  list(
    ok = TRUE,
    manifest = manifest,
    files = list(
      repo_dir = manifest$repo_dir,
      cli = manifest$cli_path,
      models_dir = manifest$models_dir,
      env_python = manifest$env_python
    )
  )
}

.dnmb_defensefinder_safe_contig_map <- function(genes) {
  contigs <- unique(base::as.character(genes$contig))
  contigs <- contigs[!base::is.na(contigs) & base::nzchar(contigs)]
  tibble::tibble(
    contig = contigs,
    safe_replicon = base::sprintf("DFREP%03d", seq_along(contigs))
  )
}

.dnmb_defensefinder_prepare_input <- function(genes, output_dir) {
  proteins <- .dnmb_prepare_query_proteins(genes)
  proteins <- proteins[order(proteins$contig, suppressWarnings(base::as.numeric(proteins$start)), suppressWarnings(base::as.numeric(proteins$end)), proteins$locus_tag), , drop = FALSE]
  rep_map <- .dnmb_defensefinder_safe_contig_map(proteins)
  proteins$safe_replicon <- rep_map$safe_replicon[match(proteins$contig, rep_map$contig)]
  proteins$defensefinder_id <- ave(
    seq_len(nrow(proteins)),
    proteins$safe_replicon,
    FUN = function(x) base::sprintf("%s_%06d", proteins$safe_replicon[x], seq_along(x))
  )
  map_tbl <- data.frame(
    defensefinder_id = proteins$defensefinder_id,
    locus_tag = proteins$locus_tag,
    contig = proteins$contig,
    start = proteins$start,
    end = proteins$end,
    direction = proteins$direction,
    product = proteins$product,
    stringsAsFactors = FALSE
  )
  faa_path <- base::file.path(output_dir, "defensefinder_input.faa")
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  con <- base::file(faa_path, open = "w")
  on.exit(base::close(con), add = TRUE)
  for (i in seq_len(nrow(proteins))) {
    strand <- if ("direction" %in% names(proteins)) base::as.character(proteins$direction[[i]]) else "NA"
    aa_len <- if ("amino_acid" %in% names(proteins)) proteins$amino_acid[[i]] else nchar(proteins$translation[[i]])
    gene_name <- if ("gene" %in% names(proteins)) base::as.character(proteins$gene[[i]]) else "NA"
    locus_tag <- proteins$locus_tag[[i]]
    protein_id <- if ("protein_id" %in% names(proteins)) base::as.character(proteins$protein_id[[i]]) else "NA"
    product <- if ("product" %in% names(proteins)) base::as.character(proteins$product[[i]]) else "NA"
    header <- base::paste(
      proteins$defensefinder_id[[i]],
      ifelse(base::is.na(strand) | !base::nzchar(strand), "NA", strand),
      "NA",
      "NA",
      ifelse(base::is.na(proteins$start[[i]]), "NA", as.character(proteins$start[[i]])),
      ifelse(base::is.na(proteins$end[[i]]), "NA", as.character(proteins$end[[i]])),
      ifelse(base::is.na(aa_len), "NA", as.character(aa_len)),
      ifelse(base::is.na(gene_name) | !base::nzchar(gene_name), "NA", gene_name),
      ifelse(base::is.na(locus_tag) | !base::nzchar(locus_tag), "NA", locus_tag),
      ifelse(base::is.na(protein_id) | !base::nzchar(protein_id), "NA", protein_id),
      ifelse(base::is.na(product) | !base::nzchar(product), "NA", product),
      sep = "\t"
    )
    base::writeLines(base::paste0(">", header), con = con)
    base::writeLines(proteins$translation[[i]], con = con)
  }
  map_path <- base::file.path(output_dir, "defensefinder_id_map.tsv")
  utils::write.table(map_tbl, file = map_path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(faa = faa_path, map = map_tbl, map_path = map_path)
}

dnmb_defensefinder_parse_genes <- function(path, id_map, systems_tbl = NULL) {
  if (!.dnmb_nonempty_file(path)) {
    return(data.frame())
  }
  genes <- tryCatch(
    utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
  if (!nrow(genes)) {
    return(data.frame())
  }
  genes <- dplyr::left_join(genes, id_map, by = c("hit_id" = "defensefinder_id"))
  genes <- .dnmb_defensefinder_annotate_models(genes)
  if (!base::is.null(systems_tbl) && base::is.data.frame(systems_tbl) && base::nrow(systems_tbl)) {
    genes <- dplyr::left_join(genes, systems_tbl, by = "sys_id")
    for (column_name in c("type", "subtype", "activity")) {
      left_name <- paste0(column_name, ".x")
      right_name <- paste0(column_name, ".y")
      if (left_name %in% names(genes) && right_name %in% names(genes)) {
        genes[[column_name]] <- dplyr::coalesce(genes[[right_name]], genes[[left_name]])
        genes[[left_name]] <- NULL
        genes[[right_name]] <- NULL
      }
    }
  }
  genes
}

dnmb_defensefinder_parse_systems <- function(path) {
  if (!.dnmb_nonempty_file(path)) {
    return(data.frame())
  }
  tryCatch(
    utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
}

dnmb_defensefinder_read_raw_best_solution <- function(path) {
  if (!.dnmb_nonempty_file(path)) {
    return(data.frame())
  }
  tryCatch(
    utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#"),
    error = function(e) data.frame()
  )
}

dnmb_defensefinder_collect_raw_genes <- function(raw_root) {
  raw_dirs <- c("DefenseFinder", "RM", "Cas", "AntiDefenseFinder")
  tables <- lapply(raw_dirs, function(dirname_value) {
    dnmb_defensefinder_read_raw_best_solution(base::file.path(raw_root, dirname_value, "best_solution.tsv"))
  })
  tables <- Filter(function(x) is.data.frame(x) && nrow(x), tables)
  if (!length(tables)) {
    return(data.frame())
  }
  .dnmb_defensefinder_annotate_models(dplyr::bind_rows(tables))
}

dnmb_defensefinder_build_systems <- function(defensefinder_genes) {
  if (base::is.null(defensefinder_genes) || !base::is.data.frame(defensefinder_genes) || !base::nrow(defensefinder_genes)) {
    return(data.frame(
      sys_id = character(),
      type = character(),
      subtype = character(),
      activity = character(),
      sys_beg = character(),
      sys_end = character(),
      protein_in_syst = character(),
      genes_count = integer(),
      name_of_profiles_in_sys = character(),
      stringsAsFactors = FALSE
    ))
  }
  tbl <- base::as.data.frame(defensefinder_genes, stringsAsFactors = FALSE)
  split_idx <- split(seq_len(nrow(tbl)), tbl$sys_id)
  rows <- lapply(split_idx, function(idx) {
    sub <- tbl[idx, , drop = FALSE]
    sub <- sub[order(suppressWarnings(as.numeric(sub$hit_pos))), , drop = FALSE]
    first_or_na <- function(column_name) {
      if (!column_name %in% names(sub)) {
        return(NA_character_)
      }
      value <- sub[[column_name]][[1]]
      if (is.null(value) || length(value) == 0L) {
        return(NA_character_)
      }
      as.character(value)
    }
    data.frame(
      sys_id = first_or_na("sys_id"),
      type = first_or_na("type"),
      subtype = first_or_na("subtype"),
      activity = first_or_na("activity"),
      sys_beg = first_or_na("hit_id"),
      sys_end = as.character(sub$hit_id[[nrow(sub)]]),
      protein_in_syst = base::paste(sub$hit_id, collapse = ","),
      genes_count = nrow(sub),
      name_of_profiles_in_sys = base::paste(sort(sub$gene_name), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

dnmb_defensefinder_normalize_hits <- function(genes_tbl) {
  if (base::is.null(genes_tbl) || !base::is.data.frame(genes_tbl) || !base::nrow(genes_tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  tbl <- base::as.data.frame(genes_tbl, stringsAsFactors = FALSE)
  for (column_name in c("type", "subtype", "activity", "sys_id", "gene_name", "hit_status", "sys_wholeness", "sys_score", "genes_count", "sys_beg", "sys_end", "protein_in_syst", "name_of_profiles_in_sys", "hit_i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov", "hit_begin_match", "hit_end_match", "counterpart", "used_in")) {
    if (base::is.null(tbl[[column_name]])) {
      tbl[[column_name]] <- base::rep(NA_character_, base::nrow(tbl))
    }
  }
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$locus_tag)
  tbl$support <- vapply(seq_len(nrow(tbl)), function(i) {
    parts <- c(
      if (!base::is.na(tbl$sys_id[[i]]) && base::nzchar(tbl$sys_id[[i]])) base::paste0("sys_id=", tbl$sys_id[[i]]) else NA_character_,
      if (!base::is.na(tbl$gene_name[[i]]) && base::nzchar(tbl$gene_name[[i]])) base::paste0("gene_name=", tbl$gene_name[[i]]) else NA_character_,
      if (!base::is.na(tbl$hit_status[[i]]) && base::nzchar(tbl$hit_status[[i]])) base::paste0("hit_status=", tbl$hit_status[[i]]) else NA_character_,
      if (!base::is.na(tbl$hit_i_eval[[i]]) && base::nzchar(base::as.character(tbl$hit_i_eval[[i]]))) base::paste0("i_eval=", tbl$hit_i_eval[[i]]) else NA_character_,
      if (!base::is.na(tbl$hit_score[[i]]) && base::nzchar(base::as.character(tbl$hit_score[[i]]))) base::paste0("hit_score=", tbl$hit_score[[i]]) else NA_character_,
      if (!base::is.na(tbl$hit_profile_cov[[i]]) && base::nzchar(base::as.character(tbl$hit_profile_cov[[i]]))) base::paste0("profile_cov=", tbl$hit_profile_cov[[i]]) else NA_character_,
      if (!base::is.na(tbl$hit_seq_cov[[i]]) && base::nzchar(base::as.character(tbl$hit_seq_cov[[i]]))) base::paste0("seq_cov=", tbl$hit_seq_cov[[i]]) else NA_character_,
      if (!base::is.na(tbl$name_of_profiles_in_sys[[i]]) && base::nzchar(tbl$name_of_profiles_in_sys[[i]])) base::paste0("profiles=", tbl$name_of_profiles_in_sys[[i]]) else NA_character_,
      if (!base::is.na(tbl$protein_in_syst[[i]]) && base::nzchar(tbl$protein_in_syst[[i]])) base::paste0("proteins=", tbl$protein_in_syst[[i]]) else NA_character_
    )
    parts <- parts[!base::is.na(parts)]
    if (!base::length(parts)) {
      return(NA_character_)
    }
    base::paste(parts, collapse = "; ")
  }, character(1))

  data.frame(
    query = tbl$query,
    source = "defensefinder",
    family_system = "DefenseFinder",
    family_id = base::as.character(tbl$type),
    hit_label = base::as.character(tbl$subtype),
    enzyme_role = base::as.character(tbl$activity),
    evidence_mode = base::as.character(tbl$hit_status),
    substrate_label = NA_character_,
    support = tbl$support,
    typing_eligible = base::as.character(tbl$hit_status) == "mandatory",
    system_type = base::as.character(tbl$type),
    system_subtype = base::as.character(tbl$subtype),
    system_activity = base::as.character(tbl$activity),
    system_id = base::as.character(tbl$sys_id),
    gene_name = base::as.character(tbl$gene_name),
    hit_status = base::as.character(tbl$hit_status),
    system_wholeness = suppressWarnings(base::as.numeric(tbl$sys_wholeness)),
    system_score = suppressWarnings(base::as.numeric(tbl$sys_score)),
    genes_count = suppressWarnings(base::as.integer(tbl$genes_count)),
    system_begin = base::as.character(tbl$sys_beg),
    system_end = base::as.character(tbl$sys_end),
    protein_in_system = base::as.character(tbl$protein_in_syst),
    profiles_in_system = base::as.character(tbl$name_of_profiles_in_sys),
    hit_i_eval = suppressWarnings(base::as.numeric(tbl$hit_i_eval)),
    hit_score = suppressWarnings(base::as.numeric(tbl$hit_score)),
    hit_profile_cov = suppressWarnings(base::as.numeric(tbl$hit_profile_cov)),
    hit_seq_cov = suppressWarnings(base::as.numeric(tbl$hit_seq_cov)),
    hit_begin_match = suppressWarnings(base::as.integer(tbl$hit_begin_match)),
    hit_end_match = suppressWarnings(base::as.integer(tbl$hit_end_match)),
    counterpart = base::as.character(tbl$counterpart),
    used_in = base::as.character(tbl$used_in),
    stringsAsFactors = FALSE
  )
}

.dnmb_defensefinder_hits_for_output <- function(hits) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  out$system_score <- suppressWarnings(base::as.numeric(out$system_score))
  out$system_wholeness <- suppressWarnings(base::as.numeric(out$system_wholeness))
  out$hit_score <- suppressWarnings(base::as.numeric(out$hit_score))
  out$hit_i_eval <- suppressWarnings(base::as.numeric(out$hit_i_eval))
  status_rank <- ifelse(out$hit_status == "mandatory", 2L, ifelse(out$hit_status == "accessory", 1L, 0L))
  out <- out[order(
    out$query,
    -ifelse(base::is.na(out$system_score), -Inf, out$system_score),
    -ifelse(base::is.na(out$system_wholeness), -Inf, out$system_wholeness),
    -status_rank,
    -ifelse(base::is.na(out$hit_score), -Inf, out$hit_score),
    ifelse(base::is.na(out$hit_i_eval), Inf, out$hit_i_eval)
  ), , drop = FALSE]
  out <- out[!duplicated(out$query), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_defensefinder_output_table <- function(genes, hits) {
  selected_hits <- .dnmb_defensefinder_hits_for_output(hits)
  out <- .dnmb_module_output_table(genes = genes, hits = selected_hits)
  drop_cols <- base::intersect(c("family_id", "hit_label", "enzyme_role", "evidence_mode", "substrate_label", "typing_eligible"), base::names(out))
  if (base::length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(
      c("system_type", "system_subtype", "system_activity", "system_id", "gene_name", "hit_status", "system_wholeness", "system_score", "genes_count", "system_begin", "system_end", "profiles_in_system", "hit_i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov", "support"),
      base::names(out)
    ),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "system_type", "system_subtype", "system_activity", "system_id", "gene_name", "hit_status", "system_wholeness", "system_score", "genes_count", "system_begin", "system_end", "profiles_in_system", "hit_i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov", "support"))
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_defensefinder_module <- function(genes,
                                          output_dir,
                                          version = .dnmb_defensefinder_default_version(),
                                          cache_root = NULL,
                                          install = TRUE,
                                          repo_url = .dnmb_defensefinder_default_repo_url(),
                                          asset_urls = NULL,
                                          cpu = 1L,
                                          genbank = NULL,
                                          include_antidefense = TRUE) {
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- base::file.path(output_dir, "defensefinder_module_trace.log")
  status <- .dnmb_defensefinder_empty_status()
  stale_outputs <- base::file.path(
    output_dir,
    c(
      "defensefinder_best_solution_genes.tsv",
      "defensefinder_systems.tsv",
      "defensefinder_hmmer.tsv",
      "defensefinder_id_map.tsv"
    )
  )
  base::unlink(stale_outputs[base::file.exists(stale_outputs)], force = TRUE)

  install_result <- dnmb_defensefinder_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    repo_url = repo_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), genes = data.frame(), systems = data.frame(), hmmer = data.frame(), output_table = data.frame()))
  }

  module <- dnmb_defensefinder_get_module(version = version, cache_root = cache_root, required = TRUE)
  stage_dir <- base::tempfile("dnmb-defensefinder-run-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  stage_models <- base::file.path(stage_dir, "models")
  .dnmb_defensefinder_copy_dir_contents(module$files$models_dir, stage_models)

  input <- .dnmb_defensefinder_prepare_input(genes = genes, output_dir = stage_dir)
  stage_out <- base::file.path(stage_dir, "out")
  base::dir.create(stage_out, recursive = TRUE, showWarnings = FALSE)

  cmd_args <- c(
    "run",
    input$faa,
    "-o", stage_out,
    "--db-type", "gembase",
    "-w", base::as.character(base::as.integer(cpu)[1]),
    "--models-dir", stage_models,
    "--preserve-raw"
  )
  if (base::isTRUE(include_antidefense)) {
    cmd_args <- c(cmd_args, "--antidefensefinder")
  }
  .dnmb_defensefinder_trace(trace_log, base::sprintf("[%s] %s", base::Sys.time(), .dnmb_format_command(module$files$cli, cmd_args)))
  run <- dnmb_run_external(
    module$files$cli,
    args = cmd_args,
    env = c(
      PATH = base::paste(
        c(base::dirname(dnmb_detect_binary("hmmsearch", required = TRUE)$path), strsplit(base::Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]),
        collapse = .Platform$path.sep
      )
    ),
    required = FALSE
  )
  .dnmb_defensefinder_trace_result(trace_log, run, label = "defensefinder_run")

  raw_root <- base::file.path(stage_out, "defense-finder-tmp")
  genes_tbl_raw <- dnmb_defensefinder_collect_raw_genes(raw_root)
  systems_tbl_raw <- dnmb_defensefinder_build_systems(genes_tbl_raw)
  genes_export <- base::list.files(stage_out, pattern = "_defense_finder_genes\\.tsv$", full.names = TRUE)
  systems_export <- base::list.files(stage_out, pattern = "_defense_finder_systems\\.tsv$", full.names = TRUE)
  hmmer_export <- base::list.files(stage_out, pattern = "_defense_finder_hmmer\\.tsv$", full.names = TRUE)

  genes_path <- base::file.path(output_dir, "defensefinder_best_solution_genes.tsv")
  systems_path <- base::file.path(output_dir, "defensefinder_systems.tsv")
  hmmer_path <- base::file.path(output_dir, "defensefinder_hmmer.tsv")
  if (nrow(genes_tbl_raw)) {
    utils::write.table(genes_tbl_raw, file = genes_path, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (base::length(genes_export) && base::file.exists(genes_export[[1]])) {
    base::file.copy(genes_export[[1]], genes_path, overwrite = TRUE)
  }
  if (nrow(systems_tbl_raw)) {
    utils::write.table(systems_tbl_raw, file = systems_path, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (base::length(systems_export) && base::file.exists(systems_export[[1]])) {
    base::file.copy(systems_export[[1]], systems_path, overwrite = TRUE)
  }
  if (base::length(hmmer_export) && base::file.exists(hmmer_export[[1]])) {
    base::file.copy(hmmer_export[[1]], hmmer_path, overwrite = TRUE)
  }
  for (path in c(input$map_path)) {
    if (base::file.exists(path)) {
      base::file.copy(path, base::file.path(output_dir, base::basename(path)), overwrite = TRUE)
    }
  }

  systems_tbl <- if (nrow(systems_tbl_raw)) systems_tbl_raw else dnmb_defensefinder_parse_systems(systems_path)
  genes_tbl <- dnmb_defensefinder_parse_genes(genes_path, id_map = input$map, systems_tbl = systems_tbl)
  hmmer_tbl <- if (.dnmb_nonempty_file(hmmer_path)) {
    tryCatch(
      utils::read.delim(hmmer_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) data.frame()
    )
  } else {
    data.frame()
  }
  hits <- dnmb_defensefinder_normalize_hits(genes_tbl)
  output_table <- .dnmb_defensefinder_output_table(genes = genes, hits = hits)
  exports_present <- base::length(genes_export) > 0L || base::length(systems_export) > 0L || base::length(hmmer_export) > 0L
  effective_ok <- base::isTRUE(run$ok) || nrow(genes_tbl_raw) > 0L || nrow(genes_tbl) > 0L || exports_present

  status <- dplyr::bind_rows(
    status,
    .dnmb_defensefinder_status_row(
      "defensefinder_run",
      if (nrow(genes_tbl_raw) || nrow(genes_tbl)) {
        "ok"
      } else if (effective_ok) {
        "empty"
      } else {
        "failed"
      },
      if (effective_ok) genes_path else (run$error %||% stage_out)
    )
  )

  list(
    ok = effective_ok,
    status = status,
    files = list(
      trace_log = trace_log,
      genes = base::file.path(output_dir, base::basename(genes_path)),
      systems = base::file.path(output_dir, base::basename(systems_path)),
      hmmer = base::file.path(output_dir, base::basename(hmmer_path)),
      id_map = base::file.path(output_dir, base::basename(input$map_path))
    ),
    genes = genes_tbl,
    systems = systems_tbl,
    hmmer = hmmer_tbl,
    hits = hits,
    output_table = output_table
  )
}
