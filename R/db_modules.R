`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) {
    y
  } else {
    x
  }
}

.dnmb_db_cache_root <- function(cache_root = NULL, create = FALSE) {
  root <- cache_root
  if (is.null(root) || !nzchar(trimws(root))) {
    # Check DNMB_CACHE_ROOT env var (for Docker / shared cache)
    env_root <- Sys.getenv("DNMB_CACHE_ROOT", unset = "")
    if (nzchar(env_root)) {
      root <- env_root
    } else {
      shared_root <- path.expand("~/.dnmb-cache")
      if (dir.exists(shared_root) || isTRUE(create)) {
        root <- shared_root
      } else {
        root <- tools::R_user_dir("DNMB", which = "cache")
      }
    }
  }
  root <- path.expand(as.character(root)[1])
  root <- file.path(root, "db_modules")

  if (isTRUE(create)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }

  if (dir.exists(root)) {
    normalizePath(root, winslash = "/", mustWork = FALSE)
  } else {
    root
  }
}

#' Get the DNMB database home directory
#'
#' @param create Logical; create the directory when `TRUE`.
#' @param cache_root Optional override path. When `NULL`, DNMB resolves the
#'   cache root in this order: `DNMB_CACHE_ROOT`, `~/.dnmb-cache` (preferred
#'   shared cache for local and Docker runs), then `tools::R_user_dir("DNMB",
#'   which = "cache")`.
#'
#' @return Character scalar path to the DNMB database home.
#' @export
dnmb_db_home <- function(create = FALSE, cache_root = NULL) {
  .dnmb_db_cache_root(cache_root = cache_root, create = create)
}

.dnmb_db_validate_key <- function(x, label) {
  value <- trimws(as.character(x)[1])
  if (is.na(value) || !nzchar(value)) {
    stop("`", label, "` must be a non-empty string.", call. = FALSE)
  }
  if (grepl("[/\\\\]", value) || grepl("(^|[.])[.]($|[.])", value, perl = TRUE)) {
    stop("`", label, "` must not contain path separators or parent-directory markers.", call. = FALSE)
  }
  value
}

.dnmb_db_module_dir <- function(module, version, cache_root = NULL, create = FALSE) {
  module <- .dnmb_db_validate_key(module, "module")
  version <- .dnmb_db_validate_key(version, "version")
  module_dir <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = create), module, version)

  if (isTRUE(create)) {
    dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (dir.exists(module_dir)) {
    normalizePath(module_dir, winslash = "/", mustWork = FALSE)
  } else {
    module_dir
  }
}

.dnmb_db_module_root <- function(module, cache_root = NULL, create = FALSE) {
  module <- .dnmb_db_validate_key(module, "module")
  root <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = create), module)
  if (isTRUE(create)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  if (dir.exists(root)) {
    normalizePath(root, winslash = "/", mustWork = FALSE)
  } else {
    root
  }
}

.dnmb_db_prune_module_versions <- function(module,
                                           keep_versions,
                                           cache_root = NULL,
                                           preserve = c("cache"),
                                           verbose = TRUE) {
  module_root <- .dnmb_db_module_root(module, cache_root = cache_root, create = FALSE)
  if (!dir.exists(module_root)) {
    return(character())
  }

  keep_versions <- trimws(as.character(keep_versions))
  preserve <- trimws(as.character(preserve))
  keep <- unique(c(keep_versions, preserve))
  keep <- keep[!is.na(keep) & nzchar(keep)]

  entries <- list.files(module_root, full.names = TRUE, no.. = TRUE)
  if (!length(entries)) {
    return(character())
  }

  is_dir <- file.info(entries)$isdir %in% TRUE
  dir_entries <- entries[is_dir]
  if (!length(dir_entries)) {
    return(character())
  }

  dir_names <- basename(dir_entries)
  prune_paths <- dir_entries[!(dir_names %in% keep)]
  if (!length(prune_paths)) {
    return(character())
  }

  unlink(prune_paths, recursive = TRUE, force = TRUE)
  removed <- prune_paths[!file.exists(prune_paths)]
  if (isTRUE(verbose) && length(removed)) {
    message(
      "[DNMB] Pruned old cache versions for ",
      module,
      ": ",
      paste(basename(removed), collapse = ", ")
    )
  }
  invisible(removed)
}

.dnmb_db_autoprune_default_versions <- function(module,
                                                version,
                                                default_version,
                                                cache_root = NULL,
                                                preserve = c("cache"),
                                                verbose = TRUE) {
  version <- trimws(as.character(version)[1])
  default_version <- trimws(as.character(default_version)[1])
  if (is.na(version) || !nzchar(version) || is.na(default_version) || !nzchar(default_version)) {
    return(character())
  }
  if (!identical(version, default_version)) {
    return(character())
  }
  .dnmb_db_prune_module_versions(
    module = module,
    keep_versions = default_version,
    cache_root = cache_root,
    preserve = preserve,
    verbose = verbose
  )
}

.dnmb_db_manifest_path <- function(module, version, cache_root = NULL) {
  file.path(.dnmb_db_module_dir(module, version, cache_root = cache_root, create = FALSE), "manifest.rds")
}

dnmb_db_write_manifest <- function(module,
                                   version,
                                   manifest = list(),
                                   cache_root = NULL,
                                   overwrite = FALSE) {
  if (!is.list(manifest)) {
    stop("`manifest` must be a list.", call. = FALSE)
  }

  module <- .dnmb_db_validate_key(module, "module")
  version <- .dnmb_db_validate_key(version, "version")
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  manifest_path <- .dnmb_db_manifest_path(module, version, cache_root = cache_root)

  if (file.exists(manifest_path) && !isTRUE(overwrite)) {
    stop("Manifest already exists for module `", module, "` version `", version, "`.", call. = FALSE)
  }

  payload <- c(
    list(
      module = module,
      version = version,
      module_dir = module_dir,
      manifest_path = manifest_path,
      written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    ),
    manifest
  )
  class(payload) <- c("dnmb_db_manifest", "list")
  saveRDS(payload, file = manifest_path)
  invisible(manifest_path)
}

dnmb_db_read_manifest <- function(module, version, cache_root = NULL, required = FALSE) {
  module <- .dnmb_db_validate_key(module, "module")
  version <- .dnmb_db_validate_key(version, "version")
  manifest_path <- .dnmb_db_manifest_path(module, version, cache_root = cache_root)

  if (!file.exists(manifest_path)) {
    if (isTRUE(required)) {
      stop("Manifest not found for module `", module, "` version `", version, "`.", call. = FALSE)
    }
    return(NULL)
  }

  manifest <- readRDS(manifest_path)
  if (is.null(manifest$module)) {
    manifest$module <- module
  }
  if (is.null(manifest$version)) {
    manifest$version <- version
  }
  # Always use the current module_dir (may differ from manifest if cache moved or Docker)
  current_module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = FALSE)
  stored_module_dir <- manifest$module_dir
  manifest$module_dir <- current_module_dir
  manifest$manifest_path <- manifest_path
  manifest <- .dnmb_db_rebase_manifest_paths(manifest, stored_module_dir, current_module_dir)
  class(manifest) <- unique(c("dnmb_db_manifest", class(manifest)))
  manifest
}

.dnmb_db_rebase_manifest_paths <- function(manifest, stored_root, current_root) {
  if (is.null(stored_root) || !nzchar(stored_root) || identical(stored_root, current_root)) {
    return(manifest)
  }
  stored_prefix <- sub("/+$", "", stored_root)
  current_prefix <- sub("/+$", "", current_root)
  rebase_one <- function(value) {
    if (!is.character(value) || !length(value)) return(value)
    hit <- startsWith(value, paste0(stored_prefix, "/")) | value == stored_prefix
    if (!any(hit, na.rm = TRUE)) return(value)
    value[hit] <- paste0(current_prefix, substring(value[hit], nchar(stored_prefix) + 1L))
    value
  }
  skip <- c("module", "version", "module_dir", "manifest_path", "resource_base_url", "repo_url")
  for (nm in setdiff(names(manifest), skip)) {
    val <- manifest[[nm]]
    if (is.character(val) && length(val)) {
      manifest[[nm]] <- rebase_one(val)
    }
  }
  manifest
}

#' Check DB freshness and warn if stale
#'
#' @param module Module name.
#' @param version Module version.
#' @param cache_root Optional cache root.
#' @param max_age_days Days before warning (default: 90).
#' @param verbose Print messages.
#' @return Logical; TRUE if fresh, FALSE if stale.
#' Check DB freshness against remote server
#'
#' Compares local manifest against remote source to detect updates.
#' @param module Module name.
#' @param version Module version.
#' @param cache_root Optional cache root.
#' @param verbose Print messages.
#' @return Logical; TRUE if up-to-date, FALSE if update available.
dnmb_db_check_freshness <- function(module, version, cache_root = NULL, verbose = TRUE) {
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)

  # If no manifest, check if module dir exists with files (e.g., eggnog, interproscan)
  if (is.null(manifest)) {
    mod_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = FALSE)
    if (dir.exists(mod_dir) && length(list.files(mod_dir)) > 0L) {
      # Use newest file modification date as proxy
      files <- list.files(mod_dir, full.names = TRUE, recursive = FALSE)
      newest <- max(file.info(files)$mtime, na.rm = TRUE)
      db_date <- as.Date(newest)
      age_days <- as.integer(Sys.Date() - db_date)
      manifest <- list(version = version, written_at = format(newest, "%Y-%m-%dT%H:%M:%SZ"))
    } else {
      if (isTRUE(verbose)) message("[DNMB] ", module, ": not installed — will be set up on first run")
      return(TRUE)
    }
  }

  db_date <- tryCatch(as.Date(sub("T.*", "", manifest$written_at)), error = function(e) NA)
  age_days <- if (!is.na(db_date)) as.integer(Sys.Date() - db_date) else NA

  # Check remote for updates (module-specific)
  remote_info <- tryCatch(.dnmb_db_check_remote_version(module, manifest), error = function(e) NULL)

  if (isTRUE(verbose)) {
    local_ver <- manifest$version %||% version
    date_str <- if (!is.na(db_date)) format(db_date, "%Y-%m-%d") else "unknown"
    ver_label <- .dnmb_db_manifest_version_label(manifest)
    ver_suffix <- if (nzchar(ver_label)) paste0(" [", ver_label, "]") else ""
    if (!is.null(remote_info) && isTRUE(remote_info$update_available)) {
      message("[DNMB] ", module, " (", local_ver, ", ", date_str, ")", ver_suffix,
              " — update available (", remote_info$remote_version, ")")
    } else if (!is.null(remote_info) && !isTRUE(remote_info$update_available)) {
      message("[DNMB] ", module, " (", local_ver, ", ", date_str, ")", ver_suffix, " — latest")
    } else {
      message("[DNMB] ", module, " (", local_ver, ", ", date_str, ")", ver_suffix)
    }
  }

  is.null(remote_info) || !isTRUE(remote_info$update_available)
}

#' Build a short label summarizing tool / data versions stored in a
#' module manifest, for display in [DNMB] status messages.
#' Returns "" when no informative field is present.
.dnmb_db_manifest_version_label <- function(manifest) {
  if (!is.list(manifest)) return("")
  fields <- c(
    db_version              = "db",
    db_release              = "db",
    models_version          = "models",
    casfinder_version       = "casfinder",
    casfinder_model_version = "casfinder_models",
    emapper_version         = "emapper",
    resolved_release_version = "release",
    dbcan_release_date      = "dbcan",
    cazydb_release_date     = "cazydb",
    model                   = "model",
    mode                    = "mode"
  )
  parts <- character()
  for (key in names(fields)) {
    v <- manifest[[key]]
    if (is.null(v) || !length(v)) next
    val <- as.character(v)[1]
    if (is.na(val) || !nzchar(val)) next
    parts <- c(parts, paste0(fields[[key]], "=", val))
  }
  paste(parts, collapse = ", ")
}

.dnmb_db_check_remote_version <- function(module, manifest) {
  check_fn <- switch(module,
    dbcan = .dnmb_db_remote_check_dbcan,
    merops = .dnmb_db_remote_check_url_changed,
    pazy = .dnmb_db_remote_check_url_changed,
    defensefinder = .dnmb_db_remote_check_defensefinder,
    interproscan = .dnmb_db_remote_check_interproscan,
    eggnog = .dnmb_db_remote_check_eggnog,
    gapmind = .dnmb_db_remote_check_gapmind,
    clean = .dnmb_db_remote_check_clean,
    prophage = .dnmb_db_remote_check_prophage,
    NULL
  )
  if (is.null(check_fn)) return(NULL)
  tryCatch(check_fn(manifest), error = function(e) NULL)
}

.dnmb_db_remote_check_dbcan <- function(manifest) {
  # Scrape dbCAN release page for latest version
  html <- tryCatch(readLines("https://pro.unl.edu/dbCAN2/download/Databases/", warn = FALSE), error = function(e) NULL)
  if (is.null(html)) return(NULL)
  ver_match <- regmatches(paste(html, collapse = ""), regexpr("V[0-9]+", paste(html, collapse = "")))
  if (!length(ver_match)) return(NULL)
  remote_ver <- ver_match[[1]]
  local_ver <- manifest$version %||% ""
  list(remote_version = remote_ver, update_available = !identical(toupper(local_ver), toupper(remote_ver)))
}

.dnmb_db_remote_check_url_changed <- function(manifest) {
  # Compare stored ETag/content-length against remote HEAD
  if (is.null(manifest$remote_asset_state)) return(NULL)
  saved <- manifest$remote_asset_state
  url <- saved$url %||% NULL
  if (is.null(url) || !nzchar(url)) return(NULL)
  headers <- tryCatch({
    run <- dnmb_run_external("curl", args = c("-sI", "-L", "-k", url), required = FALSE)
    if (!isTRUE(run$ok)) return(NULL)
    run$stdout
  }, error = function(e) NULL)
  if (is.null(headers)) return(NULL)
  remote_etag <- sub(".*ETag:\\s*", "", grep("ETag:", headers, value = TRUE, ignore.case = TRUE)[1])
  remote_length <- sub(".*Content-Length:\\s*", "", grep("Content-Length:", headers, value = TRUE, ignore.case = TRUE)[1])
  local_etag <- saved$etag %||% ""
  local_length <- as.character(saved$content_length %||% "")
  changed <- (!is.na(remote_etag) && nzchar(remote_etag) && remote_etag != local_etag) ||
             (!is.na(remote_length) && nzchar(remote_length) && remote_length != local_length)
  list(remote_version = if (changed) "newer version" else "same", update_available = changed)
}

.dnmb_db_remote_check_eggnog <- function(manifest) {
  # Check eggnog-mapper PyPI for latest version
  json <- tryCatch({
    con <- url("https://pypi.org/pypi/eggnog-mapper/json")
    on.exit(close(con))
    jsonlite::fromJSON(readLines(con, warn = FALSE))
  }, error = function(e) NULL)
  if (is.null(json)) return(NULL)
  remote_ver <- json$info$version
  local_ver <- manifest$emapper_version %||% manifest$version %||% ""
  list(remote_version = remote_ver, update_available = !identical(local_ver, remote_ver))
}

.dnmb_db_remote_check_gapmind <- function(manifest) {
  # GapMind is hosted at genomics.lbl.gov — check if files changed
  url <- "https://papers.genomics.lbl.gov/tmp/path.aa/steps.db"
  headers <- tryCatch({
    run <- dnmb_run_external("curl", args = c("-sI", "-L", "-k", url), required = FALSE)
    if (isTRUE(run$ok)) run$stdout else NULL
  }, error = function(e) NULL)
  if (is.null(headers)) return(NULL)
  remote_date <- sub(".*Last-Modified:\\s*", "", grep("Last-Modified:", headers, value = TRUE, ignore.case = TRUE)[1])
  local_date <- manifest$written_at %||% ""
  if (is.na(remote_date) || !nzchar(remote_date)) return(NULL)
  list(remote_version = trimws(remote_date), update_available = FALSE) # informational only
}

.dnmb_db_remote_check_clean <- function(manifest) {
  # CLEAN uses a fixed pretrained model — check GitHub for new releases
  json <- tryCatch({
    con <- url("https://api.github.com/repos/tttianhao/CLEAN/releases/latest")
    on.exit(close(con))
    jsonlite::fromJSON(readLines(con, warn = FALSE))
  }, error = function(e) NULL)
  if (is.null(json)) return(NULL)
  remote_tag <- json$tag_name %||% json$name %||% ""
  local_ver <- manifest$version %||% ""
  list(remote_version = remote_tag, update_available = nzchar(remote_tag) && remote_tag != local_ver)
}

.dnmb_db_remote_check_prophage <- function(manifest) {
  # PhiSpy — check GitHub releases
  json <- tryCatch({
    con <- url("https://api.github.com/repos/linsalrob/PhiSpy/releases/latest")
    on.exit(close(con))
    jsonlite::fromJSON(readLines(con, warn = FALSE))
  }, error = function(e) NULL)
  if (is.null(json)) return(NULL)
  remote_tag <- json$tag_name %||% json$name %||% ""
  list(remote_version = remote_tag, update_available = FALSE) # informational
}

.dnmb_db_remote_check_defensefinder <- function(manifest) {
  # Check both tool and models repos
  tool_ver <- tryCatch({
    con <- url("https://api.github.com/repos/mdmparis/defense-finder/releases/latest")
    on.exit(close(con))
    j <- jsonlite::fromJSON(readLines(con, warn = FALSE))
    j$tag_name %||% ""
  }, error = function(e) "")
  model_ver <- tryCatch({
    con <- url("https://api.github.com/repos/mdmparis/defense-finder-models/releases/latest")
    on.exit(close(con))
    j <- jsonlite::fromJSON(readLines(con, warn = FALSE))
    j$tag_name %||% ""
  }, error = function(e) "")
  remote_info <- paste0("tool:", tool_ver, " models:", model_ver)
  local_models <- manifest$models_version %||% ""
  updated <- (nzchar(model_ver) && model_ver != local_models) || (nzchar(tool_ver) && !grepl(gsub("^v", "", tool_ver), manifest$written_at %||% "", fixed = TRUE))
  list(remote_version = remote_info, update_available = updated)
}

.dnmb_db_remote_check_interproscan <- function(manifest) {
  # Check GitHub for latest InterProScan5 release
  json <- tryCatch({
    con <- url("https://api.github.com/repos/ebi-pf-team/interproscan/releases/latest")
    on.exit(close(con))
    jsonlite::fromJSON(readLines(con, warn = FALSE))
  }, error = function(e) NULL)
  if (is.null(json)) return(NULL)
  remote_tag <- json$tag_name %||% json$name %||% ""
  local_ver <- manifest$version %||% ""
  list(remote_version = remote_tag, update_available = nzchar(remote_tag) && remote_tag != local_ver)
}

dnmb_db_list_registry <- function(cache_root = NULL) {
  registry_root <- .dnmb_db_cache_root(cache_root = cache_root, create = FALSE)
  if (!dir.exists(registry_root)) {
    return(tibble::tibble(
      module = character(),
      version = character(),
      module_dir = character(),
      manifest_path = character(),
      written_at = character()
    ))
  }

  manifest_paths <- list.files(
    registry_root,
    pattern = "^manifest\\.rds$",
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )

  rows <- lapply(manifest_paths, function(path) {
    manifest <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(manifest)) {
      return(NULL)
    }

    tibble::tibble(
      module = as.character(manifest$module %||% NA_character_),
      version = as.character(manifest$version %||% NA_character_),
      module_dir = as.character(manifest$module_dir %||% dirname(path)),
      manifest_path = normalizePath(path, winslash = "/", mustWork = FALSE),
      written_at = as.character(manifest$written_at %||% NA_character_)
    )
  })

  dplyr::bind_rows(rows)
}

dnmb_detect_binary <- function(binary, required = FALSE) {
  binary <- .dnmb_db_validate_key(binary, "binary")
  path <- Sys.which(binary)
  found <- nzchar(path)
  message <- if (found) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    paste0("Binary not found in PATH: ", binary)
  }

  result <- list(
    binary = binary,
    path = if (found) normalizePath(path, winslash = "/", mustWork = FALSE) else "",
    found = found,
    message = message
  )
  class(result) <- c("dnmb_binary_detection", "list")

  if (!found && isTRUE(required)) {
    stop(message, call. = FALSE)
  }

  result
}

.dnmb_nonempty_lines <- function(x) {
  x <- trimws(as.character(x))
  x[nzchar(x)]
}

.dnmb_parse_tool_version <- function(x) {
  x <- .dnmb_nonempty_lines(x)
  if (!length(x)) {
    return(NA_character_)
  }

  match <- regexpr("[0-9]+(?:\\.[0-9]+)+(?:\\+)?", x, perl = TRUE)
  index <- which(match > 0L)[1]
  if (is.na(index)) {
    return(NA_character_)
  }

  regmatches(x[[index]], regexpr("[0-9]+(?:\\.[0-9]+)+(?:\\+)?", x[[index]], perl = TRUE))
}

.dnmb_version_at_least <- function(version, target) {
  version <- as.character(version)[1]
  target <- as.character(target)[1]
  if (is.na(version) || !nzchar(version) || is.na(target) || !nzchar(target)) {
    return(FALSE)
  }

  utils::compareVersion(sub("\\+$", "", version), sub("\\+$", "", target)) >= 0L
}

.dnmb_format_command <- function(command, args = character()) {
  paste(c(as.character(command)[1], as.character(args)), collapse = " ")
}

.dnmb_compact_output <- function(x, max_lines = 6L) {
  x <- .dnmb_nonempty_lines(x)
  if (!length(x)) {
    return("")
  }
  x <- x[seq_len(min(length(x), as.integer(max_lines)[1]))]
  paste(x, collapse = " | ")
}

.dnmb_makeblastdb_capabilities <- function() {
  detection <- dnmb_detect_binary("makeblastdb", required = FALSE)
  result <- list(
    binary = "makeblastdb",
    path = detection$path,
    found = isTRUE(detection$found),
    version = NA_character_,
    version_source = NA_character_,
    version_raw = character(),
    supports_blastdb_version = FALSE,
    supports_parse_seqids = FALSE,
    help_source = NA_character_,
    help_raw = character()
  )
  class(result) <- c("dnmb_makeblastdb_capabilities", "list")

  if (!isTRUE(result$found)) {
    return(result)
  }

  version_run <- dnmb_run_external("makeblastdb", args = "-version", required = FALSE)
  version_lines <- .dnmb_nonempty_lines(c(version_run$stdout, version_run$stderr))
  result$version_raw <- version_lines
  version <- .dnmb_parse_tool_version(version_lines)
  if (!is.na(version) && nzchar(version)) {
    result$version <- version
    result$version_source <- "makeblastdb -version"
  }

  help_run <- dnmb_run_external("makeblastdb", args = "-help", required = FALSE)
  help_lines <- .dnmb_nonempty_lines(c(help_run$stdout, help_run$stderr))
  result$help_raw <- help_lines
  if (length(help_lines)) {
    help_blob <- paste(help_lines, collapse = "\n")
    result$help_source <- "makeblastdb -help"
    result$supports_blastdb_version <- grepl("(^|\\s)-blastdb_version(\\s|$)", help_blob, perl = TRUE)
    result$supports_parse_seqids <- grepl("(^|\\s)-parse_seqids(\\s|$)", help_blob, perl = TRUE)
  }

  if (!isTRUE(result$supports_blastdb_version) && .dnmb_version_at_least(result$version, "2.10.0")) {
    result$supports_blastdb_version <- TRUE
  }
  if (!isTRUE(result$supports_parse_seqids) && .dnmb_version_at_least(result$version, "2.2.0")) {
    result$supports_parse_seqids <- TRUE
  }

  result
}

dnmb_run_external <- function(command,
                              args = character(),
                              required = FALSE,
                              wd = NULL,
                              env = character(),
                              stream_stderr = FALSE) {
  command <- as.character(command)[1]
  if (is.na(command) || !nzchar(trimws(command))) {
    stop("`command` must be a non-empty string.", call. = FALSE)
  }

  args <- as.character(args)
  if (!is.null(wd)) {
    wd <- path.expand(as.character(wd)[1])
    if (!dir.exists(wd)) {
      stop("`wd` does not exist: ", wd, call. = FALSE)
    }
  }

  resolved <- if (grepl("[/\\\\]", command)) {
    path.expand(command)
  } else {
    dnmb_detect_binary(command, required = required)$path
  }

  if (!nzchar(resolved)) {
    result <- list(
      command = command,
      resolved_command = "",
      args = args,
      wd = wd,
      status = NA_integer_,
      stdout = character(),
      stderr = character(),
      ok = FALSE,
      error = paste0("Binary not found in PATH: ", command)
    )
    class(result) <- c("dnmb_external_command_result", "list")
    return(result)
  }

  stdout_path <- tempfile("dnmb-stdout-")
  stderr_path <- tempfile("dnmb-stderr-")
  on.exit(unlink(c(stdout_path, stderr_path), force = TRUE), add = TRUE)

  old_wd <- NULL
  if (!is.null(wd)) {
    old_wd <- getwd()
    setwd(wd)
    on.exit(setwd(old_wd), add = TRUE)
  }

  shell_command <- paste(
    c(shQuote(resolved, type = "sh"), vapply(args, shQuote, character(1), type = "sh")),
    collapse = " "
  )

  if (length(env)) {
    env_names <- names(env)
    if (is.null(env_names) || any(!nzchar(env_names))) {
      stop("`env` must be a named character vector.", call. = FALSE)
    }
    old_env <- Sys.getenv(env_names, unset = NA_character_)
    do.call(Sys.setenv, as.list(env))
    on.exit({
      for (i in seq_along(env_names)) {
        if (is.na(old_env[[i]])) {
          Sys.unsetenv(env_names[[i]])
        } else {
          do.call(Sys.setenv, stats::setNames(list(old_env[[i]]), env_names[[i]]))
        }
      }
    }, add = TRUE)
  }

  if (isTRUE(stream_stderr)) {
    # Stream stderr to console (for long-running tools like emapper)
    status <- tryCatch(
      system(
        paste0(
          shell_command,
          " 1>", shQuote(stdout_path, type = "sh"),
          " 2>&1 | tee ", shQuote(stderr_path, type = "sh"), " >&2"
        ),
        intern = FALSE,
        wait = TRUE
      ),
      error = function(e) e
    )
  } else {
    status <- tryCatch(
      system(
        paste0(
          shell_command,
          " 1>", shQuote(stdout_path, type = "sh"),
          " 2>", shQuote(stderr_path, type = "sh")
        ),
        intern = FALSE,
        wait = TRUE
      ),
      error = function(e) e
    )
  }

  stdout <- if (file.exists(stdout_path)) readLines(stdout_path, warn = FALSE) else character()
  stderr <- if (file.exists(stderr_path)) readLines(stderr_path, warn = FALSE) else character()
  failed <- inherits(status, "error") || !identical(as.integer(status), 0L)
  error <- NULL
  exit_status <- if (inherits(status, "error")) NA_integer_ else as.integer(status)

  if (inherits(status, "error")) {
    error <- conditionMessage(status)
  } else if (!identical(exit_status, 0L)) {
    detail <- c(stderr, stdout)
    detail <- detail[nzchar(detail)]
    error <- paste0(
      "Command failed [", exit_status, "]: ",
      paste(c(command, args), collapse = " "),
      if (length(detail)) paste0("\n", paste(detail, collapse = "\n")) else ""
    )
  }

  result <- list(
    command = command,
    resolved_command = resolved,
    args = args,
    wd = wd,
    status = exit_status,
    stdout = stdout,
    stderr = stderr,
    ok = !failed,
    error = error
  )
  class(result) <- c("dnmb_external_command_result", "list")

  if (!result$ok && isTRUE(required)) {
    stop(result$error, call. = FALSE)
  }

  result
}

.dnmb_nonempty_file <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  file.exists(path) && isTRUE(file.info(path)$size > 0)
}

.dnmb_copy_local_asset <- function(source, dest) {
  source <- path.expand(source)
  if (!file.exists(source)) {
    return(list(ok = FALSE, method = "local_copy", error = paste0("Local asset not found: ", source)))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(source, dest, overwrite = TRUE)
  list(
    ok = isTRUE(ok) && file.exists(dest),
    method = "local_copy",
    error = if (isTRUE(ok) && file.exists(dest)) NULL else paste0("Failed to copy local asset: ", source)
  )
}

.dnmb_download_asset <- function(url, dest, insecure = FALSE) {
  if (!is.character(url) || !length(url) || is.na(url[[1]]) || !nzchar(trimws(url[[1]]))) {
    return(list(ok = FALSE, method = "invalid", error = "Asset URL must be a non-empty string."))
  }

  url <- trimws(url[[1]])
  if (file.exists(url)) {
    return(.dnmb_copy_local_asset(url, dest))
  }
  if (startsWith(url, "file://")) {
    return(.dnmb_copy_local_asset(sub("^file://", "", url), dest))
  }

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  curl_detection <- dnmb_detect_binary("curl", required = FALSE)
  if (isTRUE(curl_detection$found)) {
    curl_args <- c(if (isTRUE(insecure)) "-k", "-L", "-o", dest, url)
    curl_run <- dnmb_run_external("curl", curl_args, required = FALSE)
    if (isTRUE(curl_run$ok) && .dnmb_nonempty_file(dest)) {
      return(list(ok = TRUE, method = "curl", error = NULL, command = curl_run))
    }
  }

  wget_detection <- dnmb_detect_binary("wget", required = FALSE)
  if (isTRUE(wget_detection$found)) {
    wget_args <- c(if (isTRUE(insecure)) "--no-check-certificate", "-O", dest, url)
    wget_run <- dnmb_run_external("wget", wget_args, required = FALSE)
    if (isTRUE(wget_run$ok) && .dnmb_nonempty_file(dest)) {
      return(list(ok = TRUE, method = "wget", error = NULL, command = wget_run))
    }
  }

  utils_ok <- tryCatch({
    utils::download.file(url = url, destfile = dest, quiet = TRUE, mode = "wb")
    .dnmb_nonempty_file(dest)
  }, error = function(error) FALSE)
  if (isTRUE(utils_ok)) {
    return(list(ok = TRUE, method = "download.file", error = NULL))
  }

  list(ok = FALSE, method = "download", error = paste0("Failed to download asset: ", url))
}

.dnmb_parse_http_headers <- function(lines) {
  lines <- as.character(lines)
  if (!length(lines)) {
    return(character())
  }

  split_points <- c(0L, which(grepl("^HTTP/", lines)), length(lines) + 1L)
  if (length(split_points) < 2L) {
    block <- lines
  } else {
    starts <- split_points[split_points < length(lines) + 1L]
    ends <- c(starts[-1L] - 1L, length(lines))
    block <- lines[starts[[length(starts)]]:ends[[length(ends)]]]
  }

  headers <- block[grepl("^[^:]+:\\s*", block)]
  if (!length(headers)) {
    return(character())
  }

  keys <- tolower(sub(":.*$", "", headers))
  values <- trimws(sub("^[^:]+:\\s*", "", headers))
  stats::setNames(values, keys)
}

.dnmb_remote_asset_metadata <- function(url, insecure = FALSE) {
  url <- trimws(as.character(url)[1])
  result <- list(
    url = url,
    ok = FALSE,
    last_modified = NA_character_,
    content_length = NA_real_,
    etag = NA_character_,
    headers = character(),
    error = NULL
  )
  class(result) <- c("dnmb_remote_asset_metadata", "list")

  if (is.na(url) || !nzchar(url) || file.exists(url) || startsWith(url, "file://")) {
    result$error <- "Remote metadata lookup applies only to remote URLs."
    return(result)
  }

  curl_detection <- dnmb_detect_binary("curl", required = FALSE)
  if (!isTRUE(curl_detection$found)) {
    wget_detection <- dnmb_detect_binary("wget", required = FALSE)
    if (!isTRUE(wget_detection$found)) {
      result$error <- "Neither curl nor wget is available."
      return(result)
    }
  }

  headers <- character()
  run <- NULL

  if (isTRUE(curl_detection$found)) {
    curl_args <- c(if (isTRUE(insecure)) "-k", "-I", "-L", "-sS", url)
    run <- dnmb_run_external("curl", curl_args, required = FALSE)
    headers <- .dnmb_parse_http_headers(c(run$stdout, run$stderr))
  }

  if (!length(headers)) {
    wget_detection <- dnmb_detect_binary("wget", required = FALSE)
    if (isTRUE(wget_detection$found)) {
      wget_args <- c(
        if (isTRUE(insecure)) "--no-check-certificate",
        "--server-response",
        "--spider",
        url
      )
      run <- dnmb_run_external("wget", wget_args, required = FALSE)
      headers <- .dnmb_parse_http_headers(c(run$stdout, run$stderr))
    }
  }

  result$headers <- headers
  if (!length(headers)) {
    result$error <- if (!is.null(run)) run$error %||% "Failed to retrieve remote asset metadata." else "Failed to retrieve remote asset metadata."
    return(result)
  }

  result$ok <- TRUE
  result$last_modified <- if ("last-modified" %in% names(headers)) unname(headers["last-modified"][[1]]) else NA_character_
  content_length <- suppressWarnings(as.numeric(if ("content-length" %in% names(headers)) unname(headers["content-length"][[1]]) else NA_character_))
  result$content_length <- if (is.na(content_length)) NA_real_ else content_length
  result$etag <- if ("etag" %in% names(headers)) unname(headers["etag"][[1]]) else NA_character_
  result
}

dnmb_resolve_query_faa <- function(genbank = NULL,
                                   output_dir = NULL,
                                   fallback_filename = "query_proteins.faa") {
  candidates <- character()

  if (!is.null(genbank) && nzchar(trimws(as.character(genbank)[1])) && file.exists(genbank)) {
    genbank <- normalizePath(genbank, winslash = "/", mustWork = TRUE)
    gb_dir <- dirname(genbank)
    gb_base <- tools::file_path_sans_ext(basename(genbank))
    candidates <- c(
      candidates,
      file.path(gb_dir, paste0(gb_base, ".faa")),
      file.path(gb_dir, paste0(sub("\\.1$", "", gb_base), ".faa"))
    )
  }

  if (!is.null(output_dir) && nzchar(trimws(as.character(output_dir)[1]))) {
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
    candidates <- c(
      candidates,
      file.path(output_dir, fallback_filename),
      file.path(dirname(output_dir), fallback_filename)
    )
  }

  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    return(NULL)
  }

  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}
