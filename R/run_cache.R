.dnmb_file_signature <- function(paths) {
  if (is.null(paths)) {
    return(NULL)
  }

  paths <- path.expand(as.character(paths))
  paths <- paths[nzchar(paths)]
  paths <- unique(paths[file.exists(paths)])
  if (!length(paths)) {
    return(NULL)
  }

  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  info <- file.info(paths)
  files <- data.frame(
    file = basename(paths),
    path = unname(paths),
    size = unname(as.numeric(info$size)),
    md5 = unname(as.character(tools::md5sum(paths))),
    stringsAsFactors = FALSE
  )
  files <- files[order(files$md5, files$size, files$file, files$path), , drop = FALSE]
  rownames(files) <- NULL

  list(
    signature_version = 1L,
    file_count = nrow(files),
    files = files
  )
}

.dnmb_file_signature_key <- function(signature) {
  if (is.null(signature) || is.null(signature$files)) {
    return("")
  }

  files <- as.data.frame(signature$files, stringsAsFactors = FALSE)
  if (!nrow(files) || !all(c("md5", "size") %in% names(files))) {
    return("")
  }

  sort_file <- if ("file" %in% names(files)) as.character(files$file) else rep("", nrow(files))
  sort_path <- if ("path" %in% names(files)) as.character(files$path) else rep("", nrow(files))
  files <- files[order(files$md5, files$size, sort_file, sort_path), , drop = FALSE]

  paste(files$md5, files$size, collapse = "|")
}

.dnmb_same_file_signature <- function(current, previous) {
  current_key <- .dnmb_file_signature_key(current)
  previous_key <- .dnmb_file_signature_key(previous)

  nzchar(current_key) && identical(current_key, previous_key)
}

.dnmb_detect_genbank_inputs <- function(path = getwd()) {
  path <- path.expand(as.character(path)[1])

  candidates <-
    if (file.exists(path) && !dir.exists(path)) {
      if (grepl("\\.(gbk|gb|gbff)$", path, ignore.case = TRUE)) {
        path
      } else {
        character(0)
      }
    } else if (dir.exists(path)) {
      list.files(
        path,
        pattern = "\\.(gbk|gb|gbff)$",
        full.names = TRUE,
        ignore.case = TRUE
      )
    } else {
      character(0)
    }

  if (!length(candidates)) {
    return(character(0))
  }

  candidates <- candidates[!grepl("^~\\$", basename(candidates))]
  if (!length(candidates)) {
    return(character(0))
  }

  candidates <- normalizePath(candidates, winslash = "/", mustWork = TRUE)
  candidates[order(basename(candidates), candidates)]
}

.dnmb_genbank_input_signature <- function(path = getwd()) {
  .dnmb_file_signature(.dnmb_detect_genbank_inputs(path))
}

.dnmb_db_manifest_identity <- function(module, version, cache_root = NULL) {
  manifest <- tryCatch(
    dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE),
    error = function(e) NULL
  )

  if (is.null(manifest)) {
    return(list(
      module = module,
      version = version,
      installed = FALSE
    ))
  }

  fields <- unclass(manifest)
  fields$module <- module
  fields$version <- version
  fields$installed <- if (is.null(fields$install_ok)) TRUE else isTRUE(fields$install_ok)

  drop_names <- c(
    "module_dir",
    "manifest_path",
    "repo_dir",
    "path_dir",
    "env_python",
    "cli",
    "cli_path",
    "models_dir",
    "app_dir",
    "pretrained_dir",
    "db_dir"
  )
  drop_mask <- names(fields) %in% drop_names | grepl("(_dir|_path)$", names(fields))
  fields[drop_mask] <- NULL
  fields
}

.dnmb_collect_module_db_signatures <- function(module_aliases,
                                               module_version = NULL,
                                               module_cache_root = NULL,
                                               module_Prophage_backend = .dnmb_run_default_prophage_backend()) {
  if (!length(module_aliases)) {
    return(list())
  }

  signatures <- list()
  requested_version <- if (!is.null(module_version) && nzchar(trimws(as.character(module_version)[1]))) {
    trimws(as.character(module_version)[1])
  } else {
    NULL
  }

  current_version <- requested_version %||% "current"

  if ("EggNOG" %in% module_aliases) {
    signatures$EggNOG <- .dnmb_db_manifest_identity("eggnog", "data", cache_root = module_cache_root)
  }
  if ("CLEAN" %in% module_aliases) {
    signatures$CLEAN <- .dnmb_db_manifest_identity("clean", requested_version %||% "split100", cache_root = module_cache_root)
  }
  if ("DefenseFinder" %in% module_aliases) {
    signatures$DefenseFinder <- .dnmb_db_manifest_identity("defensefinder", current_version, cache_root = module_cache_root)
  }
  if ("PADLOC" %in% module_aliases) {
    signatures$PADLOC <- .dnmb_db_manifest_identity("padloc", current_version, cache_root = module_cache_root)
  }
  if ("DefensePredictor" %in% module_aliases) {
    signatures$DefensePredictor <- .dnmb_db_manifest_identity("defensepredictor", current_version, cache_root = module_cache_root)
  }
  if ("GapMind" %in% module_aliases) {
    signatures$GapMind <- list(
      aa = .dnmb_db_manifest_identity("gapmind", "aa", cache_root = module_cache_root),
      carbon = .dnmb_db_manifest_identity("gapmind", "carbon", cache_root = module_cache_root)
    )
  }
  if ("MEROPS" %in% module_aliases) {
    signatures$MEROPS <- .dnmb_db_manifest_identity("merops", current_version, cache_root = module_cache_root)
  }
  if ("dbCAN" %in% module_aliases) {
    signatures$dbCAN <- .dnmb_db_manifest_identity("dbcan", current_version, cache_root = module_cache_root)
  }
  if ("PAZy" %in% module_aliases) {
    signatures$PAZy <- .dnmb_db_manifest_identity("pazy", current_version, cache_root = module_cache_root)
  }
  if ("Prophage" %in% module_aliases) {
    signatures$Prophage <- .dnmb_db_manifest_identity("prophage", module_Prophage_backend %||% "phispy", cache_root = module_cache_root)
  }
  if ("ISelement" %in% module_aliases) {
    signatures$ISelement <- list(installed = TRUE, module = "iselement", version = current_version)
  }
  if ("REBASEfinder" %in% module_aliases) {
    signatures$REBASEfinder <- list(installed = TRUE, module = "rebasefinder", version = "embedded")
  }

  signatures
}

.dnmb_module_stage_signature <- function(genbank_signature,
                                         module_aliases,
                                         module_version = NULL,
                                         module_cache_root = NULL,
                                         module_install = TRUE,
                                         module_base_url = NULL,
                                         module_asset_urls = NULL,
                                         module_Prophage_backend = .dnmb_run_default_prophage_backend(),
                                         iselement_analysis_depth = "full",
                                         iselement_related_genbanks = NULL,
                                         iselement_related_metadata = NULL,
                                         iselement_auto_discover_related = TRUE,
                                         iselement_max_related = 5L) {
  aliases <- sort(unique(as.character(module_aliases)))
  aliases <- aliases[nzchar(aliases)]

  list(
    signature_version = 1L,
    genbank_signature = genbank_signature,
    module_aliases = aliases,
    requested = list(
      module_version = if (is.null(module_version)) NULL else as.character(module_version)[1],
      module_install = isTRUE(module_install),
      module_base_url = if (is.null(module_base_url)) NULL else as.character(module_base_url)[1],
      module_asset_urls = module_asset_urls,
      module_Prophage_backend = as.character(module_Prophage_backend)[1],
      iselement_analysis_depth = as.character(iselement_analysis_depth)[1],
      iselement_auto_discover_related = isTRUE(iselement_auto_discover_related),
      iselement_max_related = as.integer(iselement_max_related)[1]
    ),
    related_inputs = list(
      iselement_related_genbanks = .dnmb_file_signature(iselement_related_genbanks),
      iselement_related_metadata = .dnmb_file_signature(iselement_related_metadata)
    ),
    db_state = .dnmb_collect_module_db_signatures(
      module_aliases = aliases,
      module_version = module_version,
      module_cache_root = module_cache_root,
      module_Prophage_backend = module_Prophage_backend
    )
  )
}

.dnmb_module_stage_cache_path <- function(wd = getwd()) {
  file.path(wd, ".dnmb_module_stage_cache.rds")
}

.dnmb_read_module_stage_cache <- function(wd = getwd()) {
  cache_path <- .dnmb_module_stage_cache_path(wd)
  if (!file.exists(cache_path)) {
    return(NULL)
  }

  tryCatch(readRDS(cache_path), error = function(e) NULL)
}

.dnmb_write_module_stage_cache <- function(wd = getwd(),
                                           signature,
                                           module_results) {
  if (is.null(signature) || is.null(module_results)) {
    return(invisible(NULL))
  }

  cache_path <- .dnmb_module_stage_cache_path(wd)
  payload <- list(
    written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    signature = signature,
    module_results = module_results
  )
  saveRDS(payload, cache_path)
  invisible(cache_path)
}

.dnmb_module_stage_cache_status <- function(wd = getwd(), signature = NULL) {
  cache <- .dnmb_read_module_stage_cache(wd)
  if (is.null(cache) || is.null(cache$signature) || is.null(cache$module_results)) {
    return(list(reusable = FALSE, reason = "missing_cache", cache = cache))
  }

  if (is.null(signature) || !identical(signature, cache$signature)) {
    return(list(reusable = FALSE, reason = "signature_changed", cache = cache))
  }

  list(reusable = TRUE, reason = "matching_signature", cache = cache)
}

.dnmb_eggnog_state <- function(wd = getwd()) {
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)
  result_files <- list.files(
    wd,
    pattern = "\\.emapper\\.annotations\\.xlsx$",
    full.names = TRUE
  )
  result_files <- result_files[!grepl("(^|/|\\\\)~\\$", basename(result_files))]

  list(
    wd = wd,
    metadata_path = file.path(wd, ".dnmb_eggnog_input_signature.rds"),
    result_files = result_files
  )
}

.dnmb_eggnog_has_parseable_results <- function(wd = getwd()) {
  state <- .dnmb_eggnog_state(wd)
  if (!length(state$result_files)) {
    return(FALSE)
  }

  any(vapply(state$result_files, function(path) {
    header <- tryCatch(
      names(openxlsx::read.xlsx(path, startRow = 3, rows = 3, colNames = TRUE)),
      error = function(e) character()
    )
    all(c("query", "COG_category") %in% header)
  }, logical(1)))
}

.dnmb_read_eggnog_signature <- function(wd = getwd()) {
  metadata_path <- .dnmb_eggnog_state(wd)$metadata_path
  if (!file.exists(metadata_path)) {
    return(NULL)
  }

  tryCatch(readRDS(metadata_path), error = function(e) NULL)
}

.dnmb_write_eggnog_signature <- function(wd = getwd(), genbank_signature) {
  state <- .dnmb_eggnog_state(wd)
  if (is.null(genbank_signature) || !length(state$result_files)) {
    return(invisible(NULL))
  }

  payload <- list(
    input_mode = "genbank_workdir",
    written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    genbank_signature = genbank_signature,
    result_files = basename(state$result_files)
  )
  saveRDS(payload, state$metadata_path)
  invisible(state$metadata_path)
}

.dnmb_eggnog_reuse_status <- function(wd = getwd(),
                                      genbank_signature = NULL,
                                      allow_external_without_metadata = FALSE) {
  state <- .dnmb_eggnog_state(wd)
  metadata <- .dnmb_read_eggnog_signature(wd)
  has_results <- length(state$result_files) > 0L
  parseable <- .dnmb_eggnog_has_parseable_results(wd)

  if (!is.null(metadata)) {
    saved_signature <- metadata$genbank_signature
    matches_input <- .dnmb_same_file_signature(genbank_signature, saved_signature)

    return(list(
      reusable = isTRUE(matches_input) && has_results && parseable,
      has_metadata = TRUE,
      matches_input = matches_input,
      result_files = state$result_files,
      metadata_path = state$metadata_path,
      reason = if (isTRUE(matches_input)) {
        if (has_results && parseable) "matching_genbank" else if (has_results) "unparseable_results" else "results_missing"
      } else {
        "input_changed"
      }
    ))
  }

  if (isTRUE(allow_external_without_metadata) && has_results && parseable) {
    return(list(
      reusable = TRUE,
      has_metadata = FALSE,
      matches_input = NA,
      result_files = state$result_files,
      metadata_path = state$metadata_path,
      reason = "external_results_without_metadata"
    ))
  }

  list(
    reusable = FALSE,
    has_metadata = FALSE,
    matches_input = NA,
    result_files = state$result_files,
    metadata_path = state$metadata_path,
    reason = if (has_results && !parseable) "unparseable_results" else if (has_results) "results_without_metadata" else "missing_results"
  )
}

.dnmb_interproscan_state <- function(wd = getwd()) {
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)
  module_output_dir <- file.path(wd, "dnmb_module_interproscan")
  legacy_output_dir <- file.path(wd, "dnmb_interproscan")
  output_dir <- if (dir.exists(module_output_dir) || !dir.exists(legacy_output_dir)) module_output_dir else legacy_output_dir

  list(
    wd = wd,
    output_dir = output_dir,
    module_output_dir = module_output_dir,
    legacy_output_dir = legacy_output_dir,
    metadata_path = file.path(output_dir, ".input_signature.rds"),
    output_tsv = file.path(output_dir, "interproscan_results.tsv"),
    output_tsv_sites = file.path(output_dir, "interproscan_results.tsv.sites"),
    root_tsv = file.path(wd, "interproscan_results.tsv"),
    root_tsv_sites = file.path(wd, "interproscan_results.tsv.sites")
  )
}

.dnmb_read_interproscan_signature <- function(wd = getwd()) {
  state <- .dnmb_interproscan_state(wd)
  for (metadata_path in c(file.path(state$module_output_dir, ".input_signature.rds"),
                          file.path(state$legacy_output_dir, ".input_signature.rds"))) {
    if (file.exists(metadata_path)) {
      parsed <- tryCatch(readRDS(metadata_path), error = function(e) NULL)
      if (!is.null(parsed)) return(parsed)
    }
  }
  NULL
}

.dnmb_write_interproscan_signature <- function(output_dir, genbank_signature) {
  if (is.null(output_dir) || !nzchar(trimws(output_dir)) || is.null(genbank_signature)) {
    return(invisible(NULL))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata_path <- file.path(output_dir, ".input_signature.rds")
  payload <- list(
    input_mode = "genbank_workdir",
    written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    genbank_signature = genbank_signature
  )
  saveRDS(payload, metadata_path)
  invisible(metadata_path)
}

.dnmb_interproscan_reuse_status <- function(wd = getwd(),
                                            genbank_signature = NULL,
                                            allow_external_without_metadata = TRUE) {
  state <- .dnmb_interproscan_state(wd)
  has_output_results <- file.exists(state$output_tsv)
  has_module_results <- file.exists(file.path(state$module_output_dir, "interproscan_results.tsv"))
  has_legacy_results <- file.exists(file.path(state$legacy_output_dir, "interproscan_results.tsv"))
  has_root_results <- file.exists(state$root_tsv)
  metadata <- .dnmb_read_interproscan_signature(wd)

  if (!is.null(metadata)) {
    saved_signature <- metadata$genbank_signature
    matches_input <- .dnmb_same_file_signature(genbank_signature, saved_signature)

    if (matches_input && has_output_results) {
      return(list(
        reusable = TRUE,
        has_metadata = TRUE,
        matches_input = TRUE,
        source_dir = state$output_dir,
        reason = "matching_genbank"
      ))
    }

    if (matches_input && has_root_results) {
      return(list(
        reusable = TRUE,
        has_metadata = TRUE,
        matches_input = TRUE,
        source_dir = state$wd,
        reason = "matching_genbank_root_copy"
      ))
    }

    return(list(
      reusable = FALSE,
      has_metadata = TRUE,
      matches_input = matches_input,
      source_dir = NULL,
      reason = if (matches_input) "results_missing" else "input_changed"
    ))
  }

  if (isTRUE(allow_external_without_metadata) && (has_output_results || has_module_results || has_legacy_results)) {
    source_dir <- if (has_module_results) state$module_output_dir else if (has_legacy_results) state$legacy_output_dir else state$output_dir
    return(list(
      reusable = TRUE,
      has_metadata = FALSE,
      matches_input = NA,
      source_dir = source_dir,
      reason = "external_output_results"
    ))
  }

  if (isTRUE(allow_external_without_metadata) && has_root_results) {
    return(list(
      reusable = TRUE,
      has_metadata = FALSE,
      matches_input = NA,
      source_dir = state$wd,
      reason = "external_root_results"
    ))
  }

  list(
    reusable = FALSE,
    has_metadata = FALSE,
    matches_input = NA,
    source_dir = NULL,
    reason = if (has_output_results || has_module_results || has_legacy_results || has_root_results) "results_without_metadata" else "missing_results"
  )
}
