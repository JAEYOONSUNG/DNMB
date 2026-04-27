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
    "written_at",
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

.dnmb_rebasefinder_cache_identity <- function(cache_root = NULL) {
  cache_dir <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = FALSE), "rebasefinder", "cache")
  files <- c(
    file.path(cache_dir, "rebase_data.rds"),
    file.path(cache_dir, "REBASE_protein_seqs.txt"),
    file.path(cache_dir, "rebase_db.fasta"),
    file.path(cache_dir, "rebase_bairoch_lookup.rds")
  )
  files <- files[file.exists(files)]
  info <- if (length(files)) file.info(files) else data.frame()
  file_state <- if (length(files)) {
    data.frame(
      file = basename(files),
      size = unname(as.numeric(info$size)),
      mtime = unname(format(info$mtime, "%Y-%m-%d %H:%M:%S %z")),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(file = character(), size = numeric(), mtime = character())
  }
  file_state <- file_state[order(file_state$file), , drop = FALSE]
  rownames(file_state) <- NULL
  list(
    installed = file.exists(file.path(cache_dir, "rebase_data.rds")),
    module = "rebasefinder",
    version = "embedded",
    cache_dir = normalizePath(cache_dir, winslash = "/", mustWork = FALSE),
    files = file_state
  )
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
  if ("dbAPIS" %in% module_aliases) {
    signatures$dbAPIS <- .dnmb_db_manifest_identity("dbapis", current_version, cache_root = module_cache_root)
  }
  if ("AcrFinder" %in% module_aliases) {
    signatures$AcrFinder <- .dnmb_db_manifest_identity("acrfinder", current_version, cache_root = module_cache_root)
  }
  if ("Promotech" %in% module_aliases) {
    signatures$Promotech <- .dnmb_db_manifest_identity("promotech", current_version, cache_root = module_cache_root)
  }
  if ("mRNAcal" %in% module_aliases) {
    signatures$mRNAcal <- list(installed = TRUE, module = "mrnacal", version = "embedded")
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
  if ("PhiSpy" %in% module_aliases) {
    signatures$PhiSpy <- .dnmb_db_manifest_identity("prophage", "phispy", cache_root = module_cache_root)
  }
  if ("VirSorter2" %in% module_aliases) {
    signatures$VirSorter2 <- .dnmb_db_manifest_identity("prophage", "virsorter2", cache_root = module_cache_root)
  }
  if ("PIDE" %in% module_aliases) {
    signatures$PIDE <- .dnmb_db_manifest_identity("prophage", "pide", cache_root = module_cache_root)
  }
  if ("ISelement" %in% module_aliases) {
    signatures$ISelement <- list(installed = TRUE, module = "iselement", version = current_version)
  }
  if ("REBASEfinder" %in% module_aliases) {
    signatures$REBASEfinder <- .dnmb_rebasefinder_cache_identity(cache_root = module_cache_root)
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
                                         module_DefenseFinder_antidefense = TRUE,
                                         module_Prophage_backend = .dnmb_run_default_prophage_backend(),
                                         promotech_predictions = NULL,
                                         promotech_model = "RF-HOT",
                                         promotech_threshold = 0.5,
                                         promotech_max_distance = 300L,
                                         promotech_test_samples = NULL,
                                         promotech_python = "python3",
                                         promotech_download_model = TRUE,
                                         promotech_model_base_url = .dnmb_promotech_default_model_base_url(),
                                         mrnacal_upstream = 60L,
                                         mrnacal_downstream = 60L,
                                         mrnacal_rnafold_path = NULL,
                                         mrnacal_require_rnafold = TRUE,
                                         mrnacal_sd_seed = NULL,
                                         mrnacal_top_folds = 12L,
                                         iselement_analysis_depth = "full",
                                         iselement_related_genbanks = NULL,
                                         iselement_related_metadata = NULL,
                                         iselement_auto_discover_related = TRUE,
                                         iselement_max_related = 5L) {
  aliases <- sort(unique(as.character(module_aliases)))
  aliases <- aliases[nzchar(aliases)]
  promotech_asset_urls <- if (is.list(module_asset_urls) && "Promotech" %in% names(module_asset_urls)) {
    module_asset_urls[["Promotech"]]
  } else {
    module_asset_urls
  }
  promotech_asset_predictions <- character()
  if (is.list(promotech_asset_urls)) {
    prediction_names <- intersect(names(promotech_asset_urls), c("predictions", "predictions_path", "prediction_csv"))
    if (length(prediction_names)) {
      promotech_asset_predictions <- unlist(promotech_asset_urls[prediction_names], use.names = FALSE)
    }
  }

  list(
    signature_version = 1L,
    genbank_signature = genbank_signature,
    module_aliases = aliases,
    requested = list(
      module_version = if (is.null(module_version)) NULL else as.character(module_version)[1],
      module_install = isTRUE(module_install),
      module_base_url = if (is.null(module_base_url)) NULL else as.character(module_base_url)[1],
      module_asset_urls = module_asset_urls,
      module_DefenseFinder_antidefense = isTRUE(module_DefenseFinder_antidefense),
      module_Prophage_backend = as.character(module_Prophage_backend)[1],
      promotech_model = as.character(promotech_model)[1],
      promotech_threshold = as.numeric(promotech_threshold)[1],
      promotech_max_distance = as.integer(promotech_max_distance)[1],
      promotech_test_samples = if (is.null(promotech_test_samples)) NULL else as.integer(promotech_test_samples)[1],
      promotech_python = as.character(promotech_python)[1],
      promotech_download_model = isTRUE(promotech_download_model),
      promotech_model_base_url = as.character(promotech_model_base_url)[1],
      mrnacal_upstream = as.integer(mrnacal_upstream)[1],
      mrnacal_downstream = as.integer(mrnacal_downstream)[1],
      mrnacal_rnafold_path = if (is.null(mrnacal_rnafold_path)) NULL else as.character(mrnacal_rnafold_path)[1],
      mrnacal_require_rnafold = isTRUE(mrnacal_require_rnafold),
      mrnacal_sd_seed = if (is.null(mrnacal_sd_seed)) NULL else as.character(mrnacal_sd_seed)[1],
      mrnacal_top_folds = as.integer(mrnacal_top_folds)[1],
      mrnacal_algorithm = if ("mRNAcal" %in% aliases) "rnaplfold_local_tir_ncs30_tircore_intsd_caitai_v7" else NULL,
      iselement_analysis_depth = as.character(iselement_analysis_depth)[1],
      iselement_auto_discover_related = isTRUE(iselement_auto_discover_related),
      iselement_max_related = as.integer(iselement_max_related)[1]
    ),
    related_inputs = list(
      promotech_predictions = .dnmb_file_signature(c(promotech_predictions, promotech_asset_predictions)),
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

.dnmb_module_stage_expected_artifacts <- function(alias, wd = getwd()) {
  alias <- as.character(alias)[1]
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)

  switch(
    alias,
    EggNOG = list(
      mode = "any",
      paths = c(
        file.path(wd, "dnmb_module_eggnog", "eggnog_out.emapper.annotations.xlsx"),
        file.path(wd, "dnmb_module_eggnog", "eggnog_out.emapper.annotations")
      )
    ),
    CLEAN = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_clean", "clean_module_trace.log"))
    ),
    DefenseFinder = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_defensefinder", "defensefinder_systems.tsv"))
    ),
    dbAPIS = list(
      mode = "any",
      paths = c(file.path(wd, "dnmb_module_dbapis", "dbapis_hits.tsv"))
    ),
    AcrFinder = list(
      mode = "any",
      paths = c(file.path(wd, "dnmb_module_acrfinder", "acrfinder_homology.tsv"))
    ),
    Promotech = list(
      mode = "all",
      paths = c(
        file.path(wd, "dnmb_module_promotech", "promotech_predictions.tsv"),
        file.path(wd, "dnmb_module_promotech", "promotech_mapped_hits.tsv"),
        file.path(wd, "dnmb_module_promotech", "promotech_promoter_feature_for_gb")
      )
    ),
    mRNAcal = list(
      mode = "all",
      paths = c(
        file.path(wd, "dnmb_module_mrnacal", "mrnacal_translation_efficiency.tsv")
      )
    ),
    PADLOC = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_padloc", "padloc_query_proteins_padloc.csv"))
    ),
    DefensePredictor = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_defensepredictor", "defense_predictor_output.csv"))
    ),
    REBASEfinder = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_rebasefinder", "R-M_REBASE_analysis.xlsx"))
    ),
    GapMind = list(
      mode = "all",
      paths = c(
        file.path(wd, "dnmb_module_gapmindaa", "aa.sum.steps"),
        file.path(wd, "dnmb_module_gapmindcarbon", "aa.sum.steps")
      )
    ),
    MEROPS = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_merops", "merops_blastp.tsv"))
    ),
    dbCAN = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_dbcan", "run_dbcan", "overview.tsv"))
    ),
    PAZy = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_pazy", "pazy_merged.tsv"))
    ),
    ISelement = list(
      mode = "all",
      paths = c(file.path(wd, "dnmb_module_iselement", "iselement_elements.tsv"))
    ),
    PhiSpy = list(
      mode = "all",
      paths = c(
        file.path(wd, "dnmb_module_phispy", "prophage_coordinates.tsv"),
        file.path(wd, "dnmb_module_phispy", "prophage_information.tsv")
      )
    ),
    VirSorter2 = list(
      mode = "any",
      paths = c(
        file.path(wd, "dnmb_module_virsorter2", "virsorter2_boundary.tsv"),
        file.path(wd, "dnmb_module_virsorter2", "virsorter2_score.tsv")
      )
    ),
    PIDE = list(
      mode = "any",
      paths = c(
        file.path(wd, "dnmb_module_pide", "pide_cluster.csv"),
        file.path(wd, "dnmb_module_pide", "pide_virus_contig.txt")
      )
    ),
    list(mode = "none", paths = character())
  )
}

.dnmb_module_stage_result_keys <- function(alias) {
  alias <- as.character(alias)[1]
  if (identical(alias, "GapMind")) {
    return(c("GapMindAA", "GapMindCarbon"))
  }
  alias
}

.dnmb_module_stage_artifacts_exist <- function(alias, wd = getwd()) {
  spec <- .dnmb_module_stage_expected_artifacts(alias, wd = wd)
  paths <- spec$paths %||% character()
  if (!length(paths)) {
    return(TRUE)
  }

  exists_vec <- file.exists(paths)
  if (identical(spec$mode, "any")) {
    any(exists_vec)
  } else {
    all(exists_vec)
  }
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

  requested_aliases <- sort(unique(as.character(signature$module_aliases)))
  reusable_aliases <- .dnmb_module_stage_reusable_aliases(cache = cache, signature = signature, wd = wd)
  if (!setequal(requested_aliases, reusable_aliases)) {
    return(list(reusable = FALSE, reason = "missing_artifacts", cache = cache))
  }

  list(reusable = TRUE, reason = "matching_signature", cache = cache)
}

.dnmb_module_stage_reusable_aliases <- function(cache, signature = NULL, wd = getwd()) {
  if (is.null(cache) || is.null(cache$signature) || is.null(cache$module_results) || is.null(signature)) {
    return(character())
  }

  cached_sig <- cache$signature
  if (!.dnmb_same_file_signature(signature$genbank_signature, cached_sig$genbank_signature)) {
    return(character())
  }

  if (!identical(signature$requested, cached_sig$requested)) {
    return(character())
  }

  if (!identical(signature$related_inputs, cached_sig$related_inputs)) {
    return(character())
  }

  requested_aliases <- sort(unique(as.character(signature$module_aliases)))
  cached_aliases <- sort(unique(as.character(cached_sig$module_aliases)))
  overlap <- intersect(requested_aliases, cached_aliases)
  if (!length(overlap)) {
    return(character())
  }

  keep <- vapply(overlap, function(alias) {
    result_keys <- .dnmb_module_stage_result_keys(alias)
    has_result <- all(vapply(result_keys, function(key) {
      key %in% names(cache$module_results) && !is.null(cache$module_results[[key]])
    }, logical(1)))
    same_db <- identical(signature$db_state[[alias]], cached_sig$db_state[[alias]])
    has_artifacts <- .dnmb_module_stage_artifacts_exist(alias, wd = wd)
    isTRUE(has_result) && isTRUE(same_db) && isTRUE(has_artifacts)
  }, logical(1))

  overlap[keep]
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
