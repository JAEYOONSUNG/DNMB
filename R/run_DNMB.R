#' Result DNMB table
#'
#' @param module_dbCAN Logical, whether to run and append dbCAN module results.
#' @param module_MEROPS Logical, whether to run and append MEROPS module results.
#' @param module_CLEAN Logical, whether to run and append CLEAN module results.
#' @param module_PAZy Logical, whether to run and append PAZy module results.
#' @param module_GapMind Logical, whether to run and append GapMind module results.
#' @param module_DefenseFinder Logical, whether to run and append DefenseFinder module results.
#' @param module_PADLOC Logical, whether to run and append PADLOC module results.
#' @param module_DefensePredictor Logical, whether to run and append DefensePredictor module results.
#' @param module_ISelement Logical, whether to run and append IS element module results.
#' @param module_Prophage Logical, whether to run and append prophage module results.
#' @param module_EggNOG Logical, whether to run and append eggnog-mapper module results.
#' @param module_InterProScan Logical, whether to run InterProScan on protein
#'   sequences. Requires InterProScan installed (Linux/Docker only).
#' @param module_Prophage_backend Character backend for prophage detection.
#' @param module_version Optional module version string passed to
#'   `run_module_set()`.
#' @param module_cache_root Optional cache root override passed to
#'   `run_module_set()`.
#' @param module_install Logical; install module assets automatically when
#'   missing.
#' @param module_base_url Optional module asset base URL override.
#' @param module_asset_urls Optional named module asset overrides.
#' @param module_cpu Integer thread count for module external tools.
#' @param module_results Optional precomputed module runs, such as the output
#'   from `run_module_set()`.
#' @param interproscan_applications Character vector of InterProScan analyses
#'   (e.g., \code{c("Pfam", "TIGRFAM")}). \code{NULL} runs all.
#' @param interproscan_path Path to \code{interproscan.sh}. \code{NULL}
#'   auto-detects via \code{INTERPROSCAN_HOME} or PATH.
#' @param clean_previous Logical. If \code{TRUE}, remove prior run artifacts
#'   before starting. Cached module results and InterProScan outputs are kept
#'   when their saved input/database signatures still match the current run.
#' @return Invisibly returns the final `genbank_table`.
#' @export

.dnmb_default_cpu <- function() {
  cores <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA)
  if (is.na(cores) || cores < 1L) cores <- 1L
  as.integer(cores)
}

.dnmb_run_default_prophage_backend <- function() {
  backend_fn <- get0(".dnmb_prophage_default_backend", mode = "function", inherits = TRUE)
  if (is.function(backend_fn)) {
    return(backend_fn())
  }
  "phispy"
}

.dnmb_attach_runtime_packages <- function() {
  pkgs <- c(
    "tidyverse",
    "reshape2",
    "grid",
    "circlize",
    "ComplexHeatmap",
    "Biostrings",
    "openxlsx",
    "seqinr",
    "Peptides",
    "ggplotify",
    "venneuler"
  )

  missing <- character()
  for (pkg in pkgs) {
    ok <- suppressPackageStartupMessages(
      require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
    if (!isTRUE(ok)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing)) {
    stop(
      "DNMB runtime packages are missing: ",
      paste(unique(missing), collapse = ", "),
      ". Rebuild the Docker image or reinstall the package dependencies.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

run_DNMB <- function(
    module_dbCAN = TRUE,
    module_MEROPS = TRUE,
    module_CLEAN = TRUE,
    module_PAZy = TRUE,
    module_GapMind = TRUE,
    module_DefenseFinder = TRUE,
    module_PADLOC = TRUE,
    module_DefensePredictor = TRUE,
    module_REBASEfinder = TRUE,
    module_ISelement = TRUE,
    module_Prophage = TRUE,
    module_EggNOG = TRUE,
    module_InterProScan = TRUE,
    module_Prophage_backend = .dnmb_run_default_prophage_backend(),
    module_version = NULL,
    module_cache_root = NULL,
    module_install = TRUE,
    module_base_url = NULL,
    module_asset_urls = NULL,
    module_cpu = .dnmb_default_cpu(),
    module_results = NULL,
    iselement_analysis_depth = "full",
    iselement_related_genbanks = NULL,
    iselement_related_metadata = NULL,
    iselement_auto_discover_related = TRUE,
    iselement_max_related = 5L,
    interproscan_applications = NULL,
    interproscan_path = NULL,
    clean_previous = TRUE
) {
  .dnmb_attach_runtime_packages()

  genbank_signature <- .dnmb_genbank_input_signature(getwd())
  eggnog_external_status <- .dnmb_eggnog_reuse_status(
    wd = getwd(),
    genbank_signature = genbank_signature,
    allow_external_without_metadata = TRUE
  )
  requested_module_EggNOG <- isTRUE(module_EggNOG) && !isTRUE(eggnog_external_status$reusable)
  requested_module_aliases <- dnmb_enabled_module_aliases(
    module_dbCAN = module_dbCAN,
    module_MEROPS = module_MEROPS,
    module_CLEAN = module_CLEAN,
    module_PAZy = module_PAZy,
    module_GapMind = module_GapMind,
    module_DefenseFinder = module_DefenseFinder,
    module_PADLOC = module_PADLOC,
    module_DefensePredictor = module_DefensePredictor,
    module_REBASEfinder = module_REBASEfinder,
    module_ISelement = module_ISelement,
    module_Prophage = module_Prophage,
    module_EggNOG = requested_module_EggNOG,
    module_Prophage_backend = module_Prophage_backend
  )
  requested_module_stage_signature <- .dnmb_module_stage_signature(
    genbank_signature = genbank_signature,
    module_aliases = requested_module_aliases,
    module_version = module_version,
    module_cache_root = module_cache_root,
    module_install = module_install,
    module_base_url = module_base_url,
    module_asset_urls = module_asset_urls,
    module_Prophage_backend = module_Prophage_backend,
    iselement_analysis_depth = iselement_analysis_depth,
    iselement_related_genbanks = iselement_related_genbanks,
    iselement_related_metadata = iselement_related_metadata,
    iselement_auto_discover_related = iselement_auto_discover_related,
    iselement_max_related = iselement_max_related
  )

  if (isTRUE(clean_previous)) {
    message("[DNMB] Cleaning previous run artifacts...")
    dnmb_clean_previous_run(
      genbank_signature = genbank_signature,
      module_stage_signature = requested_module_stage_signature
    )
  }

  Genbank_organizer()
  message("Step1. genbank_organizer function executed")

  Codon_usage_calculator()
  message("Step2. Codon usage caculator function executed")

  tRNA_anticodon_counter()
  message("Step3. tRNA anticodon count function executed")

  genbank_fna_extractor()
  message("Step4. Genbank fna extractor function executed")

  RBS_extractor(plot_RBS = TRUE)
  message("Step5. RBS extractor function executed")

  CRISPR_array_extractor()
  message("Step6. CRISPR array extractor function executed")

  codon_usage_tRNA_heatmap_generator()
  Figure_generator()
  message("Step7. plot function executed")

  # ── EggNOG: auto-detect external results, skip module if found ──
  eggnog_external_found <- FALSE
  if (isTRUE(eggnog_external_status$reusable)) {
    EggNOG_annotations()
    .dnmb_write_eggnog_signature(getwd(), genbank_signature = genbank_signature)
    message("Step8. EggNOG external results detected and parsed.")
    eggnog_external_found <- TRUE

    EggNOG_table <- get0("EggNOG_table", envir = .GlobalEnv, inherits = FALSE)
    genbank_table <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)

    if (!is.null(EggNOG_table) && !is.null(genbank_table) &&
        "locus_tag" %in% colnames(genbank_table) && "query" %in% colnames(EggNOG_table)) {
      genbank_table <- merge(genbank_table, EggNOG_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)
      assign("genbank_table", genbank_table, envir = .GlobalEnv)
    }
  }

  # ── InterProScan: auto-detect external results OR run module ──
  if (exists("InterProScan_table", envir = .GlobalEnv, inherits = FALSE)) {
    rm("InterProScan_table", envir = .GlobalEnv)
  }
  if (exists("InterProScan_site", envir = .GlobalEnv, inherits = FALSE)) {
    rm("InterProScan_site", envir = .GlobalEnv)
  }
  interpro_external_found <- FALSE
  interpro_status <- .dnmb_interproscan_reuse_status(
    wd = getwd(),
    genbank_signature = genbank_signature,
    allow_external_without_metadata = TRUE
  )

  if (isTRUE(interpro_status$reusable)) {
    InterProScan_annotations(InterProScan_dir = interpro_status$source_dir)
    message("Step9. InterProScan existing results detected and parsed.")
    interpro_external_found <- TRUE

    InterProScan_table <- get0("InterProScan_table", envir = .GlobalEnv, inherits = FALSE)
    genbank_table <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)

    if (!is.null(InterProScan_table) && !is.null(genbank_table) &&
        "locus_tag" %in% colnames(genbank_table) && "query" %in% colnames(InterProScan_table)) {
      genbank_table <- merge(genbank_table, InterProScan_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)
      assign("genbank_table", genbank_table, envir = .GlobalEnv)
    }
  } else if (isTRUE(module_InterProScan)) {
    if (isTRUE(interpro_status$has_metadata) && identical(interpro_status$reason, "input_changed")) {
      message("[DNMB] Stored InterProScan results do not match the current GenBank input; rerunning.")
    }

    message("Step9. Running InterProScan...")
    tryCatch({
      ipr_output_dir <- .dnmb_interproscan_output_dir(getwd())
      ipr_result <- run_interproscan(
        output_dir = ipr_output_dir,
        applications = interproscan_applications,
        cpu = module_cpu,
        interproscan_path = interproscan_path,
        verbose = TRUE
      )
      .dnmb_write_interproscan_signature(
        output_dir = ipr_output_dir,
        genbank_signature = genbank_signature
      )
      # Results stay inside the module directory; parse in place rather
      # than duplicating the tsv (and optional .tsv.sites) to the working
      # directory, which would double disk usage for no consumer benefit.
      InterProScan_annotations(InterProScan_dir = dirname(ipr_result$tsv))
      interpro_external_found <- TRUE

      InterProScan_table <- get0("InterProScan_table", envir = .GlobalEnv, inherits = FALSE)
      genbank_table <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)
      if (!is.null(InterProScan_table) && !is.null(genbank_table) &&
          "locus_tag" %in% colnames(genbank_table) && "query" %in% colnames(InterProScan_table)) {
        genbank_table <- merge(genbank_table, InterProScan_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)
        assign("genbank_table", genbank_table, envir = .GlobalEnv)
      }
      message("Step9. InterProScan completed and merged.")
    }, error = function(e) {
      message("Step9. InterProScan failed: ", conditionMessage(e))
    })
  }

  # Skip EggNOG module if external results were already parsed
  if (eggnog_external_found && isTRUE(module_EggNOG)) {
    message("[DNMB] EggNOG external results already loaded — skipping module execution.")
    module_EggNOG <- FALSE
  }

  # Save base results before modules (in case modules take long)
  DNMB_table(genbank_table = get("genbank_table", envir = .GlobalEnv))
  message("Step9b. Base DNMB table saved.")

  module_aliases <- dnmb_enabled_module_aliases(
    module_dbCAN = module_dbCAN,
    module_MEROPS = module_MEROPS,
    module_CLEAN = module_CLEAN,
    module_PAZy = module_PAZy,
    module_GapMind = module_GapMind,
      module_DefenseFinder = module_DefenseFinder,
      module_PADLOC = module_PADLOC,
      module_DefensePredictor = module_DefensePredictor,
      module_REBASEfinder = module_REBASEfinder,
      module_ISelement = module_ISelement,
      module_Prophage = module_Prophage,
      module_EggNOG = module_EggNOG,
      module_Prophage_backend = module_Prophage_backend
    )
  module_stage_signature <- .dnmb_module_stage_signature(
    genbank_signature = genbank_signature,
    module_aliases = module_aliases,
    module_version = module_version,
    module_cache_root = module_cache_root,
    module_install = module_install,
    module_base_url = module_base_url,
    module_asset_urls = module_asset_urls,
    module_Prophage_backend = module_Prophage_backend,
    iselement_analysis_depth = iselement_analysis_depth,
    iselement_related_genbanks = iselement_related_genbanks,
    iselement_related_metadata = iselement_related_metadata,
    iselement_auto_discover_related = iselement_auto_discover_related,
    iselement_max_related = iselement_max_related
  )

  if (length(module_aliases) > 0L || !is.null(module_results)) {
    genbank_table <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)

    if (is.null(genbank_table)) {
      stop("genbank_table not found in the environment after Genbank_organizer().", call. = FALSE)
    }

    module_stage <- tryCatch({
      resolved_module_results <- module_results

      if (is.null(resolved_module_results) && length(module_aliases) > 0L) {
        module_cache_status <- .dnmb_module_stage_cache_status(
          wd = getwd(),
          signature = module_stage_signature
        )

        if (isTRUE(module_cache_status$reusable)) {
          message("[DNMB] Reusing cached module outputs for matching GenBank and DB state.")
          resolved_module_results <- module_cache_status$cache$module_results
        } else {
          cached_results <- module_cache_status$cache$module_results %||% list()
          reusable_aliases <- .dnmb_module_stage_reusable_aliases(
            cache = module_cache_status$cache,
            signature = module_stage_signature
          )
          reused_results <- if (length(reusable_aliases)) cached_results[reusable_aliases] else list()
          missing_aliases <- setdiff(module_aliases, reusable_aliases)

          if (length(reusable_aliases)) {
            message("[DNMB] Reusing existing module outputs for: ", paste(reusable_aliases, collapse = ", "))
          }
          if (length(missing_aliases)) {
            message("[DNMB] Running remaining modules: ", paste(missing_aliases, collapse = ", "))
          }

          newly_resolved_results <- if (length(missing_aliases)) {
            dnmb_resolve_module_results(
              module_aliases = missing_aliases,
              module_results = NULL,
              module_Prophage_backend = module_Prophage_backend,
              module_version = module_version,
              module_cache_root = module_cache_root,
              module_install = module_install,
              module_base_url = module_base_url,
              module_asset_urls = module_asset_urls,
              module_cpu = module_cpu,
              iselement_analysis_depth = iselement_analysis_depth,
              iselement_related_genbanks = iselement_related_genbanks,
              iselement_related_metadata = iselement_related_metadata,
              iselement_auto_discover_related = iselement_auto_discover_related,
              iselement_max_related = iselement_max_related
            )
          } else {
            list()
          }

          resolved_module_results <- c(reused_results, newly_resolved_results)
          .dnmb_write_module_stage_cache(
            wd = getwd(),
            signature = module_stage_signature,
            module_results = resolved_module_results
          )
        }
      }

      genbank_table <- append_module_results(genbank_table, resolved_module_results)
      assign("genbank_table", genbank_table, envir = .GlobalEnv)
      message("Step10. Module outputs appended onto genbank_table.")

      render_module_plots_fn <- get0("dnmb_render_module_plots", mode = "function", inherits = TRUE)
      if (!is.null(render_module_plots_fn)) {
        module_plots <- tryCatch(
          render_module_plots_fn(genbank_table = genbank_table, output_dir = getwd(), cache_root = module_cache_root),
          error = function(e) e
        )
        if (inherits(module_plots, "error")) {
          message("Step11. Module plots skipped: ", conditionMessage(module_plots))
        } else if (length(module_plots) > 0L) {
          message("Step11. Module plots generated in ", file.path(getwd(), "visualizations"))
        }
      }
      TRUE
    }, error = function(e) e)
    if (inherits(module_stage, "error")) {
      message("Step10. Module stage skipped: ", conditionMessage(module_stage))
    }
  }

  DNMB_table(genbank_table = get("genbank_table", envir = .GlobalEnv))
  message("Analysis end: DNMB function executed")

  invisible(get("genbank_table", envir = .GlobalEnv))
}

dnmb_enabled_module_aliases <- function(module_dbCAN = FALSE,
                                        module_MEROPS = FALSE,
                                        module_CLEAN = FALSE,
                                        module_PAZy = FALSE,
                                        module_GapMind = FALSE,
                                        module_DefenseFinder = FALSE,
                                        module_PADLOC = FALSE,
                                        module_DefensePredictor = FALSE,
                                        module_REBASEfinder = FALSE,
                                        module_ISelement = FALSE,
                                        module_Prophage = FALSE,
                                        module_EggNOG = FALSE,
                                        module_Prophage_backend = .dnmb_run_default_prophage_backend()) {
  module_flags <- c(
    EggNOG = isTRUE(module_EggNOG),
    CLEAN = isTRUE(module_CLEAN),
    DefenseFinder = isTRUE(module_DefenseFinder),
    PADLOC = isTRUE(module_PADLOC),
    DefensePredictor = isTRUE(module_DefensePredictor),
    REBASEfinder = isTRUE(module_REBASEfinder),
    GapMind = isTRUE(module_GapMind),
    MEROPS = isTRUE(module_MEROPS),
    dbCAN = isTRUE(module_dbCAN),
    PAZy = isTRUE(module_PAZy),
    ISelement = isTRUE(module_ISelement),
    Prophage = isTRUE(module_Prophage)
  )

  names(module_flags)[module_flags]
}

dnmb_resolve_module_results <- function(module_aliases,
                                        module_results = NULL,
                                        module_Prophage_backend = .dnmb_run_default_prophage_backend(),
                                        module_version = NULL,
                                        module_cache_root = NULL,
                                        module_install = TRUE,
                                        module_base_url = NULL,
                                        module_asset_urls = NULL,
                                        module_cpu = .dnmb_default_cpu(),
                                        iselement_analysis_depth = "full",
                                        iselement_related_genbanks = NULL,
                                        iselement_related_metadata = NULL,
                                        iselement_auto_discover_related = TRUE,
                                        iselement_max_related = 5L) {
  if (!is.null(module_results)) {
    return(module_results)
  }

  if (length(module_aliases) == 0L) {
    return(NULL)
  }

  run_module_set_fn <- get0("run_module_set", mode = "function", inherits = TRUE)

  if (is.null(run_module_set_fn)) {
    stop(
      "run_module_set() is not available but module aliases were enabled.",
      call. = FALSE
    )
  }

  all_module_flags <- c(
    "dbCAN",
    "MEROPS",
    "CLEAN",
    "PAZy",
    "GapMind",
    "DefenseFinder",
    "PADLOC",
    "DefensePredictor",
    "REBASEfinder",
    "ISelement",
    "Prophage",
    "EggNOG"
  )
  module_flags <- as.list(stats::setNames(rep(FALSE, length(all_module_flags)), paste0("module_", all_module_flags)))
  if (length(module_aliases) > 0L) {
    module_flags[paste0("module_", module_aliases)] <- TRUE
  }
  module_flags$module_version <- module_version
  module_flags$module_cache_root <- module_cache_root
  module_flags$module_install <- module_install
  module_flags$module_base_url <- module_base_url
  module_flags$module_asset_urls <- module_asset_urls
  module_flags$module_cpu <- module_cpu
  module_flags$module_Prophage_backend <- module_Prophage_backend
  module_flags$iselement_analysis_depth <- iselement_analysis_depth
  module_flags$iselement_related_genbanks <- iselement_related_genbanks
  module_flags$iselement_related_metadata <- iselement_related_metadata
  module_flags$iselement_auto_discover_related <- iselement_auto_discover_related
  module_flags$iselement_max_related <- iselement_max_related
  do.call(run_module_set_fn, module_flags)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Clean previous run artifacts
# ═══════════════════════════════════════════════════════════════════════════════

#' Remove previous DNMB run artifacts from the working directory
#'
#' Deletes intermediate files produced by a previous \code{run_DNMB()} call so
#' the next analysis starts from a clean slate.
#' GenBank input files (\code{*.gbff}, \code{*.gb}, \code{*.gbk}) and the
#' DB cache directory (\code{/opt/dnmb/cache} or \code{~/.dnmb-cache}) are
#' never touched.
#'
#' @param dry_run Logical. If \code{TRUE}, only prints what would be deleted
#'   without actually removing anything.
#' @param keep_excel Logical. If \code{TRUE} (the default), files matching
#'   \code{*_total.xlsx} are preserved.
#' @param genbank_signature Optional current GenBank fingerprint used to decide
#'   whether prior cached outputs still match the active input.
#' @param module_stage_signature Optional module-stage signature used to decide
#'   whether previous module outputs can be preserved.
#' @return Invisibly returns a character vector of paths that were (or would be)
#'   deleted.
#' @export
dnmb_clean_previous_run <- function(dry_run = FALSE,
                                    keep_excel = TRUE,
                                    genbank_signature = NULL,
                                    module_stage_signature = NULL) {
  wd <- getwd()
  deleted <- character(0)

  # --- helper: collect paths that exist ---
  collect <- function(paths) {
    paths[file.exists(paths)]
  }

  module_cache_status <- .dnmb_module_stage_cache_status(
    wd = wd,
    signature = module_stage_signature
  )
  interpro_status <- .dnmb_interproscan_reuse_status(
    wd = wd,
    genbank_signature = genbank_signature,
    allow_external_without_metadata = FALSE
  )
  eggnog_status <- .dnmb_eggnog_reuse_status(
    wd = wd,
    genbank_signature = genbank_signature,
    allow_external_without_metadata = FALSE
  )

  # 1. dnmb_module_* directories
  module_dirs <- list.dirs(wd, full.names = TRUE, recursive = FALSE)
  module_dirs <- module_dirs[grepl("^dnmb_module_", basename(module_dirs))]
  if (isTRUE(module_cache_status$reusable)) {
    message("[DNMB] Preserving module outputs for matching GenBank and DB state.")
  } else {
    if (length(module_dirs)) {
      message("[DNMB] Existing module output folders found but stage cache is stale or missing; keeping folders and allowing rerun to overwrite/update as needed.")
    }
  }

  # 2. InterProScan output directory
  interpro_state <- .dnmb_interproscan_state(wd)
  if (isTRUE(interpro_status$reusable)) {
    message("[DNMB] Preserving InterProScan outputs for matching GenBank input.")
  } else {
    for (cand in unique(c(interpro_state$module_output_dir, interpro_state$legacy_output_dir))) {
      if (dir.exists(cand)) {
        deleted <- c(deleted, cand)
      }
    }
    deleted <- c(
      deleted,
      collect(c(
        interpro_state$root_tsv,
        interpro_state$root_tsv_sites
      ))
    )
    if (isTRUE(interpro_status$has_metadata) && identical(interpro_status$reason, "input_changed")) {
      message("[DNMB] Removing stale InterProScan outputs because the GenBank input changed.")
    }
  }

  # 3. visualizations/ directory
  viz_dir <- file.path(wd, "visualizations")
  if (dir.exists(viz_dir)) deleted <- c(deleted, viz_dir)

  # 4. eggnog-mapper annotations
  if (isTRUE(eggnog_status$reusable)) {
    message("[DNMB] Preserving EggNOG external results for matching GenBank input.")
  } else {
    deleted <- c(deleted, eggnog_status$result_files)
    deleted <- c(deleted, collect(eggnog_status$metadata_path))
    if (isTRUE(eggnog_status$has_metadata) && identical(eggnog_status$reason, "input_changed")) {
      message("[DNMB] Removing stale EggNOG external results because the GenBank input changed.")
    }
  }

  # 5. CDS FASTA files
  cds_fasta <- list.files(wd, pattern = "_CDS_fasta\\.faa$",
                          full.names = TRUE)
  deleted <- c(deleted, cds_fasta)

  # 6. Optionally keep *_total.xlsx
  if (!isTRUE(keep_excel)) {
    xlsx_files <- list.files(wd, pattern = "_total\\.xlsx$",
                             full.names = TRUE)
    deleted <- c(deleted, xlsx_files)
  }

  deleted <- unique(deleted)

  if (length(deleted) == 0L) {
    message("[DNMB] No previous run artifacts found.")
    return(invisible(character(0)))
  }

  if (isTRUE(dry_run)) {
    message("[DNMB] Dry run -- the following would be deleted:")
    for (p in deleted) message("  ", p)
    return(invisible(deleted))
  }

  # Actually remove
  n_removed <- 0L
  for (p in deleted) {
    ok <- tryCatch({
      if (dir.exists(p)) {
        unlink(p, recursive = TRUE, force = TRUE)
      } else {
        file.remove(p)
      }
      TRUE
    }, error = function(e) {
      message("[DNMB] Could not remove: ", p, " (", conditionMessage(e), ")")
      FALSE
    })
    if (isTRUE(ok)) n_removed <- n_removed + 1L
  }

  message("[DNMB] Cleaned ", n_removed, "/", length(deleted), " artifacts.")
  invisible(deleted)
}
