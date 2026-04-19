#' Comparative DefenseFinder Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, collects each genome's
#' \code{dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv}, and
#' renders a single ComplexHeatmap of defense system counts. Subtype names
#' come from the canonical \code{model_fqn} path so they are stable across
#' runs regardless of how \code{sys_id} prefixes happen to be composed.
#'
#' @param data_root Directory containing one subfolder per genome (each
#'   subfolder may hold DefenseFinder results directly, or one level deeper
#'   under a wrapper such as \code{fullrun_*/}).
#' @param output_dir Where to save the PDF. Defaults to
#'   \code{file.path(data_root, "comparative")}.
#' @param output_file PDF file name.
#' @param color_palette Two-color gradient for heatmap dots.
#' @param bar_color Diversity barplot fill color.
#' @param line_width,line_col Grid line styling.
#' @param auto_run_missing If TRUE, folders with a GenBank record but no
#'   DefenseFinder results are run through \code{run_module_set(db = "DefenseFinder")}.
#' @param module_cache_root,module_install,module_cpu Passed through to
#'   \code{run_module_set()} when auto-running.
#' @param verbose Print progress.
#'
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_defensefinder <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_DefenseFinder_Heatmap.pdf",
    color_palette = c("white", "#330066"),
    bar_color = "#4C1C7E",
    line_width = 1,
    line_col = "grey80",
    auto_run_missing = TRUE,
    module_cache_root = NULL,
    module_install = TRUE,
    module_cpu = NULL,
    verbose = TRUE
) {
  stopifnot(is.character(data_root), length(data_root) == 1L, dir.exists(data_root))
  if (is.null(output_dir)) output_dir <- file.path(data_root, "comparative")

  marker_rel <- file.path("dnmb_module_defensefinder", "defensefinder_best_solution_genes.tsv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "DefenseFinder",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_defensefinder(data_root, marker_rel = marker_rel,
                                                       verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No DefenseFinder results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "DefenseFinder",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- DefenseFinder-specific helpers ----

.dnmb_comparative_subtype_from_model_fqn <- function(model_fqn) {
  # model_fqn looks like "defense-finder-models/DefenseFinder/BREX/BREX_III"
  # or "defense-finder-models/RM/RM/RM_Type_III" — the last path component
  # is the authoritative subtype name, used regardless of whatever sys_id
  # prefix scheme a run happens to use.
  #
  # DefenseFinder upstream occasionally ships sibling model variants (e.g.
  # RM_Type_IV and RM_Type_IV_1) that cover the same biological subtype via
  # different HMM profile families. Collapse the trailing Arabic-digit
  # suffix so variants roll up into one column; Roman numerals
  # (RM_Type_I, Wadjet_II, …) are unaffected.
  name <- basename(as.character(model_fqn))
  sub("_\\d+$", "", name)
}

.dnmb_comparative_collect_defensefinder <- function(data_root, marker_rel, verbose = TRUE) {
  # We key off best_solution_genes.tsv because it carries model_fqn — the
  # canonical subtype path — which is stable across runs. The sys_id column
  # alone is unreliable: some runs prepend "DFREP###_", others emit a bare
  # leading underscore, and regex surgery over sys_id produces ghost
  # columns like "_RM_Type_I".
  genomes <- .dnmb_comparative_discover_genome_dirs(data_root, marker = marker_rel)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    genes_tbl <- tryCatch(
      utils::read.delim(g$marker_path, header = TRUE, stringsAsFactors = FALSE,
                        na.strings = c("", "NA"), check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(genes_tbl)) next
    # Always register the genome in organism_map: the marker file exists
    # so the module has been run; a zero-hit result must render as an
    # empty row rather than being filtered out.
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    if (!nrow(genes_tbl) || !all(c("sys_id", "model_fqn") %in% names(genes_tbl))) {
      if (verbose) message("  ", genome_id, " — 0 systems (analyzed, empty)")
      next
    }
    systems_unique <- genes_tbl[!duplicated(genes_tbl$sys_id),
                                c("sys_id", "model_fqn"), drop = FALSE]
    systems_unique$subtype <- .dnmb_comparative_subtype_from_model_fqn(systems_unique$model_fqn)
    systems_unique <- systems_unique[nzchar(systems_unique$subtype) &
                                     !is.na(systems_unique$subtype), , drop = FALSE]
    if (!nrow(systems_unique)) {
      if (verbose) message("  ", genome_id, " — 0 systems (analyzed, empty)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = systems_unique$subtype,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", nrow(systems_unique), " systems")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
