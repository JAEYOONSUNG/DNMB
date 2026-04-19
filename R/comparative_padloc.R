#' Comparative PADLOC Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, collects each genome's
#' \code{dnmb_module_padloc/padloc_query_proteins_padloc.csv}, and renders a
#' single ComplexHeatmap of PADLOC defense-system counts. Uses the same
#' layout as \code{dnmb_plot_comparative_defensefinder()}.
#'
#' Each PADLOC hit row carries a \code{system} (system family/subtype) and a
#' \code{system.number} (system instance on that genome). One occurrence is
#' counted per unique (genome, system, system.number) triple so multi-gene
#' systems are tallied once.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_padloc <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_PADLOC_Heatmap.pdf",
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

  marker_rel <- file.path("dnmb_module_padloc", "padloc_query_proteins_padloc.csv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "PADLOC",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_padloc(data_root, marker_rel = marker_rel,
                                                verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No PADLOC results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "PADLOC",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- PADLOC-specific collector ----

.dnmb_comparative_collect_padloc <- function(data_root, marker_rel, verbose = TRUE) {
  # Use gbff dirs as the universe so genomes PADLOC never ran on still
  # render as empty rows. Only genomes with the marker file contribute
  # actual hits.
  genomes <- .dnmb_comparative_discover_gbff_dirs(data_root)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    # PADLOC output may sit directly under g$dir or one level deeper
    # (e.g. g$dir/fullrun_*/). Probe both layouts.
    marker_path <- file.path(g$dir, marker_rel)
    if (!file.exists(marker_path)) {
      inner <- list.files(g$dir, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
      inner <- inner[utils::file_test("-d", inner)]
      candidates <- file.path(inner, marker_rel)
      hit_idx <- which(file.exists(candidates))
      marker_path <- if (length(hit_idx)) candidates[hit_idx[1]] else NA_character_
    }
    if (is.na(marker_path) || !file.exists(marker_path)) {
      if (verbose) message("  ", genome_id, " — 0 systems (not analyzed)")
      next
    }
    hits <- tryCatch(
      utils::read.csv(marker_path, stringsAsFactors = FALSE,
                      na.strings = c("", "NA"), check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(hits)) next
    if (!nrow(hits) || !all(c("system", "system.number") %in% names(hits))) {
      if (verbose) message("  ", genome_id, " — 0 systems (analyzed, empty)")
      next
    }
    key <- paste(hits$system, hits$system.number, sep = "::")
    systems_unique <- hits[!duplicated(key), c("system", "system.number"), drop = FALSE]
    systems_unique$subtype <- as.character(systems_unique$system)
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
