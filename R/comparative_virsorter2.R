#' Comparative VirSorter2 Prophage Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, collects each genome's
#' \code{dnmb_module_virsorter2/virsorter2_boundary.tsv} (and the companion
#' \code{virsorter2_score.tsv} when available), and renders a single
#' ComplexHeatmap of VirSorter2 calls. The \code{max_score_group} field is
#' used as the subtype axis (dsDNAphage, ssDNA, NCLDV, lavidaviridae, RNA)
#' so related viral taxa cluster into stable columns across runs.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_virsorter2 <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_VirSorter2_Heatmap.pdf",
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

  marker_rel <- file.path("dnmb_module_virsorter2", "virsorter2_boundary.tsv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "VirSorter2",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_virsorter2(data_root, marker_rel = marker_rel,
                                                   verbose = verbose)

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "VirSorter2",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

.dnmb_comparative_collect_virsorter2 <- function(data_root, marker_rel, verbose = TRUE) {
  genomes <- .dnmb_comparative_discover_gbff_dirs(data_root)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    marker_path <- file.path(g$dir, marker_rel)
    if (!file.exists(marker_path)) {
      if (verbose) message("  ", genome_id, " - 0 regions (not analyzed)")
      next
    }
    score_path <- file.path(dirname(marker_path), "virsorter2_score.tsv")
    regions <- tryCatch(
      .dnmb_prophage_standardize_virsorter2_regions(
        .dnmb_prophage_parse_virsorter2_boundary(marker_path),
        if (file.exists(score_path)) .dnmb_prophage_parse_virsorter2_scores(score_path) else data.frame()
      ),
      error = function(e) NULL
    )
    if (is.null(regions) || !nrow(regions)) {
      if (verbose) message("  ", genome_id, " - 0 regions (analyzed, empty)")
      next
    }
    subtype <- as.character(regions$prophage_group)
    subtype[is.na(subtype) | !nzchar(subtype)] <- "unclassified"
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = subtype,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " - ", nrow(regions), " regions")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
