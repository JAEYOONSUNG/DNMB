#' Comparative PhiSpy Prophage Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, collects each genome's
#' \code{dnmb_module_phispy/prophage_coordinates.tsv}, and renders a single
#' ComplexHeatmap of prophage counts bucketed by region size. Uses the same
#' layout as \code{dnmb_plot_comparative_defensefinder()}.
#'
#' PhiSpy does not categorise prophages natively, so size (in bp) is used as
#' the subtype axis: Small (<10 kb), Medium (10-30 kb), Large (30-60 kb),
#' Very Large (>=60 kb).
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_phispy <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_PhiSpy_Heatmap.pdf",
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

  marker_rel <- file.path("dnmb_module_phispy", "prophage_coordinates.tsv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "PhiSpy",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_phispy(data_root, marker_rel = marker_rel,
                                                verbose = verbose)

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "PhiSpy",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

.dnmb_comparative_collect_phispy <- function(data_root, marker_rel, verbose = TRUE) {
  genomes <- .dnmb_comparative_discover_gbff_dirs(data_root)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    marker_path <- file.path(g$dir, marker_rel)
    if (!file.exists(marker_path)) {
      if (verbose) message("  ", genome_id, " - 0 prophages (not analyzed)")
      next
    }
    coord <- tryCatch(
      utils::read.delim(marker_path, header = TRUE, stringsAsFactors = FALSE,
                        na.strings = c("", "NA"), check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(coord) || !nrow(coord) || ncol(coord) < 4L) {
      if (verbose) message("  ", genome_id, " - 0 prophages (analyzed, empty)")
      next
    }
    start_bp <- suppressWarnings(as.numeric(coord[[3]]))
    end_bp <- suppressWarnings(as.numeric(coord[[4]]))
    length_bp <- abs(end_bp - start_bp) + 1
    length_bp <- length_bp[is.finite(length_bp) & length_bp > 0]
    if (!length(length_bp)) {
      if (verbose) message("  ", genome_id, " - 0 prophages (analyzed, empty)")
      next
    }
    subtype <- .dnmb_comparative_prophage_size_bucket(length_bp)
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = subtype,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " - ", length(length_bp), " prophages")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}

.dnmb_comparative_prophage_size_bucket <- function(length_bp) {
  levels_ord <- c("Small (<10kb)", "Medium (10-30kb)", "Large (30-60kb)", "Very Large (>=60kb)")
  bucket <- ifelse(length_bp < 10000L, levels_ord[1],
            ifelse(length_bp < 30000L, levels_ord[2],
            ifelse(length_bp < 60000L, levels_ord[3], levels_ord[4])))
  factor(bucket, levels = levels_ord)
}
