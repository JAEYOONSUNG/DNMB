#' Comparative CGC (CAZyme Gene Cluster) Signature Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, reads each genome's
#' \code{dnmb_module_dbcan/run_dbcan/cgc_standard_out_summary.tsv}, builds
#' a compound signature label for every CGC from which signature classes
#' are present (\code{CAZyme}, \code{TC}, \code{TF}, \code{STP},
#' \code{Sulfatase}, \code{Peptidase}), and renders a ComplexHeatmap of
#' signature-mix counts. CAZyme is always present by definition, so the
#' label is always led by \code{CAZyme}.
#'
#' Uses the same layout as \code{dnmb_plot_comparative_defensefinder()}.
#' Shares the dbCAN module auto-run with
#' \code{dnmb_plot_comparative_dbcan()}.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_cgc <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_CGC_Heatmap.pdf",
    color_palette = c("white", "#006064"),
    bar_color = "#0097A7",
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

  marker_rel <- file.path("dnmb_module_dbcan", "run_dbcan",
                          "cgc_standard_out_summary.tsv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "dbCAN",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_cgc(data_root, marker_rel = marker_rel,
                                             verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No CGC results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "CAZyme Gene Clusters (CGC)",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- CGC-specific collector ----

.dnmb_comparative_cgc_label <- function(row) {
  # Build a compound label like "CAZyme+TC+TF" from whichever signature
  # categories are present (> 0) in this CGC row. CAZymes is always > 0
  # for a CGC summary entry so the order is CAZyme, TC, TF, STP,
  # Sulfatase, Peptidase.
  parts <- c()
  if (isTRUE(row$CAZymes > 0))   parts <- c(parts, "CAZyme")
  if (isTRUE(row$TC > 0))        parts <- c(parts, "TC")
  if (isTRUE(row$TF > 0))        parts <- c(parts, "TF")
  if (isTRUE(row$STP > 0))       parts <- c(parts, "STP")
  if (isTRUE(row$Sulfatase > 0)) parts <- c(parts, "Sulfatase")
  if (isTRUE(row$Peptidase > 0)) parts <- c(parts, "Peptidase")
  if (!length(parts)) return(NA_character_)
  paste(parts, collapse = "+")
}

.dnmb_comparative_collect_cgc <- function(data_root, marker_rel, verbose = TRUE) {
  genomes <- .dnmb_comparative_discover_gbff_dirs(data_root)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    marker_path <- file.path(g$dir, marker_rel)
    if (!file.exists(marker_path)) {
      inner <- list.files(g$dir, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
      inner <- inner[utils::file_test("-d", inner)]
      candidates <- file.path(inner, marker_rel)
      hit_idx <- which(file.exists(candidates))
      marker_path <- if (length(hit_idx)) candidates[hit_idx[1]] else NA_character_
    }
    if (is.na(marker_path) || !file.exists(marker_path)) {
      if (verbose) message("  ", genome_id, " — 0 CGCs (not analyzed)")
      next
    }
    tab <- tryCatch(
      utils::read.delim(marker_path, header = TRUE, stringsAsFactors = FALSE,
                        na.strings = c("", "NA"), check.names = FALSE,
                        quote = "", comment.char = ""),
      error = function(e) NULL
    )
    needed <- c("CAZymes", "TC", "TF", "STP", "Sulfatase", "Peptidase")
    if (is.null(tab) || !nrow(tab) || !all(needed %in% names(tab))) {
      if (verbose) message("  ", genome_id, " — 0 CGCs (analyzed, empty)")
      next
    }
    for (col in needed) tab[[col]] <- suppressWarnings(as.integer(tab[[col]]))
    labels <- vapply(seq_len(nrow(tab)),
                     function(i) .dnmb_comparative_cgc_label(tab[i, , drop = FALSE]),
                     character(1))
    labels <- labels[!is.na(labels) & nzchar(labels)]
    if (!length(labels)) {
      if (verbose) message("  ", genome_id, " — 0 CGCs (analyzed, empty)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = labels,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", length(labels), " CGCs")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
