#' Comparative dbCAN CAZy-Class Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, reads each genome's
#' \code{dnmb_module_dbcan/run_dbcan/overview.tsv}, extracts the CAZy
#' class (\code{GH}, \code{GT}, \code{PL}, \code{CE}, \code{AA}, \code{CBM})
#' from the recommended call for every row with at least one tool
#' prediction, and renders a ComplexHeatmap of per-class hit counts.
#'
#' Uses the same layout as \code{dnmb_plot_comparative_defensefinder()}.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_dbcan <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_dbCAN_Heatmap.pdf",
    color_palette = c("white", "#1B5E20"),
    bar_color = "#43A047",
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

  marker_rel <- file.path("dnmb_module_dbcan", "run_dbcan", "overview.tsv")

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

  collected <- .dnmb_comparative_collect_dbcan(data_root, marker_rel = marker_rel,
                                               verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No dbCAN results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "dbCAN",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- dbCAN-specific collector ----

.dnmb_comparative_dbcan_class <- function(call) {
  # Calls look like "GH13_3", "GT2", "CE1(24-257)", "GH18_e66|CBM5_e25",
  # "CBM50_e705|CBM50_e699|GH18_e157". Take the first "|"-separated token,
  # strip any "(...)" / "_..." / "+..." suffix, and retain the leading
  # alphabetic prefix as the CAZy class.
  x <- as.character(call)
  x[is.na(x) | x == "-" | !nzchar(x)] <- NA_character_
  first <- vapply(strsplit(x, "[|+]"), function(v) if (length(v)) v[[1]] else NA_character_,
                  character(1))
  first <- sub("\\(.*$", "", first)
  first <- sub("_.*$", "", first)
  class_ <- sub("^([A-Za-z]+).*$", "\\1", first)
  class_[!nzchar(class_)] <- NA_character_
  class_
}

.dnmb_comparative_collect_dbcan <- function(data_root, marker_rel, verbose = TRUE) {
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
      if (verbose) message("  ", genome_id, " — 0 CAZymes (not analyzed)")
      next
    }
    tab <- tryCatch(
      utils::read.delim(marker_path, header = TRUE, stringsAsFactors = FALSE,
                        na.strings = c("", "NA", "-"), check.names = FALSE,
                        quote = "", comment.char = ""),
      error = function(e) NULL
    )
    if (is.null(tab) || !nrow(tab)) {
      if (verbose) message("  ", genome_id, " — 0 CAZymes (analyzed, empty)")
      next
    }
    # Prefer Recommend Results, then dbCAN_sub, then dbCAN_hmm, then DIAMOND.
    pick <- rep(NA_character_, nrow(tab))
    for (col in c("Recommend Results", "dbCAN_sub", "dbCAN_hmm", "DIAMOND")) {
      if (!col %in% names(tab)) next
      val <- as.character(tab[[col]])
      take <- is.na(pick) & !is.na(val) & nzchar(val) & val != "-"
      pick[take] <- val[take]
    }
    klass <- .dnmb_comparative_dbcan_class(pick)
    klass <- klass[!is.na(klass) & nzchar(klass)]
    if (!length(klass)) {
      if (verbose) message("  ", genome_id, " — 0 CAZymes (analyzed, empty)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = klass,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", length(klass), " CAZymes")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
