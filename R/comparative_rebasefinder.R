#' Comparative REBASEfinder Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, reads each genome's
#' organized GenBank table (\code{*_total.xlsx} sheet \code{1.GenBank_table}),
#' and renders a ComplexHeatmap of R-M system family (Type I/II/III/IV)
#' counts. Uses the same layout as
#' \code{dnmb_plot_comparative_defensefinder()}.
#'
#' The canonical REBASE family classification lives in the
#' \code{REBASEfinder_family_id} column on the organized GenBank table
#' (populated by \code{Output_combiner} after a REBASEfinder run). Each
#' hit row contributes one count in the (genome, family) cell.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_rebasefinder <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_REBASEfinder_Heatmap.pdf",
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

  marker_rel <- file.path("dnmb_module_rebasefinder", "blast_results.txt")

  if (isTRUE(auto_run_missing)) {
    # Readiness for REBASEfinder is the module-native analysis xlsx.
    # blast_results.txt alone isn't enough — we need the post-processed
    # R-M_REBASE_analysis.xlsx that carries the rm_type (Type I/II/III/IV)
    # column the heatmap renders from.
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "REBASEfinder",
      module_marker_rel = marker_rel,
      ready_check = function(dir) {
        !is.na(.dnmb_comparative_rebasefinder_xlsx(dir))
      },
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_rebasefinder(data_root, verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No REBASEfinder family assignments found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "REBASEfinder",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- REBASEfinder-specific collector ----

# Locate the module-native analysis xlsx. REBASEfinder writes it as
# <g$dir>/dnmb_module_rebasefinder/R-M_REBASE_analysis.xlsx, but some
# layouts wrap that under <g$dir>/fullrun_*/. Probe both.
.dnmb_comparative_rebasefinder_xlsx <- function(genome_dir) {
  rel <- file.path("dnmb_module_rebasefinder", "R-M_REBASE_analysis.xlsx")
  direct <- file.path(genome_dir, rel)
  if (file.exists(direct)) return(direct)
  inner <- list.files(genome_dir, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  inner <- inner[utils::file_test("-d", inner)]
  cand <- file.path(inner, rel)
  hit <- which(file.exists(cand))
  if (length(hit)) cand[hit[1]] else NA_character_
}

.dnmb_comparative_collect_rebasefinder <- function(data_root, verbose = TRUE) {
  # Key off the module-native R-M_REBASE_analysis.xlsx (produced by
  # run_module_set(db = "REBASEfinder")), because that file carries
  # the per-locus `rm_type` Type I/II/III/IV classification already.
  # That removes any dependency on Output_combiner having merged a
  # REBASEfinder_family_id column into the organized *_total.xlsx.
  genomes <- .dnmb_comparative_discover_gbff_dirs(data_root)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    # Register every gbff-bearing folder so genomes where REBASEfinder
    # never ran still render as empty rows.
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    xlsx_path <- .dnmb_comparative_rebasefinder_xlsx(g$dir)
    if (is.na(xlsx_path)) {
      if (verbose) message("  ", genome_id, " — 0 R-M hits (not analyzed)")
      next
    }
    tab <- tryCatch(
      openxlsx::read.xlsx(xlsx_path, sheet = "RM_Comprehensive"),
      error = function(e) NULL
    )
    if (is.null(tab) || !all(c("rm_type", "passed_blast_filter") %in% names(tab))) {
      if (verbose) message("  ", genome_id, " — 0 R-M hits (malformed xlsx)")
      next
    }
    pass <- as.logical(tab$passed_blast_filter)
    keep <- !is.na(pass) & pass & !is.na(tab$rm_type) & nzchar(as.character(tab$rm_type))
    if (!any(keep)) {
      if (verbose) message("  ", genome_id, " — 0 R-M hits (analyzed, empty)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = as.character(tab$rm_type)[keep],
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", sum(keep), " R-M hits")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
