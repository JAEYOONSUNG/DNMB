#' Comparative DefensePredictor Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, reads each genome's
#' \code{dnmb_module_defensepredictor/defense_predictor_output.csv}, and
#' renders a ComplexHeatmap of DefensePredictor hit counts grouped by
#' defense-protein category. Uses the same layout as
#' \code{dnmb_plot_comparative_defensefinder()}.
#'
#' DefensePredictor scores every CDS with a continuous \code{mean_log_odds}
#' defense-likeness value; there is no intrinsic system-level subtype. For
#' the comparative view we treat the per-hit product annotation (the
#' \code{name} column — e.g. "type I restriction endonuclease subunit R",
#' "DUF4145 domain-containing protein") as the defense category, after
#' trimming generic suffixes ("domain-containing protein", "family
#' protein"). This puts functionally-related hits into the same column so
#' users can compare which defense flavors each genome carries.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @param threshold Minimum \code{mean_log_odds} for a protein to count as
#'   a DefensePredictor hit. Defaults to the package-wide threshold (4.0).
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_defensepredictor <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_DefensePredictor_Heatmap.pdf",
    color_palette = c("white", "#330066"),
    bar_color = "#4C1C7E",
    line_width = 1,
    line_col = "grey80",
    threshold = .dnmb_defensepredictor_default_threshold(),
    auto_run_missing = TRUE,
    module_cache_root = NULL,
    module_install = TRUE,
    module_cpu = NULL,
    verbose = TRUE
) {
  stopifnot(is.character(data_root), length(data_root) == 1L, dir.exists(data_root))
  if (is.null(output_dir)) output_dir <- file.path(data_root, "comparative")

  marker_rel <- file.path("dnmb_module_defensepredictor", "defense_predictor_output.csv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "DefensePredictor",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_defensepredictor(
    data_root, marker_rel = marker_rel,
    threshold = as.numeric(threshold)[1], verbose = verbose
  )
  if (!length(collected$organism)) {
    if (verbose) message("No DefensePredictor results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "DefensePredictor",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- DefensePredictor-specific helpers ----

# No family rollup — DefensePredictor emits raw NCBI product names,
# and the defense-system boundary can't be inferred from the product
# string alone. Just trim generic scaffolding ("... domain-containing
# protein", "... family protein") and collapse whitespace so variants
# of the same annotation land in one column; leave everything else as
# the curator wrote it.
.dnmb_comparative_defensepredictor_category <- function(name) {
  x <- as.character(name)
  x[is.na(x) | !nzchar(x)] <- "hypothetical protein"
  x <- sub("\\s+domain-containing protein$", "", x, ignore.case = TRUE, perl = TRUE)
  x <- sub("\\s+family protein$", "", x, ignore.case = TRUE, perl = TRUE)
  x <- trimws(gsub("\\s+", " ", x))
  x
}

.dnmb_comparative_collect_defensepredictor <- function(data_root, marker_rel,
                                                       threshold, verbose = TRUE) {
  genomes <- .dnmb_comparative_discover_genome_dirs(data_root, marker = marker_rel)
  systems_all <- list()
  organism_map <- character()
  for (g in genomes) {
    hits <- tryCatch(
      utils::read.csv(g$marker_path, stringsAsFactors = FALSE,
                      na.strings = c("", "NA"), check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(hits)) next
    # Always register: the marker file exists, so DefensePredictor ran.
    # Zero-hit genomes must render as empty rows.
    genome_id <- .dnmb_comparative_genome_id(g$id)
    organism_map[genome_id] <- .dnmb_comparative_parse_organism(g$dir)
    if (!nrow(hits) || !"mean_log_odds" %in% names(hits)) {
      if (verbose) message("  ", genome_id, " — 0 DP hits (analyzed, empty)")
      next
    }
    score <- suppressWarnings(as.numeric(hits$mean_log_odds))
    keep <- !is.na(score) & score >= threshold
    if (!any(keep)) {
      if (verbose) message("  ", genome_id, " — 0 DP hits (below threshold)")
      next
    }
    name_col <- if ("name" %in% names(hits)) hits$name else rep(NA_character_, nrow(hits))
    category <- .dnmb_comparative_defensepredictor_category(name_col[keep])
    if (!length(category) || all(!nzchar(category))) {
      if (verbose) message("  ", genome_id, " — 0 DP hits (unclassified)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = category,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", length(category),
                         " DP hits across ", length(unique(category)), " categories")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
