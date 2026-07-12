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
    output_file = "Comparative_dbCAN_class_Heatmap.pdf",
    color_palette = c("white", "#E65100"),
    bar_color = "#FB8C00",
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
                                               level = "class", verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No dbCAN results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "dbCAN (class)",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

#' Comparative dbCAN CAZy-Family Heatmap Across Genomes
#'
#' Finer-grained counterpart to \code{dnmb_plot_comparative_dbcan()}: each
#' recommended call is kept at family resolution (e.g. \code{GH13},
#' \code{GT2}, \code{CBM50}) instead of being rolled up to its class
#' letter. Subfamily suffixes such as \code{GH13_3} are collapsed so
#' \code{GH13_3} and \code{GH13_5} share a single column.
#'
#' Shares the auto-run logic and overview-reading collector with
#' \code{dnmb_plot_comparative_dbcan()}; only the per-call aggregation
#' differs.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_dbcan_family <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_dbCAN_family_Heatmap.pdf",
    color_palette = c("white", "#E65100"),
    bar_color = "#FB8C00",
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
                                               level = "family", verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No dbCAN results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "dbCAN (family)",
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

.dnmb_comparative_dbcan_token <- function(call, level = c("class", "family")) {
  # A dbCAN call may contain several domains. Count each unique gene-family
  # pair instead of silently keeping only the first token.
  level <- match.arg(level)
  x <- as.character(call)
  token_list <- lapply(x, function(value) {
    families <- .dnmb_dbcan_family_tokens(value)
    if (!length(families)) return(character())
    if (level == "class") {
      unique(sub("^([A-Za-z]+).*$", "\\1", families))
    } else {
      unique(sub("^([A-Za-z]+[0-9]+).*$", "\\1", families))
    }
  })
  unname(unlist(token_list, use.names = FALSE))
}

.dnmb_comparative_collect_dbcan <- function(data_root, marker_rel,
                                            level = c("class", "family"),
                                            verbose = TRUE) {
  level <- match.arg(level)
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
    tab <- tryCatch(dnmb_dbcan_parse_overview(marker_path), error = function(e) NULL)
    if (is.null(tab) || !nrow(tab)) {
      if (verbose) message("  ", genome_id, " — 0 CAZymes (analyzed, empty)")
      next
    }
    tokens <- .dnmb_comparative_dbcan_token(tab$dbcan_all_families, level = level)
    tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
    if (!length(tokens)) {
      if (verbose) message("  ", genome_id, " — 0 CAZymes (analyzed, empty)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = tokens,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", length(tokens), " CAZymes")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
