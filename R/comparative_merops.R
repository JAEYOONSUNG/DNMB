#' Comparative MEROPS Peptidase-Family Heatmap Across Genomes
#'
#' Walks a directory of per-genome analysis folders, reads each genome's
#' \code{dnmb_module_merops/merops_blastp.tsv} (DIAMOND blastp tabular,
#' no header), picks the top hit per query, extracts the MEROPS peptidase
#' family ID (e.g. \code{C26}, \code{S8}) from the subject description
#' block \code{[C26.A01]#C26#...}, and renders a ComplexHeatmap of
#' per-family hit counts using the same layout as
#' \code{dnmb_plot_comparative_defensefinder()}.
#'
#' @inheritParams dnmb_plot_comparative_defensefinder
#' @return Invisibly, a list with \code{pdf}, \code{xlsx}, \code{matrix},
#'   and \code{genomes}.
#' @export
dnmb_plot_comparative_merops <- function(
    data_root,
    output_dir = NULL,
    output_file = "Comparative_MEROPS_Heatmap.pdf",
    color_palette = c("white", "#B71C1C"),
    bar_color = "#E53935",
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

  marker_rel <- file.path("dnmb_module_merops", "merops_blastp.tsv")

  if (isTRUE(auto_run_missing)) {
    .dnmb_comparative_autorun_module(
      data_root,
      module_db = "MEROPS",
      module_marker_rel = marker_rel,
      verbose = verbose,
      module_cache_root = module_cache_root,
      module_install = module_install,
      module_cpu = module_cpu
    )
  }

  collected <- .dnmb_comparative_collect_merops(data_root, marker_rel = marker_rel,
                                                verbose = verbose)
  if (!nrow(collected$systems)) {
    if (verbose) message("No MEROPS results found under ", data_root)
    return(invisible(NULL))
  }

  .dnmb_comparative_render_heatmap(
    long_df = collected$systems,
    organism_map = collected$organism,
    title = "MEROPS",
    output_dir = output_dir,
    output_file = output_file,
    color_palette = color_palette,
    bar_color = bar_color,
    line_width = line_width,
    line_col = line_col,
    verbose = verbose
  )
}

# ---- MEROPS-specific collector ----

.dnmb_comparative_merops_family <- function(subject_desc) {
  # Subject descriptions look like
  #   "MER0039704 - ywpE/BSU36340 ({Bacillus subtilis}) (...) [C60.A01]#C60A#..."
  # Extract family code (C60 / S8 / M20 / I04 / ...) from the bracketed
  # "[<family>.<subfamily><member>]" token. Fall back to the "#C60A#"
  # block when the bracket is absent, stripping the trailing subfamily
  # letter so C60A rolls into C60.
  s <- as.character(subject_desc)
  fam <- rep(NA_character_, length(s))

  pat_br <- "\\[[A-Z]\\d+\\."
  hit1 <- regexpr(pat_br, s) > 0
  if (any(hit1)) {
    raw <- regmatches(s[hit1], regexpr(pat_br, s[hit1]))
    fam[hit1] <- sub("^\\[", "", sub("\\.$", "", raw))
  }

  pat_hash <- "#[A-Z]\\d+[A-Za-z]?#"
  hit2 <- is.na(fam) & regexpr(pat_hash, s) > 0
  if (any(hit2)) {
    raw <- regmatches(s[hit2], regexpr(pat_hash, s[hit2]))
    fam[hit2] <- sub("[A-Za-z]$", "", gsub("#", "", raw))
  }
  fam
}

.dnmb_comparative_collect_merops <- function(data_root, marker_rel, verbose = TRUE) {
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
      if (verbose) message("  ", genome_id, " — 0 peptidases (not analyzed)")
      next
    }
    tab <- tryCatch(
      utils::read.delim(marker_path, header = FALSE, stringsAsFactors = FALSE,
                        quote = "", comment.char = "", check.names = FALSE,
                        na.strings = c("", "NA")),
      error = function(e) NULL
    )
    if (is.null(tab) || !nrow(tab)) {
      if (verbose) message("  ", genome_id, " — 0 peptidases (analyzed, empty)")
      next
    }
    qcol <- tab[[1]]
    scol <- if (ncol(tab) >= 13L) tab[[13]] else tab[[ncol(tab)]]
    # Top hit per query (tabular BLAST output is already sorted by e-value).
    keep_idx <- !duplicated(qcol)
    fam <- .dnmb_comparative_merops_family(scol[keep_idx])
    fam <- fam[!is.na(fam) & nzchar(fam)]
    if (!length(fam)) {
      if (verbose) message("  ", genome_id, " — 0 peptidases (analyzed, empty)")
      next
    }
    systems_all[[length(systems_all) + 1L]] <- data.frame(
      File = genome_id,
      subtype = fam,
      stringsAsFactors = FALSE
    )
    if (verbose) message("  ", genome_id, " — ", length(fam), " peptidases")
  }
  list(
    systems = if (length(systems_all)) do.call(rbind, systems_all) else data.frame(),
    organism = organism_map
  )
}
