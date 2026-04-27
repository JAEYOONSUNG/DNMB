#' Export a slim per-gene translation efficiency summary
#'
#' Emits an Excel workbook with a single Translation_Efficiency sheet of the
#' interpretable per-gene metrics (tir_score, percentile, codon efficiency,
#' CAI, tAI, expression_band) plus an Interpretation sheet describing how
#' to read each column. A companion focused PDF plots the four primary
#' metric distributions, a 2D scatter, and the top / bottom 20 genes by
#' expression_band-rank.
#'
#' @param results mRNAcal results table from `dnmb_run_mrnacal_module`
#'   (or the merged genbank_table with `mRNAcal_*` prefixes; the function
#'   accepts both forms).
#' @param gene_metadata Optional data.frame with locus_tag + gene/product/
#'   protein_id used to enrich the slim sheet. Falls back to the global
#'   `genbank_table`.
#' @param output_dir Output directory; defaults to
#'   `mrnacal_translation_efficiency_export/` in the current working
#'   directory.
#' @param organism Optional label.
#' @param top_n Number of top / bottom genes to list in the PDF table.
#'
#' @return List with `slim_table`, file paths in `files`.
#'
#' @export
dnmb_export_translation_efficiency <- function(results,
                                               gene_metadata = NULL,
                                               output_dir = NULL,
                                               organism = NULL,
                                               top_n = 20L) {
  results <- .dnmb_validate_load_results(results)
  if (!base::nrow(results)) {
    stop("mRNAcal results are empty.", call. = FALSE)
  }
  if (base::is.null(output_dir)) {
    output_dir <- base::file.path(base::getwd(), "mrnacal_translation_efficiency_export")
  }
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (base::is.null(gene_metadata)) {
    gene_metadata <- base::get0("genbank_table", envir = base::.GlobalEnv, inherits = FALSE)
  }

  slim <- .dnmb_te_build_slim(results, gene_metadata)

  xlsx_path <- base::file.path(output_dir, "translation_efficiency_summary.xlsx")
  tsv_path <- base::file.path(output_dir, "translation_efficiency_summary.tsv")
  readr::write_tsv(slim, tsv_path)
  ok_xlsx <- tryCatch({
    .dnmb_te_write_xlsx(slim, xlsx_path, organism = organism)
    TRUE
  }, error = function(e) {
    base::warning("Excel write failed (", conditionMessage(e), "); TSV is available.", call. = FALSE)
    FALSE
  })

  pdf_path <- base::file.path(output_dir, "translation_efficiency_distribution.pdf")
  ok_pdf <- tryCatch(
    .dnmb_te_distribution_plot(slim, pdf_path, organism = organism, top_n = top_n),
    error = function(e) {
      base::warning("Distribution PDF failed: ", conditionMessage(e), call. = FALSE)
      FALSE
    }
  )

  files <- list(tsv = tsv_path)
  if (base::isTRUE(ok_xlsx)) files$xlsx <- xlsx_path
  if (base::isTRUE(ok_pdf)) files$pdf <- pdf_path
  list(slim_table = slim, files = files, n = base::nrow(slim))
}

.dnmb_te_build_slim <- function(results, gene_metadata) {
  results <- base::as.data.frame(results, stringsAsFactors = FALSE)
  pull_num <- function(col) {
    if (col %in% base::names(results)) suppressWarnings(base::as.numeric(results[[col]]))
    else base::rep(NA_real_, base::nrow(results))
  }
  pull_chr <- function(col) {
    if (col %in% base::names(results)) base::as.character(results[[col]])
    else base::rep(NA_character_, base::nrow(results))
  }

  out <- data.frame(
    locus_tag = pull_chr("locus_tag"),
    tir_score = base::round(pull_num("tir_score"), 2),
    tir_score_percentile = base::round(pull_num("tir_score_percentile"), 1),
    codon_efficiency_score = base::round(pull_num("codon_efficiency_score"), 2),
    cai = base::round(pull_num("cai"), 3),
    tai = base::round(pull_num("tai"), 3),
    expression_band = pull_chr("expression_band"),
    start_codon = pull_chr("start_codon"),
    rbs_motif = pull_chr("rbs_motif"),
    rbs_spacer = base::as.integer(pull_num("rbs_spacer")),
    rbs_score = base::round(pull_num("rbs_score"), 1),
    duplex_score = base::round(pull_num("duplex_score"), 1),
    accessibility_score = base::round(pull_num("accessibility_score"), 1),
    ncs_mfe_dg = base::round(pull_num("ncs_mfe_dg"), 2),
    internal_sd_count = base::as.integer(pull_num("internal_sd_count")),
    stringsAsFactors = FALSE
  )

  # Enrich with gene metadata (gene name, product, protein_id) when available
  if (base::is.data.frame(gene_metadata) && base::nrow(gene_metadata) &&
      "locus_tag" %in% base::names(gene_metadata)) {
    add_cols <- base::intersect(c("gene", "product", "protein_id", "contig"), base::names(gene_metadata))
    if (base::length(add_cols)) {
      meta <- base::as.data.frame(
        gene_metadata[, c("locus_tag", add_cols), drop = FALSE],
        stringsAsFactors = FALSE
      )
      meta$locus_tag <- base::as.character(meta$locus_tag)
      meta <- meta[!base::duplicated(meta$locus_tag), , drop = FALSE]
      idx <- base::match(out$locus_tag, meta$locus_tag)
      for (col in add_cols) {
        out[[col]] <- meta[[col]][idx]
      }
      first_cols <- base::intersect(c("locus_tag", "gene", "product", "protein_id", "contig"), base::names(out))
      out <- out[, c(first_cols, base::setdiff(base::names(out), first_cols)), drop = FALSE]
    }
  }

  out <- out[base::order(-out$tir_score), , drop = FALSE]
  base::row.names(out) <- NULL
  out
}

.dnmb_te_write_xlsx <- function(slim, path, organism = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx package required for Excel output.", call. = FALSE)
  }
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Translation_Efficiency")
  openxlsx::writeData(wb, "Translation_Efficiency", slim, headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#E5E7EB"))

  band_col <- base::which(base::names(slim) == "expression_band")
  if (base::length(band_col)) {
    band_styles <- list(
      very_high = openxlsx::createStyle(fgFill = "#15803D", fontColour = "#FFFFFF", textDecoration = "bold"),
      high = openxlsx::createStyle(fgFill = "#22C55E"),
      moderate = openxlsx::createStyle(fgFill = "#FACC15"),
      low = openxlsx::createStyle(fgFill = "#F97316"),
      very_low = openxlsx::createStyle(fgFill = "#B91C1C", fontColour = "#FFFFFF", textDecoration = "bold")
    )
    for (band in base::names(band_styles)) {
      rows <- base::which(slim$expression_band == band)
      if (base::length(rows)) {
        openxlsx::addStyle(wb, "Translation_Efficiency",
                           style = band_styles[[band]],
                           rows = rows + 1L,
                           cols = band_col,
                           gridExpand = TRUE)
      }
    }
  }
  openxlsx::setColWidths(wb, "Translation_Efficiency",
                         cols = base::seq_len(base::ncol(slim)),
                         widths = "auto")
  openxlsx::freezePane(wb, "Translation_Efficiency", firstRow = TRUE, firstCol = TRUE)

  # Interpretation sheet
  openxlsx::addWorksheet(wb, "Interpretation")
  interp <- data.frame(
    column = c(
      "tir_score",
      "tir_score_percentile",
      "codon_efficiency_score",
      "cai",
      "tai",
      "expression_band",
      "start_codon",
      "rbs_motif",
      "rbs_spacer",
      "rbs_score",
      "duplex_score",
      "accessibility_score",
      "ncs_mfe_dg",
      "internal_sd_count"
    ),
    range = c(
      "0-100", "0-100", "0-100", "0-1", "0-1",
      "very_low / low / moderate / high / very_high",
      "ATG / GTG / TTG / CTG / other",
      "0-9 nt fragment of SD-like motif",
      "3-15 nt", "0-100", "0-100", "0-100",
      "kcal/mol (negative = more folded)",
      "0+"
    ),
    interpretation = c(
      "Translation INITIATION composite. Combines RBS strength, anti-SD duplex, local RNAplfold accessibility, upstream A/U enhancer, start codon, early lysine, tir_core A/G, fold MFE; subtracts internal SD penalty. Higher = ribosome initiates more easily.",
      "Within-genome rank of tir_score. 90+ = top 10% initiation efficiency in this organism.",
      "Translation ELONGATION composite = mean(CAI_score, tAI_score). Cross-validated against PaxDB on 3 species: Spearman 0.34-0.58 with measured protein abundance. SINGLE BEST predictor of steady-state protein level.",
      "Codon Adaptation Index (Sharp & Li 1987) using ribosomal proteins as reference. 0.7+ = strongly optimized codon usage.",
      "tRNA Adaptation Index (dos Reis 2004) using genome tRNA gene copy number. 0.4+ = well-matched to tRNA pool.",
      "Integrated 5-level classifier from quartile sums of tir_score and codon_efficiency_score. very_high = top expression candidate, very_low = likely poorly expressed.",
      "Annotated start codon. ATG strongly preferred in highly expressed genes; GTG/TTG in regulated/low-abundance genes.",
      "Best-matching SD-like fragment found in -35..-3 upstream window.",
      "Distance between SD motif and start codon. Optimal 7-9 nt; Gaussian-weighted.",
      "Score derived from SD motif quality and spacer optimality (0-100).",
      "Anti-SD/SD duplex free energy converted to a 0-100 score (more negative dG = higher score).",
      "RNAplfold local unpaired probability across RBS, start, TIR(-18..+10), standby, and downstream(+1..+30) windows, weighted-meaned then x100.",
      "RNAfold MFE of the +1..+30 window only (Kudla 2009). Closer to 0 = unfolded near AUG = better initiation.",
      "Number of internal Shine-Dalgarno-like motifs found in the first 60 nt of CDS. >0 may slow translation (Li et al. 2012)."
    ),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, "Interpretation", interp, headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#E5E7EB"))
  openxlsx::setColWidths(wb, "Interpretation", cols = 1, widths = 22)
  openxlsx::setColWidths(wb, "Interpretation", cols = 2, widths = 30)
  openxlsx::setColWidths(wb, "Interpretation", cols = 3, widths = 100)
  openxlsx::freezePane(wb, "Interpretation", firstRow = TRUE)

  # README sheet
  openxlsx::addWorksheet(wb, "README")
  readme <- data.frame(
    section = c("Source", "Methodology", "Validation", "How to use", "Caveats"),
    content = c(
      base::sprintf("DNMB mRNAcal module%s.",
                    if (!base::is.null(organism)) base::paste0(" — ", organism) else ""),
      "Each gene's translation efficiency is decomposed into INITIATION (tir_score) and ELONGATION (codon_efficiency_score) signals. They are intentionally orthogonal because they capture different rate-limiting steps.",
      "Cross-validated against PaxDB integrated abundance on E. coli K-12 (n=4078 CDS, codon_efficiency Spearman 0.580), B. subtilis 168 (n=4004, 0.404), and C. jejuni NCTC 11168 (n=768, 0.384). The codon_efficiency_score is the single best abundance predictor.",
      "Sort the Translation_Efficiency sheet by codon_efficiency_score for steady-state abundance ranking, by tir_score for initiation efficiency ranking, or by expression_band for a coarse 5-level overall classification.",
      "tir_score on its own does NOT predict abundance well (Spearman ~0.05 in E. coli/B. subtilis); it captures initiation rate, which is one of many factors. For overall abundance, use codon_efficiency_score or expression_band."
    ),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, "README", readme, headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#E5E7EB"))
  openxlsx::setColWidths(wb, "README", cols = 1, widths = 14)
  openxlsx::setColWidths(wb, "README", cols = 2, widths = 130)

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  TRUE
}

.dnmb_te_distribution_plot <- function(slim, path, organism = NULL, top_n = 20L) {
  if (!base::nrow(slim)) return(FALSE)
  expr_cols <- c(very_high = "#15803D", high = "#22C55E", moderate = "#FACC15",
                 low = "#F97316", very_low = "#B91C1C")
  slim$expr_band <- factor(slim$expression_band,
                           levels = c("very_high", "high", "moderate", "low", "very_low"))
  org_label <- if (!base::is.null(organism)) base::paste0(" — ", organism) else ""

  hist_panel <- function(col, title, xlab) {
    d <- data.frame(value = slim[[col]], band = slim$expr_band)
    d <- d[!base::is.na(d$value), , drop = FALSE]
    if (!base::nrow(d)) return(NULL)
    rng <- base::range(d$value, na.rm = TRUE) + c(-1, 1) * base::diff(base::range(d$value, na.rm = TRUE)) * 0.04
    bw <- base::max(base::diff(base::range(d$value, na.rm = TRUE)) / 40, 0.005)
    ggplot2::ggplot(d, ggplot2::aes(x = .data$value, fill = .data$band)) +
      ggplot2::geom_histogram(binwidth = bw, color = "white", linewidth = 0.15, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = expr_cols, drop = FALSE, na.value = "#9CA3AF") +
      ggplot2::coord_cartesian(xlim = rng) +
      ggplot2::labs(title = title, x = xlab, y = "Genes", fill = "expression_band") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"))
  }

  p1 <- hist_panel("tir_score", "Translation INITIATION (tir_score)", "tir_score (0-100)")
  p2 <- hist_panel("codon_efficiency_score", "Translation ELONGATION (codon_efficiency_score)", "codon_efficiency_score (0-100)")
  p3 <- hist_panel("cai", "CAI (Sharp & Li 1987)", "Codon Adaptation Index (0-1)")
  p4 <- hist_panel("tai", "tAI (dos Reis 2004)", "tRNA Adaptation Index (0-1)")

  scatter_d <- slim[!base::is.na(slim$tir_score) & !base::is.na(slim$codon_efficiency_score), , drop = FALSE]
  p_scatter <- ggplot2::ggplot(scatter_d, ggplot2::aes(x = .data$tir_score, y = .data$codon_efficiency_score, color = .data$expr_band)) +
    ggplot2::geom_point(alpha = 0.55, size = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = expr_cols, drop = FALSE) +
    ggplot2::labs(title = "Initiation vs Elongation", x = "tir_score (initiation)", y = "codon_efficiency_score (elongation)", color = "expression_band") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))

  # Top / bottom 20 tables
  top_idx <- base::order(-slim$codon_efficiency_score, na.last = NA)
  top_tbl <- slim[utils::head(top_idx, top_n), , drop = FALSE]
  bot_tbl <- slim[utils::tail(top_idx, top_n), , drop = FALSE]
  bot_tbl <- bot_tbl[base::order(bot_tbl$codon_efficiency_score), , drop = FALSE]

  table_panel <- function(d, title) {
    cols <- base::intersect(c("locus_tag", "gene", "product",
                              "tir_score", "codon_efficiency_score", "expression_band"),
                            base::names(d))
    if (!base::length(cols)) return(NULL)
    d2 <- d[, cols, drop = FALSE]
    if ("product" %in% base::names(d2)) {
      d2$product <- base::ifelse(
        base::is.na(d2$product), "",
        base::substr(base::as.character(d2$product), 1L, 38L)
      )
    }
    gridExtra::tableGrob(d2, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 7))
  }
  ok_grid <- requireNamespace("gridExtra", quietly = TRUE)
  if (ok_grid) {
    g_top <- table_panel(top_tbl, "Top by codon_efficiency")
    g_bot <- table_panel(bot_tbl, "Bottom by codon_efficiency")
  } else {
    g_top <- ggplot2::ggplot() + ggplot2::labs(title = "(install gridExtra to render top genes)")
    g_bot <- g_top
  }

  hist_grid <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, rel_heights = c(1, 1))
  table_grid <- cowplot::plot_grid(g_top, g_bot, ncol = 2)

  full <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
      base::paste0("Translation Efficiency Summary", org_label,
                   " — n=", base::nrow(slim), " CDS"),
      fontface = "bold", size = 12, hjust = 0.5
    ),
    hist_grid,
    p_scatter,
    table_grid,
    ncol = 1,
    rel_heights = c(0.07, 1.4, 1.0, 1.6)
  )
  ggplot2::ggsave(path, full, width = 12, height = 16, bg = "white", limitsize = FALSE)
  TRUE
}
