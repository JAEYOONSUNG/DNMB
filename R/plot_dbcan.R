.dnmb_dbcan_plot_gene_track <- function(tbl) {
  tbl <- tbl[order(tbl$dbcan_cgc_id, tbl$start, tbl$end), , drop = FALSE]
  tbl$gene_type <- ifelse(is.na(tbl$dbcan_cgc_gene_type) | !nzchar(tbl$dbcan_cgc_gene_type), "other", tbl$dbcan_cgc_gene_type)
  tbl$locus_label <- as.character(tbl$locus_tag)
  tbl$family_label <- ifelse(is.na(tbl$dbCAN_family_id) | !nzchar(tbl$dbCAN_family_id), "", as.character(tbl$dbCAN_family_id))
  tbl$label <- ifelse(
    tbl$gene_type == "CAZyme" & nzchar(tbl$family_label),
    paste0(tbl$locus_label, "\n", tbl$family_label),
    tbl$locus_label
  )
  cgc_meta <- tbl[!duplicated(tbl$dbcan_cgc_id), c("dbcan_cgc_id", "dbcan_pul_substrate"), drop = FALSE]
  cgc_meta$facet_label <- ifelse(
    is.na(cgc_meta$dbcan_pul_substrate) | !nzchar(cgc_meta$dbcan_pul_substrate),
    cgc_meta$dbcan_cgc_id,
    paste0(cgc_meta$dbcan_cgc_id, " | substrate: ", cgc_meta$dbcan_pul_substrate)
  )
  tbl$facet_label <- cgc_meta$facet_label[match(tbl$dbcan_cgc_id, cgc_meta$dbcan_cgc_id)]
  type_palette <- c(
    CAZyme = "#0F766E",
    TC = "#D97706",
    TF = "#7C3AED",
    STP = "#DC2626",
    other = "#9CA3AF"
  )
  tbl$y <- 1
  backbone <- tbl |>
    dplyr::group_by(.data$facet_label) |>
    dplyr::summarise(start = min(.data$start), end = max(.data$end), .groups = "drop")
  ggplot2::ggplot(tbl) +
    ggplot2::geom_segment(
      data = backbone,
      ggplot2::aes(x = .data$start, xend = .data$end, y = 1, yend = 1),
      linewidth = 1.1,
      color = "grey80",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = .data$y - 0.18,
        ymax = .data$y + 0.18,
        fill = .data$gene_type
      ),
      color = "grey20",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = (.data$start + .data$end) / 2,
        y = .data$y + 0.42,
        label = .data$locus_label
      ),
      size = 2.8,
      angle = 90,
      vjust = 0,
      hjust = 0.5,
      color = "black",
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = tbl[tbl$gene_type == "CAZyme" & nzchar(tbl$family_label), , drop = FALSE],
      ggplot2::aes(
        x = (.data$start + .data$end) / 2,
        y = .data$y - 0.34,
        label = .data$family_label
      ),
      size = 2.5,
      angle = 90,
      vjust = 1,
      hjust = 0.5,
      color = "#0F766E",
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~facet_label, scales = "free_x", ncol = 1) +
    ggplot2::scale_fill_manual(values = type_palette, drop = FALSE) +
    ggplot2::labs(
      title = "dbCAN CGC/PUL Summary",
      x = "Genome coordinate (bp)",
      y = NULL,
      fill = "Gene type"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(1.1, "cm"),
      plot.margin = ggplot2::margin(12, 12, 18, 12)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

.dnmb_plot_dbcan_module <- function(genbank_table, output_dir) {
  tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
  cgc_id_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_id", "dbCAN_cgc_id", "dbcan_cgc_id"))
  cgc_type_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_gene_type", "dbCAN_cgc_gene_type", "dbcan_cgc_gene_type"))
  cgc_family_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_protein_family", "dbCAN_cgc_protein_family", "dbcan_cgc_protein_family"))
  substrate_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_pul_substrate", "dbCAN_pul_substrate", "dbcan_pul_substrate"))
  family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))

  required <- c("locus_tag", "contig", "start", "end", "direction")
  if (!all(required %in% names(tbl)) || is.null(cgc_id_col) || is.null(cgc_type_col)) {
    return(NULL)
  }

  tbl$dbcan_cgc_id <- as.character(tbl[[cgc_id_col]])
  tbl$dbcan_cgc_gene_type <- as.character(tbl[[cgc_type_col]])
  tbl$dbcan_cgc_protein_family <- if (!is.null(cgc_family_col)) as.character(tbl[[cgc_family_col]]) else NA_character_
  tbl$dbcan_pul_substrate <- if (!is.null(substrate_col)) as.character(tbl[[substrate_col]]) else NA_character_
  tbl$dbCAN_family_id <- if (!is.null(family_col)) as.character(tbl[[family_col]]) else NA_character_
  tbl <- tbl[!is.na(tbl$dbcan_cgc_id) & nzchar(tbl$dbcan_cgc_id), , drop = FALSE]
  if (!nrow(tbl)) {
    return(NULL)
  }

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  p <- .dnmb_dbcan_plot_gene_track(tbl)
  pdf_path <- file.path(plot_dir, "dbcan_cgc_overview.pdf")
  .dnmb_module_plot_save(p, pdf_path = pdf_path, width = 14, height = max(5, 3.2 * length(unique(tbl$dbcan_cgc_id))))
  list(pdf = pdf_path)
}
#' Internal dbCAN plotting helpers
#'
#' Plot-construction routines for the dbCAN overview figures rendered by DNMB.
#'
#' @name dnmb_internal_plot_dbcan
#' @keywords internal
#' @noRd
NULL

