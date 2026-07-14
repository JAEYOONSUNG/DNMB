.dnmb_padloc_plot_data <- function(genbank_table, output_dir) {
  replicons <- .dnmb_module_replicon_plot_data(genbank_table, output_dir = output_dir)
  if (is.null(replicons)) {
    return(NULL)
  }
  genes <- replicons$genes
  module_dir <- file.path(output_dir, "dnmb_module_padloc")
  csv_paths <- if (dir.exists(module_dir)) {
    list.files(module_dir, pattern = "_padloc\\.csv$", full.names = TRUE)
  } else {
    character()
  }
  map_paths <- if (dir.exists(module_dir)) {
    list.files(module_dir, pattern = "_id_map\\.tsv$", full.names = TRUE)
  } else {
    character()
  }

  hits <- NULL
  if (length(csv_paths) && length(map_paths)) {
    raw <- tryCatch(
      utils::read.csv(csv_paths[[1]], stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    id_map <- tryCatch(
      utils::read.delim(map_paths[[1]], stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    required_raw <- c(
      "system.number", "seqid", "system", "target.name", "protein.name",
      "domain.iE.value", "target.coverage", "hmm.coverage", "start", "end"
    )
    if (!is.null(raw) && nrow(raw) && !is.null(id_map) && nrow(id_map) &&
        all(required_raw %in% names(raw)) &&
        all(c("query_id", "locus_tag") %in% names(id_map))) {
      map_idx <- match(as.character(raw[["target.name"]]), as.character(id_map$query_id))
      locus_tag <- as.character(id_map$locus_tag[map_idx])
      gene_idx <- match(locus_tag, as.character(genes$locus_tag))
      raw_replicon <- as.character(raw$seqid)
      if ("contig_id" %in% names(id_map)) {
        missing_replicon <- is.na(raw_replicon) | !nzchar(raw_replicon)
        raw_replicon[missing_replicon] <- as.character(id_map$contig_id[map_idx[missing_replicon]])
      }
      hits <- data.frame(
        replicon_id = raw_replicon,
        start = suppressWarnings(as.numeric(raw$start)),
        end = suppressWarnings(as.numeric(raw$end)),
        locus_tag = locus_tag,
        gene = if ("gene" %in% names(genes)) as.character(genes$gene[gene_idx]) else NA_character_,
        product = if ("product" %in% names(genes)) as.character(genes$product[gene_idx]) else NA_character_,
        PADLOC_system_number = suppressWarnings(as.integer(raw[["system.number"]])),
        PADLOC_system = as.character(raw$system),
        PADLOC_protein_name = as.character(raw[["protein.name"]]),
        PADLOC_hmm_name = if ("hmm.name" %in% names(raw)) as.character(raw[["hmm.name"]]) else NA_character_,
        PADLOC_domain_ievalue = suppressWarnings(as.numeric(raw[["domain.iE.value"]])),
        PADLOC_target_coverage = suppressWarnings(as.numeric(raw[["target.coverage"]])),
        PADLOC_hmm_coverage = suppressWarnings(as.numeric(raw[["hmm.coverage"]])),
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.null(hits) || !nrow(hits)) {
    required_combined <- c("PADLOC_system", "PADLOC_system_number")
    if (!all(required_combined %in% names(genes))) {
      return(NULL)
    }
    keep <- !is.na(genes$PADLOC_system) & nzchar(as.character(genes$PADLOC_system))
    hits <- genes[keep, , drop = FALSE]
    if (!nrow(hits)) {
      return(NULL)
    }
    hits$PADLOC_system_number <- suppressWarnings(as.integer(hits$PADLOC_system_number))
    hits$PADLOC_target_coverage <- suppressWarnings(as.numeric(hits$PADLOC_target_coverage))
    hits$PADLOC_hmm_coverage <- suppressWarnings(as.numeric(hits$PADLOC_hmm_coverage))
    hits$PADLOC_domain_ievalue <- if ("PADLOC_domain_ievalue" %in% names(hits)) {
      suppressWarnings(as.numeric(hits$PADLOC_domain_ievalue))
    } else {
      NA_real_
    }
    if (!"PADLOC_hmm_name" %in% names(hits)) hits$PADLOC_hmm_name <- NA_character_
    if (!"PADLOC_protein_name" %in% names(hits)) hits$PADLOC_protein_name <- NA_character_
  }

  hits <- hits[
    !is.na(hits$PADLOC_system) & nzchar(hits$PADLOC_system) &
      !is.na(hits$locus_tag) & nzchar(hits$locus_tag) &
      is.finite(hits$start) & is.finite(hits$end),
    , drop = FALSE
  ]
  if (!nrow(hits)) {
    return(NULL)
  }
  hits$midpoint <- (hits$start + hits$end) / 2
  hits$replicon_order <- replicons$contigs$replicon_order[
    match(hits$replicon_id, replicons$contigs$replicon_id)
  ]
  hits$gene[is.na(hits$gene)] <- ""
  hits$product[is.na(hits$product) | !nzchar(hits$product)] <- "hypothetical protein"
  list(hits = hits, contigs = replicons$contigs)
}

.dnmb_plot_padloc_module <- function(genbank_table, output_dir) {
  add_panel_header <- function(plot_obj, label, title, x = 0.0005, y = 0.975,
                               size = 12, plot_y = 0.035, plot_h = 0.90) {
    cowplot::ggdraw() +
      cowplot::draw_plot(plot_obj, x = 0, y = plot_y, width = 1, height = plot_h) +
      cowplot::draw_label(
        paste0(label, "  ", title), x = x, y = y,
        hjust = 0, vjust = 1, fontface = "bold", size = size,
        fontfamily = .dnmb_plot_font_family()
      )
  }

  plot_data <- .dnmb_padloc_plot_data(genbank_table, output_dir)
  if (is.null(plot_data)) {
    return(NULL)
  }
  hits <- plot_data$hits
  contigs <- plot_data$contigs

  system_summary <- hits |>
    dplyr::group_by(.data$PADLOC_system) |>
    dplyr::summarise(
      n_memberships = dplyr::n(),
      n_genes = dplyr::n_distinct(.data$locus_tag),
      n_loci = dplyr::n_distinct(.data$PADLOC_system_number),
      best_target_cov = max(.data$PADLOC_target_coverage, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_loci), dplyr::desc(.data$n_memberships), .data$PADLOC_system)
  system_summary$best_target_cov[!is.finite(system_summary$best_target_cov)] <- NA_real_
  system_levels <- as.character(system_summary$PADLOC_system)
  palette <- stats::setNames(
    grDevices::colorRampPalette(c(
      "#F8C35A", "#F2A93B", "#F08F20", "#EB7C12",
      "#E16A0D", "#D25709", "#BE4306"
    ))(length(system_levels)),
    system_levels
  )
  system_summary$PADLOC_system <- factor(
    system_summary$PADLOC_system,
    levels = rev(system_levels)
  )
  inventory_detail <- paste0(
    system_summary$n_memberships, " memberships / ",
    system_summary$n_genes, " genes | best target cov ",
    ifelse(is.na(system_summary$best_target_cov), "NA", sprintf("%.2f", system_summary$best_target_cov))
  )

  base_theme <- ggplot2::theme_bw(
    base_size = 10.5, base_family = .dnmb_plot_font_family()
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(family = .dnmb_plot_font_family())
    )

  p_inventory <- ggplot2::ggplot(
    system_summary,
    ggplot2::aes(x = .data$n_loci, y = .data$PADLOC_system)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$PADLOC_system),
      width = 0.62, color = "grey35", linewidth = 0.24,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$n_loci + 0.08, label = inventory_detail),
      hjust = 0, size = 2.45, color = "grey25",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(max(system_summary$n_loci)),
      expand = ggplot2::expansion(mult = c(0.01, 0.62))
    ) +
    ggplot2::labs(x = "Confirmed system loci", y = NULL) +
    base_theme +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7.4),
      legend.position = "none"
    )

  layout_genes <- hits |>
    dplyr::arrange(.data$replicon_id, .data$start, .data$PADLOC_domain_ievalue) |>
    dplyr::group_by(.data$replicon_id, .data$locus_tag) |>
    dplyr::summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      midpoint = (start + end) / 2,
      gene = dplyr::first(.data$gene),
      product = dplyr::first(.data$product),
      PADLOC_system = dplyr::first(.data$PADLOC_system),
      system_memberships = paste(unique(.data$PADLOC_system), collapse = " + "),
      priority = max(.data$PADLOC_target_coverage, na.rm = TRUE),
      .groups = "drop"
    )
  layout_genes$priority[!is.finite(layout_genes$priority)] <- 0
  layout_genes$feature_label <- .dnmb_module_feature_label(
    layout_genes$gene, layout_genes$locus_tag
  )
  layout_genes <- .dnmb_module_pack_replicon_labels(
    layout_genes, contigs,
    label_col = "feature_label", priority_col = "priority"
  )
  layout_genes$label_y <- 0.50 + 0.12 * layout_genes$label_tier
  max_label_tier <- max(layout_genes$label_tier, na.rm = TRUE)
  layout_ymax <- 0.50 + 0.12 * max_label_tier + 0.12
  contigs$contig_facet <- factor(contigs$facet_label, levels = contigs$facet_label)
  layout_genes$contig_facet <- factor(
    contigs$facet_label[match(layout_genes$replicon_id, contigs$replicon_id)],
    levels = contigs$facet_label
  )
  empty_contigs <- contigs[!contigs$replicon_id %in% layout_genes$replicon_id, , drop = FALSE]
  empty_contigs$no_call_label <- rep("0 confirmed PADLOC systems", nrow(empty_contigs))

  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = 0.18, yend = 0.18),
      linewidth = 1.05, color = "grey80"
    ) +
    ggplot2::geom_rect(
      data = layout_genes,
      ggplot2::aes(
        xmin = .data$start, xmax = .data$end,
        ymin = 0.145, ymax = 0.215, fill = .data$PADLOC_system
      ),
      color = "grey25", linewidth = 0.20, alpha = 0.95
    ) +
    ggplot2::geom_segment(
      data = layout_genes,
      ggplot2::aes(
        x = .data$midpoint, xend = .data$label_x,
        y = 0.23, yend = .data$label_y - 0.04
      ),
      linewidth = 0.23, color = "grey58", alpha = 0.82
    ) +
    ggplot2::geom_text(
      data = layout_genes,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$feature_label),
      size = 2.12, lineheight = 0.88, color = "grey20",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::geom_text(
      data = empty_contigs,
      ggplot2::aes(x = .data$length_bp / 2, y = 0.50, label = .data$no_call_label),
      size = 2.55, color = "grey48", fontface = "italic",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::facet_wrap(
      ~contig_facet, ncol = 1, scales = "free_x", strip.position = "top"
    ) +
    ggplot2::scale_fill_manual(
      values = palette, breaks = system_levels, name = "PADLOC system"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0.012, 0.032))
    ) +
    ggplot2::coord_cartesian(ylim = c(0.08, layout_ymax), clip = "on") +
    ggplot2::labs(x = "Replicon coordinate (bp)", y = NULL) +
    base_theme +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold", size = 8.5, margin = ggplot2::margin(t = 1, b = 1)
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey97", colour = "grey85", linewidth = 0.3
      ),
      panel.spacing.y = grid::unit(0.10, "lines"),
      legend.position = "none"
    )

  hits$replicon_short <- contigs$replicon_short[
    match(hits$replicon_id, contigs$replicon_id)
  ]
  hits$feature_id <- ifelse(
    !is.na(hits$gene) & nzchar(trimws(hits$gene)),
    paste0(trimws(hits$gene), " | ", hits$locus_tag),
    hits$locus_tag
  )
  hits$call_label <- paste0(
    hits$feature_id, " · ", hits$PADLOC_protein_name,
    " [", hits$PADLOC_system, " #", hits$PADLOC_system_number, "]"
  )
  hits <- hits[
    order(
      -hits$PADLOC_target_coverage, -hits$PADLOC_hmm_coverage,
      hits$PADLOC_domain_ievalue, hits$locus_tag, hits$PADLOC_system
    ),
    , drop = FALSE
  ]
  hits$call_label <- factor(hits$call_label, levels = rev(hits$call_label))
  product_short <- ifelse(
    nchar(hits$product) > 43,
    paste0(substr(hits$product, 1, 40), "..."),
    hits$product
  )
  ie_label <- ifelse(
    is.finite(hits$PADLOC_domain_ievalue),
    format(signif(hits$PADLOC_domain_ievalue, 2), scientific = TRUE),
    "NA"
  )
  hits$evidence_label <- paste0(
    product_short,
    " | HMM cov ", sprintf("%.2f", hits$PADLOC_hmm_coverage),
    " | i-E ", ie_label
  )
  annotation_x <- 1.05
  table_xmax <- 2.52

  p_hits <- ggplot2::ggplot(
    hits,
    ggplot2::aes(x = .data$PADLOC_target_coverage, y = .data$call_label)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$PADLOC_system),
      width = 0.64, color = "grey35", linewidth = 0.22
    ) +
    ggplot2::geom_vline(
      xintercept = annotation_x - 0.025, linewidth = 0.35, color = "grey82"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = annotation_x, label = .data$evidence_label),
      hjust = 0, size = 1.76, color = "grey27",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      limits = c(0, table_xmax),
      expand = ggplot2::expansion(mult = c(0, 0.006))
    ) +
    ggplot2::labs(
      x = "PADLOC target coverage | all final system-protein memberships",
      y = NULL
    ) +
    base_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 5.7, color = "grey20"),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.title.x = ggplot2::element_text(size = 8.2),
      plot.margin = ggplot2::margin(5, 8, 4, 4)
    )

  legend_plot <- p_layout +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = 8.6),
      legend.text = ggplot2::element_text(size = 7.5),
      legend.key.width = grid::unit(0.82, "lines"),
      legend.key.height = grid::unit(0.72, "lines"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  common_legend <- .dnmb_with_plot_pdf_device(
    cowplot::get_legend(legend_plot), width = 10, height = 1
  )
  legend_panel <- cowplot::ggdraw() + cowplot::draw_grob(
    common_legend, x = 0.5, y = 0.5, width = 0.99, height = 0.96,
    hjust = 0.5, vjust = 0.5
  )

  inventory_h <- max(1.65, 0.18 * nrow(system_summary) + 0.55)
  layout_h <- max(
    2.05,
    0.58 * nrow(contigs) + 0.17 * (max_label_tier + 1) + 0.65
  )
  calls_h <- max(3.2, 0.115 * nrow(hits) + 0.95)
  legend_h <- 0.68
  total_h <- inventory_h + layout_h + calls_h + legend_h
  composite <- .dnmb_with_plot_pdf_device(
    cowplot::plot_grid(
      add_panel_header(
        p_inventory, "A", "PADLOC confirmed system inventory",
        plot_y = 0.04, plot_h = 0.88
      ),
      add_panel_header(
        p_layout, "B", "PADLOC genes across chromosome and plasmids",
        plot_y = 0.03, plot_h = 0.90
      ),
      add_panel_header(
        p_hits, "C", "All confirmed PADLOC protein memberships",
        plot_y = 0.025, plot_h = 0.925
      ),
      legend_panel,
      ncol = 1,
      rel_heights = c(inventory_h, layout_h, calls_h, legend_h)
    ),
    width = 10, height = total_h
  )
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "PADLOC_overview.pdf")
  .dnmb_module_plot_save(
    composite, pdf_path, width = 10,
    height = total_h
  )
  list(
    pdf = pdf_path,
    confirmed_memberships = nrow(hits),
    confirmed_genes = dplyr::n_distinct(hits$locus_tag),
    confirmed_loci = dplyr::n_distinct(paste(hits$replicon_id, hits$PADLOC_system_number)),
    replicons_shown = nrow(contigs)
  )
}
