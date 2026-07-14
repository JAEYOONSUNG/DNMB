.dnmb_antidefense_module_prefix <- function(module_name) {
  module_name <- trimws(as.character(module_name)[1])
  switch(
    module_name,
    dbAPIS = "dbAPIS",
    AcrFinder = "AcrFinder",
    stop("Unsupported anti-defense module: ", module_name, call. = FALSE)
  )
}

.dnmb_antidefense_module_file_stub <- function(module_name) {
  module_name <- trimws(as.character(module_name)[1])
  switch(
    module_name,
    dbAPIS = "dbAPIS",
    AcrFinder = "AcrFinder",
    stop("Unsupported anti-defense module: ", module_name, call. = FALSE)
  )
}

.dnmb_antidefense_module_dir <- function(output_dir, module_name) {
  stub <- switch(
    trimws(as.character(module_name)[1]),
    dbAPIS = "dbapis",
    AcrFinder = "acrfinder",
    stop("Unsupported anti-defense module: ", module_name, call. = FALSE)
  )
  file.path(output_dir, paste0("dnmb_module_", stub))
}

.dnmb_antidefense_panel_header <- function(plot_obj,
                                           label,
                                           title,
                                           x = 0.0005,
                                           y = 0.968,
                                           size = 12,
                                           plot_y = 0.045,
                                           plot_h = 0.88) {
  cowplot::ggdraw() +
    cowplot::draw_plot(plot_obj, x = 0, y = plot_y, width = 1, height = plot_h) +
    cowplot::draw_label(
      paste0(label, "  ", title),
      x = x,
      y = y,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = size,
      fontfamily = .dnmb_plot_font_family()
    )
}

.dnmb_antidefense_display_palette <- function(values,
                                              colors = c(
                                                "#F2A93B",
                                                "#F08F20",
                                                "#EB7C12",
                                                "#E16A0D",
                                                "#D25709",
                                                "#BE4306"
                                              )) {
  values <- unique(stats::na.omit(as.character(values)))
  if (!length(values)) {
    return(stats::setNames(character(), character()))
  }
  stats::setNames(grDevices::colorRampPalette(colors)(length(values)), values)
}

.dnmb_antidefense_target_label <- function(x, max_len = 24L) {
  x <- .dnmb_integrated_defense_clean_text(x)
  x <- gsub("NAD\\+ reconstitution pathway \\(NARP\\)", "NARP", x, ignore.case = TRUE)
  x <- gsub("^restriction[- ]modification \\(RM\\)$", "Restriction-modification", x, ignore.case = TRUE)
  x <- gsub("pyrimidine cyclase system for antiphage resistance \\(Pycsar\\)", "Pycsar", x, ignore.case = TRUE)
  x <- gsub("CRISPR[-–]Cas evasion by DNA repair", "CRISPR repair evasion", x, ignore.case = TRUE)
  .dnmb_integrated_defense_short_label(x, max_len = max_len)
}

.dnmb_antidefense_module_has_outputs <- function(output_dir, module_name) {
  module_dir <- .dnmb_antidefense_module_dir(output_dir, module_name)
  if (!dir.exists(module_dir)) {
    return(FALSE)
  }
  any(file.exists(list.files(module_dir, recursive = TRUE, full.names = TRUE)))
}

.dnmb_antidefense_module_hit_table <- function(genbank_table, module_name) {
  prefix <- .dnmb_antidefense_module_prefix(module_name)
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  family_col <- paste0(prefix, "_family_id")
  if (!family_col %in% names(tbl)) {
    return(data.frame())
  }

  hit_label_col <- paste0(prefix, "_hit_label")
  support_col <- paste0(prefix, "_support")
  evidence_col <- paste0(prefix, "_evidence_mode")
  role_col <- paste0(prefix, "_enzyme_role")
  subtype_col <- paste0(prefix, "_defense_type")
  clan_col <- paste0(prefix, "_clan_defense_type")
  pident_col <- paste0(prefix, "_pident")
  evalue_col <- paste0(prefix, "_i_evalue")
  score_col <- paste0(prefix, "_hit_score")
  hmm_cov_col <- paste0(prefix, "_hmm_coverage")
  query_cov_col <- paste0(prefix, "_query_coverage")

  tbl$module_family_id <- as.character(tbl[[family_col]])
  tbl$module_hit_label <- if (hit_label_col %in% names(tbl)) as.character(tbl[[hit_label_col]]) else NA_character_
  tbl$module_support <- if (support_col %in% names(tbl)) as.character(tbl[[support_col]]) else NA_character_
  tbl$module_evidence <- if (evidence_col %in% names(tbl)) as.character(tbl[[evidence_col]]) else NA_character_
  tbl$module_role <- if (role_col %in% names(tbl)) as.character(tbl[[role_col]]) else NA_character_
  tbl$module_subtype <- if (subtype_col %in% names(tbl)) as.character(tbl[[subtype_col]]) else NA_character_
  tbl$module_clan <- if (clan_col %in% names(tbl)) as.character(tbl[[clan_col]]) else NA_character_
  tbl$module_pident <- if (pident_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[pident_col]])) else NA_real_
  tbl$module_evalue <- if (evalue_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[evalue_col]])) else NA_real_
  tbl$module_hit_score <- if (score_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[score_col]])) else NA_real_
  tbl$module_hmm_coverage <- if (hmm_cov_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[hmm_cov_col]])) else NA_real_
  tbl$module_query_coverage <- if (query_cov_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[query_cov_col]])) else NA_real_

  keep <- !is.na(tbl$module_family_id) & nzchar(tbl$module_family_id)
  tbl <- tbl[keep, , drop = FALSE]
  if (!nrow(tbl)) {
    return(tbl)
  }

  tbl$module_display <- dplyr::coalesce(
    ifelse(!is.na(tbl$module_hit_label) & nzchar(tbl$module_hit_label), tbl$module_hit_label, NA_character_),
    tbl$module_family_id
  )
  tbl$module_group <- dplyr::coalesce(
    ifelse(!is.na(tbl$module_subtype) & nzchar(tbl$module_subtype), tbl$module_subtype, NA_character_),
    ifelse(!is.na(tbl$module_clan) & nzchar(tbl$module_clan), tbl$module_clan, NA_character_),
    ifelse(!is.na(tbl$module_role) & nzchar(tbl$module_role), tbl$module_role, NA_character_),
    tbl$module_display
  )
  tbl
}

.dnmb_antidefense_status_plot <- function(module_name, output_dir, status_text) {
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(.dnmb_antidefense_module_file_stub(module_name), "_overview.pdf"))
  p <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.68, label = paste0(module_name, " anti-defense overview"), fontface = "bold", size = 6) +
    ggplot2::annotate("text", x = 0.5, y = 0.42, label = status_text, size = 4.8, color = "grey30") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")
  .dnmb_module_plot_save(p, pdf_path, width = 10, height = 4.8)
  list(pdf = pdf_path)
}

.dnmb_plot_antidefense_module <- function(genbank_table, output_dir, module_name) {
  tbl <- .dnmb_antidefense_module_hit_table(genbank_table, module_name)
  has_outputs <- .dnmb_antidefense_module_has_outputs(output_dir, module_name)
  if (!nrow(tbl)) {
    if (!has_outputs) {
      return(NULL)
    }
    return(.dnmb_antidefense_status_plot(
      module_name = module_name,
      output_dir = output_dir,
      status_text = "Module completed, but no anti-defense hits were detected for this genome."
    ))
  }

  palette <- .dnmb_defensefinder_palette(tbl$module_group)
  summary_tbl <- tbl |>
    dplyr::group_by(.data$module_group, .data$module_display) |>
    dplyr::summarise(
      n_hits = dplyr::n(),
      contigs = dplyr::n_distinct(.data$contig),
      best_pident = suppressWarnings(max(.data$module_pident, na.rm = TRUE)),
      best_evalue = suppressWarnings(min(.data$module_evalue, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_hits), .data$module_group, .data$module_display)

  summary_tbl$best_pident[!is.finite(summary_tbl$best_pident)] <- NA_real_
  summary_tbl$best_evalue[!is.finite(summary_tbl$best_evalue)] <- NA_real_
  summary_tbl$module_display <- factor(summary_tbl$module_display, levels = rev(unique(summary_tbl$module_display)))

  p_inventory <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = .data$n_hits, y = .data$module_display)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$module_group),
      width = 0.62, color = "grey35", linewidth = 0.2, alpha = 0.92
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$n_hits + 0.05,
        label = ifelse(
          is.na(.data$best_pident),
          paste0(.data$contigs, " contigs"),
          paste0(.data$contigs, " contigs | best pident ", format(round(.data$best_pident, 1), nsmall = 1))
        )
      ),
      hjust = 0, size = 2.7, color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.35))) +
    ggplot2::labs(
      title = paste0("A   ", module_name, " inventory"),
      x = "Hit count",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title.position = "plot"
    )

  layout_tbl <- tbl
  contig_levels <- unique(as.character(layout_tbl$contig))
  layout_tbl$contig <- factor(layout_tbl$contig, levels = rev(contig_levels))
  layout_tbl$start <- suppressWarnings(as.numeric(layout_tbl$start))
  layout_tbl$end <- suppressWarnings(as.numeric(layout_tbl$end))
  layout_tbl$mid <- (layout_tbl$start + layout_tbl$end) / 2
  layout_tbl$label_short <- .dnmb_integrated_defense_short_label(layout_tbl$module_display, max_len = 16L)

  p_layout <- ggplot2::ggplot(layout_tbl, ggplot2::aes(y = .data$contig)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$start, xend = .data$end, yend = .data$contig, color = .data$module_group),
      linewidth = 2.4, lineend = "round", alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$mid, color = .data$module_group),
      size = 1.8, alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(x = .data$mid, label = .data$label_short, color = .data$module_group),
      size = 2.6,
      min.segment.length = 0,
      box.padding = 0.18,
      point.padding = 0.12,
      max.overlaps = 100
    ) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(
      title = paste0("B   ", module_name, " genome layout"),
      x = "Genomic position (bp)",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title.position = "plot",
      legend.position = "bottom"
    )

  top_tbl <- summary_tbl |>
    dplyr::mutate(
      best_evalue_label = ifelse(is.na(.data$best_evalue), "NA", format(signif(.data$best_evalue, 3), scientific = TRUE)),
      label = paste0(
        .data$module_display, " [", .data$module_group, "]",
        "\nHits: ", .data$n_hits,
        " | Contigs: ", .data$contigs,
        ifelse(is.na(.data$best_pident), "", paste0(" | Best pident: ", format(round(.data$best_pident, 1), nsmall = 1))),
        "\nBest e-value: ", .data$best_evalue_label
      )
    ) |>
    dplyr::slice_head(n = 12L)

  top_tbl$rank <- factor(seq_len(nrow(top_tbl)), levels = rev(seq_len(nrow(top_tbl))))
  p_hits <- ggplot2::ggplot(top_tbl, ggplot2::aes(y = .data$rank)) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, label = .data$label),
      hjust = 0, vjust = 0.5, size = 3.0, lineheight = 1.08
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(title = paste0("C   Top ", module_name, " calls")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(.dnmb_antidefense_module_file_stub(module_name), "_overview.pdf"))
  composite <- cowplot::plot_grid(
    p_inventory,
    p_layout,
    p_hits,
    ncol = 1,
    rel_heights = c(0.95, 1.45, 1.15),
    align = "v"
  )
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = 13.2)
  list(pdf = pdf_path)
}

.dnmb_dbapis_plot_data <- function(genbank_table, output_dir) {
  replicons <- .dnmb_module_replicon_plot_data(genbank_table, output_dir = output_dir)
  if (is.null(replicons)) {
    return(NULL)
  }
  genes <- replicons$genes
  raw_path <- file.path(output_dir, "dnmb_module_dbapis", "dbapis_hits.tsv")
  raw <- if (file.exists(raw_path)) {
    tryCatch(
      utils::read.delim(raw_path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  required_raw <- c(
    "query", "family_id", "i_evalue", "score", "hmm_coverage",
    "query_coverage"
  )
  hits <- NULL
  if (!is.null(raw) && nrow(raw) && all(required_raw %in% names(raw))) {
    gene_idx <- match(as.character(raw$query), as.character(genes$locus_tag))
    hits <- data.frame(
      replicon_id = as.character(genes$replicon_id[gene_idx]),
      start = suppressWarnings(as.numeric(genes$start[gene_idx])),
      end = suppressWarnings(as.numeric(genes$end[gene_idx])),
      locus_tag = as.character(genes$locus_tag[gene_idx]),
      gene = if ("gene" %in% names(genes)) as.character(genes$gene[gene_idx]) else NA_character_,
      product = if ("product" %in% names(genes)) as.character(genes$product[gene_idx]) else NA_character_,
      module_family_id = as.character(raw$family_id),
      module_hit_label = if ("hit_label" %in% names(raw)) as.character(raw$hit_label) else NA_character_,
      module_subtype = if ("defense_type" %in% names(raw)) as.character(raw$defense_type) else NA_character_,
      module_clan_id = if ("clan_id" %in% names(raw)) as.character(raw$clan_id) else NA_character_,
      module_clan = if ("clan_defense_type" %in% names(raw)) as.character(raw$clan_defense_type) else NA_character_,
      module_description = if ("description" %in% names(raw)) as.character(raw$description) else NA_character_,
      module_evalue = suppressWarnings(as.numeric(raw$i_evalue)),
      module_hit_score = suppressWarnings(as.numeric(raw$score)),
      module_hmm_coverage = suppressWarnings(as.numeric(raw$hmm_coverage)),
      module_query_coverage = suppressWarnings(as.numeric(raw$query_coverage)),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(hits) || !nrow(hits)) {
    hits <- .dnmb_antidefense_module_hit_table(genbank_table, "dbAPIS")
    if (!nrow(hits)) {
      return(NULL)
    }
    contig_number <- if ("contig_number" %in% names(hits)) {
      suppressWarnings(as.integer(hits$contig_number))
    } else {
      as.integer(factor(as.character(hits$contig)))
    }
    hits$replicon_id <- sprintf("DNMB_CONTIG_%03d", contig_number)
    hits$module_clan_id <- if ("dbAPIS_clan_id" %in% names(hits)) as.character(hits$dbAPIS_clan_id) else NA_character_
    hits$module_description <- if ("dbAPIS_description" %in% names(hits)) as.character(hits$dbAPIS_description) else NA_character_
  }

  hits$start <- suppressWarnings(as.numeric(hits$start))
  hits$end <- suppressWarnings(as.numeric(hits$end))
  hits$module_evalue <- suppressWarnings(as.numeric(hits$module_evalue))
  hits$module_hit_score <- suppressWarnings(as.numeric(hits$module_hit_score))
  hits$module_hmm_coverage <- suppressWarnings(as.numeric(hits$module_hmm_coverage))
  hits$module_query_coverage <- suppressWarnings(as.numeric(hits$module_query_coverage))
  hits <- hits[
    !is.na(hits$locus_tag) & nzchar(as.character(hits$locus_tag)) &
      !is.na(hits$module_family_id) & nzchar(as.character(hits$module_family_id)) &
      is.finite(hits$start) & is.finite(hits$end),
    , drop = FALSE
  ]
  if (!nrow(hits)) {
    return(NULL)
  }
  hits$midpoint <- (hits$start + hits$end) / 2
  hits$gene[is.na(hits$gene)] <- ""
  hits$product[is.na(hits$product) | !nzchar(hits$product)] <- "hypothetical protein"
  clean_nonempty <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | !nzchar(x)] <- NA_character_
    x
  }
  hit_label <- clean_nonempty(hits$module_hit_label)
  clan <- clean_nonempty(hits$module_clan)
  subtype <- clean_nonempty(hits$module_subtype)
  hits$module_display <- ifelse(
    !is.na(hit_label),
    paste0(hits$module_family_id, " (", hit_label, ")"),
    hits$module_family_id
  )
  hits$module_group <- dplyr::coalesce(
    clan, subtype, rep("Unclassified anti-defense", nrow(hits))
  )
  list(hits = hits, contigs = replicons$contigs)
}

.dnmb_plot_dbapis_module <- function(genbank_table, output_dir) {
  plot_data <- .dnmb_dbapis_plot_data(genbank_table, output_dir)
  has_outputs <- .dnmb_antidefense_module_has_outputs(output_dir, "dbAPIS")
  if (is.null(plot_data)) {
    if (!has_outputs) {
      return(NULL)
    }
    return(.dnmb_antidefense_status_plot(
      module_name = "dbAPIS",
      output_dir = output_dir,
      status_text = "Module completed, but no anti-defense hits were detected for this genome."
    ))
  }
  tbl <- plot_data$hits
  contigs <- plot_data$contigs

  summary_tbl <- tbl |>
    dplyr::group_by(.data$module_group) |>
    dplyr::summarise(
      n_profiles = dplyr::n(),
      n_genes = dplyr::n_distinct(.data$locus_tag),
      n_families = dplyr::n_distinct(.data$module_family_id),
      best_evalue = min(.data$module_evalue, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_profiles), .data$module_group)
  summary_tbl$best_evalue[!is.finite(summary_tbl$best_evalue)] <- NA_real_
  summary_tbl$group_display <- .dnmb_antidefense_target_label(
    summary_tbl$module_group, max_len = 34L
  )
  group_levels <- as.character(summary_tbl$module_group)
  display_levels <- as.character(summary_tbl$group_display)
  summary_tbl$group_display <- factor(
    summary_tbl$group_display, levels = rev(display_levels)
  )
  palette <- .dnmb_antidefense_display_palette(
    group_levels,
    colors = c("#DDF6D2", "#BFEAAB", "#8FD477", "#59BB50", "#2F9D31")
  )
  legend_labels <- stats::setNames(display_levels, group_levels)
  inventory_detail <- paste0(
    summary_tbl$n_profiles, " retained profiles / ",
    summary_tbl$n_genes, " genes / ", summary_tbl$n_families, " families | best i-E ",
    ifelse(
      is.na(summary_tbl$best_evalue), "NA",
      format(signif(summary_tbl$best_evalue, 2), scientific = TRUE)
    )
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
    summary_tbl,
    ggplot2::aes(x = .data$n_profiles, y = .data$group_display)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$module_group),
      width = 0.62, color = "grey35", linewidth = 0.24,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$n_profiles + 0.18, label = inventory_detail),
      hjust = 0, size = 2.45, color = "grey25",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.65))
    ) +
    ggplot2::labs(x = "Retained family-profile matches", y = NULL) +
    base_theme +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7.6),
      legend.position = "none"
    )

  layout_genes <- tbl |>
    dplyr::arrange(.data$replicon_id, .data$start, .data$module_evalue) |>
    dplyr::group_by(.data$replicon_id, .data$locus_tag) |>
    dplyr::summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      midpoint = (start + end) / 2,
      gene = dplyr::first(.data$gene),
      product = dplyr::first(.data$product),
      module_group = dplyr::first(.data$module_group),
      n_profiles = dplyr::n(),
      priority = max(-log10(.data$module_evalue), na.rm = TRUE),
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
  empty_contigs$no_call_label <- rep("0 retained dbAPIS matches", nrow(empty_contigs))

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
        ymin = 0.145, ymax = 0.215, fill = .data$module_group
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
      size = 2.15, lineheight = 0.88, color = "grey20",
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
      values = palette, breaks = group_levels,
      labels = unname(legend_labels[group_levels]), name = "Defense target"
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

  tbl$feature_id <- ifelse(
    !is.na(tbl$gene) & nzchar(trimws(tbl$gene)),
    paste0(trimws(tbl$gene), " | ", tbl$locus_tag),
    tbl$locus_tag
  )
  tbl$replicon_short <- contigs$replicon_short[
    match(tbl$replicon_id, contigs$replicon_id)
  ]
  replicon_prefix <- ifelse(
    grepl("^Plasmid", tbl$replicon_short),
    paste0("[", tbl$replicon_short, "] "),
    ""
  )
  tbl$profile_label <- paste0(
    replicon_prefix, tbl$feature_id, " · ", tbl$module_display
  )
  tbl <- tbl[
    order(tbl$module_evalue, -tbl$module_hmm_coverage, tbl$locus_tag, tbl$module_family_id),
    , drop = FALSE
  ]
  tbl$profile_label <- factor(tbl$profile_label, levels = rev(tbl$profile_label))
  product_short <- ifelse(
    nchar(tbl$product) > 35,
    paste0(substr(tbl$product, 1, 32), "..."),
    tbl$product
  )
  ie_label <- ifelse(
    is.finite(tbl$module_evalue),
    format(signif(tbl$module_evalue, 2), scientific = TRUE),
    "NA"
  )
  group_short <- .dnmb_antidefense_target_label(tbl$module_group, max_len = 22L)
  tbl$evidence_label <- paste0(
    product_short, " | ", group_short,
    " | query cov ", sprintf("%.2f", tbl$module_query_coverage),
    " | i-E ", ie_label,
    " | score ", sprintf("%.1f", tbl$module_hit_score)
  )
  annotation_x <- 1.05
  table_xmax <- 2.72
  p_hits <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = .data$module_hmm_coverage, y = .data$profile_label)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$module_group),
      width = 0.64, color = "grey35", linewidth = 0.22
    ) +
    ggplot2::geom_vline(
      xintercept = 0.35, linetype = "dashed",
      linewidth = 0.55, color = "#333333"
    ) +
    ggplot2::geom_vline(
      xintercept = annotation_x - 0.025, linewidth = 0.35, color = "grey82"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = annotation_x, label = .data$evidence_label),
      hjust = 0, size = 1.67, color = "grey27",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.35, 0.50, 0.75, 1.00),
      limits = c(0, table_xmax),
      expand = ggplot2::expansion(mult = c(0, 0.006))
    ) +
    ggplot2::labs(
      x = "HMM coverage | retained threshold: HMM coverage ≥ 0.35 and i-E ≤ 1e-5",
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
      legend.title = ggplot2::element_text(face = "bold", size = 8.8),
      legend.text = ggplot2::element_text(size = 7.7),
      legend.key.width = grid::unit(0.90, "lines"),
      legend.key.height = grid::unit(0.76, "lines"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))
  common_legend <- .dnmb_with_plot_pdf_device(
    cowplot::get_legend(legend_plot), width = 10, height = 1
  )
  legend_panel <- cowplot::ggdraw() + cowplot::draw_grob(
    common_legend, x = 0.5, y = 0.5, width = 0.99, height = 0.96,
    hjust = 0.5, vjust = 0.5
  )

  inventory_h <- max(1.45, 0.24 * nrow(summary_tbl) + 0.62)
  layout_h <- max(
    2.05,
    0.58 * nrow(contigs) + 0.17 * (max_label_tier + 1) + 0.60
  )
  calls_h <- max(3.2, 0.125 * nrow(tbl) + 0.95)
  legend_h <- 0.52
  total_h <- inventory_h + layout_h + calls_h + legend_h
  composite <- .dnmb_with_plot_pdf_device(
    cowplot::plot_grid(
      .dnmb_antidefense_panel_header(
        p_inventory, "A", "dbAPIS retained-evidence inventory",
        y = 0.975, plot_y = 0.04, plot_h = 0.88
      ),
      .dnmb_antidefense_panel_header(
        p_layout, "B", "dbAPIS genes across chromosome and plasmids",
        y = 0.975, plot_y = 0.03, plot_h = 0.90
      ),
      .dnmb_antidefense_panel_header(
        p_hits, "C", "All retained dbAPIS family-profile matches",
        y = 0.975, plot_y = 0.025, plot_h = 0.925
      ),
      legend_panel,
      ncol = 1,
      rel_heights = c(inventory_h, layout_h, calls_h, legend_h)
    ),
    width = 10, height = total_h
  )
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "dbAPIS_overview.pdf")
  .dnmb_module_plot_save(
    composite, pdf_path, width = 10,
    height = total_h
  )
  list(
    pdf = pdf_path,
    retained_profiles = nrow(tbl),
    retained_genes = dplyr::n_distinct(tbl$locus_tag),
    retained_families = dplyr::n_distinct(tbl$module_family_id),
    replicons_shown = nrow(contigs)
  )
}

.dnmb_plot_acrfinder_module <- function(genbank_table, output_dir) {
  .dnmb_plot_antidefense_module(genbank_table, output_dir, module_name = "AcrFinder")
}
