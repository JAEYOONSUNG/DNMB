#' Plot DNMB mobileome distribution
#'
#' Build ggplot summaries from a DNMB mobileome output directory so the merged
#' calls can be visually compared against the sequence-supported
#' ISEScan-like detections and the predicted targetable regions.
#'
#' @param output_dir DNMB mobileome result directory.
#' @param genbank Optional GenBank path. If `NULL`, the function tries to infer
#'   it from `mobileome_summary.tsv`.
#' @param prefix Prefix used for the generated plot filenames.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#'
#' @return An invisible list containing the plot objects and saved file paths.
#' @export
plot_DNMB_mobileome_distribution <- function(
  output_dir,
  genbank = NULL,
  prefix = "mobileome_distribution",
  width = 15,
  height = 10,
  min_contig_fraction = 0.12
) {
  if (missing(output_dir) || is.null(output_dir) || !dir.exists(output_dir)) {
    stop("`output_dir` must be an existing DNMB mobileome result directory.")
  }

  merged_path <- file.path(output_dir, "is_elements.tsv")
  sequence_path <- file.path(output_dir, "sequence_is_elements.tsv")
  target_sites_path <- file.path(output_dir, "is_target_sites.tsv")
  targetable_path <- file.path(output_dir, "is_targetable_regions.tsv")
  target_models_path <- file.path(output_dir, "is_target_models.tsv")
  summary_path <- file.path(output_dir, "mobileome_summary.tsv")

  if (!file.exists(merged_path)) {
    stop("Missing required file: ", merged_path)
  }
  if (!file.exists(sequence_path)) {
    stop("Missing required file: ", sequence_path)
  }

  merged <- utils::read.delim(merged_path, check.names = FALSE)
  sequence <- utils::read.delim(sequence_path, check.names = FALSE)
  target_sites <- if (file.exists(target_sites_path)) utils::read.delim(target_sites_path, check.names = FALSE) else data.frame()
  targetable <- if (file.exists(targetable_path)) utils::read.delim(targetable_path, check.names = FALSE) else data.frame()
  target_models <- if (file.exists(target_models_path)) utils::read.delim(target_models_path, check.names = FALSE) else data.frame()
  merged <- .dnmb_classify_element_evidence(merged)

  if (is.null(genbank) && file.exists(summary_path)) {
    summary_tbl <- utils::read.delim(summary_path, check.names = FALSE)
    if (all(c("metric", "value") %in% colnames(summary_tbl))) {
      genbank <- summary_tbl$value[summary_tbl$metric == "source_genbank"][1]
      if (!length(genbank) || is.na(genbank) || !nzchar(genbank) || !file.exists(genbank)) {
        genbank <- NULL
      }
    }
  }

  contigs <- .dnmb_distribution_contigs(
    genbank = genbank,
    merged = merged,
    sequence = sequence,
    target_sites = target_sites,
    targetable = targetable
  )

  merged_track <- .dnmb_prepare_distribution_track(
    df = merged,
    contigs = contigs,
    track = "Merged calls",
    family_col = "element_family",
    start_col = "start",
    end_col = "end",
    extra_cols = c(
      "evidence_class",
      "primary_support_group",
      "integrated_evidence_label",
      "annotation_supported",
      "sequence_supported",
      "native_tir_supported",
      "native_tsd_supported"
    )
  )
  sequence_track <- .dnmb_prepare_distribution_track(
    df = sequence,
    contigs = contigs,
    track = "Sequence-supported",
    family_col = "element_family",
    start_col = "start",
    end_col = "end",
    extra_cols = c("tir_found", "tsd_found")
  )
  targetable_track <- if (nrow(targetable)) {
    .dnmb_prepare_distribution_track(
      df = targetable,
      contigs = contigs,
      track = "Targetable regions",
      family_col = "target_family",
      start_col = "region_start",
      end_col = "region_end"
    )
  } else {
    data.frame()
  }

  plot_df <- dplyr::bind_rows(merged_track, sequence_track, targetable_track)
  family_levels <- plot_df %>%
    dplyr::count(.data$family, sort = TRUE) %>%
    dplyr::pull(.data$family)
  family_palette <- stats::setNames(grDevices::hcl.colors(length(family_levels), "Dark 3"), rev(family_levels))
  recognition_lookup <- .dnmb_build_recognition_lookup(
    target_models = target_models,
    sequence_elements = sequence,
    family_levels = family_levels
  )
  contigs <- contigs %>%
    dplyr::mutate(
      display_bp = pmax(
        .data$sequence_length_bp,
        max(.data$sequence_length_bp, na.rm = TRUE) * min_contig_fraction
      )
    )
  contig_widths <- contigs$display_bp / sum(contigs$display_bp)

  unified_body <- .dnmb_build_unified_mobileome_body(
    merged_track = merged_track,
    sequence_track = sequence_track,
    targetable_track = targetable_track,
    contigs = contigs,
    recognition_lookup = recognition_lookup,
    family_levels = rev(family_levels),
    family_palette = family_palette
  )

  title_grob <- grid::textGrob(
    label = "DNMB Mobileome Distribution",
    x = 0,
    hjust = 0,
    gp = grid::gpar(fontsize = 16, fontface = "bold")
  )
  subtitle_grob <- grid::textGrob(
    label = paste0(
      "Single panel = integrated merged calls plus targetable-region underlay\n",
      "Wide translucent family-colored band = targetable region | hollow ring = annotation support | filled dot = sequence support | black outer ring = native TSD/TIR support\n",
      "Merged = ", nrow(merged_track),
      " | Sequence-supported raw calls = ", nrow(sequence_track),
      if (nrow(targetable_track)) paste0(" | Targetable regions = ", nrow(targetable_track)) else ""
    ),
    x = 0,
    hjust = 0,
    gp = grid::gpar(fontsize = 11)
  )
  legend_grob <- .dnmb_build_distribution_legend_grob()
  overview_plot <- gridExtra::arrangeGrob(
    grobs = list(title_grob, subtitle_grob, legend_grob, unified_body),
    ncol = 1,
    heights = grid::unit.c(
      grid::unit(0.35, "in"),
      grid::unit(0.42, "in"),
      grid::unit(0.65, "in"),
      grid::unit(1, "null")
    )
  )
  overview_core <- overview_plot
  summary_side_plot <- NULL

  count_df <- dplyr::bind_rows(
    merged_track %>% dplyr::transmute(track = "Merged calls", family = .data$family),
    sequence_track %>% dplyr::transmute(track = "Sequence-supported", family = .data$family)
  ) %>%
    dplyr::count(.data$track, .data$family, name = "n")
  count_df$family <- factor(count_df$family, levels = rev(family_levels))

  count_plot <- ggplot2::ggplot(count_df, ggplot2::aes(x = .data$n, y = .data$family, fill = .data$track)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$n),
      position = ggplot2::position_dodge(width = 0.9),
      hjust = -0.15,
      size = 3
    ) +
    ggplot2::labs(
      title = "IS Family Counts",
      subtitle = "Compare merged DNMB calls with sequence-supported ISEScan-like calls",
      x = "Count",
      y = "IS family",
      fill = "Track"
    ) +
    ggplot2::scale_fill_manual(values = c("Merged calls" = "#475569", "Sequence-supported" = "#F59E0B")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  files <- list(
    overview_png = file.path(output_dir, paste0(prefix, "_overview.png")),
    overview_pdf = file.path(output_dir, paste0(prefix, "_overview.pdf")),
    counts_png = file.path(output_dir, paste0(prefix, "_family_counts.png")),
    counts_pdf = file.path(output_dir, paste0(prefix, "_family_counts.pdf"))
  )

  ggplot2::ggsave(files$overview_png, overview_plot, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(files$overview_pdf, overview_plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(files$counts_png, count_plot, width = width * 0.8, height = height * 0.8, dpi = 300, bg = "white")
  ggplot2::ggsave(files$counts_pdf, count_plot, width = width * 0.8, height = height * 0.8, bg = "white")

  invisible(list(
    overview_plot = overview_plot,
    overview_core = overview_core,
    summary_side_plot = summary_side_plot,
    count_plot = count_plot,
    files = files
  ))
}

.dnmb_distribution_contigs <- function(genbank = NULL, merged, sequence, target_sites = data.frame(), targetable = data.frame()) {
  contigs <- NULL
  if (!is.null(genbank) && file.exists(genbank)) {
    parsed <- .dnmb_parse_genbank_features(genbank)
    contigs <- parsed$metadata %>%
      dplyr::select(contig, contig_number, sequence_length_bp, definition)
  }

  if (is.null(contigs) || !nrow(contigs)) {
    all_positions <- dplyr::bind_rows(
      merged %>% dplyr::transmute(contig = .data$contig, end = .data$end),
      sequence %>% dplyr::transmute(contig = .data$contig, end = .data$end),
      if (nrow(target_sites)) target_sites %>% dplyr::transmute(contig = .data$contig, end = .data$site_end) else data.frame(),
      if (nrow(targetable)) targetable %>% dplyr::transmute(contig = .data$contig, end = .data$region_end) else data.frame()
    )
    contigs <- all_positions %>%
      dplyr::group_by(.data$contig) %>%
      dplyr::summarise(sequence_length_bp = max(.data$end, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(contig_number = dplyr::row_number(), definition = NA_character_)
  }

  contigs %>%
    dplyr::arrange(.data$contig_number, .data$contig) %>%
    dplyr::mutate(
      contig_label = vapply(
        seq_len(dplyr::n()),
        function(i) .dnmb_pretty_contig_label(contig[[i]], definition[[i]], contig_number[[i]]),
        character(1)
      )
    )
}

.dnmb_prepare_distribution_track <- function(df, contigs, track, family_col, start_col, end_col, extra_cols = NULL) {
  if (is.null(df) || !nrow(df)) {
    return(data.frame())
  }

  selected_cols <- unique(c("contig", family_col, start_col, end_col, extra_cols))
  selected_cols <- selected_cols[selected_cols %in% colnames(df)]
  track_df <- tibble::as_tibble(df[, selected_cols, drop = FALSE])
  colnames(track_df)[match(family_col, colnames(track_df))] <- "family"
  colnames(track_df)[match(start_col, colnames(track_df))] <- "start"
  colnames(track_df)[match(end_col, colnames(track_df))] <- "end"

  track_df <- track_df %>%
    dplyr::mutate(
      family = dplyr::na_if(as.character(.data$family), ""),
      family = dplyr::coalesce(.data$family, "Unclassified")
    ) %>%
    dplyr::mutate(
      contig = factor(.data$contig, levels = contigs$contig),
      track = track
    )

  keep_cols <- c("track", "contig", "family", "start", "end")
  keep_cols <- c(keep_cols, setdiff(colnames(track_df), keep_cols))
  track_df <- track_df %>%
    dplyr::select(dplyr::all_of(keep_cols))

  track_df
}

.dnmb_build_unified_mobileome_body <- function(
  merged_track,
  sequence_track,
  targetable_track,
  contigs,
  recognition_lookup,
  family_levels,
  family_palette
) {
  family_index <- stats::setNames(seq_along(family_levels), family_levels)
  n_family <- length(family_levels)

  main_region_width <- 78
  side_region_width <- 22
  contig_gap <- 1.2
  right_gap <- 1.4
  side_inner_gap <- 1.2
  contig_layout <- contigs %>%
    dplyr::mutate(
      width_prop = .data$display_bp / sum(.data$display_bp, na.rm = TRUE)
    )
  usable_main_width <- main_region_width - contig_gap * (nrow(contig_layout) - 1)
  contig_layout$panel_width <- contig_layout$width_prop * usable_main_width
  contig_layout$x0 <- cumsum(c(0, head(contig_layout$panel_width + contig_gap, -1)))
  contig_layout$x1 <- contig_layout$x0 + contig_layout$panel_width

  recognition_slots <- max(
    4L,
    vapply(recognition_lookup$recognition_seqs, function(x) {
      if (is.null(x) || !length(x)) return(0L)
      lens <- unique(nchar(x))
      lens <- lens[is.finite(lens) & lens > 0]
      if (!length(lens)) 0L else lens[[1]]
    }, integer(1)),
    na.rm = TRUE
  )
  side_x0 <- main_region_width + right_gap
  rec_share <- 0.95 / (0.95 + 1.95)
  recognition_width <- side_region_width * rec_share
  count_width <- side_region_width - recognition_width - side_inner_gap
  recog_x0 <- side_x0
  recog_x1 <- recog_x0 + recognition_width
  count_x0 <- recog_x1 + side_inner_gap
  count_x1 <- count_x0 + count_width
  count_centers <- seq(count_x0 + 0.7, count_x1 - 0.7, length.out = 5)
  count_x1 <- max(count_centers) + 1.0
  total_xmax <- max(main_region_width, count_x1)

  layout_map <- contig_layout %>%
    dplyr::select(contig, sequence_length_bp, x0, x1)

  map_track <- function(df) {
    if (is.null(df) || !nrow(df)) {
      return(df)
    }
    out <- df %>%
      dplyr::left_join(layout_map, by = "contig") %>%
      dplyr::mutate(
        family_y = unname(family_index[as.character(.data$family)]),
        panel_width = .data$x1 - .data$x0,
        start_prop = dplyr::if_else(.data$sequence_length_bp > 1, (.data$start - 1) / (.data$sequence_length_bp - 1), 0.5),
        end_prop = dplyr::if_else(.data$sequence_length_bp > 1, (.data$end - 1) / (.data$sequence_length_bp - 1), 0.5),
        start_x = .data$x0 + .data$panel_width * .data$start_prop,
        end_x = .data$x0 + .data$panel_width * .data$end_prop,
        mid_x = (.data$start_x + .data$end_x) / 2
      )
    out
  }

  merged_draw <- map_track(merged_track)
  target_draw <- map_track(targetable_track)
  if (nrow(merged_draw)) {
    merged_draw <- merged_draw %>%
      dplyr::mutate(
        annotation_supported = dplyr::coalesce(as.logical(.data$annotation_supported), FALSE),
        sequence_supported = dplyr::coalesce(as.logical(.data$sequence_supported), FALSE),
        native_supported = dplyr::coalesce(as.logical(.data$native_tir_supported), FALSE) |
          dplyr::coalesce(as.logical(.data$native_tsd_supported), FALSE)
      )
  }

  count_long <- merged_draw %>%
    dplyr::group_by(.data$family) %>%
    dplyr::summarise(
      TOTAL = dplyr::n(),
      ANN = sum(.data$annotation_supported, na.rm = TRUE),
      SEQ = sum(.data$sequence_supported, na.rm = TRUE),
      OVL = sum(.data$annotation_supported & .data$sequence_supported, na.rm = TRUE),
      NAT = sum(.data$native_supported, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c("TOTAL", "ANN", "SEQ", "OVL", "NAT"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      family_y = unname(family_index[as.character(.data$family)]),
      metric = factor(.data$metric, levels = c("TOTAL", "ANN", "SEQ", "OVL", "NAT"))
    )
  metric_layout <- tibble::tibble(
    metric = factor(c("TOTAL", "ANN", "SEQ", "OVL", "NAT"), levels = c("TOTAL", "ANN", "SEQ", "OVL", "NAT")),
    x_center = count_centers
  )
  metric_max <- count_long %>%
    dplyr::group_by(.data$metric) %>%
    dplyr::summarise(metric_max = max(.data$value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(metric_max = pmax(.data$metric_max, 1))
  count_draw <- count_long %>%
    dplyr::left_join(metric_layout, by = "metric") %>%
    dplyr::left_join(metric_max, by = "metric") %>%
    dplyr::mutate(
      x_bg_start = .data$x_center - 0.48,
      x_bg_end = .data$x_center + 0.48,
      x_start = .data$x_bg_start,
      x_end = .data$x_start + (0.96 * .data$value / .data$metric_max),
      x_label = .data$x_center + 0.60,
      label = dplyr::if_else(.data$value > 0, as.character(.data$value), NA_character_)
    )

  tick_df <- dplyr::bind_rows(lapply(seq_len(nrow(contig_layout)), function(i) {
    row <- contig_layout[i, , drop = FALSE]
    brks <- unique(round(pretty(c(1, row$sequence_length_bp[[1]]), n = 3)))
    brks <- brks[brks >= 1 & brks <= row$sequence_length_bp[[1]]]
    if (!length(brks)) {
      return(tibble::tibble())
    }
    x <- row$x0[[1]] + (brks - 1) / max(1, row$sequence_length_bp[[1]] - 1) * (row$x1[[1]] - row$x0[[1]])
    tibble::tibble(
      x = x,
      label = format(brks, big.mark = ",", scientific = FALSE),
      contig_mid = (row$x0[[1]] + row$x1[[1]]) / 2,
      contig_label = row$contig_label[[1]]
    )
  }))
  contig_label_df <- contig_layout %>%
    dplyr::transmute(
      x = (.data$x0 + .data$x1) / 2,
      label = .data$contig_label
    )
  count_tick_df <- tibble::tibble(
    x = count_centers,
    label = c("TOTAL", "ANN", "SEQ", "OVL", "NAT")
  )

  recog_tbl <- recognition_lookup %>%
    dplyr::filter(.data$family %in% family_levels) %>%
    dplyr::mutate(family_y = unname(family_index[as.character(.data$family)]))
  recog_layers <- .dnmb_unified_recognition_layers(
    recog_tbl = recog_tbl,
    recognition_x0 = recog_x0,
    recognition_x1 = recog_x1,
    slot_count = recognition_slots
  )

  metric_colors <- c(
    "TOTAL" = "#0F172A",
    "ANN" = "#F59E0B",
    "SEQ" = "#2563EB",
    "OVL" = "#0EA5A4",
    "NAT" = "#111827"
  )
  panel_border_color <- "grey55"
  panel_border_lwd <- 0.7

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = contig_layout$x0, xmax = contig_layout$x1, ymin = 0.5, ymax = n_family + 0.5, fill = "grey98", color = panel_border_color, linewidth = panel_border_lwd) +
    ggplot2::annotate("rect", xmin = recog_x0, xmax = recog_x1, ymin = 0.5, ymax = n_family + 0.5, fill = "white", color = panel_border_color, linewidth = panel_border_lwd) +
    ggplot2::annotate("rect", xmin = count_x0, xmax = count_x1, ymin = 0.5, ymax = n_family + 0.5, fill = "white", color = panel_border_color, linewidth = panel_border_lwd) +
    ggplot2::geom_hline(yintercept = seq_len(n_family), linewidth = 0.25, color = "grey90") +
    ggplot2::geom_vline(
      xintercept = unlist(lapply(seq_len(nrow(contig_layout)), function(i) {
        row <- contig_layout[i, , drop = FALSE]
        seq(row$x0[[1]], row$x1[[1]], length.out = 4)
      })),
      linewidth = 0.2,
      color = "grey92"
    ) +
    ggplot2::geom_vline(
      xintercept = seq(recog_x0 + 0.5, recog_x1 - 0.5, by = 1),
      linewidth = 0.2,
      color = "grey93"
    ) +
    ggplot2::geom_vline(
      xintercept = count_centers,
      linewidth = 0.2,
      color = "grey93"
    )

  if (nrow(target_draw)) {
    p <- p +
      ggplot2::geom_segment(
        data = target_draw,
        ggplot2::aes(x = .data$start_x, xend = .data$end_x, y = .data$family_y, yend = .data$family_y, color = .data$family),
        linewidth = 5.0,
        alpha = 0.22,
        lineend = "round",
        show.legend = FALSE
      )
  }

  if (nrow(merged_draw)) {
    p <- p +
      ggplot2::geom_point(
        data = merged_draw %>% dplyr::filter(.data$annotation_supported),
        ggplot2::aes(x = .data$mid_x, y = .data$family_y, color = .data$family),
        shape = 21,
        fill = "white",
        stroke = 0.95,
        size = 3.0,
        alpha = 0.95,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = merged_draw %>% dplyr::filter(.data$sequence_supported),
        ggplot2::aes(x = .data$mid_x, y = .data$family_y, color = .data$family),
        shape = 16,
        size = 1.9,
        alpha = 0.95,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = merged_draw %>% dplyr::filter(.data$native_supported),
        ggplot2::aes(x = .data$mid_x, y = .data$family_y),
        shape = 21,
        fill = NA,
        color = "black",
        stroke = 0.9,
        size = 3.6,
        show.legend = FALSE
      )
  }

  if (nrow(count_draw)) {
    p <- p +
      ggplot2::geom_segment(
        data = count_draw,
        ggplot2::aes(x = .data$x_bg_start, xend = .data$x_bg_end, y = .data$family_y, yend = .data$family_y),
        linewidth = 4.6,
        lineend = "butt",
        color = "grey92",
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = count_draw,
        ggplot2::aes(x = .data$x_start, xend = .data$x_end, y = .data$family_y, yend = .data$family_y, color = .data$metric),
        linewidth = 4.6,
        lineend = "butt",
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = count_draw,
        ggplot2::aes(x = .data$x_label, y = .data$family_y, label = .data$label),
        hjust = 0,
        size = 2.7,
        show.legend = FALSE
      )
  }

  p <- p +
    ggplot2::geom_text(
      data = contig_label_df,
      ggplot2::aes(x = .data$x, y = n_family + 1.05, label = .data$label),
      fontface = "bold",
      size = 3.0,
      inherit.aes = FALSE
    ) +
    ggplot2::annotate("text", x = (recog_x0 + recog_x1) / 2, y = n_family + 1.05, label = "Recognition", fontface = "bold", size = 3.0) +
    ggplot2::annotate("text", x = (count_x0 + count_x1) / 2, y = n_family + 1.05, label = "Integrated evidence counts", fontface = "bold", size = 3.0) +
    ggplot2::geom_segment(
      data = tick_df,
      ggplot2::aes(x = .data$x, xend = .data$x, y = 0.5, yend = 0.32),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      data = count_tick_df,
      ggplot2::aes(x = .data$x, xend = .data$x, y = 0.5, yend = 0.32),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      data = tick_df,
      ggplot2::aes(x = .data$x, y = 0.12, label = .data$label),
      inherit.aes = FALSE,
      size = 2.3
    ) +
    ggplot2::geom_text(
      data = count_tick_df,
      ggplot2::aes(x = .data$x, y = 0.12, label = .data$label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 2.2
    ) +
    ggplot2::annotate("text", x = mean(range(contig_layout$x0, contig_layout$x1)), y = -0.35, label = "Contig coordinate", size = 4) +
    ggplot2::scale_color_manual(values = c(family_palette, metric_colors), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(min(contig_layout$x0) - 0.5, total_xmax + 0.8), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(
      limits = c(-0.5, n_family + 1.4),
      breaks = seq_along(family_levels),
      labels = family_levels,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(5.5, 8, 16, 5.5)
    ) +
    ggplot2::ylab("IS family") +
    ggplot2::coord_cartesian(clip = "off")

  for (layer in recog_layers) {
    p <- p + layer
  }

  p
}

.dnmb_unified_recognition_layers <- function(recog_tbl, recognition_x0, recognition_x1, slot_count) {
  if (is.null(recog_tbl) || !nrow(recog_tbl)) {
    return(list())
  }
  slot_width <- (recognition_x1 - recognition_x0) / slot_count
  slot_centers <- recognition_x0 + ((seq_len(slot_count) - 0.5) * slot_width)
  layers <- list()
  idx <- 1L
  for (i in seq_len(nrow(recog_tbl))) {
    row <- recog_tbl[i, , drop = FALSE]
    seqs <- row$recognition_seqs[[1]]
    label <- row$recognition_label[[1]]
    motif_len <- .dnmb_recognition_motif_len(seqs = seqs, label = label, slot_count = slot_count)
    start_slot <- floor((slot_count - motif_len) / 2) + 1L
    end_slot <- start_slot + motif_len - 1L
    xmin <- slot_centers[[start_slot]] - (slot_width * 0.46)
    xmax <- slot_centers[[end_slot]] + (slot_width * 0.46)
    layers[[idx]] <- ggplot2::annotation_custom(
      grob = .dnmb_recognition_logo_grob(seqs = seqs, label = label),
      xmin = xmin,
      xmax = xmax,
      ymin = row$family_y[[1]] - 0.34,
      ymax = row$family_y[[1]] + 0.34
    )
    idx <- idx + 1L
  }
  layers
}

.dnmb_build_contig_panel <- function(
  track_df,
  targetable_df = NULL,
  contig_name,
  contig_length,
  family_levels,
  family_palette,
  show_y = FALSE,
  show_x = TRUE,
  panel_title = NULL
) {
  panel_df <- track_df %>%
    dplyr::filter(.data$contig == .env$contig_name) %>%
    dplyr::mutate(
      family = factor(.data$family, levels = family_levels),
      midpoint = (.data$start + .data$end) / 2
    )
  target_df <- if (!is.null(targetable_df) && nrow(targetable_df)) {
    targetable_df %>%
      dplyr::filter(.data$contig == .env$contig_name) %>%
      dplyr::mutate(family = factor(.data$family, levels = family_levels))
  } else {
    data.frame()
  }

  track_name <- if ("track" %in% colnames(panel_df) && nrow(panel_df)) unique(as.character(panel_df$track))[1] else NA_character_
  has_merged_support <- all(c("annotation_supported", "sequence_supported") %in% colnames(panel_df))
  has_sequence_native <- any(c("tir_found", "tsd_found") %in% colnames(panel_df))

  if (has_merged_support) {
    panel_df <- panel_df %>%
      dplyr::mutate(
        annotation_supported = dplyr::coalesce(as.logical(.data$annotation_supported), FALSE),
        sequence_supported = dplyr::coalesce(as.logical(.data$sequence_supported), FALSE),
        native_supported = dplyr::coalesce(as.logical(.data$native_tir_supported), FALSE) |
          dplyr::coalesce(as.logical(.data$native_tsd_supported), FALSE)
      )
  } else if (has_sequence_native) {
    panel_df <- panel_df %>%
      dplyr::mutate(
        native_supported = dplyr::coalesce(as.logical(.data$tir_found), FALSE) |
          dplyr::coalesce(as.logical(.data$tsd_found), FALSE)
      )
  }

  p <- ggplot2::ggplot(panel_df)
  if (identical(track_name, "Merged calls") && has_merged_support) {
    p <- p +
      {
        if (nrow(target_df)) {
          ggplot2::geom_segment(
            data = target_df,
            ggplot2::aes(
              x = .data$start,
              xend = .data$end,
              y = .data$family,
              yend = .data$family,
              color = .data$family
            ),
            inherit.aes = FALSE,
            alpha = 0.22,
            linewidth = 5.2,
            lineend = "round"
          )
        }
      } +
      ggplot2::geom_point(
        data = panel_df %>% dplyr::filter(.data$annotation_supported),
        ggplot2::aes(x = .data$midpoint, y = .data$family, color = .data$family),
        shape = 21,
        fill = "white",
        stroke = 0.95,
        size = 3.0,
        alpha = 0.95
      ) +
      ggplot2::geom_point(
        data = panel_df %>% dplyr::filter(.data$sequence_supported),
        ggplot2::aes(x = .data$midpoint, y = .data$family, color = .data$family),
        shape = 16,
        size = 1.95,
        alpha = 0.95,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = panel_df %>% dplyr::filter(.data$native_supported),
        ggplot2::aes(x = .data$midpoint, y = .data$family),
        shape = 21,
        fill = NA,
        color = "black",
        stroke = 0.9,
        size = 3.6
      )
  } else if (identical(track_name, "Sequence-supported")) {
    if (!"native_supported" %in% colnames(panel_df)) {
      panel_df$native_supported <- FALSE
    }
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$start,
          xend = .data$end,
          y = .data$family,
          yend = .data$family,
          color = .data$family
        ),
        alpha = 0.15,
        linewidth = 0.6,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$midpoint, y = .data$family, color = .data$family),
        shape = 16,
        size = 2.1,
        alpha = 0.92,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = panel_df %>% dplyr::filter(.data$native_supported),
        ggplot2::aes(x = .data$midpoint, y = .data$family),
        shape = 21,
        fill = NA,
        color = "black",
        stroke = 0.85,
        size = 3.5
      )
  } else {
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$start,
          xend = .data$end,
          y = .data$family,
          yend = .data$family,
          color = .data$family
        ),
        linewidth = 2,
        alpha = 0.85,
        lineend = "round"
      )
  }

  p <- p +
    ggplot2::scale_y_discrete(
      limits = family_levels,
      drop = FALSE,
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::scale_color_manual(values = family_palette, guide = "none") +
    ggplot2::scale_fill_manual(values = family_palette, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(1, max(1, contig_length)),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = panel_title,
      x = if (isTRUE(show_x)) "Contig coordinate" else NULL,
      y = if (isTRUE(show_y)) "IS family" else NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "grey90"),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10),
      axis.text.y = if (isTRUE(show_y)) ggplot2::element_text(size = 9) else ggplot2::element_blank(),
      axis.ticks.y = if (isTRUE(show_y)) ggplot2::element_line() else ggplot2::element_blank(),
      axis.title.y = if (isTRUE(show_y)) ggplot2::element_text() else ggplot2::element_blank(),
      axis.text.x = if (isTRUE(show_x)) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
      axis.ticks.x = if (isTRUE(show_x)) ggplot2::element_line() else ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  if (!nrow(panel_df)) {
    p <- p + ggplot2::geom_blank(ggplot2::aes(x = 1, y = factor(family_levels, levels = family_levels)[1]))
  }

  p
}

.dnmb_build_track_summary_panel <- function(
  track_name,
  merged_track,
  sequence_track,
  targetable_track,
  target_sites,
  family_levels,
  recognition_lookup,
  family_palette,
  show_x = FALSE,
  show_headers = FALSE
) {
  if (identical(track_name, "Merged calls")) {
    return(.dnmb_build_integrated_side_panel(
      merged_track = merged_track,
      recognition_lookup = recognition_lookup,
      family_levels = family_levels
    ))
  }

  if (identical(track_name, "Merged calls")) {
    summary_wide <- merged_track %>%
      dplyr::group_by(.data$family) %>%
      dplyr::summarise(
        TOTAL = dplyr::n(),
        ANN = sum(dplyr::coalesce(.data$annotation_supported, FALSE), na.rm = TRUE),
        SEQ = sum(dplyr::coalesce(.data$sequence_supported, FALSE), na.rm = TRUE),
        OVL = sum(dplyr::coalesce(.data$annotation_supported, FALSE) & dplyr::coalesce(.data$sequence_supported, FALSE), na.rm = TRUE),
        NAT = sum(dplyr::coalesce(.data$native_tir_supported, FALSE) | dplyr::coalesce(.data$native_tsd_supported, FALSE), na.rm = TRUE),
        .groups = "drop"
      )

    summary_df <- summary_wide %>%
      tidyr::pivot_longer(
        cols = c("TOTAL", "ANN", "SEQ", "OVL", "NAT"),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        family = factor(.data$family, levels = family_levels),
        metric = factor(.data$metric, levels = c("TOTAL", "ANN", "SEQ", "OVL", "NAT"))
      )
    title <- paste0("Integrated evidence counts (", length(family_levels), " families)")
    metric_colors <- c(
      "TOTAL" = "#0F172A",
      "ANN" = "#F59E0B",
      "SEQ" = "#2563EB",
      "OVL" = "#0EA5A4",
      "NAT" = "#111827"
    )
    metric_layout <- tibble::tibble(
      metric = factor(c("TOTAL", "ANN", "SEQ", "OVL", "NAT"), levels = c("TOTAL", "ANN", "SEQ", "OVL", "NAT")),
      x_center = c(0.55, 1.75, 2.95, 4.15, 5.35)
    )
    metric_max <- summary_df %>%
      dplyr::group_by(.data$metric) %>%
      dplyr::summarise(metric_max = max(.data$value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(metric_max = pmax(.data$metric_max, 1))
    summary_df <- summary_df %>%
      dplyr::left_join(metric_layout, by = "metric") %>%
      dplyr::left_join(metric_max, by = "metric") %>%
      dplyr::mutate(
        x_start = .data$x_center - 0.36,
        x_end = .data$x_start + (0.72 * .data$value / .data$metric_max),
        x_label = .data$x_center + 0.41
      )

    separator_df <- tibble::tibble(x = c(1.15, 2.35, 3.55, 4.75))
    header_breaks <- metric_layout$x_center
    header_labels <- as.character(metric_layout$metric)

    bar_plot <- ggplot2::ggplot(summary_df, ggplot2::aes(y = .data$family)) +
      ggplot2::geom_vline(
        data = separator_df,
        ggplot2::aes(xintercept = .data$x),
        inherit.aes = FALSE,
        linewidth = 0.25,
        color = "grey88"
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$x_start,
          xend = .data$x_end,
          yend = .data$family,
          color = .data$metric
        ),
        linewidth = 4.2,
        lineend = "butt",
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = .data$x_label, label = .data$value),
        hjust = 0,
        size = 2.55,
        check_overlap = TRUE
      ) +
      ggplot2::labs(
        title = if (isTRUE(show_headers)) title else NULL,
        x = NULL,
        y = NULL
      ) +
      ggplot2::scale_color_manual(values = metric_colors, guide = "none") +
      ggplot2::scale_y_discrete(
        limits = family_levels,
        drop = FALSE,
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0.1, 5.95),
        breaks = header_breaks,
        labels = header_labels,
        position = "top",
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::theme_bw(base_size = 10.5) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "grey90"),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x.top = ggplot2::element_text(size = 8.2, face = "bold"),
        axis.ticks.x.top = ggplot2::element_blank(),
        axis.text.x.bottom = ggplot2::element_blank(),
        axis.ticks.x.bottom = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(5.5, 12, 5.5, 5.5)
      ) +
      ggplot2::coord_cartesian(clip = "off")
  } else {
    region_counts <- if (nrow(targetable_track)) {
      targetable_track %>%
        dplyr::count(.data$family, name = "value") %>%
        dplyr::mutate(
          family = factor(.data$family, levels = family_levels),
          metric = "Targetable regions"
        )
    } else {
      tibble::tibble(family = factor(character(), levels = family_levels), value = integer(), metric = character())
    }
    site_counts <- if (nrow(target_sites)) {
      target_sites %>%
        dplyr::mutate(family = dplyr::coalesce(dplyr::na_if(as.character(.data$target_family), ""), "Unclassified")) %>%
        dplyr::count(.data$family, name = "value") %>%
        dplyr::mutate(
          family = factor(.data$family, levels = family_levels),
          metric = "Targetable sites"
        )
    } else {
      tibble::tibble(family = factor(character(), levels = family_levels), value = integer(), metric = character())
    }
    summary_df <- dplyr::bind_rows(region_counts, site_counts) %>%
      dplyr::arrange(.data$family, .data$metric)
    title <- paste0("Targetable summary (", length(family_levels), " families)")
    bar_plot <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$value, y = .data$family, fill = .data$family))
    bar_plot <- bar_plot +
      ggplot2::geom_col(ggplot2::aes(alpha = .data$metric), position = ggplot2::position_dodge(width = 0.8), width = 0.7, show.legend = FALSE) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$value),
        position = ggplot2::position_dodge(width = 0.8),
        hjust = -0.15,
        size = 2.7,
        check_overlap = TRUE
      ) +
      ggplot2::labs(
        title = if (isTRUE(show_headers)) title else NULL,
        x = NULL,
        y = NULL
      ) +
      ggplot2::scale_y_discrete(
        limits = family_levels,
        drop = FALSE,
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
      ) +
      ggplot2::scale_fill_manual(values = family_palette[family_levels], guide = "none") +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.02, 0.10))
      ) +
      ggplot2::scale_alpha_manual(
        values = c("Targetable regions" = 0.45, "Targetable sites" = 1),
        guide = "none"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "grey90"),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = if (isTRUE(show_x)) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
        axis.ticks.x = if (isTRUE(show_x)) ggplot2::element_line() else ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.position = "none",
        plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5)
      ) +
      ggplot2::coord_cartesian(clip = "off")
  }

  if (!nrow(summary_df)) {
    summary_df <- tibble::tibble(family = factor(character(), levels = family_levels), value = integer(), metric = character())
  }

  if (!exists("bar_plot")) {
    bar_plot <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$value, y = .data$family, fill = .data$family)) +
      ggplot2::theme_void()
  }

  if (identical(track_name, "Merged calls") && !nrow(summary_df)) {
    bar_plot <- bar_plot +
      ggplot2::geom_blank(ggplot2::aes(x = 0, y = factor(family_levels, levels = family_levels)[1]))
  }

  recog_df <- recognition_lookup %>%
    dplyr::filter(.data$family %in% family_levels) %>%
    dplyr::mutate(family = factor(.data$family, levels = family_levels))
  if (!nrow(recog_df)) {
    return(bar_plot)
  }
  recog_stack <- .dnmb_build_recognition_logo_stack(
    recog_df = recog_df,
    family_levels = family_levels,
    family_palette = family_palette,
    show_title = show_headers
  )

  cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_grob(recog_stack),
    bar_plot,
    nrow = 1,
    rel_widths = c(1.15, 1.65),
    align = "v",
    axis = "lr"
  )
}

.dnmb_build_integrated_side_panel <- function(merged_track, recognition_lookup, family_levels) {
  metric_colors <- c(
    "TOTAL" = "#0F172A",
    "ANN" = "#F59E0B",
    "SEQ" = "#2563EB",
    "OVL" = "#0EA5A4",
    "NAT" = "#111827"
  )
  nt_colors <- c(
    "A" = "#16A34A",
    "C" = "#2563EB",
    "G" = "#F59E0B",
    "T" = "#DC2626"
  )

  family_df <- tibble::tibble(
    family = factor(family_levels, levels = family_levels)
  )

  slot_count <- max(
    4L,
    vapply(recognition_lookup$recognition_seqs, function(x) {
      if (is.null(x) || !length(x)) {
        return(0L)
      }
      lens <- unique(nchar(x))
      lens <- lens[is.finite(lens) & lens > 0]
      if (!length(lens)) 0L else lens[[1]]
    }, integer(1)),
    na.rm = TRUE
  )

  recog_tbl <- recognition_lookup %>%
    dplyr::filter(.data$family %in% family_levels) %>%
    dplyr::mutate(
      family = factor(.data$family, levels = family_levels)
    )

  count_wide <- merged_track %>%
    dplyr::group_by(.data$family) %>%
    dplyr::summarise(
      TOTAL = dplyr::n(),
      ANN = sum(dplyr::coalesce(.data$annotation_supported, FALSE), na.rm = TRUE),
      SEQ = sum(dplyr::coalesce(.data$sequence_supported, FALSE), na.rm = TRUE),
      OVL = sum(dplyr::coalesce(.data$annotation_supported, FALSE) & dplyr::coalesce(.data$sequence_supported, FALSE), na.rm = TRUE),
      NAT = sum(dplyr::coalesce(.data$native_tir_supported, FALSE) | dplyr::coalesce(.data$native_tsd_supported, FALSE), na.rm = TRUE),
      .groups = "drop"
    )
  count_long <- count_wide %>%
    tidyr::pivot_longer(
      cols = c("TOTAL", "ANN", "SEQ", "OVL", "NAT"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      family = factor(.data$family, levels = family_levels),
      metric = factor(.data$metric, levels = c("TOTAL", "ANN", "SEQ", "OVL", "NAT"))
    )

  metric_layout <- tibble::tibble(
    metric = factor(c("TOTAL", "ANN", "SEQ", "OVL", "NAT"), levels = c("TOTAL", "ANN", "SEQ", "OVL", "NAT")),
    x_center = slot_count + c(2.0, 3.2, 4.4, 5.6, 6.8)
  )
  metric_max <- count_long %>%
    dplyr::group_by(.data$metric) %>%
    dplyr::summarise(metric_max = max(.data$value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(metric_max = pmax(.data$metric_max, 1))
  count_draw <- count_long %>%
    dplyr::left_join(metric_layout, by = "metric") %>%
    dplyr::left_join(metric_max, by = "metric") %>%
    dplyr::mutate(
      x_start = .data$x_center - 0.32,
      x_end = .data$x_start + (0.64 * .data$value / .data$metric_max),
      x_label = .data$x_center + 0.42
    )

  count_box_xmin <- slot_count + 1.35
  count_box_xmax <- max(metric_layout$x_center) + 0.75
  recognition_box_xmin <- 0.5
  recognition_box_xmax <- slot_count + 0.5

  header_breaks <- c((recognition_box_xmin + recognition_box_xmax) / 2, metric_layout$x_center)
  header_labels <- c("Recognition", as.character(metric_layout$metric))

  p <- ggplot2::ggplot(family_df, ggplot2::aes(y = .data$family)) +
    ggplot2::annotate("rect", xmin = recognition_box_xmin, xmax = recognition_box_xmax, ymin = -Inf, ymax = Inf, fill = "grey99", color = "grey55", linewidth = 0.6) +
    ggplot2::annotate("rect", xmin = count_box_xmin, xmax = count_box_xmax, ymin = -Inf, ymax = Inf, fill = "white", color = "grey55", linewidth = 0.6) +
    ggplot2::geom_vline(
      xintercept = seq(recognition_box_xmin + 0.5, recognition_box_xmax - 0.5, by = 1),
      linewidth = 0.2,
      color = "grey92"
    ) +
    ggplot2::geom_vline(
      data = metric_layout,
      ggplot2::aes(xintercept = .data$x_center),
      inherit.aes = FALSE,
      linewidth = 0.2,
      color = "grey92"
    ) +
    ggplot2::geom_segment(
      data = count_draw,
      ggplot2::aes(x = .data$x_start, xend = .data$x_end, y = .data$family, yend = .data$family, color = .data$metric),
      inherit.aes = FALSE,
      linewidth = 4.0,
      lineend = "butt",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = count_draw,
      ggplot2::aes(x = .data$x_label, y = .data$family, label = .data$value),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.5
    ) +
    ggplot2::scale_color_manual(values = c(metric_colors, nt_colors, "LABEL" = "#475569"), guide = "none") +
    ggplot2::scale_y_discrete(
      limits = family_levels,
      drop = FALSE,
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.25, count_box_xmax + 0.7),
      breaks = header_breaks,
      labels = header_labels,
      position = "top",
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "grey90"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(size = 8.2, face = "bold"),
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.text.x.bottom = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(5.5, 8, 5.5, 5.5)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  recog_layers <- .dnmb_recognition_annotation_layers(
    recog_tbl = recog_tbl,
    family_levels = family_levels,
    slot_count = slot_count
  )
  for (layer in recog_layers) {
    p <- p + layer
  }

  p
}

.dnmb_recognition_annotation_layers <- function(recog_tbl, family_levels, slot_count) {
  layers <- list()
  layer_idx <- 1L

  for (i in seq_len(nrow(recog_tbl))) {
    row <- recog_tbl[i, , drop = FALSE]
    family <- as.character(row$family[[1]])
    y_center <- match(family, family_levels)
    if (is.na(y_center)) {
      next
    }

    seqs <- row$recognition_seqs[[1]]
    label <- row$recognition_label[[1]]
    motif_len <- .dnmb_recognition_motif_len(seqs = seqs, label = label, slot_count = slot_count)
    start_slot <- floor((slot_count - motif_len) / 2) + 1L
    end_slot <- start_slot + motif_len - 1L
    xmin <- start_slot - 0.42
    xmax <- end_slot + 0.42
    ymin <- y_center - 0.36
    ymax <- y_center + 0.36

    grob <- .dnmb_recognition_logo_grob(seqs = seqs, label = label)
    layers[[layer_idx]] <- ggplot2::annotation_custom(
      grob = grob,
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    )
    layer_idx <- layer_idx + 1L
  }

  layers
}

.dnmb_recognition_logo_grob <- function(seqs, label) {
  logo_seqs <- .dnmb_recognition_logo_sequences(seqs = seqs, label = label)

  if (length(logo_seqs) >= 1L) {
    p <- ggseqlogo::ggseqlogo(logo_seqs, method = "prob", seq_type = "dna") +
      ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5),
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
    return(ggplot2::ggplotGrob(p))
  }

  motif <- .dnmb_recognition_fallback_string(valid_seqs = logo_seqs, label = label)
  grid::textGrob(
    label = motif,
    x = 0.5,
    hjust = 0.5,
    gp = grid::gpar(fontsize = 7, col = "#475569", fontfamily = "Courier")
  )
}

.dnmb_recognition_logo_sequences <- function(seqs, label, min_reps = 6L) {
  valid_seqs <- character()
  if (!is.null(seqs) && length(seqs)) {
    valid_seqs <- toupper(seqs)
    valid_seqs <- valid_seqs[grepl("^[ACGT]+$", valid_seqs)]
    lens <- unique(nchar(valid_seqs))
    if (length(lens) > 1L) {
      valid_seqs <- valid_seqs[nchar(valid_seqs) == min(lens)]
    }
  }

  if (!length(valid_seqs)) {
    motif <- .dnmb_recognition_fallback_string(valid_seqs = character(), label = label)
    if (!is.na(motif) && nzchar(motif) && grepl("^[ACGT]+$", motif)) {
      valid_seqs <- motif
    }
  }

  if (!length(valid_seqs)) {
    return(character())
  }

  if (length(valid_seqs) == 1L) {
    return(rep(valid_seqs[[1]], max(2L, as.integer(min_reps))))
  }

  valid_seqs
}

.dnmb_recognition_fallback_string <- function(valid_seqs, label) {
  if (length(valid_seqs) == 1L) {
    return(valid_seqs[[1]])
  }

  if (!is.na(label) && nzchar(label)) {
    first_line <- strsplit(label, "\n", fixed = TRUE)[[1]][1]
    token <- gsub("[^A-Za-z0-9/+_-]", "", first_line)
    if (nzchar(token)) {
      return(token)
    }
  }

  "-"
}

.dnmb_pretty_contig_label <- function(contig, definition = NA_character_, contig_number = NA_integer_) {
  if (identical(contig, "annot")) {
    return("Chromosome")
  }
  contig
}

.dnmb_build_recognition_lookup <- function(target_models, sequence_elements, family_levels) {
  base_tbl <- tibble::tibble(family = family_levels)
  if (is.null(target_models) || !nrow(target_models)) {
    return(base_tbl %>% dplyr::mutate(recognition_label = "-", recognition_seqs = vector("list", dplyr::n())))
  }

  tir_tbl <- if (!is.null(sequence_elements) && nrow(sequence_elements)) {
    sequence_elements %>%
      dplyr::filter(isTRUE(.data$tir_found), !is.na(.data$tir_len_bp)) %>%
      dplyr::group_by(.data$element_family) %>%
      dplyr::summarise(
        tir_label = paste0("IR", suppressWarnings(max(.data$tir_len_bp, na.rm = TRUE)), "bp"),
        .groups = "drop"
      ) %>%
      dplyr::rename(family = "element_family")
  } else {
    tibble::tibble(family = character(), tir_label = character())
  }

  model_tbl <- target_models %>%
    dplyr::transmute(
      family = .data$family,
      recognition_label = dplyr::case_when(
        .data$model_type == "empirical_tsd" ~ paste0(dplyr::coalesce(.data$model_motif, "-"), "\nTSD", dplyr::coalesce(as.character(.data$dominant_tsd_len), "?")),
        .data$model_type == "motif_set" ~ paste0(vapply(.data$model_motif, .dnmb_trim_motif_label, character(1)), "\nset"),
        .data$model_type == "context_only" ~ paste0("context", ifelse(is.na(.data$dominant_tsd_len), "", paste0("\nTSD", .data$dominant_tsd_len))),
        TRUE ~ "-"
      ),
      recognition_seqs = lapply(
        seq_len(dplyr::n()),
        function(i) .dnmb_recognition_sequences_from_model(
          observed_examples = .data$observed_tsd_examples[[i]],
          dominant_tsd_len = .data$dominant_tsd_len[[i]],
          model_motif = .data$model_motif[[i]]
        )
      )
    ) %>%
    dplyr::left_join(tir_tbl, by = "family") %>%
    dplyr::mutate(
      recognition_label = dplyr::case_when(
        !is.na(.data$tir_label) & .data$recognition_label != "-" ~ paste(.data$recognition_label, .data$tir_label, sep = "\n"),
        TRUE ~ .data$recognition_label
      )
    ) %>%
    dplyr::select("family", "recognition_label", "recognition_seqs")

  base_tbl %>%
    dplyr::left_join(model_tbl, by = "family") %>%
    dplyr::mutate(
      recognition_label = dplyr::coalesce(.data$recognition_label, "-"),
      recognition_seqs = lapply(.data$recognition_seqs, function(x) if (is.null(x)) character() else x)
    )
}

.dnmb_trim_motif_label <- function(x, max_items = 2L) {
  if (is.na(x) || !nzchar(x)) {
    return("-")
  }
  vals <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  vals <- vals[nzchar(vals)]
  if (!length(vals)) {
    return("-")
  }
  vals <- vals[seq_len(min(length(vals), as.integer(max_items)))]
  paste(vals, collapse = "/")
}

.dnmb_recognition_sequences_from_model <- function(observed_examples = NA_character_, dominant_tsd_len = NA_integer_, model_motif = NA_character_) {
  seqs <- character()

  if (!is.na(observed_examples) && nzchar(observed_examples)) {
    tokens <- trimws(strsplit(observed_examples, ";", fixed = TRUE)[[1]])
    parsed <- lapply(tokens, function(tok) {
      m <- regexec("^([A-Za-z]+)\\((\\d+)\\)$", tok)
      hit <- regmatches(tok, m)[[1]]
      if (length(hit) == 3L) {
        seq <- toupper(hit[2])
        n <- as.integer(hit[3])
        if (!is.na(dominant_tsd_len) && dominant_tsd_len > 0 && nchar(seq) != dominant_tsd_len) {
          return(character())
        }
        # Preserve observed copy counts for the seqlogo instead of collapsing to unique.
        return(rep(seq, max(1L, min(n, 8L))))
      }
      character()
    })
    seqs <- unlist(parsed, use.names = FALSE)
  }

  if (!length(seqs) && !is.na(model_motif) && nzchar(model_motif)) {
    motifs <- trimws(strsplit(model_motif, ";", fixed = TRUE)[[1]])
    motifs <- toupper(motifs[nzchar(motifs)])
    motifs <- motifs[grepl("^[ACGT]+$", motifs)]
    if (!is.na(dominant_tsd_len) && dominant_tsd_len > 0) {
      motifs <- motifs[nchar(motifs) == dominant_tsd_len]
    }
    seqs <- motifs
  }

  if (length(unique(nchar(seqs))) > 1L) {
    target_len <- suppressWarnings(as.integer(stats::na.omit(dominant_tsd_len)[1]))
    if (is.na(target_len)) {
      target_len <- min(nchar(seqs))
    }
    seqs <- seqs[nchar(seqs) == target_len]
  }

  seqs
}

.dnmb_build_recognition_logo_stack <- function(recog_df, family_levels, family_palette, show_title = FALSE) {
  display_families <- rev(family_levels)
  seq_lengths <- vapply(recog_df$recognition_seqs, function(x) {
    if (is.null(x) || !length(x)) {
      return(0L)
    }
    vals <- unique(nchar(x))
    vals <- vals[is.finite(vals) & vals > 0]
    if (!length(vals)) 0L else vals[[1]]
  }, integer(1))
  slot_count <- max(4L, seq_lengths, na.rm = TRUE)

  title_grob <- if (isTRUE(show_title)) {
    grid::textGrob(
      "Recognition",
      x = 0.5,
      hjust = 0.5,
      gp = grid::gpar(fontsize = 10, fontface = "bold")
    )
  } else {
    grid::nullGrob()
  }

  logo_grobs <- lapply(display_families, function(fam) {
    row <- recog_df %>% dplyr::filter(.data$family == .env$fam)
    seqs <- if (nrow(row)) row$recognition_seqs[[1]] else character()
    label <- if (nrow(row)) row$recognition_label[[1]] else "-"
    color <- family_palette[[fam]]
    motif_len <- .dnmb_recognition_motif_len(seqs = seqs, label = label, slot_count = slot_count)

    if (length(seqs) >= 2L && length(unique(nchar(seqs))) == 1L && unique(nchar(seqs)) > 0L) {
      p <- ggseqlogo::ggseqlogo(seqs, method = "prob", seq_type = "dna") +
        ggplot2::theme_void(base_size = 9) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5),
          plot.background = ggplot2::element_rect(fill = "white", color = NA)
        )
      .dnmb_center_grob_in_slots(
        grob = ggplot2::ggplotGrob(p),
        slot_count = slot_count,
        motif_len = motif_len,
        draw_slot_guides = TRUE
      )
    } else {
      .dnmb_center_grob_in_slots(
        grob = grid::textGrob(
          label = .dnmb_compact_recognition_label(label),
          x = 0.5,
          hjust = 0.5,
          gp = grid::gpar(fontsize = 7, col = color, fontfamily = "Courier")
        ),
        slot_count = slot_count,
        motif_len = motif_len,
        draw_slot_guides = TRUE
      )
    }
  })

  gridExtra::arrangeGrob(
    grobs = c(list(title_grob), logo_grobs),
    ncol = 1,
    heights = grid::unit.c(
      grid::unit(if (isTRUE(show_title)) 0.22 else 0.01, "in"),
      grid::unit(rep(1, length(logo_grobs)), "null")
    )
  )
}

.dnmb_recognition_motif_len <- function(seqs, label, slot_count) {
  if (!is.null(seqs) && length(seqs)) {
    lens <- unique(nchar(seqs))
    lens <- lens[is.finite(lens) & lens > 0]
    if (length(lens)) {
      return(max(1L, min(slot_count, lens[[1]])))
    }
  }

  if (!is.na(label) && nzchar(label)) {
    tsd_len <- stringr::str_match(label, "TSD([0-9]+)")[, 2]
    tsd_len <- suppressWarnings(as.integer(tsd_len))
    if (length(tsd_len) && !is.na(tsd_len[[1]]) && tsd_len[[1]] > 0) {
      return(max(1L, min(slot_count, tsd_len[[1]])))
    }
  }

  max(2L, min(slot_count, 4L))
}

.dnmb_center_grob_in_slots <- function(grob, slot_count, motif_len, draw_slot_guides = FALSE) {
  motif_len <- max(1L, min(as.integer(motif_len), as.integer(slot_count)))
  start_slot <- floor((slot_count - motif_len) / 2) + 1L
  end_slot <- start_slot + motif_len - 1L

  tbl <- gtable::gtable(
    widths = grid::unit(rep(1, slot_count), "null"),
    heights = grid::unit(1, "null")
  )

  if (isTRUE(draw_slot_guides)) {
    guide_grobs <- lapply(seq_len(slot_count), function(i) {
      grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey93", lwd = 0.35))
    })
    tbl <- gtable::gtable_add_grob(
      tbl,
      grobs = guide_grobs,
      t = rep(1L, slot_count),
      l = seq_len(slot_count),
      b = rep(1L, slot_count),
      r = seq_len(slot_count),
      z = 0
    )
  }

  gtable::gtable_add_grob(
    tbl,
    grobs = grob,
    t = 1L,
    l = start_slot,
    b = 1L,
    r = end_slot,
    clip = "off",
    z = 1
  )
}

.dnmb_compact_recognition_label <- function(label) {
  if (is.na(label) || !nzchar(label)) {
    return("-")
  }
  parts <- strsplit(label, "\n", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    return("-")
  }
  if (length(parts) > 2L) {
    parts <- parts[1:2]
  }
  paste(parts, collapse = "\n")
}

.dnmb_build_distribution_legend_grob <- function() {
  support_group_colors <- c(
    "NATIVE" = "#111827",
    "OVERLAP" = "#0EA5A4",
    "ANN" = "#F59E0B",
    "SEQ" = "#2563EB"
  )

  p <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0.5, xend = 1.35, y = 1.0, yend = 1.0, linewidth = 6, alpha = 0.22, color = "#94A3B8", lineend = "round") +
    ggplot2::annotate("point", x = 2.0, y = 1.0, shape = 21, size = 3.4, stroke = 1.0, color = "#EC4899", fill = "white") +
    ggplot2::annotate("point", x = 3.25, y = 1.0, shape = 16, size = 2.4, color = "#2563EB") +
    ggplot2::annotate("point", x = 4.45, y = 1.0, shape = 21, size = 4.0, stroke = 0.95, color = "black", fill = NA) +
    ggplot2::annotate("text", x = 0.95, y = 0.45, label = "Targetable region band", size = 3.1) +
    ggplot2::annotate("text", x = 2.0, y = 0.45, label = "Annotation call", size = 3.1) +
    ggplot2::annotate("text", x = 3.25, y = 0.45, label = "Sequence call", size = 3.1) +
    ggplot2::annotate("text", x = 4.45, y = 0.45, label = "Native TSD/TIR", size = 3.1) +
    ggplot2::annotate("text", x = 5.35, y = 1.0, hjust = 0, label = "Bar groups:", fontface = "bold", size = 3.1) +
    ggplot2::annotate("rect", xmin = 6.35, xmax = 6.67, ymin = 0.84, ymax = 1.16, fill = support_group_colors[["NATIVE"]], color = NA) +
    ggplot2::annotate("text", x = 6.75, y = 1.0, hjust = 0, label = "NATIVE", size = 3.1) +
    ggplot2::annotate("rect", xmin = 7.5, xmax = 7.82, ymin = 0.84, ymax = 1.16, fill = support_group_colors[["OVERLAP"]], color = NA) +
    ggplot2::annotate("text", x = 7.9, y = 1.0, hjust = 0, label = "OVERLAP", size = 3.1) +
    ggplot2::annotate("rect", xmin = 8.85, xmax = 9.17, ymin = 0.84, ymax = 1.16, fill = support_group_colors[["ANN"]], color = NA) +
    ggplot2::annotate("text", x = 9.25, y = 1.0, hjust = 0, label = "ANN", size = 3.1) +
    ggplot2::annotate("rect", xmin = 9.75, xmax = 10.07, ymin = 0.84, ymax = 1.16, fill = support_group_colors[["SEQ"]], color = NA) +
    ggplot2::annotate("text", x = 10.15, y = 1.0, hjust = 0, label = "SEQ", size = 3.1) +
    ggplot2::coord_cartesian(xlim = c(0, 11.2), ylim = c(0.1, 1.5), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  ggplot2::ggplotGrob(p)
}
