.dnmb_gene_arrow_category_prophage <- function(product) {
  product <- tolower(as.character(product))
  dplyr::case_when(
    grepl("capsid|head protein|portal protein|terminase|major capsid|minor capsid|scaffold|prohead", product) ~ "Head/packaging",
    grepl("tail|spike|tube|baseplate|fiber|sheath|neck", product) ~ "Tail",
    grepl("helicase|primase|polymerase|dna pol|dnab|replicase|ssb|single.strand.binding", product) ~ "DNA replication",
    grepl("anti.crispr|anti.restriction|methyltransferase|methylase|ardb|ocr\\b", product) ~ "Anti-defense",
    grepl("integrase|excisionase", product) ~ "Integration",
    grepl("recombinase|bet protein|erf|sak|single.strand.anneal", product) ~ "Recombination",
    grepl("holin|endolysin|lysis|lysin|spanin", product) ~ "Lysis",
    grepl("regulator|repressor|cro|anti-termin", product) ~ "Regulation",
    grepl("hypothetical protein", product) ~ "Hypothetical",
    TRUE ~ "Other"
  )
}

.dnmb_gene_arrow_palette_prophage <- function() {
  c(
    "Head/packaging" = "#5B8FF9",
    "Tail" = "#5AD8A6",
    "DNA replication" = "#FF9F40",
    "Anti-defense" = "#FF4D4F",
    "Integration" = "#F6BD16",
    "Recombination" = "#FF6B8A",
    "Lysis" = "#E8684A",
    "Regulation" = "#6DC8EC",
    "Hypothetical" = "#D9D9D9",
    "Other" = "#9270CA"
  )
}

.dnmb_overview_window_plot <- function(contig_lengths,
                                       windows,
                                       fill_col,
                                       title,
                                       subtitle = NULL,
                                       label_col = NULL,
                                       palette = NULL,
                                       track_half_height = NULL,
                                       track_linewidth = 3.2,
                                       window_half_height = NULL,
                                       window_linewidth = NULL,
                                       window_geom = c("segment", "rect"),
                                       window_alpha = 0.92,
                                       window_border_color = "grey35",
                                       window_border_linewidth = 0.25) {
  contig_lengths <- as.data.frame(contig_lengths, stringsAsFactors = FALSE)
  windows <- as.data.frame(windows, stringsAsFactors = FALSE)
  if (!nrow(contig_lengths) || !nrow(windows)) {
    return(NULL)
  }
  contig_lengths$track <- 1
  windows$track <- 1
  window_geom <- match.arg(window_geom)
  if (is.null(window_linewidth) || is.na(window_linewidth)) {
    window_linewidth <- track_linewidth * 0.62
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contig_lengths,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = .data$track, yend = .data$track),
      linewidth = track_linewidth,
      color = "grey90",
      lineend = "round"
    ) +
    ggplot2::facet_wrap(~contig, scales = "free_x", ncol = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Genome coordinate (bp)",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.82, 1.23),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    )
  if (identical(window_geom, "segment")) {
    p <- p + ggplot2::geom_segment(
      data = windows,
      ggplot2::aes(
        x = .data$start,
        xend = .data$end,
        y = .data$track,
        yend = .data$track,
        color = .data[[fill_col]]
      ),
      linewidth = window_linewidth,
      alpha = window_alpha,
      lineend = "butt"
    )
  } else {
    p <- p + ggplot2::geom_rect(
      data = windows,
      ggplot2::aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = .data$track - window_half_height,
        ymax = .data$track + window_half_height,
        fill = .data[[fill_col]]
      ),
      color = window_border_color,
      linewidth = window_border_linewidth,
      alpha = window_alpha
    )
  }
  if (!is.null(label_col) && label_col %in% names(windows)) {
    p <- p + ggplot2::geom_text(
      data = windows,
      ggplot2::aes(x = (.data$start + .data$end) / 2, y = 1.145, label = .data[[label_col]]),
      size = 2.6,
      show.legend = FALSE
    ) +
      ggplot2::coord_cartesian(clip = "off")
  }
  if (!is.null(palette)) {
    if (identical(window_geom, "segment")) {
      p <- p + ggplot2::scale_color_manual(values = palette)
    } else {
      p <- p + ggplot2::scale_fill_manual(values = palette)
    }
  }
  p
}

.dnmb_arc_storyboard_spec <- function() {
  list(
    start_angle = pi / 2 + 0.34,
    end_angle = pi / 2 - 0.34,
    center_x = 0,
    center_y = -14.1,
    r_mid = 15.2,
    body_thickness = 0.46,
    head_thickness = 0.98,
    label_r = 16.55,
    overview_y = 2.55,
    overview_half_height = 0.08,
    xlim = c(-6.6, 6.6),
    ylim = c(-0.45, 2.95)
  )
}

.dnmb_arc_xy <- function(angle, radius, center_x = 0, center_y = 0) {
  data.frame(
    x = center_x + radius * cos(angle),
    y = center_y + radius * sin(angle),
    stringsAsFactors = FALSE
  )
}

.dnmb_arc_map_angle <- function(position,
                                xmin,
                                xmax,
                                start_angle = .dnmb_arc_storyboard_spec()$start_angle,
                                end_angle = .dnmb_arc_storyboard_spec()$end_angle) {
  position <- suppressWarnings(as.numeric(position))
  xmin <- suppressWarnings(as.numeric(xmin)[1])
  xmax <- suppressWarnings(as.numeric(xmax)[1])
  if (is.na(xmin) || is.na(xmax) || xmax <= xmin) {
    return(rep((start_angle + end_angle) / 2, length(position)))
  }
  frac <- pmax(0, pmin(1, (position - xmin) / (xmax - xmin)))
  start_angle + frac * (end_angle - start_angle)
}

.dnmb_arc_band_polygon <- function(xmin,
                                   xmax,
                                   start_angle = .dnmb_arc_storyboard_spec()$start_angle,
                                   end_angle = .dnmb_arc_storyboard_spec()$end_angle,
                                   r_inner = .dnmb_arc_storyboard_spec()$r_mid - .dnmb_arc_storyboard_spec()$head_thickness / 2,
                                   r_outer = .dnmb_arc_storyboard_spec()$r_mid + .dnmb_arc_storyboard_spec()$head_thickness / 2,
                                   center_x = .dnmb_arc_storyboard_spec()$center_x,
                                   center_y = .dnmb_arc_storyboard_spec()$center_y,
                                   n = 240L) {
  angles_outer <- seq(.dnmb_arc_map_angle(xmin, xmin, xmax, start_angle, end_angle), .dnmb_arc_map_angle(xmax, xmin, xmax, start_angle, end_angle), length.out = n)
  angles_inner <- seq(.dnmb_arc_map_angle(xmax, xmin, xmax, start_angle, end_angle), .dnmb_arc_map_angle(xmin, xmin, xmax, start_angle, end_angle), length.out = n)
  outer_xy <- .dnmb_arc_xy(angles_outer, r_outer, center_x = center_x, center_y = center_y)
  inner_xy <- .dnmb_arc_xy(angles_inner, r_inner, center_x = center_x, center_y = center_y)
  data.frame(
    x = c(outer_xy$x, inner_xy$x),
    y = c(outer_xy$y, inner_xy$y),
    stringsAsFactors = FALSE
  )
}

.dnmb_arc_gene_polygon <- function(start,
                                   end,
                                   direction,
                                   xmin,
                                   xmax,
                                   start_angle = .dnmb_arc_storyboard_spec()$start_angle,
                                   end_angle = .dnmb_arc_storyboard_spec()$end_angle,
                                   r_mid = .dnmb_arc_storyboard_spec()$r_mid,
                                   body_thickness = .dnmb_arc_storyboard_spec()$body_thickness,
                                   head_thickness = .dnmb_arc_storyboard_spec()$head_thickness,
                                   center_x = .dnmb_arc_storyboard_spec()$center_x,
                                   center_y = .dnmb_arc_storyboard_spec()$center_y,
                                   n = 110L,
                                   head_frac = 0.42) {
  start <- suppressWarnings(as.numeric(start)[1])
  end <- suppressWarnings(as.numeric(end)[1])
  if (is.na(start) || is.na(end)) {
    return(data.frame())
  }
  if (end < start) {
    tmp <- start
    start <- end
    end <- tmp
  }
  gene_len <- max(1, end - start)
  zoom_span <- max(1, xmax - xmin)
  fixed_head_len <- max(140, zoom_span * 0.07)
  head_len <- min(fixed_head_len, gene_len * 0.82)
  is_forward <- as.character(direction)[1] %in% c("+", "plus", "1")
  r_body_outer <- r_mid + body_thickness / 2
  r_body_inner <- r_mid - body_thickness / 2
  r_head_outer <- r_mid + head_thickness / 2
  r_head_inner <- r_mid - head_thickness / 2
  if (is_forward) {
    body_end <- max(start, end - head_len)
    angle_start <- .dnmb_arc_map_angle(start, xmin, xmax, start_angle, end_angle)
    angle_body_end <- .dnmb_arc_map_angle(body_end, xmin, xmax, start_angle, end_angle)
    tip_angle <- .dnmb_arc_map_angle(end, xmin, xmax, start_angle, end_angle)
    if (body_end <= start) {
      body_outer_xy <- .dnmb_arc_xy(angle_start, r_head_outer, center_x = center_x, center_y = center_y)
      body_inner_xy <- .dnmb_arc_xy(angle_start, r_head_inner, center_x = center_x, center_y = center_y)
    } else {
      body_outer_angles <- seq(angle_start, angle_body_end, length.out = n)
      body_inner_angles <- seq(angle_body_end, angle_start, length.out = n)
      body_outer_xy <- .dnmb_arc_xy(body_outer_angles, r_body_outer, center_x = center_x, center_y = center_y)
      body_inner_xy <- .dnmb_arc_xy(body_inner_angles, r_body_inner, center_x = center_x, center_y = center_y)
    }
    head_outer_xy <- .dnmb_arc_xy(angle_body_end, r_head_outer, center_x = center_x, center_y = center_y)
    head_inner_xy <- .dnmb_arc_xy(angle_body_end, r_head_inner, center_x = center_x, center_y = center_y)
    tip_xy <- .dnmb_arc_xy(tip_angle, r_mid, center_x = center_x, center_y = center_y)
    coords <- rbind(
      body_outer_xy,
      head_outer_xy,
      tip_xy,
      head_inner_xy,
      body_inner_xy
    )
  } else {
    body_start <- min(end, start + head_len)
    angle_body_start <- .dnmb_arc_map_angle(body_start, xmin, xmax, start_angle, end_angle)
    angle_end <- .dnmb_arc_map_angle(end, xmin, xmax, start_angle, end_angle)
    tip_angle <- .dnmb_arc_map_angle(start, xmin, xmax, start_angle, end_angle)
    if (body_start >= end) {
      body_outer_xy <- .dnmb_arc_xy(angle_end, r_head_outer, center_x = center_x, center_y = center_y)
      body_inner_xy <- .dnmb_arc_xy(angle_end, r_head_inner, center_x = center_x, center_y = center_y)
    } else {
      body_outer_angles <- seq(angle_body_start, angle_end, length.out = n)
      body_inner_angles <- seq(angle_end, angle_body_start, length.out = n)
      body_outer_xy <- .dnmb_arc_xy(body_outer_angles, r_body_outer, center_x = center_x, center_y = center_y)
      body_inner_xy <- .dnmb_arc_xy(body_inner_angles, r_body_inner, center_x = center_x, center_y = center_y)
    }
    head_outer_xy <- .dnmb_arc_xy(angle_body_start, r_head_outer, center_x = center_x, center_y = center_y)
    head_inner_xy <- .dnmb_arc_xy(angle_body_start, r_head_inner, center_x = center_x, center_y = center_y)
    tip_xy <- .dnmb_arc_xy(tip_angle, r_mid, center_x = center_x, center_y = center_y)
    coords <- rbind(
      head_outer_xy,
      body_outer_xy,
      body_inner_xy,
      head_inner_xy,
      tip_xy
    )
  }
  coords
}

.dnmb_arc_gene_plot <- function(tbl,
                                fill_col,
                                label_col = NULL,
                                title = NULL,
                                subtitle = NULL,
                                palette = NULL) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    return(NULL)
  }
  xmin <- min(tbl$start, na.rm = TRUE)
  xmax <- max(tbl$end, na.rm = TRUE)
  bg <- .dnmb_arc_band_polygon(xmin, xmax)
  poly_list <- lapply(seq_len(nrow(tbl)), function(i) {
    poly <- .dnmb_arc_gene_polygon(tbl$start[[i]], tbl$end[[i]], tbl$direction[[i]], xmin, xmax)
    if (!nrow(poly)) {
      return(NULL)
    }
    poly$feature_id <- i
    poly$fill_value <- as.character(tbl[[fill_col]][[i]])
    poly
  })
  poly_tbl <- dplyr::bind_rows(Filter(Negate(is.null), poly_list))
  poly_edge_tbl <- poly_tbl
  label_tbl <- NULL
  if (!is.null(label_col) && label_col %in% names(tbl)) {
    mid_angle <- .dnmb_arc_map_angle((tbl$start + tbl$end) / 2, xmin, xmax)
    label_tbl <- data.frame(
      x = 1.03 * cos(mid_angle),
      y = 1.03 * sin(mid_angle),
      angle = mid_angle * 180 / pi - 90,
      label = as.character(tbl[[label_col]]),
      stringsAsFactors = FALSE
    )
    label_tbl <- label_tbl[!is.na(label_tbl$label) & nzchar(label_tbl$label), , drop = FALSE]
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = bg, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey98", color = "grey75", linewidth = 0.4) +
    ggplot2::geom_polygon(data = poly_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value), color = "grey25", linewidth = 0.25) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      legend.position = "none",
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  if (!is.null(label_tbl) && nrow(label_tbl)) {
    p <- p + ggplot2::geom_text(
      data = label_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle),
      size = 2.4,
      hjust = 0.5,
      vjust = 0.5
    )
  }
  p
}

.dnmb_linear_map_x <- function(position, xmin, xmax, x_min_plot = -1, x_max_plot = 1) {
  position <- suppressWarnings(as.numeric(position))
  xmin <- suppressWarnings(as.numeric(xmin)[1])
  xmax <- suppressWarnings(as.numeric(xmax)[1])
  if (is.na(xmin) || is.na(xmax) || xmax <= xmin) {
    return(rep((x_min_plot + x_max_plot) / 2, length(position)))
  }
  frac <- pmax(0, pmin(1, (position - xmin) / (xmax - xmin)))
  x_min_plot + frac * (x_max_plot - x_min_plot)
}

.dnmb_overview_window_coords <- function(zoom_start,
                                         zoom_end,
                                         xmin,
                                         xmax,
                                         min_width = 0.05) {
  x_left <- .dnmb_linear_map_x(zoom_start, xmin, xmax)
  x_right <- .dnmb_linear_map_x(zoom_end, xmin, xmax)
  if (is.na(x_left) || is.na(x_right)) {
    return(c(-min_width / 2, min_width / 2))
  }
  if (x_right < x_left) {
    tmp <- x_left
    x_left <- x_right
    x_right <- tmp
  }
  width <- x_right - x_left
  if (width < min_width) {
    mid <- (x_left + x_right) / 2
    x_left <- mid - min_width / 2
    x_right <- mid + min_width / 2
  }
  c(x_left, x_right)
}

.dnmb_linear_gene_polygon <- function(start,
                                      end,
                                      direction,
                                      xmin,
                                      xmax,
                                      x_min_plot = -1,
                                      x_max_plot = 1,
                                      y_center = 1.55,
                                      height = 0.12,
                                      arrow_frac = 0.25) {
  start <- suppressWarnings(as.numeric(start)[1])
  end <- suppressWarnings(as.numeric(end)[1])
  if (is.na(start) || is.na(end)) {
    return(data.frame())
  }
  if (end < start) {
    tmp <- start
    start <- end
    end <- tmp
  }
  x_start <- .dnmb_linear_map_x(start, xmin, xmax, x_min_plot = x_min_plot, x_max_plot = x_max_plot)
  x_end <- .dnmb_linear_map_x(end, xmin, xmax, x_min_plot = x_min_plot, x_max_plot = x_max_plot)
  if (x_end < x_start) {
    tmp <- x_start
    x_start <- x_end
    x_end <- tmp
  }
  width <- max(0.01, x_end - x_start)
  arrow_len <- min(width * arrow_frac, 0.08)
  y_top <- y_center + height / 2
  y_bottom <- y_center - height / 2
  is_forward <- as.character(direction)[1] %in% c("+", "plus", "1")
  if (is_forward) {
    body_end <- max(x_start, x_end - arrow_len)
    poly <- data.frame(
      x = c(x_start, body_end, x_end, body_end, x_start),
      y = c(y_bottom, y_bottom, y_center, y_top, y_top)
    )
  } else {
    body_start <- min(x_end, x_start + arrow_len)
    poly <- data.frame(
      x = c(body_start, x_end, x_end, body_start, x_start),
      y = c(y_bottom, y_bottom, y_top, y_top, y_center)
    )
  }
  poly
}

.dnmb_zoom_connector_polygon <- function(zoom_start,
                                         zoom_end,
                                         xmin,
                                         xmax,
                                         arc_r_outer,
                                         arc_start_angle,
                                         arc_end_angle,
                                         arc_center_x = 0,
                                         arc_center_y = 0,
                                         y_overview,
                                         y_half_height = 0.14) {
  x_window <- .dnmb_overview_window_coords(zoom_start, zoom_end, xmin, xmax)
  x_left <- x_window[[1]]
  x_right <- x_window[[2]]
  ang_left <- .dnmb_arc_map_angle(zoom_start, xmin, xmax, start_angle = arc_start_angle, end_angle = arc_end_angle)
  ang_right <- .dnmb_arc_map_angle(zoom_end, xmin, xmax, start_angle = arc_start_angle, end_angle = arc_end_angle)
  left_xy <- .dnmb_arc_xy(ang_left, arc_r_outer, center_x = arc_center_x, center_y = arc_center_y)
  right_xy <- .dnmb_arc_xy(ang_right, arc_r_outer, center_x = arc_center_x, center_y = arc_center_y)
  data.frame(
    x = c(x_left, x_right, right_xy$x, left_xy$x),
    y = c(y_overview - y_half_height, y_overview - y_half_height, right_xy$y, left_xy$y)
  )
}

.dnmb_storyboard_plot <- function(tbl,
                                  contig_length,
                                  fill_col,
                                  label_col = NULL,
                                  title = NULL,
                                  subtitle = NULL,
                                  palette = NULL) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    return(NULL)
  }
  spec <- .dnmb_arc_storyboard_spec()
  xmin <- 0
  xmax <- suppressWarnings(as.numeric(contig_length)[1])
  if (is.na(xmax) || xmax <= 0) {
    xmax <- max(tbl$end, na.rm = TRUE)
  }
  zoom_start_raw <- min(tbl$start, na.rm = TRUE)
  zoom_end_raw <- max(tbl$end, na.rm = TRUE)
  zoom_pad <- max((zoom_end_raw - zoom_start_raw) * 0.18, 250)
  zoom_start <- max(xmin, zoom_start_raw - zoom_pad)
  zoom_end <- min(xmax, zoom_end_raw + zoom_pad)
  bg <- .dnmb_arc_band_polygon(
    zoom_start,
    zoom_end,
    spec$start_angle,
    spec$end_angle,
    r_inner = spec$r_mid - spec$head_thickness / 2,
    r_outer = spec$r_mid + spec$head_thickness / 2,
    center_x = spec$center_x,
    center_y = spec$center_y
  )
  connector <- .dnmb_zoom_connector_polygon(
    zoom_start,
    zoom_end,
    xmin,
    xmax,
    spec$r_mid + spec$head_thickness / 2,
    spec$start_angle,
    spec$end_angle,
    arc_center_x = spec$center_x,
    arc_center_y = spec$center_y,
    y_overview = spec$overview_y,
    y_half_height = spec$overview_half_height
  )
  x_window <- .dnmb_overview_window_coords(zoom_start, zoom_end, xmin, xmax)
  poly_list <- lapply(seq_len(nrow(tbl)), function(i) {
    poly <- .dnmb_arc_gene_polygon(
      tbl$start[[i]],
      tbl$end[[i]],
      tbl$direction[[i]],
      zoom_start,
      zoom_end,
      spec$start_angle,
      spec$end_angle,
      r_mid = spec$r_mid,
      body_thickness = spec$body_thickness,
      head_thickness = spec$head_thickness,
      center_x = spec$center_x,
      center_y = spec$center_y
    )
    if (!nrow(poly)) {
      return(NULL)
    }
    poly$feature_id <- i
    poly$fill_value <- as.character(tbl[[fill_col]][[i]])
    poly
  })
  arc_tbl <- dplyr::bind_rows(Filter(Negate(is.null), poly_list))
  overview_gene_polys <- lapply(seq_len(nrow(tbl)), function(i) {
    poly <- .dnmb_linear_gene_polygon(tbl$start[[i]], tbl$end[[i]], tbl$direction[[i]], xmin, xmax)
    if (!nrow(poly)) {
      return(NULL)
    }
    poly$feature_id <- i
    poly$fill_value <- as.character(tbl[[fill_col]][[i]])
    poly
  })
  overview_tbl <- dplyr::bind_rows(Filter(Negate(is.null), overview_gene_polys))
  label_tbl <- NULL
  if (!is.null(label_col) && label_col %in% names(tbl)) {
    mid_angle <- .dnmb_arc_map_angle((tbl$start + tbl$end) / 2, zoom_start, zoom_end, spec$start_angle, spec$end_angle)
    label_xy <- .dnmb_arc_xy(mid_angle, spec$label_r, center_x = spec$center_x, center_y = spec$center_y)
    label_tbl <- data.frame(
      x = label_xy$x,
      y = label_xy$y,
      angle = mid_angle * 180 / pi - 90,
      label = as.character(tbl[[label_col]]),
      stringsAsFactors = FALSE
    )
    label_tbl <- label_tbl[!is.na(label_tbl$label) & nzchar(label_tbl$label), , drop = FALSE]
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = connector, ggplot2::aes(x = .data$x, y = .data$y), fill = "#CBD5E1", alpha = 0.38, color = NA) +
    ggplot2::geom_segment(ggplot2::aes(x = -1, xend = 1, y = spec$overview_y, yend = spec$overview_y), linewidth = 4.6, color = "grey90", lineend = "round") +
    ggplot2::geom_polygon(data = overview_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value), color = "grey25", linewidth = 0.2) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x_window[[1]], xmax = x_window[[2]], ymin = spec$overview_y - spec$overview_half_height, ymax = spec$overview_y + spec$overview_half_height), fill = NA, color = "#475569", linewidth = 0.55, linetype = "dashed") +
    ggplot2::geom_polygon(data = bg, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey98", color = "grey75", linewidth = 0.4) +
    ggplot2::geom_polygon(data = arc_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value), color = "grey25", linewidth = 0.28, linejoin = "mitre") +
    ggplot2::coord_equal(xlim = spec$xlim, ylim = spec$ylim, clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 6, 4, 6)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  if (!is.null(label_tbl) && nrow(label_tbl)) {
    p <- p + ggplot2::geom_text(
      data = label_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle),
      size = 2.4
    )
  }
  p
}

.dnmb_gene_arrow_plot <- function(tbl,
                                  y_col,
                                  fill_col,
                                  title,
                                  subtitle = NULL,
                                  label_col = NULL,
                                  facet_col = NULL,
                                  palette = NULL) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  tbl$forward <- tbl$direction %in% c("+", "plus", "1", 1)
  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      xmin = .data$start,
      xmax = .data$end,
      y = .data[[y_col]],
      fill = .data[[fill_col]],
      forward = .data$forward
    )
  ) +
    gggenes::geom_gene_arrow(
      arrowhead_height = grid::unit(6.4, "mm"),
      arrowhead_width = grid::unit(3.4, "mm"),
      color = "grey25"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Genome coordinate (bp)",
      y = NULL,
      fill = fill_col
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
  if (!is.null(label_col) && label_col %in% names(tbl)) {
    p <- p + ggplot2::geom_text(
      data = tbl[!is.na(tbl[[label_col]]) & nzchar(tbl[[label_col]]), , drop = FALSE],
      ggplot2::aes(x = (.data$start + .data$end) / 2, y = .data[[y_col]], label = .data[[label_col]]),
      size = 2.5,
      vjust = -0.8,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
      ggplot2::coord_cartesian(clip = "off")
  }
  if (!is.null(facet_col) && facet_col %in% names(tbl)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_x", ncol = 1)
  }
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  p
}

.dnmb_compact_window_track <- function(contig_length,
                                       contig,
                                       start,
                                       end,
                                       fill = "#94A3B8",
                                       title = NULL,
                                       subtitle = NULL) {
  contig_length <- suppressWarnings(as.numeric(contig_length)[1])
  if (is.na(contig_length) || contig_length <= 0) {
    contig_length <- max(suppressWarnings(as.numeric(end)[1]), 1)
  }
  contig_tbl <- data.frame(
    contig = as.character(contig),
    length_bp = contig_length,
    track = 1,
    stringsAsFactors = FALSE
  )
  win_tbl <- data.frame(
    contig = as.character(contig),
    start = suppressWarnings(as.numeric(start)),
    end = suppressWarnings(as.numeric(end)),
    category = "window",
    track = 1,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contig_tbl,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = .data$track, yend = .data$track),
      linewidth = 3.2,
      color = "grey90",
      lineend = "round"
    ) +
    ggplot2::geom_rect(
      data = win_tbl,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = .data$track - 0.12, ymax = .data$track + 0.12),
      fill = fill,
      color = "grey25",
      linewidth = 0.25,
      alpha = 0.96
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      plot.margin = ggplot2::margin(3, 4, 1, 4)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

.dnmb_prophage_storyboard_plot <- function(sub_tbl,
                                           contig_length,
                                           region_id,
                                           palette,
                                           show_legend = TRUE,
                                           region_class = "complete") {
  sub_tbl <- as.data.frame(sub_tbl, stringsAsFactors = FALSE)
  if (!nrow(sub_tbl)) {
    return(NULL)
  }
  region_start_col <- .dnmb_pick_column(sub_tbl, c("Prophage_prophage_start", "prophage_start"))
  region_end_col <- .dnmb_pick_column(sub_tbl, c("Prophage_prophage_end", "prophage_end"))
  zoom_start_raw <- if (!is.null(region_start_col) && any(!is.na(sub_tbl[[region_start_col]]))) {
    suppressWarnings(min(as.numeric(sub_tbl[[region_start_col]]), na.rm = TRUE))
  } else {
    min(sub_tbl$start, na.rm = TRUE)
  }
  zoom_end_raw <- if (!is.null(region_end_col) && any(!is.na(sub_tbl[[region_end_col]]))) {
    suppressWarnings(max(as.numeric(sub_tbl[[region_end_col]]), na.rm = TRUE))
  } else {
    max(sub_tbl$end, na.rm = TRUE)
  }
  zoom_start <- max(0, zoom_start_raw)
  zoom_end <- min(contig_length, zoom_end_raw)

  plot_x_left <- -5.72
  plot_x_right <- 5.72
  overview_left <- -4.9
  overview_right <- 4.9
  line_y <- 2.24
  line_half <- 0.045
  highlight_half <- 0.018
  lower_left <- plot_x_left
  lower_right <- plot_x_right
  gene_center_y <- 1.80
  gene_height <- 0.08 * 0.60
  gene_top_y <- gene_center_y + gene_height / 2
  gene_bottom_y <- gene_center_y - gene_height / 2
  gene_boundary_height <- gene_height * 0.10
  lower_top <- gene_top_y + 0.028
  lower_bottom <- gene_bottom_y - 0.012
  gene_left <- lower_left
  gene_right <- lower_right
  funnel_left <- gene_left + 0.22
  funnel_right <- gene_right - 0.22
  compact_prophage_label <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("^phage\\s+", "", x, ignore.case = TRUE)
    x <- gsub("\\bfamily protein\\b", "protein", x, ignore.case = TRUE)
    x <- gsub("\\bhelix-turn-helix\\b", "HTH", x, ignore.case = TRUE)
    x <- gsub("\\s+", " ", x)
    vapply(
      x,
      function(one_label) paste(base::strwrap(one_label, width = 34), collapse = "\n"),
      character(1)
    )
  }

  highlight_x <- .dnmb_linear_map_x(
    c(zoom_start_raw, zoom_end_raw),
    xmin = 0,
    xmax = contig_length,
    x_min_plot = overview_left,
    x_max_plot = overview_right
  )
  highlight_pad <- min(0.05, max(0, (highlight_x[[2]] - highlight_x[[1]]) * 0.12))
  highlight_x_draw <- c(highlight_x[[1]] + highlight_pad, highlight_x[[2]] - highlight_pad)
  if (!is.finite(highlight_x_draw[[1]]) || !is.finite(highlight_x_draw[[2]]) || highlight_x_draw[[2]] <= highlight_x_draw[[1]]) {
    highlight_x_draw <- highlight_x
  }
  funnel_tbl <- data.frame(
    x = c(highlight_x_draw[[1]], highlight_x_draw[[2]], funnel_right, funnel_left),
    y = c(line_y - highlight_half, line_y - highlight_half, lower_top, lower_top),
    stringsAsFactors = FALSE
  )
  funnel_slices <- lapply(seq_len(14L), function(i) {
    t0 <- (i - 1) / 14
    t1 <- i / 14
    x_left0 <- highlight_x_draw[[1]] + (funnel_left - highlight_x_draw[[1]]) * t0
    x_right0 <- highlight_x_draw[[2]] + (funnel_right - highlight_x_draw[[2]]) * t0
    x_left1 <- highlight_x_draw[[1]] + (funnel_left - highlight_x_draw[[1]]) * t1
    x_right1 <- highlight_x_draw[[2]] + (funnel_right - highlight_x_draw[[2]]) * t1
    y0 <- (line_y - highlight_half) + (lower_top - (line_y - highlight_half)) * t0
    y1 <- (line_y - highlight_half) + (lower_top - (line_y - highlight_half)) * t1
    data.frame(
      x = c(x_left0, x_right0, x_right1, x_left1),
      y = c(y0, y0, y1, y1),
      alpha = 0.32 - (0.32 - 0.04) * ((t0 + t1) / 2),
      stringsAsFactors = FALSE
    )
  })
  lower_box_tbl <- data.frame(
    x = c(lower_left, lower_right, lower_right, lower_left),
    y = c(lower_bottom, lower_bottom, lower_top, lower_top),
    stringsAsFactors = FALSE
  )
  # Region class styling: partial = dashed border + pale, complete = solid
  is_partial <- identical(region_class, "partial")
  region_border_lty <- if (is_partial) "dashed" else "solid"
  region_border_col <- if (is_partial) "#D1D5DB" else "#94A3B8"
  region_highlight_fill <- if (is_partial) "#E5E7EB" else "#94A3B8"
  region_gene_alpha <- if (is_partial) 0.5 else 1.0

  poly_list <- lapply(seq_len(nrow(sub_tbl)), function(i) {
    poly <- .dnmb_linear_gene_polygon(
      start = sub_tbl$start[[i]],
      end = sub_tbl$end[[i]],
      direction = sub_tbl$direction[[i]],
      xmin = zoom_start,
      xmax = zoom_end,
      x_min_plot = gene_left,
      x_max_plot = gene_right,
      y_center = gene_center_y,
      height = gene_height,
      arrow_frac = 0.22
    )
    if (!nrow(poly)) {
      return(NULL)
    }
    poly$feature_id <- i
    poly$fill_value <- as.character(sub_tbl$category[[i]])
    poly
  })
  poly_tbl <- dplyr::bind_rows(Filter(Negate(is.null), poly_list))

  in_gene_label_tbl <- sub_tbl
  in_gene_label_tbl$xmin_plot <- .dnmb_linear_map_x(
    pmin(as.numeric(in_gene_label_tbl$start), as.numeric(in_gene_label_tbl$end)),
    xmin = zoom_start,
    xmax = zoom_end,
    x_min_plot = gene_left,
    x_max_plot = gene_right
  )
  in_gene_label_tbl$xmax_plot <- .dnmb_linear_map_x(
    pmax(as.numeric(in_gene_label_tbl$start), as.numeric(in_gene_label_tbl$end)),
    xmin = zoom_start,
    xmax = zoom_end,
    x_min_plot = gene_left,
    x_max_plot = gene_right
  )
  in_gene_label_tbl$x <- (in_gene_label_tbl$xmin_plot + in_gene_label_tbl$xmax_plot) / 2
  in_gene_label_tbl$plot_width <- abs(in_gene_label_tbl$xmax_plot - in_gene_label_tbl$xmin_plot)
  old_locus_label <- if ("old_locus_tag" %in% names(in_gene_label_tbl)) trimws(as.character(in_gene_label_tbl$old_locus_tag)) else rep("", nrow(in_gene_label_tbl))
  locus_label <- if ("locus_tag" %in% names(in_gene_label_tbl)) trimws(as.character(in_gene_label_tbl$locus_tag)) else rep("", nrow(in_gene_label_tbl))
  base_locus_label <- ifelse(nzchar(locus_label), locus_label, old_locus_label)
  wrapped_locus_label <- sub("_", "_\n", base_locus_label, fixed = TRUE)
  fallback_label <- if ("gene" %in% names(in_gene_label_tbl)) trimws(as.character(in_gene_label_tbl$gene)) else rep("", nrow(in_gene_label_tbl))
  in_gene_label_tbl$label <- ifelse(nzchar(base_locus_label), wrapped_locus_label, fallback_label)
  max_line_chars <- vapply(
    strsplit(in_gene_label_tbl$label, "\n", fixed = TRUE),
    function(parts) {
      parts <- parts[nzchar(parts)]
      if (!length(parts)) return(0L)
      max(nchar(parts), na.rm = TRUE)
    },
    integer(1)
  )
  line_count <- pmax(1L, lengths(strsplit(in_gene_label_tbl$label, "\n", fixed = TRUE)))
  label_threshold <- pmax(0.50, max_line_chars * 0.082 + (line_count - 1) * 0.09)
  in_gene_label_tbl$fill_hex <- unname(palette[as.character(in_gene_label_tbl$category)])
  in_gene_label_tbl$fill_hex[is.na(in_gene_label_tbl$fill_hex)] <- "#9270CA"
  in_gene_label_tbl$label_color <- vapply(in_gene_label_tbl$fill_hex, function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    luminance <- 0.2126 * rgb[1] + 0.7152 * rgb[2] + 0.0722 * rgb[3]
    if (is.finite(luminance) && luminance > 0.62) "#1F2937" else "white"
  }, character(1))
  in_gene_label_tbl <- in_gene_label_tbl[
    !is.na(in_gene_label_tbl$label) &
      nzchar(in_gene_label_tbl$label) &
      is.finite(in_gene_label_tbl$plot_width) &
      in_gene_label_tbl$plot_width >= label_threshold,
    ,
    drop = FALSE
  ]
  in_gene_label_tbl$y <- if (nrow(in_gene_label_tbl)) gene_center_y else numeric(0)

  core_categories <- c("Tail", "Head/packaging", "Integration", "Lysis", "Regulation", "DNA replication", "Anti-defense", "Recombination")
  label_tbl <- sub_tbl[sub_tbl$category %in% core_categories, , drop = FALSE]
  if (nrow(label_tbl)) {
    label_palette <- palette
    label_palette[names(label_palette) == "Hypothetical"] <- "#6B7280"
    label_tbl$label_color <- unname(label_palette[as.character(label_tbl$category)])
    label_tbl$label_color[is.na(label_tbl$label_color)] <- "#374151"
    label_tbl$annotation_label <- compact_prophage_label(label_tbl$product)
    label_tbl$display_label <- paste0(label_tbl$locus_tag, "\n", label_tbl$annotation_label)
    label_tbl$x <- .dnmb_linear_map_x(
      (label_tbl$start + label_tbl$end) / 2,
      xmin = zoom_start,
      xmax = zoom_end,
      x_min_plot = gene_left,
      x_max_plot = gene_right
    )
    label_tbl$y <- gene_top_y + 0.05
    label_tbl$priority <- dplyr::case_when(
      grepl("terminase|portal|capsid", tolower(label_tbl$product)) ~ 100,
      grepl("integrase|holin|lysin|endolysin", tolower(label_tbl$product)) ~ 90,
      grepl("tail|head|spike", tolower(label_tbl$product)) ~ 70,
      TRUE ~ 30
    ) + scales::rescale(abs(as.numeric(label_tbl$end) - as.numeric(label_tbl$start)), to = c(0, 20))
    label_tbl <- label_tbl[order(-label_tbl$priority), , drop = FALSE]
    label_tbl <- head(label_tbl, 8L)
    label_tbl$y <- gene_top_y + 0.05
  }

  cluster_tick_tbl <- data.frame(
    x = c(gene_left, gene_right),
    y = c(gene_bottom_y - gene_boundary_height - 0.004, gene_bottom_y - gene_boundary_height - 0.004),
    label = c(.dnmb_fmt_bp_exact(zoom_start_raw), .dnmb_fmt_bp_exact(zoom_end_raw)),
    hjust = c(0, 1),
    stringsAsFactors = FALSE
  )
  lane_label <- .dnmb_prophage_compact_genome_label(sub_tbl$contig[[1]])
  if (is.na(lane_label) || !nzchar(lane_label)) {
    lane_label <- as.character(sub_tbl$contig[[1]])
  }
  panel_ymin <- max(0.98, min(cluster_tick_tbl$y) - 0.035)

  plot_obj <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      ggplot2::aes(x = overview_left, xend = overview_right, y = line_y, yend = line_y),
      linewidth = 4.4,
      color = "grey90",
      lineend = "round"
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = highlight_x_draw[[1]],
        xmax = highlight_x_draw[[2]],
        ymin = line_y - highlight_half,
        ymax = line_y + highlight_half
      ),
      fill = region_highlight_fill,
      color = NA,
      alpha = 0.96
    ) +
    ggplot2::geom_polygon(
      data = lower_box_tbl,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = grDevices::adjustcolor("white", alpha.f = 0.9),
      color = NA
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = gene_left,
        xmax = gene_right,
        ymin = gene_bottom_y - gene_boundary_height,
        ymax = gene_bottom_y
      ),
      fill = "grey45",
      color = NA
    ) +
    ggplot2::geom_polygon(
      data = poly_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value),
      color = NA,
      alpha = region_gene_alpha
    ) +
    ggplot2::labs(
      title = paste0("Prophage ", region_id,
                     if (is_partial) " (partial/questionable)" else "",
                     " (", .dnmb_fmt_bp_label(zoom_end_raw - zoom_start_raw + 1), ")"),
      x = NULL,
      y = NULL,
      fill = "Category"
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 11.5, margin = ggplot2::margin(l = 18, b = 1.5)),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 4, r = 8, b = 12, l = 12),
      legend.position = if (show_legend) "bottom" else "none",
      legend.margin = ggplot2::margin(t = 0, b = 0),
      legend.box.margin = ggplot2::margin(t = 0, b = 0),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.height = grid::unit(10, "pt"),
      legend.key.width = grid::unit(16, "pt"),
      legend.title = ggplot2::element_text(margin = ggplot2::margin(r = 8))
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(plot_x_left, plot_x_right), ylim = c(panel_ymin, 2.38), clip = "off") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::guides(fill = ggplot2::guide_legend(title.position = "left", nrow = 1, byrow = TRUE)) +
    ggplot2::geom_text(
      data = cluster_tick_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, hjust = .data$hjust),
      size = 3.0,
      vjust = 1,
      color = "grey30",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = overview_left,
      y = line_y + 0.060,
      label = lane_label,
      hjust = 0,
      vjust = 0.5,
      size = 2.9,
      color = "#6B7280"
    )
  plot_obj <- plot_obj +
    ggplot2::geom_polygon(
      data = poly_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id),
      fill = NA,
      color = if (is_partial) "grey60" else "grey25",
      linewidth = 0.14,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = in_gene_label_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, color = .data$label_color),
      size = 1.55,
      lineheight = 0.84,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = cluster_tick_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, hjust = .data$hjust),
      size = 3.0,
      vjust = 1,
      color = "grey30",
      show.legend = FALSE,
      inherit.aes = FALSE
    )
  for (slice_tbl in funnel_slices) {
    plot_obj <- plot_obj +
      ggplot2::geom_polygon(
        data = slice_tbl,
        ggplot2::aes(x = .data$x, y = .data$y),
        fill = grDevices::adjustcolor("#CBD5E1", alpha.f = slice_tbl$alpha[[1]]),
        color = NA,
        inherit.aes = FALSE
      )
  }
  if (nrow(label_tbl)) {
    plot_obj <- plot_obj +
      ggrepel::geom_text_repel(
        data = label_tbl,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$display_label, color = .data$label_color),
        size = 2.3,
        direction = "y",
        nudge_y = 0.12,
        segment.size = 0.2,
        segment.color = "grey60",
        box.padding = 0.12,
        point.padding = 0.08,
        min.segment.length = 0,
        max.overlaps = 12,
        fontface = "bold",
        lineheight = 0.92,
        show.legend = FALSE,
        inherit.aes = FALSE
      )
  }
  attr(plot_obj, "panel_ymin") <- panel_ymin
  attr(plot_obj, "panel_ymax") <- 2.38
  attr(plot_obj, "track_top") <- gene_top_y
  attr(plot_obj, "track_bottom") <- gene_bottom_y - gene_boundary_height
  attr(plot_obj, "track_center") <- gene_center_y
  plot_obj
}

.dnmb_circular_window_plot <- function(contig_length,
                                       window_start,
                                       window_end,
                                       title = NULL,
                                       subtitle = NULL,
                                       fill = "#F87171") {
  contig_length <- suppressWarnings(as.numeric(contig_length)[1])
  if (is.na(contig_length) || contig_length <= 0) {
    return(NULL)
  }
  window_start <- max(0, suppressWarnings(as.numeric(window_start)[1]))
  window_end <- min(contig_length, suppressWarnings(as.numeric(window_end)[1]))
  if (is.na(window_start) || is.na(window_end) || window_end <= window_start) {
    return(NULL)
  }
  start_angle <- pi / 2
  end_angle <- pi / 2 - 2 * pi
  center_x <- 0
  center_y <- 1.5
  r_mid <- 1.06
  body_half <- 0.08
  bg <- .dnmb_arc_band_polygon(
    0,
    contig_length,
    start_angle = start_angle,
    end_angle = end_angle,
    r_inner = r_mid - body_half,
    r_outer = r_mid + body_half,
    center_x = center_x,
    center_y = center_y,
    n = 360L
  )
  win <- .dnmb_arc_band_polygon(
    window_start,
    window_end,
    start_angle = start_angle,
    end_angle = end_angle,
    r_inner = r_mid - body_half,
    r_outer = r_mid + body_half,
    center_x = center_x,
    center_y = center_y,
    n = 120L
  )
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = bg, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey95", color = "grey75", linewidth = 0.35) +
    ggplot2::geom_polygon(data = win, ggplot2::aes(x = .data$x, y = .data$y), fill = fill, color = "grey25", linewidth = 0.3, alpha = 0.92) +
    ggplot2::annotate("text", x = 0, y = -0.02, label = if (is.null(title)) "" else title, size = 5.4, fontface = "bold") +
    ggplot2::annotate("text", x = 0, y = -0.28, label = if (is.null(subtitle)) "" else subtitle, size = 3.6) +
    ggplot2::coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.45, 1.35), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
}

.dnmb_arc_zoom_plot <- function(tbl,
                                fill_col,
                                label_col = NULL,
                                title = NULL,
                                subtitle = NULL,
                                palette = NULL) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    return(NULL)
  }
  spec <- .dnmb_arc_storyboard_spec()
  zoom_start <- min(tbl$start, na.rm = TRUE)
  zoom_end <- max(tbl$end, na.rm = TRUE)
  bg <- .dnmb_arc_band_polygon(
    zoom_start,
    zoom_end,
    spec$start_angle,
    spec$end_angle,
    r_inner = spec$r_mid - spec$head_thickness / 2,
    r_outer = spec$r_mid + spec$head_thickness / 2,
    center_x = spec$center_x,
    center_y = spec$center_y
  )
  poly_list <- lapply(seq_len(nrow(tbl)), function(i) {
    poly <- .dnmb_arc_gene_polygon(
      tbl$start[[i]],
      tbl$end[[i]],
      tbl$direction[[i]],
      zoom_start,
      zoom_end,
      spec$start_angle,
      spec$end_angle,
      r_mid = spec$r_mid,
      body_thickness = spec$body_thickness,
      head_thickness = spec$head_thickness,
      center_x = spec$center_x,
      center_y = spec$center_y
    )
    if (!nrow(poly)) {
      return(NULL)
    }
    poly$feature_id <- i
    poly$fill_value <- as.character(tbl[[fill_col]][[i]])
    poly
  })
  poly_tbl <- dplyr::bind_rows(Filter(Negate(is.null), poly_list))
  label_tbl <- NULL
  if (!is.null(label_col) && label_col %in% names(tbl)) {
    labeled <- tbl[!is.na(tbl[[label_col]]) & nzchar(tbl[[label_col]]), , drop = FALSE]
    if (nrow(labeled)) {
      mid_angle <- .dnmb_arc_map_angle((labeled$start + labeled$end) / 2, zoom_start, zoom_end, spec$start_angle, spec$end_angle)
      label_xy <- .dnmb_arc_xy(mid_angle, spec$label_r, center_x = spec$center_x, center_y = spec$center_y)
      label_tbl <- data.frame(
        x = label_xy$x,
        y = label_xy$y,
        angle = mid_angle * 180 / pi - 90,
        label = as.character(labeled[[label_col]]),
        stringsAsFactors = FALSE
      )
    }
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = bg, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey98", color = "grey75", linewidth = 0.38) +
    ggplot2::geom_polygon(data = poly_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value), color = "grey25", linewidth = 0.28, linejoin = "mitre") +
    ggplot2::coord_equal(xlim = spec$xlim, ylim = c(-0.2, 2.6), clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      plot.margin = ggplot2::margin(0, 4, 0, 4)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  if (!is.null(label_tbl) && nrow(label_tbl)) {
    p <- p + ggplot2::geom_text(
      data = label_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle),
      size = 2.7,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.92
    )
  }
  p
}

.dnmb_arc_tick_segments <- function(positions,
                                    xmin,
                                    xmax,
                                    start_angle,
                                    end_angle,
                                    r_inner,
                                    r_outer,
                                    center_x,
                                    center_y) {
  positions <- suppressWarnings(as.numeric(positions))
  positions <- positions[!is.na(positions)]
  if (!length(positions)) {
    return(data.frame())
  }
  ang <- .dnmb_arc_map_angle(positions, xmin, xmax, start_angle, end_angle)
  inner_xy <- .dnmb_arc_xy(ang, r_inner, center_x = center_x, center_y = center_y)
  outer_xy <- .dnmb_arc_xy(ang, r_outer, center_x = center_x, center_y = center_y)
  data.frame(
    x = inner_xy$x,
    y = inner_xy$y,
    xend = outer_xy$x,
    yend = outer_xy$y,
    stringsAsFactors = FALSE
  )
}

.dnmb_arc_path_segments <- function(start,
                                    end,
                                    xmin,
                                    xmax,
                                    start_angle,
                                    end_angle,
                                    radius,
                                    center_x,
                                    center_y,
                                    n = 120L) {
  ang <- seq(
    .dnmb_arc_map_angle(start, xmin, xmax, start_angle, end_angle),
    .dnmb_arc_map_angle(end, xmin, xmax, start_angle, end_angle),
    length.out = n
  )
  xy <- .dnmb_arc_xy(ang, radius, center_x = center_x, center_y = center_y)
  data.frame(
    x = xy$x,
    y = xy$y,
    stringsAsFactors = FALSE
  )
}

.dnmb_arc_angle_path <- function(angle_start,
                                 angle_end,
                                 radius,
                                 center_x = 0,
                                 center_y = 0,
                                 n = 120L) {
  ang <- seq(angle_start, angle_end, length.out = n)
  xy <- .dnmb_arc_xy(ang, radius, center_x = center_x, center_y = center_y)
  data.frame(
    x = xy$x,
    y = xy$y,
    stringsAsFactors = FALSE
  )
}

.dnmb_funnel_polygon <- function(neck_angle,
                                 neck_radius,
                                 neck_span,
                                 mouth_angle,
                                 mouth_radius,
                                 mouth_span,
                                 center_x = 0,
                                 center_y = 0,
                                 n = 90L,
                                 neck_hold = 0.58) {
  t <- seq(0, 1, length.out = n)
  smooth <- pmax(0, pmin(1, (t - neck_hold) / (1 - neck_hold)))
  smooth <- smooth * smooth * (3 - 2 * smooth)
  mid_angle <- neck_angle + smooth * (mouth_angle - neck_angle)
  radii <- neck_radius + t * (mouth_radius - neck_radius)
  half_span <- (neck_span / 2) + smooth * ((mouth_span / 2) - (neck_span / 2))
  left_ang <- mid_angle + half_span
  right_ang <- mid_angle - half_span
  left_xy <- .dnmb_arc_xy(left_ang, radii, center_x = center_x, center_y = center_y)
  right_xy <- .dnmb_arc_xy(rev(right_ang), rev(radii), center_x = center_x, center_y = center_y)
  data.frame(
    x = c(left_xy$x, right_xy$x),
    y = c(left_xy$y, right_xy$y),
    stringsAsFactors = FALSE
  )
}

.dnmb_funnel_outline <- function(neck_angle,
                                 neck_radius,
                                 neck_span,
                                 mouth_angle,
                                 mouth_radius,
                                 mouth_span,
                                 center_x = 0,
                                 center_y = 0,
                                 n = 140L,
                                 neck_hold = 0.80) {
  t <- seq(0, 1, length.out = n)
  smooth <- pmax(0, pmin(1, (t - neck_hold) / (1 - neck_hold)))
  smooth <- smooth * smooth * (3 - 2 * smooth)

  neck_center <- .dnmb_arc_xy(neck_angle, neck_radius, center_x = center_x, center_y = center_y)
  mouth_center <- .dnmb_arc_xy(mouth_angle, mouth_radius, center_x = center_x, center_y = center_y)
  neck_left <- .dnmb_arc_xy(neck_angle + neck_span / 2, neck_radius, center_x = center_x, center_y = center_y)
  neck_right <- .dnmb_arc_xy(neck_angle - neck_span / 2, neck_radius, center_x = center_x, center_y = center_y)
  mouth_left <- .dnmb_arc_xy(mouth_angle + mouth_span / 2, mouth_radius, center_x = center_x, center_y = center_y)
  mouth_right <- .dnmb_arc_xy(mouth_angle - mouth_span / 2, mouth_radius, center_x = center_x, center_y = center_y)

  dx <- mouth_center$x - neck_center$x
  dy <- mouth_center$y - neck_center$y
  norm <- sqrt(dx^2 + dy^2)
  if (is.na(norm) || norm <= 0) {
    norm <- 1
  }
  ux <- dx / norm
  uy <- dy / norm
  vx <- -uy
  vy <- ux

  stem_x <- neck_center$x + t * dx
  stem_y <- neck_center$y + t * dy
  neck_half_width <- sqrt((neck_left$x - neck_center$x)^2 + (neck_left$y - neck_center$y)^2)
  mouth_half_width <- sqrt((mouth_left$x - mouth_center$x)^2 + (mouth_left$y - mouth_center$y)^2)
  half_width <- neck_half_width + smooth * (mouth_half_width - neck_half_width)

  center_xy <- data.frame(x = stem_x, y = stem_y, stringsAsFactors = FALSE)
  left_xy <- data.frame(
    x = stem_x + vx * half_width,
    y = stem_y + vy * half_width,
    stringsAsFactors = FALSE
  )
  right_xy <- data.frame(
    x = stem_x - vx * half_width,
    y = stem_y - vy * half_width,
    stringsAsFactors = FALSE
  )
  rbind(
    data.frame(x = center_xy$x, y = center_xy$y, part = "center", point_id = seq_along(t), stringsAsFactors = FALSE),
    data.frame(x = left_xy$x, y = left_xy$y, part = "left", point_id = seq_along(t), stringsAsFactors = FALSE),
    data.frame(x = right_xy$x, y = right_xy$y, part = "right", point_id = seq_along(t), stringsAsFactors = FALSE)
  )
}

.dnmb_sector_layout <- function(contig_lengths,
                                start_angle = pi / 2,
                                gap_angle = 0.08) {
  contig_lengths <- as.data.frame(contig_lengths, stringsAsFactors = FALSE)
  contig_lengths <- contig_lengths[!is.na(contig_lengths$length_bp) & contig_lengths$length_bp > 0, , drop = FALSE]
  if (!nrow(contig_lengths)) {
    return(data.frame())
  }
  total_span <- 2 * pi - gap_angle * nrow(contig_lengths)
  lengths <- suppressWarnings(as.numeric(contig_lengths$length_bp))
  spans <- total_span * lengths / sum(lengths)
  theta_start <- numeric(length(spans))
  theta_end <- numeric(length(spans))
  cursor <- start_angle
  for (i in seq_along(spans)) {
    theta_start[[i]] <- cursor
    theta_end[[i]] <- cursor - spans[[i]]
    cursor <- theta_end[[i]] - gap_angle
  }
  data.frame(
    contig = as.character(contig_lengths$contig),
    length_bp = lengths,
    theta_start = theta_start,
    theta_end = theta_end,
    span = spans,
    stringsAsFactors = FALSE
  )
}

.dnmb_sector_position_angle <- function(contig,
                                        position,
                                        sector_layout) {
  idx <- match(as.character(contig), sector_layout$contig)
  pos <- suppressWarnings(as.numeric(position))
  out <- rep(NA_real_, length(pos))
  ok <- !is.na(idx) & !is.na(pos)
  if (!any(ok)) {
    return(out)
  }
  lengths <- sector_layout$length_bp[idx[ok]]
  frac <- pmax(0, pmin(1, pos[ok] / lengths))
  out[ok] <- sector_layout$theta_start[idx[ok]] + frac * (sector_layout$theta_end[idx[ok]] - sector_layout$theta_start[idx[ok]])
  out
}

.dnmb_label_angle_facing <- function(angle_deg) {
  angle_deg <- suppressWarnings(as.numeric(angle_deg))
  normalized <- ((angle_deg + 180) %% 360) - 180
  needs_flip <- normalized < -90 | normalized > 90
  display_angle <- ifelse(needs_flip, normalized + 180, normalized)
  display_angle <- ifelse(display_angle > 180, display_angle - 360, display_angle)
  display_angle <- ifelse(display_angle < -180, display_angle + 360, display_angle)
  hjust <- rep(0.5, length(display_angle))
  data.frame(
    angle = display_angle,
    hjust = hjust,
    stringsAsFactors = FALSE
  )
}

.dnmb_label_radial_facing <- function(angle_deg) {
  angle_deg <- suppressWarnings(as.numeric(angle_deg))
  normalized <- ((angle_deg + 180) %% 360) - 180
  needs_flip <- normalized < -90 | normalized > 90
  display_angle <- ifelse(needs_flip, normalized + 180, normalized)
  display_angle <- ifelse(display_angle > 180, display_angle - 360, display_angle)
  display_angle <- ifelse(display_angle < -180, display_angle + 360, display_angle)
  data.frame(
    angle = display_angle,
    hjust = ifelse(needs_flip, 1, 0),
    stringsAsFactors = FALSE
  )
}

.dnmb_repel_label_angles <- function(angle_rad, label_text, base_sep = 0.035) {
  angle_rad <- suppressWarnings(as.numeric(angle_rad))
  label_text <- as.character(label_text)
  if (!length(angle_rad)) {
    return(angle_rad)
  }
  ord <- order(angle_rad)
  sorted_angle <- angle_rad[ord]
  sorted_text <- label_text[ord]
  adjusted <- sorted_angle
  if (length(adjusted) >= 2L) {
    for (i in 2:length(adjusted)) {
      prev_width <- nchar(sorted_text[[i - 1]])
      curr_width <- nchar(sorted_text[[i]])
      min_sep <- max(base_sep, min(0.085, 0.020 + (prev_width + curr_width) * 0.0012))
      if ((adjusted[[i]] - adjusted[[i - 1]]) < min_sep) {
        adjusted[[i]] <- adjusted[[i - 1]] + min_sep
      }
    }
  }
  shift <- mean(sorted_angle, na.rm = TRUE) - mean(adjusted, na.rm = TRUE)
  adjusted <- adjusted + shift
  out <- angle_rad
  out[ord] <- adjusted
  out
}

.dnmb_repel_label_angles_bounded <- function(angle_rad,
                                             label_text,
                                             lower_bound,
                                             upper_bound,
                                             base_sep = 0.035) {
  angle_rad <- suppressWarnings(as.numeric(angle_rad))
  lower_bound <- suppressWarnings(as.numeric(lower_bound))
  upper_bound <- suppressWarnings(as.numeric(upper_bound))
  label_text <- as.character(label_text)
  if (!length(angle_rad)) {
    return(angle_rad)
  }
  ord <- order(angle_rad)
  a <- angle_rad[ord]
  lo <- lower_bound[ord]
  hi <- upper_bound[ord]
  txt <- label_text[ord]
  lo2 <- pmin(lo, hi)
  hi2 <- pmax(lo, hi)
  adjusted <- pmax(pmin(a, hi2), lo2)
  if (length(adjusted) >= 2L) {
    for (pass in seq_len(2L)) {
      for (i in 2:length(adjusted)) {
        prev_width <- nchar(txt[[i - 1]])
        curr_width <- nchar(txt[[i]])
        min_sep <- max(base_sep, min(0.09, 0.018 + (prev_width + curr_width) * 0.0012))
        if ((adjusted[[i]] - adjusted[[i - 1]]) < min_sep) {
          adjusted[[i]] <- adjusted[[i - 1]] + min_sep
        }
      }
      adjusted <- pmin(adjusted, hi2)
      for (i in seq(length(adjusted) - 1L, 1L, by = -1L)) {
        next_width <- nchar(txt[[i + 1]])
        curr_width <- nchar(txt[[i]])
        min_sep <- max(base_sep, min(0.09, 0.018 + (next_width + curr_width) * 0.0012))
        if ((adjusted[[i + 1]] - adjusted[[i]]) < min_sep) {
          adjusted[[i]] <- adjusted[[i + 1]] - min_sep
        }
      }
      adjusted <- pmax(adjusted, lo2)
    }
  }
  out <- angle_rad
  out[ord] <- adjusted
  out
}

.dnmb_assign_panel_layout <- function(angle_mid,
                                      span_rad,
                                      ring_radii = c(2.45, 3.15, 3.85),
                                      gap_rad = 0.06) {
  angle_mid <- suppressWarnings(as.numeric(angle_mid))
  angle_mid <- atan2(sin(angle_mid), cos(angle_mid))
  span_rad <- suppressWarnings(as.numeric(span_rad))
  ord <- order(angle_mid, decreasing = TRUE)
  assigned <- data.frame(
    display_angle = angle_mid,
    span_rad = span_rad,
    ring_id = rep(1L, length(angle_mid)),
    stringsAsFactors = FALSE
  )
  assigned <- assigned[ord, , drop = FALSE]
  ring_items <- vector("list", length(ring_radii))
  angle_dist <- function(a, b) {
    abs(atan2(sin(a - b), cos(a - b)))
  }
  for (i in seq_len(nrow(assigned))) {
    ang <- assigned$display_angle[[i]]
    span <- assigned$span_rad[[i]]
    picked <- NA_integer_
    for (ring in seq_along(ring_radii)) {
      if (!length(ring_items[[ring]])) {
        picked <- ring
        break
      }
      d <- vapply(ring_items[[ring]], function(one) {
        angle_dist(ang, one$mid) - ((span + one$span) / 2 + gap_rad)
      }, numeric(1))
      if (all(d >= 0)) {
        picked <- ring
        break
      }
    }
    if (is.na(picked)) {
      ring_load <- vapply(ring_items, length, integer(1))
      picked <- which.min(ring_load)
    }
    assigned$ring_id[[i]] <- picked
    ring_items[[picked]][[length(ring_items[[picked]]) + 1L]] <- list(mid = ang, span = span)
  }
  back <- assigned
  back$orig_index <- ord
  back <- back[order(back$orig_index), , drop = FALSE]
  back$r_mid <- ring_radii[back$ring_id]
  back
}

.dnmb_arc_context_zoom_plot <- function(tbl,
                                        context_start,
                                        context_end,
                                        highlight_start,
                                        highlight_end,
                                        fill_col,
                                        label_col = NULL,
                                        title = NULL,
                                        subtitle = NULL,
                                        palette = NULL,
                                        highlight_color = "#F87171") {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    return(NULL)
  }
  spec <- .dnmb_arc_storyboard_spec()
  zoom_start <- min(tbl$start, na.rm = TRUE)
  zoom_end <- max(tbl$end, na.rm = TRUE)

  context_r_mid <- spec$r_mid - 1.72
  context_half <- 0.13
  tick_inner <- context_r_mid - 0.05
  tick_outer <- context_r_mid + 0.05

  context_bg <- .dnmb_arc_band_polygon(
    context_start,
    context_end,
    spec$start_angle,
    spec$end_angle,
    r_inner = context_r_mid - context_half,
    r_outer = context_r_mid + context_half,
    center_x = spec$center_x,
    center_y = spec$center_y,
    n = 220L
  )
  context_hi <- .dnmb_arc_band_polygon(
    highlight_start,
    highlight_end,
    spec$start_angle,
    spec$end_angle,
    r_inner = context_r_mid - context_half,
    r_outer = context_r_mid + context_half,
    center_x = spec$center_x,
    center_y = spec$center_y,
    n = 120L
  )
  context_ticks <- .dnmb_arc_tick_segments(
    positions = tbl$midpoint,
    xmin = context_start,
    xmax = context_end,
    start_angle = spec$start_angle,
    end_angle = spec$end_angle,
    r_inner = tick_inner,
    r_outer = tick_outer,
    center_x = spec$center_x,
    center_y = spec$center_y
  )
  defense_ticks <- .dnmb_arc_tick_segments(
    positions = tbl$midpoint[tbl[[fill_col]] != "neighbor"],
    xmin = context_start,
    xmax = context_end,
    start_angle = spec$start_angle,
    end_angle = spec$end_angle,
    r_inner = tick_inner - 0.01,
    r_outer = tick_outer + 0.01,
    center_x = spec$center_x,
    center_y = spec$center_y
  )
  context_outline <- .dnmb_arc_path_segments(
    start = highlight_start,
    end = highlight_end,
    xmin = context_start,
    xmax = context_end,
    start_angle = spec$start_angle,
    end_angle = spec$end_angle,
    radius = context_r_mid + context_half + 0.015,
    center_x = spec$center_x,
    center_y = spec$center_y,
    n = 90L
  )

  poly_list <- lapply(seq_len(nrow(tbl)), function(i) {
    poly <- .dnmb_arc_gene_polygon(
      tbl$start[[i]],
      tbl$end[[i]],
      tbl$direction[[i]],
      zoom_start,
      zoom_end,
      spec$start_angle,
      spec$end_angle,
      r_mid = spec$r_mid,
      body_thickness = spec$body_thickness,
      head_thickness = spec$head_thickness,
      center_x = spec$center_x,
      center_y = spec$center_y
    )
    if (!nrow(poly)) {
      return(NULL)
    }
    poly$feature_id <- i
    poly$fill_value <- as.character(tbl[[fill_col]][[i]])
    poly
  })
  outer_tbl <- dplyr::bind_rows(Filter(Negate(is.null), poly_list))
  zoom_backbone <- .dnmb_arc_path_segments(
    start = zoom_start,
    end = zoom_end,
    xmin = zoom_start,
    xmax = zoom_end,
    start_angle = spec$start_angle,
    end_angle = spec$end_angle,
    radius = spec$r_mid,
    center_x = spec$center_x,
    center_y = spec$center_y,
    n = 140L
  )
  endpoint_angles <- .dnmb_arc_map_angle(
    c(zoom_start, zoom_end),
    zoom_start,
    zoom_end,
    spec$start_angle,
    spec$end_angle
  )
  endpoint_xy <- .dnmb_arc_xy(
    endpoint_angles,
    spec$r_mid + spec$head_thickness / 2 + 0.16,
    center_x = spec$center_x,
    center_y = spec$center_y
  )
  endpoint_tbl <- data.frame(
    x = endpoint_xy$x,
    y = endpoint_xy$y,
    label = c(.dnmb_fmt_bp_exact(zoom_start), .dnmb_fmt_bp_exact(zoom_end)),
    hjust = c(1, 0),
    stringsAsFactors = FALSE
  )

  label_tbl <- NULL
  if (!is.null(label_col) && label_col %in% names(tbl)) {
    labeled <- tbl[!is.na(tbl[[label_col]]) & nzchar(tbl[[label_col]]), , drop = FALSE]
    if (nrow(labeled)) {
      mid_angle <- .dnmb_arc_map_angle((labeled$start + labeled$end) / 2, zoom_start, zoom_end, spec$start_angle, spec$end_angle)
      label_xy <- .dnmb_arc_xy(mid_angle, spec$label_r, center_x = spec$center_x, center_y = spec$center_y)
      label_tbl <- data.frame(
        x = label_xy$x,
        y = label_xy$y,
        angle = mid_angle * 180 / pi - 90,
        label = as.character(labeled[[label_col]]),
        stringsAsFactors = FALSE
      )
    }
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = context_bg, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey97", color = NA) +
    ggplot2::geom_segment(data = context_ticks, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), linewidth = 0.18, color = "#CBD5E1") +
    ggplot2::geom_segment(data = defense_ticks, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), linewidth = 0.26, color = highlight_color) +
    ggplot2::geom_path(data = context_outline, ggplot2::aes(x = .data$x, y = .data$y), linewidth = 1.0, color = highlight_color, lineend = "round") +
    ggplot2::geom_path(data = zoom_backbone, ggplot2::aes(x = .data$x, y = .data$y), linewidth = 0.34, color = "grey45", lineend = "round") +
    ggplot2::geom_polygon(data = outer_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value), color = "grey25", linewidth = 0.28, linejoin = "mitre") +
    ggplot2::coord_equal(xlim = spec$xlim, ylim = c(-0.25, 2.7), clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      plot.margin = ggplot2::margin(0, 4, 0, 4)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  if (!is.null(label_tbl) && nrow(label_tbl)) {
    p <- p + ggplot2::geom_text(
      data = label_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle),
      size = 2.8
    )
  }
  p <- p + ggplot2::geom_text(
    data = endpoint_tbl,
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, hjust = .data$hjust),
    size = 2.5,
    vjust = 0.5,
    color = "grey25",
    inherit.aes = FALSE
  )
  p
}
