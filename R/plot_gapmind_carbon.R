.dnmb_gapmind_carbon_cazy_nodes <- function() {
  nodes <- .dnmb_cct_metabolism_nodes()
  # Remap types: backbone/ppp/entry_intermediate -> precursor (old convention)
  old_type <- ifelse(nodes$type %in% c("backbone", "ppp", "entry_intermediate"),
                     "precursor", nodes$type)
  data.frame(id = nodes$id, x = nodes$x, y = nodes$y,
             label = nodes$label, type = old_type,
             stringsAsFactors = FALSE)
}

# Step-level detection status for carbon pathways from genbank_table + aa.sum.steps
.dnmb_gapmind_carbon_step_status <- function(genbank_table, output_dir = NULL) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  path_col <- "GapMindCarbon_pathway_id"; step_col <- "GapMindCarbon_step_id"
  conf_col <- "GapMindCarbon_confidence"; bp_col <- "GapMindCarbon_on_best_path"
  lt_col <- "locus_tag"
  found <- data.frame(
    pathway_id = character(),
    step_id = character(),
    confidence = character(),
    locus_tag = character(),
    stringsAsFactors = FALSE
  )
  if (path_col %in% base::names(tbl)) {
    tbl <- tbl[!is.na(tbl[[path_col]]) & base::nzchar(tbl[[path_col]]), , drop = FALSE]
    if (base::nrow(tbl)) {
      has_bp <- bp_col %in% names(tbl)
      bp_tbl <- if (has_bp) tbl[!is.na(tbl[[bp_col]]) & tbl[[bp_col]] == TRUE, , drop = FALSE] else tbl
      has_lt <- lt_col %in% names(bp_tbl)
      found <- bp_tbl |>
        dplyr::group_by(pathway_id = .data[[path_col]], step_id = .data[[step_col]]) |>
        dplyr::summarise(
          confidence = dplyr::case_when(
            any(.data[[conf_col]] == "high") ~ "high",
            any(.data[[conf_col]] == "medium") ~ "medium", TRUE ~ "low"),
          locus_tag = if (has_lt) dplyr::first(.data[[lt_col]]) else NA_character_,
          .groups = "drop") |>
        as.data.frame(stringsAsFactors = FALSE)
    }
  }
  # Supplement from aa.sum.steps in carbon module directory
  steps_file <- NULL
  if (!is.null(output_dir)) {
    for (f in c(file.path(output_dir, "dnmb_module_gapmind_carbon", "aa.sum.steps"),
                file.path(output_dir, "dnmb_module_gapmindcarbon", "aa.sum.steps"),
                file.path(output_dir, "aa.sum.steps"))) {
      if (file.exists(f)) { steps_file <- f; break }
    }
  }
  if (!is.null(steps_file)) {
    st <- utils::read.delim(steps_file, stringsAsFactors = FALSE)
    if (all(c("onBestPath","pathway","step","score") %in% names(st))) {
      bp_st <- st[st$onBestPath == 1 & !duplicated(st[,c("pathway","step")]), , drop = FALSE]
      bp_st$confidence <- dplyr::case_when(bp_st$score>=2~"high", bp_st$score>=1~"medium", TRUE~"low")
      lt_s <- if ("locusId" %in% names(bp_st)) bp_st$locusId else NA_character_
      extra <- data.frame(pathway_id=bp_st$pathway, step_id=bp_st$step,
                          confidence=bp_st$confidence, locus_tag=lt_s, stringsAsFactors=FALSE)
      ekey <- paste0(extra$pathway_id,"::",extra$step_id)
      fkey <- paste0(found$pathway_id,"::",found$step_id)
      found <- rbind(found, extra[!ekey %in% fkey, , drop=FALSE])
    }
  }
  if (!base::nrow(found)) return(NULL)
  found
}

# Per-pathway completion stats for carbon sources
.dnmb_gapmind_carbon_pathway_stats <- function(step_status, output_dir = NULL) {
  if (is.null(step_status) || nrow(step_status) == 0) return(NULL)
  # Read total steps from aa.sum.steps
  steps_file <- NULL
  if (!is.null(output_dir)) {
    for (f in c(file.path(output_dir, "dnmb_module_gapmind_carbon", "aa.sum.steps"),
                file.path(output_dir, "dnmb_module_gapmindcarbon", "aa.sum.steps"),
                file.path(output_dir, "aa.sum.steps"))) {
      if (file.exists(f)) { steps_file <- f; break }
    }
  }
  # Get per-pathway total best-path steps and found counts
  if (!is.null(steps_file)) {
    st <- utils::read.delim(steps_file, stringsAsFactors = FALSE)
    if (all(c("onBestPath","pathway","step","score") %in% names(st))) {
      bp_st <- st[st$onBestPath == 1 & !duplicated(st[,c("pathway","step")]), , drop = FALSE]
      pw_totals <- stats::aggregate(step ~ pathway, data = bp_st, FUN = length)
      names(pw_totals) <- c("pathway_id", "n_total")
      pw_found <- stats::aggregate(
        step ~ pathway,
        data = bp_st[bp_st$score >= 1, , drop = FALSE],
        FUN = length)
      names(pw_found) <- c("pathway_id", "n_found")
      pstats <- merge(pw_totals, pw_found, by = "pathway_id", all.x = TRUE)
      pstats$n_found[is.na(pstats$n_found)] <- 0L
      pstats$fraction <- pstats$n_found / pstats$n_total
      # Determine best confidence per pathway
      conf_rank <- c(none = 0, low = 1, medium = 2, high = 3)
      bp_st$confidence <- dplyr::case_when(bp_st$score>=2~"high", bp_st$score>=1~"medium", TRUE~"none")
      pw_conf <- stats::aggregate(confidence ~ pathway, data = bp_st, FUN = function(x) {
        names(which.max(conf_rank[x]))
      })
      names(pw_conf) <- c("pathway_id", "best_confidence")
      pstats <- merge(pstats, pw_conf, by = "pathway_id", all.x = TRUE)
      return(pstats)
    }
  }
  # Fallback: from step_status only
  pw_found <- stats::aggregate(step_id ~ pathway_id, data = step_status, FUN = length)
  names(pw_found) <- c("pathway_id", "n_found")
  pw_found$n_total <- pw_found$n_found  # unknown total, assume all found
  pw_found$fraction <- 1.0
  pw_found
}

# Extract CAZy family summary from genbank_table
.dnmb_gapmind_carbon_cazy_summary <- function(genbank_table) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))
  if (is.null(family_col)) return(NULL)
  tbl <- tbl[!is.na(tbl[[family_col]]) & base::nzchar(tbl[[family_col]]), , drop = FALSE]
  if (!base::nrow(tbl)) return(NULL)
  tbl$family <- as.character(tbl[[family_col]])
  # Extract class prefix (GH, GT, CE, AA, PL, CBM)
  tbl$class <- sub("^(GH|GT|CE|AA|PL|CBM).*", "\\1", tbl$family)
  counts <- stats::aggregate(family ~ class, data = tbl, FUN = function(x) length(x))
  names(counts) <- c("class", "count")
  # Unique families per class
  fam_list <- stats::aggregate(family ~ class, data = tbl, FUN = function(x) {
    paste(sort(unique(x)), collapse = ", ")
  })
  names(fam_list) <- c("class", "families")
  merge(counts, fam_list, by = "class")
}

# Connection lines between carbon source nodes and backbone entry points
# Now delegates to the comprehensive .dnmb_cct_entry_edges()
.dnmb_gapmind_carbon_connections <- function(nodes) {
  .dnmb_cct_entry_edges(nodes)
}

# Main carbon + CAZy pathway map function (SNFG metabolism frame version)
.dnmb_plot_gapmind_carbon_cazy_pathway_map <- function(genbank_table, output_dir,
                                                        file_stub = "cazy_carbon_metabolism_frame") {
  step_status <- .dnmb_gapmind_carbon_step_status(genbank_table, output_dir = output_dir)
  if (is.null(step_status)) return(NULL)

  # Use new SNFG-annotated node table
  styled_nodes <- .dnmb_cct_nodes_styled()
  node_x   <- stats::setNames(styled_nodes$x, styled_nodes$id)
  node_y   <- stats::setNames(styled_nodes$y, styled_nodes$id)
  pstats   <- .dnmb_gapmind_carbon_pathway_stats(step_status, output_dir = output_dir)
  cazy_sum <- .dnmb_gapmind_carbon_cazy_summary(genbank_table)

  # Backbone + PPP internal edges
  bb_edges <- .dnmb_cct_backbone_edges(styled_nodes)
  # Carbon source -> intermediate -> backbone edges
  conn <- .dnmb_cct_entry_edges(styled_nodes)

  # --- Node subsets ---
  backbone_nodes <- styled_nodes[styled_nodes$type %in% c("backbone", "ppp"), , drop = FALSE]
  entry_inter    <- styled_nodes[styled_nodes$type == "entry_intermediate", , drop = FALSE]
  carbon_src     <- styled_nodes[styled_nodes$type == "carbon_source", , drop = FALSE]
  cazy_bridges   <- styled_nodes[styled_nodes$type == "cazy_bridge", , drop = FALSE]

  # --- Pathway category colors (derived from sugar_type) ---
  sugar_to_cat <- c(
    glucose = "monosaccharide", galactose = "monosaccharide",
    mannose = "monosaccharide", fructose = "monosaccharide",
    xylose = "monosaccharide", arabinose = "monosaccharide",
    ribose = "monosaccharide", rhamnose = "monosaccharide",
    fucose = "monosaccharide", GlcNAc = "sugar_deriv",
    GlcA = "sugar_deriv", GalA = "sugar_deriv",
    glucosamine = "sugar_deriv", NAG = "sugar_deriv",
    glycerol = "sugar_alcohol", mannitol = "sugar_alcohol",
    sorbitol = "sugar_alcohol", myoinositol = "sugar_alcohol",
    xylitol = "sugar_alcohol", organic_acid = "organic_acid",
    amino_acid = "amino_acid", nucleoside = "nucleoside",
    aromatic = "aromatic", generic = "other",
    trehalose = "disaccharide", sucrose = "disaccharide",
    maltose = "disaccharide", cellobiose = "disaccharide",
    lactose = "disaccharide", phospho_sugar = "intermediate",
    intermediate = "intermediate"
  )
  cat_colors <- c(
    monosaccharide = "#E6550D", disaccharide = "#31A354",
    sugar_deriv = "#756BB1", sugar_alcohol = "#3182BD",
    organic_acid = "#E7298A", amino_acid = "#D6604D",
    nucleoside = "#8C6BB1", aromatic = "#636363",
    intermediate = "#999999", other = "#AAAAAA"
  )

  # --- Assign confidence + fill to carbon source nodes ---
  conf_colors <- c(high = "#2CA25F", medium = "#FEC44F", low = "#F03B20", none = "#CCCCCC")
  carbon_src$confidence <- "none"
  carbon_src$fill_color <- conf_colors["none"]
  carbon_src$ratio_text <- ""
  carbon_src$cat <- vapply(carbon_src$sugar_type, function(st) {
    ct <- sugar_to_cat[st]; if (is.na(ct)) "other" else ct
  }, character(1))
  carbon_src$cat_color <- cat_colors[carbon_src$cat]
  carbon_src$cat_color[is.na(carbon_src$cat_color)] <- "#AAAAAA"
  for (i in seq_len(nrow(carbon_src))) {
    pid <- carbon_src$id[i]
    if (!is.null(pstats) && pid %in% pstats$pathway_id) {
      ps <- pstats[pstats$pathway_id == pid, , drop = FALSE]
      carbon_src$ratio_text[i] <- paste0(ps$n_found[1], "/", ps$n_total[1])
      frac <- ps$fraction[1]
      if (frac >= 0.75) {
        carbon_src$fill_color[i] <- conf_colors["high"]
        carbon_src$confidence[i] <- "high"
      } else if (frac >= 0.5) {
        carbon_src$fill_color[i] <- conf_colors["medium"]
        carbon_src$confidence[i] <- "medium"
      } else if (frac > 0) {
        carbon_src$fill_color[i] <- conf_colors["low"]
        carbon_src$confidence[i] <- "low"
      }
    }
  }
  carbon_src$confidence <- factor(carbon_src$confidence,
                                   levels = c("high", "medium", "low", "none"))

  # --- Connection edge colors ---
  if (!is.null(conn) && nrow(conn) > 0) {
    conn$edge_color <- "#DDDDDD"
    conn$edge_alpha <- 0.3
    for (i in seq_len(nrow(conn))) {
      src_id <- conn$from[i]
      idx <- match(src_id, carbon_src$id)
      if (!is.na(idx)) {
        conf <- as.character(carbon_src$confidence[idx])
        conn$edge_color[i] <- conf_colors[conf]
        conn$edge_alpha[i] <- ifelse(conf == "none", 0.15, 0.45)
      }
    }
  }

  # --- Backbone edge data ---
  glycolysis_edges <- bb_edges[bb_edges$pathway == "glycolysis", , drop = FALSE]
  ppp_edges        <- bb_edges[bb_edges$pathway == "ppp", , drop = FALSE]
  ed_edges         <- bb_edges[bb_edges$pathway == "ed", , drop = FALSE]

  # --- CAZy class summary ---
  cazy_label_df <- data.frame(
    id = c("CAZy_GH", "CAZy_CE", "CAZy_GT", "CAZy_AA"),
    class = c("GH", "CE", "GT", "AA"), stringsAsFactors = FALSE)
  cazy_label_df$detail <- ""
  cazy_label_df$count <- 0L
  if (!is.null(cazy_sum)) {
    for (i in seq_len(nrow(cazy_label_df))) {
      cls <- cazy_label_df$class[i]
      idx <- match(cls, cazy_sum$class)
      if (!is.na(idx)) {
        cazy_label_df$count[i] <- cazy_sum$count[idx]
        cazy_label_df$detail[i] <- cazy_sum$families[idx]
      }
    }
  }
  cazy_label_df$x <- node_x[cazy_label_df$id]
  cazy_label_df$y <- node_y[cazy_label_df$id]
  cazy_label_df$full_label <- ifelse(
    cazy_label_df$count > 0,
    paste0(cazy_label_df$class, " (", cazy_label_df$count, ")\n", cazy_label_df$detail),
    paste0(cazy_label_df$class, "\n(none)"))

  # --- CAZy -> disaccharide + polysaccharide targets ---
  cazy_tgt <- c("sucrose", "maltose", "cellobiose", "trehalose", "lactose")
  cazy_tgt <- cazy_tgt[cazy_tgt %in% names(node_x)]
  cazy_links <- if (length(cazy_tgt) > 0) {
    data.frame(
      x = rep(node_x["CAZy_GH"], length(cazy_tgt)),
      y = rep(node_y["CAZy_GH"], length(cazy_tgt)),
      xend = node_x[cazy_tgt], yend = node_y[cazy_tgt],
      stringsAsFactors = FALSE)
  } else NULL

  # ==========================================================
  # BUILD THE PLOT
  # ==========================================================
  p <- ggplot2::ggplot()

  # L0: Cell boundary — rounded rectangle
  cell_rect <- data.frame(
    xmin = 3.5, xmax = 13.0, ymin = 6.5, ymax = 23.0, stringsAsFactors = FALSE)
  p <- p +
    ggplot2::geom_rect(
      data = cell_rect,
      ggplot2::aes(xmin = .data$xmin, ymin = .data$ymin,
                   xmax = .data$xmax, ymax = .data$ymax),
      fill = "#FAFCFA", color = "#CCCCCC", linewidth = 0.6, linetype = "solid") +
    ggplot2::annotate("text", x = 12.8, y = 23.2, label = "Cytoplasm",
                      hjust = 1, vjust = 0, size = 2.2, fontface = "italic",
                      color = "#AAAAAA")

  # L1: Glycolysis backbone (thick, dark)
  if (nrow(glycolysis_edges) > 0) {
    p <- p + ggplot2::geom_segment(
      data = glycolysis_edges,
      ggplot2::aes(x = .data$from_x, y = .data$from_y,
                   xend = .data$to_x, yend = .data$to_y),
      color = "#555555", linewidth = 2.2, lineend = "round",
      arrow = ggplot2::arrow(length = grid::unit(0.05, "inches"), type = "closed"))
  }

  # L1b: PPP branch (teal)
  if (nrow(ppp_edges) > 0) {
    p <- p + ggplot2::geom_segment(
      data = ppp_edges,
      ggplot2::aes(x = .data$from_x, y = .data$from_y,
                   xend = .data$to_x, yend = .data$to_y),
      color = "#2B8CBE", linewidth = 1.2, lineend = "round", alpha = 0.7,
      arrow = ggplot2::arrow(length = grid::unit(0.04, "inches"), type = "closed"))
  }

  # L1c: Entner-Doudoroff branch (orange)
  if (nrow(ed_edges) > 0) {
    p <- p + ggplot2::geom_segment(
      data = ed_edges,
      ggplot2::aes(x = .data$from_x, y = .data$from_y,
                   xend = .data$to_x, yend = .data$to_y),
      color = "#D95F02", linewidth = 1.0, lineend = "round", alpha = 0.7,
      arrow = ggplot2::arrow(length = grid::unit(0.04, "inches"), type = "closed"))
  }

  # L2: Entry connection edges
  if (!is.null(conn) && nrow(conn) > 0) {
    for (i in seq_len(nrow(conn))) {
      row <- conn[i, , drop = FALSE]
      p <- p + ggplot2::geom_segment(
        data = row,
        ggplot2::aes(x = .data$from_x, y = .data$from_y,
                     xend = .data$to_x, yend = .data$to_y),
        color = row$edge_color, linewidth = 0.35, alpha = row$edge_alpha,
        lineend = "round",
        arrow = ggplot2::arrow(length = grid::unit(0.04, "inches"), type = "closed"))
    }
  }

  # L3: CAZy bridge -> target links (dashed purple)
  if (!is.null(cazy_links) && nrow(cazy_links) > 0) {
    p <- p + ggplot2::geom_segment(
      data = cazy_links,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#756BB1", linewidth = 0.4, linetype = "dashed", alpha = 0.5)
  }

  # L4: Backbone + PPP + entry intermediate nodes — SNFG shapes
  for (nset in list(backbone_nodes, entry_inter)) {
    if (nrow(nset) > 0) {
      for (i in seq_len(nrow(nset))) {
        nd <- nset[i, , drop = FALSE]
        p <- p + ggplot2::geom_point(
          data = nd, ggplot2::aes(x = .data$x, y = .data$y),
          shape = nd$snfg_shape, size = 4.5,
          fill = nd$snfg_fill, color = nd$snfg_color, stroke = 0.6)
      }
    }
  }
  # Backbone labels — right side for backbone, adaptive for entry intermediates
  bb_only <- backbone_nodes[backbone_nodes$type == "backbone", , drop = FALSE]
  ppp_only <- backbone_nodes[backbone_nodes$type == "ppp", , drop = FALSE]
  # Backbone: label to the right
  p <- p + ggplot2::geom_text(
    data = bb_only,
    ggplot2::aes(x = .data$x + 0.5, y = .data$y, label = .data$label),
    size = 1.8, fontface = "bold", color = "#222222", hjust = 0)
  # PPP: label to the right
  if (nrow(ppp_only) > 0) {
    p <- p + ggplot2::geom_text(
      data = ppp_only,
      ggplot2::aes(x = .data$x + 0.5, y = .data$y, label = .data$label),
      size = 1.6, fontface = "bold", color = "#2B8CBE", hjust = 0)
  }
  # Entry intermediates: label to the left (they sit between sources and backbone)
  if (nrow(entry_inter) > 0) {
    p <- p + ggplot2::geom_text(
      data = entry_inter,
      ggplot2::aes(x = .data$x, y = .data$y + 0.3, label = .data$label),
      size = 1.4, fontface = "bold", color = "#555555", hjust = 0.5)
  }

  # L5: Carbon source nodes — SNFG shape + confidence overlay
  for (i in seq_len(nrow(carbon_src))) {
    nd <- carbon_src[i, , drop = FALSE]
    p <- p + ggplot2::geom_point(
      data = nd, ggplot2::aes(x = .data$x, y = .data$y),
      shape = nd$snfg_shape, size = 4.5,
      fill = "white", color = nd$snfg_color, stroke = 0.5)
  }
  # Pie fill for pathway fraction
  if (!is.null(pstats)) {
    for (i in seq_len(nrow(carbon_src))) {
      pid <- carbon_src$id[i]
      if (pid %in% pstats$pathway_id) {
        ps <- pstats[pstats$pathway_id == pid, , drop = FALSE]
        frac <- ps$fraction[1]
        if (frac > 0) {
          pie <- .dnmb_pie_polygon(carbon_src$x[i], carbon_src$y[i], 0.15, frac)
          pie$group <- pid
          p <- p + ggplot2::geom_polygon(
            data = pie,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
            fill = carbon_src$fill_color[i], color = NA, alpha = 0.7)
        }
      }
    }
  }
  # Confidence border ring
  p <- p + ggplot2::geom_point(
    data = carbon_src,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$confidence),
    shape = 21, size = 4.5, fill = NA, stroke = 0.7, show.legend = TRUE)

  # L6: Carbon source labels — positioned based on x relative to backbone
  # Left-side sources: label to the left; right-side: label to the right
  carbon_src$lab_x <- ifelse(carbon_src$x < 7.0,
                              carbon_src$x - 0.4,
                              carbon_src$x + 0.4)
  carbon_src$lab_hjust <- ifelse(carbon_src$x < 7.0, 1, 0)
  for (i in seq_len(nrow(carbon_src))) {
    cs <- carbon_src[i, , drop = FALSE]
    p <- p + ggplot2::geom_label(
      data = cs,
      ggplot2::aes(x = .data$lab_x, y = .data$y, label = .data$label),
      size = 1.5, fontface = "bold", color = "white",
      fill = cs$cat_color, hjust = cs$lab_hjust,
      label.padding = grid::unit(0.06, "lines"), label.size = 0)
  }
  # Ratio text below node
  has_ratio <- nzchar(carbon_src$ratio_text)
  if (any(has_ratio)) {
    p <- p + ggplot2::geom_text(
      data = carbon_src[has_ratio, , drop = FALSE],
      ggplot2::aes(x = .data$x, y = .data$y - 0.2, label = .data$ratio_text),
      size = 1.0, color = "#444444", fontface = "bold")
  }

  # L7: CAZy bridge nodes — compact horizontal row
  p <- p +
    ggplot2::geom_point(
      data = cazy_bridges,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 23, size = 6, fill = "#E8E0F0", color = "#756BB1", stroke = 0.6)
  # CAZy labels: class name + count above, families below
  cazy_label_df$short_label <- ifelse(
    cazy_label_df$count > 0,
    paste0(cazy_label_df$class, " (", cazy_label_df$count, ")"),
    cazy_label_df$class)
  # Truncate family detail to max 30 chars
  cazy_label_df$short_detail <- ifelse(
    nchar(cazy_label_df$detail) > 30,
    paste0(substr(cazy_label_df$detail, 1, 28), "..."),
    cazy_label_df$detail)
  p <- p +
    ggplot2::geom_text(
      data = cazy_label_df,
      ggplot2::aes(x = .data$x, y = .data$y + 0.35, label = .data$short_label),
      size = 1.3, color = "#756BB1", fontface = "bold") +
    ggplot2::geom_text(
      data = cazy_label_df[nzchar(cazy_label_df$short_detail), , drop = FALSE],
      ggplot2::aes(x = .data$x, y = .data$y - 0.35, label = .data$short_detail),
      size = 0.9, color = "#9E88C7", fontface = "italic")

  # L8: Pathway labels
  p <- p +
    ggplot2::annotate("text", x = 7.0, y = 23.2, label = "GLYCOLYSIS",
                      size = 2.8, fontface = "bold", color = "#555555") +
    ggplot2::annotate("text", x = 10.5, y = 21.2, label = "PPP",
                      size = 2.2, fontface = "bold", color = "#2B8CBE") +
    ggplot2::annotate("text", x = 2.5, y = 14.8, label = "ED pathway",
                      size = 1.8, fontface = "bold", color = "#D95F02")

  # L9: Scale + legends + theme
  p <- p +
    ggplot2::scale_color_manual(
      name = "Pathway\nCompleteness",
      values = c(high = "#2CA25F", medium = "#FEC44F", low = "#F03B20", none = "#CCCCCC"),
      labels = c(high = "High (>=75%)", medium = "Medium (50-74%)",
                 low = "Low (<50%)", none = "Not found"),
      drop = FALSE)

  # SNFG legend — two-column layout at bottom-right
  snfg_items <- data.frame(
    sugar_type = c("glucose", "galactose", "mannose", "fructose",
                   "GlcNAc", "xylose", "arabinose", "fucose",
                   "rhamnose", "GlcA", "GalA", "ribose"),
    label = c("Glucose (blue circle)", "Galactose (yellow circle)",
              "Mannose (green circle)", "Fructose (lt-green circle)",
              "GlcNAc (blue square)", "Xylose (orange star)",
              "Arabinose (green star)", "Fucose (red triangle)",
              "Rhamnose (green triangle)", "GlcA (blue diamond)",
              "GalA (yellow diamond)", "Ribose (purple circle)"),
    stringsAsFactors = FALSE)
  snfg_s <- .dnmb_snfg_lookup(snfg_items$sugar_type)
  n_snfg <- nrow(snfg_items)
  n_col1 <- ceiling(n_snfg / 2)
  snfg_legend_x1 <- -1.5
  snfg_legend_x2 <- 4.5
  snfg_legend_y <- 5.5
  snfg_items$x <- ifelse(seq_len(n_snfg) <= n_col1, snfg_legend_x1, snfg_legend_x2)
  snfg_items$y <- vapply(seq_len(n_snfg), function(i) {
    if (i <= n_col1) snfg_legend_y - (i - 1) * 0.45
    else snfg_legend_y - (i - n_col1 - 1) * 0.45
  }, numeric(1))
  snfg_items$shape <- snfg_s$snfg_shape
  snfg_items$fill  <- snfg_s$snfg_fill
  snfg_items$color <- snfg_s$snfg_color
  for (i in seq_len(nrow(snfg_items))) {
    si <- snfg_items[i, , drop = FALSE]
    p <- p + ggplot2::geom_point(
      data = si, ggplot2::aes(x = .data$x, y = .data$y),
      shape = si$shape, size = 2.5, fill = si$fill, color = si$color, stroke = 0.5,
      show.legend = FALSE)
  }
  p <- p +
    ggplot2::geom_text(
      data = snfg_items,
      ggplot2::aes(x = .data$x + 0.25, y = .data$y, label = .data$label),
      size = 1.4, hjust = 0, color = "#444444") +
    ggplot2::annotate("text", x = snfg_legend_x1, y = snfg_legend_y + 0.5,
                      hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
                      color = "#333333", label = "SNFG Sugar Symbol Notation")

  # Category legend (below SNFG, two-column)
  legend_cats <- c("Monosaccharide", "Disaccharide", "Sugar derivative",
                   "Sugar alcohol", "Organic acid", "Amino acid",
                   "Nucleoside-derived", "Aromatic")
  legend_cols <- c("#E6550D", "#31A354", "#756BB1", "#3182BD",
                   "#E7298A", "#D6604D", "#8C6BB1", "#636363")
  n_cat1 <- ceiling(length(legend_cats) / 2)
  cat_legend_y <- snfg_legend_y - ceiling(n_snfg / 2) * 0.45 - 1.0
  n_cats <- length(legend_cats)
  cat_df <- data.frame(
    x = ifelse(seq_len(n_cats) <= n_cat1, snfg_legend_x1, snfg_legend_x2),
    y = vapply(seq_len(n_cats), function(i) {
      if (i <= n_cat1) cat_legend_y - (i - 1) * 0.45
      else cat_legend_y - (i - n_cat1 - 1) * 0.45
    }, numeric(1)),
    label = legend_cats, color = legend_cols, stringsAsFactors = FALSE)
  p <- p +
    ggplot2::geom_point(
      data = cat_df, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 15, size = 2.5, color = cat_df$color, show.legend = FALSE) +
    ggplot2::geom_text(
      data = cat_df, ggplot2::aes(x = .data$x + 0.25, y = .data$y, label = .data$label),
      size = 1.5, hjust = 0, color = "#333333") +
    ggplot2::annotate("text", x = snfg_legend_x1, y = cat_legend_y + 0.5,
                      hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
                      color = "#333333", label = "Carbon Source Category")

  p <- p +
    ggplot2::coord_fixed(ratio = 1,
      xlim = c(-2.0, 14.5), ylim = c(cat_legend_y - 2.5, 25.0), clip = "off") +
    ggplot2::labs(
      title = "Carbon Utilization & CAZy Metabolism Frame",
      subtitle = paste0(
        "Glycolysis + PPP + ED backbone | SNFG sugar notation | ",
        "Fill = pathway completeness | Arrow = metabolic entry")) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0,
                                          margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0, color = "#666666",
                                             margin = ggplot2::margin(b = 5)),
      legend.position = c(0.85, 0.12),
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(size = 7, face = "bold"),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.width = grid::unit(12, "pt"),
      legend.key.height = grid::unit(5, "pt"),
      plot.margin = ggplot2::margin(5, 5, 5, 5))

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  .dnmb_module_plot_save(p, pdf_path, width = 16, height = 20)
  list(pdf = pdf_path)
}
