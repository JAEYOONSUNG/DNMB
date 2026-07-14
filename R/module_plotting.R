.dnmb_module_plot_dir <- function(output_dir = getwd()) {
  dir <- file.path(output_dir, "visualizations")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(dir, winslash = "/", mustWork = FALSE)
}

.dnmb_module_plot_save <- function(plot, pdf_path, width = 12, height = 6,
                                   min_dim = 4, max_dim = 45,
                                   device = NULL) {
  ar <- width / height
  width  <- max(min_dim, min(max_dim, width))
  height <- max(min_dim, min(max_dim, height))
  # Restore aspect ratio after clamping: shrink the larger side
  if (width / height != ar) {
    if (ar >= 1) {
      height <- min(max_dim, max(min_dim, width / ar))
    } else {
      width <- min(max_dim, max(min_dim, height * ar))
    }
  }
  save_args <- list(
    filename = pdf_path, plot = plot,
    width = width, height = height, bg = "white",
    device = .dnmb_plot_pdf_device
  )
  # `device` remains in the internal signature for compatibility. All module
  # PDFs intentionally use the package device so no caller can bypass the
  # embedded, editable Arial contract.
  do.call(ggplot2::ggsave, save_args)
}

.dnmb_module_plot_save_multi <- function(plots, pdf_path, labels = "AUTO", ncol = 1, rel_heights = NULL, width = 12, height = 8) {
  composite <- cowplot::plot_grid(plotlist = plots, labels = labels, ncol = ncol, rel_heights = rel_heights)
  ggplot2::ggsave(pdf_path, composite, width = width, height = height, bg = "white",
                  device = .dnmb_plot_pdf_device)
  invisible(composite)
}

.dnmb_overview_title <- function(label, title) {
  label <- trimws(as.character(label)[1])
  title <- trimws(as.character(title)[1])
  if (is.na(label) || !nzchar(label)) {
    return(title)
  }
  if (is.na(title) || !nzchar(title)) {
    return(label)
  }
  paste0(label, "   ", title)
}

.dnmb_overview_title_theme <- function(size = 11) {
  ggplot2::theme(
    text = ggplot2::element_text(family = .dnmb_plot_font_family()),
    plot.title = ggplot2::element_text(
      family = .dnmb_plot_font_family(),
      face = "bold",
      size = size,
      hjust = 0,
      vjust = 1,
      margin = ggplot2::margin(t = 0, r = 0, b = 2, l = 2)
    ),
    plot.title.position = "panel"
  )
}

.dnmb_overview_tag_plot <- function(plot, label, title = NULL, size = 11) {
  if (is.null(plot)) {
    return(NULL)
  }
  if (is.null(title)) {
    title <- plot$labels$title %||% ""
  }
  plot +
    ggplot2::labs(title = .dnmb_overview_title(label, title)) +
    .dnmb_overview_title_theme(size = size)
}

.dnmb_inline_discrete_legend_plot <- function(title, colors, max_per_row = 5L) {
  if (is.null(colors) || !length(colors)) {
    return(NULL)
  }
  labels <- names(colors)
  if (is.null(labels) || !length(labels)) {
    labels <- as.character(seq_along(colors))
  }
  n <- length(labels)
  n_rows <- ceiling(n / max_per_row)
  items_per_row <- ceiling(n / n_rows)
  title_width <- max(0.82, nchar(title) * 0.09)
  key_width <- 0.18
  key_gap <- 0.055
  item_gap <- 0.14
  label_widths <- pmax(0.56, nchar(labels) * 0.072)
  item_widths <- key_width + key_gap + label_widths
  slot_width <- max(item_widths) + item_gap
  outer_pad <- 0.55
  row_height <- 0.32
  row_gap <- 0.18
  total_height <- n_rows * row_height + (n_rows - 1) * row_gap + 0.36
  row_items <- split(seq_len(n), ceiling(seq_len(n) / items_per_row))
  row_max_width <- max(vapply(row_items, function(idx) {
    length(idx) * slot_width - item_gap
  }, numeric(1)))
  inner_width <- title_width + 0.24 + row_max_width
  total_width <- inner_width + 2 * outer_pad
  title_x <- outer_pad + title_width / 2
  key_left <- outer_pad + title_width + 0.24
  item_tbl_list <- list()
  for (r in seq_along(row_items)) {
    idx <- row_items[[r]]
    row_y_center <- total_height - 0.18 - (r - 1) * (row_height + row_gap)
    row_width <- length(idx) * slot_width - item_gap
    row_left <- key_left + (row_max_width - row_width) / 2
    for (j in seq_along(idx)) {
      i <- idx[[j]]
      left <- row_left + (j - 1) * slot_width
      item_tbl_list[[length(item_tbl_list) + 1L]] <- data.frame(
        label = labels[[i]],
        fill = unname(colors[[i]]),
        key_xmin = left,
        key_xmax = left + key_width,
        key_ymin = row_y_center - 0.12,
        key_ymax = row_y_center + 0.12,
        label_x = left + key_width + key_gap,
        label_y = row_y_center,
        stringsAsFactors = FALSE
      )
    }
  }
  item_tbl <- dplyr::bind_rows(item_tbl_list)
  title_y <- total_height - 0.18
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = title_x, y = title_y, label = title, hjust = 0.5, vjust = 0.5, size = 4.1) +
    ggplot2::geom_rect(
      data = item_tbl,
      ggplot2::aes(xmin = .data$key_xmin, xmax = .data$key_xmax, ymin = .data$key_ymin, ymax = .data$key_ymax, fill = .data$fill),
      color = "grey35",
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = item_tbl,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$label),
      hjust = 0,
      vjust = 0.5,
      size = 3.6,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, total_width), ylim = c(0, total_height), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
}

.dnmb_fmt_bp_label <- function(bp) {
  bp <- suppressWarnings(as.numeric(bp)[1])
  if (is.na(bp) || bp <= 0) {
    return(NA_character_)
  }
  if (bp >= 1e6) {
    return(paste0(format(round(bp / 1e6, 2), trim = TRUE, nsmall = 2), " Mb"))
  }
  paste0(format(round(bp / 1e3, 1), trim = TRUE, nsmall = 1), " kb")
}

.dnmb_fmt_bp_exact <- function(bp) {
  bp <- suppressWarnings(as.numeric(bp)[1])
  if (is.na(bp)) {
    return(NA_character_)
  }
  paste0(format(round(bp), big.mark = ",", scientific = FALSE, trim = TRUE), " bp")
}

.dnmb_pretty_contig_label <- function(definition = NA_character_, accession = NA_character_, fallback = NA_character_) {
  label <- definition
  if (is.na(label) || !nzchar(label)) {
    label <- fallback
  }
  label <- trimws(gsub("\\s+", " ", label))
  label <- sub(",\\s*complete genome\\.?$", "", label, ignore.case = TRUE)
  label <- sub(",\\s*complete sequence\\.?$", "", label, ignore.case = TRUE)
  label <- sub("\\s+chromosome,?.*$", "", label, ignore.case = TRUE)
  label <- sub("\\s+plasmid,?.*$", "", label, ignore.case = TRUE)
  label <- sub("\\s+contig,?.*$", "", label, ignore.case = TRUE)
  label <- sub("\\.$", "", label)
  if (!is.na(accession) && nzchar(accession)) {
    return(paste0(label, " | ", accession))
  }
  label
}

.dnmb_pick_column <- function(tbl, candidates) {
  found <- candidates[candidates %in% names(tbl)]
  if (!length(found)) {
    return(NULL)
  }
  found[[1]]
}

.dnmb_find_gbff_for_plot <- function(output_dir = getwd()) {
  candidates <- c(
    list.files(output_dir, pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE, ignore.case = TRUE),
    list.files(getwd(), pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE, ignore.case = TRUE)
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates)) {
    return(NULL)
  }
  candidates[[1]]
}

.dnmb_parse_gbff_records <- function(gbff_path) {
  if (is.null(gbff_path) || !file.exists(gbff_path)) {
    return(data.frame())
  }

  lines <- readLines(gbff_path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) {
    return(data.frame())
  }

  trim_field <- function(x) gsub("\\s+", " ", trimws(x))
  records <- list()
  current <- NULL
  i <- 1L
  while (i <= length(lines)) {
    line <- lines[[i]]
    if (startsWith(line, "LOCUS")) {
      if (!is.null(current)) {
        records[[length(records) + 1L]] <- current
      }
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      current <- list(
        order = length(records) + 1L,
        accession = if (length(parts) >= 2L) parts[[2]] else NA_character_,
        definition = NA_character_,
        length_bp = NA_real_
      )
      i <- i + 1L
      next
    }

    if (is.null(current)) {
      i <- i + 1L
      next
    }

    if (startsWith(line, "DEFINITION")) {
      definition <- sub("^DEFINITION\\s+", "", line)
      j <- i + 1L
      while (j <= length(lines) && grepl("^\\s{12}\\S", lines[[j]])) {
        definition <- paste(definition, trim_field(lines[[j]]))
        j <- j + 1L
      }
      current$definition <- trim_field(definition)
      i <- j
      next
    }

    if (grepl("^\\s+source\\s+\\d+\\.\\.\\d+", line)) {
      bounds <- sub("^\\s+source\\s+(\\d+)\\.\\.(\\d+).*$", "\\1;\\2", line)
      parts <- strsplit(bounds, ";", fixed = TRUE)[[1]]
      current$length_bp <- suppressWarnings(as.numeric(parts[[2]]))
      i <- i + 1L
      next
    }

    i <- i + 1L
  }

  if (!is.null(current)) {
    records[[length(records) + 1L]] <- current
  }
  if (!length(records)) {
    return(data.frame())
  }

  out <- do.call(rbind, lapply(records, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

.dnmb_contig_lengths_for_plot <- function(genbank_table, output_dir = getwd()) {
  tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
  tbl$start <- suppressWarnings(as.numeric(tbl$start))
  tbl$end <- suppressWarnings(as.numeric(tbl$end))

  contigs <- unique(tbl[, intersect(c("contig_number", "contig"), names(tbl)), drop = FALSE])
  if (!nrow(contigs) || !"contig" %in% names(contigs)) {
    return(data.frame())
  }
  if ("contig_number" %in% names(contigs)) {
    contigs$contig_number <- suppressWarnings(as.numeric(contigs$contig_number))
    contigs <- contigs[order(contigs$contig_number, contigs$contig), , drop = FALSE]
  } else {
    contigs <- contigs[order(contigs$contig), , drop = FALSE]
    contigs$contig_number <- seq_len(nrow(contigs))
  }
  # Deduplicate by contig name to avoid duplicate factor levels in circlize
  contigs <- contigs[!duplicated(contigs$contig), , drop = FALSE]
  rownames(contigs) <- NULL

  max_end <- tbl |>
    dplyr::group_by(.data$contig) |>
    dplyr::summarise(max_end = max(.data$end, na.rm = TRUE), .groups = "drop")
  contigs$max_end <- max_end$max_end[match(contigs$contig, max_end$contig)]

  # Map per-replicon lengths by replicon order (contig_number), never by row
  # position. Both the global `contig_length` table and the parsed gbff records
  # are ordered by replicon, so a positional take (seq_len(nrow(contigs)))
  # misassigns lengths whenever `contigs` is a subset that does not start at
  # replicon 1 and run contiguously — e.g. a module table filtered to systems
  # sitting on the chromosome plus a high-numbered plasmid. That misassignment
  # can hand a contig a length shorter than its gene coordinates, collapsing
  # downstream zoom windows to zero rows.
  contig_idx <- suppressWarnings(as.integer(contigs$contig_number))

  contig_length_obj <- get0("contig_length", envir = .GlobalEnv, inherits = FALSE)
  if (is.data.frame(contig_length_obj) && nrow(contig_length_obj) >= nrow(contigs)) {
    length_col <- contig_length_obj[[1]]
    if (!is.null(length_col)) {
      pick <- ifelse(!is.na(contig_idx) & contig_idx >= 1L & contig_idx <= length(length_col),
                     contig_idx, NA_integer_)
      contigs$length_bp <- suppressWarnings(as.numeric(length_col[pick]))
      contigs$length_source <- "global_contig_length"
    }
  }

  if (!"length_bp" %in% names(contigs) || all(is.na(contigs$length_bp))) {
    gbff_records <- .dnmb_parse_gbff_records(.dnmb_find_gbff_for_plot(output_dir))
    if (nrow(gbff_records) >= nrow(contigs)) {
      g_idx <- match(contig_idx, gbff_records$order)
      contigs$length_bp <- gbff_records$length_bp[g_idx]
      contigs$accession <- gbff_records$accession[g_idx]
      contigs$gbff_definition <- gbff_records$definition[g_idx]
      contigs$length_source <- "gbff_source"
    }
  }

  if (!"length_bp" %in% names(contigs)) {
    contigs$length_bp <- NA_real_
  }
  missing_length <- is.na(contigs$length_bp) | contigs$length_bp <= 0
  if (any(missing_length)) {
    contigs$length_bp[missing_length] <- contigs$max_end[missing_length]
    if (!"length_source" %in% names(contigs)) {
      contigs$length_source <- NA_character_
    }
    contigs$length_source[missing_length] <- "max_gene_end"
  }

  if (!"accession" %in% names(contigs)) {
    contigs$accession <- NA_character_
  }
  contigs$sector_label <- vapply(
    seq_len(nrow(contigs)),
    function(i) .dnmb_pretty_contig_label(
      definition = if ("gbff_definition" %in% names(contigs)) contigs$gbff_definition[[i]] else NA_character_,
      accession = contigs$accession[[i]],
      fallback = contigs$contig[[i]]
    ),
    character(1)
  )
  contigs
}

.dnmb_module_replicon_plot_data <- function(genbank_table, output_dir = getwd()) {
  genes <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
  required <- c("start", "end", "locus_tag")
  if (!nrow(genes) || !all(required %in% names(genes))) {
    return(NULL)
  }

  genes$start <- suppressWarnings(as.numeric(genes$start))
  genes$end <- suppressWarnings(as.numeric(genes$end))
  genes$midpoint <- (genes$start + genes$end) / 2
  if ("contig_number" %in% names(genes)) {
    contig_number <- suppressWarnings(as.integer(genes$contig_number))
  } else if ("contig" %in% names(genes)) {
    contig_number <- as.integer(factor(as.character(genes$contig), levels = unique(as.character(genes$contig))))
  } else {
    contig_number <- rep.int(1L, nrow(genes))
  }
  if (anyNA(contig_number)) {
    replacement <- as.integer(factor(
      if ("contig" %in% names(genes)) as.character(genes$contig) else seq_len(nrow(genes)),
      levels = unique(if ("contig" %in% names(genes)) as.character(genes$contig) else seq_len(nrow(genes)))
    ))
    contig_number[is.na(contig_number)] <- replacement[is.na(contig_number)]
  }
  genes$replicon_order <- contig_number
  genes$replicon_id <- sprintf("DNMB_CONTIG_%03d", contig_number)

  contig_cols <- intersect(c("replicon_id", "replicon_order", "contig"), names(genes))
  contigs <- genes[!duplicated(genes$replicon_id), contig_cols, drop = FALSE]
  contigs <- contigs[order(contigs$replicon_order, contigs$replicon_id), , drop = FALSE]
  max_end <- tapply(genes$end, genes$replicon_id, max, na.rm = TRUE)

  gbff_records <- .dnmb_parse_gbff_records(.dnmb_find_gbff_for_plot(output_dir))
  if (nrow(gbff_records)) {
    idx <- match(contigs$replicon_order, gbff_records$order)
    contigs$length_bp <- suppressWarnings(as.numeric(gbff_records$length_bp[idx]))
    contigs$accession <- as.character(gbff_records$accession[idx])
    contigs$definition <- as.character(gbff_records$definition[idx])
  } else {
    contigs$length_bp <- NA_real_
    contigs$accession <- NA_character_
    contigs$definition <- NA_character_
  }
  missing_length <- !is.finite(contigs$length_bp) | contigs$length_bp <= 0
  contigs$length_bp[missing_length] <- unname(max_end[contigs$replicon_id[missing_length]])

  is_plasmid <- grepl("\\bplasmid\\b", contigs$definition, ignore.case = TRUE)
  plasmid_name <- sub(
    ".*\\bplasmid\\s+([^,.;]+).*$", "\\1", contigs$definition,
    ignore.case = TRUE
  )
  contigs$replicon_short <- ifelse(
    is_plasmid,
    paste("Plasmid", plasmid_name),
    ifelse(contigs$replicon_order == 1L, "Chromosome", paste("Replicon", contigs$replicon_order))
  )
  accession_label <- ifelse(
    !is.na(contigs$accession) & nzchar(contigs$accession),
    paste0(" | ", contigs$accession),
    ""
  )
  contigs$sector_label <- paste0(contigs$replicon_short, accession_label)
  contigs$facet_label <- paste0(
    contigs$sector_label, " (", scales::label_comma()(contigs$length_bp), " bp)"
  )
  rownames(contigs) <- NULL
  list(genes = genes, contigs = contigs)
}

.dnmb_module_feature_label <- function(gene, locus_tag) {
  gene <- trimws(as.character(gene))
  locus_tag <- as.character(locus_tag)
  locus_suffix <- sub("^.*_", "", locus_tag)
  ifelse(
    !is.na(gene) & nzchar(gene),
    paste0(gene, " (", locus_suffix, ")"),
    locus_tag
  )
}

.dnmb_module_pack_replicon_labels <- function(tbl, contigs,
                                                label_col = "feature_label",
                                                priority_col = NULL,
                                                panel_width_in = 7.8,
                                                font_size_pt = 6.2) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    tbl$label_tier <- integer()
    tbl$label_x <- numeric()
    tbl$label_half_width <- numeric()
    return(tbl)
  }
  if (!all(c("replicon_id", "midpoint", label_col) %in% names(tbl))) {
    stop("Replicon label packing requires replicon_id, midpoint, and label text.", call. = FALSE)
  }
  priority <- if (!is.null(priority_col) && priority_col %in% names(tbl)) {
    suppressWarnings(as.numeric(tbl[[priority_col]]))
  } else {
    rep.int(0, nrow(tbl))
  }
  priority[!is.finite(priority)] <- 0
  tbl$label_tier <- 0L
  tbl$label_x <- suppressWarnings(as.numeric(tbl$midpoint))
  tbl$label_half_width <- 0

  for (replicon_id in unique(as.character(tbl$replicon_id))) {
    ii <- which(as.character(tbl$replicon_id) == replicon_id)
    contig_len <- contigs$length_bp[match(replicon_id, contigs$replicon_id)]
    if (!length(contig_len) || !is.finite(contig_len) || contig_len <= 0) {
      contig_len <- max(tbl$end[ii], na.rm = TRUE)
    }
    label_chars <- vapply(
      strsplit(as.character(tbl[[label_col]][ii]), "\n", fixed = TRUE),
      function(x) max(nchar(x), 1L),
      numeric(1)
    )
    width_fraction <- (
      label_chars * font_size_pt * 0.52 / 72 / panel_width_in
    ) + 0.012
    half_width <- pmin(0.24, width_fraction / 2) * contig_len
    label_x <- pmin(pmax(tbl$midpoint[ii], half_width), contig_len - half_width)
    left <- label_x - half_width
    right <- label_x + half_width
    ord <- order(label_x, -priority[ii])
    tier_ends <- numeric()
    assigned <- integer(length(ii))
    for (jj in ord) {
      available <- which(left[jj] >= tier_ends)
      tier <- if (length(available)) available[[1]] else length(tier_ends) + 1L
      if (tier > length(tier_ends)) {
        tier_ends <- c(tier_ends, -Inf)
      }
      assigned[jj] <- tier - 1L
      tier_ends[tier] <- right[jj] + 0.006 * contig_len
    }
    tbl$label_tier[ii] <- assigned
    tbl$label_x[ii] <- label_x
    tbl$label_half_width[ii] <- half_width
  }
  tbl
}

.dnmb_contig_ordered_table <- function(genbank_table, required = c("contig", "start", "end")) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  if (!base::all(required %in% base::names(tbl))) {
    return(data.frame())
  }
  tbl$start <- suppressWarnings(base::as.numeric(tbl$start))
  tbl$end <- suppressWarnings(base::as.numeric(tbl$end))
  tbl$midpoint <- (tbl$start + tbl$end) / 2
  contigs <- .dnmb_contig_lengths_for_plot(tbl)
  if (base::nrow(contigs)) {
    tbl$contig <- factor(tbl$contig, levels = contigs$contig)
  }
  tbl
}

# ---------- GapMind AA Pathway Map — iPath-style metabolic network ----------

# Metabolite nodes — systematic textbook layout
# Glycolysis backbone vertical at x=1.0, pathways branch rightward, step ~ 0.5 units
dnmb_render_module_plots <- function(genbank_table, output_dir = getwd(), cache_root = NULL) {
  plots <- list()
  dropped <- list()

  run_plot <- function(name, fn) {
    reason <- NULL
    res <- tryCatch(
      fn(),
      error = function(e) {
        reason <<- paste0("error: ", conditionMessage(e))
        NULL
      }
    )
    if (is.list(res)) {
      plots[[name]] <<- res
    } else {
      if (is.null(reason)) {
        reason <- "no data (plot function returned NULL)"
      }
      dropped[[name]] <<- reason
    }
  }

  run_plot("GapMindAA",        function() .dnmb_plot_gapmind_aa_pathway_map(genbank_table, output_dir = output_dir))
  run_plot("CAZyTransport", function() .dnmb_plot_cazy_carbon_transport_map(
    genbank_table,
    output_dir = output_dir,
    file_stub = "CAZy_overview"
  ))
  run_plot("DefenseFinder",    function() .dnmb_plot_defensefinder_module(genbank_table, output_dir = output_dir))
  run_plot("dbAPIS",           function() .dnmb_plot_dbapis_module(genbank_table, output_dir = output_dir))
  run_plot("AcrFinder",        function() .dnmb_plot_acrfinder_module(genbank_table, output_dir = output_dir))
  run_plot("mRNAcal",          function() .dnmb_plot_mrnacal_module(genbank_table, output_dir = output_dir))
  run_plot("PADLOC",           function() .dnmb_plot_padloc_module(genbank_table, output_dir = output_dir))
  run_plot("DefensePredictor", function() .dnmb_plot_defensepredictor_module(genbank_table, output_dir = output_dir))
  run_plot("ISelement",        function() .dnmb_plot_iselement_module(genbank_table, output_dir = output_dir))
  if ("Prophage_prophage_id" %in% names(genbank_table)) {
    run_plot("Prophage",       function() .dnmb_plot_prophage_module(genbank_table, output_dir = output_dir, cache_root = cache_root))
  }
  run_plot("dbCAN",            function() .dnmb_plot_dbcan_module(genbank_table, output_dir = output_dir))
  run_plot("MEROPS",           function() .dnmb_plot_merops_module(genbank_table, output_dir = output_dir))
  run_plot("PAZy",             function() .dnmb_plot_pazy_module(genbank_table, output_dir = output_dir, cache_root = cache_root))
  run_plot("REBASEfinder",     function() .dnmb_plot_rebasefinder_module(genbank_table, output_dir = output_dir, cache_root = cache_root))
  run_plot("DefenseOverview",      function() .dnmb_plot_integrated_defense_module(genbank_table, output_dir = output_dir, defensefinder_activity = "Defense"))
  run_plot("AntiDefenseOverview",  function() .dnmb_plot_integrated_defense_module(genbank_table, output_dir = output_dir, defensefinder_activity = "Anti-defense"))

  if (length(dropped)) {
    message("[DNMB plots] ", length(dropped), " module plot(s) skipped:")
    for (name in names(dropped)) {
      message("  - ", name, ": ", dropped[[name]])
    }
  }

  plots
}

# REBASEfinder overview — delegates to R/plot_rebasefinder.R.
# Do not fall back to DefenseViz's native RM_system_dotplot.pdf here; that is
# the old-style plot and should never overwrite DNMB's REBASE_overview.pdf.
.dnmb_plot_rebasefinder_module <- function(genbank_table, output_dir, cache_root = NULL) {
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  expected_pdf <- file.path(plot_dir, "REBASE_overview.pdf")
  had_previous <- file.exists(expected_pdf)
  previous_mtime <- if (had_previous) file.info(expected_pdf)$mtime else as.POSIXct(NA)
  result <- tryCatch(
    .dnmb_plot_rebasefinder_overview(genbank_table, output_dir, cache_root = cache_root),
    error = function(e) NULL
  )
  if (!is.null(result) && file.exists(expected_pdf) && file.info(expected_pdf)$size > 500) {
    return(list(pdf = expected_pdf))
  }
  if (file.exists(expected_pdf)) {
    current_mtime <- file.info(expected_pdf)$mtime
    newly_written <- !had_previous || (!is.na(previous_mtime) && !is.na(current_mtime) && current_mtime > previous_mtime)
    if (isTRUE(newly_written) && file.info(expected_pdf)$size > 500) {
      return(list(pdf = expected_pdf))
    }
    if (had_previous) {
      unlink(expected_pdf, force = TRUE)
    }
  }
  NULL
}
