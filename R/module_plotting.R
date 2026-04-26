.dnmb_module_plot_dir <- function(output_dir = getwd()) {
  dir <- file.path(output_dir, "visualizations")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(dir, winslash = "/", mustWork = FALSE)
}

.dnmb_module_plot_save <- function(plot, pdf_path, width = 12, height = 6,
                                   min_dim = 4, max_dim = 45) {
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
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, bg = "white")
}

.dnmb_module_plot_save_multi <- function(plots, pdf_path, labels = "AUTO", ncol = 1, rel_heights = NULL, width = 12, height = 8) {
  composite <- cowplot::plot_grid(plotlist = plots, labels = labels, ncol = ncol, rel_heights = rel_heights)
  ggplot2::ggsave(pdf_path, composite, width = width, height = height, bg = "white")
  invisible(composite)
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

  contig_length_obj <- get0("contig_length", envir = .GlobalEnv, inherits = FALSE)
  if (is.data.frame(contig_length_obj) && nrow(contig_length_obj) >= nrow(contigs)) {
    length_col <- contig_length_obj[[1]]
    if (!is.null(length_col)) {
      contigs$length_bp <- suppressWarnings(as.numeric(length_col[seq_len(nrow(contigs))]))
      contigs$length_source <- "global_contig_length"
    }
  }

  if (!"length_bp" %in% names(contigs) || all(is.na(contigs$length_bp))) {
    gbff_records <- .dnmb_parse_gbff_records(.dnmb_find_gbff_for_plot(output_dir))
    if (nrow(gbff_records) >= nrow(contigs)) {
      contigs$length_bp <- gbff_records$length_bp[seq_len(nrow(contigs))]
      contigs$accession <- gbff_records$accession[seq_len(nrow(contigs))]
      contigs$gbff_definition <- gbff_records$definition[seq_len(nrow(contigs))]
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
  run_plot("CAZyTransport",    function() .dnmb_plot_cazy_carbon_transport_map(genbank_table, output_dir = output_dir))
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
