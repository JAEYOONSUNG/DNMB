#' Shared machinery for comparative per-module heatmaps.
#'
#' All `dnmb_plot_comparative_*()` plots share the same visual grammar:
#' rows = genomes grouped by Genus/Species, columns = module subtypes,
#' dot size / gradient = occurrence count, right-side Diversity bar,
#' Count(n) legend centered under the bar, Genus/Species legend under
#' the rightmost strain labels. The only per-module difference is the
#' collector that produces the long `File`/`subtype` frame (and the
#' title / output filename / auto-run module alias).
#'
#' @keywords internal
#' @noRd

.dnmb_comparative_discover_genome_dirs <- function(data_root, marker) {
  subs <- list.files(data_root, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  subs <- subs[utils::file_test("-d", subs)]
  subs <- subs[!grepl("^(comparative|scripts|visualizations)$", basename(subs), ignore.case = TRUE)]
  genomes <- list()
  for (d in subs) {
    if (file.exists(file.path(d, marker))) {
      genomes[[length(genomes) + 1L]] <- list(id = basename(d), dir = d,
                                               marker_path = file.path(d, marker))
      next
    }
    inner <- list.files(d, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
    inner <- inner[utils::file_test("-d", inner)]
    for (dd in inner) {
      if (file.exists(file.path(dd, marker))) {
        genomes[[length(genomes) + 1L]] <- list(id = basename(d), dir = dd,
                                                 marker_path = file.path(dd, marker))
        break
      }
    }
  }
  genomes
}

.dnmb_comparative_discover_gbff_dirs <- function(data_root) {
  subs <- list.files(data_root, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  subs <- subs[utils::file_test("-d", subs)]
  subs <- subs[!grepl("^(comparative|scripts|visualizations)$", basename(subs), ignore.case = TRUE)]
  out <- list()
  for (d in subs) {
    direct <- list.files(d, pattern = "\\.(gbff|gbk|gb)$",
                         full.names = TRUE, ignore.case = TRUE)
    if (length(direct)) {
      out[[length(out) + 1L]] <- list(id = basename(d), dir = d, gbff = direct[[1]])
      next
    }
    inner <- list.files(d, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
    inner <- inner[utils::file_test("-d", inner)]
    for (dd in inner) {
      nested <- list.files(dd, pattern = "\\.(gbff|gbk|gb)$",
                           full.names = TRUE, ignore.case = TRUE)
      if (length(nested)) {
        out[[length(out) + 1L]] <- list(id = basename(d), dir = dd, gbff = nested[[1]])
        break
      }
    }
  }
  out
}

.dnmb_comparative_parse_organism <- function(analysis_dir) {
  gbff <- .dnmb_find_gbff_for_plot(analysis_dir)
  if (is.null(gbff) || !file.exists(gbff)) return(NA_character_)
  con <- base::file(gbff, "r")
  on.exit(close(con))
  in_record <- FALSE
  while (TRUE) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return(NA_character_)
    if (startsWith(line, "LOCUS")) {
      in_record <- TRUE
      next
    }
    if (in_record && grepl("^\\s*ORGANISM\\s+", line)) {
      return(trimws(sub("^\\s*ORGANISM\\s+", "", line)))
    }
    if (in_record && startsWith(line, "FEATURES")) return(NA_character_)
  }
}

.dnmb_comparative_genome_id <- function(dir_basename) {
  id <- sub("_done$", "", dir_basename)
  sub("\\.\\d+$", "", id)  # GCF_xxx.1 -> GCF_xxx
}

.dnmb_comparative_load_genbank_table <- function(analysis_dir) {
  xlsx_candidates <- list.files(analysis_dir, pattern = "_total\\.xlsx$", full.names = TRUE)
  if (!length(xlsx_candidates)) return(NULL)
  genes <- tryCatch(
    as.data.frame(
      readxl::read_excel(xlsx_candidates[[1]], sheet = "1.GenBank_table", guess_max = 10000),
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(genes) || !nrow(genes)) NULL else genes
}

# Generic auto-runner: for each folder with a genbank but missing the
# module's canonical marker file, invoke run_module_set() with the named
# database to produce it. Prepares the genbank_table from any existing
# *_total.xlsx and falls back to running Genbank_organizer() in place.
.dnmb_comparative_autorun_module <- function(data_root,
                                              module_db,
                                              module_marker_rel,
                                              ready_check = NULL,
                                              verbose = TRUE,
                                              module_cache_root = NULL,
                                              module_install = TRUE,
                                              module_cpu = NULL) {
  genomes <- .dnmb_comparative_discover_gbff_dirs(data_root)
  if (!length(genomes)) return(invisible(FALSE))
  cpu <- if (is.null(module_cpu)) .dnmb_default_cpu() else as.integer(module_cpu)[1]
  ran_any <- FALSE
  for (g in genomes) {
    ready <- if (!is.null(ready_check)) isTRUE(ready_check(g$dir))
             else file.exists(file.path(g$dir, module_marker_rel))
    if (ready) next
    if (verbose) message("[comparative] ", module_db, " missing for ", g$id, " — running module.")
    genes <- .dnmb_comparative_load_genbank_table(g$dir)
    if (is.null(genes)) {
      old_wd <- getwd()
      setwd(g$dir)
      tryCatch(
        Genbank_organizer(),
        error = function(e) {
          if (verbose) message("  Genbank_organizer() failed: ", conditionMessage(e))
        },
        finally = setwd(old_wd)
      )
      genes <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)
    }
    if (is.null(genes) || !is.data.frame(genes) || !nrow(genes)) {
      if (verbose) message("  unable to prepare gene table — skipping ", g$id)
      next
    }
    module_dir <- file.path(g$dir, paste0("dnmb_module_", tolower(module_db)))
    dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      run_module_set(
        db = module_db,
        genbank_table = genes,
        genbank = g$gbff,
        output_dir = module_dir,
        module_cache_root = module_cache_root,
        module_install = module_install,
        module_cpu = cpu,
        verbose = verbose
      ),
      error = function(e) {
        if (verbose) message("  ", module_db, " run failed: ", conditionMessage(e))
      }
    )
    ran_any <- ran_any || file.exists(file.path(g$dir, module_marker_rel))
  }
  invisible(ran_any)
}

# Given a long_df (File, subtype) and a File -> organism map, render the
# full comparative heatmap (pdf + xlsx) using the shared layout. Every
# module wrapper funnels through this.
.dnmb_comparative_render_heatmap <- function(long_df,
                                             organism_map,
                                             title,
                                             output_dir,
                                             output_file,
                                             color_palette = c("white", "#330066"),
                                             bar_color = "#4C1C7E",
                                             line_width = 1,
                                             line_col = "grey80",
                                             verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # organism_map is the source of truth for which genomes were analyzed,
  # regardless of whether each produced any hits. Analyzed-but-empty
  # genomes must render as zero-count rows so the user can distinguish
  # "ran and found nothing" from "not yet run".
  universe <- if (length(organism_map)) names(organism_map) else character(0)
  has_hits <- is.data.frame(long_df) && nrow(long_df) > 0
  if (has_hits) stopifnot(all(c("File", "subtype") %in% names(long_df)))

  if (!length(universe) && !has_hits) return(invisible(NULL))

  if (has_hits) {
    subtype_levels <- sort(unique(long_df$subtype))
    mat_df <- stats::aggregate(
      rep(1L, nrow(long_df)),
      by = list(File = long_df$File, subtype = long_df$subtype),
      FUN = sum
    )
    names(mat_df)[3] <- "count"
    mat_df$subtype <- factor(mat_df$subtype, levels = subtype_levels)
    wide_df <- stats::reshape(mat_df, idvar = "File", timevar = "subtype",
                              direction = "wide", sep = "__")
    names(wide_df) <- sub("^count__", "", names(wide_df))
    wide_df[is.na(wide_df)] <- 0L
    for (s in subtype_levels) {
      if (!s %in% names(wide_df)) wide_df[[s]] <- 0L
    }
    wide_df <- wide_df[, c("File", subtype_levels), drop = FALSE]
  } else {
    subtype_levels <- character(0)
    wide_df <- data.frame(File = character(0), stringsAsFactors = FALSE)
  }

  # Pad the full analyzed universe: every genome in organism_map gets a
  # row, even if no subtypes were detected. Empty rows render as a blank
  # heatmap line + zero-length diversity bar.
  missing_files <- setdiff(universe, wide_df$File)
  if (length(missing_files)) {
    add_df <- data.frame(File = missing_files, stringsAsFactors = FALSE)
    for (s in subtype_levels) add_df[[s]] <- 0L
    wide_df <- rbind(wide_df, add_df)
  }
  if (!length(subtype_levels)) {
    if (verbose) message("No subtype hits across any analyzed genome — nothing to plot.")
    return(invisible(NULL))
  }

  organism <- organism_map[wide_df$File]
  organism[is.na(organism) | !nzchar(organism)] <- wide_df$File[is.na(organism) | !nzchar(organism)]
  genus <- vapply(strsplit(organism, "\\s+"),
                  function(x) if (length(x)) x[[1]] else NA_character_, character(1))
  species <- vapply(strsplit(organism, "\\s+"),
                    function(x) if (length(x) >= 2L) x[[2]] else NA_character_, character(1))

  wide_df$SOURCE <- organism
  wide_df$Genus <- factor(genus)
  wide_df$Species <- factor(species)
  wide_df <- wide_df[order(wide_df$Genus, wide_df$Species, wide_df$File), , drop = FALSE]

  counts_only <- as.matrix(wide_df[, subtype_levels, drop = FALSE])
  storage.mode(counts_only) <- "numeric"
  diversity <- rowSums(counts_only > 0)
  diversity_ratio <- diversity / max(1L, length(subtype_levels))

  xlsx_path <- file.path(output_dir, sub("\\.pdf$", ".xlsx", output_file))
  summary_df <- data.frame(
    File = wide_df$File,
    SOURCE = wide_df$SOURCE,
    Genus = as.character(wide_df$Genus),
    Species = as.character(wide_df$Species),
    wide_df[, subtype_levels, drop = FALSE],
    diversity = diversity,
    diversity_ratio = diversity_ratio,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  tryCatch(
    openxlsx::write.xlsx(summary_df, xlsx_path, rowNames = FALSE, colNames = TRUE, overwrite = TRUE),
    error = function(e) if (verbose) message("Could not write xlsx: ", conditionMessage(e))
  )

  mat <- counts_only
  rownames(mat) <- wide_df$SOURCE

  unique_genus <- levels(wide_df$Genus)
  n_genus <- length(unique_genus)
  hue_start <- 15
  hues <- (hue_start + (seq_len(n_genus) - 1L) * (360 / n_genus)) %% 360
  genus_colors <- grDevices::hcl(h = hues, c = 72, l = 58, fixup = TRUE)
  names(genus_colors) <- unique_genus

  species_colors <- list()
  for (gi in seq_along(unique_genus)) {
    g <- unique_genus[[gi]]
    g_species <- unique(as.character(wide_df$Species[wide_df$Genus == g & !is.na(wide_df$Species)]))
    n_sp <- length(g_species)
    if (n_sp) {
      sp_l <- if (n_sp == 1L) 85 else seq(90, 72, length.out = n_sp)
      sp_c <- if (n_sp == 1L) 30 else seq(22, 42, length.out = n_sp)
      pal <- grDevices::hcl(h = hues[[gi]], c = sp_c, l = sp_l, fixup = TRUE)
      names(pal) <- g_species
      species_colors <- c(species_colors, pal)
    }
  }
  species_colors <- unlist(species_colors)
  species_colors <- species_colors[!is.na(names(species_colors)) & nzchar(names(species_colors))]

  row_anno <- data.frame(Genus = wide_df$Genus, Species = wide_df$Species)
  heatmap_colors <- list(Genus = genus_colors, Species = species_colors)

  transparent_col_fun <- circlize::colorRamp2(
    c(min(mat, na.rm = TRUE), max(1, max(mat, na.rm = TRUE))),
    c("transparent", "transparent")
  )
  gradient_col_fun <- circlize::colorRamp2(
    c(min(mat, na.rm = TRUE), max(1, max(mat, na.rm = TRUE))),
    color_palette
  )
  max_size <- max(1, max(mat, na.rm = TRUE))

  diversity_mat <- matrix(diversity, ncol = 1L)
  rownames(diversity_mat) <- wide_df$SOURCE
  bar_anno <- ComplexHeatmap::rowAnnotation(
    Diversity = ComplexHeatmap::anno_barplot(
      diversity_mat,
      width = grid::unit(3, "cm"),
      gp = grid::gpar(fill = bar_color, col = NA),
      border = TRUE,
      axis_param = list(gp = grid::gpar(fontsize = 12))
    ),
    annotation_name_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  make_italics <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
  make_plains <- function(x) as.expression(lapply(x, function(y) bquote(plain(.(y)))))
  "%c%" <- function(x, y) parse(text = paste(x, "~", y, sep = " "))
  first_two <- paste(
    vapply(strsplit(rownames(mat), "\\s+"),
           function(x) if (length(x)) x[[1]] else "", character(1)),
    vapply(strsplit(rownames(mat), "\\s+"),
           function(x) if (length(x) >= 2L) x[[2]] else "", character(1)),
    sep = " "
  )
  remainder <- trimws(substring(rownames(mat), nchar(first_two) + 1))
  row_labels <- make_italics(first_two) %c% make_plains(remainder)

  count_min <- 0
  count_max <- max(1, max(mat, na.rm = TRUE))
  count_breaks <- seq(count_min, count_max, length.out = length(color_palette))
  count_legend <- ComplexHeatmap::Legend(
    title = "Count (n)",
    col_fun = circlize::colorRamp2(breaks = count_breaks, colors = color_palette),
    direction = "vertical",
    title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 12)
  )
  genus_legend <- ComplexHeatmap::Legend(
    labels = names(genus_colors),
    title = "Genus",
    type = "points", pch = 16,
    background = NA,
    legend_gp = grid::gpar(col = genus_colors),
    title_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  species_legend <- if (length(species_colors)) {
    ComplexHeatmap::Legend(
      labels = names(species_colors),
      title = "Species",
      type = "points", pch = 16,
      background = NA,
      legend_gp = grid::gpar(col = species_colors),
      title_gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  } else NULL
  combined_legend <- if (!is.null(species_legend)) {
    ComplexHeatmap::packLegend(genus_legend, species_legend,
                               direction = "horizontal", gap = grid::unit(8, "mm"))
  } else genus_legend

  cell_size_mm    <- 8
  body_width_mm   <- ncol(mat) * cell_size_mm
  body_height_mm  <- nrow(mat) * cell_size_mm
  left_strip_mm   <- 14
  bar_width_mm    <- 30
  row_label_mm    <- max(90, 10 + max(nchar(rownames(mat)), na.rm = TRUE) * 2.8)
  col_label_mm    <- max(45, 10 + max(nchar(colnames(mat)), na.rm = TRUE) * 2.8)
  title_mm        <- 14
  legend_block_mm <- 55
  side_pad_mm     <- 10

  pdf_width_in  <- (side_pad_mm + left_strip_mm + body_width_mm +
                      bar_width_mm + row_label_mm + side_pad_mm) / 25.4
  pdf_height_in <- (title_mm + body_height_mm + col_label_mm +
                      legend_block_mm + side_pad_mm) / 25.4

  pdf_path <- file.path(output_dir, output_file)
  grDevices::pdf(pdf_path, width = pdf_width_in, height = pdf_height_in)

  dot_hm <- ComplexHeatmap::Heatmap(
    mat,
    name = "subtype",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_dend = FALSE, show_row_names = TRUE, border = TRUE,
    width  = grid::unit(body_width_mm, "mm"),
    height = grid::unit(body_height_mm, "mm"),
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 12),
    row_labels = row_labels,
    left_annotation = ComplexHeatmap::rowAnnotation(
      df = row_anno, col = heatmap_colors,
      annotation_name_gp = grid::gpar(fontsize = 12, fontface = "italic")
    ),
    right_annotation = bar_anno,
    col = transparent_col_fun
  )

  ComplexHeatmap::draw(
    dot_hm,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    padding = grid::unit(c(legend_block_mm, side_pad_mm, title_mm, side_pad_mm), "mm")
  )

  ComplexHeatmap::decorate_heatmap_body("subtype", {
    # Grid lines first so circles render on top of them.
    for (i in seq_len(nrow(mat))) {
      grid::grid.lines(x = c(0, 1),
                       y = grid::unit((nrow(mat) - i + 0.5) / nrow(mat), "npc"),
                       gp = grid::gpar(lty = 1, lwd = line_width, col = line_col))
    }
    for (j in seq_len(ncol(mat))) {
      grid::grid.lines(x = grid::unit((j - 0.5) / ncol(mat), "npc"),
                       y = c(0, 1),
                       gp = grid::gpar(lty = 1, lwd = line_width, col = line_col))
    }
    # Circles (sized by count) with the count written inside.
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        v <- mat[i, j]
        if (!is.finite(v) || v <= 0) next
        cx <- grid::unit((j - 0.5) / ncol(mat), "npc")
        cy <- grid::unit((nrow(mat) - i + 0.5) / nrow(mat), "npc")
        r_mm <- sqrt(v / max_size) * 4
        grid::grid.circle(cx, cy, grid::unit(r_mm, "mm"),
                          gp = grid::gpar(fill = gradient_col_fun(v), col = NA))
        grid::grid.text(
          as.character(v), cx, cy,
          gp = grid::gpar(
            fontsize = max(6, min(10, r_mm * 2.2)),
            fontface = "bold",
            col = if (v >= max_size * 0.5) "white" else "black"
          )
        )
      }
    }
  })

  ComplexHeatmap::decorate_annotation("Diversity", {
    grid::pushViewport(grid::viewport(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0, "npc") - grid::unit(20, "mm"),
      width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
      just = c("center", "top"), clip = "off"
    ))
    ComplexHeatmap::draw(count_legend,
                         x = grid::unit(0.5, "npc"), y = grid::unit(1, "npc"),
                         just = c("center", "top"))
    grid::popViewport()
  })

  row_label_x_offset <- grid::unit(body_width_mm + bar_width_mm + 14, "mm")
  ComplexHeatmap::decorate_heatmap_body("subtype", {
    grid::pushViewport(grid::viewport(
      x = row_label_x_offset,
      y = grid::unit(0, "npc") - grid::unit(6, "mm"),
      width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
      just = c("left", "top"), clip = "off"
    ))
    ComplexHeatmap::draw(combined_legend,
                         x = grid::unit(0, "npc"), y = grid::unit(1, "npc"),
                         just = c("left", "top"))
    grid::popViewport()
  })

  ComplexHeatmap::decorate_annotation("Genus", {
    grid::grid.text(
      title,
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc") + grid::unit(6, "mm"),
      just = c("left", "bottom"),
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )
  })

  grDevices::dev.off()
  if (verbose) message("Heatmap saved: ", pdf_path)

  invisible(list(pdf = pdf_path, xlsx = xlsx_path, matrix = mat, genomes = wide_df$File))
}
