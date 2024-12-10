#' Generate a Heatmap for DefenseFinder Results
#'
#' This function generates a heatmap visualization from DefenseFinder results, including annotations
#' for diversity metrics and customizable aesthetics such as color gradients and barplot styles.
#'
#' @param analysis_dir Character. The directory containing DefenseFinder result files. Default is the current working directory.
#' @param output_file Character. The path to save the generated heatmap as a PDF. Default is `"DefenseFinder_Heatmap.pdf"`.
#' @param color_palette Character vector. A two-color gradient palette for the heatmap cells. Default is `c("white", "#330066")`.
#' @param bar_color Character. The color of the barplot showing diversity metrics. Default is `"#4C1C7E"`.
#' @param gradient_length Integer. Number of steps in the color gradient for the heatmap. Default is `100`.
#' @param line_width Numeric. Width of grid lines in the heatmap. Default is `1`.
#' @param line_col Character. Color of grid lines in the heatmap. Default is `"grey80"`.
#' @param debug Logical. Whether to display debug messages during function execution. Default is `TRUE`.
#'
#' @return A `ComplexHeatmap` object representing the heatmap visualization. The heatmap is also saved to a PDF file.
#'
#' @details
#' This function processes DefenseFinder results to generate a heatmap visualization. It links GenBank features
#' from the `Genome_summary_df` dataset and annotates rows with genus and species colors. The barplot annotation
#' shows diversity metrics per sample.
#'
#' Customization options include:
#' - Heatmap cell color gradient (`color_palette`).
#' - Barplot fill color (`bar_color`).
#' - Font sizes and line aesthetics (`line_width`, `line_col`).
#'
#' The function requires `ComplexHeatmap`, `circlize`, `grid`, `reshape2`, and `tidyverse` packages.
#'
#' @examples
#' # Generate a heatmap using default settings:
#' DefenseFinder_Heatmap()
#'
#' # Customize the color palette and barplot color:
#' DefenseFinder_Heatmap(
#'     color_palette = c("white", "#FF5733"),
#'     bar_color = "#33A1FF",
#'     debug = FALSE
#' )
#'
#' @import ComplexHeatmap circlize grid reshape2 tidyverse
#' @export

DefenseFinder_Heatmap <- function(
    analysis_dir = getwd(),
    output_file = "DefenseFinder_Heatmap.pdf",
    color_palette = c("white", "#330066"),
    bar_color = "#4C1C7E",
    gradient_length = 100,
    line_width = 1,
    line_col = "grey80",
    debug = TRUE
) {
  # Helper functions
  "%c%" <- function(x, y) parse(text = paste(x, "~", y, sep = " "))

  make_italics <- function(x) {
    as.expression(lapply(x, function(y) bquote(italic(.(y)))))
  }

  make_plains <- function(x) {
    as.expression(lapply(x, function(y) bquote(plain(.(y)))))
  }

  make_symmetric <- function(mat) {
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    return(mat)
  }

  # Check for Genome_summary_df in the environment
  if (!exists("Genome_summary_df", envir = .GlobalEnv)) {
    if (debug) message("Genome_summary_df not found in the environment. Generating using Genbank_info_combiner...")
    Genome_summary_df <- Genbank_info_combiner()
    assign("Genome_summary_df", Genome_summary_df, envir = .GlobalEnv)
  } else {
    if (debug) message("Genome_summary_df found in the environment.")
    Genome_summary_df <- get("Genome_summary_df", envir = .GlobalEnv)
  }

  # Internal genus color palette
  base_genus_colors <- c(
    "#ff4947", "#f67e7d", "#f69c7d", "#f6b07d", "#AF967D", "#af7a46","#E63946",
    "#dda46b", "#f2c9a1", "#f8c779", "#ffc661", "#ffb532", "#f59e03","#FFE000",
    "#f5b247", "#f8c924", "#f3e75b", "#fff098", "#c6dda6", "#bad780","#1D3557",
    "#919f7f", "#95b46a", "#A3BE8C", "#80d7c6", "#4fc6d0", "#4faad0",
    "#5c9ce4", "#337cd6", "#5565ca", "#5E81AC", "#AD8CAE", "#bc97ab",
    "#bb84a1", "#cb72a1", "#cc5293", "#a04876", "#484860", "#184860",
    "#c0c0c0", "#E66101", "#FDB863", "#CA0020", "#F4A582", "#545863",
    "#4DAC26", "#B8E186", "#0571B0", "#92C5DE", "#B2ABD2", "#5E3C99",
    "#80CDC1", "#018571", "#DFC27D", "#BABABA", "#404040", "#0571B0"
  )

  analysis_list <- list.files()[grep("finder_systems\\.tsv$",list.files())]
  defense_finder <- vector("list", length = length(analysis_list))

  for (i in 1:length(analysis_list)){
    #dbCAN[[analysis_list[i]]][[dbCAN_list[j]]] <- read_tsv(file=paste(analysis_list[i], dbCAN_list[j], sep="_")) # when using prefix
    file_path <- paste0(getwd(), "/", analysis_list[i])
    defense_finder[[i]] <- readr::read_delim(file=file_path) # when using no prefix
  }

  if (any(grepl("^(GCA_|GCF_)", analysis_list))){
    names(defense_finder) <- analysis_list %>% qdap::beg2char(.,char="\\.", noc =2) # naming check
  } else{
    names(defense_finder) <- analysis_list %>% qdap::beg2char(.,char="\\.", noc =2) # naming check
  }

  # data.frame
  defense_finder_df <- defense_finder %>% tibble::enframe(name="File") %>% unnest(cols=dplyr::everything())

  # calculate diversity
  subtype_list <- defense_finder_df %>% dplyr::select(`subtype`) %>% dplyr::distinct() %>% pull()

  # re-organize
  defense_finder_df <- defense_finder_df %>% reshape2::dcast(., formula = File ~ `subtype`, fun.aggregate = length, value.var = "subtype")

  # matrix
  defense_finder_df <- defense_finder_df %>% dplyr::mutate(diversity=rowSums(select(., -File) != 0))
  defense_finder_df <- defense_finder_df %>% dplyr::mutate(diversity_ratio = diversity / length(subtype_list))

  # database link to genbank feature
  defense_finder_df <- merge(x=defense_finder_df, y=Genome_summary_df, by.x="File", by.y="File", all=TRUE)
  defense_finder_df <- defense_finder_df %>% dplyr::select('File', 'SOURCE', subtype_list, "diversity", "diversity_ratio")

  defense_finder_df <- defense_finder_df %>% dplyr::mutate(
    Genus = stringr::word(SOURCE, 1) %>% as.factor(),
    Species = stringr::word(SOURCE, 2) %>% as.factor()
  )

  # arrange column sequence
  defense_finder_df <- arrange.vars(defense_finder_df, c("Genus"=which(colnames(defense_finder_df) == "SOURCE")+1,
                                                         "Species"=which(colnames(defense_finder_df) == "SOURCE")+2)
  )
  defense_finder_df <- defense_finder_df %>% dplyr::arrange(Genus, Species)

  # save to xlsx
  file_name <- paste("Table. DefenseFinder",length(defense_finder), "total.xlsx", sep="_")
  write.xlsx(defense_finder_df, file_name, rownames=FALSE,colnames=TRUE)


  # visualization
  mat <- as.matrix(defense_finder_df %>% dplyr::select(sort(subtype_list)))
  rownames(mat) <- defense_finder_df$SOURCE

  diversity <- as.matrix(defense_finder_df$diversity)
  rownames(diversity) <- defense_finder_df$SOURCE

  row_anno <- data.frame(Genus = defense_finder_df$Genus, Species = defense_finder_df$Species)
  row_anno <- as.data.frame(lapply(row_anno, factor))

  defense_finder_df$Genus%>% unique()
  defense_finder_df$Species%>% unique()

  # Generate Genus colors dynamically (Random picking)
  unique_genus <- levels(defense_finder_df$Genus)
  if (length(unique_genus) <= length(base_genus_colors)) {
    genus_colors <- sample(base_genus_colors, length(unique_genus))
  } else {
    extra_colors_needed <- length(unique_genus) - length(base_genus_colors)
    random_extra_colors <- grDevices::colors()[sample(grDevices::colors(), extra_colors_needed)]
    genus_colors <- c(base_genus_colors, random_extra_colors)
    genus_colors <- sample(genus_colors, length(unique_genus))
  }
  names(genus_colors) <- unique_genus
  genus_colors <- unlist(genus_colors)

  if (debug) {
    message("Assigned Genus colors:")
    print(genus_colors)
  }

  # Generate species-specific colors grouped by Genus
  species_colors <- list()

  for (genus in unique_genus) {
    genus_species <- defense_finder_df %>%
      dplyr::filter(Genus == genus & !is.na(Species)) %>% # NA 제거
      dplyr::pull(Species) %>%
      unique()

    if (length(genus_species) > 0) {
      species_palette <- grDevices::colorRampPalette(c(genus_colors[genus], "white"))(length(genus_species))
      names(species_palette) <- genus_species
      species_colors <- c(species_colors, species_palette)
    }
  }

  species_colors <- unlist(species_colors)

  if (debug) {
    message("Species colors before cleanup:")
    print(species_colors)
  }

  # Ensure NA and invalid names are handled
  species_colors <- species_colors[!is.na(names(species_colors)) & names(species_colors) != ""]

  # Validate that species_colors contains only valid colors
  if (any(is.na(species_colors)) || any(!grepl("^#", species_colors))) {
    stop("Invalid color values detected in species_colors after cleanup.")
  }

  heatmap_colors <- list(
    "Genus" = genus_colors,
    "Species" = species_colors
  )

  if (debug) {
    message("Final species colors:")
    print(species_colors)
  }

  # Set up heatmap gradient for identity values
  heatmap_gradient <- grDevices::colorRampPalette(color_palette)(gradient_length)

  # ComplexHeatmap
  mat <- as.matrix(defense_finder_df %>% select(sort(subtype_list)))
  rownames(mat) <- defense_finder_df$SOURCE

  max_size <- max(mat, na.rm = TRUE)

  # Dot plot generation
  transparent_col_fun <- circlize::colorRamp2(c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)), c("transparent", "transparent"))
  gradient_col_fun <- circlize::colorRamp2(c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)), color_palette)

  bar_anno <- rowAnnotation(Diversity = anno_barplot(diversity,
                                                     width = unit(3, "cm"),
                                                     gp = gpar(fill = bar_color, col = NA),
                                                     border = TRUE,
                                                     axis_param = list(gp = gpar(fontsize = 12))),
                                                     annotation_name_gp = gpar(fontsize = 14, fontface = "bold"))





  # Set up heatmap gradient for identity values
  heatmap_gradient <- grDevices::colorRampPalette(color_palette)(gradient_length)

  # Generate labels with italics and plain text
  row_labels <- paste(
    stringr::word(rownames(mat), 1), # Extract the first word
    stringr::word(rownames(mat), 2), # Extract the second word
    sep = " "
  ) %>% make_italics() %c%
    (rownames(mat) %>%
       stringr::str_remove("^\\S+\\s+\\S+\\s*") %>% # Remove first two words and their spaces
       stringr::str_squish() %>% # Remove any extra spaces
       make_plains())

  # Create legends
  genus_legend <- ComplexHeatmap::Legend(
    labels = names(genus_colors),
    title = "Genus",
    type = "points",
    pch = 16,
    legend_gp = grid::gpar(col = genus_colors),
    title_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  species_legend <- ComplexHeatmap::Legend(
    labels = names(species_colors),
    title = "Species",
    type = "points",
    pch = 16,
    legend_gp = grid::gpar(col = species_colors),
    title_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  # dot_hm not here
  count_min <- min(mat, na.rm = TRUE)
  count_max <- max(mat, na.rm = TRUE)
  count_breaks <- seq(count_min, count_max, length.out = length(color_palette))

  count_legend <- ComplexHeatmap::Legend(
    title = "Count (n)",
    col_fun = circlize::colorRamp2(
      breaks = count_breaks,
      colors = color_palette
    ),
    direction = "vertical",
    title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 12)
  )

  combined_legend <- ComplexHeatmap::packLegend(
    genus_legend,
    species_legend,
    gap = grid::unit(5, "mm")
  )

  pdf(output_file, width = 14, height = 12)

  dot_hm <- ComplexHeatmap::Heatmap(mat, name = "subtype",
                                    cluster_rows = FALSE,
                                    cluster_columns=FALSE,
                                    show_row_dend = FALSE,
                                    show_row_names = TRUE,
                                    border = TRUE,
                                    row_names_gp = gpar(fontsize = 12),
                                    column_names_gp = gpar(fontsize = 12),
                                    row_labels =row_labels,
                                    left_annotation = rowAnnotation(
                                      df = row_anno,
                                      col = heatmap_colors,
                                      annotation_name_gp = gpar(fontsize = 12, fontface = "italic")
                                    ),
                                    right_annotation = bar_anno,
                                    col = transparent_col_fun,
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      size <- sqrt(mat[i, j] / max_size) * unit(5, "mm")
                                      numeric_size <- sqrt(mat[i, j] / max_size)*4
                                      circle_size <- unit(numeric_size, "mm")
                                      color <- gradient_col_fun(mat[i, j])
                                      if(mat[i, j] > 0) {
                                        grid.circle(x, y, circle_size, gp = gpar(fill = color, col = NA))
                                      }
                                    })

  {
    row_idx <- 1:nrow(mat) # target
    btm <- 1 - (row_idx / nrow(mat))
    top <- btm + (1/nrow(mat))

    col_idx <- 1:nrow(mat) # target
    left <- 1 - (col_idx / ncol(mat))
    right <- left + (1/ncol(mat))

    # Add boxes around our rows
    line_col <- 'grey80'
    line_width <- 1
  }

  ComplexHeatmap::draw(dot_hm,
                       show_heatmap_legend = FALSE,
                       show_annotation_legend = FALSE
  )

  decorate_heatmap_body("subtype", {
    # Draw horizontal lines at row centers
    for (i in seq_len(nrow(mat))) {
      grid.lines(
        x = c(0, 1),  # Left to right within the heatmap body
        y = unit((nrow(mat) - i + 0.5) / nrow(mat), "npc"),  # Center of each row
        gp = gpar(lty = 1, lwd = line_width, col = line_col)
      )
    }

    # Draw vertical lines at column centers
    for (j in seq_len(ncol(mat))) {
      grid.lines(
        x = unit((j - 0.5) / ncol(mat), "npc"),  # Center of each column
        y = c(0, 1),  # Bottom to top within the heatmap body
        gp = gpar(lty = 1, lwd = line_width, col = line_col)
      )
    }
  })

  grid::pushViewport(grid::viewport(x = unit(0.81, "npc"), y = unit(0.635, "npc"), width = 0.1, height = 1.0, just = c("left", "top")))
  ComplexHeatmap::draw(count_legend)
  grid::popViewport()

  ComplexHeatmap::draw(combined_legend, x = unit(0.9, "npc"), y = unit(0.1835, "npc"), just = c("left", "top"))

  grid::grid.text(
    "DefenseFinder", x = 0.01, y = 1.01, just = c("left", "top"),
    gp = grid::gpar(fontsize = 16, fontface = "bold")
  )
  # Close the PDF device
  dev.off()
  if (debug) message("Heatmap successfully saved to: ", output_file)

  # Return the heatmap object for use in RStudio
  return(dot_hm)
}
