#' Generate a heatmap for codon usage and tRNA copy number.
#'
#' This function creates a heatmap that visualizes codon usage and the corresponding tRNA copy number distribution.
#'
#' @param codon_usage A data frame containing codon usage data.
#' @param tRNA_distribution A data frame containing tRNA copy number data.
#' @param output_pdf A string specifying the name of the output PDF file.
#' @param width A numeric value specifying the width of the output PDF in inches.
#' @param height A numeric value specifying the height of the output PDF in inches.
#' @param save_pdf A logical value indicating whether to save the plot as a PDF (default is TRUE).
#' @return A heatmap figure combining codon usage and tRNA copy number information.
#' @export

.dnmb_draw_aa_block_labels <- function(aa_blocks = NULL, fontsize = 8) {
  if (is.null(aa_blocks)) {
    if (!exists("Codon_usage_AA_blocks", envir = .GlobalEnv)) {
      return(invisible(NULL))
    }
    aa_blocks <- get("Codon_usage_AA_blocks", envir = .GlobalEnv)
  }

  if (is.null(aa_blocks) || !nrow(aa_blocks)) {
    return(invisible(NULL))
  }

  ComplexHeatmap::decorate_annotation("AA", {
    for (idx in seq_len(nrow(aa_blocks))) {
      grid::grid.text(label = aa_blocks$AA[idx],
                      x = grid::unit(0.5, "npc"),
                      y = grid::unit(aa_blocks$center_npc[idx], "npc"),
                      rot = 90,
                      gp = grid::gpar(col = aa_blocks$text_col[idx],
                                      fontsize = fontsize,
                                      fontface = "bold"))
    }
  })

  invisible(NULL)
}

codon_usage_tRNA_heatmap_generator <- function(codon_usage = NULL,
                                               tRNA_distribution = NULL,
                                               output_pdf = "codon_usage_tRNA_heatmap.pdf",
                                               width = 10,
                                               height = 10,
                                               save_pdf = FALSE) {
  # Load necessary libraries
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)

  # Fetch codon_usage and tRNA_distribution from the environment if they are NULL
  if (is.null(codon_usage)) {
    if (exists("codon_usage", envir = .GlobalEnv)) {
      codon_usage <- get("codon_usage", envir = .GlobalEnv)
      message("codon_usage fetched from the global environment.")
    } else {
      stop("codon_usage is not provided and not found in the global environment.")
    }
  }

  if (is.null(tRNA_distribution)) {
    if (exists("tRNA_distribution", envir = .GlobalEnv)) {
      tRNA_distribution <- get("tRNA_distribution", envir = .GlobalEnv)
      message("tRNA_distribution fetched from the global environment.")
    } else {
      stop("tRNA_distribution is not provided and not found in the global environment.")
    }
  }

  # Define AA and chemical properties
  AA <- c('S', 'T', 'C', 'Y', 'N', 'Q', 'G', 'A', 'V', 'L', 'I', 'F', 'W', 'M', 'P', 'H', 'K', 'R', 'D', 'E', 'O', 'U', 'B', 'Z', 'X', 'J')
  Chemical_properties <- c('polar', 'polar', 'polar', 'polar', 'polar', 'polar', 'nonpolar', 'nonpolar', 'nonpolar', 'nonpolar', 'nonpolar', 'nonpolar', 'nonpolar', 'nonpolar', 'nonpolar', 'basic', 'basic', 'basic', 'acidic', 'acidic', '_', '_', '_', '_', '_', "_")
  aa_chemical_prop. <- data.frame("AA"=AA[1:20], "Chemical_property"=Chemical_properties[1:20])

  aa_palette <- c("*" = "#cbd2d6",
                  "S" = "#f5cadb",
                  "T" = "#eea1be",
                  "C" = "#e678a1",
                  "Y" = "#df4e85",
                  "N" = "#d62668",
                  "Q" = "#ac1f54",
                  "G" = "#d0f0e8",
                  "A" = "#aae5d7",
                  "V" = "#85d9c5",
                  "L" = "#5fceb4",
                  "I" = "#3bc1a1",
                  "F" = "#309c82",
                  "W" = "#247663",
                  "M" = "#185143",
                  "P" = "#0d2b24",
                  "H" = "#8fb5cf",
                  "K" = "#6c9ec1",
                  "R" = "#3c6d8f",
                  "D" = "#f1dd73",
                  "E" = "#ebce31")

  .dnmb_label_contrast <- function(cols) {
    cols <- as.character(cols)
    vapply(cols, function(col) {
      rgb <- grDevices::col2rgb(col) / 255
      luminance <- (0.299 * rgb[1, 1]) + (0.587 * rgb[2, 1]) + (0.114 * rgb[3, 1])
      if (luminance > 0.65) "#1F1F1F" else "white"
    }, character(1))
  }

  # Merge codon usage with chemical properties
  codon_usage <- merge(x = codon_usage, y = aa_chemical_prop., by = "AA", all.x = TRUE)
  annotation_table <- codon_usage %>% dplyr::select(-eff, -freq, -RSCU)
  annotation_table <- annotation_table %>% group_by(Chemical_property) %>% arrange(Chemical_property, AA) %>% ungroup() %>% as.data.frame()
  rownames(annotation_table) <- annotation_table$codon

  custom_codon_order <- annotation_table %>% pull(codon)
  aa_vec <- as.character(annotation_table$AA)
  aa_runs <- rle(aa_vec)
  aa_block_starts <- cumsum(c(1L, head(aa_runs$lengths, -1L)))
  aa_block_ends <- cumsum(aa_runs$lengths)
  aa_blocks <- data.frame(
    AA = aa_runs$values,
    start = aa_block_starts,
    end = aa_block_ends,
    stringsAsFactors = FALSE
  )
  aa_blocks$center_npc <- 1 - (((aa_blocks$start + aa_blocks$end) / 2) - 0.5) / nrow(annotation_table)
  aa_blocks$fill <- unname(aa_palette[aa_blocks$AA])
  aa_blocks$text_col <- .dnmb_label_contrast(aa_blocks$fill)
  assign("Codon_usage_AA_blocks", aa_blocks, envir = .GlobalEnv)

  # Define custom AA order
  custom_AA_order <- c("E", "D", "R", "K", "H", "P", "M", "W", "F", "I", "L", "V", "A", "G", "Q", "N", "Y", "C", "T", "S", "*")
  annotation_table$AA <- factor(annotation_table$AA, levels = custom_AA_order)

  # Rearrange codon usage plot
  codon_usage_plot <- codon_usage %>% dplyr::select(codon, RSCU)
  rownames(codon_usage_plot) <- codon_usage$codon
  codon_usage_plot <- codon_usage_plot %>% dplyr::select(-codon)
  codon_usage_plot <- codon_usage_plot[custom_codon_order, , drop = FALSE]

  # Rearrange tRNA distribution plot
  tRNA_distribution_plot <- tRNA_distribution %>% dplyr::select(codon, codon_count)
  rownames(tRNA_distribution_plot) <- tRNA_distribution_plot$codon
  tRNA_distribution_plot <- tRNA_distribution_plot %>% dplyr::select(-codon)
  tRNA_distribution_plot <- tRNA_distribution_plot[custom_codon_order, , drop = FALSE]

  palettelength <- 100

  # Define custom color palette with white for value 1
  custom_colors <- circlize::colorRamp2(c(min(as.matrix(codon_usage_plot)), 1, max(as.matrix(codon_usage_plot))),
                                        c("#337cd6", "white", "#FE3D62"))

  # Create the heatmap
  my_heatmap <-  ComplexHeatmap::Heatmap(as.matrix(codon_usage_plot),
                                         col = custom_colors,
                                         name = "Z-score",
                                         heatmap_height = unit(4.0, "mm")*nrow(codon_usage_plot),
                                         heatmap_width = unit(0.8, "mm")*nrow(codon_usage_plot),
                                         cluster_rows = FALSE,
                                         row_order = custom_codon_order,
                                         show_row_names = TRUE,
                                         show_column_names = TRUE,
                                         left_annotation = rowAnnotation(
                                           df = annotation_table %>% dplyr::select(Chemical_property, AA),
                                           col = list("Chemical_property" = c("acidic" = "#ede9d0",
                                                                              "basic" = "#d0d4ed",
                                                                              "nonpolar" = "#d0edd9",
                                                                              "polar" = "#edd0e4",
                                                                              "NA" = "#ECEFF1"),
                                                      "AA" = aa_palette
                                           ),
                                           annotation_legend_param = list(
                                             Chemical_property = list(title = "Chemical\nproperty")
                                           )
                                         ),
                                         right_annotation = rowAnnotation(
                                           "No. of\ntRNA" = anno_barplot(as.matrix(tRNA_distribution_plot)[custom_codon_order, , drop = FALSE],
                                                                  gp = gpar(fill = "#595959", col = NA),
                                                                  border = FALSE)
                                         ),
                                         column_title = "Codon usage and tRNA distribution",
                                         column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         column_title_side = "top"
  )

  # Conditionally save the heatmap as a PDF if save_pdf is TRUE
  if (save_pdf) {
    pdf(output_pdf, width = width, height = height)
    ComplexHeatmap::draw(my_heatmap)
    .dnmb_draw_aa_block_labels(fontsize = 8)
    dev.off()
    message("Heatmap saved as ", output_pdf)
  } else {
    message("PDF not saved. To save, set save_pdf = TRUE.")
  }

  # Assign the result to 'Codon_usage_tRNA_plot' in the global environment
  assign("Codon_usage_tRNA_plot", my_heatmap, envir = .GlobalEnv)
  message("The heatmap has been saved to the R environment variable 'Codon_usage_tRNA_plot'.")
}
