#' Generate a complex plot with heatmap, seqlogo, and histogram
#'
#' @param output_pdf The output path for the PDF file.
#' @param pdf_width The width of the output PDF file (default is 12).
#' @param pdf_height The height of the output PDF file (default is 12).
#' @return A PDF file with the arranged plots and labels.
#' @export
Figure_generator <- function(output_pdf = "arranged_plot_with_labels.pdf",
                             pdf_width = 12, pdf_height = 12) {

  library(gridExtra)
  library(ComplexHeatmap)
  library(grid)

  # Check if Codon_usage_tRNA_plot exists in the environment
  if (exists("Codon_usage_tRNA_plot", envir = .GlobalEnv)) {
    Codon_usage_tRNA_plot <- get("Codon_usage_tRNA_plot", envir = .GlobalEnv)
  } else {
    stop("Codon_usage_tRNA_plot not found in the environment.")
  }

  # Check if RBS_seqlogo exists in the environment
  if (exists("RBS_seqlogo", envir = .GlobalEnv)) {
    RBS_seqlogo <- get("RBS_seqlogo", envir = .GlobalEnv)
  } else {
    stop("RBS_seqlogo not found in the environment.")
  }

  # Check if spacer_histogram exists in the environment
  if (exists("spacer_histogram", envir = .GlobalEnv)) {
    spacer_histogram <- get("spacer_histogram", envir = .GlobalEnv)
  } else {
    stop("spacer_histogram not found in the environment.")
  }

  # Reduce plot margins if using ggplot objects
  RBS_seqlogo <- RBS_seqlogo + theme(plot.margin = unit(c(3, 0, 3, 0), "cm"))
  spacer_histogram <- spacer_histogram + theme(plot.margin = unit(c(3, 0, 3, 0), "cm"))

  # Adjust padding in ComplexHeatmap
  heatmap_grob <- grid.grabExpr(draw(Codon_usage_tRNA_plot, padding = unit(c(1, 0, 1, 0), "cm")))

  # Create the layout with grid.arrange, ensuring widths and heights are compatible with layout_matrix
  arranged_plot <- grid.arrange(
    heatmap_grob,
    RBS_seqlogo,
    spacer_histogram,
    layout_matrix = rbind(c(1, 2),
                          c(1, 2),
                          c(1, 3),
                          c(1, 3))
  )

  # Open a new PDF device
  pdf(output_pdf, width = pdf_width, height = pdf_height)

  # Draw the plot
  grid.draw(arranged_plot)

  # Add labels (a, b, c) to the plots
  grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize = 32, fontface = "bold"))
  grid.text("b", x = 0.48, y = 0.98, gp = gpar(fontsize = 32, fontface = "bold"))
  grid.text("c", x = 0.48, y = 0.48, gp = gpar(fontsize = 32, fontface = "bold"))

  # Close the PDF device
  dev.off()

  # Optionally, return a message indicating completion
  message("Figure generation completed, and saved to: ", output_pdf)
}
