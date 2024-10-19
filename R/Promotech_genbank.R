
#' Promotech to sanpgene #
#'
#' @param df A data frame containing the data from eggNOG-mapper (http://eggnog-mapper.embl.de).
#' @param EggNOG_dir A string specifying the directory where the eggNOG-mapper output files are located. If NULL, the current working directory is used.
#' @param save_output A logical value indicating whether to save the output to an Excel file in XLSX format. Defaults to FALSE.
#' @return A list containing summary information of the processed data.
#' @export
#'

# attach promotech promoter prediction result to genbank for the snapgene visualization
# processing my promotech_output
Promoter_prediction_to_genbank <- function(xlsx_dir = getwd(), output_file = "promoter_feature_for_gb") {
  library(dplyr)
  library(stringr)

  # Get the list of genome prediction files (csv format)
  xlsx_file <- list.files(xlsx_dir, pattern = "genome_predictions\\.csv$", full.names = TRUE)

  if (length(xlsx_file) == 0) {
    stop("No genome_predictions.csv file found in the directory.")
  }

  # Read the genome predictions CSV file
  genome_predictions <- read.csv(xlsx_file, sep = "\t")

  # Add promoter position format for SnapGene (depending on strand direction)
  genome_predictions <- genome_predictions %>%
    mutate(position = case_when(
      strand == "-" ~ paste0("complement(", start, "..", end, ")"),
      strand == "+" ~ paste0(start, "..", end)
    ))

  # Generate the SnapGene feature format for promoters
  Feature <- paste0(
    str_pad("Promoter", 21, side = "right"), genome_predictions$position, "\n",
    str_pad("", 21, side = "right"), paste0("/gene=", genome_predictions$score)
  )

  # Write the feature to a file
  write.table(Feature, file = output_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

  message(paste("Promoter features written to", output_file))
}
# for snapgene format
#primer_bind     complement(3617..3634)
#/note="M13 rev"
#/note="common sequencing primer, one of multiple similar variants"
#/note="color: #a020f0; direction: LEFT"
#misc_feature    2591..2806
#/note="pRplsWT"
#/note="color: #ffcc99; direction: RIGHT"

# for primers
#primer_bind     2611..2627
#/label=James
#/note="'hello'"
#/note="color: black; sequence: gtttttgcgccgcccgg"

#primer_bind     complement(2611..2627)
#/label=James
#/note="'hello'"
#/note="color: black; sequence: gtttttgcgccgcccgg"

