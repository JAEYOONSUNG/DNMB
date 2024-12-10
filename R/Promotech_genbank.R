#' Promotech to sanpgene #
#'
#' @param df A data frame containing the data from eggNOG-mapper (http://eggnog-mapper.embl.de).
#' @param EggNOG_dir A string specifying the directory where the eggNOG-mapper output files are located. If NULL, the current working directory is used.
#' @param save_output A logical value indicating whether to save the output to an Excel file in XLSX format. Defaults to FALSE.
#' @return A list containing summary information of the processed data.
#' @export
#'

Promoter_to_genbank <- function(promotech = NULL, genbank = NULL) {
  # Load necessary packages
  library(dplyr)
  library(stringr)
  library(tools)  # For file path operations

  # Get the current working directory
  current_dir <- getwd()

  # Check and search for files ending with promotech.csv
  if (is.null(promotech)) {
    promotech_files <- list.files(current_dir, pattern = "promotech\\.csv$", full.names = TRUE)
    if (length(promotech_files) == 0) {
      stop("No files ending with promotech.csv found in the current directory.")
    } else if (length(promotech_files) > 1) {
      stop("Multiple files ending with promotech.csv found. Please specify the file using the promotech parameter.")
    } else {
      promotech <- promotech_files[1]
    }
  } else {
    if (!file.exists(promotech)) {
      stop(paste("The specified promotech file does not exist:", promotech))
    }
  }

  # Check and search for GenBank files
  if (is.null(genbank)) {
    gb_files <- list.files(current_dir, pattern = "\\.(gb|gbk|gbff)$", full.names = TRUE)
    if (length(gb_files) == 0) {
      stop("No GenBank files ending with .gb, .gbk, or .gbff found in the current directory.")
    } else if (length(gb_files) > 1) {
      stop("Multiple GenBank files found. Please specify the file using the genbank parameter.")
    } else {
      genbank <- gb_files[1]
    }
  } else {
    if (!file.exists(genbank)) {
      stop(paste("The specified GenBank file does not exist:", genbank))
    }
  }

  # Read the promotech.csv file
  promotech_data <- read.csv(promotech, sep = "\t")

  # Read the GenBank file
  gb_content <- readLines(genbank)

  # Find the FEATURES section
  features_start <- which(grepl("^FEATURES", gb_content))
  if (length(features_start) == 0) {
    stop("Cannot find the FEATURES section in the GenBank file.")
  }

  # Find the start position of the ORIGIN section
  origin_start <- which(grepl("^(ORIGIN|BASE COUNT|CONTIG)", gb_content))
  if (length(origin_start) == 0) {
    origin_start <- length(gb_content) + 1
  }

  # Extract existing FEATURES
  existing_features <- gb_content[(features_start + 1):(origin_start - 1)]

  # Prepare promoter features with score in the /label qualifier
  promotech_data <- promotech_data %>%
    mutate(position = case_when(
      strand == "-" ~ paste0("complement(", start, "..", end, ")"),
      strand == "+" ~ paste0(start, "..", end)
    )) %>%
    mutate(feature = paste0(
      "     promoter        ", position, "\n",
      "                     /label=\"Promoter (score: ", score, ")\""
    ))

  promoter_features <- promotech_data$feature

  # Combine existing FEATURES with new promoter FEATURES
  new_features <- c(existing_features, promoter_features)

  # Reconstruct the modified GenBank file content
  new_gb_content <- c(
    gb_content[1:features_start],
    new_features,
    gb_content[origin_start:length(gb_content)]
  )

  # Create output file name by adding '_promotech' before the extension
  base_name <- tools::file_path_sans_ext(genbank)
  extension <- tools::file_ext(genbank)
  output_file <- paste0(base_name, "_promotech.", extension)

  # Save the modified GenBank content to the output file
  writeLines(new_gb_content, output_file)

  message(paste("Modified GenBank file has been saved as", output_file))
}
