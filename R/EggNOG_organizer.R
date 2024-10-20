#' EggNOG emapper table #
#'
#' @param df A data frame containing the data from eggNOG-mapper (http://eggnog-mapper.embl.de).
#' @param EggNOG_dir A string specifying the directory where the eggNOG-mapper output files are located. If NULL, the current working directory is used.
#' @param save_output A logical value indicating whether to save the output to an Excel file in XLSX format. Defaults to FALSE.
#' @return A list containing summary information of the processed data.
#' @export
#'

EggNOG_annotations <- function(EggNOG_dir = NULL, save_output = FALSE) {
  library(dplyr)
  library(stringr)
  library(openxlsx)

  # If EggNOG_dir is NULL, use the current working directory
  if (is.null(EggNOG_dir)) {
    EggNOG_dir <- getwd()
  }

  # Get the list of .emapper.annotations.xlsx files in the directory
  EggNOG_file <- list.files(EggNOG_dir, pattern = "\\.emapper.annotations\\.xlsx$")
  EggNOG_file <- EggNOG_file[!grepl("^~\\$", EggNOG_file)]

  # Read each file and convert it into a data frame
  eggNOG_list <- lapply(EggNOG_file, function(file) {
    Egg_temp <- openxlsx::read.xlsx(
      file.path(EggNOG_dir, file),
      startRow = 3,    # Adjust startRow to correctly read column names
      colNames = TRUE  # Ensure column names are read
    )
    # Remove the last 3 rows if necessary
    if (nrow(Egg_temp) > 3) {
      Egg_temp <- Egg_temp[1:(nrow(Egg_temp) - 3), ]
    }
    Egg_temp
  })

  # Combine all data frames into one
  eggNOG <- bind_rows(eggNOG_list)

  # Check if 'query' column exists
  if ("query" %in% colnames(eggNOG)) {
    # Clean up the 'query' column
    eggNOG <- eggNOG %>%
      mutate(
        query = if_else(str_detect(query, "gnl\\|"), gsub("gnl\\|", "", query), query),
        query = if_else(str_detect(query, "\\|"), gsub("\\|", ":", query), query),
        query = if_else(str_detect(query, "lcl:(AC|NC|BG|NT|NW|NZ)_([a-zA-Z]+)?[0-9]+\\.[0-9]+_prot_"), gsub("lcl:(AC|NC|BG|NT|NW|NZ)_([a-zA-Z]+)?[0-9]+\\.[0-9]+_prot_", "", query), query),
        query = if_else(str_detect(query, "_[0-9]{1,4}$"), gsub("_[0-9]{1,4}$", "", query), query),
        query = if_else(str_detect(query, "^extdb:"), gsub("^extdb:", "", query), query)
      ) %>%
      filter(!str_detect(query, "^[0-9]+$"))
  } else {
    # If 'query' column doesn't exist, proceed without cleanup
    message("The 'query' column was not found in the eggNOG data frame. Proceeding without 'query' cleanup.")
  }

  # Create 'COG_category_for_plot' column and assign color codes and legends
  eggNOG <- eggNOG %>%
    mutate(COG_category_for_plot = str_sub(COG_category, 1, 1)) %>%
    mutate(COG_color = case_when(
      COG_category_for_plot == "J" ~ '#ff0000',
      COG_category_for_plot == "A" ~ '#c2af58',
      COG_category_for_plot == "K" ~ '#ff9900',
      COG_category_for_plot == "L" ~ '#ffff00',
      COG_category_for_plot == "B" ~ '#ffc600',
      COG_category_for_plot == "D" ~ '#99ff00',
      COG_category_for_plot == "Y" ~ '#493126',
      COG_category_for_plot == "V" ~ '#ff008a',
      COG_category_for_plot == "T" ~ '#0000ff',
      COG_category_for_plot == "M" ~ '#9ec928',
      COG_category_for_plot == "N" ~ '#006633',
      COG_category_for_plot == "Z" ~ '#660099',
      COG_category_for_plot == "W" ~ '#336699',
      COG_category_for_plot == "U" ~ '#33cc99',
      COG_category_for_plot == "O" ~ '#00ffff',
      COG_category_for_plot == "C" ~ '#9900ff',
      COG_category_for_plot == "G" ~ '#805642',
      COG_category_for_plot == "E" ~ '#ff00ff',
      COG_category_for_plot == "F" ~ '#99334d',
      COG_category_for_plot == "H" ~ '#727dcc',
      COG_category_for_plot == "I" ~ '#5c5a1b',
      COG_category_for_plot == "P" ~ '#0099ff',
      COG_category_for_plot == "Q" ~ '#ffcc99',
      COG_category_for_plot == "R" ~ '#ff9999',
      COG_category_for_plot == "S" ~ '#d6aadf',
      TRUE ~ NA_character_
    )) %>%
    mutate(COG_legend = case_when(
      COG_category_for_plot == "J" ~ '[J] Translation, ribosomal structure and biogenesis',
      COG_category_for_plot == "A" ~ '[A] RNA processing and modification',
      COG_category_for_plot == "K" ~ '[K] Transcription',
      COG_category_for_plot == "L" ~ '[L] Replication, recombination and repair',
      COG_category_for_plot == "B" ~ '[B] Chromatin structure and dynamics',
      COG_category_for_plot == "D" ~ '[D] Cell cycle control, cell division, chromosome partitioning',
      COG_category_for_plot == "Y" ~ '[Y] Nuclear structure',
      COG_category_for_plot == "V" ~ '[V] Defense mechanisms',
      COG_category_for_plot == "T" ~ '[T] Signal transduction mechanisms',
      COG_category_for_plot == "M" ~ '[M] Cell wall/membrane/envelope biogenesis',
      COG_category_for_plot == "N" ~ '[N] Cell motility',
      COG_category_for_plot == "Z" ~ '[Z] Cytoskeleton',
      COG_category_for_plot == "W" ~ '[W] Extracellular structures',
      COG_category_for_plot == "U" ~ '[U] Intracellular trafficking, secretion, and vesicular transport',
      COG_category_for_plot == "O" ~ '[O] Posttranslational modification, protein turnover, chaperones',
      COG_category_for_plot == "C" ~ '[C] Energy production and conversion',
      COG_category_for_plot == "G" ~ '[G] Carbohydrate transport and metabolism',
      COG_category_for_plot == "E" ~ '[E] Amino acid transport and metabolism',
      COG_category_for_plot == "F" ~ '[F] Nucleotide transport and metabolism',
      COG_category_for_plot == "H" ~ '[H] Coenzyme transport and metabolism',
      COG_category_for_plot == "I" ~ '[I] Lipid transport and metabolism',
      COG_category_for_plot == "P" ~ '[P] Inorganic ion transport and metabolism',
      COG_category_for_plot == "Q" ~ '[Q] Secondary metabolites biosynthesis, transport and catabolism',
      COG_category_for_plot == "R" ~ '[R] General function prediction only',
      COG_category_for_plot == "S" ~ '[S] Function unknown',
      TRUE ~ NA_character_
    ))

  # Save the output as 'EggNOG_table.xlsx' if save_output is TRUE
  if (save_output) {
    output_file <- file.path(EggNOG_dir, "EggNOG_table.xlsx")
    openxlsx::write.xlsx(eggNOG, output_file, rowNames = FALSE)
    message(paste("Results have been saved to:", output_file))
  }

  # Assign the result to 'EggNOG_table' in the global environment
  assign("EggNOG_table", eggNOG, envir = .GlobalEnv)
  message("The result has been saved to the R environment variable 'EggNOG_table'.")

  # Return the eggNOG data frame
  return(tibble::tibble(eggNOG))
}
