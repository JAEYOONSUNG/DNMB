#' genbank feature table
#'
#' This is a simple example function.
#'
#' @param gb_dir A genbank file containing nucleotide such as ".gb", "gbk", ".gbff". the directory where the genbank files are located. If NULL, the current working directory is used.
#' @return A data frame from genbank file.
#' @export
#'


gb_info <- function(gb_dir = NULL) {
  # Load necessary libraries
  library(dplyr)
  library(readr)
  library(tidyr)
  library(openxlsx)
  library(plyr)

  # Initialize empty data frame for storing genome table
  Genome_summary <- data.frame()

  # Read the specified GenBank file

  # If InterProScan_dir is NULL, use the current working directory
  if (is.null(gb_dir)) {
    gb_dir <- getwd()
  }

  # Get the list of .tsv files in the directory
  gb_files <- list.files(
    gb_dir,
    pattern = "\\.gbk|gb|gbff$",
    full.names = FALSE
  )
  # Exclude temporary files starting with '~$'
  gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  # Check if InterPro search files exist
  if (length(gb_files) == 0) {
    stop("No genbank files found.")
  }

  origin <- readr::read_fwf(gb_files, fwf_widths(c(10, Inf)))

  # Extract portion before "FEATURES"
  origin <- origin %>%
    slice(., 1:as.numeric(grep("FEATURES", origin$X1, ignore.case = FALSE)[1]) - 1) %>%
    sapply(function(x) str_squish(x)) %>%
    as.data.frame()

  # Process genome info table
  Genome_info_table <- origin %>%
    slice(., 1:as.numeric(grep("##", origin$X2, ignore.case = FALSE)[1]) - 1) %>%
    mutate_all(~replace(., is.na(.), NA)) %>%
    tidyr::fill(X1) %>%
    group_by(X1) %>%
    dplyr::summarise(X2 = paste(X2, collapse = " ")) %>%
    mutate(X2 = ifelse(X2 == "NA", "", X2))

  # Process genome annotation table
  Genome_annotation_table <- origin %>%
    slice(as.numeric(grep("^##Genome-Annotation-Data-START##", .$X2, ignore.case = FALSE)[1] + 1):
            as.numeric(grep("^##Genome-Annotation-Data-END##", .$X2, ignore.case = FALSE)[1] - 1)) %>%
    separate(col = "X2", into = c("X1", "X2"), sep = " :: ", fill = "left", extra = "merge") %>%
    mutate_all(~replace(., is.na(.), NA)) %>%
    tidyr::fill(X1)   %>%
    group_by(X1) %>%
    dplyr::summarise(X2 = paste(X2, collapse = " "))

  # Combine info and annotation tables
  Genome_summary <- rbind(Genome_info_table, Genome_annotation_table) %>% as.data.frame()

  # Remove empty first row if necessary
  Genome_summary <- Genome_summary %>%
    slice(if (is.na(Genome_summary[1, 1]) | Genome_summary[1, 1] == "") 2:n() else 1:n())

  # Final processing on the genome table
  rownames(Genome_summary) <- Genome_summary$X1
  Genome_summary <- Genome_summary %>% dplyr::select(-c(X1))
  Genome_summary <- rbind("File" = gsub("\\.gb(k|ff)?$", "", basename(gb_files)), Genome_summary)
  Genome_summary <- t(Genome_summary) %>% as.data.frame()
  rownames(Genome_summary) <- NULL

  # Assign the table to the global environment (optional)
  assign("Genome_summary", Genome_summary, envir = .GlobalEnv)
}
