#' Result DNMB table
#'
#' @param gb_table A data frame containing the data.
#' @param InterPro_search Logical, whether to include p-values for means.
#' @param InterPro_site Logical, whether to include p-values for medians.
#' @param codon Logical, whether to use chi-square test for categorical variables.
#' @return A data frame containing the baseline table.
#' @export
#'

DNMB_table <- function(
    genbank_table = TRUE,
    EggNOG_table = FALSE,
    InterProScan_table = FALSE,
    InterProScan_site = FALSE,
    codon_usage = TRUE,
    tRNA_anticodon = TRUE,
    tRNA_distribution = TRUE,
    RBS_table = TRUE,
    CRISPR_by_spacer = TRUE,
    save_dir = NULL
) {
  library(openxlsx)

  # Set the save directory to the current working directory if not provided
  if (is.null(save_dir)) {
    save_dir <- getwd()
  }

  # Check for genbank_table
  if (genbank_table == TRUE) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      genbank_table <- get("genbank_table", envir = .GlobalEnv)
    } else {
      stop("genbank_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for InterProScan_table
  if (EggNOG_table == TRUE) {
    if (exists("EggNOG_table", envir = .GlobalEnv)) {
      InterPro_table <- get("EggNOG_table", envir = .GlobalEnv)
    } else {
      stop("EggNOG_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for codon_usage_table
  if (codon_usage == TRUE) {
    if (exists("codon_usage", envir = .GlobalEnv)) {
      codon_usage <- get("codon_usage", envir = .GlobalEnv)
    } else {
      stop("codon_usage is set to TRUE but not found in the environment.")
    }
  }

  # Check for tRNA_anticodon_table
  if (tRNA_anticodon == TRUE) {
    if (exists("tRNA_anticodon", envir = .GlobalEnv)) {
      tRNA_anticodon <- get("tRNA_anticodon", envir = .GlobalEnv)
    } else {
      stop("tRNA_anticodon is set to TRUE but not found in the environment.")
    }
  }

  # Check for tRNA_distribution_table
  if (tRNA_distribution == TRUE) {
    if (exists("tRNA_distribution", envir = .GlobalEnv)) {
      tRNA_distribution <- get("tRNA_distribution", envir = .GlobalEnv)
    } else {
      stop("tRNA_distribution is set to TRUE but not found in the environment.")
    }
  }

  # Check for RBS_table
  if (RBS_table == TRUE) {
    if (exists("RBS_table", envir = .GlobalEnv)) {
      RBS_table <- get("RBS_table", envir = .GlobalEnv)
    } else {
      stop("RBS_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for CRISPR_by_spacer
  if (CRISPR_by_spacer == TRUE) {
    if (exists("CRISPR_by_spacer", envir = .GlobalEnv)) {
      CRISPR_by_spacer <- get("CRISPR_by_spacer", envir = .GlobalEnv)
    } else {
      stop("CRISPR_by_spacer is set to TRUE but not found in the environment.")
    }
  }

  # Check for InterProScan_site
  if (InterProScan_site == TRUE) {
    if (exists("InterProScan_site", envir = .GlobalEnv)) {
      InterProScan_site <- get("InterProScan_site", envir = .GlobalEnv)
      message("InterProScan_site found in the environment.")
    } else {
      message("InterProScan_site is set to TRUE but not found in the environment. Skipping the InterProScan_site sheet.")
      InterProScan_site <- NULL
    }
  }

  # Create a new workbook
  wb <- createWorkbook()

  # Add worksheets
  addWorksheet(wb, "1.GenBank_table")
  addWorksheet(wb, "2.Codon_usage")
  addWorksheet(wb, "3.tRNA_anticodon")
  addWorksheet(wb, "4.tRNA_distribution")
  addWorksheet(wb, "5.RBS")
  addWorksheet(wb, "6.CRISPR_table")

  # Write data to worksheets
  writeData(wb, "1.GenBank_table", genbank_table, startRow = 1, startCol = 1)
  writeData(wb, "2.Codon_usage", codon_usage, startRow = 1, startCol = 1)
  writeData(wb, "3.tRNA_anticodon", tRNA_anticodon, startRow = 1, startCol = 1)
  writeData(wb, "4.tRNA_distribution", tRNA_distribution, startRow = 1, startCol = 1)
  writeData(wb, "5.RBS", RBS_table, startRow = 1, startCol = 1)
  writeData(wb, "6.CRISPR_table", CRISPR_by_spacer, startRow = 1, startCol = 1)

  # Adding the worksheet only if InterProScan_site exists
  if (!is.null(InterProScan_site)) {
    addWorksheet(wb, "7.InterPro_site")
    writeData(wb, "7.InterPro_site", InterProScan_site, startRow = 1, startCol = 1)
    message("InterProScan_site sheet added.")
  } else {
    message("InterProScan_site sheet skipped.")
  }

  # Construct the file name
    gb_dir <- getwd()
    gb_files <- list.files(gb_dir,pattern = "\\.gbk$|\\.gb$|\\.gbff$",full.names = FALSE)
  # Exclude temporary files starting with '~$'
    gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  file_name <- paste0(qdap::beg2char(gb_files, "."), "_total.xlsx")
  save_path <- file.path(save_dir, file_name)

  # Save the workbook
  saveWorkbook(wb, file = save_path, overwrite = TRUE)

  message(paste("Workbook saved as", save_path))
}
