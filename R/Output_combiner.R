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

  # Check for EggNOG_table
  if (EggNOG_table == TRUE) {
    if (exists("EggNOG_table", envir = .GlobalEnv)) {
      InterPro_table <- get("EggNOG_table", envir = .GlobalEnv)
    } else {
      stop("EggNOG_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for codon_usage
  if (codon_usage == TRUE) {
    if (exists("codon_usage", envir = .GlobalEnv)) {
      codon_usage <- get("codon_usage", envir = .GlobalEnv)
    } else {
      stop("codon_usage is set to TRUE but not found in the environment.")
    }
  }

  # Check for tRNA_anticodon
  if (tRNA_anticodon == TRUE) {
    if (exists("tRNA_anticodon", envir = .GlobalEnv)) {
      tRNA_anticodon <- get("tRNA_anticodon", envir = .GlobalEnv)
    } else {
      stop("tRNA_anticodon is set to TRUE but not found in the environment.")
    }
  }

  # Check for tRNA_distribution
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
      warning("CRISPR_by_spacer is set to TRUE but not found in the environment. Creating an empty sheet.")
      CRISPR_by_spacer <- NULL
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
  wb <- openxlsx::createWorkbook()

  # Add worksheets
  openxlsx::addWorksheet(wb, "1.GenBank_table")
  openxlsx::addWorksheet(wb, "2.Codon_usage")
  openxlsx::addWorksheet(wb, "3.tRNA_anticodon")
  openxlsx::addWorksheet(wb, "4.tRNA_distribution")
  openxlsx::addWorksheet(wb, "5.RBS")
  openxlsx::addWorksheet(wb, "6.CRISPR_table")

  no_wrap_style <- openxlsx::createStyle(wrapText = FALSE)

  # Write data to worksheets
  openxlsx::writeData(wb, "1.GenBank_table", genbank_table, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "2.Codon_usage", codon_usage, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "3.tRNA_anticodon", tRNA_anticodon, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "4.tRNA_distribution", tRNA_distribution, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "5.RBS", RBS_table, startRow = 1, startCol = 1)

  # Handle CRISPR table (empty or not)
  if (is.null(CRISPR_by_spacer)) {
    openxlsx::writeData(wb, "6.CRISPR_table", data.frame(), startRow = 1, startCol = 1)
  } else {
    openxlsx::writeData(wb, "6.CRISPR_table", CRISPR_by_spacer, startRow = 1, startCol = 1)
  }

  openxlsx::addStyle(wb, sheet = "1.GenBank_table", style = no_wrap_style, rows = 1:nrow(genbank_table) + 1, cols = 1:ncol(genbank_table), gridExpand = TRUE)

  # Adding the worksheet only if InterProScan_site exists
  if (!is.null(InterProScan_site)) {
    openxlsx::addWorksheet(wb, "7.InterPro_site")
    openxlsx::writeData(wb, "7.InterPro_site", InterProScan_site, startRow = 1, startCol = 1)
    message("InterProScan_site sheet added.")
  }

  # Construct the file name
  gb_dir <- getwd()
  gb_files <- list.files(gb_dir, pattern = "\\.gbk$|\\.gb$|\\.gbff$", full.names = FALSE)
  gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  file_name <- paste0(qdap::beg2char(gb_files, "."), "_total.xlsx")
  save_path <- file.path(save_dir, file_name)

  # Save the workbook
  openxlsx::saveWorkbook(wb, file = save_path, overwrite = TRUE)

  message(paste("Workbook saved as", save_path))
}
