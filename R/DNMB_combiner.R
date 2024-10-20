#' Result DNMB table
#'
#' @param gb_table A data frame containing the data.
#' @param InterPro_search Logical, whether to include p-values for means.
#' @param InterPro_site Logical, whether to include p-values for medians.
#' @param codon Logical, whether to use chi-square test for categorical variables.
#' @return A data frame containing the baseline table.
#' @export
#'

run_DNMB_combiner <- function(genbank_table = NULL, EggNOG_table = NULL, InterProScan_table = NULL) {
  library(openxlsx)

    # Check if genbank_table is NULL and exists in the environment
  if (is.null(genbank_table)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      genbank_table <- get("genbank_table", envir = .GlobalEnv)
      message("genbank_table fetched from the global environment.")
    } else {
      stop("Error: 'genbank_table' not found in the environment.")
    }
  }

  # Step 7: Call EggNOG_annotations and merge with genbank_table
  EggNOG_annotations()
  message("Step7. EggNOG_annotations function executed")

  # Check if EggNOG_table is NULL and exists in the environment
  if (is.null(EggNOG_table)) {
    if (exists("EggNOG_table", envir = .GlobalEnv)) {
      EggNOG_table <- get("EggNOG_table", envir = .GlobalEnv)
      message("EggNOG_table fetched from the global environment.")
    } else {
      stop("Error: 'EggNOG_table' not found in the environment.")
    }
  }

  # Check if the necessary columns exist in the tables before merging
  if (!"locus_tag" %in% colnames(genbank_table)) stop("Error: 'locus_tag' not found in genbank_table")
  if (!"query" %in% colnames(EggNOG_table)) stop("Error: 'query' not found in EggNOG_table")

  genbank_table <- merge(genbank_table, EggNOG_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)


  # Step 8: Call InterProScan_annotations and merge with genbank_table
  InterProScan_annotations()
  message("Step8. InterProScan_annotations function executed")

    # Check if InterProScan_table is NULL and exists in the environment
  if (is.null(InterProScan_table)) {
    if (exists("InterProScan_table", envir = .GlobalEnv)) {
      InterProScan_table <- get("InterProScan_table", envir = .GlobalEnv)
      message("InterProScan_table fetched from the global environment.")
    } else {
      stop("Error: 'InterProScan_table' not found in the environment.")
    }
  }


  # Check if the necessary columns exist in the tables before merging
  if (!"locus_tag" %in% colnames(genbank_table)) stop("Error: 'locus_tag' not found in genbank_table")
  if (!"query" %in% colnames(InterProScan_table)) stop("Error: 'query' not found in InterProScan_table")

  genbank_table <- merge(genbank_table, InterProScan_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)

  # Assign the result to 'genbank_table' in the global environment
  assign("genbank_table", genbank_table, envir = .GlobalEnv)
  message("The result has been saved to the R environment variable 'genbank_table'.")

  # Optionally, return a message indicating completion
  message("All functions executed sequentially.")
}
