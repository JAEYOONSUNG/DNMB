#' Result DNMB table
#'
#' @param gb_table A data frame containing the data.
#' @param InterPro_search Logical, whether to include p-values for means.
#' @param InterPro_site Logical, whether to include p-values for medians.
#' @param codon Logical, whether to use chi-square test for categorical variables.
#' @return A data frame containing the baseline table.
#' @export
#'

run_DNMB <- function() {
  # Call Function Genbank_organizer
  Genbank_organizer()
  message("Step1. genbank_organizer function executed")

  # Call Function codon_usage_calculator
  Codon_usage_calculator()
  message("Step2. Codon usage caculator function executed")

  # Call Count_tRNA_anticodon Function
  tRNA_anticodon_counter()
  message("Step3. tRNA anticodon count function executed")

  # Call Genbank fna extractor Function
  genbank_fna_extractor()
  message("Step4. Genbank fna extractor function executed")

  # Call RBS extractor Function
  RBS_extractor(plot_RBS = TRUE)
  message("Step5. RBS extractor function executed")

  # Call RBS extractor Function
  CRISPR_array_extractor()
  message("Step6. CRISPR array extractor function executed")

  # Call DNMB tRNA_RBS plot
  codon_usage_tRNA_heatmap_generator()
  Figure_generator()
  message("Step7. plot function executed")

  # Step 7: Call EggNOG_annotations and merge with genbank_table
  # Check if 'emapper.annotations' exists in the current working directory
  if (any(grepl("\\.emapper\\.annotations\\.", list.files(getwd())))) {
    # Step 7: Call EggNOG_annotations and merge with genbank_table
    EggNOG_annotations()
    message("Step8. EggNOG_annotations function executed")

    # Check if EggNOG_table is NULL and exists in the environment
    if (is.null(EggNOG_table)) {
      if (exists("EggNOG_table", envir = .GlobalEnv)) {
        EggNOG_table <- get("EggNOG_table", envir = .GlobalEnv)
        message("EggNOG_table fetched from the global environment.")
      } else {
        message("EggNOG_table not found in the environment. Skipping this step.")
        return(invisible(NULL))  # Skip this step and continue
      }
    }

    # Check if the necessary columns exist in the tables before merging
    if (!"locus_tag" %in% colnames(genbank_table)) {
      message("Error: 'locus_tag' not found in genbank_table. Skipping the merge step.")
      return(invisible(NULL))
    }
    if (!"query" %in% colnames(EggNOG_table)) {
      message("Error: 'query' not found in EggNOG_table. Skipping the merge step.")
      return(invisible(NULL))
    }

    # Merge the genbank_table with EggNOG_table
    genbank_table <- merge(genbank_table, EggNOG_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)

    # Assign the result to 'genbank_table' in the global environment
    assign("genbank_table", genbank_table, envir = .GlobalEnv)
    message("The result has been saved to the R environment variable 'genbank_table'.")
  } else {
    message("No 'emapper.annotations' file found in the working directory. Skipping EggNOG_annotations.")
  }

  # Step 8: Call InterProScan_annotations and merge with genbank_table
  # Check if '.tsv' and '.tsv.sites' files exist in the current working directory
  if (any(grepl("\\.tsv$", list.files(getwd()))) || any(grepl("\\.tsv.sites$", list.files(getwd())))) {

    InterProScan_annotations()
    message("Step9. InterProScan_annotations function executed")

    # Check if InterProScan_table is NULL and exists in the environment
    if (is.null(InterProScan_table)) {
      if (exists("InterProScan_table", envir = .GlobalEnv)) {
        InterProScan_table <- get("InterProScan_table", envir = .GlobalEnv)
        message("InterProScan_table fetched from the global environment.")
      } else {
        message("InterProScan_table not found in the environment. Skipping this step.")
        return(invisible(NULL))  # Skip this step and continue
      }
    }

    # Check if the necessary columns exist in the tables before merging
    if (!"locus_tag" %in% colnames(genbank_table)) {
      message("Error: 'locus_tag' not found in genbank_table. Skipping the merge step.")
      return(invisible(NULL))
    }
    if (!"query" %in% colnames(InterProScan_table)) {
      message("Error: 'query' not found in InterProScan_table. Skipping the merge step.")
      return(invisible(NULL))
    }

    # Merge the genbank_table with InterProScan_table
    genbank_table <- merge(genbank_table, InterProScan_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)

    # Assign the result to 'genbank_table' in the global environment
    assign("genbank_table", genbank_table, envir = .GlobalEnv)
    message("The result has been saved to the R environment variable 'genbank_table'.")
  } else {
    message("No '.tsv' or '.tsv.sites' files found in the working directory. Skipping InterProScan_annotations.")
  }

  DNMB_table()
  message("Analysis end: DNMB function executed")
}
