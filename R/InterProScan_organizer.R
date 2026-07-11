#' InterProScan table
#'
#' This is a simple example function.
#'
#' @param df A data frame containing the data from the InterProScan (tsv format).
#' @param InterProScan_dir A string specifying the directory where the InterProScan output files are located. If NULL, the current working directory is used.
#' @return A data frame containing the baseline table.
#' @export
#'

# Load InterProScan
InterProScan_annotations <- function(InterProScan_dir = NULL) {
  # Load necessary libraries
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))

  # If InterProScan_dir is NULL, use the current working directory
  if (is.null(InterProScan_dir)) {
    InterProScan_dir <- getwd()
  }

  # Get the list of .tsv files in the directory
  InterPro_search_files <- list.files(
    InterProScan_dir,
    pattern = "\\.tsv$",
    full.names = FALSE
  )
  # Exclude temporary files starting with '~$'
  InterPro_search_files <- InterPro_search_files[!grepl("^~\\$", InterPro_search_files)]

  # Check if InterPro search files exist
  if (length(InterPro_search_files) == 0) {
    stop("No InterProScan files found.")
  }

  # Read each InterPro search file and convert to data frame
  InterPro_search_list <- lapply(InterPro_search_files, function(file) {
    read.csv(
      file.path(InterProScan_dir, file),
      header = FALSE,
      quote = "",
      sep = "\t",
      stringsAsFactors = FALSE
    )
  })
  InterPro_search <- bind_rows(InterPro_search_list)

  # InterProScan emits 13 base columns and appends GO/pathway columns when
  # --goterms and --pathways are enabled. Preserve those optional values when
  # present instead of replacing V14/V15 during standardization.
  if (ncol(InterPro_search) == 13L) {
    InterPro_search$V14 <- "-"
    InterPro_search$V15 <- "-"
  } else if (ncol(InterPro_search) == 14L) {
    InterPro_search$V15 <- "-"
  } else if (ncol(InterPro_search) != 15L) {
    stop("Unexpected number of columns in InterPro search files.")
  }
  InterPro_search$V14[is.na(InterPro_search$V14) | !nzchar(InterPro_search$V14)] <- "-"
  InterPro_search$V15[is.na(InterPro_search$V15) | !nzchar(InterPro_search$V15)] <- "-"
  InterPro_search <- InterPro_search[!is.na(InterPro_search$V3), , drop = FALSE]

  colnames(InterPro_search) <- c(
    "query",
    "Sequence MD5 digest",
    "Sequence length",
    "Analysis",
    "Signature accession",
    "Signature description",
    "Start location",
    "Stop location",
    "Score",
    "Status",
    "Date",
    "InterPro annotations - accession",
    "InterPro annotations - description",
    "GO annotations",
    "Pathways annotations"
  )

  # Load InterProScan_Site
  # Get the list of .tsv.sites files in the directory
  InterPro_site_files <- list.files(
    InterProScan_dir,
    pattern = "\\.tsv\\.sites$",
    full.names = FALSE
  )
  # Exclude temporary files starting with '~$'
  InterPro_site_files <- InterPro_site_files[!grepl("^~\\$", InterPro_site_files)]

  # Check if InterPro site files exist and read them
  if (length(InterPro_site_files) > 0) {
    InterPro_site_list <- lapply(InterPro_site_files, function(file) {
      read.csv(
        file.path(InterProScan_dir, file),
        header = FALSE,
        quote = "",
        sep = "\t",
        stringsAsFactors = FALSE
      )
    })
    InterPro_site <- bind_rows(InterPro_site_list)

    # Assign column names to InterPro site data
    colnames(InterPro_site) <- c(
      "query",
      "Sequence MD5 digest",
      "Sequence length",
      "Identifier",
      "Signature ac",
      "rpsblast-site start",
      "rpsblast-site end",
      "numLocations",
      "site-location residue",
      "site-location residue start",
      "site-location residue end",
      "rpsblast-site description"
    )
  } else {
    InterPro_site <- NULL
  }

  # Expose the parsed site-level table to callers (Output_combiner writes
  # the `7.InterPro_site` sheet from .GlobalEnv$InterProScan_site). Without
  # this assign() the parsed .tsv.sites data would be dropped at function
  # return, leaving the xlsx sheet empty even when a non-empty sites file
  # exists on disk.
  if (exists("InterProScan_site", envir = .GlobalEnv, inherits = FALSE)) {
    rm("InterProScan_site", envir = .GlobalEnv)
  }
  if (!is.null(InterPro_site) && is.data.frame(InterPro_site) && nrow(InterPro_site) > 0) {
    assign("InterProScan_site", InterPro_site, envir = .GlobalEnv)
  }

  # Process InterPro search data
  InterPro <- reshape2::melt(
    InterPro_search,
    id.vars = c("query", "Analysis"),
    measure.vars = c(
      "Signature accession",
      "Signature description",
      "Score",
      "GO annotations",
      "Pathways annotations"
    )
  )

  InterPro <- reshape2::dcast(
    InterPro,
    formula = query ~ variable + Analysis,
    fun.aggregate = toString
  )

  # Reorder columns: make 'query' the first column
  InterPro <- data.table::setcolorder(InterPro, order(sub('.*_', '', names(InterPro))))
  InterPro <- arrange.vars(InterPro, c("query"=1))

  # Clean up the 'query' column
  InterPro$query <- InterPro$query %>%
    str_replace_all("gnl\\|", "") %>%
    str_replace_all("\\|", ":") %>%
    str_replace_all("extdb:", "")

  # Assign the processed data frame to 'InterPro_table' in the R environment
  assign("InterProScan_table", InterPro, envir = .GlobalEnv)
  message("The processed data frame has been saved to the R environment variable 'InterProScan_table'.")

  # Return the summary information
  return(tibble::tibble(InterPro))
}
