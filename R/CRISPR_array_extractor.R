#' CRISPR array extractor
#'
#' @param df A data frame containing the data from the InterProScan (tsv format).
#' @param gb_dir A string specifying the directory where the InterProScan output files are located. If NULL, the current working directory is used.
#' @return A data frame containing the baseline table.
#' @export
#'


CRISPR_array_extractor <- function(gb_dir = NULL, save_output = FALSE) {
  library(readr)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(stringr)
  library(splitstackshape)
  library(openxlsx)

  # If gb_dir is NULL, use the current working directory
  if (is.null(gb_dir)) {
    gb_dir <- getwd()
  }

  # Get the list of .tsv files in the directory
  gb_files <- list.files(
    gb_dir,
    pattern = "\\.gbk$|\\.gb$|\\.gbff$",
    full.names = FALSE
  )
  # Exclude temporary files starting with '~$'
  gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  # Read GenBank file
  origin <- read_fwf(gb_files, readr::fwf_widths(c(12, Inf)), show_col_types = FALSE)

  # Delete nucleotide sequences
  origin <- origin %>% dplyr::slice(., 1:as.numeric(grep("ORIGIN", origin$X1, ignore.case = FALSE)[1]) - 1)
  origin <- origin %>% sapply(function(x) str_squish(x)) %>% as.data.frame()

  # Report table
  report <- origin %>% dplyr::slice(., 1:as.numeric(grep("FEATURES", origin$X1, ignore.case = FALSE)[1]) - 1)
  report <- report %>% sapply(function(x) str_squish(x)) %>% as.data.frame()

  # Feature table
  feature <- origin %>%
    dplyr::slice(as.numeric(grep("FEATURES", origin$X1, ignore.case = FALSE)[1]):nrow(.))

  # Now feature is created, and you can reference it in the next step
  feature <- feature %>%
    dplyr::slice(as.numeric(grep("gene", feature$X1, ignore.case = FALSE)[1]):nrow(.)) %>%
    tidyr::fill(X1)

  # Feature info
  feature_info <- rle(feature$X1)
  feature_info <- data.frame(unclass(feature_info)) %>%
    mutate("end" = cumsum(lengths), "start" = end - lengths + 1) %>%
    dplyr::arrange(values, start, end, lengths)

  # Extract CRISPR features
  CRISPR <- feature_info %>% dplyr::filter(str_detect(values, "repeat"))

  if (nrow(CRISPR) != 0) {
    CRISPR <- CRISPR %>% mutate("discriminator" = 1:nrow(CRISPR))
    temp_rows <- unlist(lapply(1:nrow(CRISPR), function(i) CRISPR$start[i]:CRISPR$end[i]))
    CRISPR_table <- feature[temp_rows,] %>%
      mutate("discriminator" = rep(CRISPR$discriminator, CRISPR$lengths)) %>%
      tidyr::separate(., col = "X2", into = c("category", "value"), sep = "\\=") %>%
      tidyr::separate(., col = "category", into = c("category", "value2"), sep = " ", fill = "right", extra = "merge") %>%
      mutate("category" = gsub("/", "", category)) %>%
      dplyr::select(-X1) %>%
      replace(is.na(.), "") %>%
      reshape2::melt(., id = c("discriminator", "category"), na.rm = TRUE) %>%
      dcast(formula = discriminator ~ category, fun.aggregate = toString) %>%
      mutate_all(~str_replace_all(., " ", "")) %>%
      mutate_all(~str_replace_all(., "^,{1,}|,{1,}$", "")) %>%
      mutate_all(~str_replace_all(., "\"", "")) %>%
      tidyr::separate(region, into = c("repeat_start", "repeat_end"), sep = "\\.\\.")

    # Combine CRISPR tables for each contig
    for (i in 1:nrow(contig_length)) {
      assign(paste("CRISPR_table", i, sep = "_"), CRISPR_table %>% mutate("sequence" = str_sub(get(paste("contig", i, "seq", sep = "_")), repeat_start, repeat_end)))
      assign(paste("CRISPR_table", i, sep = "_"), get(paste("CRISPR_table", i, sep = "_")) %>% mutate("sequence_length" = nchar(sequence)))
    }

    if (exists(paste0("CRISPR_table_", 1))) {
      CRISPR_table <- CRISPR_table_1
      for (i in (1:nrow(contig_length))[-1]) {
        if (exists(paste0("CRISPR_table_", i))) {
          CRISPR_table <- rbind.fill(CRISPR_table, get(paste("CRISPR_table", i, sep = "_")))
        }
      }
    }

    CRISPR_table <- CRISPR_table %>%
      mutate(rpt_unit_seq = toupper(rpt_unit_seq)) %>%
      mutate_at("sequence", ~na_if(., '')) %>%
      dplyr::filter(!is.na(sequence))

    # Process spacers
    CRISPR_spacer <- lapply(1:nrow(CRISPR_table), function(i) strsplit(CRISPR_table$sequence[i], split = CRISPR_table$rpt_unit_seq[i]) %>% setNames(CRISPR_table$rpt_unit_seq[i]))
    CRISPR_by_spacer <- data.frame()

    for (i in 1:length(CRISPR_spacer)) {
      discriminator <- names(CRISPR_spacer[[i]])
      by_spacer <- cbind(discriminator, CRISPR_spacer[[i]] %>% as.data.frame() %>% setNames("spacer"))
      CRISPR_by_spacer <- plyr::rbind.fill(CRISPR_by_spacer, by_spacer)
    }

    CRISPR_by_spacer <- CRISPR_by_spacer %>%
      mutate_all(na_if, "") %>%
      drop_na(spacer) %>%
      mutate("length" = nchar(spacer))

    CRISPR_by_spacer_dist <- CRISPR_by_spacer %>%
      distinct(spacer, .keep_all = TRUE) %>%
      mutate("N20" = str_sub(spacer, start = -20), "length_N20" = nchar(N20))

    # Save results to Excel and FASTA if save_output is TRUE
    if (save_output) {
      write.xlsx(x = CRISPR_by_spacer_dist, file = file.path(gb_dir, "CRISPR_by_spacer.xlsx"), row.names = FALSE)
      CRISPR_by_spacer_dist %>%
        dplyr::select(discriminator, spacer) %>%
        mutate("discriminator" = paste0(">", discriminator)) %>%
        write.table(file.path(gb_dir, "CRISPR_spacer.fasta"), row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)

      CRISPR_by_spacer_dist %>%
        dplyr::select(discriminator, N20) %>%
        mutate("discriminator" = paste0(">", discriminator)) %>%
        write.table(file.path(gb_dir, "CRISPR_spacer_N20.fasta"), row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)
    }

    # Assign the result to 'CRISPR array' in the global environment
    assign("CRISPR_by_spacer", CRISPR_by_spacer_dist, envir = .GlobalEnv)
    message("The result has been saved to the R environment variable 'CRISPR array'.")

  } else {
    message("No CRISPR sequences found.")
  }
}
