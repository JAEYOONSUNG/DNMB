#' codon usage table
#'
#' @param df A data frame containing the data from the genbank_organizer function in DNMB package.
#' @return A data frame containing the codon usage information.
#' @export
#'

Codon_usage_calculator <- function(target = NULL, save_output = NULL) {
  library(dplyr)
  library(Biostrings)
  library(seqinr)
  library(stringr)

  # If target is NULL, check if 'genbank_table' exists in the global environment
  if (is.null(target)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      target <- get("genbank_table", envir = .GlobalEnv)
      message("Using 'genbank_table' from the global environment as target.")
    } else {
      stop("You need to provide a 'target' or ensure 'genbank_table' exists in the global environment.")
    }
  }


  # Convert factor columns to character if needed
  target <- target %>% dplyr::mutate_if(is.factor, as.character)

  # Filter rows that have a valid protein_id
  codon_CDS <- target %>% dplyr::filter(!protein_id == "")

  # Initialize codon count
  codon_count <- seqinr::uco(s2c((codon_CDS$rearranged_nt_seq)[1]), frame = 0, index = c("eff"), as.data.frame = TRUE, NA.rscu = NA)
  codon <- codon_count[,1:2]  # Extract codon and amino acid columns

  # Loop over all sequences and accumulate codon usage
  for(i in 1:nrow(codon_CDS)) {
    codon_count <- codon_count + seqinr::uco(s2c((codon_CDS$rearranged_nt_seq)[i]), frame = 0, index = c("eff"), as.data.frame = TRUE, NA.rscu = NA)
  }

  # Add frequency to codon table
  codon <- cbind(codon, eff = codon_count[, 3])
  codon <- codon %>% mutate(freq = eff / sum(eff))

  # Calculate RSCU (Relative Synonymous Codon Usage)
  codon_rscu <- codon %>%
    group_by(AA) %>%
    reframe(codon = codon, RSCU = n() * eff / sum(eff)) %>%
    ungroup() %>%
    dplyr::arrange(AA, codon) %>%
    select(-c(AA))  # Remove AA column to prevent duplicates

  # Merge RSCU back into codon table
  codon <- merge(x = codon, y = codon_rscu, by = "codon", all.x = TRUE, sort = FALSE) %>%
    dplyr::arrange(AA, codon) %>%
    dplyr::mutate(codon = tolower(codon))  # Convert codons to uppercase

  # Sort codons by order of standard genetic code
  codon <- codon %>% arrange(match(codon, toupper(seqinr::words())))

  # Group and order by amino acid
  codon_table <- as.data.frame(matrix(c(seqinr::words(), seqinr::translate(sapply(seqinr::words(), s2c))),
                                      ncol = 2, byrow = FALSE))
  colnames(codon_table) <- c("codon", "AA")
  codon_table %>% group_by(AA) %>% dplyr::arrange(codon_table, .by_group = FALSE)
  codon_order <- codon_table %>% group_by(AA) %>% top_n(1) %>% dplyr::select(AA) %>% pull() %>% as.vector()
  codon <- codon %>% dplyr::rename(aaa = AA)
  codon <- codon %>% mutate("AA" = seqinr::a(aaa))
  codon <- merge(codon_table, codon %>% select(-AA), by.x = "codon", by.y = "codon", all.x = TRUE)
  codon <- codon %>% dplyr::arrange(match(AA, codon_order))

  # Save the output as 'codon_usage.xlsx' if save_output is TRUE
  if (!is.null(save_output) && save_output) {
    output_file <- file.path(getwd(), paste0("codon_usage",".xlsx"))
    openxlsx::write.xlsx(codon, output_file, rowNames = FALSE)
    message(paste("Results have been saved to:", output_file))
  }

  # Assign the table to the global environment (optional)
  assign("codon_usage", codon, envir = .GlobalEnv)

  # Return the summary information
  return(tibble::tibble(codon))
}
