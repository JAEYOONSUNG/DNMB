#' Count tRNA anticodon
#'
#' @param df A data frame containing the data from the genbank.
#' @param gb_dir A string specifying the directory where the genbank output files are located. If NULL, the current working directory is used.
#' @return A data frame containing the tRNA coding region and distribution.
#' @export
#'

tRNA_anticodon_counter <- function(target = NULL, save_output = NULL) {
  library(dplyr)
  library(stringr)
  library(seqinr)
  library(Biostrings)

  # If target is NULL, check if 'genbank_table' exists in the global environment
  if (is.null(target)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      target <- get("genbank_table", envir = .GlobalEnv)
      message("Using 'genbank_table' from the global environment as target.")
    } else {
      stop("You need to provide a 'target' or ensure 'genbank_table' exists in the global environment.")
    }
  }

  # Filter for rows with anticodon information
  tRNA_anticodon <- target %>% dplyr::filter(anticodon != "")

  # Extract anticodon location and clean up the column
  tRNA_anticodon <- tRNA_anticodon %>%
    mutate(
      anticodon_location = stringr::str_extract(anticodon, "pos\\:.{1,}[0-9]{1,}\\.\\.[0-9]{1,}") %>% gsub("pos:", "", .),
      anticodon_location = dplyr::case_when(
        grepl("complement", anticodon_location) ~ paste0(anticodon_location, ")"),
        TRUE ~ anticodon_location
      ),
      anticodon = str_extract(anticodon_location, "[0-9]{1,}\\.\\.[0-9]{1,}"),

      # Extract anticodon from the nt_seq and ensure it is converted to character
      anticodon_position = as.numeric(str_extract(anticodon, "^[0-9]+")) - as.numeric(start),
      )

  tRNA_anticodon <- tRNA_anticodon %>%
    dplyr::mutate(
      anticodon = mapply(function(nt_seq, pos, anticodon_location) {
        # Ensure that anticodon_position is valid and within range
        if (!is.na(pos) && pos > 0 && pos <= (nchar(nt_seq) - 2)) {
          anticodon_seq <- substr(nt_seq, pos + 1, pos + 3)

          # Check if 'complement' is in anticodon_location, and reverse complement if needed
          if (grepl("complement", anticodon_location)) {
            return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(anticodon_seq))))
          } else {
            return(anticodon_seq)
          }
        } else {
          return(NA)  # Return NA if position is invalid
        }
      }, nt_seq, anticodon_position, anticodon_location) # Add anticodon_location as a third argument
    )

  # Ensure anticodon is not NA and is a valid DNA string
  tRNA_anticodon <- tRNA_anticodon %>%
    dplyr::filter(!is.na(anticodon) & anticodon != "") %>%  # Filter out NAs or empty anticodon values
    dplyr::mutate(anticodon = as.character(anticodon))  # Ensure anticodon is treated as a character

  tRNA_anticodon <- tRNA_anticodon %>%
    mutate(
      tRNA_codon = dplyr::case_when(
        grepl("^complement\\([0-9]{1,}\\.\\.[0-9]{1,}", anticodon_location) ~ sapply(anticodon, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))),
        grepl("^[0-9]{1,}\\.\\.[0-9]{1,}", anticodon_location) ~ sapply(anticodon, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))
      )
    )

  tRNA_anticodon <- tRNA_anticodon %>%
    mutate(
      anticodon_location = gsub("complement\\(", "", anticodon_location)  # "complement(" trimming
    )

  tRNA_anticodon <- tRNA_anticodon %>% dplyr::select(locus_tag, start, end, direction, product, anticodon, tRNA_codon, nt_seq, anticodon_location, anticodon_position)

  tRNA_anticodon <- tRNA_anticodon %>%
    mutate(
      anticodon = tolower(anticodon),
      tRNA_codon = tolower(tRNA_codon),
      AA = sapply(tRNA_codon, function(codon) {
        amino_acid <- seqinr::translate(s2c(codon))
        return(amino_acid)
      })
    )

  # Create the tRNA distribution table and match it with the codon table
  # tRNA_anticodon: The existing tRNA codon table (the data frame you already have)
  # Group by AA and tRNA_codon to calculate codon count

  tRNA_codon_count <- tRNA_anticodon %>%
    dplyr::group_by(AA, tRNA_codon) %>%
    dplyr::summarise(codon_count = n(), .groups = 'drop')  # Calculate the count of each tRNA codon

  # Load the standard codon table provided by the seqinr package
  codon_table <- as.data.frame(matrix(c(seqinr::words(), seqinr::translate(sapply(seqinr::words(), s2c))),
                                      ncol = 2, byrow = FALSE))
  colnames(codon_table) <- c("codon", "AA")
  codon_table %>% group_by(AA) %>% dplyr::arrange(codon_table, .by_group = FALSE)
  codon_order <- codon_table %>% group_by(AA) %>% top_n(1) %>% dplyr::select(AA) %>% pull() %>% as.vector()

  # Match the codon table with the tRNA distribution
  # The codon_table has codon and amino acid information, match it with the tRNA distribution
  tRNA_distribution <- merge(codon_table, tRNA_codon_count %>% dplyr::select(-AA), by.x = "codon", by.y = "tRNA_codon", all.x = TRUE)
  tRNA_distribution$codon_count[is.na(tRNA_distribution$codon_count)] <- 0  # Set unmatched codon counts to 0
  tRNA_distribution <- tRNA_distribution %>% dplyr::arrange(match(AA, codon_order))
  # Save the output as 'codon_usage.xlsx' if save_output is TRUE
  if (!is.null(save_output) && save_output) {
    output_file <- file.path(getwd(), paste0("tRNA_anticodon",".xlsx"))
    openxlsx::write.xlsx(tRNA_anticodon, output_file, rowNames = FALSE)
    message(paste("Results have been saved to:", output_file))
  }

  # Assign the final result to the environment (optional)
  assign("tRNA_anticodon", tRNA_anticodon, envir = .GlobalEnv)
  assign("tRNA_distribution", tRNA_distribution, envir = .GlobalEnv)
}
