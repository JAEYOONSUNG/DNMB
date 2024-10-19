#' Extract the RBS (Ribosome Binding Site) region.
#'
#' This function extracts the Ribosome Binding Site (RBS) region from the target data frame.
#'
#' @param target A data frame containing the target data for analysis.
#' @param dependent A string specifying the dependent variable for the analysis.
#' @param plot_RBS A logical value indicating whether to generate a plot of the RBS region.
#' @param plot_path A string specifying the file path to save the RBS plot (if `plot_RBS` is TRUE).
#' @return A data frame containing the extracted RBS regions.
#' @export

RBS_extractor <- function(target = NULL, plot_RBS = FALSE, save_plot = FALSE, plot_path = NULL) {
  library(dplyr)
  library(Biostrings)
  library(ggplot2)
  library(ggseqlogo)
  library(gridExtra)

  # If target is NULL, check if 'genbank_table' exists in the global environment
  if (is.null(target)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      target <- get("genbank_table", envir = .GlobalEnv)
      message("Using 'genbank_table' from the global environment as target.")
    } else {
      stop("You need to provide a 'target' or ensure 'genbank_table' exists in the global environment.")
    }
  }

  # Extract RBS sequences based on 16S ribosomal RNA
  contig_names <- ls(envir = .GlobalEnv, pattern = "^contig_[0-9]{1,}_seq$")

  # 리스트의 모든 데이터를 rbind로 결합
  contig_list <- lapply(contig_names, function(x) get(x, envir = .GlobalEnv))
  contig_list_rc <- lapply(contig_names, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(get(x, envir = .GlobalEnv)))))

  RBS <- genbank_table %>%
    filter(product == "16S ribosomal RNA") %>%
    select(rearranged_nt_seq) %>%
    mutate(RBS = str_sub(rearranged_nt_seq, start = -9))

  # Reverse complement for RBS sequences
  RBS$RBS <- sapply(RBS$RBS, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))

  # Find the most frequent RBS sequence
  most_frequent_RBS <- toupper(max(RBS$RBS))

  # Search for RBS on the top strand
  strand_plus_RBS <- vmatchPattern(
    most_frequent_RBS,
    toupper(contig_list),
    max.mismatch = 2,
    min.mismatch = 0,
    with.indels = FALSE,
    fixed = TRUE
  ) %>% as.data.frame() %>%
    mutate(strand = "+")

  # Process matches and extract regions
  temp_df <- data.frame()
  for (i in 1:length(contig_list)) {
    temp_RBS_df <- strand_plus_RBS %>%
      filter(group == i) %>%
      mutate(RBS = str_sub(contig_list[i], start, end),
             TR = str_sub(contig_list[i], start, end + 100))
    temp_df <- rbind.fill(temp_df, temp_RBS_df)
  }
  strand_plus_RBS <- temp_df

  # Set range and extract RBS hits
  range <- "{1,10}"
  strand_plus_RBS <- strand_plus_RBS %>%
    mutate(RBS_hit = stringr::str_extract(TR, paste0(RBS, ".", range, "ATG"))) %>%
    mutate(spacer = nchar(RBS_hit) - nchar("ATG") - nchar(RBS)) %>%
    mutate(ORFstart = start + nchar(RBS) + spacer) %>%
    mutate(ORFmatch = ifelse(ORFstart %in% genbank_table$start, "match", "unmatch"))

  # Search for RBS on the bottom strand
  strand_minus_RBS <- vmatchPattern(
    most_frequent_RBS,
    toupper(contig_list_rc),
    max.mismatch = 2,
    min.mismatch = 0,
    with.indels = FALSE,
    fixed = TRUE
  ) %>% as.data.frame() %>%
    mutate(strand = "-")

  # Process matches for the bottom strand
  temp_df <- data.frame()
  for (i in 1:length(contig_list_rc)) {
    temp_RBS_df <- strand_minus_RBS %>%
      filter(group == i) %>%
      mutate(RBS = str_sub(contig_list_rc[i], start, end),
             TR = str_sub(contig_list_rc[i], start, end + 100))
    temp_df <- rbind.fill(temp_df, temp_RBS_df)
  }
  strand_minus_RBS <- temp_df

  # Extract hits and calculate spacers for the bottom strand
  strand_minus_RBS <- strand_minus_RBS %>%
    mutate(RBS_hit = stringr::str_extract(TR, paste0(RBS, ".", range, "ATG"))) %>%
    mutate(spacer = nchar(RBS_hit) - nchar("ATG") - nchar(RBS)) %>%
    mutate(ORFstart = start + nchar(RBS) + spacer)

  # Adjust ORFstart for reverse strand
  temp_df <- data.frame()
  for (i in 1:length(contig_list_rc)) {
    temp_RBS_df <- strand_minus_RBS %>%
      filter(group == i) %>%
      mutate(ORFstart = nchar(contig_list_rc[i]) - as.numeric(start + nchar(RBS) + spacer) + 1)
    temp_df <- rbind.fill(temp_df, temp_RBS_df)
  }
  strand_minus_RBS <- temp_df

  # Match ORFs for the reverse strand
  strand_minus_RBS <- strand_minus_RBS %>%
    mutate(ORFmatch = ifelse(ORFstart %in% (genbank_table %>% filter(direction == "-") %>% select(end) %>% pull), "match", "unmatch"))

  # Combine both strands
  putative_RBS <- rbind.fill(strand_plus_RBS, strand_minus_RBS)
  putative_matched_RBS <- putative_RBS %>% filter(ORFmatch == "match")

  # Assign the result to 'Putative RBS table' in the global environment
  assign("RBS_table", putative_RBS, envir = .GlobalEnv)
  assign("matched_RBS_table", putative_matched_RBS, envir = .GlobalEnv)
  message("The result has been saved to the R environment variable 'RBS_table'.")


  # Plot RBS sequence logo and spacer length histogram if requested
  if (plot_RBS) {
    # Sequence logo plot
    cs1 <- make_col_scheme(chars = c('A', 'T', 'C', 'G'), groups = c('A', 'T', 'G', 'C'),
                           cols = c('#76B56A', '#D06461', '#3B81A1', '#EFCA70'))
    seqlogo_plot <- ggseqlogo(putative_matched_RBS$RBS, col_scheme = cs1)+
      theme_bw() +
      theme(
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(angle = 90, size = 30),
        axis.text.x = element_text(color = "grey20", size = 30),
        axis.text.y = element_text(color = "grey20", size = 30),
        plot.title = element_text(size = 30)
      )

    # Spacer length histogram plot
    spacer_histogram <- ggplot(putative_matched_RBS, aes(x = spacer)) +
      geom_histogram(aes(y = ..density..), color = "black", fill = "#EBE5DE", binwidth = 1, alpha = 1.0) +
      scale_x_continuous(breaks = seq(1, 10, by = 1)) +
      theme_bw() +
      theme(
        panel.background = element_blank(),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(angle = 90, size = 30),
        axis.text.x = element_text(color = "grey20", size = 30),
        axis.text.y = element_text(color = "grey20", size = 30),
        plot.title = element_text(size = 30)
      ) +
      xlab("Spacer") +
      ylab("Density") +
      ggtitle("RBS to AUG")

    # Assign the plots to 'RBS seqlogo and spacer histogram' in the global environment
    assign("RBS_seqlogo", seqlogo_plot, envir = .GlobalEnv)
    message("The heatmap has been saved to the R environment variable 'RBS_seqlogo'.")

    assign("spacer_histogram", spacer_histogram, envir = .GlobalEnv)
    message("The heatmap has been saved to the R environment variable 'spacer_histogram'.")


    # Combine plots
    grid.arrange(seqlogo_plot, spacer_histogram, nrow = 2)

    # Save plot if requested
    if (save_plot) {
      if (is.null(plot_path)) {
        plot_path <- file.path(getwd(), "RBS_plots.pdf")
      }
      ggsave(plot_path, plot = grid.arrange(seqlogo_plot, spacer_histogram, nrow = 2))
      message(paste("Plot saved at:", plot_path))
    }
  }
}
