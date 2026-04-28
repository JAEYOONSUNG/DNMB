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

.dnmb_empty_rbs_table <- function() {
  data.frame(
    group = integer(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = character(),
    RBS = character(),
    TR = character(),
    RBS_hit = character(),
    spacer = integer(),
    ORFstart = integer(),
    ORFmatch = character(),
    stringsAsFactors = FALSE
  )
}

.dnmb_empty_rbs_plot <- function(title, subtitle) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5)
    )
}

.dnmb_rrna_product_mask <- function(product) {
  product <- tolower(trimws(as.character(product)))
  grepl("16s ribosomal rna|small subunit ribosomal rna|ssu rrna", product, perl = TRUE)
}

.dnmb_rbs_start_codons <- function(translation_domain = "bacteria") {
  if (identical(translation_domain, "archaea")) {
    return(c("ATG", "GTG", "TTG"))
  }
  "ATG"
}

.dnmb_assign_empty_rbs_outputs <- function(plot_RBS = FALSE, reason = "No 16S ribosomal RNA annotation detected.") {
  empty_rbs <- .dnmb_empty_rbs_table()
  assign("RBS_table", empty_rbs, envir = .GlobalEnv)
  assign("matched_RBS_table", empty_rbs, envir = .GlobalEnv)
  message(reason)
  message("The result has been saved to the R environment variable 'RBS_table'.")

  if (isTRUE(plot_RBS)) {
    seqlogo_plot <- .dnmb_empty_rbs_plot("RBS motif", reason)
    spacer_histogram <- .dnmb_empty_rbs_plot("RBS to AUG", reason)
    assign("RBS_seqlogo", seqlogo_plot, envir = .GlobalEnv)
    assign("spacer_histogram", spacer_histogram, envir = .GlobalEnv)
    message("Placeholder RBS plots have been saved to the global environment.")
  }

  invisible(empty_rbs)
}

.dnmb_archaea_upstream_window <- function(row, contig_list, flank = 25L) {
  contig_idx <- suppressWarnings(as.integer(row[["contig_number"]]))
  if (is.na(contig_idx) || contig_idx < 1L || contig_idx > length(contig_list)) {
    return(NULL)
  }
  seq_top <- as.character(contig_list[[contig_idx]])
  seq_len <- nchar(seq_top)
  start <- suppressWarnings(as.integer(row[["start"]]))
  end <- suppressWarnings(as.integer(row[["end"]]))
  strand <- as.character(row[["direction"]])
  rearranged <- as.character(row[["rearranged_nt_seq"]] %||% "")
  start_codon <- toupper(substr(rearranged, 1, 3))

  if (!nzchar(start_codon) || !(start_codon %in% .dnmb_rbs_start_codons("archaea"))) {
    return(NULL)
  }

  if (identical(strand, "+")) {
    left <- max(1L, start - flank)
    right <- start - 1L
    if (right < left) {
      return(NULL)
    }
    upstream <- substr(seq_top, left, right)
    return(list(
      upstream = toupper(upstream),
      start_codon = start_codon,
      gene_anchor = start,
      strand = strand,
      contig_idx = contig_idx,
      flank = nchar(upstream),
      genomic_left = left,
      genomic_right = right
    ))
  }

  if (identical(strand, "-")) {
    left <- end + 1L
    right <- min(seq_len, end + flank)
    if (right < left) {
      return(NULL)
    }
    genomic <- substr(seq_top, left, right)
    upstream <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(genomic)))
    return(list(
      upstream = toupper(upstream),
      start_codon = start_codon,
      gene_anchor = end,
      strand = strand,
      contig_idx = contig_idx,
      flank = nchar(upstream),
      genomic_left = left,
      genomic_right = right
    ))
  }

  NULL
}

.dnmb_archaea_candidate_kmers <- function(rrna_seed, min_len = 5L, max_len = 9L) {
  rrna_seed <- toupper(trimws(as.character(rrna_seed)[1]))
  if (is.na(rrna_seed) || !nzchar(rrna_seed)) {
    return(character())
  }
  n <- nchar(rrna_seed)
  kmers <- character()
  for (k in seq.int(min_len, min(max_len, n))) {
    if (k > n) {
      next
    }
    for (i in seq_len(n - k + 1L)) {
      kmers <- c(kmers, substr(rrna_seed, i, i + k - 1L))
    }
  }
  unique(kmers[nzchar(kmers)])
}

.dnmb_archaea_hamming <- function(a, b) {
  aa <- strsplit(a, "", fixed = TRUE)[[1]]
  bb <- strsplit(b, "", fixed = TRUE)[[1]]
  sum(aa != bb)
}

.dnmb_archaea_best_sd_hit <- function(upstream_seq, kmers) {
  upstream_seq <- toupper(trimws(as.character(upstream_seq)[1]))
  if (is.na(upstream_seq) || !nzchar(upstream_seq) || !length(kmers)) {
    return(NULL)
  }
  best <- NULL
  up_len <- nchar(upstream_seq)
  for (motif in kmers) {
    k <- nchar(motif)
    if (k > up_len) {
      next
    }
    for (pos in seq_len(up_len - k + 1L)) {
      hit <- substr(upstream_seq, pos, pos + k - 1L)
      mismatches <- .dnmb_archaea_hamming(hit, motif)
      spacer <- up_len - (pos + k - 1L)
      if (spacer < 0L || spacer > 20L) {
        next
      }
      score <- (k * 3) - (mismatches * 4)
      candidate <- list(
        motif = motif,
        hit = hit,
        pos = pos,
        len = k,
        mismatches = mismatches,
        spacer = spacer,
        score = score
      )
      if (is.null(best) ||
          candidate$score > best$score ||
          (identical(candidate$score, best$score) && candidate$len > best$len) ||
          (identical(candidate$score, best$score) && identical(candidate$len, best$len) && candidate$mismatches < best$mismatches)) {
        best <- candidate
      }
    }
  }
  best
}

.dnmb_archaea_plot_window9 <- function(upstream_seq, best_hit) {
  upstream_seq <- toupper(trimws(as.character(upstream_seq)[1]))
  if (is.null(best_hit) || is.na(upstream_seq) || !nzchar(upstream_seq)) {
    return(NA_character_)
  }
  up_len <- nchar(upstream_seq)
  hit_end <- best_hit$pos + best_hit$len - 1L
  plot_start <- max(1L, hit_end - 8L)
  plot_end <- min(up_len, plot_start + 8L)
  if ((plot_end - plot_start + 1L) < 9L) {
    plot_start <- max(1L, plot_end - 8L)
  }
  out <- substr(upstream_seq, plot_start, plot_start + 8L)
  if (nchar(out) != 9L) {
    return(NA_character_)
  }
  out
}

.dnmb_archaea_estimate_spacer_peak <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl) || !nrow(tbl)) {
    return(NA_integer_)
  }
  work <- tbl
  work$score <- suppressWarnings(as.numeric(work$score))
  work$spacer <- suppressWarnings(as.integer(work$spacer))
  work <- work[!is.na(work$score) & !is.na(work$spacer), , drop = FALSE]
  if (!nrow(work)) {
    return(NA_integer_)
  }
  work$weight <- pmax(0.5, work$score - min(work$score, na.rm = TRUE) + 1)
  spacer_tbl <- work |>
    dplyr::group_by(.data$spacer) |>
    dplyr::summarise(weight = sum(.data$weight, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$weight), dplyr::desc(.data$n), .data$spacer)
  suppressWarnings(as.integer(spacer_tbl$spacer[[1]]))
}

.dnmb_archaea_status_plots <- function(status_text) {
  list(
    seqlogo_plot = .dnmb_empty_rbs_plot("SD-like motif", status_text),
    spacer_histogram = .dnmb_empty_rbs_plot("Initiation architecture", status_text)
  )
}

.dnmb_archaea_rbs_extractor <- function(target, contig_list, plot_RBS = FALSE) {
  rrna_tbl <- target[.dnmb_rrna_product_mask(target$product), , drop = FALSE]
  rrna_tbl <- rrna_tbl[!is.na(rrna_tbl$rearranged_nt_seq) & nzchar(as.character(rrna_tbl$rearranged_nt_seq)), , drop = FALSE]
  if (!nrow(rrna_tbl)) {
    empty <- .dnmb_empty_rbs_table()
    assign("RBS_table", empty, envir = .GlobalEnv)
    assign("matched_RBS_table", empty, envir = .GlobalEnv)
    if (isTRUE(plot_RBS)) {
      plots <- .dnmb_archaea_status_plots("No 16S/SSU ribosomal RNA annotation detected; skipping archaeal SD-like analysis.")
      assign("RBS_seqlogo", plots$seqlogo_plot, envir = .GlobalEnv)
      assign("spacer_histogram", plots$spacer_histogram, envir = .GlobalEnv)
    }
    message("No 16S/SSU ribosomal RNA annotation detected; skipping archaeal SD-like analysis.")
    message("The result has been saved to the R environment variable 'RBS_table'.")
    return(invisible(empty))
  }

  rrna_tail <- unique(vapply(rrna_tbl$rearranged_nt_seq, function(x) {
    seq <- toupper(as.character(x))
    tail <- substr(seq, max(1L, nchar(seq) - 11L), nchar(seq))
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(tail)))
  }, character(1)))
  kmers <- unique(unlist(lapply(rrna_tail, .dnmb_archaea_candidate_kmers), use.names = FALSE))
  kmers <- kmers[nchar(kmers) >= 5L]

  gene_tbl <- target[
    !is.na(target$start) &
      !is.na(target$end) &
      !is.na(target$direction) &
      !is.na(target$rearranged_nt_seq) &
      nzchar(as.character(target$rearranged_nt_seq)),
    , drop = FALSE
  ]
  gene_tbl <- gene_tbl[!grepl("ribosomal RNA", gene_tbl$product %||% "", ignore.case = TRUE), , drop = FALSE]

  rows <- list()
  for (i in seq_len(nrow(gene_tbl))) {
    row <- gene_tbl[i, , drop = FALSE]
    win <- .dnmb_archaea_upstream_window(row, contig_list = contig_list, flank = 25L)
    if (is.null(win) || !nzchar(win$upstream)) {
      next
    }
    best <- .dnmb_archaea_best_sd_hit(win$upstream, kmers)
    if (is.null(best) || best$len < 5L) {
      next
    }

    if (identical(win$strand, "+")) {
      motif_start <- win$genomic_left + best$pos - 1L
      motif_end <- motif_start + best$len - 1L
    } else {
      motif_end <- win$genomic_right - best$pos + 1L
      motif_start <- motif_end - best$len + 1L
    }

    rows[[length(rows) + 1L]] <- data.frame(
      group = suppressWarnings(as.integer(row$contig_number)),
      start = motif_start,
      end = motif_end,
      width = best$len,
      strand = as.character(row$direction),
      RBS = best$hit,
      TR = paste0(win$upstream, win$start_codon),
      RBS_hit = paste0(best$hit, strrep("N", best$spacer), win$start_codon),
      spacer = best$spacer,
      ORFstart = if (identical(win$strand, "+")) suppressWarnings(as.integer(row$start)) else suppressWarnings(as.integer(row$end)),
      ORFmatch = "match",
      locus_tag = as.character(row$locus_tag),
      start_codon = win$start_codon,
      score = best$score,
      mismatches = best$mismatches,
      motif_seed = best$motif,
      RBS_plot9 = .dnmb_archaea_plot_window9(win$upstream, best),
      stringsAsFactors = FALSE
    )
  }

  putative_RBS <- if (length(rows)) dplyr::bind_rows(rows) else .dnmb_empty_rbs_table()
  matched_RBS_table <- putative_RBS

  assign("RBS_table", putative_RBS, envir = .GlobalEnv)
  assign("matched_RBS_table", matched_RBS_table, envir = .GlobalEnv)
  message("The result has been saved to the R environment variable 'RBS_table'.")

  if (isTRUE(plot_RBS)) {
    if (!nrow(matched_RBS_table)) {
      plots <- .dnmb_archaea_status_plots("No archaeal SD-like candidates were detected.")
      assign("RBS_seqlogo", plots$seqlogo_plot, envir = .GlobalEnv)
      assign("spacer_histogram", plots$spacer_histogram, envir = .GlobalEnv)
      message("Placeholder archaeal initiation plots have been saved to the global environment.")
      return(invisible(putative_RBS))
    }

    score_cutoff <- max(stats::quantile(matched_RBS_table$score, probs = 0.85, na.rm = TRUE, names = FALSE), 13)
    strong_pool <- matched_RBS_table[
      matched_RBS_table$score >= score_cutoff &
        matched_RBS_table$width >= 6L &
        matched_RBS_table$mismatches <= 1L,
      , drop = FALSE
    ]
    if (!nrow(strong_pool) && nrow(matched_RBS_table)) {
      strong_pool <- matched_RBS_table[order(-matched_RBS_table$score, matched_RBS_table$mismatches), , drop = FALSE]
      strong_pool <- utils::head(strong_pool, min(20L, nrow(strong_pool)))
      strong_pool <- strong_pool[strong_pool$score >= 10, , drop = FALSE]
    }

    if (nrow(strong_pool) < 4L) {
      plots <- .dnmb_archaea_status_plots("No strong archaeal SD-like subset was detected; this genome may be weak-SD or leaderless-rich.")
      assign("RBS_seqlogo", plots$seqlogo_plot, envir = .GlobalEnv)
      assign("spacer_histogram", plots$spacer_histogram, envir = .GlobalEnv)
      message("Strong archaeal SD-like subset was not detected; status plots were saved instead.")
      return(invisible(putative_RBS))
    }

    spacer_peak <- .dnmb_archaea_estimate_spacer_peak(strong_pool)
    if (is.na(spacer_peak)) {
      spacer_peak <- 6L
    }
    strong <- strong_pool[abs(strong_pool$spacer - spacer_peak) <= 2L, , drop = FALSE]
    if (nrow(strong) < 4L) {
      strong <- strong_pool[abs(strong_pool$spacer - spacer_peak) <= 3L, , drop = FALSE]
    }
    if (nrow(strong) < 4L) {
      strong <- strong_pool
    }

    strong$RBS_plot9 <- as.character(strong$RBS_plot9)
    strong_logo <- strong[!is.na(strong$RBS_plot9) & nchar(strong$RBS_plot9) == 9L, , drop = FALSE]
    if (nrow(strong_logo) < 4L) {
      plots <- .dnmb_archaea_status_plots("Strong archaeal SD-like candidates were found, but motif lengths were too heterogeneous for a stable sequence logo.")
      assign("RBS_seqlogo", plots$seqlogo_plot, envir = .GlobalEnv)
      assign("spacer_histogram", plots$spacer_histogram, envir = .GlobalEnv)
      message("Archaeal SD-like candidates were heterogeneous; status plots were saved instead.")
      return(invisible(putative_RBS))
    }

    cs1 <- ggseqlogo::make_col_scheme(
      chars = c("A", "T", "C", "G"),
      groups = c("A", "T", "G", "C"),
      cols = c("#76B56A", "#D06461", "#3B81A1", "#EFCA70")
    )
    seqlogo_plot <- ggseqlogo::ggseqlogo(strong_logo$RBS_plot9, col_scheme = cs1) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "RBS motif") +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 30),
        axis.title.x = ggplot2::element_text(size = 30),
        axis.title.y = ggplot2::element_text(angle = 90, size = 30),
        axis.text.x = ggplot2::element_text(color = "grey20", size = 30),
        axis.text.y = ggplot2::element_text(color = "grey20", size = 30),
        plot.title = ggplot2::element_text(size = 30)
      )

    spacer_histogram <- ggplot2::ggplot(strong, ggplot2::aes(x = .data$spacer)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(density)),
        color = "black", fill = "#EBE5DE", binwidth = 1, alpha = 1
      ) +
      ggplot2::scale_x_continuous(breaks = seq(1, 10, by = 1), limits = c(1, 10)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 30),
        axis.title.y = ggplot2::element_text(angle = 90, size = 30),
        axis.text.x = ggplot2::element_text(color = "grey20", size = 30),
        axis.text.y = ggplot2::element_text(color = "grey20", size = 30),
        plot.title = ggplot2::element_text(size = 30)
      ) +
      ggplot2::xlab("Spacer") +
      ggplot2::ylab("Density") +
      ggplot2::ggtitle("RBS to AUG")

    assign("RBS_seqlogo", seqlogo_plot, envir = .GlobalEnv)
    assign("spacer_histogram", spacer_histogram, envir = .GlobalEnv)
    message("Archaeal SD-like plots have been saved to the global environment.")
  }

  invisible(putative_RBS)
}

RBS_extractor <- function(target = NULL,
                          plot_RBS = FALSE,
                          save_plot = FALSE,
                          plot_path = NULL,
                          translation_domain = NULL,
                          gb_path = NULL) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggseqlogo))
  suppressPackageStartupMessages(library(gridExtra))

  # If target is NULL, check if 'genbank_table' exists in the global environment
  if (is.null(target)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      target <- get("genbank_table", envir = .GlobalEnv)
      message("Using 'genbank_table' from the global environment as target.")
    } else {
      stop("You need to provide a 'target' or ensure 'genbank_table' exists in the global environment.")
    }
  }

  translation_domain <- translation_domain %||% .dnmb_detect_translation_domain(target = target, gb_path = gb_path)
  if (identical(translation_domain, "archaea")) {
    contig_names <- ls(envir = .GlobalEnv, pattern = "^contig_[0-9]{1,}_seq$")
    if (!length(contig_names)) {
      return(.dnmb_assign_empty_rbs_outputs(
        plot_RBS = plot_RBS,
        reason = "No contig sequences found in the global environment for archaeal initiation analysis."
      ))
    }
    contig_list <- lapply(contig_names, function(x) get(x, envir = .GlobalEnv))
    return(.dnmb_archaea_rbs_extractor(target = target, contig_list = contig_list, plot_RBS = plot_RBS))
  }

  # Extract RBS sequences based on 16S ribosomal RNA
  contig_names <- ls(envir = .GlobalEnv, pattern = "^contig_[0-9]{1,}_seq$")

  if (!length(contig_names)) {
    return(.dnmb_assign_empty_rbs_outputs(
      plot_RBS = plot_RBS,
      reason = "No contig sequences found in the global environment for RBS extraction."
    ))
  }

  # 리스트의 모든 데이터를 rbind로 결합
  contig_list <- lapply(contig_names, function(x) get(x, envir = .GlobalEnv))
  contig_list_rc <- lapply(contig_names, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(get(x, envir = .GlobalEnv)))))

  RBS <- target %>%
    filter(.dnmb_rrna_product_mask(product)) %>%
    select(rearranged_nt_seq) %>%
    mutate(RBS = stringr::str_sub(rearranged_nt_seq, start = -9))

  if (!nrow(RBS)) {
    return(.dnmb_assign_empty_rbs_outputs(
      plot_RBS = plot_RBS,
      reason = "No 16S ribosomal RNA annotation detected; skipping RBS extraction."
    ))
  }

  # Reverse complement for RBS sequences
  RBS$RBS <- sapply(RBS$RBS, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))
  RBS <- RBS[!is.na(RBS$RBS) & nzchar(RBS$RBS), , drop = FALSE]

  if (!nrow(RBS)) {
    return(.dnmb_assign_empty_rbs_outputs(
      plot_RBS = plot_RBS,
      reason = "No valid 16S-derived RBS seed sequence detected; skipping RBS extraction."
    ))
  }

  # Find the most frequent RBS sequence
  most_frequent_RBS <- names(sort(table(toupper(RBS$RBS)), decreasing = TRUE))[1]

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
      mutate(RBS = stringr::str_sub(contig_list[i], start, end),
             TR = stringr::str_sub(contig_list[i], start, end + 100))
    temp_df <- plyr::rbind.fill(temp_df, temp_RBS_df)
  }
  strand_plus_RBS <- temp_df

  # Set range and extract RBS hits
  range <- "{1,10}"
  strand_plus_RBS <- strand_plus_RBS %>%
    mutate(RBS_hit = stringr::str_extract(TR, paste0(RBS, ".", range, "ATG"))) %>%
    mutate(spacer = nchar(RBS_hit) - nchar("ATG") - nchar(RBS)) %>%
    mutate(ORFstart = start + nchar(RBS) + spacer) %>%
    mutate(ORFmatch = ifelse(ORFstart %in% target$start, "match", "unmatch"))

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
      mutate(RBS = stringr::str_sub(contig_list_rc[i], start, end),
             TR = stringr::str_sub(contig_list_rc[i], start, end + 100))
    temp_df <- plyr::rbind.fill(temp_df, temp_RBS_df)
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
    temp_df <- plyr::rbind.fill(temp_df, temp_RBS_df)
  }
  strand_minus_RBS <- temp_df

  # Match ORFs for the reverse strand
  strand_minus_RBS <- strand_minus_RBS %>%
    mutate(ORFmatch = ifelse(ORFstart %in% (target %>% filter(direction == "-") %>% select(end) %>% pull), "match", "unmatch"))

  # Combine both strands
  putative_RBS <- plyr::rbind.fill(strand_plus_RBS, strand_minus_RBS)
  putative_matched_RBS <- putative_RBS %>% filter(ORFmatch == "match")

  # Assign the result to 'Putative RBS table' in the global environment
  assign("RBS_table", putative_RBS, envir = .GlobalEnv)
  assign("matched_RBS_table", putative_matched_RBS, envir = .GlobalEnv)
  message("The result has been saved to the R environment variable 'RBS_table'.")


  # Plot RBS sequence logo and spacer length histogram if requested
  if (plot_RBS) {
    if (!nrow(putative_matched_RBS)) {
      seqlogo_plot <- .dnmb_empty_rbs_plot("RBS motif", "No 16S-linked RBS matches were detected.")
      spacer_histogram <- .dnmb_empty_rbs_plot("RBS to AUG", "No 16S-linked RBS matches were detected.")
      assign("RBS_seqlogo", seqlogo_plot, envir = .GlobalEnv)
      assign("spacer_histogram", spacer_histogram, envir = .GlobalEnv)
      message("Placeholder RBS plots have been saved to the global environment.")
      return(invisible(putative_RBS))
    }

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


    # Save plot if requested
    if (save_plot) {
      if (is.null(plot_path)) {
        plot_path <- file.path(getwd(), "RBS_plots.pdf")
      }
      combined_plot <- arrangeGrob(seqlogo_plot, spacer_histogram, nrow = 2)
      ggsave(plot_path, plot = combined_plot)
      message(paste("Plot saved at:", plot_path))
    }
  }
}
