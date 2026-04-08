# Mobileome_insertion_verify.R — Cross-genome IS insertion verification
#
# For each IS element in the focal genome, BLAST flanking sequences against
# related genomes to verify insertions and recover validated TSD sequences.
# Adapted from is110coverage blast_verify approach.

# ---------------------------------------------------------------------------
# Main entry: verify all IS insertions against related genomes
# ---------------------------------------------------------------------------
.dnmb_verify_is_insertions <- function(
  elements,
  focal_metadata,
  related_genbanks,
  output_dir,
  flank_bp = 500L,
  evalue = 1e-5,
  min_identity = 70,
  max_empty_gap = 50L,
  verbose = TRUE
) {
  blastn <- Sys.which("blastn")
  makeblastdb <- Sys.which("makeblastdb")
  if (!nzchar(blastn) || !nzchar(makeblastdb)) {
    .dnmb_iselement_verbose(verbose, "BLAST+ not found. Skipping insertion verification.")
    return(.dnmb_empty_verification())
  }
  if (!nrow(elements) || !length(related_genbanks)) {
    return(.dnmb_empty_verification())
  }

  verify_dir <- file.path(output_dir, "insertion_verify")
  dir.create(verify_dir, recursive = TRUE, showWarnings = FALSE)

  # Build focal sequence map
  focal_seq_map <- stats::setNames(focal_metadata$sequence, focal_metadata$contig)

  all_results <- list()
  all_tsd <- list()

  for (ref_idx in seq_along(related_genbanks)) {
    ref_path <- related_genbanks[ref_idx]
    if (!file.exists(ref_path)) next
    .dnmb_iselement_verbose(verbose, paste0("Verifying against: ", basename(ref_path)))

    # Parse reference genome sequences
    ref_seqs <- tryCatch({
      ref_parsed <- .dnmb_parse_genbank_features(ref_path)
      stats::setNames(ref_parsed$metadata$sequence, ref_parsed$metadata$contig)
    }, error = function(e) character())
    if (!length(ref_seqs)) next

    # Build BLAST database for this reference (use tempdir to avoid path-with-spaces issues)
    db_dir <- tempdir()
    db_prefix <- file.path(db_dir, paste0("dnmb_ref_", ref_idx))
    ref_fasta <- paste0(db_prefix, ".fasta")
    .dnmb_write_named_fasta(ref_seqs, ref_fasta)
    system2(makeblastdb, c("-in", ref_fasta, "-dbtype", "nucl", "-out", db_prefix),
            stdout = FALSE, stderr = FALSE)

    # Verify each IS element
    for (el_idx in seq_len(nrow(elements))) {
      el <- elements[el_idx, , drop = FALSE]
      ctg <- el$contig
      focal_seq <- focal_seq_map[[ctg]]
      if (is.null(focal_seq) || !nzchar(focal_seq)) next

      el_start <- el$start; el_end <- el$end
      el_len <- el_end - el_start + 1L

      # Extract flanking sequences
      left_start <- max(1L, el_start - flank_bp)
      left_end <- el_start - 1L
      right_start <- el_end + 1L
      right_end <- min(nchar(focal_seq), el_end + flank_bp)

      left_flank <- substr(focal_seq, left_start, left_end)
      right_flank <- substr(focal_seq, right_start, right_end)

      if (nchar(left_flank) < 50 || nchar(right_flank) < 50) next

      # BLAST each flank
      left_hits <- .dnmb_blast_flank(
        flank_seq = left_flank, flank_name = paste0("el", el_idx, "_L"),
        db_prefix = db_prefix, evalue = evalue, min_identity = min_identity
      )
      right_hits <- .dnmb_blast_flank(
        flank_seq = right_flank, flank_name = paste0("el", el_idx, "_R"),
        db_prefix = db_prefix, evalue = evalue, min_identity = min_identity
      )

      # Classify insertion
      cls <- .dnmb_classify_insertion(left_hits, right_hits, el_len, max_empty_gap)

      result_row <- tibble::tibble(
        element_id = el$element_id,
        element_family = el$element_family,
        contig = ctg,
        start = el_start,
        end = el_end,
        element_length = el_len,
        ref_file = basename(ref_path),
        ref_contig = cls$ref_contig,
        status = cls$status,
        gap_bp = cls$gap_bp,
        left_pident = cls$left_pident,
        right_pident = cls$right_pident,
        left_ref_end = cls$left_ref_end,
        right_ref_start = cls$right_ref_start
      )
      all_results[[length(all_results) + 1]] <- result_row

      # If empty_site found, extract validated TSD from reference gap
      if (cls$status == "empty_site" && !is.na(cls$gap_bp) && cls$gap_bp >= 0 && cls$gap_bp <= 20) {
        tsd_seq <- .dnmb_extract_gap_tsd(
          ref_seqs = ref_seqs,
          ref_contig = cls$ref_contig,
          left_ref_end = cls$left_ref_end,
          right_ref_start = cls$right_ref_start
        )
        if (!is.na(tsd_seq) && nzchar(tsd_seq)) {
          # Normalize to IS coding direction: if IS is on minus strand, revcomp the TSD
          el_strand <- if ("strand" %in% names(el)) el$strand else "+"
          if (!is.na(el_strand) && el_strand == "-") {
            tsd_seq <- .dnmb_revcomp_tsd(tsd_seq)
          }
          all_tsd[[length(all_tsd) + 1]] <- tibble::tibble(
            element_id = el$element_id,
            family = el$element_family,
            tsd_seq = tsd_seq,
            tsd_len_bp = nchar(tsd_seq),
            source = "comparative_verified",
            ref_file = basename(ref_path),
            ref_contig = cls$ref_contig,
            left_pident = cls$left_pident,
            right_pident = cls$right_pident
          )
        }
      }
    }
  }

  verification <- dplyr::bind_rows(all_results)
  validated_tsd <- dplyr::bind_rows(all_tsd)

  # Write outputs
  if (nrow(verification)) {
    utils::write.table(verification, file.path(verify_dir, "insertion_verification.tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE)
  }
  if (nrow(validated_tsd)) {
    utils::write.table(validated_tsd, file.path(verify_dir, "validated_tsd.tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE)
  }

  .dnmb_iselement_verbose(verbose, paste0(
    "Insertion verification: ", nrow(verification), " checks, ",
    sum(verification$status == "empty_site", na.rm = TRUE), " empty sites, ",
    nrow(validated_tsd), " validated TSDs"
  ))

  list(
    verification = verification,
    validated_tsd = validated_tsd
  )
}

.dnmb_empty_verification <- function() {
  list(
    verification = tibble::tibble(
      element_id = character(), element_family = character(),
      contig = character(), start = integer(), end = integer(),
      element_length = integer(), ref_file = character(),
      ref_contig = character(), status = character(),
      gap_bp = integer(), left_pident = numeric(), right_pident = numeric(),
      left_ref_end = integer(), right_ref_start = integer()
    ),
    validated_tsd = tibble::tibble(
      element_id = character(), family = character(),
      tsd_seq = character(), tsd_len_bp = integer(),
      source = character(), ref_file = character(),
      ref_contig = character(), left_pident = numeric(), right_pident = numeric()
    )
  )
}

# ---------------------------------------------------------------------------
# BLAST a single flank sequence
# ---------------------------------------------------------------------------
.dnmb_blast_flank <- function(flank_seq, flank_name, db_prefix, evalue = 1e-5, min_identity = 70) {
  query_file <- tempfile(fileext = ".fasta")
  on.exit(unlink(query_file), add = TRUE)
  writeLines(c(paste0(">", flank_name), flank_seq), query_file)

  out_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(out_file), add = TRUE)

  system2("blastn", c(
    "-query", query_file,
    "-db", db_prefix,
    "-outfmt", "6",
    "-evalue", as.character(evalue),
    "-perc_identity", as.character(min_identity),
    "-max_target_seqs", "5",
    "-out", out_file
  ), stdout = FALSE, stderr = FALSE)

  if (!file.exists(out_file) || file.size(out_file) < 5) {
    return(data.frame(sseqid = character(), pident = numeric(), length = integer(),
                      sstart = integer(), send = integer(), bitscore = numeric(),
                      stringsAsFactors = FALSE))
  }

  tryCatch({
    hits <- utils::read.delim(out_file, header = FALSE, check.names = FALSE)
    names(hits) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore")
    hits
  }, error = function(e) {
    data.frame(sseqid = character(), pident = numeric(), length = integer(),
               sstart = integer(), send = integer(), bitscore = numeric(),
               stringsAsFactors = FALSE)
  })
}

# ---------------------------------------------------------------------------
# Classify insertion by comparing flank positions in reference
# ---------------------------------------------------------------------------
.dnmb_classify_insertion <- function(left_hits, right_hits, is_element_length, max_empty_gap = 50L) {
  result <- list(status = "absent", gap_bp = NA_integer_, ref_contig = NA_character_,
                 left_pident = NA_real_, right_pident = NA_real_,
                 left_ref_end = NA_integer_, right_ref_start = NA_integer_)

  no_left <- is.null(left_hits) || !nrow(left_hits)
  no_right <- is.null(right_hits) || !nrow(right_hits)

  if (no_left && no_right) return(result)

  if (no_left || no_right) {
    result$status <- "partial"
    if (!no_left) {
      best <- left_hits[which.max(left_hits$bitscore), , drop = FALSE]
      result$left_pident <- best$pident; result$ref_contig <- best$sseqid
    } else {
      best <- right_hits[which.max(right_hits$bitscore), , drop = FALSE]
      result$right_pident <- best$pident; result$ref_contig <- best$sseqid
    }
    return(result)
  }

  # Both flanks have hits — find best pair on same contig
  shared_contigs <- intersect(left_hits$sseqid, right_hits$sseqid)
  if (!length(shared_contigs)) {
    result$status <- "split"
    return(result)
  }

  best_gap <- Inf
  best_pair <- NULL
  for (ctg in shared_contigs) {
    l_sub <- left_hits[left_hits$sseqid == ctg, , drop = FALSE]
    r_sub <- right_hits[right_hits$sseqid == ctg, , drop = FALSE]
    for (li in seq_len(nrow(l_sub))) {
      for (ri in seq_len(nrow(r_sub))) {
        l_end <- max(l_sub$sstart[li], l_sub$send[li])
        r_start <- min(r_sub$sstart[ri], r_sub$send[ri])
        gap <- r_start - l_end - 1L
        if (abs(gap) < abs(best_gap)) {
          best_gap <- gap
          best_pair <- list(
            left = l_sub[li, , drop = FALSE],
            right = r_sub[ri, , drop = FALSE],
            contig = ctg, gap = gap,
            left_ref_end = l_end,
            right_ref_start = r_start
          )
        }
      }
    }
  }

  if (is.null(best_pair)) {
    result$status <- "rearranged"
    return(result)
  }

  result$gap_bp <- as.integer(best_pair$gap)
  result$ref_contig <- best_pair$contig
  result$left_pident <- best_pair$left$pident
  result$right_pident <- best_pair$right$pident
  result$left_ref_end <- best_pair$left_ref_end
  result$right_ref_start <- best_pair$right_ref_start

  abs_gap <- abs(best_pair$gap)
  if (abs_gap <= max_empty_gap) {
    result$status <- "empty_site"
  } else if (abs_gap >= is_element_length * 0.8 && abs_gap <= is_element_length * 1.2) {
    result$status <- "filled_site"
  } else {
    result$status <- "rearranged"
  }

  result
}

# ---------------------------------------------------------------------------
# Extract TSD from reference genome gap between flanks
# ---------------------------------------------------------------------------
.dnmb_extract_gap_tsd <- function(ref_seqs, ref_contig, left_ref_end, right_ref_start) {
  if (is.na(ref_contig) || is.na(left_ref_end) || is.na(right_ref_start)) return(NA_character_)
  seq <- ref_seqs[[ref_contig]]
  if (is.null(seq) || !nzchar(seq)) return(NA_character_)

  gap_start <- left_ref_end + 1L
  gap_end <- right_ref_start - 1L
  if (gap_start > gap_end || gap_end > nchar(seq)) return(NA_character_)

  tsd <- toupper(substr(seq, gap_start, gap_end))
  if (!grepl("^[ACGT]+$", tsd)) return(NA_character_)
  tsd
}

# ---------------------------------------------------------------------------
# Reverse complement for TSD normalization to IS coding direction
# ---------------------------------------------------------------------------
.dnmb_revcomp_tsd <- function(seq) {
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  chars <- rev(strsplit(toupper(seq), "")[[1]])
  paste(comp[chars], collapse = "")
}

# ---------------------------------------------------------------------------
# Write named character vector as FASTA
# ---------------------------------------------------------------------------
.dnmb_write_named_fasta <- function(seqs, path) {
  con <- file(path, "w")
  on.exit(close(con))
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm), con)
    s <- seqs[[nm]]
    n <- nchar(s)
    starts <- seq(1, n, 70)
    for (st in starts) writeLines(substr(s, st, min(st + 69, n)), con)
  }
}
