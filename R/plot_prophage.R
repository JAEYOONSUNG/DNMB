.dnmb_prophage_reference_cache_dir <- function(output_dir, region_id) {
  dir <- file.path(output_dir, "dnmb_module_prophage", "_cache", paste0("prophage_ref_", region_id))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(dir, winslash = "/", mustWork = FALSE)
}

.dnmb_normalize_contig_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\bcomplete genome\\b", "complete", x)
  x <- gsub("\\bcomplete sequence\\b", "complete", x)
  x <- gsub("\\s+", " ", x)
  x
}

.dnmb_write_single_fasta <- function(header, sequence, path) {
  con <- base::file(path, "w")
  on.exit(base::close(con), add = TRUE)
  base::writeLines(paste0(">", header), con = con)
  base::writeLines(sequence, con = con)
  invisible(path)
}

.dnmb_prophage_ensure_viral_db <- function(cache_root = NULL) {
  db_dir <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = TRUE), "prophage", "viral_refdb")
  db_file <- file.path(db_dir, "ref_viruses_rep_genomes.nin")
  if (file.exists(db_file)) return(file.path(db_dir, "ref_viruses_rep_genomes"))
  dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
  message("[Prophage] Downloading NCBI viral reference genome DB to ", db_dir, "...")
  run <- dnmb_run_external("update_blastdb.pl", args = c("--decompress", "ref_viruses_rep_genomes"), wd = db_dir, required = FALSE)
  if (!isTRUE(run$ok)) {
    # Fallback: manual download
    for (ext in c("ndb","nhr","nin","njs","nog","nos","not","nsq","ntf","nto")) {
      url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/db/ref_viruses_rep_genomes.", ext)
      dnmb_run_external("curl", args = c("-L", "-k", "-o", file.path(db_dir, paste0("ref_viruses_rep_genomes.", ext)), url), required = FALSE)
    }
  }
  if (file.exists(db_file)) file.path(db_dir, "ref_viruses_rep_genomes") else NULL
}

.dnmb_prophage_remote_reference_hits <- function(query_fasta, cache_dir, top_n = 5L, cache_root = NULL) {
  out_tsv <- file.path(cache_dir, "reference_hits.tsv")
  if (!file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
    # Try local DB first
    local_db <- .dnmb_prophage_ensure_viral_db(cache_root = cache_root)
    if (!is.null(local_db)) {
      args <- c(
        "-task", "megablast",
        "-db", local_db,
        "-query", query_fasta,
        "-max_target_seqs", "30",
        "-evalue", "1e-10",
        "-outfmt", "6 saccver pident length qcovs bitscore stitle",
        "-out", out_tsv
      )
    } else {
      # Fallback: remote BLAST
      args <- c(
        "-task", "megablast",
        "-db", "ref_viruses_rep_genomes",
        "-remote",
        "-query", query_fasta,
        "-max_target_seqs", "30",
        "-evalue", "1e-10",
        "-outfmt", "6 saccver pident length qcovs bitscore stitle",
        "-out", out_tsv
      )
    }
    run <- dnmb_run_external("blastn", args = args, required = FALSE)
    if (!isTRUE(run$ok)) {
      return(data.frame())
    }
  }
  if (!file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
    return(data.frame())
  }
  tbl <- utils::read.table(out_tsv, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    return(data.frame())
  }
  names(tbl) <- c("accession", "pident", "length", "qcovs", "bitscore", "title")
  tbl$accession <- as.character(tbl$accession)
  tbl$title <- as.character(tbl$title)
  tbl$pident <- suppressWarnings(as.numeric(tbl$pident))
  tbl$qcovs <- suppressWarnings(as.numeric(tbl$qcovs))
  tbl$bitscore <- suppressWarnings(as.numeric(tbl$bitscore))
  tbl <- tbl[!duplicated(tbl$accession), , drop = FALSE]
  tbl <- tbl[order(-tbl$qcovs, -tbl$bitscore, -tbl$pident), , drop = FALSE]
  head(tbl, max(1L, as.integer(top_n)[1]))
}

.dnmb_fetch_genbank_text <- function(accession, dest_path) {
  if (file.exists(dest_path) && file.info(dest_path)$size > 0) {
    return(dest_path)
  }
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",
    accession,
    "&rettype=gbwithparts&retmode=text"
  )
  run <- dnmb_run_external("curl", args = c("-L", "-k", url, "-o", dest_path), required = FALSE)
  if (!isTRUE(run$ok) || !file.exists(dest_path) || file.info(dest_path)$size == 0) {
    return("")
  }
  dest_path
}

.dnmb_parse_single_genbank_record <- function(path) {
  if (!file.exists(path)) {
    return(list(sequence = "", cds = data.frame()))
  }
  lines <- readLines(path, warn = FALSE)
  origin_idx <- grep("^ORIGIN", lines)
  end_idx <- grep("^//", lines)
  sequence <- ""
  if (length(origin_idx) && length(end_idx) && end_idx[[1]] > origin_idx[[1]]) {
    sequence <- paste(lines[(origin_idx[[1]] + 1L):(end_idx[[1]] - 1L)], collapse = "")
    sequence <- gsub("[^ACGTacgt]", "", sequence)
    sequence <- toupper(sequence)
  }
  feat_idx <- grep("^FEATURES", lines)
  cds_tbl <- data.frame()
  if (length(feat_idx) && length(origin_idx) && origin_idx[[1]] > feat_idx[[1]]) {
    feat_lines <- lines[(feat_idx[[1]] + 1L):(origin_idx[[1]] - 1L)]
    entries <- list()
    i <- 1L
    while (i <= length(feat_lines)) {
      line <- feat_lines[[i]]
      if (grepl("^     CDS\\s+", line)) {
        location <- trimws(sub("^     CDS\\s+", "", line))
        j <- i + 1L
        qualifiers <- character()
        while (j <= length(feat_lines) && grepl("^                     ", feat_lines[[j]])) {
          qualifiers <- c(qualifiers, trimws(feat_lines[[j]]))
          j <- j + 1L
        }
        starts <- suppressWarnings(as.numeric(unlist(regmatches(location, gregexpr("[0-9]+", location)))))
        if (length(starts)) {
          gene <- NA_character_
          locus_tag <- NA_character_
          product <- NA_character_
          if (length(qualifiers)) {
            gene <- sub('^/gene=\"([^\"]+)\".*$', '\\1', qualifiers[grepl('^/gene=\"', qualifiers)][1])
            locus_tag <- sub('^/locus_tag=\"([^\"]+)\".*$', '\\1', qualifiers[grepl('^/locus_tag=\"', qualifiers)][1])
            product <- sub('^/product=\"([^\"]+)\".*$', '\\1', qualifiers[grepl('^/product=\"', qualifiers)][1])
          }
          entries[[length(entries) + 1L]] <- data.frame(
            start = min(starts),
            end = max(starts),
            direction = if (grepl("complement", location)) "-" else "+",
            gene = ifelse(is.na(gene), NA_character_, gene),
            locus_tag = ifelse(is.na(locus_tag), NA_character_, locus_tag),
            product = ifelse(is.na(product), NA_character_, product),
            stringsAsFactors = FALSE
          )
        }
        i <- j
        next
      }
      i <- i + 1L
    }
    cds_tbl <- dplyr::bind_rows(entries)
  }
  list(sequence = sequence, cds = cds_tbl)
}

.dnmb_prophage_pairwise_alignment <- function(query_fasta, subject_fasta, out_tsv) {
  if (!file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
    run <- dnmb_run_external(
      "blastn",
      args = c(
        "-task", "blastn",
        "-query", query_fasta,
        "-subject", subject_fasta,
        "-dust", "no",
        "-max_hsps", "200",
        "-outfmt", "6 qstart qend sstart send pident length mismatch gaps qseq sseq sstrand",
        "-out", out_tsv
      ),
      required = FALSE
    )
    if (!isTRUE(run$ok)) {
      return(data.frame())
    }
  }
  if (!file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
    return(data.frame())
  }
  tbl <- utils::read.table(out_tsv, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    return(data.frame())
  }
  names(tbl) <- c("qstart", "qend", "sstart", "send", "pident", "length", "mismatch", "gaps", "qseq", "sseq", "sstrand")
  tbl$qstart <- as.numeric(tbl$qstart)
  tbl$qend <- as.numeric(tbl$qend)
  tbl$sstart <- as.numeric(tbl$sstart)
  tbl$send <- as.numeric(tbl$send)
  tbl$pident <- as.numeric(tbl$pident)
  tbl$length <- as.numeric(tbl$length)
  tbl$mismatch <- as.numeric(tbl$mismatch)
  tbl$gaps <- as.numeric(tbl$gaps)
  tbl[order(tbl$qstart, tbl$qend), , drop = FALSE]
}

.dnmb_prophage_alignment_blocks <- function(aln_tbl) {
  if (!is.data.frame(aln_tbl) || !nrow(aln_tbl)) {
    return(data.frame())
  }
  blocks <- data.frame(
    start = pmin(aln_tbl$qstart, aln_tbl$qend),
    end = pmax(aln_tbl$qstart, aln_tbl$qend),
    stringsAsFactors = FALSE
  )
  blocks <- blocks[order(blocks$start, blocks$end), , drop = FALSE]
  merged <- list()
  cur_start <- blocks$start[[1]]
  cur_end <- blocks$end[[1]]
  if (nrow(blocks) > 1L) {
    for (i in 2:nrow(blocks)) {
      if (blocks$start[[i]] <= cur_end + 1L) {
        cur_end <- max(cur_end, blocks$end[[i]])
      } else {
        merged[[length(merged) + 1L]] <- data.frame(start = cur_start, end = cur_end)
        cur_start <- blocks$start[[i]]
        cur_end <- blocks$end[[i]]
      }
    }
  }
  merged[[length(merged) + 1L]] <- data.frame(start = cur_start, end = cur_end)
  dplyr::bind_rows(merged)
}

.dnmb_prophage_alignment_regions <- function(aln_tbl, gap_tol = 250L) {
  if (!is.data.frame(aln_tbl) || !nrow(aln_tbl)) {
    return(data.frame())
  }
  segs <- data.frame(
    start = pmin(aln_tbl$qstart, aln_tbl$qend),
    end = pmax(aln_tbl$qstart, aln_tbl$qend),
    length = suppressWarnings(as.numeric(aln_tbl$length)),
    pident = suppressWarnings(as.numeric(aln_tbl$pident)),
    stringsAsFactors = FALSE
  )
  segs <- segs[order(segs$start, segs$end), , drop = FALSE]
  merged <- list()
  cur_start <- segs$start[[1]]
  cur_end <- segs$end[[1]]
  cur_weight <- ifelse(is.na(segs$length[[1]]), max(1, cur_end - cur_start + 1), segs$length[[1]])
  cur_sum <- cur_weight * segs$pident[[1]]
  if (nrow(segs) > 1L) {
    for (i in 2:nrow(segs)) {
      next_start <- segs$start[[i]]
      next_end <- segs$end[[i]]
      next_weight <- ifelse(is.na(segs$length[[i]]), max(1, next_end - next_start + 1), segs$length[[i]])
      next_sum <- next_weight * segs$pident[[i]]
      if (next_start <= cur_end + gap_tol) {
        cur_end <- max(cur_end, next_end)
        cur_weight <- cur_weight + next_weight
        cur_sum <- cur_sum + next_sum
      } else {
        merged[[length(merged) + 1L]] <- data.frame(
          start = cur_start,
          end = cur_end,
          pident = cur_sum / cur_weight,
          stringsAsFactors = FALSE
        )
        cur_start <- next_start
        cur_end <- next_end
        cur_weight <- next_weight
        cur_sum <- next_sum
      }
    }
  }
  merged[[length(merged) + 1L]] <- data.frame(
    start = cur_start,
    end = cur_end,
    pident = cur_sum / cur_weight,
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(merged)
}

.dnmb_prophage_snv_positions <- function(aln_tbl) {
  if (!is.data.frame(aln_tbl) || !nrow(aln_tbl)) {
    return(data.frame())
  }
  out <- list()
  for (i in seq_len(nrow(aln_tbl))) {
    qpos <- aln_tbl$qstart[[i]]
    qseq <- strsplit(aln_tbl$qseq[[i]], "", fixed = TRUE)[[1]]
    sseq <- strsplit(aln_tbl$sseq[[i]], "", fixed = TRUE)[[1]]
    rows <- list()
    for (j in seq_along(qseq)) {
      q <- qseq[[j]]
      s <- sseq[[j]]
      if (q != "-") {
        current_qpos <- qpos
        qpos <- qpos + 1L
      } else {
        current_qpos <- qpos
      }
      if (q == "-" || s == "-") {
        next
      }
      if (!identical(q, s)) {
        rows[[length(rows) + 1L]] <- data.frame(position = current_qpos, stringsAsFactors = FALSE)
      }
    }
    if (length(rows)) {
      out[[length(out) + 1L]] <- dplyr::bind_rows(rows)
    }
  }
  dplyr::bind_rows(out)
}

.dnmb_prophage_shared_block_bands <- function(block_tbl, min_rows = 2L) {
  if (!is.data.frame(block_tbl) || !nrow(block_tbl)) {
    return(data.frame())
  }
  x_points <- sort(unique(c(block_tbl$start, block_tbl$end + 1)))
  if (length(x_points) < 2L) {
    return(data.frame())
  }
  segs <- vector("list", length(x_points) - 1L)
  for (i in seq_len(length(x_points) - 1L)) {
    seg_start <- x_points[[i]]
    seg_end <- x_points[[i + 1L]] - 1L
    covered_rows <- unique(block_tbl$row[block_tbl$start <= seg_end & block_tbl$end >= seg_start])
    covered_rows <- covered_rows[!is.na(covered_rows)]
    if (length(covered_rows) >= min_rows) {
      segs[[i]] <- data.frame(
        start = seg_start,
        end = seg_end,
        shared_n = length(covered_rows),
        stringsAsFactors = FALSE
      )
    }
  }
  seg_tbl <- dplyr::bind_rows(segs)
  if (!nrow(seg_tbl)) {
    return(seg_tbl)
  }
  seg_tbl <- seg_tbl[order(seg_tbl$start, seg_tbl$end), , drop = FALSE]
  merged <- list()
  cur_start <- seg_tbl$start[[1]]
  cur_end <- seg_tbl$end[[1]]
  cur_n <- seg_tbl$shared_n[[1]]
  if (nrow(seg_tbl) > 1L) {
    for (i in 2:nrow(seg_tbl)) {
      if (seg_tbl$start[[i]] <= cur_end + 1L && seg_tbl$shared_n[[i]] == cur_n) {
        cur_end <- max(cur_end, seg_tbl$end[[i]])
      } else {
        merged[[length(merged) + 1L]] <- data.frame(start = cur_start, end = cur_end, shared_n = cur_n)
        cur_start <- seg_tbl$start[[i]]
        cur_end <- seg_tbl$end[[i]]
        cur_n <- seg_tbl$shared_n[[i]]
      }
    }
  }
  merged[[length(merged) + 1L]] <- data.frame(start = cur_start, end = cur_end, shared_n = cur_n)
  out <- dplyr::bind_rows(merged)
  if (nrow(out)) {
    out$label <- paste0(out$shared_n, " refs")
  }
  out
}

.dnmb_prophage_reference_display_name <- function(title) {
  title <- trimws(as.character(title)[1])
  if (is.na(title) || !nzchar(title)) {
    return(NA_character_)
  }
  title <- sub(",?\\s*complete genome\\.?$", "", title, ignore.case = TRUE)
  title <- sub(",?\\s*DNA\\.?$", "", title, ignore.case = TRUE)
  title <- gsub("\\s+", " ", title)
  trimws(title)
}

.dnmb_prophage_reference_short_name <- function(title) {
  title <- .dnmb_prophage_reference_display_name(title)
  if (is.na(title) || !nzchar(title)) {
    return(NA_character_)
  }
  title <- sub("^Geobacillus virus\\s+", "", title, ignore.case = TRUE)
  title <- sub("^Geobacillus phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^Bacillus phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^Phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^Deep-sea thermophilic phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^phi\\s+", "", title, ignore.case = TRUE)
  title <- sub("^phage\\s+", "", title, ignore.case = TRUE)
  trimws(title)
}

.dnmb_prophage_alignment_identity <- function(aln_tbl) {
  if (!is.data.frame(aln_tbl) || !nrow(aln_tbl)) {
    return(NA_real_)
  }
  weights <- suppressWarnings(as.numeric(aln_tbl$length))
  weights[is.na(weights) | weights <= 0] <- 1
  stats::weighted.mean(suppressWarnings(as.numeric(aln_tbl$pident)), w = weights, na.rm = TRUE)
}

.dnmb_prophage_alignment_coverage <- function(region_tbl, query_len) {
  if (!is.data.frame(region_tbl) || !nrow(region_tbl) || is.na(query_len) || query_len <= 0) {
    return(NA_real_)
  }
  100 * sum(region_tbl$end - region_tbl$start + 1, na.rm = TRUE) / query_len
}

.dnmb_prophage_query_segments <- function(query_genes, query_len) {
  if (!is.data.frame(query_genes) || !nrow(query_genes) || is.na(query_len) || query_len <= 0) {
    return(data.frame())
  }
  genes <- query_genes
  genes$seg_start <- pmin(as.numeric(genes$plot_start), as.numeric(genes$plot_end))
  genes$seg_end <- pmax(as.numeric(genes$plot_start), as.numeric(genes$plot_end))
  genes <- genes[order(genes$seg_start, genes$seg_end), , drop = FALSE]
  segs <- list()
  cursor <- 1L
  seg_id <- 1L
  for (i in seq_len(nrow(genes))) {
    st <- max(1L, as.integer(genes$seg_start[[i]]))
    en <- min(as.integer(query_len), as.integer(genes$seg_end[[i]]))
    if (st > cursor) {
      segs[[length(segs) + 1L]] <- data.frame(
        segment_id = seg_id,
        segment_type = "intergenic",
        start = cursor,
        end = st - 1L,
        stringsAsFactors = FALSE
      )
      seg_id <- seg_id + 1L
    }
    segs[[length(segs) + 1L]] <- data.frame(
      segment_id = seg_id,
      segment_type = "gene",
      start = st,
      end = en,
      locus_tag = genes$locus_tag[[i]],
      stringsAsFactors = FALSE
    )
    seg_id <- seg_id + 1L
    cursor <- en + 1L
  }
  if (cursor <= query_len) {
    segs[[length(segs) + 1L]] <- data.frame(
      segment_id = seg_id,
      segment_type = "intergenic",
      start = cursor,
      end = query_len,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(segs)
}

.dnmb_reverse_complement_dna <- function(seq) {
  seq <- toupper(gsub("[^ACGT]", "N", as.character(seq)[1]))
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  paste(rev(chartr("ACGTN", "TGCAN", chars)), collapse = "")
}

.dnmb_read_fasta_records <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) {
    return(list())
  }
  idx <- grep("^>", lines)
  if (!length(idx)) {
    return(list())
  }
  ends <- c(idx[-1] - 1L, length(lines))
  out <- list()
  for (i in seq_along(idx)) {
    header <- sub("^>", "", lines[[idx[[i]]]])
    seq <- paste(lines[(idx[[i]] + 1L):ends[[i]]], collapse = "")
    out[[header]] <- seq
  }
  out
}

.dnmb_translate_codon <- local({
  tbl <- c(
    TTT="F",TTC="F",TTA="L",TTG="L",TCT="S",TCC="S",TCA="S",TCG="S",
    TAT="Y",TAC="Y",TAA="*",TAG="*",TGT="C",TGC="C",TGA="*",TGG="W",
    CTT="L",CTC="L",CTA="L",CTG="L",CCT="P",CCC="P",CCA="P",CCG="P",
    CAT="H",CAC="H",CAA="Q",CAG="Q",CGT="R",CGC="R",CGA="R",CGG="R",
    ATT="I",ATC="I",ATA="I",ATG="M",ACT="T",ACC="T",ACA="T",ACG="T",
    AAT="N",AAC="N",AAA="K",AAG="K",AGT="S",AGC="S",AGA="R",AGG="R",
    GTT="V",GTC="V",GTA="V",GTG="V",GCT="A",GCC="A",GCA="A",GCG="A",
    GAT="D",GAC="D",GAA="E",GAG="E",GGT="G",GGC="G",GGA="G",GGG="G"
  )
  function(codon) {
    codon <- toupper(gsub("[^ACGT]", "N", as.character(codon)[1]))
    if (nchar(codon) != 3L || grepl("N", codon, fixed = TRUE)) {
      return(NA_character_)
    }
    unname(tbl[[codon]])
  }
})

.dnmb_prophage_reference_interval_from_blast <- function(ref_sequence, aln_tbl) {
  if (!nzchar(ref_sequence) || !is.data.frame(aln_tbl) || !nrow(aln_tbl)) {
    return(list(sequence = "", ref_start = NA_real_, ref_end = NA_real_, strand = "+"))
  }
  ref_start <- min(pmin(aln_tbl$sstart, aln_tbl$send), na.rm = TRUE)
  ref_end <- max(pmax(aln_tbl$sstart, aln_tbl$send), na.rm = TRUE)
  strand_votes <- ifelse(aln_tbl$send >= aln_tbl$sstart, "+", "-")
  strand <- names(sort(table(strand_votes), decreasing = TRUE))[1]
  subseq <- substr(ref_sequence, ref_start, ref_end)
  if (identical(strand, "-")) {
    subseq <- .dnmb_reverse_complement_dna(subseq)
  }
  list(sequence = subseq, ref_start = ref_start, ref_end = ref_end, strand = strand)
}

.dnmb_prophage_mafft_pairwise_alignment <- function(query_seq, ref_seq, cache_prefix) {
  aln_fa <- paste0(cache_prefix, "_mafft.fa")
  in_fa <- paste0(cache_prefix, "_mafft_input.fa")
  if (!file.exists(aln_fa) || file.info(aln_fa)$size == 0) {
    .dnmb_write_single_fasta("query", query_seq, in_fa)
    cat("\n", file = in_fa, append = TRUE)
    cat(">reference\n", ref_seq, "\n", file = in_fa, append = TRUE, sep = "")
    run <- dnmb_run_external("mafft", args = c("--quiet", "--auto", in_fa), required = FALSE)
    if (!isTRUE(run$ok) || is.null(run$stdout) || !nzchar(run$stdout)) {
      return(NULL)
    }
    writeLines(run$stdout, aln_fa)
  }
  recs <- .dnmb_read_fasta_records(aln_fa)
  if (!all(c("query", "reference") %in% names(recs))) {
    return(NULL)
  }
  list(query = recs[["query"]], reference = recs[["reference"]], path = aln_fa)
}

.dnmb_prophage_mutation_profile <- function(query_seq, query_genes, query_aln, ref_aln, gap_tol = 150L, window_bp = 200L) {
  q_chars <- strsplit(query_aln, "", fixed = TRUE)[[1]]
  r_chars <- strsplit(ref_aln, "", fixed = TRUE)[[1]]
  q_len <- nchar(query_seq)
  q_base <- rep(NA_character_, q_len)
  r_base <- rep(NA_character_, q_len)
  ref_pos <- rep(NA_real_, q_len)
  qpos <- 0L
  rpos <- 0L
  for (i in seq_along(q_chars)) {
    if (!identical(r_chars[[i]], "-")) {
      rpos <- rpos + 1L
    }
    if (!identical(q_chars[[i]], "-")) {
      qpos <- qpos + 1L
      q_base[[qpos]] <- toupper(q_chars[[i]])
      r_base[[qpos]] <- toupper(r_chars[[i]])
      ref_pos[[qpos]] <- if (!identical(r_chars[[i]], "-")) rpos else NA_real_
    }
  }
  aligned <- !is.na(r_base) & r_base != "-"
  identical_base <- aligned & q_base == r_base
  win_bp <- max(50L, as.integer(window_bp)[1])
  win_starts <- seq.int(1L, q_len, by = win_bp)
  window_rows <- lapply(win_starts, function(st) {
    en <- min(q_len, st + win_bp - 1L)
    idx <- st:en
    aligned_n <- sum(aligned[idx], na.rm = TRUE)
    if (!aligned_n) {
      return(NULL)
    }
    data.frame(
      start = st,
      end = en,
      pident = 100 * mean(identical_base[idx][aligned[idx]], na.rm = TRUE),
      aligned_frac = 100 * mean(aligned[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  window_tbl <- dplyr::bind_rows(window_rows)
  segs <- list()
  i <- 1L
  while (i <= q_len) {
    if (!aligned[[i]]) {
      i <- i + 1L
      next
    }
    start <- i
    end <- i
    while (end < q_len && (aligned[[end + 1L]] || (end + 1L - end <= gap_tol && !aligned[[end + 1L]] && end + 1L < q_len && aligned[[min(q_len, end + 2L)]]))) {
      if (aligned[[end + 1L]]) {
        end <- end + 1L
      } else {
        break
      }
    }
    idx <- start:end
    segs[[length(segs) + 1L]] <- data.frame(
      start = start,
      end = end,
      pident = 100 * mean(identical_base[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    i <- end + 1L
  }
  region_tbl <- dplyr::bind_rows(segs)
  gene_idx <- rep(NA_integer_, q_len)
  if (nrow(query_genes)) {
    for (gi in seq_len(nrow(query_genes))) {
      st <- min(query_genes$plot_start[[gi]], query_genes$plot_end[[gi]])
      en <- max(query_genes$plot_start[[gi]], query_genes$plot_end[[gi]])
      st <- max(1L, as.integer(st)); en <- min(q_len, as.integer(en))
      gene_idx[st:en] <- gi
    }
  }
  muts <- list()
  mismatch_pos <- which(aligned & !identical_base)
  for (pos in mismatch_pos) {
    cls <- "Noncoding"
    gi <- gene_idx[[pos]]
    if (!is.na(gi) && gi > 0L) {
      gene <- query_genes[gi, , drop = FALSE]
      gst <- min(gene$plot_start[[1]], gene$plot_end[[1]])
      gen <- max(gene$plot_start[[1]], gene$plot_end[[1]])
      dir <- as.character(gene$direction[[1]])
      if (identical(dir, "-")) {
        offset <- gen - pos
        codon_end <- gen - 3L * floor(offset / 3L)
        codon_pos <- (codon_end - 2L):codon_end
      } else {
        offset <- pos - gst
        codon_start <- gst + 3L * floor(offset / 3L)
        codon_pos <- codon_start:(codon_start + 2L)
      }
      if (all(codon_pos >= 1L & codon_pos <= q_len)) {
        qcodon_bases <- substring(query_seq, codon_pos, codon_pos)
        rcodon_bases <- qcodon_bases
        ref_codon <- r_base[codon_pos]
        if (all(!is.na(ref_codon)) && all(ref_codon != "-") && all(grepl("^[ACGT]$", ref_codon))) {
          rcodon_bases <- ref_codon
          qcodon <- paste(qcodon_bases, collapse = "")
          rcodon <- paste(rcodon_bases, collapse = "")
          if (identical(dir, "-")) {
            qcodon <- .dnmb_reverse_complement_dna(qcodon)
            rcodon <- .dnmb_reverse_complement_dna(rcodon)
          }
          aa_q <- .dnmb_translate_codon(qcodon)
          aa_r <- .dnmb_translate_codon(rcodon)
          if (!is.na(aa_q) && !is.na(aa_r)) {
            cls <- if (identical(aa_q, aa_r)) "Synonymous" else "Nonsynonymous"
          } else {
            cls <- "Nonsynonymous"
          }
        } else {
          cls <- "Nonsynonymous"
        }
      }
    }
    muts[[length(muts) + 1L]] <- data.frame(position = pos, class = cls, stringsAsFactors = FALSE)
  }
  list(
    regions = region_tbl,
    windows = window_tbl,
    mutations = dplyr::bind_rows(muts),
    identity = 100 * mean(identical_base[aligned], na.rm = TRUE),
    coverage = 100 * mean(aligned, na.rm = TRUE),
    aligned = aligned,
    identical = identical_base,
    ref_pos = ref_pos
  )
}

.dnmb_prophage_segment_identity_profile <- function(segment_tbl, aligned, identical_base, ref_pos) {
  if (!is.data.frame(segment_tbl) || !nrow(segment_tbl)) {
    return(data.frame())
  }
  out <- lapply(seq_len(nrow(segment_tbl)), function(i) {
    st <- max(1L, as.integer(segment_tbl$start[[i]]))
    en <- max(st, as.integer(segment_tbl$end[[i]]))
    idx <- st:en
    aligned_frac <- 100 * mean(aligned[idx], na.rm = TRUE)
    pident <- if (sum(aligned[idx], na.rm = TRUE) > 0) {
      100 * mean(identical_base[idx][aligned[idx]], na.rm = TRUE)
    } else {
      NA_real_
    }
    mapped_ref <- ref_pos[idx]
    mapped_ref <- mapped_ref[is.finite(mapped_ref)]
    data.frame(
      segment_id = segment_tbl$segment_id[[i]],
      segment_type = segment_tbl$segment_type[[i]],
      query_start = st,
      query_end = en,
      start = if (length(mapped_ref)) min(mapped_ref) else NA_real_,
      end = if (length(mapped_ref)) max(mapped_ref) else NA_real_,
      query_span = en - st + 1,
      aligned_frac = aligned_frac,
      pident = pident,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

.dnmb_prophage_merge_synteny_blocks <- function(segment_tbl, gap_tol = 400L, identity_tol = 18) {
  if (!is.data.frame(segment_tbl) || !nrow(segment_tbl)) {
    return(data.frame())
  }
  tbl <- segment_tbl[order(segment_tbl$row, segment_tbl$start, segment_tbl$end), , drop = FALSE]
  out <- list()
  for (row_id in unique(tbl$row)) {
    sub_tbl <- tbl[tbl$row == row_id, , drop = FALSE]
    keep <- !is.na(sub_tbl$pident) & !is.na(sub_tbl$aligned_frac) & sub_tbl$aligned_frac > 0
    sub_tbl <- sub_tbl[keep, , drop = FALSE]
    if (!nrow(sub_tbl)) next
    cur_start <- sub_tbl$start[[1]]
    cur_end <- sub_tbl$end[[1]]
    cur_pid <- sub_tbl$pident[[1]]
    cur_cov <- sub_tbl$aligned_frac[[1]]
    cur_type <- sub_tbl$segment_type[[1]]
    seg_ids <- c(sub_tbl$segment_id[[1]])
    if (nrow(sub_tbl) > 1L) {
      for (i in 2:nrow(sub_tbl)) {
        same_type <- identical(cur_type, sub_tbl$segment_type[[i]])
        close_gap <- sub_tbl$start[[i]] <= cur_end + gap_tol
        similar_identity <- abs(sub_tbl$pident[[i]] - cur_pid) <= identity_tol
        if (same_type && close_gap && similar_identity) {
          cur_end <- max(cur_end, sub_tbl$end[[i]])
          cur_pid <- mean(c(cur_pid, sub_tbl$pident[[i]]), na.rm = TRUE)
          cur_cov <- mean(c(cur_cov, sub_tbl$aligned_frac[[i]]), na.rm = TRUE)
          seg_ids <- c(seg_ids, sub_tbl$segment_id[[i]])
        } else {
          out[[length(out) + 1L]] <- data.frame(
            row = row_id,
            start = cur_start,
            end = cur_end,
            pident = cur_pid,
            aligned_frac = cur_cov,
            segment_type = cur_type,
            segment_ids = paste(seg_ids, collapse = ","),
            stringsAsFactors = FALSE
          )
          cur_start <- sub_tbl$start[[i]]
          cur_end <- sub_tbl$end[[i]]
          cur_pid <- sub_tbl$pident[[i]]
          cur_cov <- sub_tbl$aligned_frac[[i]]
          cur_type <- sub_tbl$segment_type[[i]]
          seg_ids <- c(sub_tbl$segment_id[[i]])
        }
      }
    }
    out[[length(out) + 1L]] <- data.frame(
      row = row_id,
      start = cur_start,
      end = cur_end,
      pident = cur_pid,
      aligned_frac = cur_cov,
      segment_type = cur_type,
      segment_ids = paste(seg_ids, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(out)
}

.dnmb_prophage_chain_rows <- function(segment_tbl, row_identity_map) {
  rows <- names(row_identity_map)
  if (length(rows) <= 1L) {
    return(rows)
  }
  remaining <- rows
  first <- rows[[which.max(as.numeric(row_identity_map[rows]))]]
  order_out <- first
  remaining <- setdiff(remaining, first)
  while (length(remaining)) {
    prev <- order_out[[length(order_out)]]
    prev_tbl <- segment_tbl[segment_tbl$row == prev, , drop = FALSE]
    scores <- vapply(remaining, function(candidate) {
      cand_tbl <- segment_tbl[segment_tbl$row == candidate, , drop = FALSE]
      joined <- dplyr::inner_join(
        prev_tbl[, c("segment_id", "start", "end", "pident", "aligned_frac")],
        cand_tbl[, c("segment_id", "start", "end", "pident", "aligned_frac")],
        by = c("segment_id", "start", "end"),
        suffix = c("_prev", "_cand")
      )
      if (!nrow(joined)) {
        return(-Inf)
      }
      overlap_len <- joined$end - joined$start + 1
      pair_id <- pmin(joined$pident_prev, joined$pident_cand, na.rm = TRUE)
      pair_cov <- pmin(joined$aligned_frac_prev, joined$aligned_frac_cand, na.rm = TRUE)
      sum(overlap_len * pmax(0, pair_id - 40) * (pair_cov / 100), na.rm = TRUE)
    }, numeric(1))
    if (all(!is.finite(scores))) {
      next_row <- remaining[[which.max(as.numeric(row_identity_map[remaining]))]]
    } else {
      next_row <- remaining[[which.max(scores)]]
    }
    order_out <- c(order_out, next_row)
    remaining <- setdiff(remaining, next_row)
  }
  order_out
}

.dnmb_plot_prophage_reference_panel <- function(sub_tbl, region_id, output_dir, top_n = 5L, cache_root = NULL) {
  text_family <- "sans"
  gbff_path <- .dnmb_find_gbff_for_plot(output_dir)
  if (is.null(gbff_path) || !file.exists(gbff_path)) {
    return(NULL)
  }
  records <- .dnmb_prophage_parse_genbank_records(gbff_path)
  if (!nrow(records)) {
    return(NULL)
  }
  contig_target <- .dnmb_normalize_contig_label(sub_tbl$contig[[1]])
  contig_match <- .dnmb_normalize_contig_label(records$definition) == contig_target |
    .dnmb_normalize_contig_label(records$accession) == contig_target |
    .dnmb_normalize_contig_label(records$locus) == contig_target
  if (!any(contig_match)) {
    return(NULL)
  }
  query_sequence <- records$sequence[[which(contig_match)[1]]]
  if (is.na(query_sequence) || !nzchar(query_sequence)) {
    return(NULL)
  }
  region_start <- min(pmin(as.numeric(sub_tbl$start), as.numeric(sub_tbl$end)), na.rm = TRUE)
  region_end <- max(pmax(as.numeric(sub_tbl$start), as.numeric(sub_tbl$end)), na.rm = TRUE)
  if (!is.finite(region_start) || !is.finite(region_end) || region_end <= region_start) {
    return(NULL)
  }
  query_len <- region_end - region_start + 1
  query_region_seq <- substr(query_sequence, region_start, region_end)
  cache_dir <- .dnmb_prophage_reference_cache_dir(output_dir, region_id)
  query_fasta <- file.path(cache_dir, "query_region.fna")
  .dnmb_write_single_fasta(paste0("Prophage_", region_id), query_region_seq, query_fasta)
  hit_tbl <- .dnmb_prophage_remote_reference_hits(query_fasta, cache_dir = cache_dir, top_n = top_n, cache_root = cache_root)
  if (!nrow(hit_tbl)) {
    return(NULL)
  }

  query_genes <- sub_tbl
  query_genes$plot_start <- pmin(as.numeric(query_genes$start), as.numeric(query_genes$end)) - region_start + 1
  query_genes$plot_end <- pmax(as.numeric(query_genes$start), as.numeric(query_genes$end)) - region_start + 1
  query_genes$category <- .dnmb_gene_arrow_category_prophage(query_genes$product)
  query_segments <- .dnmb_prophage_query_segments(query_genes, query_len)
  palette <- .dnmb_gene_arrow_palette_prophage()
  seg_fill_fun <- scales::col_numeric(palette = c("#D8E5FF", "#7FB1FF", "#1D4ED8"), domain = c(40, 95))

  ref_seg_rows <- list()
  ref_mut_rows <- list()
  ref_label_rows <- list()
  row_identity_map <- numeric()
  row_coverage_map <- numeric()

  for (i in seq_len(nrow(hit_tbl))) {
    acc <- hit_tbl$accession[[i]]
    gbk_path <- .dnmb_fetch_genbank_text(acc, file.path(cache_dir, paste0(acc, ".gbk")))
    if (!nzchar(gbk_path) || !file.exists(gbk_path)) next
    ref <- .dnmb_parse_single_genbank_record(gbk_path)
    if (!nzchar(ref$sequence)) next
    ref_fasta <- file.path(cache_dir, paste0(acc, ".fna"))
    .dnmb_write_single_fasta(acc, ref$sequence, ref_fasta)
    aln_tbl <- .dnmb_prophage_pairwise_alignment(query_fasta, ref_fasta, file.path(cache_dir, paste0(acc, "_pairwise.tsv")))
    if (!nrow(aln_tbl)) next
    ref_interval <- .dnmb_prophage_reference_interval_from_blast(ref$sequence, aln_tbl)
    if (!nzchar(ref_interval$sequence)) next
    mafft_aln <- .dnmb_prophage_mafft_pairwise_alignment(query_region_seq, ref_interval$sequence, file.path(cache_dir, paste0(acc, "_stable")))
    if (is.null(mafft_aln)) next
    profile <- .dnmb_prophage_mutation_profile(
      query_seq = query_region_seq,
      query_genes = query_genes,
      query_aln = mafft_aln$query,
      ref_aln = mafft_aln$reference,
      window_bp = max(150L, floor(query_len / 90L))
    )
    seg_profile <- .dnmb_prophage_segment_identity_profile(
      segment_tbl = query_segments,
      aligned = profile$aligned,
      identical_base = profile$identical,
      ref_pos = profile$ref_pos
    )
    if (!nrow(seg_profile)) next
    # Clamp reference coordinates to query_len range
    seg_profile$start <- pmax(1, pmin(query_len, seg_profile$start))
    seg_profile$end   <- pmax(1, pmin(query_len, seg_profile$end))
    seg_profile$row <- acc
    seg_profile$fill_hex <- ifelse(
      !is.na(seg_profile$pident) & !is.na(seg_profile$aligned_frac) & seg_profile$pident >= 40 & seg_profile$aligned_frac >= 10,
      seg_fill_fun(pmin(95, pmax(40, seg_profile$pident))),
      "#FFFFFF"
    )
    seg_profile$seg_height <- ifelse(seg_profile$segment_type == "gene", 0.16, 0.07)
    ref_seg_rows[[length(ref_seg_rows) + 1L]] <- seg_profile
    if (nrow(profile$mutations)) {
      muts <- profile$mutations
      muts$row <- acc
      ref_mut_rows[[length(ref_mut_rows) + 1L]] <- muts
    }
    ref_label_rows[[length(ref_label_rows) + 1L]] <- data.frame(
      row = acc,
      label = paste0(acc, " | ", round(profile$identity), "/", round(profile$coverage), "\n", .dnmb_prophage_reference_short_name(hit_tbl$title[[i]])),
      stringsAsFactors = FALSE
    )
    row_identity_map[[acc]] <- profile$identity
    row_coverage_map[[acc]] <- profile$coverage
  }

  segment_tbl <- dplyr::bind_rows(ref_seg_rows)
  snv_tbl <- dplyr::bind_rows(ref_mut_rows)
  label_tbl <- dplyr::bind_rows(ref_label_rows)
  if (!nrow(segment_tbl)) {
    return(NULL)
  }
  # Sort by composite score (coverage * identity) so best overall match is on top
  row_score_map <- row_identity_map * row_coverage_map / 100
  row_levels <- c("Query prophage", names(sort(row_score_map, decreasing = TRUE)))
  row_spacing <- 1.55
  row_map <- stats::setNames(rev(seq_along(row_levels)) * row_spacing, row_levels)
  segment_tbl$y <- row_map[segment_tbl$row]
  snv_tbl$y <- row_map[snv_tbl$row]
  label_tbl$y <- row_map[label_tbl$row]

  query_poly_list <- lapply(seq_len(nrow(query_genes)), function(i) {
    poly <- .dnmb_linear_gene_polygon(
      start = query_genes$plot_start[[i]],
      end = query_genes$plot_end[[i]],
      direction = query_genes$direction[[i]],
      xmin = 1,
      xmax = query_len,
      x_min_plot = 1,
      x_max_plot = query_len,
      y_center = row_map[["Query prophage"]],
      height = 0.46,
      arrow_frac = 0.18
    )
    if (!nrow(poly)) return(NULL)
    poly$feature_id <- i
    poly$fill_hex <- unname(palette[query_genes$category[[i]]])
    poly
  })
  query_poly_tbl <- dplyr::bind_rows(Filter(Negate(is.null), query_poly_list))

  # Build ribbon polygons connecting adjacent rows
  ribbon_tbl <- data.frame()
  ref_rows <- row_levels[row_levels != "Query prophage"]
  if (nrow(segment_tbl) && length(ref_rows)) {
    ribbon_list <- list()
    ribbon_idx <- 1L
    for (pair_i in seq_len(length(row_levels) - 1L)) {
      upper_row <- row_levels[[pair_i]]
      lower_row <- row_levels[[pair_i + 1L]]
      lower_segs <- segment_tbl[
        segment_tbl$row == lower_row &
          segment_tbl$segment_type == "gene" &
          segment_tbl$fill_hex != "#FFFFFF" &
          !is.na(segment_tbl$pident) &
          segment_tbl$pident >= 40 &
          segment_tbl$aligned_frac >= 10,
        , drop = FALSE
      ]
      if (!nrow(lower_segs)) next
      upper_y <- unname(row_map[[upper_row]])
      lower_y <- unname(row_map[[lower_row]])
      top_y <- upper_y - if (identical(upper_row, "Query prophage")) 0.24 else 0.18
      bottom_y <- lower_y + 0.18
      # Get upper row segments for ref-to-ref ribbons
      upper_segs <- if (!identical(upper_row, "Query prophage")) {
        segment_tbl[segment_tbl$row == upper_row & segment_tbl$segment_type == "gene", , drop = FALSE]
      } else {
        NULL
      }
      for (j in seq_len(nrow(lower_segs))) {
        # Top coordinates: query position (for Query row) or ref position (for ref rows)
        qs <- lower_segs$query_start[[j]]
        qe <- lower_segs$query_end[[j]]
        if (identical(upper_row, "Query prophage")) {
          top_left <- qs
          top_right <- qe
        } else if (!is.null(upper_segs) && nrow(upper_segs)) {
          # Find matching upper segment by segment_id
          match_idx <- which(upper_segs$segment_id == lower_segs$segment_id[[j]])
          if (length(match_idx)) {
            us <- upper_segs[match_idx[1], ]
            top_left <- if (is.finite(us$start)) us$start else qs
            top_right <- if (is.finite(us$end)) us$end else qe
          } else {
            top_left <- qs
            top_right <- qe
          }
        } else {
          top_left <- qs
          top_right <- qe
        }
        # Bottom coordinates: reference position from lower row (clamped to query range)
        rs <- lower_segs$start[[j]]
        re <- lower_segs$end[[j]]
        bottom_left <- pmax(1, pmin(query_len, if (is.finite(rs)) rs else qs))
        bottom_right <- pmax(1, pmin(query_len, if (is.finite(re)) re else qe))
        top_w <- max(1, top_right - top_left + 1)
        bottom_w <- max(1, bottom_right - bottom_left + 1)
        top_taper <- min(top_w * 0.10, max(6, top_w * 0.16))
        bottom_taper <- min(bottom_w * 0.10, max(6, bottom_w * 0.16))
        pident_val <- lower_segs$pident[[j]]
        alpha_val <- pmin(0.72, pmax(0.15, scales::rescale(pident_val, to = c(0.15, 0.72), from = c(40, 100))))
        seg_mid <- (qs + qe) / 2
        gene_idx <- which(query_genes$plot_start <= seg_mid & query_genes$plot_end >= seg_mid)
        seg_category <- if (length(gene_idx)) as.character(query_genes$category[[gene_idx[1]]]) else "Other"
        seg_color <- unname(palette[seg_category])
        if (is.na(seg_color)) seg_color <- "#9270CA"
        blend_frac <- scales::rescale(pident_val, to = c(0.5, 1.0), from = c(40, 100))
        seg_rgb <- grDevices::col2rgb(seg_color) / 255
        blended_rgb <- seg_rgb * blend_frac + 1.0 * (1 - blend_frac)
        seg_color <- grDevices::rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3])
        ribbon_list[[ribbon_idx]] <- data.frame(
          x = c(top_left + top_taper, top_right - top_taper, bottom_right - bottom_taper, bottom_left + bottom_taper),
          y = c(top_y, top_y, bottom_y, bottom_y),
          group = ribbon_idx,
          fill_hex = seg_color,
          alpha_val = alpha_val,
          pident = pident_val,
          stringsAsFactors = FALSE
        )
        ribbon_idx <- ribbon_idx + 1L
      }
    }
    if (length(ribbon_list)) {
      ribbon_tbl <- dplyr::bind_rows(ribbon_list)
    }
  }

  show_snv <- nrow(snv_tbl) > 0L

  scale_len <- if (query_len >= 5000) 5000 else if (query_len >= 2000) 2000 else 1000
  legend_y_base <- -0.55
  scale_tbl <- data.frame(x = 1, xend = scale_len + 1, y = legend_y_base + 0.10, label = .dnmb_fmt_bp_label(scale_len))
  grad_fun <- scales::col_numeric(palette = c("#D8E5FF", "#7FB1FF", "#1D4ED8"), domain = c(40, 95))
  grad_vals <- seq(40, 95, length.out = 30L)
  grad_xmin <- query_len * 0.35
  grad_w <- query_len * 0.015
  identity_legend_tbl <- data.frame(
    xmin = grad_xmin + grad_w * (seq_along(grad_vals) - 1L),
    xmax = grad_xmin + grad_w * seq_along(grad_vals),
    ymin = legend_y_base,
    ymax = legend_y_base + 0.16,
    fill_hex = grad_fun(grad_vals)
  )
  identity_label_tbl <- data.frame(x = grad_xmin - query_len * 0.02, y = legend_y_base + 0.08, label = "Identity (%)", hjust = 1)
  identity_tick_tbl <- data.frame(
    x = c(identity_legend_tbl$xmin[[1]], identity_legend_tbl$xmax[[nrow(identity_legend_tbl)]]),
    y = c(legend_y_base - 0.08, legend_y_base - 0.08),
    label = c("40", "95")
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = data.frame(row = row_levels, y = row_map[row_levels]), ggplot2::aes(x = 1, xend = query_len, y = .data$y, yend = .data$y), linewidth = 4.2, color = "#E5E7EB", lineend = "round")
  if (nrow(ribbon_tbl)) {
    p <- p + ggplot2::geom_polygon(
      data = ribbon_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$fill_hex, alpha = .data$alpha_val),
      color = NA, inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::geom_polygon(data = query_poly_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_hex), color = "grey25", linewidth = 0.25) +
    ggplot2::geom_rect(data = segment_tbl, ggplot2::aes(xmin = .data$query_start, xmax = .data$query_end, ymin = .data$y - .data$seg_height, ymax = .data$y + .data$seg_height, fill = .data$fill_hex), color = NA, alpha = 0.9, inherit.aes = FALSE)
  if (show_snv && nrow(snv_tbl)) {
    p <- p + ggplot2::geom_segment(data = snv_tbl, ggplot2::aes(x = .data$position, xend = .data$position, y = .data$y - 0.22, yend = .data$y + 0.22, color = .data$class), linewidth = 0.36)
  }
  p <- p +
    ggplot2::geom_segment(data = scale_tbl, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y), linewidth = 0.5, color = "grey35") +
    ggplot2::geom_rect(data = identity_legend_tbl, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, fill = .data$fill_hex), color = NA, inherit.aes = FALSE) +
    ggplot2::geom_text(data = identity_label_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = identity_label_tbl$hjust[[1]], vjust = 0.5, size = 2.8, color = "#4B5563", family = text_family, inherit.aes = FALSE) +
    ggplot2::geom_text(data = identity_tick_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = c(0, 1), vjust = 0.5, size = 2.5, color = "#6B7280", family = text_family, inherit.aes = FALSE) +
    ggplot2::geom_text(data = scale_tbl, ggplot2::aes(x = (.data$x + .data$xend) / 2, y = .data$y - 0.18, label = .data$label), size = 3.2, vjust = 1, family = text_family) +
    ggplot2::geom_text(data = data.frame(label = "Query prophage", y = row_map[["Query prophage"]]), ggplot2::aes(x = query_len + query_len * 0.016, y = .data$y, label = .data$label), hjust = 0, vjust = 0.5, size = 3.0, family = text_family, inherit.aes = FALSE) +
    ggplot2::geom_text(data = label_tbl, ggplot2::aes(x = query_len + query_len * 0.016, y = .data$y, label = .data$label), hjust = 0, vjust = 0.5, size = 2.6, lineheight = 0.88, family = text_family, inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = c("Synonymous" = "#A5BFF9", "Nonsynonymous" = "#2563EB", "Noncoding" = "#D9468B"), guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(1, query_len), ylim = c(legend_y_base - 0.18, max(row_map) + 0.8), clip = "off") +
    ggplot2::labs(title = "Reference Phage Comparison", x = "Prophage coordinate (bp)", y = NULL) +
    ggplot2::theme_bw(base_size = 11, base_family = text_family) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, family = text_family),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 180, 14, 12)
    ) +
    ggplot2::geom_hline(yintercept = 0.6, linewidth = 0.3, color = "#D1D5DB")
  p
}
