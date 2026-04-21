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

.dnmb_prophage_capsule_polygon <- function(xmin, xmax, y_center, half_height, x_radius = NULL, n_cap = 48L) {
  xmin <- suppressWarnings(as.numeric(xmin)[1])
  xmax <- suppressWarnings(as.numeric(xmax)[1])
  y_center <- suppressWarnings(as.numeric(y_center)[1])
  half_height <- suppressWarnings(as.numeric(half_height)[1])
  x_radius <- suppressWarnings(as.numeric(x_radius)[1])
  if (!is.finite(xmin) || !is.finite(xmax) || !is.finite(y_center) || !is.finite(half_height) ||
      xmax <= xmin || half_height <= 0) {
    return(data.frame())
  }
  width <- xmax - xmin
  if (!is.finite(x_radius) || x_radius <= 0) {
    x_radius <- half_height
  }
  x_radius <- min(x_radius, width / 2)
  if (width <= 2 * x_radius + 1e-8) {
    theta <- seq(pi / 2, pi / 2 - 2 * pi, length.out = max(24L, 2L * n_cap))
    return(data.frame(
      x = (xmin + xmax) / 2 + (width / 2) * cos(theta),
      y = y_center + half_height * sin(theta),
      stringsAsFactors = FALSE
    ))
  }
  left_center <- xmin + x_radius
  right_center <- xmax - x_radius
  right_theta <- seq(pi / 2, -pi / 2, length.out = n_cap)
  left_theta <- seq(-pi / 2, pi / 2, length.out = n_cap)
  data.frame(
    x = c(
      left_center,
      right_center,
      right_center + x_radius * cos(right_theta),
      right_center,
      left_center,
      left_center - x_radius * cos(left_theta)
    ),
    y = c(
      y_center + half_height,
      y_center + half_height,
      y_center + half_height * sin(right_theta),
      y_center - half_height,
      y_center - half_height,
      y_center + half_height * sin(left_theta)
    ),
    stringsAsFactors = FALSE
  )
}

.dnmb_prophage_capsule_half_height_at_x <- function(x, lane_xmin, lane_xmax, half_height, x_radius = NULL) {
  x <- suppressWarnings(as.numeric(x))
  lane_xmin <- suppressWarnings(as.numeric(lane_xmin)[1])
  lane_xmax <- suppressWarnings(as.numeric(lane_xmax)[1])
  half_height <- suppressWarnings(as.numeric(half_height)[1])
  x_radius <- suppressWarnings(as.numeric(x_radius)[1])
  if (!length(x) || !is.finite(lane_xmin) || !is.finite(lane_xmax) || !is.finite(half_height) ||
      lane_xmax <= lane_xmin || half_height <= 0) {
    return(rep(NA_real_, length(x)))
  }
  width <- lane_xmax - lane_xmin
  if (!is.finite(x_radius) || x_radius <= 0) {
    x_radius <- half_height
  }
  x_radius <- min(x_radius, width / 2)
  left_center <- lane_xmin + x_radius
  right_center <- lane_xmax - x_radius
  out <- rep(half_height, length(x))
  left_idx <- which(x < left_center)
  right_idx <- which(x > right_center)
  if (length(left_idx)) {
    dx <- pmax(-x_radius, pmin(x_radius, x[left_idx] - left_center))
    out[left_idx] <- half_height * sqrt(pmax(0, 1 - (dx / x_radius)^2))
  }
  if (length(right_idx)) {
    dx <- pmax(-x_radius, pmin(x_radius, x[right_idx] - right_center))
    out[right_idx] <- half_height * sqrt(pmax(0, 1 - (dx / x_radius)^2))
  }
  out
}

.dnmb_prophage_clipped_bar_polygon <- function(xmin, xmax, y_center, bar_half, lane_xmin, lane_xmax, lane_half, lane_x_radius = NULL, n = 40L) {
  xmin <- suppressWarnings(as.numeric(xmin)[1])
  xmax <- suppressWarnings(as.numeric(xmax)[1])
  y_center <- suppressWarnings(as.numeric(y_center)[1])
  bar_half <- suppressWarnings(as.numeric(bar_half)[1])
  lane_xmin <- suppressWarnings(as.numeric(lane_xmin)[1])
  lane_xmax <- suppressWarnings(as.numeric(lane_xmax)[1])
  lane_half <- suppressWarnings(as.numeric(lane_half)[1])
  lane_x_radius <- suppressWarnings(as.numeric(lane_x_radius)[1])
  if (!is.finite(xmin) || !is.finite(xmax) || !is.finite(y_center) || !is.finite(bar_half) ||
      !is.finite(lane_xmin) || !is.finite(lane_xmax) || !is.finite(lane_half) ||
      xmax <= xmin || lane_xmax <= lane_xmin || bar_half <= 0 || lane_half <= 0) {
    return(data.frame())
  }
  xseq <- seq(max(xmin, lane_xmin), min(xmax, lane_xmax), length.out = max(8L, n))
  if (!length(xseq)) {
    return(data.frame())
  }
  allowed_half <- .dnmb_prophage_capsule_half_height_at_x(
    xseq,
    lane_xmin = lane_xmin,
    lane_xmax = lane_xmax,
    half_height = lane_half,
    x_radius = lane_x_radius
  )
  half_vec <- pmin(bar_half, allowed_half)
  data.frame(
    x = c(xseq, rev(xseq)),
    y = c(y_center + half_vec, rev(y_center - half_vec)),
    stringsAsFactors = FALSE
  )
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
  title <- sub(",?\\s*complete sequence\\.?$", "", title, ignore.case = TRUE)
  title <- sub(",?\\s*DNA\\.?$", "", title, ignore.case = TRUE)
  title <- gsub("\\s+", " ", title)
  trimws(title)
}

.dnmb_prophage_compact_genome_label <- function(title) {
  title <- trimws(as.character(title)[1])
  if (is.na(title) || !nzchar(title)) {
    return(NA_character_)
  }
  title <- sub(",?\\s*complete genome\\.?$", "", title, ignore.case = TRUE)
  title <- sub(",?\\s*complete sequence\\.?$", "", title, ignore.case = TRUE)
  title <- sub(",?\\s*chromosome\\b", "", title, ignore.case = TRUE)
  title <- sub(",?\\s*complete\\b", "", title, ignore.case = TRUE)
  title <- gsub("\\s+,", ",", title)
  title <- gsub(",\\s*,", ", ", title)
  title <- gsub("\\s+", " ", title)
  title <- gsub(",\\s*$", "", title)
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
  title <- sub("^Bacteriophage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^Phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^Deep-sea thermophilic phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("^phi\\s+", "", title, ignore.case = TRUE)
  title <- sub("^phage\\s+", "", title, ignore.case = TRUE)
  title <- sub("\\s+DNA$", "", title, ignore.case = TRUE)
  trimws(title)
}

.dnmb_prophage_reference_cluster_key <- function(title) {
  title <- .dnmb_prophage_reference_short_name(title)
  if (is.na(title) || !nzchar(title)) {
    return(NA_character_)
  }
  title <- tolower(title)
  title <- gsub("[^a-z0-9]+", "", title)
  trimws(title)
}

.dnmb_prophage_query_overlap_fraction <- function(tbl_a, tbl_b) {
  if (!is.data.frame(tbl_a) || !nrow(tbl_a) || !is.data.frame(tbl_b) || !nrow(tbl_b)) {
    return(0)
  }
  a <- tbl_a[, c("query_start", "query_end"), drop = FALSE]
  b <- tbl_b[, c("query_start", "query_end"), drop = FALSE]
  a$width <- a$query_end - a$query_start + 1
  b$width <- b$query_end - b$query_start + 1
  overlap <- 0
  for (i in seq_len(nrow(a))) {
    for (j in seq_len(nrow(b))) {
      st <- max(a$query_start[[i]], b$query_start[[j]])
      en <- min(a$query_end[[i]], b$query_end[[j]])
      if (is.finite(st) && is.finite(en) && en >= st) {
        overlap <- overlap + (en - st + 1)
      }
    }
  }
  denom <- min(sum(a$width, na.rm = TRUE), sum(b$width, na.rm = TRUE))
  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }
  overlap / denom
}

.dnmb_prophage_prune_reference_candidates <- function(segment_tbl, candidate_meta, max_refs = 5L) {
  if (!is.data.frame(candidate_meta) || !nrow(candidate_meta)) {
    return(character())
  }
  candidate_meta <- candidate_meta[order(-candidate_meta$score, -candidate_meta$coverage, -candidate_meta$identity), , drop = FALSE]
  n <- nrow(candidate_meta)
  if (n <= 1L) {
    return(candidate_meta$row[seq_len(min(n, max_refs))])
  }

  adjacency <- vector("list", n)
  for (i in seq_len(n)) {
    adjacency[[i]] <- i
  }

  for (i in seq_len(n - 1L)) {
    row_i <- candidate_meta$row[[i]]
    key_i <- candidate_meta$cluster_key[[i]]
    tbl_i <- segment_tbl[segment_tbl$row == row_i & segment_tbl$segment_type == "gene" & segment_tbl$fill_hex != "#FFFFFF", , drop = FALSE]
    for (j in seq.int(i + 1L, n)) {
      row_j <- candidate_meta$row[[j]]
      key_j <- candidate_meta$cluster_key[[j]]
      same_named_cluster <- !is.na(key_i) && nzchar(key_i) && !is.na(key_j) && nzchar(key_j) && identical(key_i, key_j)
      tbl_j <- segment_tbl[segment_tbl$row == row_j & segment_tbl$segment_type == "gene" & segment_tbl$fill_hex != "#FFFFFF", , drop = FALSE]
      ov_frac <- .dnmb_prophage_query_overlap_fraction(tbl_i, tbl_j)
      similar_id <- abs(candidate_meta$identity[[i]] - candidate_meta$identity[[j]]) <= 3
      similar_cov <- abs(candidate_meta$coverage[[i]] - candidate_meta$coverage[[j]]) <= 8
      similar_len <- abs(candidate_meta$row_len[[i]] - candidate_meta$row_len[[j]]) <= max(500, 0.05 * max(candidate_meta$row_len[[i]], candidate_meta$row_len[[j]]))
      same_alignment_cluster <- ov_frac >= 0.92 && similar_id && similar_cov && similar_len
      if (same_named_cluster || same_alignment_cluster) {
        adjacency[[i]] <- c(adjacency[[i]], j)
        adjacency[[j]] <- c(adjacency[[j]], i)
      }
    }
  }

  component_id <- rep(NA_integer_, n)
  component_counter <- 0L
  for (i in seq_len(n)) {
    if (!is.na(component_id[[i]])) next
    component_counter <- component_counter + 1L
    stack <- i
    while (length(stack)) {
      cur <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      if (!is.na(component_id[[cur]])) next
      component_id[[cur]] <- component_counter
      nbrs <- adjacency[[cur]]
      stack <- c(stack, nbrs[is.na(component_id[nbrs])])
    }
  }

  candidate_meta$component_id <- component_id
  reps <- lapply(split(candidate_meta, candidate_meta$component_id), function(df) {
    df <- df[order(-df$score, -df$coverage, -df$identity), , drop = FALSE]
    df[1, , drop = FALSE]
  })
  rep_meta <- dplyr::bind_rows(reps)
  rep_meta <- rep_meta[order(-rep_meta$score, -rep_meta$coverage, -rep_meta$identity), , drop = FALSE]
  head(rep_meta$row, max_refs)
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
    muts[[length(muts) + 1L]] <- data.frame(
      position = pos,
      query_position = pos,
      ref_position = ref_pos[[pos]],
      class = cls,
      stringsAsFactors = FALSE
    )
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
    mapped_query <- idx[is.finite(ref_pos[idx])]
    data.frame(
      segment_id = segment_tbl$segment_id[[i]],
      segment_type = segment_tbl$segment_type[[i]],
      query_start = st,
      query_end = en,
      start = if (length(mapped_ref)) min(mapped_ref) else NA_real_,
      end = if (length(mapped_ref)) max(mapped_ref) else NA_real_,
      projected_query_start = if (length(mapped_query)) min(mapped_query) else st,
      projected_query_end = if (length(mapped_query)) max(mapped_query) else en,
      query_span = en - st + 1,
      aligned_frac = aligned_frac,
      pident = pident,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

.dnmb_prophage_segment_profile_from_regions <- function(segment_tbl, region_tbl, query_len) {
  if (!is.data.frame(segment_tbl) || !nrow(segment_tbl) || !is.data.frame(region_tbl) || !nrow(region_tbl)) {
    return(data.frame())
  }
  out <- lapply(seq_len(nrow(segment_tbl)), function(i) {
    st <- max(1L, as.integer(segment_tbl$start[[i]]))
    en <- min(as.integer(query_len), max(st, as.integer(segment_tbl$end[[i]])))
    overlaps <- region_tbl[region_tbl$start <= en & region_tbl$end >= st, , drop = FALSE]
    if (!nrow(overlaps)) {
      return(NULL)
    }
    ov_start <- pmax(st, overlaps$start)
    ov_end <- pmin(en, overlaps$end)
    ov_width <- pmax(0, ov_end - ov_start + 1)
    keep <- ov_width > 0
    if (!any(keep)) {
      return(NULL)
    }
    ov_start <- ov_start[keep]
    ov_end <- ov_end[keep]
    ov_width <- ov_width[keep]
    overlaps <- overlaps[keep, , drop = FALSE]
    data.frame(
      segment_id = segment_tbl$segment_id[[i]],
      segment_type = segment_tbl$segment_type[[i]],
      query_start = st,
      query_end = en,
      start = min(ov_start),
      end = max(ov_end),
      query_span = en - st + 1,
      aligned_frac = 100 * sum(ov_width) / max(1, en - st + 1),
      pident = stats::weighted.mean(overlaps$pident, w = ov_width, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(Filter(Negate(is.null), out))
}

.dnmb_prophage_segment_profile_from_alignment <- function(segment_tbl, aln_tbl, query_len) {
  if (!is.data.frame(segment_tbl) || !nrow(segment_tbl) || !is.data.frame(aln_tbl) || !nrow(aln_tbl)) {
    return(data.frame())
  }
  aln <- aln_tbl
  aln$q_lo <- pmin(aln$qstart, aln$qend)
  aln$q_hi <- pmax(aln$qstart, aln$qend)
  out <- lapply(seq_len(nrow(segment_tbl)), function(i) {
    st <- max(1L, as.integer(segment_tbl$start[[i]]))
    en <- min(as.integer(query_len), max(st, as.integer(segment_tbl$end[[i]])))
    hits <- aln[aln$q_lo <= en & aln$q_hi >= st, , drop = FALSE]
    if (!nrow(hits)) {
      return(NULL)
    }
    pieces <- lapply(seq_len(nrow(hits)), function(j) {
      h <- hits[j, , drop = FALSE]
      ov_qs <- max(st, h$q_lo[[1]])
      ov_qe <- min(en, h$q_hi[[1]])
      if (!is.finite(ov_qs) || !is.finite(ov_qe) || ov_qe < ov_qs) {
        return(NULL)
      }
      q_dir <- if (h$qend[[1]] >= h$qstart[[1]]) 1 else -1
      s_dir <- if (h$send[[1]] >= h$sstart[[1]]) 1 else -1
      offset_start <- if (q_dir == 1) ov_qs - h$qstart[[1]] else h$qstart[[1]] - ov_qe
      offset_end <- if (q_dir == 1) ov_qe - h$qstart[[1]] else h$qstart[[1]] - ov_qs
      s_a <- h$sstart[[1]] + s_dir * offset_start
      s_b <- h$sstart[[1]] + s_dir * offset_end
      data.frame(
        start = min(s_a, s_b),
        end = max(s_a, s_b),
        width = ov_qe - ov_qs + 1,
        pident = suppressWarnings(as.numeric(h$pident[[1]])),
        stringsAsFactors = FALSE
      )
    })
    pieces <- dplyr::bind_rows(Filter(Negate(is.null), pieces))
    if (!nrow(pieces)) {
      return(NULL)
    }
    data.frame(
      segment_id = segment_tbl$segment_id[[i]],
      segment_type = segment_tbl$segment_type[[i]],
      query_start = st,
      query_end = en,
      start = min(pieces$start, na.rm = TRUE),
      end = max(pieces$end, na.rm = TRUE),
      query_span = en - st + 1,
      aligned_frac = 100 * sum(pieces$width, na.rm = TRUE) / max(1, en - st + 1),
      pident = stats::weighted.mean(pieces$pident, w = pieces$width, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(Filter(Negate(is.null), out))
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

.dnmb_prophage_union_width <- function(start, end) {
  start <- suppressWarnings(as.numeric(start))
  end <- suppressWarnings(as.numeric(end))
  keep <- is.finite(start) & is.finite(end)
  if (!any(keep)) {
    return(NA_real_)
  }
  iv <- data.frame(
    start = pmin(start[keep], end[keep]),
    end = pmax(start[keep], end[keep]),
    stringsAsFactors = FALSE
  )
  iv <- iv[order(iv$start, iv$end), , drop = FALSE]
  cur_s <- iv$start[[1]]
  cur_e <- iv$end[[1]]
  total <- 0
  if (nrow(iv) > 1L) {
    for (i in 2:nrow(iv)) {
      s <- iv$start[[i]]
      e <- iv$end[[i]]
      if (s <= cur_e + 1) {
        cur_e <- max(cur_e, e)
      } else {
        total <- total + (cur_e - cur_s + 1)
        cur_s <- s
        cur_e <- e
      }
    }
  }
  total + (cur_e - cur_s + 1)
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

.dnmb_prophage_reference_shared_legend_plot <- function(max_len = 5000L) {
  max_len <- suppressWarnings(as.numeric(max_len)[1])
  if (!is.finite(max_len) || max_len <= 0) {
    max_len <- 5000
  }
  scale_len <- if (max_len >= 5000) 5000 else if (max_len >= 2000) 2000 else 1000
  grad_vals <- seq(40, 95, length.out = 30L)
  grad_fun <- scales::col_numeric(palette = c("#D8E5FF", "#7FB1FF", "#1D4ED8"), domain = c(40, 95))
  scale_tbl <- data.frame(x = 0.08, xend = 0.22, y = 0.52, label = .dnmb_fmt_bp_label(scale_len))
  identity_legend_tbl <- data.frame(
    xmin = seq(0.42, 0.72, length.out = length(grad_vals)),
    xmax = seq(0.42, 0.72, length.out = length(grad_vals)) + 0.30 / length(grad_vals),
    ymin = 0.38,
    ymax = 0.66,
    fill_hex = grad_fun(grad_vals)
  )
  identity_label_tbl <- data.frame(x = 0.42, y = 0.80, label = "Identity (%)")
  identity_tick_tbl <- data.frame(x = c(0.42, 0.72 + 0.30 / length(grad_vals)), y = c(0.26, 0.26), label = c("40", "95"))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = scale_tbl, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y), linewidth = 0.5, color = "grey35", inherit.aes = FALSE) +
    ggplot2::geom_text(data = scale_tbl, ggplot2::aes(x = (.data$x + .data$xend) / 2, y = .data$y - 0.18, label = .data$label), size = 3.0, vjust = 1, inherit.aes = FALSE) +
    ggplot2::geom_rect(data = identity_legend_tbl, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, fill = .data$fill_hex), color = NA, inherit.aes = FALSE) +
    ggplot2::geom_text(data = identity_label_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = 0, size = 2.8, color = "#4B5563", inherit.aes = FALSE) +
    ggplot2::geom_text(data = identity_tick_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = c(0, 1), size = 2.5, color = "#6B7280", inherit.aes = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 40, 4, 12))
}

.dnmb_prophage_reference_style_spec <- function(embedded = FALSE, n_ref_rows = 1L) {
  n_ref_rows <- max(1L, as.integer(n_ref_rows[1]))
  multi_ref_extra <- max(0, n_ref_rows - 1L)
  lane_fill_pad <- 0.0045 + multi_ref_extra * 0.0030
  inter_ref_ribbon_clearance <- 0.018 + multi_ref_extra * 0.012
  inter_ref_ribbon_alpha_scale <- max(0.30, 0.62 - multi_ref_extra * 0.14)
  inter_ref_taper_frac <- min(0.34, 0.20 + multi_ref_extra * 0.05)
  inter_ref_taper_min_bp <- 10
  embedded_label_gutter <- 0.016 + max(0, n_ref_rows - 1L) * 0.004
  embedded_label_margin <- 118 + max(0, n_ref_rows - 1L) * 22
  embedded_label_size <- max(2.0, 2.45 - max(0, n_ref_rows - 1L) * 0.18)
  embedded_projection_gain <- 1.00
  list(
    lane_fill_pad = lane_fill_pad,
    lane_cap_outer_multiplier = 1.9,
    lane_cap_inner_multiplier = 1.7,
    query_outline_width = 0.10,
    ref_outline_width = 0.08,
    lane_outline_width = 0.18,
    ribbon_alpha_min = 0.24,
    ribbon_alpha_max = 0.82,
    embedded_base_taper_frac = 0.02,
    embedded_identity_taper_frac = 0.10,
    embedded_taper_min_bp = 4,
    inter_ref_ribbon_clearance = inter_ref_ribbon_clearance,
    inter_ref_ribbon_alpha_scale = inter_ref_ribbon_alpha_scale,
    inter_ref_taper_frac = inter_ref_taper_frac,
    inter_ref_taper_min_bp = inter_ref_taper_min_bp,
    embedded_projection_gain = embedded_projection_gain,
    label_gutter_frac = if (embedded) embedded_label_gutter else 0.018,
    label_right_margin = if (embedded) embedded_label_margin else 165,
    label_size = if (embedded) embedded_label_size else 2.8,
    panel_bottom_pad = if (embedded) 0.02 else 0.05,
    panel_top_pad = if (embedded) 0 else 0.05
  )
}

.dnmb_prophage_project_interval_to_query_axis <- function(start, end, src_min, src_max, query_len) {
  start <- suppressWarnings(as.numeric(start)[1])
  end <- suppressWarnings(as.numeric(end)[1])
  src_min <- suppressWarnings(as.numeric(src_min)[1])
  src_max <- suppressWarnings(as.numeric(src_max)[1])
  query_len <- suppressWarnings(as.numeric(query_len)[1])
  if (!all(is.finite(c(start, end, src_min, src_max, query_len))) ||
      src_max <= src_min || query_len <= 0) {
    return(c(start, end))
  }
  st <- min(start, end)
  en <- max(start, end)
  if ((src_max - src_min) <= 1 || query_len <= 1) {
    return(c(1, max(1, query_len)))
  }
  p_st <- 1 + (st - src_min) * (query_len - 1) / (src_max - src_min)
  p_en <- 1 + (en - src_min) * (query_len - 1) / (src_max - src_min)
  c(max(1, p_st), min(query_len, p_en))
}

.dnmb_prophage_rescale_interval_to_frame <- function(start, end, src_len, dst_min, dst_max) {
  start <- suppressWarnings(as.numeric(start)[1])
  end <- suppressWarnings(as.numeric(end)[1])
  src_len <- suppressWarnings(as.numeric(src_len)[1])
  dst_min <- suppressWarnings(as.numeric(dst_min)[1])
  dst_max <- suppressWarnings(as.numeric(dst_max)[1])
  if (!all(is.finite(c(start, end, src_len, dst_min, dst_max))) || src_len <= 0 || dst_max <= dst_min) {
    return(c(start, end))
  }
  st <- min(start, end)
  en <- max(start, end)
  if (src_len <= 1) {
    return(c(dst_min, dst_max))
  }
  span <- dst_max - dst_min
  p_st <- dst_min + (st - 1) * span / (src_len - 1)
  p_en <- dst_min + (en - 1) * span / (src_len - 1)
  c(max(dst_min, p_st), min(dst_max, p_en))
}

.dnmb_prophage_rescale_interval_between_domains <- function(start, end, src_min, src_max, dst_min, dst_max) {
  start <- suppressWarnings(as.numeric(start)[1])
  end <- suppressWarnings(as.numeric(end)[1])
  src_min <- suppressWarnings(as.numeric(src_min)[1])
  src_max <- suppressWarnings(as.numeric(src_max)[1])
  dst_min <- suppressWarnings(as.numeric(dst_min)[1])
  dst_max <- suppressWarnings(as.numeric(dst_max)[1])
  if (!all(is.finite(c(start, end, src_min, src_max, dst_min, dst_max))) || src_max <= src_min || dst_max <= dst_min) {
    return(c(start, end))
  }
  st <- min(start, end)
  en <- max(start, end)
  p_st <- dst_min + (st - src_min) * (dst_max - dst_min) / (src_max - src_min)
  p_en <- dst_min + (en - src_min) * (dst_max - dst_min) / (src_max - src_min)
  c(max(dst_min, p_st), min(dst_max, p_en))
}

.dnmb_prophage_comparison_layout <- function(n_ref_rows) {
  n_ref_rows <- max(1L, as.integer(n_ref_rows[1]))
  gene_full_height <- 0.08 * 0.60
  lane_inner_half <- gene_full_height / 2
  lane_border_gap <- 0.004
  lane_outer_half <- lane_inner_half + lane_border_gap
  row_pitch <- 0.20
  ref_base_y <- 0.20
  ref_top_y <- ref_base_y + (n_ref_rows - 1L) * row_pitch
  ref_centers <- ref_top_y - (seq_len(n_ref_rows) - 1L) * row_pitch
  query_y <- ref_top_y + row_pitch
  list(
    gene_full_height = gene_full_height,
    lane_inner_half = lane_inner_half,
    lane_outer_half = lane_outer_half,
    lane_border_gap = lane_border_gap,
    row_pitch = row_pitch,
    ref_centers = ref_centers,
    query_y = query_y,
    panel_top = query_y + lane_outer_half + 0.03,
    panel_bottom = min(ref_centers) - lane_outer_half - 0.08
  )
}

.dnmb_plot_prophage_reference_panel <- function(sub_tbl, region_id, output_dir, top_n = 5L, cache_root = NULL, embedded = FALSE, show_panel_legend = FALSE) {
  text_family <- "sans"
  sub_tbl <- .dnmb_prophage_apply_detected_bounds(sub_tbl, output_dir = output_dir, region_id = region_id)
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
  region_start_col <- .dnmb_pick_column(sub_tbl, c("Prophage_prophage_start", "prophage_start"))
  region_end_col <- .dnmb_pick_column(sub_tbl, c("Prophage_prophage_end", "prophage_end"))
  region_start <- if (!is.null(region_start_col) && any(!is.na(sub_tbl[[region_start_col]]))) {
    suppressWarnings(min(as.numeric(sub_tbl[[region_start_col]]), na.rm = TRUE))
  } else {
    min(pmin(as.numeric(sub_tbl$start), as.numeric(sub_tbl$end)), na.rm = TRUE)
  }
  region_end <- if (!is.null(region_end_col) && any(!is.na(sub_tbl[[region_end_col]]))) {
    suppressWarnings(max(as.numeric(sub_tbl[[region_end_col]]), na.rm = TRUE))
  } else {
    max(pmax(as.numeric(sub_tbl$start), as.numeric(sub_tbl$end)), na.rm = TRUE)
  }
  if (!is.finite(region_start) || !is.finite(region_end) || region_end <= region_start) {
    return(NULL)
  }
  query_len <- region_end - region_start + 1
  query_region_seq <- substr(query_sequence, region_start, region_end)
  cache_dir <- .dnmb_prophage_reference_cache_dir(output_dir, region_id)
  query_fasta <- file.path(cache_dir, "query_region.fna")
  .dnmb_write_single_fasta(paste0("Prophage_", region_id), query_region_seq, query_fasta)
  hit_tbl <- .dnmb_prophage_remote_reference_hits(query_fasta, cache_dir = cache_dir, top_n = max(12L, top_n * 4L), cache_root = cache_root)
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
  row_length_map <- c("Query prophage" = query_len)
  row_effective_len_map <- c("Query prophage" = query_len)

  for (i in seq_len(nrow(hit_tbl))) {
    acc <- hit_tbl$accession[[i]]
    gbk_path <- .dnmb_fetch_genbank_text(acc, file.path(cache_dir, paste0(acc, ".gbk")))
    if (!nzchar(gbk_path) || !file.exists(gbk_path)) next
    ref <- .dnmb_parse_single_genbank_record(gbk_path)
    if (!nzchar(ref$sequence)) next
    ref_total_len <- nchar(ref$sequence)
    ref_fasta <- file.path(cache_dir, paste0(acc, ".fna"))
    .dnmb_write_single_fasta(acc, ref$sequence, ref_fasta)
    aln_tbl <- .dnmb_prophage_pairwise_alignment(query_fasta, ref_fasta, file.path(cache_dir, paste0(acc, "_pairwise.tsv")))
    if (!nrow(aln_tbl)) next
    ref_interval <- .dnmb_prophage_reference_interval_from_blast(ref$sequence, aln_tbl)
    if (!nzchar(ref_interval$sequence)) next
    mafft_aln <- .dnmb_prophage_mafft_pairwise_alignment(query_region_seq, ref_interval$sequence, file.path(cache_dir, paste0(acc, "_stable")))
    if (is.null(mafft_aln)) {
      seg_profile <- .dnmb_prophage_segment_profile_from_alignment(
        segment_tbl = query_segments,
        aln_tbl = aln_tbl,
        query_len = query_len
      )
      if (!nrow(seg_profile)) next
      region_tbl <- .dnmb_prophage_alignment_regions(aln_tbl)
      profile_identity <- .dnmb_prophage_alignment_identity(aln_tbl)
      profile_coverage <- .dnmb_prophage_alignment_coverage(region_tbl, query_len)
      muts <- data.frame()
    } else {
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
      profile_identity <- profile$identity
      profile_coverage <- profile$coverage
      muts <- profile$mutations
      if (nrow(muts)) {
        if (identical(ref_interval$strand, "-")) {
          muts$position <- ref_interval$ref_end - muts$ref_position + 1
        } else {
          muts$position <- ref_interval$ref_start + muts$ref_position - 1
        }
      }
    }
    # Convert to full reference-genome coordinates
    if (is.null(mafft_aln)) {
      seg_profile$start <- pmax(1, pmin(ref_total_len, seg_profile$start))
      seg_profile$end   <- pmax(1, pmin(ref_total_len, seg_profile$end))
    } else if (identical(ref_interval$strand, "-")) {
      full_start <- ref_interval$ref_end - seg_profile$end + 1
      full_end <- ref_interval$ref_end - seg_profile$start + 1
      seg_profile$start <- pmin(full_start, full_end)
      seg_profile$end <- pmax(full_start, full_end)
    } else {
      seg_profile$start <- ref_interval$ref_start + seg_profile$start - 1
      seg_profile$end <- ref_interval$ref_start + seg_profile$end - 1
    }
    seg_profile$start <- pmax(1, pmin(ref_total_len, seg_profile$start))
    seg_profile$end   <- pmax(1, pmin(ref_total_len, seg_profile$end))
    seg_profile$ref_window_start <- ref_interval$ref_start
    seg_profile$ref_window_end <- ref_interval$ref_end
    seg_profile$ref_total_len <- ref_total_len
    seg_profile$effective_ref_span <- max(
      1,
      .dnmb_prophage_union_width(seg_profile$start, seg_profile$end)
    )
    seg_profile$row <- acc
    seg_profile$fill_hex <- ifelse(
      !is.na(seg_profile$pident) & !is.na(seg_profile$aligned_frac) & seg_profile$pident >= 40 & seg_profile$aligned_frac >= 10,
      seg_fill_fun(pmin(95, pmax(40, seg_profile$pident))),
      "#FFFFFF"
    )
    seg_profile$seg_height <- ifelse(seg_profile$segment_type == "gene", 0.072, 0.026)
    ref_seg_rows[[length(ref_seg_rows) + 1L]] <- seg_profile
    if (nrow(muts)) {
      muts$row <- acc
      ref_mut_rows[[length(ref_mut_rows) + 1L]] <- muts
    }
    ref_short <- .dnmb_prophage_reference_short_name(hit_tbl$title[[i]])
    if (is.na(ref_short) || !nzchar(ref_short)) {
      ref_short <- acc
    }
    ref_label_rows[[length(ref_label_rows) + 1L]] <- data.frame(
      row = acc,
      label = paste0(
        ref_short,
        " (", .dnmb_fmt_bp_label(ref_total_len), ")",
        "\n",
        "Coverage ", sprintf("%.1f", profile_coverage),
        "\n",
        "Identity (%) ", sprintf("%.1f", profile_identity)
      ),
      stringsAsFactors = FALSE
    )
    row_identity_map[[acc]] <- profile_identity
    row_coverage_map[[acc]] <- profile_coverage
    row_length_map[[acc]] <- ref_total_len
    cov_scaled_len <- max(1, query_len * (profile_coverage / 100))
    row_effective_len_map[[acc]] <- min(seg_profile$effective_ref_span[[1]], cov_scaled_len)
  }

  segment_tbl <- dplyr::bind_rows(ref_seg_rows)
  snv_tbl <- dplyr::bind_rows(ref_mut_rows)
  label_tbl <- dplyr::bind_rows(ref_label_rows)
  if (!nrow(segment_tbl)) {
    return(NULL)
  }
  candidate_meta <- data.frame(
    row = names(row_identity_map),
    identity = as.numeric(row_identity_map),
    coverage = as.numeric(row_coverage_map),
    score = as.numeric(row_identity_map) * as.numeric(row_coverage_map) / 100,
    row_len = as.numeric(row_length_map[names(row_identity_map)]),
    cluster_key = vapply(names(row_identity_map), function(r) {
      idx <- which(label_tbl$row == r)
      if (!length(idx)) return(NA_character_)
      .dnmb_prophage_reference_cluster_key(label_tbl$label[[idx[1]]])
    }, character(1)),
    stringsAsFactors = FALSE
  )
  keep_rows <- .dnmb_prophage_prune_reference_candidates(segment_tbl, candidate_meta, max_refs = top_n)
  if (!length(keep_rows)) {
    return(NULL)
  }
  segment_tbl <- segment_tbl[segment_tbl$row %in% keep_rows, , drop = FALSE]
  snv_tbl <- snv_tbl[snv_tbl$row %in% keep_rows, , drop = FALSE]
  label_tbl <- label_tbl[label_tbl$row %in% keep_rows, , drop = FALSE]
  row_identity_map <- row_identity_map[keep_rows]
  row_coverage_map <- row_coverage_map[keep_rows]
  row_length_map <- c("Query prophage" = query_len, row_length_map[keep_rows])
  row_effective_len_map <- c("Query prophage" = query_len, row_effective_len_map[keep_rows])
  # Sort by composite score (coverage * identity) so best overall match is on top
  row_score_map <- row_identity_map * row_coverage_map / 100
  ref_order <- names(sort(row_score_map, decreasing = TRUE))
  row_levels <- c("Query prophage", ref_order)
  layout_spec <- .dnmb_prophage_comparison_layout(length(ref_order))
  style_spec <- .dnmb_prophage_reference_style_spec(
    embedded = embedded,
    n_ref_rows = length(ref_order)
  )
  # Panel frame = query's native bp frame [1, query_len]. Storyboard above uses
  # the same scale, so query bp-per-pixel matches storyboard exactly.
  # Each ref row is rescaled from its native bp [1, ref_native_len] into the panel
  # frame — ref segments therefore show where in the ref a region sits (ribbons
  # become diagonal because ref bp ≠ query bp even when they align), while the
  # visible row width still equals the query width instead of stretching with ref.
  row_native_len_map <- row_length_map
  row_native_len_map[["Query prophage"]] <- query_len
  panel_len <- query_len
  row_xmin_map <- stats::setNames(rep(1, length(row_levels)), row_levels)
  row_xmax_map <- stats::setNames(rep(panel_len, length(row_levels)), row_levels)
  row_map <- c("Query prophage" = layout_spec$query_y, stats::setNames(layout_spec$ref_centers, ref_order))
  draw_row_levels <- if (embedded) ref_order else row_levels

  rescale_native_to_panel <- function(pos, native_len) {
    native_len <- as.numeric(native_len)
    pos <- as.numeric(pos)
    out <- pos
    ok <- is.finite(native_len) & native_len > 1 & is.finite(panel_len) & panel_len > 1
    if (length(native_len) == 1L) {
      if (isTRUE(ok)) out <- 1 + (pos - 1) * (panel_len - 1) / (native_len - 1)
    } else {
      idx <- which(ok & is.finite(pos))
      if (length(idx)) out[idx] <- 1 + (pos[idx] - 1) * (panel_len - 1) / (native_len[idx] - 1)
    }
    out
  }

  segment_tbl$y <- row_map[segment_tbl$row]
  snv_tbl$y <- row_map[snv_tbl$row]

  # Ref segment bottom x: native ref bp rescaled so the ref's full length spans the
  # same width as the query. Diagonal ribbons emerge because query_start/query_end
  # (top) and the ref's own ref_start/ref_end position (bottom) are on different
  # percentages of each row's native length.
  seg_native_len <- as.numeric(row_native_len_map[segment_tbl$row])
  segment_tbl$plot_start <- pmax(1, pmin(panel_len, rescale_native_to_panel(segment_tbl$start, seg_native_len)))
  segment_tbl$plot_end   <- pmax(1, pmin(panel_len, rescale_native_to_panel(segment_tbl$end,   seg_native_len)))
  if (is.data.frame(label_tbl) && nrow(label_tbl) && "row" %in% names(label_tbl)) {
    label_tbl$y <- row_map[label_tbl$row]
    label_tbl$x_label <- as.numeric(row_xmax_map[label_tbl$row]) + panel_len * style_spec$label_gutter_frac
    label_tbl$label_hjust <- 0
  } else {
    label_tbl <- data.frame(row = character(), label = character(), y = numeric(), stringsAsFactors = FALSE)
  }

  lane_outer_half <- layout_spec$lane_outer_half
  lane_border_gap <- layout_spec$lane_border_gap
  lane_inner_half <- layout_spec$lane_inner_half
  lane_fill_pad <- style_spec$lane_fill_pad
  lane_cap_x_radius_outer <- function(row_len) lane_outer_half * style_spec$lane_cap_outer_multiplier
  lane_cap_x_radius_inner <- function(row_len) lane_inner_half * style_spec$lane_cap_inner_multiplier
  gene_bar_half <- lane_inner_half - lane_fill_pad
  aux_bar_half <- gene_bar_half
  snv_tick_half <- lane_inner_half - 0.004
  query_anchor_y <- unname(row_map[["Query prophage"]]) - gene_bar_half

  query_bar_tbl <- query_genes
  # Query genes are already in query-local bp [1, query_len] — that is also the
  # panel frame, so no rescale is needed.
  query_bar_tbl$xmin <- pmax(1, pmin(panel_len, as.numeric(query_bar_tbl$plot_start)))
  query_bar_tbl$xmax <- pmax(1, pmin(panel_len, as.numeric(query_bar_tbl$plot_end)))
  query_bar_tbl$y <- row_map[["Query prophage"]]
  query_bar_tbl$bar_hex <- unname(palette[query_bar_tbl$category])
  query_bar_tbl$bar_hex[is.na(query_bar_tbl$bar_hex)] <- "#9270CA"
  query_bar_tbl$bar_half <- gene_bar_half
  query_bar_tbl$ymin <- query_bar_tbl$y - query_bar_tbl$bar_half
  query_bar_tbl$ymax <- query_bar_tbl$y + query_bar_tbl$bar_half
  segment_tbl$render_half <- ifelse(segment_tbl$segment_type == "gene", gene_bar_half, aux_bar_half)
  segment_tbl$ymin <- segment_tbl$y - segment_tbl$render_half
  segment_tbl$ymax <- segment_tbl$y + segment_tbl$render_half

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
      # Get upper row segments for ref-to-ref ribbons
      upper_segs <- if (!identical(upper_row, "Query prophage")) {
        segment_tbl[segment_tbl$row == upper_row & segment_tbl$segment_type == "gene", , drop = FALSE]
      } else {
        NULL
      }
      for (j in seq_len(nrow(lower_segs))) {
        # All ribbon endpoints live on the query-anchored panel frame [1, query_len],
        # so ribbons are simple rectangles spanning the query bp range that aligns
        # between the two rows.
        qs <- pmax(1, pmin(panel_len, as.numeric(lower_segs$query_start[[j]])))
        qe <- pmax(1, pmin(panel_len, as.numeric(lower_segs$query_end[[j]])))
        if (identical(upper_row, "Query prophage")) {
          top_left <- qs
          top_right <- qe
          top_half <- gene_bar_half
        } else if (!is.null(upper_segs) && nrow(upper_segs)) {
          match_idx <- which(upper_segs$segment_id == lower_segs$segment_id[[j]])
          if (length(match_idx)) {
            us <- upper_segs[match_idx[1], ]
            top_left <- if (is.finite(us$plot_start)) us$plot_start else qs
            top_right <- if (is.finite(us$plot_end))   us$plot_end   else qe
            top_half <- if (is.finite(us$render_half)) us$render_half else gene_bar_half
          } else {
            top_left <- qs
            top_right <- qe
            top_half <- gene_bar_half
          }
        } else {
          top_left <- qs
          top_right <- qe
          top_half <- gene_bar_half
        }
        bottom_left  <- if (is.finite(lower_segs$plot_start[[j]])) lower_segs$plot_start[[j]] else qs
        bottom_right <- if (is.finite(lower_segs$plot_end[[j]]))   lower_segs$plot_end[[j]]   else qe
        bottom_half <- if (is.finite(lower_segs$render_half[[j]])) lower_segs$render_half[[j]] else aux_bar_half
        top_y <- if (embedded && identical(upper_row, "Query prophage")) query_anchor_y else upper_y - top_half
        bottom_y <- lower_y + bottom_half
        pident_val <- lower_segs$pident[[j]]
        if (!identical(upper_row, "Query prophage")) {
          top_y <- top_y - style_spec$inter_ref_ribbon_clearance
          bottom_y <- bottom_y + style_spec$inter_ref_ribbon_clearance
        }
        top_w <- max(1, top_right - top_left + 1)
        bottom_w <- max(1, bottom_right - bottom_left + 1)
        if (embedded) {
          # Embedded panels keep projected alignment, but use a small identity-aware
          # taper so ribbons remain readable instead of collapsing into vertical boxes.
          taper_frac <- style_spec$embedded_base_taper_frac +
            (1 - pmin(100, pmax(40, pident_val)) / 100) * style_spec$embedded_identity_taper_frac
          top_taper <- min(top_w * (taper_frac * 0.6), max(style_spec$embedded_taper_min_bp, top_w * 0.01))
          bottom_taper <- min(bottom_w * taper_frac, max(style_spec$embedded_taper_min_bp, bottom_w * 0.015))
        } else {
          taper_frac <- if (embedded) style_spec$inter_ref_taper_frac else 0.10
          taper_min_bp <- if (embedded) style_spec$inter_ref_taper_min_bp else 6
          top_taper <- min(top_w * taper_frac, max(taper_min_bp, top_w * 0.16))
          bottom_taper <- min(bottom_w * taper_frac, max(taper_min_bp, bottom_w * 0.16))
        }
        alpha_val <- pmin(style_spec$ribbon_alpha_max, pmax(style_spec$ribbon_alpha_min, scales::rescale(
          pident_val,
          to = c(style_spec$ribbon_alpha_min, style_spec$ribbon_alpha_max),
          from = c(40, 100)
        )))
        if (!identical(upper_row, "Query prophage")) {
          alpha_val <- alpha_val * style_spec$inter_ref_ribbon_alpha_scale
        }
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
  if (show_snv) {
    snv_palette <- c("Synonymous" = "#A5BFF9", "Nonsynonymous" = "#2563EB", "Noncoding" = "#D9468B")
    snv_tbl$color_hex <- unname(snv_palette[snv_tbl$class])
    snv_tbl$color_hex[is.na(snv_tbl$color_hex)] <- "#94A3B8"
    # SNVs are plotted at the query bp position where the mismatch lives, since the
    # panel is now query-anchored.
    snv_tbl$plot_position <- pmax(1, pmin(panel_len, as.numeric(snv_tbl$query_position)))
    snv_tbl <- dplyr::bind_rows(lapply(split(snv_tbl, snv_tbl$row), function(d) {
      snv_allowed_half <- .dnmb_prophage_capsule_half_height_at_x(
        x = d$plot_position,
        lane_xmin = 1,
        lane_xmax = panel_len,
        half_height = lane_inner_half,
        x_radius = lane_cap_x_radius_inner(panel_len)
      )
      snv_half <- pmin(snv_tick_half, pmax(0.004, snv_allowed_half - 0.003))
      d$ymin <- d$y - snv_half
      d$ymax <- d$y + snv_half
      d
    }))
  }

  lane_tbl <- data.frame(
    row = draw_row_levels,
    y = unname(row_map[draw_row_levels]),
    x_min = as.numeric(row_xmin_map[draw_row_levels]),
    x_max = as.numeric(row_xmax_map[draw_row_levels]),
    stringsAsFactors = FALSE
  )
  lane_outer_poly <- dplyr::bind_rows(lapply(seq_len(nrow(lane_tbl)), function(i) {
    row_len <- lane_tbl$x_max[[i]] - lane_tbl$x_min[[i]] + 1
    poly <- .dnmb_prophage_capsule_polygon(lane_tbl$x_min[[i]], lane_tbl$x_max[[i]], lane_tbl$y[[i]], lane_outer_half, x_radius = lane_cap_x_radius_outer(row_len))
    if (!nrow(poly)) return(NULL)
    poly$group <- paste0("lane_outer_", i)
    poly
  }))
  lane_inner_poly <- dplyr::bind_rows(lapply(seq_len(nrow(lane_tbl)), function(i) {
    row_len <- lane_tbl$x_max[[i]] - lane_tbl$x_min[[i]] + 1
    poly <- .dnmb_prophage_capsule_polygon(lane_tbl$x_min[[i]], lane_tbl$x_max[[i]], lane_tbl$y[[i]], lane_inner_half, x_radius = lane_cap_x_radius_inner(row_len))
    if (!nrow(poly)) return(NULL)
    poly$group <- paste0("lane_inner_", i)
    poly
  }))
  query_bar_poly <- dplyr::bind_rows(lapply(seq_len(nrow(query_bar_tbl)), function(i) {
    poly <- .dnmb_prophage_clipped_bar_polygon(
      xmin = query_bar_tbl$xmin[[i]],
      xmax = query_bar_tbl$xmax[[i]],
      y_center = query_bar_tbl$y[[i]],
      bar_half = query_bar_tbl$bar_half[[i]],
      lane_xmin = 1,
      lane_xmax = panel_len,
      lane_half = lane_inner_half,
      lane_x_radius = lane_cap_x_radius_inner(panel_len)
    )
    if (!nrow(poly)) return(NULL)
    poly$group <- paste0("query_bar_", i)
    poly$fill_hex <- query_bar_tbl$bar_hex[[i]]
    poly
  }))
  ref_bar_poly <- dplyr::bind_rows(lapply(seq_len(nrow(segment_tbl)), function(i) {
    poly <- .dnmb_prophage_clipped_bar_polygon(
      xmin = segment_tbl$plot_start[[i]],
      xmax = segment_tbl$plot_end[[i]],
      y_center = segment_tbl$y[[i]],
      bar_half = segment_tbl$render_half[[i]],
      lane_xmin = as.numeric(row_xmin_map[[segment_tbl$row[[i]]]]),
      lane_xmax = as.numeric(row_xmax_map[[segment_tbl$row[[i]]]]),
      lane_half = lane_inner_half,
      lane_x_radius = lane_cap_x_radius_inner(as.numeric(row_xmax_map[[segment_tbl$row[[i]]]] - row_xmin_map[[segment_tbl$row[[i]]]] + 1))
    )
    if (!nrow(poly)) return(NULL)
    poly$group <- paste0("ref_bar_", i)
    poly$fill_hex <- segment_tbl$fill_hex[[i]]
    poly
  }))

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = lane_outer_poly,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      fill = "#C7CED7",
      color = NA,
      inherit.aes = FALSE
    )
  if (nrow(ribbon_tbl)) {
    p <- p + ggplot2::geom_polygon(
      data = ribbon_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$fill_hex, alpha = .data$alpha_val),
      color = NA, inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::geom_polygon(
      data = lane_inner_poly,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      fill = "#F4F6F8",
      color = NA,
      inherit.aes = FALSE
    )
  if (!embedded) {
    p <- p +
      ggplot2::geom_polygon(
        data = query_bar_poly,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$fill_hex),
        color = NA,
        inherit.aes = FALSE
      )
  }
  p <- p +
    ggplot2::geom_polygon(
      data = ref_bar_poly,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$fill_hex),
      color = NA,
      alpha = 0.94,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  if (show_snv && nrow(snv_tbl)) {
    p <- p + ggplot2::geom_segment(
      data = snv_tbl,
      ggplot2::aes(x = .data$plot_position, xend = .data$plot_position, y = .data$ymin, yend = .data$ymax, color = .data$color_hex),
      linewidth = 0.34
    )
  }
  if (!embedded) {
    p <- p +
      ggplot2::geom_polygon(
        data = query_bar_poly,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
        fill = NA,
        color = "#AEB8C3",
        linewidth = style_spec$query_outline_width,
        inherit.aes = FALSE
      )
  }
  p <- p +
    ggplot2::geom_polygon(
      data = ref_bar_poly,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      fill = NA,
      color = "#B7C0CA",
      linewidth = style_spec$ref_outline_width,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_polygon(
      data = lane_outer_poly,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      fill = NA,
      color = "#B8C1CC",
      linewidth = style_spec$lane_outline_width,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(data = label_tbl, ggplot2::aes(x = .data$x_label, y = .data$y, label = .data$label), hjust = label_tbl$label_hjust[[1]], vjust = 0.5, size = style_spec$label_size, lineheight = 0.92, family = text_family, inherit.aes = FALSE) +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(
      xlim = c(1, panel_len),
      ylim = c(
        min(lane_tbl$y - lane_outer_half) - 0.05,
        if (embedded) query_anchor_y else max(lane_tbl$y + lane_outer_half) + style_spec$panel_top_pad
      ),
      clip = "off"
    ) +
    ggplot2::labs(title = NULL, x = "Phage coordinate (bp)", y = NULL) +
    ggplot2::theme_bw(base_size = 11, base_family = text_family) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = if (embedded) {
        ggplot2::margin(0, style_spec$label_right_margin, 2, 12)
      } else {
        ggplot2::margin(2, style_spec$label_right_margin, 14, 12)
      }
    ) +
    ggplot2::geom_hline(yintercept = 0.6, linewidth = 0.3, color = "#D1D5DB")
  panel_ymin <- min(lane_tbl$y - lane_outer_half) - style_spec$panel_bottom_pad
  panel_ymax <- max(
    if (embedded) query_anchor_y else -Inf,
    max(lane_tbl$y + lane_outer_half),
    na.rm = TRUE
  ) + style_spec$panel_top_pad
  attr(p, "embedded_ok") <- TRUE
  attr(p, "reference_rows") <- setdiff(row_levels, "Query prophage")
  attr(p, "reference_row_count") <- length(draw_row_levels)
  attr(p, "panel_ymin") <- panel_ymin
  attr(p, "panel_ymax") <- panel_ymax
  attr(p, "track_top") <- max(lane_tbl$y + lane_outer_half)
  attr(p, "track_bottom") <- min(lane_tbl$y - lane_outer_half)
  attr(p, "first_row_center") <- max(lane_tbl$y)
  attr(p, "row_pitch") <- layout_spec$row_pitch
  attr(p, "panel_y_span") <- panel_ymax - panel_ymin
  attr(p, "plot_left_margin_pt") <- 12
  attr(p, "plot_right_margin_pt") <- style_spec$label_right_margin
  p
}

# --- PhiSpy x VirSorter2 consensus helper -----------------------------------

# Returns a small data.frame of VS2 calls that overlap the PhiSpy region, or
# NULL when none match. Used to annotate PhiSpy storyboards with VS2 consensus
# info (group, score, partial flag, hallmark count) on a single summary line.
.dnmb_prophage_vs2_overlap <- function(output_dir, contig_id, region_start_bp, region_end_bp) {
  if (base::is.null(output_dir) || !base::nzchar(output_dir)) return(NULL)
  vs2_path <- base::file.path(output_dir, "dnmb_module_virsorter2", "virsorter2_boundary.tsv")
  if (!base::file.exists(vs2_path)) return(NULL)
  bd <- tryCatch(
    utils::read.delim(vs2_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (base::is.null(bd) || !base::nrow(bd)) return(NULL)
  norm_contig <- function(x) base::sub("\\.[0-9]+$", "", base::as.character(x))
  contig <- norm_contig(bd$seqname)
  query_contig <- norm_contig(contig_id)
  rstart <- suppressWarnings(base::as.numeric(bd$trim_bp_start))
  rend   <- suppressWarnings(base::as.numeric(bd$trim_bp_end))
  keep <- !base::is.na(rstart) & !base::is.na(rend) &
          contig == query_contig &
          rend   >= region_start_bp &
          rstart <= region_end_bp
  if (!base::any(keep)) return(NULL)
  bd <- bd[keep, , drop = FALSE]
  group <- bd$final_max_score_group %||% bd$group
  score <- suppressWarnings(base::as.numeric(bd$final_max_score %||% bd$pr_full))
  partial <- suppressWarnings(base::as.integer(bd$partial %||% 0L)) == 1L
  hallmarks <- suppressWarnings(base::as.integer(bd$hallmark_cnt %||% NA))
  base::data.frame(
    group = base::as.character(group),
    score = score,
    partial = partial,
    hallmarks = hallmarks,
    stringsAsFactors = FALSE
  )
}

.dnmb_prophage_vs2_summary_piece <- function(overlap) {
  if (base::is.null(overlap) || !base::nrow(overlap)) return("")
  o <- overlap[1, , drop = FALSE]
  parts <- base::paste0("VS2: ", o$group)
  if (base::is.finite(o$score)) parts <- base::paste0(parts, base::sprintf(" score=%.2f", o$score))
  if (base::is.finite(o$hallmarks) && o$hallmarks > 0) parts <- base::paste0(parts, " hm=", o$hallmarks)
  if (isTRUE(o$partial)) parts <- base::paste0(parts, " (partial)")
  parts
}

