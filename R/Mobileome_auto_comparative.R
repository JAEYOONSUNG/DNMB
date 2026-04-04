# Mobileome_auto_comparative.R — Automated related genome discovery via NCBI + ANI
#
# Parses organism from GenBank, searches NCBI for related assemblies,
# downloads top N genomes, computes ANI, and selects the best relatives.
# Users just provide a GenBank file — everything else is automatic.

# ---------------------------------------------------------------------------
# Main entry: auto-discover related genomes
# ---------------------------------------------------------------------------
.dnmb_auto_discover_related_genomes <- function(
  genbank,
  output_dir,
  max_related = 5L,
  min_ani = 90,
  assembly_level = c("Complete Genome", "Chromosome"),
  cache_dir = NULL,
  verbose = TRUE
) {
  if (is.null(cache_dir)) {
    cache_dir <- file.path(output_dir, "related_genomes_cache")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: Parse organism from GenBank
  .dnmb_iselement_verbose(verbose, "Parsing organism from GenBank file")
  organism_info <- .dnmb_parse_organism_from_genbank(genbank)
  if (is.null(organism_info$species) || !nzchar(organism_info$species)) {
    .dnmb_iselement_verbose(verbose, "Could not parse organism from GenBank. Skipping auto-comparative.")
    return(.dnmb_empty_auto_comparative())
  }
  .dnmb_iselement_verbose(verbose, paste0("Organism: ", organism_info$species,
                                           if (!is.na(organism_info$strain)) paste0(" (strain: ", organism_info$strain, ")") else ""))

  # Step 2: Search NCBI for related assemblies
  .dnmb_iselement_verbose(verbose, "Searching NCBI for related assemblies")
  assemblies <- .dnmb_search_ncbi_assemblies(
    species = organism_info$species,
    genus = organism_info$genus,
    focal_accession = organism_info$assembly_accession,
    assembly_level = assembly_level,
    max_results = max_related * 3L,
    verbose = verbose
  )
  if (!nrow(assemblies)) {
    .dnmb_iselement_verbose(verbose, "No related assemblies found on NCBI.")
    return(.dnmb_empty_auto_comparative())
  }
  .dnmb_iselement_verbose(verbose, paste0("Found ", nrow(assemblies), " candidate assemblies"))

  # Step 3: Download GenBank files
  .dnmb_iselement_verbose(verbose, "Downloading related GenBank files")
  downloaded <- .dnmb_download_ncbi_genbanks(
    assemblies = assemblies,
    cache_dir = cache_dir,
    verbose = verbose
  )
  if (!length(downloaded$paths)) {
    .dnmb_iselement_verbose(verbose, "No genomes could be downloaded.")
    return(.dnmb_empty_auto_comparative())
  }
  .dnmb_iselement_verbose(verbose, paste0("Downloaded ", length(downloaded$paths), " genomes"))

  # Step 4: Compute ANI
  .dnmb_iselement_verbose(verbose, "Computing ANI against focal genome")
  ani_results <- .dnmb_compute_ani(
    focal_genbank = genbank,
    related_paths = downloaded$paths,
    cache_dir = cache_dir,
    verbose = verbose
  )

  # Step 5: Filter and rank by ANI
  if (!nrow(ani_results)) {
    .dnmb_iselement_verbose(verbose, "ANI computation returned no results.")
    return(.dnmb_empty_auto_comparative())
  }

  selected <- ani_results %>%
    dplyr::filter(.data$ani >= min_ani) %>%
    dplyr::arrange(dplyr::desc(.data$ani)) %>%
    dplyr::slice_head(n = as.integer(max_related))

  if (!nrow(selected)) {
    .dnmb_iselement_verbose(verbose, paste0("No genomes passed ANI threshold (min_ani=", min_ani, ")"))
    return(.dnmb_empty_auto_comparative())
  }

  .dnmb_iselement_verbose(verbose, paste0(
    "Selected ", nrow(selected), " related genomes (ANI range: ",
    sprintf("%.1f", min(selected$ani)), "-", sprintf("%.1f", max(selected$ani)), "%)"
  ))

  # Build metadata table
  metadata <- selected %>%
    dplyr::left_join(downloaded$metadata, by = "genbank_path") %>%
    dplyr::transmute(
      accession_key = .data$accession,
      ani_to_focal = .data$ani,
      completeness = 100,
      assembly_level = dplyr::coalesce(.data$assembly_level, "Complete Genome"),
      organism_name = .data$organism
    )

  list(
    related_genbanks = selected$genbank_path,
    related_metadata = metadata,
    ani_results = ani_results,
    organism_info = organism_info,
    assemblies_searched = assemblies
  )
}

.dnmb_empty_auto_comparative <- function() {
  list(
    related_genbanks = character(),
    related_metadata = tibble::tibble(),
    ani_results = tibble::tibble(),
    organism_info = list(),
    assemblies_searched = tibble::tibble()
  )
}

# ---------------------------------------------------------------------------
# Parse organism name, genus, species, strain from GenBank header
# ---------------------------------------------------------------------------
.dnmb_parse_organism_from_genbank <- function(genbank) {
  lines <- readLines(genbank, n = 100, warn = FALSE)

  # Try SOURCE line first
  source_idx <- grep("^SOURCE\\s+", lines)
  organism_idx <- grep("^\\s+ORGANISM\\s+", lines)
  assembly_idx <- grep("Assembly:\\s+", lines)

  species <- NA_character_
  genus <- NA_character_
  strain <- NA_character_
  assembly_accession <- NA_character_
  taxonomy_lineage <- NA_character_

  if (length(organism_idx)) {
    organism_line <- sub("^\\s+ORGANISM\\s+", "", lines[organism_idx[1]])
    organism_line <- trimws(organism_line)
    # Parse: "Geobacillus stearothermophilus ATCC 12980"
    tokens <- strsplit(organism_line, "\\s+")[[1]]
    if (length(tokens) >= 2) {
      genus <- tokens[1]
      species <- paste(tokens[1:2], collapse = " ")
      if (length(tokens) >= 3) {
        strain <- paste(tokens[3:length(tokens)], collapse = " ")
      }
    } else if (length(tokens) == 1) {
      genus <- tokens[1]
      species <- tokens[1]
    }
    # Get taxonomy lineage from the next line(s)
    if (organism_idx[1] + 1 <= length(lines)) {
      tax_line <- lines[organism_idx[1] + 1]
      if (grepl("^\\s+[A-Z]", tax_line)) {
        taxonomy_lineage <- trimws(sub("\\.$", "", trimws(tax_line)))
      }
    }
  } else if (length(source_idx)) {
    source_line <- sub("^SOURCE\\s+", "", lines[source_idx[1]])
    source_line <- trimws(source_line)
    tokens <- strsplit(source_line, "\\s+")[[1]]
    if (length(tokens) >= 2) {
      genus <- tokens[1]
      species <- paste(tokens[1:2], collapse = " ")
      if (length(tokens) >= 3) {
        strain <- paste(tokens[3:length(tokens)], collapse = " ")
      }
    }
  } else {
    # Fallback: try DEFINITION
    def_idx <- grep("^DEFINITION\\s+", lines)
    if (length(def_idx)) {
      def_line <- sub("^DEFINITION\\s+", "", lines[def_idx[1]])
      def_line <- trimws(def_line)
      tokens <- strsplit(def_line, "\\s+")[[1]]
      if (length(tokens) >= 2) {
        genus <- tokens[1]
        species <- paste(tokens[1:2], collapse = " ")
      }
    }
  }

  if (length(assembly_idx)) {
    assembly_accession <- trimws(sub(".*Assembly:\\s+", "", lines[assembly_idx[1]]))
  }

  list(
    species = species,
    genus = genus,
    strain = strain,
    assembly_accession = assembly_accession,
    taxonomy_lineage = taxonomy_lineage
  )
}

# ---------------------------------------------------------------------------
# Search NCBI Assembly database using Entrez Direct
# ---------------------------------------------------------------------------
.dnmb_search_ncbi_assemblies <- function(
  species,
  genus,
  focal_accession = NA_character_,
  assembly_level = c("Complete Genome", "Chromosome"),
  max_results = 15L,
  verbose = TRUE
) {
  esearch <- Sys.which("esearch")
  efetch <- Sys.which("efetch")

  if (!nzchar(esearch) || !nzchar(efetch)) {
    # Try common conda installation paths
    conda_candidates <- c(
      file.path(Sys.getenv("CONDA_DIR", "/opt/conda"), "bin"),
      path.expand("~/miniforge3/bin"),
      path.expand("~/miniconda3/bin"),
      path.expand("~/anaconda3/bin"),
      "/opt/conda/bin"
    )
    for (cand_dir in conda_candidates) {
      cand_esearch <- file.path(cand_dir, "esearch")
      cand_efetch <- file.path(cand_dir, "efetch")
      if (file.exists(cand_esearch) && file.exists(cand_efetch)) {
        esearch <- cand_esearch
        efetch <- cand_efetch
        break
      }
    }
    if (!nzchar(esearch) || !file.exists(esearch) || !nzchar(efetch) || !file.exists(efetch)) {
      if (isTRUE(verbose)) {
        message("[DNMB iselement] esearch/efetch not found. Install entrez-direct: conda install -c bioconda entrez-direct")
      }
      return(tibble::tibble())
    }
  }

  # Build search query: same species, RefSeq, complete/chromosome level
  level_filter <- paste0("(", paste(paste0('"', assembly_level, '"[Assembly Level]'), collapse = " OR "), ")")
  query <- paste0('"', species, '"[Organism] AND "latest refseq"[filter] AND ', level_filter)

  tryCatch({
    # Run esearch + efetch to get assembly summaries
    cmd <- paste0(
      shQuote(esearch), " -db assembly -query ", shQuote(query),
      " | ", shQuote(efetch), " -format docsum"
    )
    raw <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    if (!length(raw)) {
      # Fallback: search by genus if species search fails
      if (isTRUE(verbose)) {
        message("[DNMB iselement] No results for species '", species, "', trying genus '", genus, "'")
      }
      query <- paste0('"', genus, '"[Organism] AND "latest refseq"[filter] AND ', level_filter)
      cmd <- paste0(
        shQuote(esearch), " -db assembly -query ", shQuote(query),
        " | ", shQuote(efetch), " -format docsum"
      )
      raw <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    }
    if (!length(raw)) return(tibble::tibble())

    xml_text <- paste(raw, collapse = "\n")
    .dnmb_parse_assembly_docsum(xml_text, focal_accession, max_results)
  }, error = function(e) {
    if (isTRUE(verbose)) {
      message("[DNMB iselement] NCBI search failed: ", conditionMessage(e))
    }
    tibble::tibble()
  })
}

.dnmb_parse_assembly_docsum <- function(xml_text, focal_accession = NA, max_results = 15L) {
  # Parse assembly accession, organism, FTP path from XML docsum
  accessions <- regmatches(xml_text, gregexpr('GCF_[0-9]+\\.[0-9]+', xml_text))[[1]]
  organisms <- regmatches(xml_text, gregexpr('<Organism>[^<]+</Organism>', xml_text))[[1]]
  organisms <- gsub('</?Organism>', '', organisms)
  ftp_paths <- regmatches(xml_text, gregexpr('<FtpPath_RefSeq>[^<]+</FtpPath_RefSeq>', xml_text))[[1]]
  ftp_paths <- gsub('</?FtpPath_RefSeq>', '', ftp_paths)
  assembly_levels <- regmatches(xml_text, gregexpr('<AssemblyStatus>[^<]+</AssemblyStatus>', xml_text))[[1]]
  assembly_levels <- gsub('</?AssemblyStatus>', '', assembly_levels)

  n <- min(length(accessions), length(ftp_paths))
  if (!n) return(tibble::tibble())

  tbl <- tibble::tibble(
    accession = accessions[seq_len(n)],
    organism = if (length(organisms) >= n) organisms[seq_len(n)] else rep(NA_character_, n),
    ftp_path = ftp_paths[seq_len(n)],
    assembly_level = if (length(assembly_levels) >= n) assembly_levels[seq_len(n)] else rep(NA_character_, n)
  )

  # Remove focal genome
  if (!is.na(focal_accession) && nzchar(focal_accession)) {
    tbl <- tbl %>% dplyr::filter(.data$accession != focal_accession)
  }

  tbl %>% dplyr::slice_head(n = as.integer(max_results))
}

# ---------------------------------------------------------------------------
# Download GenBank files from NCBI FTP
# ---------------------------------------------------------------------------
.dnmb_download_ncbi_genbanks <- function(assemblies, cache_dir, verbose = TRUE) {
  paths <- character()
  meta_rows <- list()

  for (i in seq_len(nrow(assemblies))) {
    acc <- assemblies$accession[i]
    ftp <- assemblies$ftp_path[i]
    if (is.na(ftp) || !nzchar(ftp)) next

    # Construct download URL for genomic.gbff.gz
    basename_part <- basename(ftp)
    gbff_url <- paste0(ftp, "/", basename_part, "_genomic.gbff.gz")
    # Convert FTP to HTTPS for wget/curl
    gbff_url <- sub("^ftp://", "https://", gbff_url)

    local_gz <- file.path(cache_dir, paste0(acc, ".gbff.gz"))
    local_gbff <- file.path(cache_dir, paste0(acc, ".gbff"))

    if (file.exists(local_gbff) && file.size(local_gbff) > 1000) {
      if (isTRUE(verbose)) message("[DNMB iselement] Using cached: ", acc)
      paths <- c(paths, local_gbff)
      meta_rows[[length(meta_rows) + 1]] <- tibble::tibble(
        genbank_path = local_gbff,
        accession = acc,
        organism = assemblies$organism[i],
        assembly_level = assemblies$assembly_level[i]
      )
      next
    }

    dl_ok <- tryCatch({
      if (isTRUE(verbose)) message("[DNMB iselement] Downloading: ", acc)
      utils::download.file(gbff_url, local_gz, mode = "wb", quiet = !verbose)
      # Decompress
      system2("gunzip", c("-f", local_gz), stdout = FALSE, stderr = FALSE)
      file.exists(local_gbff) && file.size(local_gbff) > 1000
    }, error = function(e) FALSE)

    if (isTRUE(dl_ok)) {
      paths <- c(paths, local_gbff)
      meta_rows[[length(meta_rows) + 1]] <- tibble::tibble(
        genbank_path = local_gbff,
        accession = acc,
        organism = assemblies$organism[i],
        assembly_level = assemblies$assembly_level[i]
      )
    } else {
      if (isTRUE(verbose)) message("[DNMB iselement] Download failed: ", acc)
    }
  }

  list(
    paths = paths,
    metadata = dplyr::bind_rows(meta_rows)
  )
}

# ---------------------------------------------------------------------------
# Compute ANI using available tool (skani > fastANI > BLAST-based fallback)
# ---------------------------------------------------------------------------
.dnmb_compute_ani <- function(focal_genbank, related_paths, cache_dir, verbose = TRUE) {
  if (!length(related_paths)) return(tibble::tibble())

  # Extract FASTA from GenBank for ANI tools
  focal_fna <- file.path(cache_dir, "focal_genome.fna")
  .dnmb_genbank_to_fasta(focal_genbank, focal_fna)

  related_fnas <- vapply(related_paths, function(p) {
    fna <- file.path(cache_dir, paste0(tools::file_path_sans_ext(basename(p)), ".fna"))
    if (!file.exists(fna) || file.size(fna) < 100) {
      .dnmb_genbank_to_fasta(p, fna)
    }
    fna
  }, character(1))

  # Try skani first (fastest), then fastANI, then BLAST fallback
  skani <- Sys.which("skani")
  fastani <- Sys.which("fastANI")

  if (nzchar(skani)) {
    return(.dnmb_run_skani(focal_fna, related_fnas, related_paths, cache_dir, verbose))
  }
  if (nzchar(fastani)) {
    return(.dnmb_run_fastani(focal_fna, related_fnas, related_paths, cache_dir, verbose))
  }

  # BLAST-based ANI fallback
  .dnmb_iselement_verbose(verbose, "No ANI tool found (skani/fastANI). Using BLAST-based ANI estimate.")
  .dnmb_run_blast_ani(focal_fna, related_fnas, related_paths, cache_dir, verbose)
}

.dnmb_genbank_to_fasta <- function(genbank, fasta_out) {
  lines <- readLines(genbank, warn = FALSE)
  locus_idx <- grep("^LOCUS\\s+", lines)
  end_idx <- grep("^//\\s*$", lines)

  con <- file(fasta_out, "w")
  on.exit(close(con))

  for (i in seq_along(locus_idx)) {
    start <- locus_idx[i]
    end_rec <- end_idx[end_idx >= start]
    end_rec <- if (length(end_rec)) end_rec[1] else length(lines)
    block <- lines[start:end_rec]

    contig <- sub("^LOCUS\\s+(\\S+).*", "\\1", block[1])
    origin_start <- grep("^ORIGIN\\s*$", block)
    if (!length(origin_start)) next

    seq_lines <- block[(origin_start[1] + 1):length(block)]
    seq_lines <- seq_lines[!grepl("^//", seq_lines)]
    seq_chars <- gsub("[^a-zA-Z]", "", seq_lines)
    sequence <- paste(seq_chars, collapse = "")
    if (nchar(sequence) < 100) next

    writeLines(paste0(">", contig), con)
    # Write in 70-char lines
    seq_len <- nchar(sequence)
    starts <- seq(1, seq_len, 70)
    for (s in starts) {
      writeLines(substr(sequence, s, min(s + 69, seq_len)), con)
    }
  }
}

.dnmb_run_skani <- function(focal_fna, related_fnas, related_paths, cache_dir, verbose) {
  skani <- Sys.which("skani")
  query_list <- file.path(cache_dir, "related_list.txt")
  writeLines(related_fnas, query_list)
  out_file <- file.path(cache_dir, "skani_results.tsv")

  cmd <- paste(
    shQuote(skani), "dist",
    "-q", shQuote(focal_fna),
    "--rl", shQuote(query_list),
    "-o", shQuote(out_file),
    "-t 2"
  )
  system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)

  if (!file.exists(out_file) || file.size(out_file) < 10) {
    return(tibble::tibble())
  }

  results <- utils::read.delim(out_file, header = TRUE, check.names = FALSE)
  fna_to_gbff <- stats::setNames(related_paths, related_fnas)

  if (nrow(results) && ncol(results) >= 3) {
    tibble::tibble(
      genbank_path = fna_to_gbff[results[[2]]],
      ani = as.numeric(results[[3]])
    ) %>% dplyr::filter(!is.na(.data$ani))
  } else {
    tibble::tibble()
  }
}

.dnmb_run_fastani <- function(focal_fna, related_fnas, related_paths, cache_dir, verbose) {
  fastani <- Sys.which("fastANI")
  query_list <- file.path(cache_dir, "related_list.txt")
  writeLines(related_fnas, query_list)
  out_file <- file.path(cache_dir, "fastani_results.tsv")

  cmd <- paste(
    shQuote(fastani),
    "-q", shQuote(focal_fna),
    "--rl", shQuote(query_list),
    "-o", shQuote(out_file),
    "-t 2"
  )
  system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)

  if (!file.exists(out_file) || file.size(out_file) < 10) {
    return(tibble::tibble())
  }

  results <- utils::read.delim(out_file, header = FALSE, check.names = FALSE)
  fna_to_gbff <- stats::setNames(related_paths, related_fnas)

  if (nrow(results) && ncol(results) >= 3) {
    tibble::tibble(
      genbank_path = fna_to_gbff[results[[2]]],
      ani = as.numeric(results[[3]])
    ) %>% dplyr::filter(!is.na(.data$ani))
  } else {
    tibble::tibble()
  }
}

.dnmb_run_blast_ani <- function(focal_fna, related_fnas, related_paths, cache_dir, verbose) {
  blastn <- Sys.which("blastn")
  makeblastdb <- Sys.which("makeblastdb")

  if (!nzchar(blastn) || !nzchar(makeblastdb)) {
    .dnmb_iselement_verbose(verbose, "blastn/makeblastdb not found. Cannot compute ANI.")
    return(tibble::tibble())
  }

  # Build focal DB
  focal_db <- file.path(cache_dir, "focal_db")
  system2(makeblastdb, c("-in", shQuote(focal_fna), "-dbtype", "nucl", "-out", shQuote(focal_db)),
          stdout = FALSE, stderr = FALSE)

  fna_to_gbff <- stats::setNames(related_paths, related_fnas)
  results <- list()

  for (fna in related_fnas) {
    out <- file.path(cache_dir, paste0(basename(fna), "_blast.tsv"))
    cmd <- paste(
      shQuote(blastn),
      "-query", shQuote(fna),
      "-db", shQuote(focal_db),
      "-outfmt '6 qseqid sseqid pident length mismatch gapopen qlen'",
      "-max_target_seqs 1",
      "-evalue 1e-10",
      "-num_threads 2",
      ">", shQuote(out), "2>/dev/null"
    )
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (file.exists(out) && file.size(out) > 0) {
      hits <- tryCatch(
        utils::read.delim(out, header = FALSE, check.names = FALSE),
        error = function(e) data.frame()
      )
      if (nrow(hits) && ncol(hits) >= 4) {
        # Weighted average identity by alignment length
        total_aligned <- sum(hits[[4]], na.rm = TRUE)
        if (total_aligned > 0) {
          weighted_id <- sum(hits[[3]] * hits[[4]], na.rm = TRUE) / total_aligned
          results[[length(results) + 1]] <- tibble::tibble(
            genbank_path = fna_to_gbff[[fna]],
            ani = weighted_id
          )
        }
      }
    }
  }

  dplyr::bind_rows(results)
}
#' Internal helpers for automatic mobileome comparative genome selection
#'
#' Utility routines used to discover, filter, and prioritize related genomes
#' for comparative mobileome analysis.
#'
#' @name dnmb_internal_mobileome_auto_comparative
#' @keywords internal
#' @noRd
NULL
