.dnmb_run_comparative_mobileome <- function(
  focal,
  related_genbanks,
  related_metadata = NULL,
  min_related_completeness = 95,
  min_related_ani = 95,
  output_dir,
  detection_mode = "annotation",
  verbose = TRUE
) {
  if (is.null(related_genbanks) || !length(related_genbanks)) {
    return(.dnmb_empty_comparative())
  }

  related_genbanks <- normalizePath(related_genbanks, mustWork = TRUE)
  comp_dir <- file.path(output_dir, "comparative")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  resources <- .dnmb_detect_comparative_resources()
  status <- resources$status

  metadata_filter <- .dnmb_filter_related_genbanks(
    related_genbanks = related_genbanks,
    related_metadata = related_metadata,
    min_related_completeness = min_related_completeness,
    min_related_ani = min_related_ani
  )
  related_genbanks <- metadata_filter$related_genbanks
  status <- dplyr::bind_rows(status, metadata_filter$status)
  if (!length(related_genbanks)) {
    status <- dplyr::bind_rows(
      status,
      tibble::tibble(component = "comparative_prefilter", status = "empty", detail = "no related genomes passed metadata filters")
    )
    return(.dnmb_empty_comparative(status = status))
  }

  .dnmb_mobileome_message(verbose, "Preparing comparative genomes")
  genomes <- c(list(focal), lapply(seq_along(related_genbanks), function(i) {
    .dnmb_prepare_comparative_genome(
      genbank = related_genbanks[[i]],
      genome_id = paste0("related_", sprintf("%02d", i)),
      detection_mode = detection_mode,
      output_dir = comp_dir,
      verbose = verbose
    )
  }))
  names(genomes) <- vapply(genomes, function(x) x$genome_id, character(1))

  if (!nzchar(resources$diamond)) {
    status <- dplyr::bind_rows(
      status,
      tibble::tibble(component = "comparative", status = "failed", detail = "diamond not found")
    )
    return(.dnmb_empty_comparative(status = status))
  }

  .dnmb_mobileome_message(verbose, "Running reciprocal orthology searches")
  orthology <- .dnmb_run_comparative_orthology(
    focal = focal,
    related = genomes[names(genomes) != focal$genome_id],
    output_dir = comp_dir,
    diamond = resources$diamond,
    verbose = verbose
  )

  .dnmb_mobileome_message(verbose, "Building comparative locus anchors")
  loci <- .dnmb_build_comparative_loci(
    focal = focal,
    related = genomes[names(genomes) != focal$genome_id],
    orthology = orthology
  )

  .dnmb_mobileome_message(verbose, "Inferring occupied-empty matrix")
  occupied_empty <- .dnmb_build_occupied_empty_matrix(
    focal = focal,
    related = genomes[names(genomes) != focal$genome_id],
    loci = loci,
    orthology = orthology
  )

  .dnmb_mobileome_message(verbose, "Scoring comparative hotspots")
  hotspots <- .dnmb_score_comparative_hotspots(
    focal = focal,
    loci = loci,
    occupied_empty = occupied_empty,
    orthology = orthology
  )

  .dnmb_mobileome_message(verbose, "Proposing event chronology")
  chronology <- .dnmb_build_event_chronology(
    focal = focal,
    hotspots = hotspots,
    occupied_empty = occupied_empty,
    orthology = orthology
  )

  locus_master <- .dnmb_build_family_locus_master(
    focal = focal,
    hotspots = hotspots,
    occupied_empty = occupied_empty,
    chronology = chronology
  )

  files <- list(
    orthology_tsv = file.path(comp_dir, "orthology_rbh.tsv"),
    loci_tsv = file.path(comp_dir, "comparative_loci.tsv"),
    occupied_empty_tsv = file.path(comp_dir, "occupied_empty_matrix.tsv"),
    hotspots_tsv = file.path(comp_dir, "comparative_hotspots.tsv"),
    chronology_tsv = file.path(comp_dir, "event_chronology.tsv"),
    locus_master_tsv = file.path(comp_dir, "family_locus_master.tsv")
  )
  .dnmb_write_tsv(orthology, files$orthology_tsv)
  .dnmb_write_tsv(loci, files$loci_tsv)
  .dnmb_write_tsv(occupied_empty, files$occupied_empty_tsv)
  .dnmb_write_tsv(hotspots, files$hotspots_tsv)
  .dnmb_write_tsv(chronology, files$chronology_tsv)
  .dnmb_write_tsv(locus_master, files$locus_master_tsv)

  status <- dplyr::bind_rows(
    status,
    tibble::tibble(
      component = c("comparative_genomes", "comparative_loci", "comparative_hotspots", "event_chronology"),
      status = "ok",
      detail = c(
        paste("related genomes:", length(related_genbanks)),
        paste("loci:", nrow(loci)),
        paste("hotspots:", nrow(hotspots)),
        paste("chronology rows:", nrow(chronology))
      )
    )
  )

  list(
    status = status,
    related_metadata = metadata_filter$metadata,
    genomes = genomes,
    orthology = orthology,
    loci = loci,
    occupied_empty = occupied_empty,
    hotspots = hotspots,
    chronology = chronology,
    locus_master = locus_master,
    files = files
  )
}

.dnmb_filter_related_genbanks <- function(
  related_genbanks,
  related_metadata = NULL,
  min_related_completeness = 95,
  min_related_ani = 95
) {
  status <- tibble::tibble(component = character(), status = character(), detail = character())
  if (is.null(related_metadata)) {
    status <- dplyr::bind_rows(
      status,
      tibble::tibble(component = "comparative_prefilter", status = "skipped", detail = "related_metadata not supplied")
    )
    return(list(related_genbanks = related_genbanks, metadata = tibble::tibble(), status = status))
  }

  meta <- .dnmb_normalize_related_metadata(related_metadata)
  if (!nrow(meta)) {
    status <- dplyr::bind_rows(
      status,
      tibble::tibble(component = "comparative_prefilter", status = "empty", detail = "related_metadata had no usable rows")
    )
    return(list(related_genbanks = related_genbanks, metadata = meta, status = status))
  }

  path_tbl <- tibble::tibble(genbank = related_genbanks) %>%
    dplyr::mutate(
      accession_key = .dnmb_accession_key_from_path(.data$genbank),
      filename_key = basename(.data$genbank)
    )

  meta_joined <- path_tbl %>%
    dplyr::left_join(meta, by = c("accession_key")) %>%
    dplyr::mutate(
      completeness_ok = dplyr::case_when(
        !is.na(.data$completeness) ~ .data$completeness >= min_related_completeness,
        TRUE ~ TRUE
      ),
      ani_ok = dplyr::case_when(
        !is.na(.data$ani_to_focal) ~ .data$ani_to_focal >= min_related_ani,
        TRUE ~ TRUE
      ),
      assembly_ok = dplyr::case_when(
        !is.na(.data$assembly_level) & nzchar(.data$assembly_level) ~ .data$assembly_level %in% c("Complete Genome", "Chromosome", "Scaffold", "Contig"),
        TRUE ~ TRUE
      ),
      keep = .data$completeness_ok & .data$ani_ok & .data$assembly_ok
    )

  filtered <- meta_joined %>% dplyr::filter(.data$keep) %>% dplyr::pull(.data$genbank)
  status <- dplyr::bind_rows(
    status,
    tibble::tibble(
      component = "comparative_prefilter",
      status = "ok",
      detail = paste("kept", length(filtered), "of", length(related_genbanks), "related genomes")
    )
  )

  list(
    related_genbanks = filtered,
    metadata = meta_joined,
    status = status
  )
}

.dnmb_normalize_related_metadata <- function(meta) {
  meta <- tibble::as_tibble(meta)
  if (!nrow(meta)) {
    return(meta)
  }

  choose_col <- function(candidates) {
    hit <- candidates[candidates %in% colnames(meta)]
    if (length(hit)) hit[[1]] else NA_character_
  }

  accession_col <- choose_col(c(
    "RefSeq assembly accession",
    "Assembly Accession",
    "assembly_accession",
    "accession",
    "Filename"
  ))
  completeness_col <- choose_col(c(
    "completeness",
    "Completeness",
    "CheckM completeness",
    "checkm_completeness",
    "completeness_score"
  ))
  ani_col <- choose_col(c(
    "ani_to_focal",
    "ANI",
    "ANI to focal",
    "fastani",
    "skani_ani"
  ))
  assembly_level_col <- choose_col(c(
    "assembly_level",
    "Assembly level",
    "Assembly_Level"
  ))
  source_col <- choose_col(c("SOURCE", "source", "Organism Name", "organism_name"))
  definition_col <- choose_col(c("DEFINITION", "definition"))

  if (is.na(accession_col)) {
    return(tibble::tibble())
  }

  meta %>%
    dplyr::transmute(
      accession_key = .dnmb_accession_key_from_path(as.character(.data[[accession_col]])),
      completeness = if (!is.na(completeness_col)) suppressWarnings(as.numeric(.data[[completeness_col]])) else NA_real_,
      ani_to_focal = if (!is.na(ani_col)) suppressWarnings(as.numeric(.data[[ani_col]])) else NA_real_,
      assembly_level = if (!is.na(assembly_level_col)) as.character(.data[[assembly_level_col]]) else NA_character_,
      source = if (!is.na(source_col)) as.character(.data[[source_col]]) else NA_character_,
      definition = if (!is.na(definition_col)) as.character(.data[[definition_col]]) else NA_character_
    ) %>%
    dplyr::distinct(.data$accession_key, .keep_all = TRUE)
}

.dnmb_accession_key_from_path <- function(x) {
  x <- basename(x)
  x <- sub("\\.(gbk|gbff|gb|fna|fa|fasta)$", "", x, ignore.case = TRUE)
  x
}

.dnmb_empty_comparative <- function(status = NULL) {
  if (is.null(status)) {
    status <- tibble::tibble(component = character(), status = character(), detail = character())
  }
  list(
    status = status,
    genomes = list(),
    orthology = tibble::tibble(),
    loci = tibble::tibble(),
    occupied_empty = tibble::tibble(),
    hotspots = tibble::tibble(),
    chronology = tibble::tibble(),
    locus_master = tibble::tibble(),
    files = list()
  )
}

.dnmb_detect_comparative_resources <- function() {
  diamond <- Sys.which("diamond")
  status <- tibble::tibble(
    component = "diamond",
    status = if (nzchar(diamond)) "found" else "missing",
    detail = diamond
  )
  list(diamond = diamond, status = status)
}

.dnmb_prepare_comparative_genome <- function(genbank, genome_id, detection_mode = "annotation", output_dir, verbose = TRUE) {
  parsed <- .dnmb_parse_genbank_features(genbank)
  genes <- .dnmb_predict_gene_essentiality(.dnmb_build_gene_table(parsed$features))
  annotation_elements <- .dnmb_detect_is_elements(parsed$features, genes)
  sequence_elements <- if (detection_mode %in% c("sequence", "hybrid")) {
    seq_dir <- file.path(output_dir, genome_id, "sequence_engine")
    dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)
    .dnmb_run_sequence_engine(
      parsed = parsed,
      output_dir = file.path(output_dir, genome_id),
      verbose = verbose
    )$sequence_elements
  } else {
    tibble::tibble()
  }
  elements <- .dnmb_merge_element_calls(
    annotation_elements = annotation_elements,
    sequence_elements = sequence_elements,
    genes = genes
  )
  proteins <- genes %>%
    dplyr::filter(!is.na(.data$translation), nzchar(.data$translation)) %>%
    dplyr::transmute(
      genome_id = genome_id,
      gene_id = .data$gene_id,
      locus_tag = .data$locus_tag,
      protein_seq = .data$translation,
      protein_label = paste0(genome_id, "|", .data$gene_id)
    )
  faa_path <- file.path(output_dir, genome_id, paste0(genome_id, ".faa"))
  dir.create(dirname(faa_path), recursive = TRUE, showWarnings = FALSE)
  .dnmb_write_protein_fasta(proteins, faa_path)

  list(
    genome_id = genome_id,
    genbank = genbank,
    metadata = parsed$metadata,
    features = parsed$features,
    genes = genes,
    elements = elements,
    proteins = proteins,
    faa = faa_path
  )
}

.dnmb_write_protein_fasta <- function(proteins, path) {
  if (!nrow(proteins)) {
    writeLines(character(), con = path)
    return(invisible(path))
  }
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(proteins))) {
    seq_clean <- gsub("[^A-Za-z*]", "", proteins$protein_seq[[i]])
    seq_clean <- toupper(seq_clean)
    seq_clean <- gsub("\\*", "", seq_clean)
    if (!nzchar(seq_clean)) next
    writeLines(paste0(">", proteins$protein_label[[i]]), con = con)
    writeLines(seq_clean, con = con)
  }
  invisible(path)
}

.dnmb_run_comparative_orthology <- function(focal, related, output_dir, diamond, verbose = TRUE) {
  if (!length(related)) {
    return(tibble::tibble())
  }

  rows <- list()
  out_dir <- file.path(output_dir, "diamond")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (genome in related) {
    if (!nrow(focal$proteins) || !nrow(genome$proteins)) {
      next
    }
    db_prefix <- file.path(out_dir, genome$genome_id)
    .dnmb_run_system_command(diamond, c("makedb", "--in", shQuote(genome$faa), "--db", shQuote(db_prefix)))
    out1 <- file.path(out_dir, paste0(focal$genome_id, "_vs_", genome$genome_id, ".tsv"))
    .dnmb_run_system_command(diamond, c(
      "blastp",
      "--query", shQuote(focal$faa),
      "--db", shQuote(db_prefix),
      "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
      "--max-target-seqs", "5",
      "--out", shQuote(out1)
    ))

    db_prefix2 <- file.path(out_dir, focal$genome_id)
    .dnmb_run_system_command(diamond, c("makedb", "--in", shQuote(focal$faa), "--db", shQuote(db_prefix2)))
    out2 <- file.path(out_dir, paste0(genome$genome_id, "_vs_", focal$genome_id, ".tsv"))
    .dnmb_run_system_command(diamond, c(
      "blastp",
      "--query", shQuote(genome$faa),
      "--db", shQuote(db_prefix2),
      "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
      "--max-target-seqs", "5",
      "--out", shQuote(out2)
    ))

    fwd <- .dnmb_parse_diamond_tsv(out1, focal$genome_id, genome$genome_id)
    rev <- .dnmb_parse_diamond_tsv(out2, genome$genome_id, focal$genome_id)
    if (!nrow(fwd) || !nrow(rev)) next

    fwd_best <- fwd %>% dplyr::arrange(.data$query_gene, .data$evalue, dplyr::desc(.data$bitscore)) %>% dplyr::group_by(.data$query_gene) %>% dplyr::slice(1) %>% dplyr::ungroup()
    rev_best <- rev %>% dplyr::arrange(.data$query_gene, .data$evalue, dplyr::desc(.data$bitscore)) %>% dplyr::group_by(.data$query_gene) %>% dplyr::slice(1) %>% dplyr::ungroup()

    rbh <- fwd_best %>%
      dplyr::inner_join(
        rev_best,
        by = c("query_gene" = "subject_gene", "subject_gene" = "query_gene"),
        suffix = c("_fwd", "_rev")
      ) %>%
      dplyr::transmute(
        focal_genome = .data$query_genome_fwd,
        related_genome = .data$subject_genome_fwd,
        focal_gene = .data$query_gene,
        related_gene = .data$subject_gene,
        pident = pmin(.data$pident_fwd, .data$pident_rev),
        bitscore = pmin(.data$bitscore_fwd, .data$bitscore_rev),
        evalue = pmax(.data$evalue_fwd, .data$evalue_rev)
      ) %>%
      dplyr::filter(.data$pident >= 30, .data$bitscore >= 50)

    rows[[length(rows) + 1L]] <- rbh
  }

  dplyr::bind_rows(rows)
}

.dnmb_parse_diamond_tsv <- function(path, query_genome, subject_genome) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(tibble::tibble())
  }
  tbl <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!nrow(tbl)) return(tibble::tibble())
  colnames(tbl) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore")
  tibble::as_tibble(tbl) %>%
    dplyr::transmute(
      query_genome = query_genome,
      subject_genome = subject_genome,
      query_gene = sub("^[^|]+\\|", "", .data$qseqid),
      subject_gene = sub("^[^|]+\\|", "", .data$sseqid),
      pident = as.numeric(.data$pident),
      bitscore = as.numeric(.data$bitscore),
      evalue = as.numeric(.data$evalue)
    )
}

.dnmb_build_comparative_loci <- function(focal, related, orthology) {
  if (!nrow(focal$elements) || !nrow(orthology)) {
    return(tibble::tibble())
  }
  focal_genes <- focal$genes %>% dplyr::arrange(.data$contig, .data$start, .data$end)
  rows <- lapply(seq_len(nrow(focal$elements)), function(i) {
    elem <- focal$elements[i, , drop = FALSE]
    left_gene <- focal_genes %>%
      dplyr::filter(.data$contig == elem$contig[[1]], .data$end < elem$start[[1]]) %>%
      dplyr::arrange(dplyr::desc(.data$end)) %>%
      dplyr::slice_head(n = 1)
    right_gene <- focal_genes %>%
      dplyr::filter(.data$contig == elem$contig[[1]], .data$start > elem$end[[1]]) %>%
      dplyr::arrange(.data$start) %>%
      dplyr::slice_head(n = 1)
    tibble::tibble(
      locus_id = paste0("LOC_", sprintf("%04d", i)),
      element_id = elem$element_id[[1]],
      family = elem$element_family[[1]],
      focal_contig = elem$contig[[1]],
      focal_start = elem$start[[1]],
      focal_end = elem$end[[1]],
      left_anchor_gene = if (nrow(left_gene)) left_gene$gene_id[[1]] else NA_character_,
      right_anchor_gene = if (nrow(right_gene)) right_gene$gene_id[[1]] else NA_character_,
      anchor_type = dplyr::case_when(
        nrow(left_gene) & nrow(right_gene) ~ "flanking_pair",
        nrow(left_gene) | nrow(right_gene) ~ "single_anchor",
        TRUE ~ "none"
      )
    )
  })
  dplyr::bind_rows(rows)
}

.dnmb_build_occupied_empty_matrix <- function(focal, related, loci, orthology) {
  if (!nrow(loci) || !length(related)) {
    return(tibble::tibble())
  }

  rows <- list()
  for (genome in related) {
    gene_tbl <- genome$genes %>% dplyr::arrange(.data$contig, .data$start, .data$end)
    elem_tbl <- genome$elements
    orth <- orthology %>% dplyr::filter(.data$related_genome == genome$genome_id)
    for (i in seq_len(nrow(loci))) {
      loc <- loci[i, , drop = FALSE]
      left_related <- orth$related_gene[match(loc$left_anchor_gene[[1]], orth$focal_gene)]
      right_related <- orth$related_gene[match(loc$right_anchor_gene[[1]], orth$focal_gene)]
      left_row <- if (!is.na(left_related)) gene_tbl %>% dplyr::filter(.data$gene_id == left_related) %>% dplyr::slice_head(n = 1) else gene_tbl[0, , drop = FALSE]
      right_row <- if (!is.na(right_related)) gene_tbl %>% dplyr::filter(.data$gene_id == right_related) %>% dplyr::slice_head(n = 1) else gene_tbl[0, , drop = FALSE]

      state <- "unknown"
      region_start <- NA_integer_
      region_end <- NA_integer_
      inserted_family <- NA_character_
      if (nrow(left_row) && nrow(right_row) && identical(left_row$contig[[1]], right_row$contig[[1]])) {
        region_start <- min(left_row$end[[1]], right_row$end[[1]]) + 1L
        region_end <- max(left_row$start[[1]], right_row$start[[1]]) - 1L
        overlapping <- elem_tbl %>%
          dplyr::filter(
            .data$contig == left_row$contig[[1]],
            .data$start <= region_end,
            .data$end >= region_start
          )
        if (!nrow(overlapping)) {
          state <- "empty"
        } else if (any(overlapping$element_family == loc$family[[1]], na.rm = TRUE)) {
          state <- "occupied_same_family"
          inserted_family <- paste(unique(stats::na.omit(overlapping$element_family)), collapse = "; ")
        } else {
          state <- "occupied_other_family"
          inserted_family <- paste(unique(stats::na.omit(overlapping$element_family)), collapse = "; ")
        }
      } else if (nrow(left_row) || nrow(right_row)) {
        state <- "partial_anchor"
      }

      rows[[length(rows) + 1L]] <- tibble::tibble(
        locus_id = loc$locus_id[[1]],
        focal_element_id = loc$element_id[[1]],
        family = loc$family[[1]],
        genome_id = genome$genome_id,
        state = state,
        left_anchor_related = left_related,
        right_anchor_related = right_related,
        region_contig = if (nrow(left_row)) left_row$contig[[1]] else if (nrow(right_row)) right_row$contig[[1]] else NA_character_,
        region_start = region_start,
        region_end = region_end,
        inserted_family = inserted_family
      )
    }
  }
  dplyr::bind_rows(rows)
}

.dnmb_score_comparative_hotspots <- function(focal, loci, occupied_empty, orthology) {
  if (!nrow(loci)) return(tibble::tibble())
  relatedness <- orthology %>%
    dplyr::group_by(.data$related_genome) %>%
    dplyr::summarise(mean_pident = mean(.data$pident, na.rm = TRUE), ortholog_count = dplyr::n(), .groups = "drop")

  loci %>%
    dplyr::left_join(
      occupied_empty %>%
        dplyr::group_by(.data$locus_id) %>%
        dplyr::summarise(
          empty_count = sum(.data$state == "empty", na.rm = TRUE),
          occupied_same_family_count = sum(.data$state == "occupied_same_family", na.rm = TRUE),
          occupied_other_family_count = sum(.data$state == "occupied_other_family", na.rm = TRUE),
          partial_anchor_count = sum(.data$state == "partial_anchor", na.rm = TRUE),
          .groups = "drop"
        ),
      by = "locus_id"
    ) %>%
    dplyr::mutate(
      empty_count = dplyr::coalesce(.data$empty_count, 0L),
      occupied_same_family_count = dplyr::coalesce(.data$occupied_same_family_count, 0L),
      occupied_other_family_count = dplyr::coalesce(.data$occupied_other_family_count, 0L),
      comparative_support_score = pmax(0, .data$empty_count * 1.0 + .data$occupied_same_family_count * 0.5 - .data$occupied_other_family_count * 0.3),
      hotspot_class = dplyr::case_when(
        .data$empty_count >= 2 ~ "recurrently_empty_hotspot",
        .data$occupied_same_family_count >= 2 ~ "shared_occupied_locus",
        TRUE ~ "weak_comparative_signal"
      )
    ) %>%
    dplyr::left_join(
      relatedness %>% dplyr::summarise(
        related_genome_count = dplyr::n(),
        mean_related_pident = mean(.data$mean_pident, na.rm = TRUE)
      ),
      by = character()
    )
}

.dnmb_build_event_chronology <- function(focal, hotspots, occupied_empty, orthology) {
  if (!nrow(hotspots)) return(tibble::tibble())
  relatedness <- orthology %>%
    dplyr::group_by(.data$related_genome) %>%
    dplyr::summarise(mean_pident = mean(.data$pident, na.rm = TRUE), .groups = "drop")

  rows <- lapply(seq_len(nrow(hotspots)), function(i) {
    hot <- hotspots[i, , drop = FALSE]
    states <- occupied_empty %>%
      dplyr::filter(.data$locus_id == hot$locus_id[[1]]) %>%
      dplyr::left_join(relatedness, by = c("genome_id" = "related_genome")) %>%
      dplyr::arrange(dplyr::desc(.data$mean_pident))

    closest_state <- if (nrow(states)) states$state[[1]] else NA_character_
    chronology_class <- dplyr::case_when(
      hot$empty_count[[1]] >= 2 & hot$occupied_same_family_count[[1]] == 0 ~ "single_step_recent_insertion_likely",
      hot$empty_count[[1]] >= 1 & hot$occupied_same_family_count[[1]] >= 1 ~ "multi_step_or_recurrent_usage_likely",
      hot$occupied_same_family_count[[1]] >= 2 ~ "ancestral_or_shared_insertion_likely",
      TRUE ~ "chronology_uncertain"
    )
    proposed_order <- dplyr::case_when(
      chronology_class == "single_step_recent_insertion_likely" ~ "empty locus in related genomes -> insertion in focal lineage",
      chronology_class == "multi_step_or_recurrent_usage_likely" ~ "pre-existing permissive hotspot -> repeated family movement across lineages -> focal occupancy",
      chronology_class == "ancestral_or_shared_insertion_likely" ~ "older insertion event before recent strain divergence",
      TRUE ~ "insufficient comparative evidence"
    )
    tibble::tibble(
      locus_id = hot$locus_id[[1]],
      focal_element_id = hot$element_id[[1]],
      family = hot$family[[1]],
      chronology_class = chronology_class,
      closest_related_state = closest_state,
      proposed_order = proposed_order,
      related_state_summary = if (nrow(states)) paste(paste(states$genome_id, states$state, sep = ":"), collapse = "; ") else NA_character_
    )
  })
  dplyr::bind_rows(rows)
}

.dnmb_build_family_locus_master <- function(focal, hotspots, occupied_empty, chronology) {
  if (!nrow(hotspots)) return(tibble::tibble())
  hotspots %>%
    dplyr::left_join(chronology, by = c("locus_id", "element_id" = "focal_element_id", "family")) %>%
    dplyr::rename(focal_element_id = .data$element_id) %>%
    dplyr::select(
      .data$locus_id,
      .data$focal_element_id,
      .data$family,
      .data$focal_contig,
      .data$focal_start,
      .data$focal_end,
      .data$left_anchor_gene,
      .data$right_anchor_gene,
      .data$anchor_type,
      .data$empty_count,
      .data$occupied_same_family_count,
      .data$occupied_other_family_count,
      .data$comparative_support_score,
      .data$hotspot_class,
      .data$chronology_class,
      .data$proposed_order,
      .data$related_state_summary
    )
}

#' Plot comparative mobileome loci
#'
#' Visualize focal and related genomes around comparative mobileome loci using
#' gggenes. The plot highlights predicted insertions/deletions/transpositions by
#' showing aligned flanking genes and IS families across genomes.
#'
#' @param comparative A comparative result object returned by the DNMB
#'   comparative layer.
#' @param output_dir Directory where plots should be written.
#' @param top_n Number of loci to plot.
#'
#' @return Invisible list with plot file paths.
#' @export
plot_DNMB_mobileome_comparative <- function(comparative, output_dir, top_n = 6) {
  if (is.null(comparative) || !length(comparative$hotspots) || !nrow(comparative$hotspots)) {
    stop("Comparative result does not contain hotspots to plot.")
  }
  plot_dir <- file.path(output_dir, "comparative")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list()
  if (nrow(comparative$occupied_empty)) {
    mat_df <- comparative$occupied_empty %>%
      dplyr::left_join(
        comparative$chronology %>% dplyr::select(.data$locus_id, .data$chronology_class),
        by = "locus_id"
      ) %>%
      dplyr::mutate(
        state = factor(.data$state, levels = c("occupied_same_family", "empty", "occupied_other_family", "partial_anchor", "unknown")),
        genome_id = factor(.data$genome_id, levels = rev(unique(.data$genome_id)))
      )
    p_overview <- ggplot2::ggplot(mat_df, ggplot2::aes(x = .data$locus_id, y = .data$genome_id, fill = .data$state)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::facet_grid(. ~ chronology_class, scales = "free_x", space = "free_x") +
      ggplot2::scale_fill_manual(values = c(
        occupied_same_family = "#2A9D8F",
        empty = "#E9C46A",
        occupied_other_family = "#E76F51",
        partial_anchor = "#A0AEC0",
        unknown = "#E5E7EB"
      )) +
      ggplot2::labs(
        title = "Comparative occupied/empty hotspot matrix",
        subtitle = "Same-family occupancy, empty sites, and other-family occupancy across related genomes",
        x = "Comparative locus",
        y = "Genome",
        fill = "State"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
    files$overview_pdf <- file.path(plot_dir, "comparative_hotspot_matrix.pdf")
    files$overview_png <- file.path(plot_dir, "comparative_hotspot_matrix.png")
    ggplot2::ggsave(files$overview_pdf, p_overview, width = 12, height = 5)
    ggplot2::ggsave(files$overview_png, p_overview, width = 12, height = 5, dpi = 300)
  }

  top_loci <- comparative$hotspots %>%
    dplyr::left_join(
      comparative$chronology %>% dplyr::select(.data$locus_id, .data$proposed_order),
      by = "locus_id"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$comparative_support_score), .data$locus_id) %>%
    dplyr::slice_head(n = as.integer(top_n))

  for (i in seq_len(nrow(top_loci))) {
    locus <- top_loci[i, , drop = FALSE]
    gene_df <- .dnmb_collect_comparative_plot_genes(comparative, locus$locus_id[[1]])
    if (!nrow(gene_df)) next
    p <- ggplot2::ggplot(gene_df, ggplot2::aes(
      xmin = .data$start,
      xmax = .data$end,
      y = .data$genome_id,
      fill = .data$family_fill,
      forward = .data$forward
    )) +
      gggenes::geom_gene_arrow(color = "grey25", alpha = 0.95) +
      ggplot2::facet_wrap(~ locus_id, scales = "free_x", ncol = 1) +
      ggplot2::geom_text(
        data = gene_df %>% dplyr::filter(.data$is_element),
        ggplot2::aes(x = (.data$start + .data$end) / 2, y = .data$genome_id, label = .data$label),
        inherit.aes = FALSE,
        size = 3,
        vjust = -0.6
      ) +
      ggplot2::labs(
        title = paste("Comparative locus", locus$locus_id[[1]], "family", locus$family[[1]]),
        subtitle = dplyr::coalesce(locus$proposed_order[[1]], ""),
        x = "Aligned locus coordinate",
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      )
    files[[paste0("locus_", locus$locus_id[[1]], "_pdf")]] <- file.path(plot_dir, paste0("comparative_", locus$locus_id[[1]], ".pdf"))
    ggplot2::ggsave(files[[paste0("locus_", locus$locus_id[[1]], "_pdf")]], p, width = 12, height = 4)
  }
  invisible(files)
}

.dnmb_collect_comparative_plot_genes <- function(comparative, locus_id) {
  hot <- comparative$hotspots %>% dplyr::filter(.data$locus_id == .env$locus_id)
  if (!nrow(hot)) return(tibble::tibble())
  genomes <- comparative$genomes
  rows <- list()
  for (genome in genomes) {
    gene_tbl <- genome$genes
    elems <- genome$elements
    left_anchor <- if (identical(genome$genome_id, comparative$genomes[[1]]$genome_id)) hot$left_anchor_gene[[1]] else {
      orth <- comparative$orthology %>% dplyr::filter(.data$related_genome == genome$genome_id, .data$focal_gene == hot$left_anchor_gene[[1]])
      if (nrow(orth)) orth$related_gene[[1]] else NA_character_
    }
    right_anchor <- if (identical(genome$genome_id, comparative$genomes[[1]]$genome_id)) hot$right_anchor_gene[[1]] else {
      orth <- comparative$orthology %>% dplyr::filter(.data$related_genome == genome$genome_id, .data$focal_gene == hot$right_anchor_gene[[1]])
      if (nrow(orth)) orth$related_gene[[1]] else NA_character_
    }
    anchor_genes <- gene_tbl %>% dplyr::filter(.data$gene_id %in% c(left_anchor, right_anchor))
    if (!nrow(anchor_genes)) next
    region_start <- max(1L, min(anchor_genes$start) - 5000L)
    region_end <- max(anchor_genes$end) + 5000L
    gene_block <- gene_tbl %>%
      dplyr::filter(.data$contig == anchor_genes$contig[[1]], .data$end >= region_start, .data$start <= region_end) %>%
      dplyr::transmute(
        locus_id = locus_id,
        genome_id = genome$genome_id,
        start = .data$start,
        end = .data$end,
        forward = .data$strand != "-",
        family_fill = "gene",
        label = dplyr::coalesce(.data$gene, .data$locus_tag, .data$gene_id),
        is_element = FALSE
      )
    elem_block <- elems %>%
      dplyr::filter(.data$contig == anchor_genes$contig[[1]], .data$end >= region_start, .data$start <= region_end) %>%
      dplyr::transmute(
        locus_id = locus_id,
        genome_id = genome$genome_id,
        start = .data$start,
        end = .data$end,
        forward = .data$strand != "-",
        family_fill = dplyr::coalesce(.data$element_family, "Unclassified"),
        label = dplyr::coalesce(.data$element_family, .data$element_id),
        is_element = TRUE
      )
    rows[[length(rows) + 1L]] <- dplyr::bind_rows(gene_block, elem_block)
  }
  dplyr::bind_rows(rows)
}
