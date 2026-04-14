.dnmb_iselement_module_name <- function() {
  "iselement"
}

.dnmb_iselement_default_version <- function() {
  "current"
}

.dnmb_iselement_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_iselement_empty_status <- function() {
  .dnmb_iselement_status_row(character(), character(), character())
}

.dnmb_iselement_normalize_hits <- function(elements) {
  if (base::is.null(elements) || !base::is.data.frame(elements) || !base::nrow(elements)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  tbl <- base::as.data.frame(elements, stringsAsFactors = FALSE)
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$locus_tag)
  tbl <- tbl[!base::is.na(tbl$query) & base::nzchar(tbl$query), , drop = FALSE]
  if (!base::nrow(tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  data.frame(
    query = tbl$query,
    source = "iselement",
    family_system = "ISfinder-like",
    family_id = base::as.character(tbl$element_family),
    hit_label = base::as.character(tbl$mobile_element_type),
    enzyme_role = base::as.character(tbl$feature_type),
    evidence_mode = base::as.character(tbl$confidence),
    substrate_label = NA_character_,
    support = base::as.character(tbl$evidence_summary),
    typing_eligible = base::is.character(tbl$confidence) & tbl$confidence %in% c("high", "medium", "low"),
    element_id = base::as.character(tbl$element_id),
    element_type = base::as.character(tbl$mobile_element_type),
    feature_type = base::as.character(tbl$feature_type),
    confidence = base::as.character(tbl$confidence),
    confidence_score = suppressWarnings(base::as.numeric(tbl$confidence_score)),
    evidence_label = base::as.character(tbl$evidence_label),
    supporting_gene_count = suppressWarnings(base::as.integer(tbl$supporting_gene_count)),
    supporting_genes = base::as.character(tbl$supporting_genes),
    supporting_products = base::as.character(tbl$supporting_products),
    stringsAsFactors = FALSE
  )
}

.dnmb_iselement_output_table <- function(genes, hits) {
  out <- .dnmb_module_output_table(genes = genes, hits = hits)
  drop_cols <- base::intersect(
    c("family_id", "hit_label", "enzyme_role", "evidence_mode", "substrate_label", "typing_eligible"),
    base::names(out)
  )
  if (base::length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(
      c("element_id", "element_type", "feature_type", "confidence", "confidence_score", "family_system", "evidence_label", "supporting_gene_count", "supporting_genes", "support"),
      base::names(out)
    ),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "element_id", "element_type", "feature_type", "confidence", "confidence_score", "family_system", "evidence_label", "supporting_gene_count", "supporting_genes", "support"))
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_iselement_module <- function(genes,
                                      output_dir,
                                      genbank = NULL,
                                      detection_mode = c("hybrid", "annotation"),
                                      prodigal_mode = c("single", "meta"),
                                      isescan_dir = NULL,
                                      isfinder_db = NULL,
                                      isfinder_fasta = NULL,
                                      analysis_depth = c("standard", "full"),
                                      related_genbanks = NULL,
                                      related_metadata = NULL,
                                      auto_discover_related = TRUE,
                                      max_related = 5L,
                                      min_related_completeness = 95,
                                      min_related_ani = 90,
                                      max_target_sites_per_family = 10000,
                                      site_scan_strategy = c("balanced", "full"),
                                      max_landing_pads = 500L,
                                      verbose = FALSE) {
  detection_mode <- base::match.arg(detection_mode)
  prodigal_mode <- base::match.arg(prodigal_mode)
  analysis_depth <- base::match.arg(analysis_depth)
  site_scan_strategy <- base::match.arg(site_scan_strategy)
  genbank <- genbank %||% .dnmb_module_detect_genbank(base::getwd())
  if (base::is.null(genbank) || !base::file.exists(genbank)) {
    return(list(
      ok = FALSE,
      status = .dnmb_iselement_status_row("iselement_genbank", "missing", "GenBank file is required for IS element analysis."),
      files = list(),
      elements = tibble::tibble(),
      hits = .dnmb_module_empty_optional_long_table(),
      output_table = data.frame()
    ))
  }

  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  parsed <- .dnmb_parse_genbank_features(genbank)
  gene_table <- .dnmb_build_gene_table(parsed$features)

  annotation_elements <- .dnmb_detect_is_elements(parsed$features, gene_table)
  sequence_engine <- .dnmb_empty_sequence_engine()
  if (identical(detection_mode, "hybrid")) {
    sequence_engine <- .dnmb_run_sequence_engine(
      parsed = parsed,
      output_dir = output_dir,
      prodigal_mode = prodigal_mode,
      isescan_dir = isescan_dir,
      isfinder_db = isfinder_db,
      isfinder_fasta = isfinder_fasta,
      use_genbank_proteins = TRUE,
      verbose = verbose
    )
  }

  # If sequence engine produced results, filter annotation-only elements:
  # Keep annotation elements only if they overlap with HMM-verified elements
  # or have high-confidence annotation (mobile_element feature type)
  seq_el <- sequence_engine$sequence_elements
  if (nrow(seq_el) && nrow(annotation_elements)) {
    .dnmb_iselement_verbose(verbose, "Filtering annotation elements against HMM results")
    keep <- vapply(seq_len(nrow(annotation_elements)), function(i) {
      ae <- annotation_elements[i, , drop = FALSE]
      # Keep if it's a directly annotated mobile_element feature
      if (!is.na(ae$evidence_label) && identical(ae$evidence_label, "annotated_mobile_element")) return(TRUE)
      # Keep if overlaps with any sequence engine element
      ae_contig <- ae$contig
      ae_start <- ae$start; ae_end <- ae$end
      if (is.na(ae_contig) || is.na(ae_start) || is.na(ae_end)) return(FALSE)
      any(seq_el$contig == ae_contig &
          seq_el$start <= ae_end &
          seq_el$end >= ae_start, na.rm = TRUE)
    }, logical(1))
    keep[is.na(keep)] <- FALSE
    n_removed <- sum(!keep)
    annotation_elements <- annotation_elements[keep, , drop = FALSE]
    if (n_removed > 0L) {
      .dnmb_iselement_verbose(verbose, paste0("Removed ", n_removed, " annotation-only false positives"))
    }
  }

  elements <- .dnmb_merge_element_calls(
    annotation_elements = annotation_elements,
    sequence_elements = seq_el,
    genes = gene_table
  )
  # Filter out "new" family (ISEScan unclassified cluster — XerC/XerD recombinases, not true IS)
  if (nrow(elements) && "element_family" %in% names(elements)) {
    n_new <- sum(elements$element_family == "new", na.rm = TRUE)
    if (n_new > 0) {
      elements <- elements[elements$element_family != "new" | is.na(elements$element_family), , drop = FALSE]
      .dnmb_iselement_verbose(verbose, paste0("Filtered ", n_new, " 'new' (unclassified) elements"))
    }
  }
  hits <- .dnmb_iselement_normalize_hits(elements)
  output_table <- .dnmb_iselement_output_table(genes = genes, hits = hits)

  element_file <- base::file.path(output_dir, "iselement_elements.tsv")
  if (base::nrow(elements)) {
    utils::write.table(elements, file = element_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  status <- dplyr::bind_rows(
    .dnmb_iselement_empty_status(),
    .dnmb_iselement_status_row("iselement_annotation", "ok", base::paste0("annotation_elements=", nrow(annotation_elements))),
    .dnmb_iselement_status_row("iselement_sequence", if (base::nrow(sequence_engine$sequence_elements)) "ok" else "skipped", base::paste0("sequence_elements=", nrow(sequence_engine$sequence_elements))),
    .dnmb_iselement_status_row("iselement_merge", "ok", base::paste0("merged_elements=", nrow(elements)))
  )

  result <- list(
    ok = TRUE,
    status = status,
    files = list(element_table = element_file),
    elements = elements,
    hits = hits,
    output_table = output_table,
    sequence_engine = sequence_engine,
    census = NULL,
    target_models = NULL,
    landing_pads = NULL,
    comparative = NULL
  )

  # --- Full analysis depth: deep IS analysis with landing pad scoring ---
  if (identical(analysis_depth, "full") && nrow(elements)) {
    .dnmb_iselement_verbose(verbose, "Scoring gene essentiality")
    gene_table <- .dnmb_predict_gene_essentiality(gene_table)

    .dnmb_iselement_verbose(verbose, "Building IS target recognition models")
    target_models <- .dnmb_build_target_models(
      elements = elements,
      metadata = parsed$metadata,
      genes = gene_table
    )

    # Collect raw TSD sequences per family for seqlogo
    .dnmb_iselement_verbose(verbose, "Collecting TSD sequences for recognition analysis")
    seq_map <- stats::setNames(parsed$metadata$sequence, parsed$metadata$contig)
    rules <- .dnmb_family_variant_rules()
    tsd_raw_list <- list()
    for (fam_name in unique(elements$element_family)) {
      if (is.na(fam_name) || !nzchar(fam_name)) next
      fam_el <- elements[elements$element_family == fam_name, , drop = FALSE]
      # Filter to contigs we have sequences for
      fam_el <- fam_el[fam_el$contig %in% names(seq_map), , drop = FALSE]
      if (!nrow(fam_el)) next
      rule <- rules[rules$family == fam_name, , drop = FALSE]
      exp_len <- if (nrow(rule)) rule$dr_expected_bp[1] else NA_integer_
      tsd_rows <- tryCatch(
        .dnmb_collect_empirical_tsd_signals(fam_el, seq_map, exp_len),
        error = function(e) tibble::tibble()
      )
      if (nrow(tsd_rows)) {
        tsd_rows$family <- fam_name
        tsd_raw_list[[length(tsd_raw_list) + 1]] <- tsd_rows
      }
    }
    tsd_raw <- dplyr::bind_rows(tsd_raw_list)
    if (nrow(tsd_raw)) {
      tsd_raw_file <- file.path(output_dir, "iselement_tsd_sequences.tsv")
      utils::write.table(tsd_raw, file = tsd_raw_file, sep = "\t", row.names = FALSE, quote = FALSE)
    }

    target_sites <- tibble::tibble()
    targetable_regions <- tibble::tibble()
    if (nrow(target_models)) {
      .dnmb_iselement_verbose(verbose, "Scanning genome for IS target sites")
      target_sites <- .dnmb_predict_target_sites(
        target_models = target_models,
        metadata = parsed$metadata,
        genes = gene_table,
        max_target_sites_per_family = max_target_sites_per_family,
        site_scan_strategy = site_scan_strategy
      )

      .dnmb_iselement_verbose(verbose, "Aggregating targetable regions")
      targetable_regions <- .dnmb_build_targetable_regions(
        target_models = target_models,
        target_sites = target_sites,
        metadata = parsed$metadata,
        genes = gene_table,
        site_scan_strategy = site_scan_strategy
      )
    }

    comparative <- NULL
    # Auto-discover related genomes if none provided
    if ((is.null(related_genbanks) || !length(related_genbanks)) && isTRUE(auto_discover_related)) {
      .dnmb_iselement_verbose(verbose, "Auto-discovering related genomes from NCBI (ANI-based)")
      auto_result <- tryCatch(
        .dnmb_auto_discover_related_genomes(
          genbank = genbank,
          output_dir = output_dir,
          max_related = max_related,
          min_ani = min_related_ani,
          verbose = verbose
        ),
        error = function(e) {
          if (isTRUE(verbose)) {
            message("[DNMB iselement] Auto-discovery failed: ", conditionMessage(e))
          }
          .dnmb_empty_auto_comparative()
        }
      )
      if (length(auto_result$related_genbanks)) {
        related_genbanks <- auto_result$related_genbanks
        related_metadata <- auto_result$related_metadata
        .dnmb_iselement_verbose(verbose, paste0(
          "Auto-discovered ", length(related_genbanks), " related genomes"
        ))
      }
    }
    if (!is.null(related_genbanks) && length(related_genbanks)) {
      .dnmb_iselement_verbose(verbose, "Running comparative mobileome analysis")
      focal_genome <- list(
        genome_id = "focal",
        genbank = genbank,
        metadata = parsed$metadata,
        features = parsed$features,
        genes = gene_table,
        elements = elements,
        proteins = gene_table %>%
          dplyr::filter(!is.na(.data$translation), nzchar(.data$translation)) %>%
          dplyr::transmute(
            genome_id = "focal",
            gene_id = .data$gene_id,
            locus_tag = .data$locus_tag,
            protein_seq = .data$translation,
            protein_label = paste0("focal|", .data$gene_id)
          ),
        faa = file.path(output_dir, "comparative", "focal", "focal.faa")
      )
      dir.create(dirname(focal_genome$faa), recursive = TRUE, showWarnings = FALSE)
      .dnmb_write_protein_fasta(focal_genome$proteins, focal_genome$faa)

      comparative <- tryCatch(
        .dnmb_run_comparative_mobileome(
          focal = focal_genome,
          related_genbanks = related_genbanks,
          related_metadata = related_metadata,
          min_related_completeness = min_related_completeness,
          min_related_ani = min_related_ani,
          output_dir = output_dir,
          detection_mode = "annotation",
          verbose = verbose
        ),
        error = function(e) {
          if (isTRUE(verbose)) {
            message("[DNMB iselement] Comparative analysis failed: ", conditionMessage(e))
          }
          NULL
        }
      )
    }

    # Cross-genome insertion verification (if related genomes available)
    insertion_verify <- .dnmb_empty_verification()
    if (!is.null(related_genbanks) && length(related_genbanks)) {
      .dnmb_iselement_verbose(verbose, "Verifying IS insertions against related genomes")
      insertion_verify <- tryCatch(
        .dnmb_verify_is_insertions(
          elements = elements,
          focal_metadata = parsed$metadata,
          related_genbanks = related_genbanks,
          output_dir = output_dir,
          flank_bp = 1000L,
          verbose = verbose
        ),
        error = function(e) {
          if (isTRUE(verbose)) message("[DNMB iselement] Insertion verify failed: ", conditionMessage(e))
          .dnmb_empty_verification()
        }
      )
      # Merge validated TSDs with boundary-detected TSDs
      if (nrow(insertion_verify$validated_tsd)) {
        validated_for_merge <- insertion_verify$validated_tsd %>%
          dplyr::select(family = .data$family, tsd_seq = .data$tsd_seq,
                        tsd_len_bp = .data$tsd_len_bp, source = .data$source)
        tsd_raw <- dplyr::bind_rows(tsd_raw, validated_for_merge)
        # Re-write combined TSD file
        tsd_raw_file <- file.path(output_dir, "iselement_tsd_sequences.tsv")
        utils::write.table(tsd_raw, file = tsd_raw_file, sep = "\t", row.names = FALSE, quote = FALSE)
        .dnmb_iselement_verbose(verbose, paste0(
          "Added ", nrow(insertion_verify$validated_tsd), " comparative-verified TSDs"
        ))
      }
    }

    .dnmb_iselement_verbose(verbose, "Scoring landing pads")
    landing_pads <- .dnmb_score_landing_pads(
      genes = gene_table,
      elements = elements,
      target_models = target_models,
      target_sites = target_sites,
      metadata = parsed$metadata,
      comparative = comparative,
      insertion_verify = insertion_verify,
      max_landing_pads = max_landing_pads
    )

    .dnmb_iselement_verbose(verbose, "Building IS element census")
    census <- .dnmb_build_is_census(
      elements = elements,
      target_models = target_models,
      landing_pads = landing_pads
    )

    # Write output files
    census_file <- file.path(output_dir, "iselement_census.tsv")
    if (nrow(census)) {
      utils::write.table(census, file = census_file, sep = "\t", row.names = FALSE, quote = FALSE)
    }

    landing_pads_file <- file.path(output_dir, "iselement_landing_pads.tsv")
    landing_pads_bed <- file.path(output_dir, "iselement_landing_pads.bed")
    if (nrow(landing_pads)) {
      # Ensure numeric columns before arithmetic
      for (.lp_col in c("region_start", "region_end", "landing_pad_score")) {
        if (.lp_col %in% names(landing_pads) && !is.numeric(landing_pads[[.lp_col]])) {
          landing_pads[[.lp_col]] <- as.numeric(as.character(landing_pads[[.lp_col]]))
        }
      }
      utils::write.table(landing_pads, file = landing_pads_file, sep = "\t", row.names = FALSE, quote = FALSE)
      bed_tbl <- landing_pads %>%
        dplyr::transmute(
          chrom = .data$contig,
          chromStart = pmax(0L, .data$region_start - 1L),
          chromEnd = .data$region_end,
          name = .data$landing_pad_id,
          score = as.integer(round(1000 * pmin(1, .data$landing_pad_score))),
          strand = "."
        )
      utils::write.table(bed_tbl, file = landing_pads_bed, sep = "\t",
                         row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    target_models_file <- file.path(output_dir, "iselement_target_models.tsv")
    if (nrow(target_models)) {
      utils::write.table(target_models, file = target_models_file, sep = "\t", row.names = FALSE, quote = FALSE)
    }

    # Enrich output_table with essentiality and landing pad info
    if (nrow(output_table) && "essentiality_class" %in% names(gene_table)) {
      ess_map <- stats::setNames(gene_table$essentiality_class, gene_table$locus_tag)
      locus_col <- if ("locus_tag" %in% names(output_table)) output_table$locus_tag else rownames(output_table)
      output_table$iselement_essentiality_class <- ess_map[locus_col]
    }
    if (nrow(landing_pads) && nrow(output_table)) {
      output_table$iselement_nearest_landing_pad_rank <- .dnmb_nearest_landing_pad_rank(output_table, landing_pads, genes)
    }

    # Update result
    result$output_table <- output_table
    result$census <- census
    result$target_models <- target_models
    result$target_sites <- target_sites
    result$targetable_regions <- targetable_regions
    result$landing_pads <- landing_pads
    result$comparative <- comparative
    result$tsd_raw <- tsd_raw
    result$insertion_verify <- insertion_verify
    result$files$census_table <- census_file
    result$files$landing_pads_table <- landing_pads_file
    result$files$landing_pads_bed <- landing_pads_bed
    result$files$target_models_table <- target_models_file

    result$status <- dplyr::bind_rows(
      result$status,
      .dnmb_iselement_status_row("iselement_target_models", "ok", paste0("models=", nrow(target_models))),
      .dnmb_iselement_status_row("iselement_target_sites", "ok", paste0("sites=", nrow(target_sites))),
      .dnmb_iselement_status_row("iselement_landing_pads", "ok", paste0("pads=", nrow(landing_pads), " high=", sum(landing_pads$landing_pad_confidence == "high", na.rm = TRUE))),
      .dnmb_iselement_status_row("iselement_census", "ok", paste0("families=", nrow(census))),
      .dnmb_iselement_status_row("iselement_comparative",
                                  if (!is.null(comparative)) "ok" else "skipped",
                                  if (!is.null(comparative)) paste0("hotspots=", nrow(comparative$hotspots %||% tibble::tibble())) else "no related genbanks supplied")
    )

    .dnmb_iselement_verbose(verbose, paste0(
      "IS element deep analysis complete: ",
      nrow(elements), " elements, ",
      nrow(census), " families, ",
      nrow(landing_pads), " landing pads (",
      sum(landing_pads$landing_pad_confidence == "high", na.rm = TRUE), " high)"
    ))
  }

  result
}

# ---------------------------------------------------------------------------
# IS Census Summary
# ---------------------------------------------------------------------------
.dnmb_build_is_census <- function(elements, target_models, landing_pads) {
  # Ensure start/end are numeric (may arrive as factor/character from TSV)
  for (col in c("start", "end")) {
    if (col %in% names(elements) && !is.numeric(elements[[col]])) {
      elements[[col]] <- as.numeric(as.character(elements[[col]]))
    }
  }
  if (!nrow(elements)) {
    return(tibble::tibble(
      family = character(), element_count = integer(),
      high_confidence_count = integer(), medium_confidence_count = integer(),
      low_confidence_count = integer(),
      annotation_supported_count = integer(), sequence_supported_count = integer(),
      mean_element_length_bp = numeric(), contig_distribution = character(),
      recognition_motif = character(), recognition_strategy = character(),
      model_confidence = character(), dominant_tsd_len_bp = integer(),
      observed_tsd_examples = character(),
      landing_pad_count_high = integer(), landing_pad_count_medium = integer()
    ))
  }

  family_col <- if ("element_family" %in% names(elements)) "element_family" else "family"
  families <- elements[[family_col]]
  families[is.na(families) | !nzchar(families)] <- "unknown"
  elements$census_family <- families

  summary_tbl <- elements %>%
    dplyr::group_by(.data$census_family) %>%
    dplyr::summarise(
      element_count = dplyr::n(),
      high_confidence_count = sum(.data$confidence == "high", na.rm = TRUE),
      medium_confidence_count = sum(.data$confidence == "medium", na.rm = TRUE),
      low_confidence_count = sum(.data$confidence == "low", na.rm = TRUE),
      annotation_supported_count = sum(grepl("annotation", .data$evidence_summary, ignore.case = TRUE), na.rm = TRUE),
      sequence_supported_count = sum(grepl("hmm|sequence|phmmer", .data$evidence_summary, ignore.case = TRUE), na.rm = TRUE),
      mean_element_length_bp = round(mean(as.numeric(as.character(.data$end)) - as.numeric(as.character(.data$start)) + 1L, na.rm = TRUE), 0),
      contig_distribution = paste(sort(unique(.data$contig)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::rename(family = .data$census_family)

  # Join target model info
  if (!is.null(target_models) && nrow(target_models)) {
    model_info <- target_models %>%
      dplyr::transmute(
        family = .data$family,
        recognition_motif = dplyr::coalesce(.data$model_motif, NA_character_),
        recognition_strategy = .data$recognition_strategy,
        model_confidence = .data$model_confidence,
        dominant_tsd_len_bp = .data$dominant_tsd_len,
        observed_tsd_examples = .data$observed_tsd_examples
      )
    summary_tbl <- summary_tbl %>%
      dplyr::left_join(model_info, by = "family")
  } else {
    summary_tbl$recognition_motif <- NA_character_
    summary_tbl$recognition_strategy <- NA_character_
    summary_tbl$model_confidence <- NA_character_
    summary_tbl$dominant_tsd_len_bp <- NA_integer_
    summary_tbl$observed_tsd_examples <- NA_character_
  }

  # Join landing pad counts by matching IS families
  if (!is.null(landing_pads) && nrow(landing_pads) &&
      "matching_is_families" %in% names(landing_pads)) {
    lp_family_counts <- .dnmb_count_landing_pads_by_family(landing_pads, unique(summary_tbl$family))
    summary_tbl <- summary_tbl %>%
      dplyr::left_join(lp_family_counts, by = "family")
  } else {
    summary_tbl$landing_pad_count_high <- 0L
    summary_tbl$landing_pad_count_medium <- 0L
  }

  summary_tbl %>%
    dplyr::arrange(dplyr::desc(.data$element_count))
}

.dnmb_count_landing_pads_by_family <- function(landing_pads, families) {
  # Count landing pads that are NOT at risk from each family (i.e., safe from that family)
  high_pads <- landing_pads %>%
    dplyr::filter(.data$landing_pad_confidence == "high")
  medium_pads <- landing_pads %>%
    dplyr::filter(.data$landing_pad_confidence == "medium")

  tibble::tibble(
    family = families,
    landing_pad_count_high = vapply(families, function(f) {
      sum(is.na(high_pads$matching_is_families) |
          !grepl(f, high_pads$matching_is_families, fixed = TRUE), na.rm = TRUE)
    }, integer(1)),
    landing_pad_count_medium = vapply(families, function(f) {
      sum(is.na(medium_pads$matching_is_families) |
          !grepl(f, medium_pads$matching_is_families, fixed = TRUE), na.rm = TRUE)
    }, integer(1))
  )
}

# ---------------------------------------------------------------------------
# Nearest landing pad rank per gene
# ---------------------------------------------------------------------------
.dnmb_nearest_landing_pad_rank <- function(output_table, landing_pads, genes) {
  if (!nrow(landing_pads) || !nrow(output_table)) {
    return(rep(NA_integer_, nrow(output_table)))
  }
  output_table <- as.data.frame(output_table, stringsAsFactors = FALSE)
  landing_pads <- as.data.frame(landing_pads, stringsAsFactors = FALSE)
  for (col in c("start", "end")) {
    if (col %in% names(output_table) && !is.numeric(output_table[[col]])) {
      output_table[[col]] <- suppressWarnings(as.numeric(output_table[[col]]))
    }
  }
  for (col in c("region_start", "region_end")) {
    if (col %in% names(landing_pads) && !is.numeric(landing_pads[[col]])) {
      landing_pads[[col]] <- suppressWarnings(as.numeric(landing_pads[[col]]))
    }
  }
  backbone_cols <- c("locus_tag", "contig", "start", "end")
  has_coords <- all(c("contig", "start", "end") %in% names(output_table))
  if (!has_coords && "locus_tag" %in% names(output_table) && nrow(genes)) {
    coord_cols <- intersect(c("locus_tag", "contig", "start", "end"), names(genes))
    coord_map <- as.data.frame(genes[, coord_cols, drop = FALSE], stringsAsFactors = FALSE)
    output_table <- dplyr::left_join(output_table, coord_map, by = "locus_tag")
  }
  if (!all(c("contig", "start", "end") %in% names(output_table))) {
    return(rep(NA_integer_, nrow(output_table)))
  }

  vapply(seq_len(nrow(output_table)), function(i) {
    ctg <- output_table$contig[[i]]
    mid <- (output_table$start[[i]] + output_table$end[[i]]) / 2
    lp_ctg <- landing_pads[landing_pads$contig == ctg, , drop = FALSE]
    if (!nrow(lp_ctg)) return(NA_integer_)
    lp_ctg$dist <- pmin(abs(lp_ctg$region_start - mid), abs(lp_ctg$region_end - mid))
    nearest_idx <- which.min(lp_ctg$dist)
    if (!length(nearest_idx)) return(NA_integer_)
    lp_id <- lp_ctg$landing_pad_id[[nearest_idx]]
    match(lp_id, landing_pads$landing_pad_id)
  }, integer(1))
}

# ---------------------------------------------------------------------------
# Verbose helper
# ---------------------------------------------------------------------------
.dnmb_iselement_verbose <- function(verbose, msg) {
  if (isTRUE(verbose)) {
    message("[DNMB iselement] ", msg)
  }
}
