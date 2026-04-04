# Mobileome_landing_pads.R — Landing pad scoring for safe genome insertion sites
#
# 2-Phase landing pad selection algorithm:
#   Phase 1: Hard filters remove biologically unsafe candidates
#   Phase 2: Hybrid multiplicative-additive scoring ranks survivors
#
# final_score = safety_gate * quality_score
#   safety_gate  = essentiality_gate * is_distance_gate * operon_gate
#   quality_score = weighted sum of 7 sub-scores (comparative conservation,
#                   metabolic neutrality, region size, IS density,
#                   redundancy buffer, transcriptional isolation, GC normality)

.dnmb_score_landing_pads <- function(
  genes,
  elements,
  target_models,
  target_sites,
  metadata,
  comparative = NULL,
  insertion_verify = NULL,
  max_landing_pads = 500L
) {
  if (!nrow(genes) || !nrow(elements)) {
    return(.dnmb_empty_landing_pads())
  }

  contigs <- unique(genes$contig)
  contigs <- contigs[!is.na(contigs) & nzchar(contigs)]
  if (!length(contigs)) {
    return(.dnmb_empty_landing_pads())
  }

  candidate_list <- lapply(contigs, function(ctg) {
    .dnmb_build_landing_pad_candidates_contig(
      contig = ctg,
      genes = genes,
      metadata = metadata
    )
  })
  candidates <- dplyr::bind_rows(candidate_list)
  if (!nrow(candidates)) {
    return(.dnmb_empty_landing_pads())
  }

  # ===================================================================
  # Phase 1: Hard filters — remove biologically unsafe candidates
  # ===================================================================
  filter_result <- .dnmb_hard_filter_landing_pads(candidates, genes, elements)
  candidates <- filter_result$passed
  filtered_out <- filter_result$filtered
  if (!nrow(candidates)) {
    return(.dnmb_empty_landing_pads())
  }

  # ===================================================================
  # Phase 2: Hybrid multiplicative-additive scoring
  # ===================================================================

  # --- Safety gates (multiplicative, range 0-1) ---
  candidates$essentiality_gate <- .dnmb_score_essentiality_gate(candidates, genes)
  candidates$is_distance_gate <- .dnmb_score_is_distance_gate(candidates, elements)
  candidates$operon_gate <- .dnmb_score_operon_gate(candidates)

  candidates$safety_gate <- (
    candidates$essentiality_gate *
    candidates$is_distance_gate *
    candidates$operon_gate
  )

  # --- Quality sub-scores (additive, range 0-1) ---

  # Sub-score 1: Comparative conservation (weight 0.25)
  candidates$comparative_conservation <- .dnmb_score_comparative_conservation(
    candidates = candidates,
    comparative = comparative,
    insertion_verify = insertion_verify,
    elements = elements
  )

  # Sub-score 2: Metabolic neutrality (weight 0.20)
  candidates$metabolic_neutrality <- .dnmb_score_metabolic_neutrality(candidates)

  # Sub-score 3: Region size quality (weight 0.15)
  candidates$region_size_score <- .dnmb_score_region_size(candidates)

  # Sub-score 4: IS density penalty (weight 0.10)
  candidates$is_distance_score <- .dnmb_score_is_distance(candidates, elements)

  # Sub-score 5: IS target site vulnerability (weight 0.10) — construct stability risk
  tsv_result <- .dnmb_score_target_site_vulnerability(candidates, target_sites, target_models)
  candidates$target_site_vulnerability <- tsv_result$scores
  candidates$matching_is_families <- tsv_result$families
  candidates$target_site_context <- tsv_result$context

  # Sub-score 6: Redundancy buffer (weight 0.10)
  candidates$redundancy_buffer <- .dnmb_score_redundancy_buffer(candidates, genes)

  # Sub-score 7: Transcriptional isolation (weight 0.10)
  candidates$transcriptional_isolation <- .dnmb_score_transcriptional_isolation(candidates)

  # Sub-score 8: GC content normality (weight 0.05)
  candidates$gc_normality <- .dnmb_score_gc_normality(candidates, metadata)

  # --- Quality score (weighted additive, 8 sub-scores) ---
  candidates$quality_score <- (
    0.20 * candidates$comparative_conservation +
    0.20 * candidates$metabolic_neutrality +
    0.15 * candidates$region_size_score +
    0.10 * candidates$is_distance_score +
    0.10 * candidates$target_site_vulnerability +
    0.10 * candidates$redundancy_buffer +
    0.10 * candidates$transcriptional_isolation +
    0.05 * candidates$gc_normality
  )

  # --- Final score = safety_gate * quality_score ---
  candidates$landing_pad_score <- candidates$safety_gate * candidates$quality_score

  # --- Confidence classification ---
  # Confidence thresholds adjusted for 2-phase algorithm (safety_gate × quality lowers scores)
  candidates$landing_pad_confidence <- dplyr::case_when(
    candidates$landing_pad_score >= 0.60 & candidates$safety_gate >= 0.7 ~ "high",
    candidates$landing_pad_score >= 0.35 ~ "medium",
    TRUE ~ "low"
  )

  # --- Genome-wide spacing: select diverse representatives ---
  candidates <- .dnmb_apply_genome_spacing(candidates, max_landing_pads)

  # --- Select and order output columns ---
  candidates %>%
    dplyr::select(
      .data$landing_pad_id,
      .data$contig,
      .data$region_start,
      .data$region_end,
      region_size_bp = .data$region_size,
      .data$landing_pad_score,
      .data$landing_pad_confidence,
      .data$safety_gate,
      .data$quality_score,
      .data$essentiality_gate,
      .data$is_distance_gate,
      .data$operon_gate,
      .data$comparative_conservation,
      .data$metabolic_neutrality,
      .data$region_size_score,
      .data$is_distance_score,
      .data$target_site_vulnerability,
      .data$target_site_context,
      .data$redundancy_buffer,
      .data$transcriptional_isolation,
      .data$gc_normality,
      .data$context_class,
      .data$gene_orientation,
      .data$left_gene_id,
      .data$left_gene_product,
      .data$left_gene_essentiality_class,
      .data$left_gene_strand,
      .data$right_gene_id,
      .data$right_gene_product,
      .data$right_gene_essentiality_class,
      .data$right_gene_strand,
      dplyr::any_of(c("is_comparative_hotspot", "comparative_occupancy_ratio", "matching_is_families"))
    )
}

# Keep backward-compatible essentiality_safety for any external callers
.dnmb_score_essentiality_safety <- function(candidates) {
  left_ess <- dplyr::coalesce(candidates$left_gene_essentiality_score, 0.05)
  right_ess <- dplyr::coalesce(candidates$right_gene_essentiality_score, 0.05)
  max_ess <- pmax(left_ess, right_ess)
  pmax(0.0, pmin(1.0, 1.0 - max_ess))
}

.dnmb_empty_landing_pads <- function() {
  tibble::tibble(
    landing_pad_id = character(),
    contig = character(),
    region_start = integer(),
    region_end = integer(),
    region_size_bp = integer(),
    landing_pad_score = numeric(),
    landing_pad_confidence = character(),
    safety_gate = numeric(),
    quality_score = numeric(),
    essentiality_gate = numeric(),
    is_distance_gate = numeric(),
    operon_gate = numeric(),
    comparative_conservation = numeric(),
    metabolic_neutrality = numeric(),
    region_size_score = numeric(),
    is_distance_score = numeric(),
    target_site_vulnerability = numeric(),
    target_site_context = character(),
    redundancy_buffer = numeric(),
    transcriptional_isolation = numeric(),
    gc_normality = numeric(),
    context_class = character(),
    gene_orientation = character(),
    left_gene_id = character(),
    left_gene_product = character(),
    left_gene_essentiality_class = character(),
    left_gene_strand = character(),
    right_gene_id = character(),
    right_gene_product = character(),
    right_gene_essentiality_class = character(),
    right_gene_strand = character()
  )
}

# ===========================================================================
# Phase 1: Hard filters — remove biologically unsafe candidates
# ===========================================================================
.dnmb_hard_filter_landing_pads <- function(candidates, genes, elements) {
  n_total <- nrow(candidates)
  reasons <- rep(NA_character_, n_total)

  # Filter 1: Region size < 100bp — too small for homologous recombination
  too_small <- candidates$region_size < 100L
  reasons[too_small & is.na(reasons)] <- "region_size_lt_100bp"

  # Filter 2: Flanking gene essentiality_class == "high" — direct disruption risk
  left_high <- dplyr::coalesce(candidates$left_gene_essentiality_class, "low") == "high"
  right_high <- dplyr::coalesce(candidates$right_gene_essentiality_class, "low") == "high"
  flanking_essential <- left_high | right_high
  reasons[flanking_essential & is.na(reasons)] <- "flanking_gene_high_essentiality"

  # Filter 3: Any gene with essentiality_class == "high" within 500bp of boundary
  # Check ALL genes, not just immediate flanking genes
  has_ess_col <- "essentiality_class" %in% names(genes)
  if (has_ess_col) {
    high_genes <- genes[dplyr::coalesce(genes$essentiality_class, "low") == "high", , drop = FALSE]
    if (nrow(high_genes)) {
      nearby_essential <- vapply(seq_len(n_total), function(i) {
        if (!is.na(reasons[i])) return(FALSE)  # already filtered
        cand_ctg <- candidates$contig[i]
        cand_start <- candidates$region_start[i]
        cand_end <- candidates$region_end[i]
        ctg_high <- high_genes[high_genes$contig == cand_ctg, , drop = FALSE]
        if (!nrow(ctg_high)) return(FALSE)
        # Check if any high-essentiality gene is within 500bp of boundary
        any(
          (ctg_high$end >= (cand_start - 500L) & ctg_high$start <= (cand_start + 500L)) |
          (ctg_high$end >= (cand_end - 500L) & ctg_high$start <= (cand_end + 500L))
        )
      }, logical(1))
      reasons[nearby_essential & is.na(reasons)] <- "high_essentiality_within_500bp"
    }
  }

  # Filter 4: Nearest IS element within 2kb — recombination hotspot risk
  if (nrow(elements)) {
    is_too_close <- vapply(seq_len(n_total), function(i) {
      if (!is.na(reasons[i])) return(FALSE)
      mid <- (candidates$region_start[i] + candidates$region_end[i]) / 2
      ctg_el <- elements[elements$contig == candidates$contig[i], , drop = FALSE]
      if (!nrow(ctg_el)) return(FALSE)
      dists <- pmin(abs(ctg_el$start - mid), abs(ctg_el$end - mid))
      min(dists, na.rm = TRUE) < 2000
    }, logical(1))
    reasons[is_too_close & is.na(reasons)] <- "is_element_within_2kb"
  }

  # Filter 5: >=3 IS elements within 10kb — genomic island / unstable region
  if (nrow(elements)) {
    is_dense <- vapply(seq_len(n_total), function(i) {
      if (!is.na(reasons[i])) return(FALSE)
      mid <- (candidates$region_start[i] + candidates$region_end[i]) / 2
      ctg_el <- elements[elements$contig == candidates$contig[i], , drop = FALSE]
      if (!nrow(ctg_el)) return(FALSE)
      n_10kb <- sum(ctg_el$start <= (mid + 10000) & ctg_el$end >= (mid - 10000))
      n_10kb >= 3L
    }, logical(1))
    reasons[is_dense & is.na(reasons)] <- "ge3_is_elements_within_10kb"
  }

  # Filter 6: Operon interior — both flanking genes same strand AND
  #            intergenic distance < 100bp — polar effect will disrupt downstream
  if ("gene_orientation" %in% names(candidates)) {
    operon_interior <- candidates$gene_orientation == "co-directional" &
                       candidates$region_size < 100L
    reasons[operon_interior & is.na(reasons)] <- "operon_interior"
  }

  is_filtered <- !is.na(reasons)
  candidates$filter_reason <- reasons

  list(
    passed = candidates[!is_filtered, , drop = FALSE],
    filtered = candidates[is_filtered, , drop = FALSE],
    n_total = n_total,
    n_passed = sum(!is_filtered),
    n_filtered = sum(is_filtered)
  )
}

# ===========================================================================
# Build candidate intergenic regions per contig
# ===========================================================================
.dnmb_build_landing_pad_candidates_contig <- function(contig, genes, metadata) {
  ctg_genes <- genes %>%
    dplyr::filter(.data$contig == .env$contig) %>%
    dplyr::arrange(.data$start, .data$end)
  if (nrow(ctg_genes) < 2L) {
    return(tibble::tibble())
  }

  gap_starts <- ctg_genes$end[-nrow(ctg_genes)] + 1L
  gap_ends <- ctg_genes$start[-1L] - 1L
  valid <- gap_ends >= gap_starts & (gap_ends - gap_starts + 1L) >= 50L
  if (!any(valid)) {
    return(tibble::tibble())
  }

  left_genes <- ctg_genes[-nrow(ctg_genes), , drop = FALSE][valid, , drop = FALSE]
  right_genes <- ctg_genes[-1L, , drop = FALSE][valid, , drop = FALSE]

  # Determine strand info for flanking genes
  left_strand <- if ("strand" %in% names(left_genes)) {
    dplyr::coalesce(left_genes$strand, "+")
  } else {
    rep("+", nrow(left_genes))
  }
  right_strand <- if ("strand" %in% names(right_genes)) {
    dplyr::coalesce(right_genes$strand, "+")
  } else {
    rep("+", nrow(right_genes))
  }

  # Classify gene orientation
  gene_orientation <- dplyr::case_when(
    left_strand == "+" & right_strand == "-" ~ "divergent",
    left_strand == "-" & right_strand == "+" ~ "convergent",
    TRUE ~ "co-directional"
  )

  tibble::tibble(
    contig = contig,
    region_start = gap_starts[valid],
    region_end = gap_ends[valid],
    region_size = gap_ends[valid] - gap_starts[valid] + 1L,
    context_class = "intergenic",
    left_gene_id = left_genes$gene_id,
    left_gene = dplyr::coalesce(left_genes$gene, ""),
    left_gene_product = dplyr::coalesce(left_genes$product, ""),
    left_gene_essentiality_score = dplyr::coalesce(left_genes$essentiality_score, 0.05),
    left_gene_essentiality_class = dplyr::coalesce(left_genes$essentiality_class, "low"),
    left_gene_strand = left_strand,
    right_gene_id = right_genes$gene_id,
    right_gene = dplyr::coalesce(right_genes$gene, ""),
    right_gene_product = dplyr::coalesce(right_genes$product, ""),
    right_gene_essentiality_score = dplyr::coalesce(right_genes$essentiality_score, 0.05),
    right_gene_essentiality_class = dplyr::coalesce(right_genes$essentiality_class, "low"),
    right_gene_strand = right_strand,
    gene_orientation = gene_orientation
  )
}

# ===========================================================================
# Safety gate: Essentiality gate (multiplicative)
# Check all genes within 2kb; score = (1 - max_nearby_essentiality)^2
# ===========================================================================
.dnmb_score_essentiality_gate <- function(candidates, genes) {
  has_ess <- "essentiality_score" %in% names(genes)
  if (!has_ess) {
    # Fall back to flanking gene scores only
    left_ess <- dplyr::coalesce(candidates$left_gene_essentiality_score, 0.05)
    right_ess <- dplyr::coalesce(candidates$right_gene_essentiality_score, 0.05)
    max_ess <- pmax(left_ess, right_ess)
    return(pmax(0.0, pmin(1.0, (1.0 - max_ess)^2)))
  }

  vapply(seq_len(nrow(candidates)), function(i) {
    cand_ctg <- candidates$contig[i]
    cand_start <- candidates$region_start[i]
    cand_end <- candidates$region_end[i]
    mid <- (cand_start + cand_end) / 2

    ctg_genes <- genes[genes$contig == cand_ctg, , drop = FALSE]
    if (!nrow(ctg_genes)) return(1.0)

    # All genes within 2kb of the candidate region boundary
    nearby <- ctg_genes[
      (ctg_genes$end >= (cand_start - 2000L) & ctg_genes$start <= (cand_end + 2000L)),
      , drop = FALSE
    ]
    if (!nrow(nearby)) return(1.0)

    max_ess <- max(dplyr::coalesce(nearby$essentiality_score, 0.05), na.rm = TRUE)
    pmax(0.0, pmin(1.0, (1.0 - max_ess)^2))
  }, numeric(1))
}

# ===========================================================================
# Safety gate: IS distance gate (multiplicative)
# Linear ramp: 0 at 2kb, 1.0 at 10kb
# ===========================================================================
.dnmb_score_is_distance_gate <- function(candidates, elements) {
  if (!nrow(elements)) {
    return(rep(1.0, nrow(candidates)))
  }

  vapply(seq_len(nrow(candidates)), function(i) {
    mid <- (candidates$region_start[i] + candidates$region_end[i]) / 2
    ctg_el <- elements[elements$contig == candidates$contig[i], , drop = FALSE]
    if (!nrow(ctg_el)) return(1.0)

    dists <- pmin(abs(ctg_el$start - mid), abs(ctg_el$end - mid))
    min_is_dist <- min(dists, na.rm = TRUE)

    # Linear ramp: 0 at <=2kb, 1.0 at >=10kb
    pmin(1.0, pmax(0.0, (min_is_dist - 2000) / 8000))
  }, numeric(1))
}

# ===========================================================================
# Safety gate: Operon gate (multiplicative)
# Based on flanking gene orientation and intergenic distance
# ===========================================================================
.dnmb_score_operon_gate <- function(candidates) {
  orientation <- if ("gene_orientation" %in% names(candidates)) {
    candidates$gene_orientation
  } else {
    rep("co-directional", nrow(candidates))
  }
  region_size <- candidates$region_size

  dplyr::case_when(
    orientation == "divergent" ~ 1.0,
    orientation == "convergent" ~ 0.8,
    # co-directional with various intergenic distances
    orientation == "co-directional" & region_size > 300L ~ 0.7,
    orientation == "co-directional" & region_size >= 100L ~ 0.5,
    orientation == "co-directional" & region_size < 100L ~ 0.3,
    TRUE ~ 0.5
  )
}

# ===========================================================================
# Quality sub-score: Metabolic neutrality (weight 0.20)
# Keyword-based assessment of flanking gene metabolic importance
# ===========================================================================
.dnmb_score_metabolic_neutrality <- function(candidates) {
  # Critical metabolic keywords -> high penalty
  critical_kw <- c(
    "glycolysis", "tca", "citrate", "pyruvate", "phosphotransferase",
    "pentose phosphate", "electron transport", "nadh", "cytochrome",
    "atp synthase", "acetyl-coa", "ribosom"
  )
  # Neutral keywords -> safe
  neutral_kw <- c("hypothetical protein", "unknown function", "uncharacterized")
  # Moderate keywords -> moderate safety
  moderate_kw <- c("regulatory", "regulator", "signaling", "signal", "motility",
                    "chemotaxis", "transcriptional repressor", "transcriptional activator")

  vapply(seq_len(nrow(candidates)), function(i) {
    left_prod <- tolower(dplyr::coalesce(candidates$left_gene_product[i], ""))
    right_prod <- tolower(dplyr::coalesce(candidates$right_gene_product[i], ""))
    combined <- paste(left_prod, right_prod)

    # Check critical first
    if (any(vapply(critical_kw, function(kw) grepl(kw, combined, fixed = TRUE), logical(1)))) {
      return(0.2)
    }
    # Check neutral (both flanking are hypothetical/unknown)
    left_neutral <- any(vapply(neutral_kw, function(kw) grepl(kw, left_prod, fixed = TRUE), logical(1)))
    right_neutral <- any(vapply(neutral_kw, function(kw) grepl(kw, right_prod, fixed = TRUE), logical(1)))
    if (left_neutral && right_neutral) {
      return(1.0)
    }
    if (left_neutral || right_neutral) {
      # One side is neutral
      other <- if (left_neutral) right_prod else left_prod
      if (any(vapply(moderate_kw, function(kw) grepl(kw, other, fixed = TRUE), logical(1)))) {
        return(0.8)
      }
      return(0.9)
    }
    # Check moderate
    if (any(vapply(moderate_kw, function(kw) grepl(kw, combined, fixed = TRUE), logical(1)))) {
      return(0.7)
    }
    # Default: unknown function category -> moderate safety
    0.6
  }, numeric(1))
}

# ===========================================================================
# Quality sub-score: Transcriptional isolation (weight 0.10)
# ===========================================================================
.dnmb_score_transcriptional_isolation <- function(candidates) {
  orientation <- if ("gene_orientation" %in% names(candidates)) {
    candidates$gene_orientation
  } else {
    rep("co-directional", nrow(candidates))
  }

  dplyr::case_when(
    orientation == "divergent" ~ 1.0,
    orientation == "convergent" ~ 0.7,
    orientation == "co-directional" ~ 0.4,
    TRUE ~ 0.4
  )
}

# ===========================================================================
# Quality sub-score: GC content normality (weight 0.05)
# Compare local GC% to genome-wide mean
# ===========================================================================
.dnmb_score_gc_normality <- function(candidates, metadata) {
  # Build sequence map from metadata
  if (is.null(metadata) || !nrow(metadata) ||
      !"sequence" %in% names(metadata) || !"contig" %in% names(metadata)) {
    return(rep(0.7, nrow(candidates)))
  }

  seq_map <- stats::setNames(metadata$sequence, metadata$contig)

  # Compute genome-wide GC%
  all_seqs <- paste(metadata$sequence, collapse = "")
  total_len <- nchar(all_seqs)
  if (total_len == 0L) return(rep(0.7, nrow(candidates)))

  gc_count_total <- nchar(gsub("[^GCgc]", "", all_seqs))
  genome_gc <- gc_count_total / total_len

  # Compute genome-wide GC% standard deviation using sliding windows
  # Use 1kb windows across the genome
  window_gcs <- numeric()
  for (ctg in names(seq_map)) {
    seq <- seq_map[[ctg]]
    seq_len_bp <- nchar(seq)
    if (seq_len_bp < 1000L) next
    n_windows <- as.integer(seq_len_bp / 1000L)
    for (w in seq_len(min(n_windows, 100L))) {  # cap at 100 windows per contig for speed
      start_pos <- sample.int(seq_len_bp - 999L, 1L)
      window_seq <- substr(seq, start_pos, start_pos + 999L)
      gc_w <- nchar(gsub("[^GCgc]", "", window_seq)) / 1000
      window_gcs <- c(window_gcs, gc_w)
    }
  }
  gc_sd <- if (length(window_gcs) > 1L) stats::sd(window_gcs) else 0.05

  # Score each candidate
  vapply(seq_len(nrow(candidates)), function(i) {
    ctg <- candidates$contig[i]
    if (!ctg %in% names(seq_map)) return(0.7)
    seq <- seq_map[[ctg]]
    start_pos <- max(1L, candidates$region_start[i])
    end_pos <- min(nchar(seq), candidates$region_end[i])
    if (end_pos <= start_pos) return(0.7)

    region_seq <- substr(seq, start_pos, end_pos)
    region_len <- nchar(region_seq)
    if (region_len < 10L) return(0.7)

    local_gc <- nchar(gsub("[^GCgc]", "", region_seq)) / region_len
    deviation <- abs(local_gc - genome_gc) / max(gc_sd, 0.01)

    # Score: 1.0 if within 1SD, 0.7 if 1-2SD, 0.3 if >2SD
    if (deviation <= 1.0) {
      1.0
    } else if (deviation <= 2.0) {
      0.7
    } else {
      0.3
    }
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Sub-score: Comparative conservation (weight 0.25) — kept from original
# ---------------------------------------------------------------------------
.dnmb_score_comparative_conservation <- function(candidates, comparative,
                                                  insertion_verify = NULL,
                                                  elements = NULL) {
  has_comparative <- !is.null(comparative) && is.list(comparative) &&
    !is.null(comparative$hotspots) && nrow(comparative$hotspots) &&
    !is.null(comparative$occupied_empty) && nrow(comparative$occupied_empty)

  has_verify <- !is.null(insertion_verify) && is.list(insertion_verify) &&
    !is.null(insertion_verify$verification) && nrow(insertion_verify$verification) &&
    !is.null(elements) && nrow(elements)

  if (!has_comparative && !has_verify) {
    candidates$is_comparative_hotspot <- NA
    candidates$comparative_occupancy_ratio <- NA_real_
    return(rep(0.5, nrow(candidates)))
  }

  hotspots <- if (has_comparative) comparative$hotspots else NULL
  occupied_empty <- if (has_comparative) comparative$occupied_empty else NULL
  verify_df <- if (has_verify) insertion_verify$verification else NULL

  scores <- vapply(seq_len(nrow(candidates)), function(i) {
    cand <- candidates[i, , drop = FALSE]
    left_id <- cand$left_gene_id
    right_id <- cand$right_gene_id

    # --- Score from comparative hotspot data (original logic) ---
    comp_score <- NA_real_
    if (has_comparative) {
      matching_loci <- hotspots %>%
        dplyr::filter(
          (.data$left_anchor_gene_id == .env$left_id & .data$right_anchor_gene_id == .env$right_id) |
          (.data$left_anchor_gene_id == .env$right_id & .data$right_anchor_gene_id == .env$left_id)
        )
      if (nrow(matching_loci)) {
        locus_ids <- matching_loci$locus_id
        oe_rows <- occupied_empty %>%
          dplyr::filter(.data$locus_id %in% .env$locus_ids)
        if (nrow(oe_rows)) {
          n_total <- nrow(oe_rows)
          n_empty <- sum(oe_rows$state == "empty", na.rm = TRUE)
          n_occupied <- sum(grepl("occupied", oe_rows$state), na.rm = TRUE)
          empty_ratio <- n_empty / n_total

          if (n_occupied > 0) {
            comp_score <- 0.3
          } else if (empty_ratio >= 0.8) {
            comp_score <- 1.0
          } else {
            comp_score <- 0.5 + 0.5 * empty_ratio
          }
        }
      }
    }

    # --- Score from insertion verification data ---
    verify_score <- NA_real_
    if (has_verify) {
      cand_mid <- (cand$region_start + cand$region_end) / 2
      nearby_el <- elements[
        elements$contig == cand$contig &
        elements$start <= (cand_mid + 5000) &
        elements$end >= (cand_mid - 5000),
        , drop = FALSE
      ]
      if (nrow(nearby_el)) {
        nearby_ids <- nearby_el$element_id
        nearby_verify <- verify_df[verify_df$element_id %in% nearby_ids, , drop = FALSE]
        # Keep only conclusive statuses
        conclusive <- nearby_verify[nearby_verify$status %in% c("empty_site", "filled_site"), , drop = FALSE]
        if (nrow(conclusive)) {
          n_empty <- sum(conclusive$status == "empty_site")
          n_filled <- sum(conclusive$status == "filled_site")
          n_total <- n_empty + n_filled
          if (n_total > 0) {
            empty_frac <- n_empty / n_total
            # All empty_site -> 1.0 (safe); all filled_site -> 0.2 (dangerous)
            verify_score <- 0.2 + 0.8 * empty_frac
          }
        }
      }
    }

    # --- Combine: prefer comparative hotspot when available, otherwise use verify ---
    if (!is.na(comp_score) && !is.na(verify_score)) {
      return(0.6 * comp_score + 0.4 * verify_score)
    } else if (!is.na(comp_score)) {
      return(comp_score)
    } else if (!is.na(verify_score)) {
      return(verify_score)
    }
    return(0.5)  # neutral fallback
  }, numeric(1))

  candidates$is_comparative_hotspot <- scores <= 0.35
  candidates$comparative_occupancy_ratio <- scores
  scores
}

# ---------------------------------------------------------------------------
# Sub-score: IS target site vulnerability (construct stability risk)
# If IS target sites exist within/near the landing pad, a transposon could
# later insert into the construct, disrupting it. Higher score = safer.
# Also annotates target site context (none/internal/flanking_one/flanking_both)
# ---------------------------------------------------------------------------
.dnmb_score_target_site_vulnerability <- function(candidates, target_sites, target_models) {
  if (is.null(target_sites) || !nrow(target_sites) ||
      is.null(target_models) || !nrow(target_models)) {
    return(list(scores = rep(1.0, nrow(candidates)),
                families = rep(NA_character_, nrow(candidates)),
                context = rep("none", nrow(candidates))))
  }

  buffer_bp <- 500L
  scores <- numeric(nrow(candidates))
  families_vec <- character(nrow(candidates))
  context_vec <- character(nrow(candidates))

  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, , drop = FALSE]
    # Sites INSIDE the landing pad region
    internal <- target_sites[
      target_sites$contig == cand$contig &
      target_sites$site_start <= cand$region_end &
      target_sites$site_end >= cand$region_start, , drop = FALSE]
    # Sites in buffer zone (within 500bp outside the region)
    left_buffer <- target_sites[
      target_sites$contig == cand$contig &
      target_sites$site_end >= (cand$region_start - buffer_bp) &
      target_sites$site_start < cand$region_start, , drop = FALSE]
    right_buffer <- target_sites[
      target_sites$contig == cand$contig &
      target_sites$site_start <= (cand$region_end + buffer_bp) &
      target_sites$site_end > cand$region_end, , drop = FALSE]

    n_internal <- nrow(internal)
    n_left <- nrow(left_buffer)
    n_right <- nrow(right_buffer)
    all_nearby <- dplyr::bind_rows(internal, left_buffer, right_buffer)
    n_families <- if (nrow(all_nearby)) length(unique(all_nearby$target_family)) else 0L

    # Risk: internal sites weighted more heavily than buffer sites
    risk <- pmin(1.0, n_families * 0.20 + n_internal * 0.05 + (n_left + n_right) * 0.015)
    scores[i] <- 1.0 - risk

    # Context annotation
    context_vec[i] <- if (n_internal > 0) {
      "internal"
    } else if (n_left > 0 && n_right > 0) {
      "flanking_both"
    } else if (n_left > 0 || n_right > 0) {
      "flanking_one"
    } else {
      "none"
    }

    families_vec[i] <- if (nrow(all_nearby)) {
      paste(sort(unique(all_nearby$target_family)), collapse = "; ")
    } else NA_character_
  }

  list(scores = scores, families = families_vec, context = context_vec)
}

# ---------------------------------------------------------------------------
# Sub-score: Region size (larger intergenic regions are better landing pads)
# ---------------------------------------------------------------------------
.dnmb_score_region_size <- function(candidates) {
  sizes <- candidates$region_size
  dplyr::case_when(
    sizes >= 1000L ~ 1.0,
    sizes >= 500L ~ 0.7 + 0.3 * (sizes - 500) / 500,
    sizes >= 300L ~ 0.5 + 0.2 * (sizes - 300) / 200,
    sizes >= 100L ~ 0.3 + 0.2 * (sizes - 100) / 200,
    TRUE ~ 0.1 + 0.2 * sizes / 100
  )
}

# ---------------------------------------------------------------------------
# Sub-score: Distance from IS elements (further = safer/more stable)
# ---------------------------------------------------------------------------
.dnmb_score_is_distance <- function(candidates, elements) {
  if (!nrow(elements)) {
    return(rep(1.0, nrow(candidates)))
  }

  vapply(seq_len(nrow(candidates)), function(i) {
    cand <- candidates[i, , drop = FALSE]
    mid <- (cand$region_start + cand$region_end) / 2
    ctg_elements <- elements %>%
      dplyr::filter(.data$contig == cand$contig)
    if (!nrow(ctg_elements)) return(1.0)

    # Count IS elements within windows
    n_5kb <- sum(ctg_elements$start <= (mid + 5000) & ctg_elements$end >= (mid - 5000))
    n_20kb <- sum(ctg_elements$start <= (mid + 20000) & ctg_elements$end >= (mid - 20000))
    n_50kb <- sum(ctg_elements$start <= (mid + 50000) & ctg_elements$end >= (mid - 50000))

    # Min distance to any IS element
    dists <- pmin(abs(ctg_elements$start - mid), abs(ctg_elements$end - mid))
    min_dist <- min(dists, na.rm = TRUE)

    # Composite: penalize proximity and density
    dist_score <- pmin(1.0, min_dist / 10000)  # 0 at IS, 1.0 at >10kb
    density_penalty <- pmin(1.0, n_5kb * 0.15 + n_20kb * 0.03 + n_50kb * 0.005)
    pmax(0.0, dist_score * (1.0 - density_penalty))
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Sub-score: Redundancy buffer
# ---------------------------------------------------------------------------
.dnmb_score_redundancy_buffer <- function(candidates, genes) {
  if (!"redundancy_count" %in% names(genes)) {
    return(rep(0.25, nrow(candidates)))
  }
  rc_map <- stats::setNames(genes$redundancy_count, genes$gene_id)

  vapply(seq_len(nrow(candidates)), function(i) {
    left_rc <- rc_map[[candidates$left_gene_id[[i]]]]
    right_rc <- rc_map[[candidates$right_gene_id[[i]]]]
    if (is.null(left_rc)) left_rc <- 1L
    if (is.null(right_rc)) right_rc <- 1L
    max_rc <- max(left_rc, right_rc, na.rm = TRUE)
    pmin(1.0, (max_rc - 1) * 0.25 + 0.25)
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Genome-wide spacing: ensure landing pads are distributed across the genome
# ---------------------------------------------------------------------------
.dnmb_apply_genome_spacing <- function(candidates, max_landing_pads) {
  if (nrow(candidates) <= max_landing_pads) {
    candidates <- candidates %>%
      dplyr::arrange(dplyr::desc(.data$landing_pad_score)) %>%
      dplyr::mutate(landing_pad_id = paste0("LP_", sprintf("%04d", dplyr::row_number())))
    return(candidates)
  }

  # Sort by score first
  candidates <- candidates %>%
    dplyr::arrange(dplyr::desc(.data$landing_pad_score))

  # Genome-wide spacing: divide genome into bins, pick best from each bin,
  # then fill remaining slots by score
  contigs <- unique(candidates$contig)
  selected_idx <- integer()

  for (ctg in contigs) {
    ctg_mask <- candidates$contig == ctg
    ctg_rows <- which(ctg_mask)
    if (!length(ctg_rows)) next

    ctg_len <- max(candidates$region_end[ctg_rows], na.rm = TRUE)
    # Create bins of ~50kb each
    n_bins <- max(10L, as.integer(ctg_len / 50000))
    bin_size <- ctg_len / n_bins
    bins <- ceiling(candidates$region_start[ctg_rows] / bin_size)

    # Pick best candidate from each bin
    for (b in sort(unique(bins))) {
      bin_rows <- ctg_rows[bins == b]
      selected_idx <- c(selected_idx, bin_rows[1])  # already sorted by score
    }
  }

  # Remove duplicates
  selected_idx <- unique(selected_idx)

  # Fill remaining slots with next-best by score, maintaining minimum spacing
  remaining <- setdiff(seq_len(nrow(candidates)), selected_idx)
  min_spacing <- 5000L  # minimum 5kb between selected landing pads

  for (r in remaining) {
    if (length(selected_idx) >= max_landing_pads) break
    cand_mid <- (candidates$region_start[r] + candidates$region_end[r]) / 2
    cand_ctg <- candidates$contig[r]
    sel_in_ctg <- selected_idx[candidates$contig[selected_idx] == cand_ctg]
    if (length(sel_in_ctg)) {
      sel_mids <- (candidates$region_start[sel_in_ctg] + candidates$region_end[sel_in_ctg]) / 2
      if (min(abs(sel_mids - cand_mid)) < min_spacing) next
    }
    selected_idx <- c(selected_idx, r)
  }

  candidates <- candidates[selected_idx, , drop = FALSE] %>%
    dplyr::arrange(dplyr::desc(.data$landing_pad_score)) %>%
    dplyr::slice_head(n = as.integer(max_landing_pads)) %>%
    dplyr::mutate(landing_pad_id = paste0("LP_", sprintf("%04d", dplyr::row_number())))

  candidates
}
#' Internal helpers for mobileome landing-pad discovery
#'
#' Scoring and classification routines for identifying candidate landing pads
#' and targetable regions in mobileome analyses.
#'
#' @name dnmb_internal_mobileome_landing_pads
#' @keywords internal
#' @noRd
NULL
