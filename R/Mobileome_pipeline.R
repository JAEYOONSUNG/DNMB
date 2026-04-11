#' Run the DNMB mobileome impact pipeline
#'
#' This pipeline accepts a single GenBank file and performs a hybrid
#' mobile-element analysis focused on IS-like elements and their possible local
#' genome consequences. The workflow parses annotated features, optionally adds
#' sequence-driven evidence inspired by ISEScan (ORF calling plus pHMM/singleton
#' reference searches and family-aware TIR/TSD heuristics), merges annotation
#' and sequence evidence, predicts nearby mutation-prone genes, infers
#' deletion/inversion/transposition scenarios that could be mediated by element
#' pairs, and summarizes plausible phenotype changes from affected genes.
#'
#' The predictions are annotation-driven and are intended as a prioritization
#' layer for downstream manual review or experimental validation.
#'
#' @param genbank Path to a single `.gb`, `.gbk`, or `.gbff` file.
#' @param output_dir Directory where result tables will be written. If `NULL`,
#'   a directory named `<genbank_basename>_dnmb_mobileome` is created beside the
#'   input file.
#' @param related_genbanks Optional vector of related-strain GenBank paths used
#'   for comparative hotspot and chronology inference.
#' @param related_metadata Optional data frame describing related genomes. When
#'   supplied, it should contain a GenBank path column plus `completeness` and
#'   `ani_to_focal` columns so comparative genomes can be filtered before
#'   hotspot inference.
#' @param min_related_completeness Minimum completeness percentage required for
#'   comparative genomes when `related_metadata` is supplied.
#' @param min_related_ani Minimum ANI to the focal genome required for
#'   comparative genomes when `related_metadata` is supplied.
#' @param related_detection_mode Detection mode used for related genomes in the
#'   comparative layer. Defaults to `annotation` for speed.
#' @param detection_mode Detection mode. `hybrid` combines annotation and
#'   sequence-driven calls, `sequence` uses only the sequence engine, and
#'   `annotation` uses only annotation-derived calls.
#' @param prodigal_mode Prodigal mode used by the sequence engine, either
#'   `single` or `meta`.
#' @param isescan_dir Optional path to an ISEScan installation directory. If
#'   `NULL`, the directory is inferred from `isescan.py` in `PATH`.
#' @param isfinder_db Optional BLAST database prefix built from curated ISfinder
#'   references or a compatible custom IS reference set.
#' @param isfinder_fasta Optional nucleotide FASTA used to create a BLAST
#'   database for reference-based family refinement.
#' @param flank_bp Number of bases on each side of an element used when calling
#'   nearby mutation-prone genes.
#' @param pair_max_dist_bp Maximum distance between two candidate elements on
#'   the same contig when generating pair-mediated structural-variant scenarios.
#' @param export_excel Logical; if `TRUE`, also write an Excel workbook.
#' @param save_rds Logical; if `TRUE`, save the full result object as an RDS.
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @param max_target_sites_per_family Maximum number of predicted target sites
#'   retained per IS family after recognition-sequence and viability scoring.
#' @param max_variant_targets_per_element Maximum number of target-site variants
#'   retained per source element when building the variant catalog.
#' @param site_scan_strategy Target-site scan strategy. `balanced` keeps
#'   genome-native, higher-confidence models for speed; `full` scans all
#'   available models.
#'
#' @return An invisible list containing parsed features, genes, IS/mobile
#'   elements, mutation points, structural-variant scenarios, phenotype
#'   predictions, and output file paths.
#' @export
`%>%` <- dplyr::`%>%`

run_DNMB_mobileome <- function(
  genbank,
  output_dir = NULL,
  related_genbanks = NULL,
  related_metadata = NULL,
  min_related_completeness = 95,
  min_related_ani = 95,
  related_detection_mode = c("annotation", "hybrid", "sequence"),
  detection_mode = c("hybrid", "sequence", "annotation"),
  prodigal_mode = c("single", "meta"),
  isescan_dir = NULL,
  isfinder_db = NULL,
  isfinder_fasta = NULL,
  flank_bp = 10000,
  pair_max_dist_bp = 200000,
  export_excel = TRUE,
  save_rds = TRUE,
  verbose = TRUE,
  max_target_sites_per_family = 2000,
  max_variant_targets_per_element = 50,
  site_scan_strategy = c("balanced", "full")
) {
  if (missing(genbank) || is.null(genbank) || !nzchar(genbank)) {
    stop("`genbank` must be a path to a single GenBank file.")
  }
  if (!file.exists(genbank)) {
    stop("GenBank file does not exist: ", genbank)
  }
  detection_mode <- match.arg(detection_mode)
  related_detection_mode <- match.arg(related_detection_mode)
  prodigal_mode <- match.arg(prodigal_mode)
  site_scan_strategy <- match.arg(site_scan_strategy)
  min_related_completeness <- as.numeric(min_related_completeness)
  min_related_ani <- as.numeric(min_related_ani)

  flank_bp <- as.integer(flank_bp)
  pair_max_dist_bp <- as.integer(pair_max_dist_bp)
  if (is.na(flank_bp) || flank_bp < 0) {
    stop("`flank_bp` must be a non-negative integer.")
  }
  if (is.na(pair_max_dist_bp) || pair_max_dist_bp < 1) {
    stop("`pair_max_dist_bp` must be a positive integer.")
  }

  genbank <- normalizePath(genbank, mustWork = TRUE)
  if (is.null(output_dir)) {
    base_name <- tools::file_path_sans_ext(basename(genbank))
    output_dir <- file.path(dirname(genbank), paste0(base_name, "_dnmb_mobileome"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  .dnmb_mobileome_message(verbose, "Parsing GenBank features")
  parsed <- .dnmb_parse_genbank_features(genbank)
  genes <- .dnmb_build_gene_table(parsed$features)
  .dnmb_mobileome_message(verbose, "Scoring essentiality and redundancy")
  genes <- .dnmb_predict_gene_essentiality(genes)

  annotation_elements <- .dnmb_empty_elements()
  if (detection_mode %in% c("hybrid", "annotation")) {
    .dnmb_mobileome_message(verbose, "Detecting annotation-derived candidate IS/mobile elements")
    annotation_elements <- .dnmb_detect_is_elements(parsed$features, genes)
  }

  sequence_engine <- .dnmb_empty_sequence_engine()
  if (detection_mode %in% c("hybrid", "sequence")) {
    .dnmb_mobileome_message(verbose, "Running sequence-based IS engine")
    sequence_engine <- .dnmb_run_sequence_engine(
      parsed = parsed,
      output_dir = output_dir,
      prodigal_mode = prodigal_mode,
      isescan_dir = isescan_dir,
      isfinder_db = isfinder_db,
      isfinder_fasta = isfinder_fasta,
      verbose = verbose
    )
  } else {
    sequence_engine$status <- tibble::tibble(
      component = "sequence_engine",
      status = "skipped",
      detail = "detection_mode=annotation"
    )
  }

  .dnmb_mobileome_message(verbose, "Merging IS evidence layers")
  elements <- .dnmb_merge_element_calls(
    annotation_elements = annotation_elements,
    sequence_elements = sequence_engine$sequence_elements,
    genes = genes
  )

  .dnmb_mobileome_message(verbose, "Predicting local mutation points")
  mutation_points <- .dnmb_predict_mutation_points(elements, genes, flank_bp = flank_bp)

  .dnmb_mobileome_message(verbose, "Generating structural-variant scenarios")
  scenarios <- .dnmb_predict_variant_scenarios(
    elements,
    genes,
    pair_max_dist_bp = pair_max_dist_bp
  )

  .dnmb_mobileome_message(verbose, "Building target-site recognition models")
  target_models <- .dnmb_build_target_models(
    elements = elements,
    metadata = parsed$metadata,
    genes = genes
  )

  .dnmb_mobileome_message(verbose, "Scanning targetable genome sites")
  target_sites <- .dnmb_predict_target_sites(
    target_models = target_models,
    metadata = parsed$metadata,
    genes = genes,
    max_target_sites_per_family = max_target_sites_per_family,
    site_scan_strategy = site_scan_strategy
  )

  .dnmb_mobileome_message(verbose, "Aggregating family-wise targetable regions")
  targetable_regions <- .dnmb_build_targetable_regions(
    target_models = target_models,
    target_sites = target_sites,
    metadata = parsed$metadata,
    genes = genes,
    site_scan_strategy = site_scan_strategy
  )

  .dnmb_mobileome_message(verbose, "Enumerating genome-variant possibilities")
  variant_catalog <- .dnmb_build_variant_catalog(
    elements = elements,
    target_sites = target_sites,
    scenarios = scenarios,
    genes = genes,
    max_variant_targets_per_element = max_variant_targets_per_element
  )

  .dnmb_mobileome_message(verbose, "Building master integration table")
  master_table <- .dnmb_build_master_table(
    elements = elements,
    sequence_elements = sequence_engine$sequence_elements,
    target_models = target_models,
    targetable_regions = targetable_regions,
    target_sites = target_sites,
    variant_catalog = variant_catalog
  )

  comparative <- .dnmb_empty_comparative()
  if (!is.null(related_genbanks) && length(related_genbanks)) {
    .dnmb_mobileome_message(verbose, "Running comparative mobileome layer")
    focal_genome <- list(
      genome_id = "focal",
      genbank = genbank,
      metadata = parsed$metadata,
      features = parsed$features,
      genes = genes,
      elements = elements,
      proteins = genes %>%
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

    comparative <- .dnmb_run_comparative_mobileome(
      focal = focal_genome,
      related_genbanks = related_genbanks,
      related_metadata = related_metadata,
      min_related_completeness = min_related_completeness,
      min_related_ani = min_related_ani,
      output_dir = output_dir,
      detection_mode = related_detection_mode,
      verbose = verbose
    )
  }

  compact_table <- .dnmb_build_compact_summary_table(
    elements = elements,
    sequence_elements = sequence_engine$sequence_elements,
    target_models = target_models,
    targetable_regions = targetable_regions,
    variant_catalog = variant_catalog,
    comparative = comparative
  )

  .dnmb_mobileome_message(verbose, "Summarizing phenotype consequences")
  phenotype_predictions <- .dnmb_build_phenotype_predictions(
    elements = elements,
    mutation_points = mutation_points,
    scenarios = scenarios,
    variant_catalog = variant_catalog
  )

  summary_table <- tibble::tibble(
    metric = c(
      "source_genbank",
      "analysis_timestamp",
      "detection_mode",
      "related_genome_count",
      "site_scan_strategy",
      "contig_count",
      "feature_count",
      "gene_count",
      "annotation_element_count",
      "sequence_element_count",
      "merged_element_count",
      "reference_hit_count",
      "high_essential_gene_count",
      "target_model_count",
      "target_site_count",
      "targetable_region_count",
      "variant_catalog_count",
      "master_table_count",
      "compact_table_count",
      "comparative_locus_count",
      "comparative_hotspot_count",
      "comparative_chronology_count",
      "mutation_point_count",
      "structural_variant_scenario_count",
      "phenotype_prediction_count"
    ),
    value = c(
      genbank,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      detection_mode,
      as.character(if (is.null(related_genbanks)) 0 else length(related_genbanks)),
      site_scan_strategy,
      as.character(length(unique(parsed$features$contig))),
      as.character(nrow(parsed$features)),
      as.character(nrow(genes)),
      as.character(nrow(annotation_elements)),
      as.character(nrow(sequence_engine$sequence_elements)),
      as.character(nrow(elements)),
      as.character(nrow(sequence_engine$reference_hits)),
      as.character(sum(genes$essentiality_class == "high", na.rm = TRUE)),
      as.character(nrow(target_models)),
      as.character(nrow(target_sites)),
      as.character(nrow(targetable_regions)),
      as.character(nrow(variant_catalog)),
      as.character(nrow(master_table)),
      as.character(nrow(compact_table)),
      as.character(nrow(comparative$loci)),
      as.character(nrow(comparative$hotspots)),
      as.character(nrow(comparative$chronology)),
      as.character(nrow(mutation_points)),
      as.character(nrow(scenarios)),
      as.character(nrow(phenotype_predictions))
    )
  )

  summary_lines <- c(
    paste("DNMB mobileome analysis for:", basename(genbank)),
    paste("Source GenBank:", genbank),
    paste("Output directory:", output_dir),
    paste("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    paste("Contigs:", length(unique(parsed$features$contig))),
    paste("Parsed features:", nrow(parsed$features)),
    paste("Gene-like records:", nrow(genes)),
    paste("Detection mode:", detection_mode),
    paste("Related genomes:", if (is.null(related_genbanks)) 0 else length(related_genbanks)),
    paste("Target-site scan strategy:", site_scan_strategy),
    paste("Annotation-derived elements:", nrow(annotation_elements)),
    paste("Sequence-derived elements:", nrow(sequence_engine$sequence_elements)),
    paste("Merged IS/mobile elements:", nrow(elements)),
    paste("Reference-supported calls:", nrow(sequence_engine$reference_hits)),
    paste("High-essentiality genes:", sum(genes$essentiality_class == "high", na.rm = TRUE)),
    paste("Target-site models:", nrow(target_models)),
    paste("Predicted target sites:", nrow(target_sites)),
    paste("Targetable regions:", nrow(targetable_regions)),
    paste("Variant catalog rows:", nrow(variant_catalog)),
    paste("Master table rows:", nrow(master_table)),
    paste("Compact table rows:", nrow(compact_table)),
    paste("Comparative loci:", nrow(comparative$loci)),
    paste("Comparative hotspots:", nrow(comparative$hotspots)),
    paste("Chronology rows:", nrow(comparative$chronology)),
    paste("Predicted mutation-point rows:", nrow(mutation_points)),
    paste("Predicted structural-variant scenarios:", nrow(scenarios)),
    paste("Phenotype prediction rows:", nrow(phenotype_predictions)),
    "",
    "Interpretation notes:",
    "- Hybrid mode combines annotation evidence with ISEScan-inspired sequence evidence.",
    "- A single genome is sufficient for local impact hypotheses, but not for destination-site prediction after transposition.",
    "- ISfinder-style reference refinement is used only when a curated reference FASTA or BLAST DB is supplied.",
    "- Phenotype outputs indicate likely affected functions, not guaranteed direction or magnitude of change."
  )

  file_map <- list(
    summary_tsv = file.path(output_dir, "mobileome_summary.tsv"),
    summary_txt = file.path(output_dir, "mobileome_summary.txt"),
    features_tsv = file.path(output_dir, "genbank_features.tsv"),
    genes_tsv = file.path(output_dir, "gene_table.tsv"),
    essentiality_tsv = file.path(output_dir, "gene_essentiality.tsv"),
    engine_status_tsv = file.path(output_dir, "sequence_engine_status.tsv"),
    orfs_tsv = file.path(output_dir, "sequence_orfs.tsv"),
    hmmer_hits_tsv = file.path(output_dir, "isescan_hmmer_hits.tsv"),
    annotation_elements_tsv = file.path(output_dir, "annotation_is_elements.tsv"),
    sequence_elements_tsv = file.path(output_dir, "sequence_is_elements.tsv"),
    reference_hits_tsv = file.path(output_dir, "isfinder_reference_hits.tsv"),
    target_models_tsv = file.path(output_dir, "is_target_models.tsv"),
    target_sites_tsv = file.path(output_dir, "is_target_sites.tsv"),
    targetable_regions_tsv = file.path(output_dir, "is_targetable_regions.tsv"),
    targetable_regions_bed = file.path(output_dir, "is_targetable_regions.bed"),
    variant_catalog_tsv = file.path(output_dir, "genome_variant_catalog.tsv"),
    master_table_tsv = file.path(output_dir, "mobileome_master_table.tsv"),
    compact_table_tsv = file.path(output_dir, "mobileome_compact_table.tsv"),
    comparative_status_tsv = file.path(output_dir, "comparative_status.tsv"),
    comparative_loci_tsv = file.path(output_dir, "comparative", "comparative_loci.tsv"),
    occupied_empty_tsv = file.path(output_dir, "comparative", "occupied_empty_matrix.tsv"),
    comparative_hotspots_tsv = file.path(output_dir, "comparative", "comparative_hotspots.tsv"),
    chronology_tsv = file.path(output_dir, "comparative", "event_chronology.tsv"),
    family_locus_master_tsv = file.path(output_dir, "comparative", "family_locus_master.tsv"),
    elements_tsv = file.path(output_dir, "is_elements.tsv"),
    mutation_points_tsv = file.path(output_dir, "is_mutation_points.tsv"),
    scenarios_tsv = file.path(output_dir, "is_variant_scenarios.tsv"),
    phenotypes_tsv = file.path(output_dir, "phenotype_predictions.tsv")
  )

  .dnmb_write_tsv(summary_table, file_map$summary_tsv)
  writeLines(summary_lines, con = file_map$summary_txt)
  .dnmb_write_tsv(parsed$features, file_map$features_tsv)
  .dnmb_write_tsv(genes, file_map$genes_tsv)
  .dnmb_write_tsv(
    genes %>%
      dplyr::select(
        gene_id,
        contig,
        locus_tag,
        gene,
        product,
        essentiality_score,
        essentiality_class,
        redundancy_count,
        essentiality_reasons
      ),
    file_map$essentiality_tsv
  )
  .dnmb_write_tsv(sequence_engine$status, file_map$engine_status_tsv)
  .dnmb_write_tsv(sequence_engine$orfs, file_map$orfs_tsv)
  .dnmb_write_tsv(sequence_engine$hmmer_hits, file_map$hmmer_hits_tsv)
  .dnmb_write_tsv(annotation_elements, file_map$annotation_elements_tsv)
  .dnmb_write_tsv(sequence_engine$sequence_elements, file_map$sequence_elements_tsv)
  .dnmb_write_tsv(sequence_engine$reference_hits, file_map$reference_hits_tsv)
  .dnmb_write_tsv(target_models, file_map$target_models_tsv)
  .dnmb_write_tsv(target_sites, file_map$target_sites_tsv)
  .dnmb_write_tsv(targetable_regions, file_map$targetable_regions_tsv)
  .dnmb_write_bed(
    targetable_regions %>%
      dplyr::transmute(
        chrom = .data$contig,
        chromStart = pmax(0L, .data$region_start - 1L),
        chromEnd = .data$region_end,
        name = paste(.data$target_family, .data$region_type, sep = "|"),
        score = as.integer(round(1000 * pmin(1, .data$max_target_site_score))),
        strand = "."
      ),
    file_map$targetable_regions_bed
  )
  .dnmb_write_tsv(variant_catalog, file_map$variant_catalog_tsv)
  .dnmb_write_tsv(master_table, file_map$master_table_tsv)
  .dnmb_write_tsv(compact_table, file_map$compact_table_tsv)
  .dnmb_write_tsv(comparative$status, file_map$comparative_status_tsv)
  .dnmb_write_tsv(elements, file_map$elements_tsv)
  .dnmb_write_tsv(mutation_points, file_map$mutation_points_tsv)
  .dnmb_write_tsv(scenarios, file_map$scenarios_tsv)
  .dnmb_write_tsv(phenotype_predictions, file_map$phenotypes_tsv)

  if (isTRUE(export_excel)) {
    excel_path <- file.path(output_dir, "dnmb_mobileome_report.xlsx")
    .dnmb_write_mobileome_workbook(
      excel_path = excel_path,
      summary_table = summary_table,
      features = parsed$features,
      genes = genes,
      engine_status = sequence_engine$status,
      orfs = sequence_engine$orfs,
      hmmer_hits = sequence_engine$hmmer_hits,
      annotation_elements = annotation_elements,
      sequence_elements = sequence_engine$sequence_elements,
      reference_hits = sequence_engine$reference_hits,
      target_models = target_models,
      target_sites = target_sites,
      targetable_regions = targetable_regions,
      variant_catalog = variant_catalog,
      master_table = master_table,
      compact_table = compact_table,
      comparative = comparative,
      elements = elements,
      mutation_points = mutation_points,
      scenarios = scenarios,
      phenotype_predictions = phenotype_predictions
    )
    file_map$excel <- excel_path
  }

  result <- list(
    source_genbank = genbank,
    output_dir = output_dir,
    metadata = parsed$metadata,
    features = parsed$features,
    genes = genes,
    annotation_elements = annotation_elements,
    sequence_engine = sequence_engine,
    target_models = target_models,
    target_sites = target_sites,
    targetable_regions = targetable_regions,
    variant_catalog = variant_catalog,
    master_table = master_table,
    compact_table = compact_table,
    comparative = comparative,
    elements = elements,
    mutation_points = mutation_points,
    scenarios = scenarios,
    phenotype_predictions = phenotype_predictions,
    files = file_map
  )

  if (isTRUE(save_rds)) {
    rds_path <- file.path(output_dir, "dnmb_mobileome_result.rds")
    saveRDS(result, file = rds_path)
    result$files$rds <- rds_path
  }

  .dnmb_mobileome_message(verbose, "Mobileome pipeline finished")
  invisible(result)
}

.dnmb_mobileome_message <- function(verbose, text) {
  if (isTRUE(verbose)) {
    message("[DNMB mobileome] ", text)
  }
}

.dnmb_write_tsv <- function(x, path) {
  utils::write.table(
    x = x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = ""
  )
}

.dnmb_write_bed <- function(x, path) {
  utils::write.table(
    x = x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

.dnmb_write_mobileome_workbook <- function(
  excel_path,
  summary_table,
  features,
  genes,
  engine_status,
  orfs,
  hmmer_hits,
  annotation_elements,
  sequence_elements,
  reference_hits,
  target_models,
  target_sites,
  targetable_regions,
  variant_catalog,
  master_table,
  compact_table,
  comparative,
  elements,
  mutation_points,
  scenarios,
  phenotype_predictions
) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", summary_table)
  openxlsx::addWorksheet(wb, "Features")
  openxlsx::writeData(wb, "Features", features)
  openxlsx::addWorksheet(wb, "Genes")
  openxlsx::writeData(wb, "Genes", genes)
  openxlsx::addWorksheet(wb, "Engine_Status")
  openxlsx::writeData(wb, "Engine_Status", engine_status)
  openxlsx::addWorksheet(wb, "Sequence_ORFs")
  openxlsx::writeData(wb, "Sequence_ORFs", orfs)
  openxlsx::addWorksheet(wb, "HMMER_Hits")
  openxlsx::writeData(wb, "HMMER_Hits", hmmer_hits)
  openxlsx::addWorksheet(wb, "Annotation_IS")
  openxlsx::writeData(wb, "Annotation_IS", annotation_elements)
  openxlsx::addWorksheet(wb, "Sequence_IS")
  openxlsx::writeData(wb, "Sequence_IS", sequence_elements)
  openxlsx::addWorksheet(wb, "Reference_Hits")
  openxlsx::writeData(wb, "Reference_Hits", reference_hits)
  openxlsx::addWorksheet(wb, "Target_Models")
  openxlsx::writeData(wb, "Target_Models", target_models)
  openxlsx::addWorksheet(wb, "Target_Sites")
  openxlsx::writeData(wb, "Target_Sites", target_sites)
  openxlsx::addWorksheet(wb, "Targetable_Regions")
  openxlsx::writeData(wb, "Targetable_Regions", targetable_regions)
  openxlsx::addWorksheet(wb, "Variant_Catalog")
  openxlsx::writeData(wb, "Variant_Catalog", variant_catalog)
  openxlsx::addWorksheet(wb, "Master_Table")
  openxlsx::writeData(wb, "Master_Table", master_table)
  openxlsx::addWorksheet(wb, "Compact_Table")
  openxlsx::writeData(wb, "Compact_Table", compact_table)
  openxlsx::addWorksheet(wb, "Comparative_Status")
  openxlsx::writeData(wb, "Comparative_Status", comparative$status)
  if (nrow(comparative$loci)) {
    openxlsx::addWorksheet(wb, "Comparative_Loci")
    openxlsx::writeData(wb, "Comparative_Loci", comparative$loci)
  }
  if (nrow(comparative$occupied_empty)) {
    openxlsx::addWorksheet(wb, "Occupied_Empty")
    openxlsx::writeData(wb, "Occupied_Empty", comparative$occupied_empty)
  }
  if (nrow(comparative$hotspots)) {
    openxlsx::addWorksheet(wb, "Comparative_Hotspots")
    openxlsx::writeData(wb, "Comparative_Hotspots", comparative$hotspots)
  }
  if (nrow(comparative$chronology)) {
    openxlsx::addWorksheet(wb, "Chronology")
    openxlsx::writeData(wb, "Chronology", comparative$chronology)
  }
  if (nrow(comparative$locus_master)) {
    openxlsx::addWorksheet(wb, "Family_Locus_Master")
    openxlsx::writeData(wb, "Family_Locus_Master", comparative$locus_master)
  }
  openxlsx::addWorksheet(wb, "IS_Elements")
  openxlsx::writeData(wb, "IS_Elements", elements)
  openxlsx::addWorksheet(wb, "Mutation_Points")
  openxlsx::writeData(wb, "Mutation_Points", mutation_points)
  openxlsx::addWorksheet(wb, "SV_Scenarios")
  openxlsx::writeData(wb, "SV_Scenarios", scenarios)
  openxlsx::addWorksheet(wb, "Phenotypes")
  openxlsx::writeData(wb, "Phenotypes", phenotype_predictions)
  openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
}

.dnmb_parse_genbank_features <- function(genbank) {
  lines <- readLines(genbank, warn = FALSE)
  locus_idx <- grep("^LOCUS\\s+", lines)
  if (!length(locus_idx)) {
    stop("No LOCUS records were found in the GenBank file.")
  }

  end_idx <- grep("^//\\s*$", lines)
  features_list <- vector("list", length(locus_idx))
  metadata_list <- vector("list", length(locus_idx))

  for (i in seq_along(locus_idx)) {
    start_idx <- locus_idx[i]
    candidate_ends <- end_idx[end_idx >= start_idx]
    end_record <- if (length(candidate_ends)) candidate_ends[1] else length(lines)
    block <- lines[start_idx:end_record]

    locus_line <- block[1]
    contig <- .dnmb_extract_locus_name(locus_line)
    definition <- .dnmb_extract_multiline_field(block, "DEFINITION")
    accession <- .dnmb_extract_multiline_field(block, "ACCESSION")
    version <- .dnmb_extract_multiline_field(block, "VERSION")

    features_start <- grep("^FEATURES\\s+", block)
    origin_start <- grep("^ORIGIN\\s*$", block)
    sequence <- .dnmb_extract_origin_sequence(block)
    metadata_list[[i]] <- tibble::tibble(
      contig = contig,
      contig_number = i,
      definition = definition,
      accession = accession,
      version = version,
      sequence_length_bp = nchar(sequence),
      sequence = sequence
    )

    if (!length(features_start) || !length(origin_start) || origin_start[1] <= features_start[1]) {
      features_list[[i]] <- tibble::tibble()
      next
    }

    feature_lines <- block[(features_start[1] + 1):(origin_start[1] - 1)]
    features_list[[i]] <- .dnmb_parse_feature_block(
      feature_lines = feature_lines,
      contig = contig,
      contig_number = i,
      definition = definition,
      source_file = basename(genbank)
    )
  }

  list(
    metadata = dplyr::bind_rows(metadata_list),
    features = dplyr::bind_rows(features_list)
  )
}

.dnmb_extract_locus_name <- function(locus_line) {
  tokens <- strsplit(trimws(locus_line), "\\s+")[[1]]
  if (length(tokens) >= 2) {
    tokens[2]
  } else {
    "unknown_contig"
  }
}

.dnmb_extract_multiline_field <- function(lines, field_name) {
  idx <- grep(paste0("^", field_name, "\\s+"), lines)
  if (!length(idx)) {
    return(NA_character_)
  }

  value <- sub(paste0("^", field_name, "\\s+"), "", lines[idx[1]])
  next_idx <- idx[1] + 1
  while (next_idx <= length(lines) && grepl("^\\s{12}\\S", lines[next_idx])) {
    value <- paste(value, trimws(lines[next_idx]))
    next_idx <- next_idx + 1
  }
  stringr::str_squish(value)
}

.dnmb_extract_origin_sequence <- function(block) {
  origin_start <- grep("^ORIGIN\\s*$", block)
  if (!length(origin_start)) {
    return("")
  }
  end_idx <- grep("^//\\s*$", block)
  end_idx <- if (length(end_idx)) end_idx[1] else length(block)
  seq_lines <- block[(origin_start[1] + 1):(end_idx - 1)]
  seq_text <- paste(seq_lines, collapse = "")
  seq_text <- gsub("[0-9[:space:]]+", "", seq_text)
  toupper(seq_text)
}

.dnmb_parse_feature_block <- function(feature_lines, contig, contig_number, definition, source_file) {
  current <- NULL
  current_qualifier <- NULL
  parsed_features <- list()

  push_current <- function(feature_obj) {
    if (is.null(feature_obj)) {
      return(NULL)
    }
    parsed_features[[length(parsed_features) + 1]] <<- .dnmb_feature_to_row(
      feature_obj = feature_obj,
      contig = contig,
      contig_number = contig_number,
      definition = definition,
      source_file = source_file
    )
    NULL
  }

  for (line in feature_lines) {
    if (grepl("^\\s{5}\\S", line)) {
      push_current(current)
      current <- list(
        feature_type = trimws(substr(line, 6, 20)),
        raw_location = trimws(substr(line, 22, nchar(line))),
        qualifiers = list()
      )
      current_qualifier <- NULL
      next
    }

    if (is.null(current)) {
      next
    }

    if (grepl("^\\s{21}/", line)) {
      qualifier_line <- trimws(line)
      qualifier_line <- sub("^/", "", qualifier_line)
      if (grepl("=", qualifier_line, fixed = TRUE)) {
        qualifier_key <- sub("=.*$", "", qualifier_line)
        qualifier_value <- sub("^[^=]+=", "", qualifier_line)
        qualifier_value <- .dnmb_strip_quotes(qualifier_value)
      } else {
        qualifier_key <- qualifier_line
        qualifier_value <- "TRUE"
      }

      if (is.null(current$qualifiers[[qualifier_key]])) {
        current$qualifiers[[qualifier_key]] <- qualifier_value
      } else {
        current$qualifiers[[qualifier_key]] <- c(
          current$qualifiers[[qualifier_key]],
          qualifier_value
        )
      }
      current_qualifier <- qualifier_key
      next
    }

    continuation <- trimws(line)
    if (!nzchar(continuation)) {
      next
    }

    if (!is.null(current_qualifier)) {
      values <- current$qualifiers[[current_qualifier]]
      values[length(values)] <- stringr::str_squish(paste(
        values[length(values)],
        .dnmb_strip_quotes(continuation)
      ))
      current$qualifiers[[current_qualifier]] <- values
    } else {
      current$raw_location <- paste0(current$raw_location, continuation)
    }
  }

  push_current(current)
  dplyr::bind_rows(parsed_features)
}

.dnmb_strip_quotes <- function(x) {
  x <- sub('^"', "", x)
  x <- sub('"$', "", x)
  x
}

.dnmb_feature_to_row <- function(feature_obj, contig, contig_number, definition, source_file) {
  coords <- .dnmb_parse_location(feature_obj$raw_location)
  qualifier_text <- .dnmb_collapse_qualifiers(feature_obj$qualifiers)
  tibble::tibble(
    contig = contig,
    contig_number = contig_number,
    definition = definition,
    source_file = source_file,
    feature_type = feature_obj$feature_type,
    raw_location = feature_obj$raw_location,
    start = coords$start,
    end = coords$end,
    strand = coords$strand,
    segment_count = coords$segment_count,
    partial = coords$partial,
    joined = coords$joined,
    locus_tag = .dnmb_first_qualifier_value(feature_obj$qualifiers, "locus_tag"),
    gene = .dnmb_first_qualifier_value(feature_obj$qualifiers, "gene"),
    product = .dnmb_first_qualifier_value(feature_obj$qualifiers, "product"),
    note = .dnmb_collapse_qualifier_value(feature_obj$qualifiers, "note"),
    mobile_element_type = .dnmb_collapse_qualifier_value(feature_obj$qualifiers, "mobile_element_type"),
    protein_id = .dnmb_first_qualifier_value(feature_obj$qualifiers, "protein_id"),
    translation = .dnmb_first_qualifier_value(feature_obj$qualifiers, "translation"),
    qualifiers_text = qualifier_text
  )
}

.dnmb_first_qualifier_value <- function(qualifiers, key) {
  values <- qualifiers[[key]]
  if (is.null(values) || !length(values)) {
    return(NA_character_)
  }
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) {
    return(NA_character_)
  }
  values[1]
}

.dnmb_collapse_qualifier_value <- function(qualifiers, key) {
  values <- qualifiers[[key]]
  if (is.null(values) || !length(values)) {
    return(NA_character_)
  }
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) {
    return(NA_character_)
  }
  paste(unique(values), collapse = "; ")
}

.dnmb_collapse_qualifiers <- function(qualifiers) {
  if (!length(qualifiers)) {
    return(NA_character_)
  }
  parts <- vapply(
    names(qualifiers),
    function(key) {
      value <- qualifiers[[key]]
      value <- value[!is.na(value) & nzchar(value)]
      if (!length(value)) {
        return(NA_character_)
      }
      paste0(key, "=", paste(unique(value), collapse = " | "))
    },
    FUN.VALUE = character(1)
  )
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (!length(parts)) {
    return(NA_character_)
  }
  paste(parts, collapse = "; ")
}

.dnmb_parse_location <- function(raw_location) {
  positions <- stringr::str_extract_all(raw_location, "\\d+")[[1]]
  positions <- as.integer(positions)
  positions <- positions[!is.na(positions)]

  if (!length(positions)) {
    return(list(
      start = NA_integer_,
      end = NA_integer_,
      strand = NA_character_,
      segment_count = 0L,
      partial = grepl("[<>]", raw_location),
      joined = grepl("join|order", raw_location)
    ))
  }

  list(
    start = min(positions),
    end = max(positions),
    strand = if (grepl("complement", raw_location)) "-" else "+",
    segment_count = length(positions) %/% 2L,
    partial = grepl("[<>]", raw_location),
    joined = grepl("join|order", raw_location)
  )
}

.dnmb_build_gene_table <- function(features) {
  gene_like_types <- c("CDS", "gene", "rRNA", "tRNA", "ncRNA", "misc_RNA", "tmRNA")
  genes <- features %>%
    dplyr::filter(.data$feature_type %in% gene_like_types) %>%
    dplyr::mutate(
      locus_tag = dplyr::na_if(.data$locus_tag, ""),
      gene_key = dplyr::if_else(
        !is.na(.data$locus_tag),
        .data$locus_tag,
        paste(.data$contig, .data$start, .data$end, .data$strand, sep = ":")
      ),
      feature_rank = dplyr::case_when(
        .data$feature_type == "CDS" ~ 1L,
        .data$feature_type == "gene" ~ 2L,
        .data$feature_type == "rRNA" ~ 3L,
        .data$feature_type == "tRNA" ~ 4L,
        TRUE ~ 5L
      )
    ) %>%
    dplyr::arrange(.data$contig, .data$start, .data$end, .data$feature_rank) %>%
    dplyr::group_by(.data$gene_key) %>%
    dplyr::summarise(
      contig = .dnmb_first_non_empty(.data$contig),
      contig_number = dplyr::first(.data$contig_number),
      definition = .dnmb_first_non_empty(.data$definition),
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      strand = .dnmb_first_non_empty(.data$strand),
      locus_tag = .dnmb_first_non_empty(.data$locus_tag),
      gene = .dnmb_first_non_empty(.data$gene),
      product = .dnmb_first_non_empty(.data$product),
      note = .dnmb_first_non_empty(.data$note),
      protein_id = .dnmb_first_non_empty(.data$protein_id),
      translation = .dnmb_first_non_empty(.data$translation),
      feature_types = paste(unique(.data$feature_type), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene_id = dplyr::if_else(
        !is.na(.data$locus_tag),
        .data$locus_tag,
        paste0("GENE_", sprintf("%05d", dplyr::row_number()))
      ),
      gene_length_bp = .data$end - .data$start + 1L,
      annotation_text = stringr::str_squish(paste(
        .data$locus_tag,
        .data$gene,
        .data$product,
        .data$note
      ))
    ) %>%
    dplyr::select(
      gene_id,
      contig,
      contig_number,
      definition,
      locus_tag,
      gene,
      product,
      note,
      protein_id,
      translation,
      feature_types,
      start,
      end,
      gene_length_bp,
      strand,
      annotation_text
    )

  if (!nrow(genes)) {
    genes <- tibble::tibble(
      gene_id = character(),
      contig = character(),
      contig_number = integer(),
      definition = character(),
      locus_tag = character(),
      gene = character(),
      product = character(),
      note = character(),
      protein_id = character(),
      translation = character(),
      feature_types = character(),
      start = integer(),
      end = integer(),
      gene_length_bp = integer(),
      strand = character(),
      annotation_text = character()
    )
  }

  genes
}

.dnmb_first_non_empty <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  x[1]
}

.dnmb_detect_is_elements <- function(features, genes, merge_gap_bp = 2000L) {
  if (!nrow(features)) {
    return(.dnmb_empty_elements())
  }

  element_patterns <- paste(
    c(
      "\\binsertion sequence\\b",
      "\\btransposase\\b",
      "\\btransposon\\b",
      "\\bmobile element\\b",
      "\\bIS\\d+[A-Za-z0-9_-]*\\b",
      "\\bIS[A-Z][A-Za-z0-9_-]*\\b",
      "\\bTn\\d+[A-Za-z0-9_-]*\\b",
      "\\btnp[[:alnum:]_-]*\\b",
      "family transposase"
    ),
    collapse = "|"
  )

  seed_table <- features %>%
    dplyr::mutate(
      evidence_text = stringr::str_to_lower(stringr::str_squish(paste(
        .data$feature_type,
        .data$gene,
        .data$product,
        .data$note,
        .data$mobile_element_type,
        .data$qualifiers_text
      ))),
      direct_mobile_feature = .data$feature_type %in% c("mobile_element", "repeat_region"),
      keyword_match = stringr::str_detect(.data$evidence_text, stringr::regex(element_patterns, ignore_case = TRUE)),
      is_candidate = (.data$feature_type == "mobile_element") |
        (.data$direct_mobile_feature & .data$keyword_match) |
        (.data$feature_type %in% c("CDS", "gene") & .data$keyword_match),
      family_text = paste(.data$mobile_element_type, .data$product, .data$note, .data$qualifiers_text)
    ) %>%
    dplyr::filter(.data$is_candidate) %>%
    dplyr::mutate(
      element_family = vapply(.data$family_text, .dnmb_extract_is_family, character(1)),
      confidence_score = dplyr::case_when(
        .data$feature_type == "mobile_element" ~ 3L,
        !is.na(.data$mobile_element_type) & nzchar(.data$mobile_element_type) ~ 3L,
        stringr::str_detect(.data$evidence_text, "insertion sequence|transposase|family transposase") ~ 2L,
        TRUE ~ 1L
      ),
      evidence_label = dplyr::case_when(
        .data$feature_type == "mobile_element" ~ "annotated_mobile_element",
        .data$feature_type == "repeat_region" ~ "annotated_repeat_region",
        TRUE ~ "mobile_gene_annotation"
      )
    ) %>%
    dplyr::filter(!is.na(.data$start), !is.na(.data$end)) %>%
    dplyr::arrange(.data$contig, .data$start, .data$end)

  # Filter out false positives from mobile_gene_annotation category
  # These are genes where IS-like patterns matched non-IS annotations
  # (e.g., "isomerase" matching IS\d+, "transporter" matching tnp)
  if (nrow(seed_table)) {
    fp_patterns <- paste(c(
      "\\bisomerase\\b", "\\bkinase\\b", "\\bsynthase\\b", "\\bsynthetase\\b",
      "\\breductase\\b", "\\bdehydrogenase\\b", "\\boxidase\\b", "\\bmutase\\b",
      "\\bligase\\b", "\\blyase\\b", "\\bhydrolase\\b", "\\bepimerase\\b",
      "\\bphosphatase\\b", "\\btransporter\\b", "\\bpermease\\b",
      "\\bribosomal protein\\b", "\\bATP synthase\\b", "\\bflagellar\\b",
      "\\bchemotaxis\\b", "\\bammonium\\b", "\\bcytochrome\\b",
      "\\btRNA\\b", "\\bphage tail\\b"
    ), collapse = "|")

    is_mobile_gene_annot <- seed_table$evidence_label == "mobile_gene_annotation"
    has_real_is_keyword <- stringr::str_detect(
      seed_table$evidence_text,
      stringr::regex("\\btransposase\\b|\\bmobile.element\\b|\\binsertion.sequence\\b|\\bTnpB\\b|\\bTnpA\\b|\\bIS\\d{2,}\\b",
                     ignore_case = TRUE)
    )
    has_fp_keyword <- stringr::str_detect(seed_table$product, stringr::regex(fp_patterns, ignore_case = TRUE))

    # Remove: mobile_gene_annotation + has false positive keyword + no real IS keyword
    is_false_positive <- is_mobile_gene_annot & has_fp_keyword & !has_real_is_keyword
    if (any(is_false_positive)) {
      seed_table <- seed_table[!is_false_positive, , drop = FALSE]
    }
  }

  if (!nrow(seed_table)) {
    return(.dnmb_empty_elements())
  }

  element_rows <- list()
  current <- seed_table[1, , drop = FALSE]

  for (idx in 2:nrow(seed_table)) {
    next_row <- seed_table[idx, , drop = FALSE]
    same_contig <- identical(current$contig[[1]], next_row$contig[[1]])
    compatible_family <- .dnmb_family_compatible(current$element_family[[1]], next_row$element_family[[1]])
    # isTRUE() guards against NA in any of next_row$start, current$end, or
    # merge_gap_bp — if coordinates are missing, treat the rows as not
    # mergeable rather than letting `if (... && NA)` abort the module.
    close_enough <- isTRUE(next_row$start[[1]] <= (current$end[[1]] + merge_gap_bp))

    if (same_contig && compatible_family && close_enough) {
      current$start[[1]] <- min(current$start[[1]], next_row$start[[1]], na.rm = TRUE)
      current$end[[1]] <- max(current$end[[1]], next_row$end[[1]], na.rm = TRUE)
      current$strand[[1]] <- .dnmb_merge_strand(current$strand[[1]], next_row$strand[[1]])
      current$confidence_score[[1]] <- max(current$confidence_score[[1]], next_row$confidence_score[[1]], na.rm = TRUE)
      current$element_family[[1]] <- .dnmb_merge_family(current$element_family[[1]], next_row$element_family[[1]])
      current$evidence_label[[1]] <- paste(unique(c(
        unlist(strsplit(current$evidence_label[[1]], "; ")),
        next_row$evidence_label[[1]]
      )), collapse = "; ")
      current$locus_tag[[1]] <- paste(unique(na.omit(c(
        unlist(strsplit(dplyr::coalesce(current$locus_tag[[1]], ""), "; ")),
        next_row$locus_tag[[1]]
      ))), collapse = "; ")
      current$product[[1]] <- paste(unique(na.omit(c(
        unlist(strsplit(dplyr::coalesce(current$product[[1]], ""), "; ")),
        next_row$product[[1]]
      ))), collapse = "; ")
      current$feature_type[[1]] <- paste(unique(c(
        unlist(strsplit(current$feature_type[[1]], "; ")),
        next_row$feature_type[[1]]
      )), collapse = "; ")
    } else {
      element_rows[[length(element_rows) + 1]] <- current
      current <- next_row
    }
  }
  element_rows[[length(element_rows) + 1]] <- current

  elements <- dplyr::bind_rows(element_rows) %>%
    dplyr::mutate(
      element_id = paste0("IS_", sprintf("%04d", dplyr::row_number())),
      size_bp = .data$end - .data$start + 1L,
      confidence = dplyr::case_when(
        .data$confidence_score >= 3L ~ "high",
        .data$confidence_score == 2L ~ "medium",
        TRUE ~ "low"
      ),
      evidence_summary = stringr::str_squish(paste(
        .data$evidence_label,
        .data$product,
        .data$mobile_element_type
      ))
    ) %>%
    dplyr::select(
      element_id,
      contig,
      contig_number,
      definition,
      start,
      end,
      size_bp,
      strand,
      element_family,
      mobile_element_type,
      confidence,
      confidence_score,
      feature_type,
      locus_tag,
      product,
      evidence_label,
      evidence_summary
    )

  if (!nrow(genes)) {
    elements$supporting_gene_count <- 0L
    elements$supporting_genes <- NA_character_
    elements$supporting_products <- NA_character_
    return(elements)
  }

  support_list <- lapply(seq_len(nrow(elements)), function(i) {
    overlaps <- genes %>%
      dplyr::filter(
        .data$contig == elements$contig[[i]],
        .data$start <= elements$end[[i]],
        .data$end >= elements$start[[i]]
      )
    if (!nrow(overlaps)) {
      return(tibble::tibble(
        supporting_gene_count = 0L,
        supporting_genes = NA_character_,
        supporting_products = NA_character_
      ))
    }
    tibble::tibble(
      supporting_gene_count = nrow(overlaps),
      supporting_genes = paste(unique(stats::na.omit(overlaps$gene_id)), collapse = "; "),
      supporting_products = paste(unique(stats::na.omit(overlaps$product)), collapse = "; ")
    )
  })

  dplyr::bind_cols(elements, dplyr::bind_rows(support_list))
}

.dnmb_empty_elements <- function() {
  tibble::tibble(
    element_id = character(),
    contig = character(),
    contig_number = integer(),
    definition = character(),
    start = integer(),
    end = integer(),
    size_bp = integer(),
    strand = character(),
    element_family = character(),
    mobile_element_type = character(),
    confidence = character(),
    confidence_score = integer(),
    feature_type = character(),
    locus_tag = character(),
    product = character(),
    evidence_label = character(),
    evidence_summary = character(),
    supporting_gene_count = integer(),
    supporting_genes = character(),
    supporting_products = character()
  )
}

.dnmb_family_compatible <- function(a, b) {
  if (is.na(a) || !nzchar(a) || is.na(b) || !nzchar(b)) {
    return(TRUE)
  }
  identical(a, b)
}

.dnmb_merge_family <- function(a, b) {
  if (!is.na(a) && nzchar(a)) {
    return(a)
  }
  if (!is.na(b) && nzchar(b)) {
    return(b)
  }
  NA_character_
}

.dnmb_merge_strand <- function(a, b) {
  if (is.na(a) || !nzchar(a)) {
    return(b)
  }
  if (is.na(b) || !nzchar(b)) {
    return(a)
  }
  if (identical(a, b)) {
    return(a)
  }
  "mixed"
}

.dnmb_known_is_families <- function() {
  c(
    "IS200/IS605", "ISAZO13", "IS1182", "IS1380", "IS1595", "IS1634",
    "ISAS1", "ISKRA4", "ISNCY", "ISH3", "IS110", "IS21", "IS256",
    "IS30", "IS481", "IS607", "IS630", "IS701", "IS982", "ISL3",
    "IS66", "IS91", "IS3", "IS4", "IS5", "IS6", "IS1"
  )
}

.dnmb_extract_is_family <- function(text) {
  if (all(is.na(text)) || !nzchar(paste(text, collapse = ""))) {
    return(NA_character_)
  }

  known <- .dnmb_known_is_families()
  patterns <- stringr::str_replace_all(known, "([/])", "\\\\\\1")
  patterns <- patterns[order(nchar(patterns), decreasing = TRUE)]
  candidate <- stringr::str_extract(text, stringr::regex(paste0("\\b(", paste(patterns, collapse = "|"), ")\\b"), ignore_case = TRUE))
  if (!is.na(candidate) && nzchar(candidate)) {
    return(toupper(candidate))
  }

  family_phrase <- stringr::str_extract(
    text,
    stringr::regex("\\bIS\\d+\\s+family\\b", ignore_case = TRUE)
  )
  if (!is.na(family_phrase) && nzchar(family_phrase)) {
    return(stringr::str_to_upper(stringr::str_replace(family_phrase, "\\s+family$", "")))
  }

  NA_character_
}

.dnmb_predict_mutation_points <- function(elements, genes, flank_bp = 10000L) {
  if (!nrow(elements) || !nrow(genes)) {
    return(.dnmb_empty_mutation_points())
  }

  edt <- data.table::as.data.table(elements)
  gdt <- data.table::as.data.table(genes)

  edt[, win_start := pmax(1L, start - as.integer(flank_bp))]
  edt[, win_end := end + as.integer(flank_bp)]
  edt[, element_start := start]
  edt[, element_end := end]

  data.table::setnames(edt, old = c("start", "end"), new = c("interval_start", "interval_end"))
  data.table::setnames(gdt, old = c("start", "end"), new = c("gene_start", "gene_end"))

  data.table::setkey(edt, contig, win_start, win_end)
  data.table::setkey(gdt, contig, gene_start, gene_end)

  joined <- data.table::foverlaps(
    gdt,
    edt,
    by.x = c("contig", "gene_start", "gene_end"),
    by.y = c("contig", "win_start", "win_end"),
    type = "any",
    nomatch = 0L
  )

  if (!nrow(joined)) {
    return(.dnmb_empty_mutation_points())
  }

  joined[, relationship := data.table::fcase(
    gene_start <= element_end & gene_end >= element_start, "direct_overlap",
    gene_end < element_start, "upstream_neighbor",
    gene_start > element_end, "downstream_neighbor",
    default = "local_neighbor"
  )]
  joined[, distance_to_element_bp := data.table::fcase(
    relationship == "direct_overlap", 0L,
    relationship == "upstream_neighbor", as.integer(element_start - gene_end),
    relationship == "downstream_neighbor", as.integer(gene_start - element_end),
    default = 0L
  )]
  joined[, mutation_priority := data.table::fcase(
    relationship == "direct_overlap", "high",
    distance_to_element_bp <= 500L, "high",
    distance_to_element_bp <= 2000L, "medium",
    default = "low"
  )]
  joined[, candidate_variant_modes := data.table::fcase(
    relationship == "direct_overlap", "insertional_disruption; excision_footprint",
    distance_to_element_bp <= 500L, "boundary_insertion; promoter_or_operon_effect",
    default = "regional_rearrangement_context"
  )]
  joined[, candidate_position_start := data.table::fcase(
    relationship == "direct_overlap", pmax(gene_start, element_start),
    relationship == "upstream_neighbor", element_start,
    relationship == "downstream_neighbor", element_end,
    default = element_start
  )]
  joined[, candidate_position_end := data.table::fcase(
    relationship == "direct_overlap", pmin(gene_end, element_end),
    relationship == "upstream_neighbor", element_start,
    relationship == "downstream_neighbor", element_end,
    default = element_end
  )]
  joined[, interval_context := data.table::fcase(
    relationship == "direct_overlap", "gene_body",
    distance_to_element_bp <= 150L, "promoter",
    default = "intergenic"
  )]
  joined[, essentiality_penalty := data.table::fcase(
    interval_context == "gene_body" & essentiality_class == "high", 0.95,
    interval_context == "gene_body" & essentiality_class == "medium", 0.55,
    interval_context == "gene_body", 0.2,
    interval_context == "promoter" & essentiality_class == "high", 0.7,
    interval_context == "promoter" & essentiality_class == "medium", 0.4,
    interval_context == "promoter", 0.15,
    interval_context == "intergenic" & essentiality_class == "high", 0.12,
    interval_context == "intergenic" & essentiality_class == "medium", 0.07,
    interval_context == "intergenic" & !is.na(essentiality_class), 0.03,
    interval_context == "intergenic", 0.05,
    default = 0.15
  )]
  joined[, expected_viability := pmax(0, pmin(1, 1 - essentiality_penalty))]
  joined[, element_confidence := confidence]
  joined[, phenotype_categories := gene_phenotype_categories]
  joined[, phenotype_prediction := gene_phenotype_prediction]

  data.table::setorder(joined, distance_to_element_bp, gene_start)

  result <- joined[, .(
    element_id,
    element_family,
    element_confidence,
    contig,
    gene_id,
    locus_tag,
    gene,
    product,
    note,
    start = gene_start,
    end = gene_end,
    strand,
    essentiality_score,
    essentiality_class,
    relationship,
    distance_to_element_bp,
    mutation_priority,
    candidate_variant_modes,
    candidate_position_start,
    candidate_position_end,
    essentiality_penalty,
    expected_viability,
    phenotype_categories,
    phenotype_prediction,
    annotation_text
  )]

  tibble::as_tibble(result)
}

.dnmb_empty_mutation_points <- function() {
  tibble::tibble(
    element_id = character(),
    element_family = character(),
    element_confidence = character(),
    contig = character(),
    gene_id = character(),
    locus_tag = character(),
    gene = character(),
    product = character(),
    note = character(),
    start = integer(),
    end = integer(),
    strand = character(),
    essentiality_score = numeric(),
    essentiality_class = character(),
    relationship = character(),
    distance_to_element_bp = integer(),
    mutation_priority = character(),
    candidate_variant_modes = character(),
    candidate_position_start = integer(),
    candidate_position_end = integer(),
    essentiality_penalty = numeric(),
    expected_viability = numeric(),
    phenotype_categories = character(),
    phenotype_prediction = character(),
    annotation_text = character()
  )
}

.dnmb_predict_variant_scenarios <- function(elements, genes, pair_max_dist_bp = 200000L) {
  if (nrow(elements) < 2L) {
    return(.dnmb_empty_scenarios())
  }

  scenario_rows <- list()
  scenario_id <- 1L

  for (i in seq_len(nrow(elements) - 1L)) {
    for (j in (i + 1L):nrow(elements)) {
      same_contig <- identical(elements$contig[[i]], elements$contig[[j]])
      if (!same_contig) {
        next
      }

      left_idx <- if (elements$start[[i]] <= elements$start[[j]]) i else j
      right_idx <- if (left_idx == i) j else i

      span_bp <- elements$start[[right_idx]] - elements$end[[left_idx]]
      if (is.na(span_bp) || span_bp < 0L || span_bp > pair_max_dist_bp) {
        next
      }

      if (!.dnmb_elements_can_pair(elements[left_idx, , drop = FALSE], elements[right_idx, , drop = FALSE], pair_max_dist_bp)) {
        next
      }

      left <- elements[left_idx, , drop = FALSE]
      right <- elements[right_idx, , drop = FALSE]
      inner_start <- left$end[[1]] + 1L
      inner_end <- right$start[[1]] - 1L
      affected_genes <- genes %>%
        dplyr::filter(
          .data$contig == left$contig[[1]],
          .data$start <= inner_end,
          .data$end >= inner_start
        )

      orientation_class <- if (identical(left$strand[[1]], right$strand[[1]]) && left$strand[[1]] != "mixed") {
        "same_orientation"
      } else {
        "opposite_or_mixed_orientation"
      }

      scenario_type <- if (orientation_class == "same_orientation") {
        "deletion_or_excision"
      } else {
        "inversion"
      }

      scenario_rows[[length(scenario_rows) + 1L]] <- .dnmb_build_scenario_row(
        scenario_id = paste0("SV_", sprintf("%04d", scenario_id)),
        scenario_type = scenario_type,
        left = left,
        right = right,
        inner_start = inner_start,
        inner_end = inner_end,
        affected_genes = affected_genes,
        orientation_class = orientation_class
      )
      scenario_id <- scenario_id + 1L

      if (inner_end > inner_start && (inner_end - inner_start + 1L) <= min(100000L, pair_max_dist_bp)) {
        scenario_rows[[length(scenario_rows) + 1L]] <- .dnmb_build_scenario_row(
          scenario_id = paste0("SV_", sprintf("%04d", scenario_id)),
          scenario_type = "composite_transposition_candidate",
          left = left,
          right = right,
          inner_start = inner_start,
          inner_end = inner_end,
          affected_genes = affected_genes,
          orientation_class = orientation_class
        )
        scenario_id <- scenario_id + 1L
      }
    }
  }

  if (!length(scenario_rows)) {
    return(.dnmb_empty_scenarios())
  }

  dplyr::bind_rows(scenario_rows) %>%
    dplyr::arrange(.data$contig, .data$interval_start, .data$interval_end, .data$scenario_type)
}

.dnmb_empty_scenarios <- function() {
  tibble::tibble(
    scenario_id = character(),
    scenario_type = character(),
    contig = character(),
    source_elements = character(),
    source_element_families = character(),
    orientation_class = character(),
    interval_start = integer(),
    interval_end = integer(),
    interval_size_bp = integer(),
    affected_gene_count = integer(),
    affected_genes = character(),
    affected_products = character(),
    high_essential_gene_count = integer(),
    medium_essential_gene_count = integer(),
    highest_essentiality_class = character(),
    phenotype_categories = character(),
    phenotype_prediction = character(),
    confidence = character(),
    confidence_score = numeric(),
    essentiality_penalty = numeric(),
    viability_score = numeric(),
    viability_class = character(),
    rationale = character()
  )
}

.dnmb_elements_can_pair <- function(left, right, pair_max_dist_bp) {
  families_known <- !is.na(left$element_family[[1]]) && nzchar(left$element_family[[1]]) &&
    !is.na(right$element_family[[1]]) && nzchar(right$element_family[[1]])
  if (families_known) {
    return(identical(left$element_family[[1]], right$element_family[[1]]))
  }

  both_highish <- left$confidence_score[[1]] >= 2L && right$confidence_score[[1]] >= 2L
  gap_bp <- right$start[[1]] - left$end[[1]]
  both_highish && gap_bp <= min(50000L, pair_max_dist_bp)
}

.dnmb_build_scenario_row <- function(
  scenario_id,
  scenario_type,
  left,
  right,
  inner_start,
  inner_end,
  affected_genes,
  orientation_class
) {
  gene_ids <- if (nrow(affected_genes)) paste(unique(stats::na.omit(affected_genes$gene_id)), collapse = "; ") else NA_character_
  products <- if (nrow(affected_genes)) paste(unique(stats::na.omit(affected_genes$product)), collapse = "; ") else NA_character_
  annotation_text <- if (nrow(affected_genes)) paste(affected_genes$annotation_text, collapse = " ") else NA_character_
  phenotype_categories <- .dnmb_collapse_categories(annotation_text)
  phenotype_prediction <- .dnmb_collapse_effects(annotation_text)
  high_essential_gene_count <- if (nrow(affected_genes)) sum(affected_genes$essentiality_class == "high", na.rm = TRUE) else 0L
  medium_essential_gene_count <- if (nrow(affected_genes)) sum(affected_genes$essentiality_class == "medium", na.rm = TRUE) else 0L
  highest_essentiality_class <- if (high_essential_gene_count > 0L) {
    "high"
  } else if (medium_essential_gene_count > 0L) {
    "medium"
  } else if (nrow(affected_genes) > 0L) {
    "low"
  } else {
    NA_character_
  }

  confidence <- if (
    !is.na(left$element_family[[1]]) &&
    !is.na(right$element_family[[1]]) &&
    identical(left$element_family[[1]], right$element_family[[1]])
  ) {
    "high"
  } else if (left$confidence_score[[1]] >= 2L && right$confidence_score[[1]] >= 2L) {
    "medium"
  } else {
    "low"
  }
  confidence_score <- dplyr::case_when(
    confidence == "high" ~ 0.85,
    confidence == "medium" ~ 0.55,
    TRUE ~ 0.3
  )

  essentiality_penalty <- dplyr::case_when(
    scenario_type == "deletion_or_excision" ~ min(1, 0.55 * high_essential_gene_count + 0.2 * medium_essential_gene_count),
    scenario_type == "composite_transposition_candidate" ~ min(1, 0.35 * high_essential_gene_count + 0.15 * medium_essential_gene_count),
    TRUE ~ min(1, 0.25 * high_essential_gene_count + 0.1 * medium_essential_gene_count)
  )
  viability_score <- max(0, min(1, confidence_score - essentiality_penalty))
  viability_class <- dplyr::case_when(
    viability_score >= 0.7 ~ "high",
    viability_score >= 0.35 ~ "medium",
    TRUE ~ "low"
  )

  rationale <- dplyr::case_when(
    scenario_type == "deletion_or_excision" ~ "Elements on the same contig and orientation can support recombination-mediated deletion or excision of the intervening interval.",
    scenario_type == "inversion" ~ "Elements on the same contig with opposite or mixed orientation can support inversion of the intervening interval.",
    TRUE ~ "Two nearby IS/mobile elements flank a gene block that could behave as a composite transposition cargo region."
  )

  tibble::tibble(
    scenario_id = scenario_id,
    scenario_type = scenario_type,
    contig = left$contig[[1]],
    source_elements = paste(left$element_id[[1]], right$element_id[[1]], sep = "; "),
    source_element_families = paste(
      dplyr::coalesce(left$element_family[[1]], "NA"),
      dplyr::coalesce(right$element_family[[1]], "NA"),
      sep = "; "
    ),
    orientation_class = orientation_class,
    interval_start = inner_start,
    interval_end = inner_end,
    interval_size_bp = max(0L, inner_end - inner_start + 1L),
    affected_gene_count = nrow(affected_genes),
    affected_genes = gene_ids,
    affected_products = products,
    high_essential_gene_count = high_essential_gene_count,
    medium_essential_gene_count = medium_essential_gene_count,
    highest_essentiality_class = highest_essentiality_class,
    phenotype_categories = phenotype_categories,
    phenotype_prediction = phenotype_prediction,
    confidence = confidence,
    confidence_score = confidence_score,
    essentiality_penalty = essentiality_penalty,
    viability_score = viability_score,
    viability_class = viability_class,
    rationale = rationale
  )
}

.dnmb_build_phenotype_predictions <- function(elements, mutation_points, scenarios, variant_catalog = NULL) {
  element_context <- if (nrow(mutation_points)) {
    mutation_points %>%
      dplyr::group_by(.data$element_id, .data$element_family, .data$element_confidence, .data$contig) %>%
      dplyr::summarise(
        context_type = "local_element_context",
        affected_gene_count = dplyr::n_distinct(.data$gene_id),
        affected_genes = paste(unique(stats::na.omit(.data$gene_id)), collapse = "; "),
        phenotype_categories = .dnmb_collapse_categories(paste(.data$annotation_text, collapse = " ")),
        phenotype_prediction = .dnmb_collapse_effects(paste(.data$annotation_text, collapse = " ")),
        supporting_rationale = "Based on genes overlapping or flanking the candidate element.",
        confidence = dplyr::first(.data$element_confidence),
        viability = dplyr::case_when(
          min(.data$expected_viability, na.rm = TRUE) >= 0.7 ~ "high",
          min(.data$expected_viability, na.rm = TRUE) >= 0.35 ~ "medium",
          TRUE ~ "low"
        ),
        .groups = "drop"
      ) %>%
      dplyr::rename(context_id = element_id)
  } else {
    tibble::tibble(
      context_id = character(),
      element_family = character(),
      element_confidence = character(),
      contig = character(),
      context_type = character(),
      affected_gene_count = integer(),
      affected_genes = character(),
      phenotype_categories = character(),
      phenotype_prediction = character(),
      supporting_rationale = character(),
      confidence = character(),
      viability = character()
    )
  }

  scenario_context <- if (nrow(scenarios)) {
    scenarios %>%
      dplyr::transmute(
        context_id = scenario_id,
        element_family = source_element_families,
        element_confidence = NA_character_,
        contig = contig,
        context_type = scenario_type,
        affected_gene_count = affected_gene_count,
        affected_genes = affected_genes,
        phenotype_categories = phenotype_categories,
        phenotype_prediction = phenotype_prediction,
        supporting_rationale = rationale,
        confidence = confidence,
        viability = viability_class
      )
  } else {
    tibble::tibble(
      context_id = character(),
      element_family = character(),
      element_confidence = character(),
      contig = character(),
      context_type = character(),
      affected_gene_count = integer(),
      affected_genes = character(),
      phenotype_categories = character(),
      phenotype_prediction = character(),
      supporting_rationale = character(),
      confidence = character(),
      viability = character()
    )
  }

  variant_context <- if (!is.null(variant_catalog) && nrow(variant_catalog)) {
    variant_catalog %>%
      dplyr::group_by(.data$variant_type, .data$source_family) %>%
      dplyr::summarise(
        context_id = paste0("variant:", dplyr::first(.data$variant_type), ":", dplyr::first(.data$source_family)),
        element_family = dplyr::first(.data$source_family),
        element_confidence = NA_character_,
        contig = paste(unique(stats::na.omit(.data$contig)), collapse = "; "),
        context_type = paste0("variant_catalog_", dplyr::first(.data$variant_type)),
        affected_gene_count = dplyr::n_distinct(stats::na.omit(.data$affected_genes)),
        affected_genes = paste(unique(stats::na.omit(.data$affected_genes)), collapse = "; "),
        phenotype_categories = .dnmb_collapse_categories(paste(.data$phenotype_prediction, collapse = " ")),
        phenotype_prediction = paste(unique(stats::na.omit(.data$phenotype_prediction)), collapse = "; "),
        supporting_rationale = "Aggregated from scored target-site and structural-variant predictions.",
        confidence = dplyr::case_when(
          mean(.data$likelihood_score, na.rm = TRUE) >= 0.7 ~ "high",
          mean(.data$likelihood_score, na.rm = TRUE) >= 0.35 ~ "medium",
          TRUE ~ "low"
        ),
        viability = dplyr::case_when(
          mean(.data$viability_score, na.rm = TRUE) >= 0.7 ~ "high",
          mean(.data$viability_score, na.rm = TRUE) >= 0.35 ~ "medium",
          TRUE ~ "low"
        ),
        .groups = "drop"
      )
  } else {
    tibble::tibble(
      context_id = character(),
      element_family = character(),
      element_confidence = character(),
      contig = character(),
      context_type = character(),
      affected_gene_count = integer(),
      affected_genes = character(),
      phenotype_categories = character(),
      phenotype_prediction = character(),
      supporting_rationale = character(),
      confidence = character(),
      viability = character()
    )
  }

  dplyr::bind_rows(element_context, scenario_context, variant_context)
}

.dnmb_collapse_categories <- function(annotation_text) {
  categories <- .dnmb_predict_phenotype_categories(annotation_text)
  if (!length(categories)) {
    return(NA_character_)
  }
  paste(categories, collapse = "; ")
}

.dnmb_collapse_effects <- function(annotation_text) {
  categories <- .dnmb_predict_phenotype_categories(annotation_text)
  if (!length(categories)) {
    return("general functional change possible, but phenotype category is not specific from the current annotation")
  }
  effects <- unique(unname(.dnmb_phenotype_effect_map()[categories]))
  paste(effects, collapse = "; ")
}

.dnmb_predict_phenotype_categories <- function(annotation_text) {
  text <- stringr::str_to_lower(paste(annotation_text, collapse = " "))
  if (!nzchar(text)) {
    return(character())
  }

  category_patterns <- .dnmb_phenotype_category_patterns()

  hits <- names(category_patterns)[vapply(
    category_patterns,
    function(pattern) stringr::str_detect(text, stringr::regex(pattern, ignore_case = TRUE)),
    FUN.VALUE = logical(1)
  )]
  hits
}

.dnmb_phenotype_category_patterns <- function() {
  list(
    antimicrobial_resistance = "resistan|efflux|beta-lactamase|multidrug|aminoglycoside|vancomycin|tetracycline|macrolide|drug exporter",
    motility_chemotaxis = "flagell|chemotaxis|motility|fli[[:alnum:]]|flg[[:alnum:]]|mota|motb",
    biofilm_attachment = "biofilm|adhesin|fimbr|pilus|pili|exopolysaccharide|capsule biosynthesis|pel|psl|alginate|curli",
    sporulation_development = "sporulat|spore|germin|spo[[:alnum:]]|sigf|sige|sigh|sigk",
    secretion_conjugation = "secretion system|type iii secretion|type iv secretion|type vi secretion|sec pathway|tat pathway|conjug",
    cell_envelope = "cell wall|peptidoglycan|murein|lipopolysaccharide|outer membrane|surface protein|capsule|porin",
    stress_response = "stress|heat shock|cold shock|oxidative|superoxide|catalase|peroxidase|dna repair protein reca|dna damage|chaperone|dnak|groel|usp",
    genome_stability = "dna repair|recombination|replication|helicase|polymerase|muts|mutl|uvr|ruv|reca|par[ab]",
    metabolism = "dehydrogenase|synthase|hydrolase|kinase|transferase|isomerase|lyase|metabolism|catabolic|biosynth|ferment|respirat",
    transport = "transporter|permease|abc transporter|mfs transporter|symporter|antiporter|channel protein|porin",
    transcriptional_regulation = "transcriptional regulator|response regulator|repressor|activator|sigma factor|sensor kinase|two-component|dna-binding",
    virulence_host_interaction = "virulence|toxin|hemolysin|invasion|pathogenic|host interaction|adhesin",
    defense_anti_phage = "restriction|methylase|crispr|cas[[:digit:]]|abortive infection|defense system|rm system",
    growth_division_membrane = "cell division|ftsz|ftsa|mreb|rodz|membrane protein|atpase|lipid|fatty acid"
  )
}

.dnmb_predict_phenotype_strings_vector <- function(annotation_texts) {
  texts <- stringr::str_to_lower(dplyr::coalesce(annotation_texts, ""))
  patterns <- .dnmb_phenotype_category_patterns()
  if (!length(texts)) {
    return(list(categories = character(), predictions = character()))
  }

  hit_matrix <- vapply(
    patterns,
    function(pattern) grepl(pattern, texts, perl = TRUE),
    FUN.VALUE = logical(length(texts))
  )
  if (is.null(dim(hit_matrix))) {
    hit_matrix <- matrix(hit_matrix, ncol = length(patterns))
    colnames(hit_matrix) <- names(patterns)
  }

  categories <- vapply(seq_len(nrow(hit_matrix)), function(i) {
    hits <- colnames(hit_matrix)[hit_matrix[i, ]]
    if (!length(hits)) NA_character_ else paste(hits, collapse = "; ")
  }, character(1))

  effect_map <- .dnmb_phenotype_effect_map()
  predictions <- vapply(seq_along(categories), function(i) {
    if (is.na(categories[[i]]) || !nzchar(categories[[i]])) {
      return("general functional change possible, but phenotype category is not specific from the current annotation")
    }
    hits <- strsplit(categories[[i]], "; ", fixed = TRUE)[[1]]
    effects <- unique(unname(effect_map[hits]))
    paste(effects, collapse = "; ")
  }, character(1))

  list(categories = categories, predictions = predictions)
}

.dnmb_phenotype_effect_map <- function() {
  c(
    antimicrobial_resistance = "altered antimicrobial susceptibility",
    motility_chemotaxis = "altered motility or chemotaxis",
    biofilm_attachment = "altered biofilm formation or surface attachment",
    sporulation_development = "altered sporulation, germination, or developmental timing",
    secretion_conjugation = "altered secretion or horizontal-transfer capacity",
    cell_envelope = "altered cell envelope structure or permeability",
    stress_response = "altered stress tolerance",
    genome_stability = "altered genome stability, DNA repair, or chromosome maintenance",
    metabolism = "altered nutrient use or metabolic flux",
    transport = "altered uptake or export of metabolites and drugs",
    transcriptional_regulation = "altered regulatory network behavior",
    virulence_host_interaction = "altered host interaction or virulence potential",
    defense_anti_phage = "altered defense against phage or foreign DNA",
    growth_division_membrane = "altered growth, morphology, or membrane homeostasis"
  )
}
