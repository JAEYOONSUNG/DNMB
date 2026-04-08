suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(tidyr)
  library(cowplot); library(circlize); library(ComplexHeatmap)
  library(grid); library(gridBase); library(gggenes)
})
dnmb_iselement_context <- get0("dnmb_iselement_context", inherits = TRUE, ifnotfound = list())
.dnmb_iselement_ctx_get <- function(name, default = NULL) {
  value <- dnmb_iselement_context[[name]]
  if (is.null(value)) default else value
}

pkg_dir <- .dnmb_iselement_ctx_get("pkg_dir", normalizePath(getwd(), winslash = "/", mustWork = FALSE))
output_dir <- .dnmb_iselement_ctx_get("output_dir", getwd())
plot_dir <- .dnmb_iselement_ctx_get("plot_dir", file.path(output_dir, "visualizations"))
genbank_path <- .dnmb_iselement_ctx_get("genbank_path", file.path(output_dir, "GCF_030369615.1.gbff"))
full_result_path <- .dnmb_iselement_ctx_get("full_result_path", file.path(output_dir, "iselement_full_result.rds"))
landing_pads_path <- .dnmb_iselement_ctx_get("landing_pads_path", file.path(output_dir, "dnmb_module_iselement", "iselement_landing_pads.tsv"))
tsd_path <- .dnmb_iselement_ctx_get("tsd_path", file.path(output_dir, "dnmb_module_iselement", "iselement_tsd_sequences.tsv"))
verify_tsd_path <- .dnmb_iselement_ctx_get("verify_tsd_path", file.path(output_dir, "dnmb_module_iselement", "insertion_verify", "validated_tsd.tsv"))
iv_verify_path <- .dnmb_iselement_ctx_get("insertion_verify_path", file.path(output_dir, "insertion_verify_result.rds"))
iv_verify_tsv <- .dnmb_iselement_ctx_get("insertion_verify_tsv_path", file.path(output_dir, "dnmb_module_iselement", "insertion_verify", "insertion_verification.tsv"))
out_file <- .dnmb_iselement_ctx_get("out_file", file.path(plot_dir, "ISelement_overview.pdf"))

# Source only if functions not already available (standalone mode)
if (!exists(".dnmb_parse_genbank_features", mode = "function")) {
  for (f in c("R/db_modules.R","R/module_api.R","R/Mobileome_pipeline.R",
              "R/Mobileome_sequence_engine.R","R/Mobileome_variant_engine.R",
              "R/Mobileome_comparative.R","R/Mobileome_landing_pads.R",
              "R/Mobileome_auto_comparative.R","R/Mobileome_insertion_verify.R",
              "R/db_module_iselement.R"))
    if (file.exists(file.path(pkg_dir, f))) source(file.path(pkg_dir, f), local = FALSE)
}

parsed <- .dnmb_iselement_ctx_get("parsed", .dnmb_parse_genbank_features(genbank_path))
genes <- .dnmb_iselement_ctx_get("genes", .dnmb_predict_gene_essentiality(.dnmb_build_gene_table(parsed$features)))
result <- .dnmb_iselement_ctx_get("full_result", readRDS(full_result_path))
elements <- .dnmb_iselement_ctx_get("elements", result$elements)
census <- .dnmb_iselement_ctx_get("census", result$census)
target_models <- .dnmb_iselement_ctx_get("target_models", result$target_models)
landing_pads <- .dnmb_iselement_ctx_get("landing_pads", read.delim(landing_pads_path, check.names=FALSE))
genome_len <- max(genes$end, na.rm=TRUE)
ctg <- .dnmb_iselement_ctx_get("main_contig", unique(genes$contig)[1])

# Use ALL contigs (chromosome + plasmid)
all_contigs <- unique(genes$contig)
ess_track <- genes
is_track <- elements %>%
  mutate(fam = ifelse(is.na(element_family)|!nzchar(element_family), "unknown", element_family))
lp_track <- landing_pads
fam_cols <- c(IS110="#7E57C2",IS982="#26A69A",IS701="#EF5350",IS630="#42A5F5",
              IS66="#FFA726",IS4="#AB47BC",IS21="#26C6DA",IS3="#EC407A",
              IS481="#5C6BC0",ISL3="#66BB6A","IS200/IS605"="#8D6E63",
              IS256="#78909C",IS5="#D4E157",IS1595="#FF7043",unknown="#BDBDBD")

# ======================================================================
# Pre-render unified Panel B: one row per IS family
# Each row: [Col1: Family info + bar | Col2: TSD SeqLogo | Col3: gggenes synteny]
# ======================================================================

# --- Load TSD data ---
tsd_raw <- if (file.exists(tsd_path)) read.delim(tsd_path, check.names = FALSE) else data.frame()

if (file.exists(verify_tsd_path)) {
  verify_tsd <- read.delim(verify_tsd_path, check.names = FALSE)
  if (nrow(verify_tsd) && !"source" %in% names(verify_tsd))
    verify_tsd$source <- "comparative_verified"
  if (nrow(verify_tsd)) {
    verify_for_merge <- verify_tsd[, intersect(names(tsd_raw), names(verify_tsd)), drop = FALSE]
    if (!"source" %in% names(tsd_raw) && nrow(tsd_raw))
      tsd_raw$source <- "boundary_duplication"
    tsd_raw <- dplyr::bind_rows(tsd_raw, verify_for_merge)
    tsd_raw <- tsd_raw[!duplicated(tsd_raw[, c("family", "tsd_seq")]), , drop = FALSE]
  }
}

# --- Load insertion verification ---
insertion_verify_df <- NULL
if (file.exists(iv_verify_path)) {
  iv_obj <- readRDS(iv_verify_path)
  if (!is.null(iv_obj$verification) && nrow(iv_obj$verification))
    insertion_verify_df <- iv_obj$verification
} else if (file.exists(iv_verify_tsv)) {
  insertion_verify_df <- read.delim(iv_verify_tsv, check.names = FALSE)
}

# --- Re-run insertion verification (if related genomes cached) ---
ref_cache_dir <- file.path(output_dir, "dnmb_module_iselement", "related_genomes_cache")
ref_genbanks <- if (dir.exists(ref_cache_dir)) {
  list.files(ref_cache_dir, pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE)
} else character(0)

# Use cached verification results if available (skip expensive re-run)
iv_cache <- iv_verify_path
iv_updated <- if (file.exists(iv_cache)) readRDS(iv_cache) else NULL
if (!is.null(iv_updated) && !is.null(iv_updated$verification) && nrow(iv_updated$verification)) {
  insertion_verify_df <- iv_updated$verification
}
iv_for_panel <- insertion_verify_df

# --- Target models for TSD lengths ---
tm <- target_models %>% filter(n_elements >= 2) %>% arrange(desc(n_elements))

# --- Family selection: top families with empty_site verification ---
fam_counts_df <- elements %>%
  filter(!is.na(element_family) & nzchar(element_family)) %>%
  count(element_family, sort = TRUE)

synteny_families <- character()
empty_site_df <- data.frame()
if (!is.null(iv_for_panel) && nrow(iv_for_panel) &&
    "left_ref_end" %in% names(iv_for_panel)) {
  empty_site_df <- iv_for_panel %>%
    filter(status == "empty_site",
           !is.na(left_ref_end), !is.na(right_ref_start),
           !is.na(element_family) & nzchar(element_family)) %>%
    group_by(element_family) %>%
    slice_max(order_by = left_pident + right_pident, n = 1) %>%
    ungroup()
  synteny_families <- intersect(fam_counts_df$element_family, empty_site_df$element_family)
}

# Show ALL families with >= 2 elements, sorted by count descending.
# Include "unknown" family as well.
all_fam_counts <- elements %>%
  mutate(element_family = ifelse(is.na(element_family) | !nzchar(element_family),
                                  "unknown", element_family)) %>%
  count(element_family, sort = TRUE) %>%
  filter(n >= 1)
show_families <- all_fam_counts$element_family
n_flanking <- 4  # Show 4 genes on each side of IS element
flank_bp <- 10000L

# --- Parse reference genomes (genes + organism info) ---
ref_gene_tables <- list()
ref_organism_info <- list()
for (rp in ref_genbanks) {
  ref_parsed_tmp <- tryCatch(.dnmb_parse_genbank_features(rp), error = function(e) NULL)
  if (!is.null(ref_parsed_tmp)) {
    ref_genes_tmp <- .dnmb_build_gene_table(ref_parsed_tmp$features)
    ref_gene_tables[[basename(rp)]] <- ref_genes_tmp
    # Parse organism from SOURCE line
    ref_lines <- tryCatch(readLines(rp, n = 30, warn = FALSE), error = function(e) character())
    source_idx <- grep("^SOURCE\\s+", ref_lines)
    accession_idx <- grep("^ACCESSION\\s+", ref_lines)
    organism <- if (length(source_idx)) trimws(sub("^SOURCE\\s+", "", ref_lines[source_idx[1]])) else basename(rp)
    accession <- if (length(accession_idx)) trimws(sub("^ACCESSION\\s+", "", ref_lines[accession_idx[1]])) else ""
    ref_organism_info[[basename(rp)]] <- list(organism = organism, accession = accession)
  }
}

# --- Helpers ---
.short_gene_label <- function(gene_name, product, max_len = 12) {
  if (!is.na(gene_name) && nzchar(gene_name)) return(gene_name)
  if (!is.na(product) && nzchar(product)) return(substr(product, 1, max_len))
  return("")
}
.get_ess_type <- function(ess_score, g_start, g_end, el_start, el_end) {
  if (!is.na(el_start) && !is.na(el_end) && g_end >= el_start && g_start <= el_end) return("IS-associated")
  if (!is.na(ess_score) && ess_score >= 0.5) return("essential")
  return("non-essential")
}

# Color palette
type_cols <- c(
  "IS_element" = "#E53935", "essential" = "#1565C0", "non-essential" = "#CFD8DC",
  "IS-associated" = "#FF8F00", "ref_gene" = "#90A4AE", "empty_site" = "white"
)
type_labels <- c(
  "IS_element" = "IS Element", "essential" = "Essential gene",
  "non-essential" = "Non-essential", "IS-associated" = "IS-associated",
  "ref_gene" = "Reference gene", "empty_site" = "Empty site"
)
conf_cols <- c(high = "#66BB6A", medium = "#42A5F5", low = "#FFC107")

Y_FOCAL <- 0.5
Y_REF   <- 0

# ======================================================================
# Pre-pass: compute unified x-range based on max IS length + flank_bp
# All families will show: [-flank_bp - max_IS/2, +flank_bp + max_IS/2]
# Fixed total view = ±10kb from IS center = 20kb total for all families
# ======================================================================
unified_x_half <- flank_bp  # exactly 10kb each side = 20kb total
unified_xlim <- c(-unified_x_half, unified_x_half)

# ======================================================================
# Build unified rows: one per family
# ======================================================================
unified_rows <- list()
# Shared x-axis scale: max element count across all shown families
max_fam_count <- max(all_fam_counts$n[all_fam_counts$element_family %in% show_families], na.rm = TRUE)

for (fam_name in show_families) {
  fam_count <- all_fam_counts$n[all_fam_counts$element_family == fam_name]
  if (is.na(fam_count) || length(fam_count) == 0) fam_count <- 0
  fam_census <- census %>% filter(family == fam_name) %>% slice(1)
  if (fam_name == "unknown") fam_census <- census %>% filter(family == fam_name | is.na(family)) %>% slice(1)
  has_empty_site <- fam_name %in% empty_site_df$element_family
  rep_verify <- if (has_empty_site) empty_site_df %>% filter(element_family == fam_name) %>% slice(1) else data.frame()
  rep_id <- if (has_empty_site) rep_verify$element_id else NA

  # ── Column 1: Family info panel ──
  # Confidence distribution bar data
  hi_n <- if (nrow(fam_census) && !is.null(fam_census$high_confidence_count)) fam_census$high_confidence_count else 0
  md_n <- if (nrow(fam_census) && !is.null(fam_census$medium_confidence_count)) fam_census$medium_confidence_count else 0
  lo_n <- if (nrow(fam_census) && !is.null(fam_census$low_confidence_count)) fam_census$low_confidence_count else 0
  # For unknown family, estimate from elements directly if census is missing

  if ((hi_n + md_n + lo_n) == 0 && fam_count > 0) {
    fam_els <- elements %>%
      mutate(element_family = ifelse(is.na(element_family) | !nzchar(element_family), "unknown", element_family)) %>%
      filter(element_family == fam_name)
    if ("confidence" %in% names(fam_els)) {
      hi_n <- sum(fam_els$confidence == "high", na.rm = TRUE)
      md_n <- sum(fam_els$confidence == "medium", na.rm = TRUE)
      lo_n <- sum(fam_els$confidence == "low", na.rm = TRUE)
    } else {
      lo_n <- fam_count
    }
  }
  conf_bar <- data.frame(
    lev = factor(c("high", "medium", "low"), levels = c("low", "medium", "high")),
    cnt = c(hi_n, md_n, lo_n), stringsAsFactors = FALSE
  )

  # Verification status for this family
  fam_verify_all <- if (!is.null(iv_for_panel) && nrow(iv_for_panel) && fam_name != "unknown") {
    iv_for_panel %>% filter(element_family == fam_name)
  } else data.frame()
  n_empty_fam <- sum(fam_verify_all$status == "empty_site", na.rm = TRUE)
  n_filled_fam <- sum(fam_verify_all$status == "filled_site", na.rm = TRUE)
  n_no_hit_fam <- sum(fam_verify_all$status == "no_hit", na.rm = TRUE)

  # Build status label and color
  if (n_empty_fam > 0) {
    status_label <- paste0("empty_site (", n_empty_fam, ")")
    status_col <- "#2E7D32"
    status_bg <- "#E8F5E9"
  } else if (n_filled_fam > 0) {
    status_label <- paste0("filled_site (", n_filled_fam, ")")
    status_col <- "#E65100"
    status_bg <- "#FFF3E0"
  } else if (n_no_hit_fam > 0) {
    status_label <- paste0("no_hit (", n_no_hit_fam, ")")
    status_col <- "#616161"
    status_bg <- "#F5F5F5"
  } else {
    status_label <- "not verified"
    status_col <- "#9E9E9E"
    status_bg <- "#FAFAFA"
  }

  # Build Column 1: horizontal stacked bar (x = count, shared scale across families)
  bar_data <- data.frame(
    family = fam_name,
    confidence = factor(c("low", "medium", "high"), levels = c("low", "medium", "high")),
    count = c(lo_n, md_n, hi_n)
  )

  # Column 1: Horizontal stacked bar with family name as rotated y-axis label
  p_info <- ggplot(bar_data, aes(y = fam_name, x = count, fill = confidence)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c(low = "#FFC107", medium = "#42A5F5", high = "#66BB6A")) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              position = position_stack(vjust = 0.5), size = 2.3, color = "white", fontface = "bold") +
    annotate("text", x = fam_count + max_fam_count * 0.02, y = fam_name,
             label = paste0("n=", fam_count), hjust = 0, size = 2.3, color = "grey30") +
    scale_x_continuous(limits = c(0, max_fam_count * 1.12), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "bold", size = 10, angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = 6),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(1, 4, 1, 2))

  # ── Column 2: TSD SeqLogo ──
  fam_tm <- tm %>% filter(family == fam_name) %>% slice(1)
  dom_len <- if (nrow(fam_tm)) fam_tm$dominant_tsd_len else NA

  fam_tsd <- if (nrow(tsd_raw)) tsd_raw[tsd_raw$family == fam_name, , drop = FALSE] else data.frame()
  p_logo <- NULL

  if (nrow(fam_tsd) && !is.na(dom_len)) {
    fam_tsd_dom <- fam_tsd[fam_tsd$tsd_len_bp == dom_len, , drop = FALSE]
    seqs <- toupper(fam_tsd_dom$tsd_seq[!is.na(fam_tsd_dom$tsd_seq) &
                                         nzchar(fam_tsd_dom$tsd_seq) &
                                         grepl("^[ACGT]+$", fam_tsd_dom$tsd_seq)])
    n_bd <- sum(fam_tsd_dom$source == "boundary_duplication", na.rm = TRUE)
    n_cv <- sum(fam_tsd_dom$source == "comparative_verified", na.rm = TRUE)

    if (length(seqs) >= 1) {
      src_tag <- if (n_cv > 0) paste0("bd=", n_bd, " cv=", n_cv) else paste0("n=", length(seqs))
      p_logo <- ggplot() +
        ggseqlogo::geom_logo(seqs, method = "prob") +
        labs(x = NULL, y = NULL) +
        annotate("text", x = (dom_len + 1) / 2, y = 1.05,
                 label = paste0("TSD ", dom_len, "bp (", src_tag, ")"),
                 size = 2.3, color = if (n_cv > 0) "#1B5E20" else "grey30",
                 fontface = "bold") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_void(base_size = 8) +
        theme(plot.background = element_rect(fill = "white", color = NA),
              plot.margin = margin(0, 2, 0, 2))
    }
  }

  if (is.null(p_logo)) {
    p_logo <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No TSD",
               size = 2.5, color = "grey50", fontface = "italic") +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            plot.margin = margin(0, 2, 0, 2))
  }

  # ── Column 3: gggenes synteny view (focal vs reference) ──
  # Pick representative element: use verified element if available, else pick largest
  if (!is.na(rep_id)) {
    rep_el <- elements %>% filter(element_id == rep_id) %>% slice(1)
  } else {
    rep_el <- elements %>%
      mutate(element_family = ifelse(is.na(element_family) | !nzchar(element_family), "unknown", element_family)) %>%
      filter(element_family == fam_name) %>%
      mutate(el_len = end - start) %>%
      slice_max(order_by = el_len, n = 1) %>%
      slice(1)
  }
  el_contig <- rep_el$contig; el_start <- rep_el$start; el_end <- rep_el$end
  el_len <- el_end - el_start + 1
  # Center on IS element so it sits at x=0; show actual position as annotation
  is_center <- (el_start + el_end) / 2

  contig_genes <- genes %>% filter(contig == el_contig) %>% arrange(start)
  p_synteny <- NULL

  if (nrow(contig_genes)) {
    # Select all genes within ±flank_bp of IS center, clamp at boundary
    view_left <- is_center - flank_bp
    view_right <- is_center + flank_bp
    flank_genes <- contig_genes[
      contig_genes$end >= view_left & contig_genes$start <= view_right, , drop = FALSE]

    # Focal rows — skip IS-overlapping genes, clamp boundary genes
    focal_rows <- list()
    for (j in seq_len(nrow(flank_genes))) {
      g <- flank_genes[j, ]
      # Skip genes that overlap with IS element
      if (g$end >= el_start && g$start <= el_end) next
      # Clamp gene coordinates to view boundary (cut as rectangle)
      g_start <- max(g$start, view_left)
      g_end <- min(g$end, view_right)
      # Don't label genes that are clamped (cut at boundary)
      is_clamped <- (g$start < view_left) || (g$end > view_right)
      focal_rows[[length(focal_rows) + 1]] <- data.frame(
        row = "focal",
        xmin = g_start - is_center, xmax = g_end - is_center,
        forward = g$strand == "+",
        gene_label = if (is_clamped) "" else dplyr::coalesce(g$locus_tag, g$gene_id),
        type = .get_ess_type(g$essentiality_score, g$start, g$end, el_start, el_end),
        y = Y_FOCAL, stringsAsFactors = FALSE)
    }
    # IS element with product name
    is_product <- if (!is.na(rep_el$supporting_products) && nzchar(rep_el$supporting_products))
      paste0(fam_name, " (", substr(rep_el$supporting_products, 1, 25), ")") else fam_name
    focal_rows[[length(focal_rows) + 1]] <- data.frame(
      row = "focal",
      xmin = el_start - is_center, xmax = el_end - is_center,
      forward = if (!is.na(rep_el$strand) && rep_el$strand == "-") FALSE else TRUE,
      gene_label = is_product, type = "IS_element",
      y = Y_FOCAL, stringsAsFactors = FALSE)
    focal_df <- do.call(rbind, focal_rows)
    is_xmin <- el_start - is_center
    is_xmax <- el_end - is_center

    # Reference rows — only if this family has empty_site verification
    ref_rows <- list()
    ref_ribbon_data <- NULL

    if (has_empty_site) {
      ref_file <- rep_verify$ref_file
      ref_contig_name <- rep_verify$ref_contig
      left_ref_end_pos <- rep_verify$left_ref_end
      right_ref_start_pos <- rep_verify$right_ref_start

      ref_genes_all <- ref_gene_tables[[ref_file]]
      ref_contig_genes <- if (!is.null(ref_genes_all)) {
        ref_genes_all %>% filter(contig == ref_contig_name) %>% arrange(start)
      } else data.frame()

      if (nrow(ref_contig_genes)) {
        ref_gap_center <- (left_ref_end_pos + right_ref_start_pos) / 2
        # Align reference LEFT boundary with focal IS LEFT boundary
        # is_xmin is focal IS left edge in relative coords
        # left_ref_end_pos is reference gap left edge in absolute coords
        # Shift reference so left_ref_end_pos maps to is_xmin
        ref_offset <- is_xmin - (left_ref_end_pos - ref_gap_center)

        # Select ref genes within ±flank_bp of gap center
        ref_view_left <- ref_gap_center - flank_bp
        ref_view_right <- ref_gap_center + flank_bp
        ref_flank <- ref_contig_genes[
          ref_contig_genes$end >= ref_view_left &
          ref_contig_genes$start <= ref_view_right, , drop = FALSE]

        # Filter out abnormally large genes (annotation artifacts spanning whole contigs)
        ref_flank <- ref_flank[(ref_flank$end - ref_flank$start + 1) < 50000, , drop = FALSE]

        for (j in seq_len(nrow(ref_flank))) {
          g <- ref_flank[j, ]
          rx_min <- (g$start - ref_gap_center) + ref_offset
          rx_max <- (g$end - ref_gap_center) + ref_offset
          # Skip reference genes that fall entirely outside unified_xlim
          if (rx_max < unified_xlim[1] || rx_min > unified_xlim[2]) next
          # Clamp coordinates to unified_xlim, don't label clamped genes
          is_clamped_ref <- (rx_min < unified_xlim[1]) || (rx_max > unified_xlim[2])
          rx_min <- max(rx_min, unified_xlim[1])
          rx_max <- min(rx_max, unified_xlim[2])
          ref_rows[[length(ref_rows) + 1]] <- data.frame(
            row = "reference",
            xmin = rx_min,
            xmax = rx_max,
            forward = g$strand == "+",
            gene_label = if (is_clamped_ref) "" else dplyr::coalesce(g$locus_tag, g$gene_id),
            type = "ref_gene",
            y = Y_REF, stringsAsFactors = FALSE)
        }

        ref_gap_xmin <- (left_ref_end_pos - ref_gap_center) + ref_offset
        ref_gap_xmax <- (right_ref_start_pos - ref_gap_center) + ref_offset

        # --- Ribbon polygons: seamless trapezoids connecting focal and ref rows ---
        # No y-offset needed — polygons span the full gap between Y_FOCAL and Y_REF
        focal_left  <- focal_df %>% filter(xmax <= is_xmin & type != "IS_element")
        focal_right <- focal_df %>% filter(xmin >= is_xmax & type != "IS_element")
        ref_df_tmp  <- do.call(rbind, ref_rows)
        ref_left  <- ref_df_tmp %>% filter(xmax <= ref_gap_xmin)
        ref_right <- ref_df_tmp %>% filter(xmin >= ref_gap_xmax)

        if (nrow(focal_left) && nrow(ref_left)) {
          focal_left_xmin <- max(min(focal_left$xmin), unified_xlim[1])
          ref_left_xmin   <- max(min(ref_left$xmin), unified_xlim[1])
          ref_ribbon_data <- list()
          ref_ribbon_data$left <- data.frame(
            x = c(focal_left_xmin, is_xmin, ref_gap_xmin, ref_left_xmin),
            y = c(Y_FOCAL, Y_FOCAL, Y_REF, Y_REF))
        }
        if (nrow(focal_right) && nrow(ref_right)) {
          focal_right_xmax <- min(max(focal_right$xmax), unified_xlim[2])
          ref_right_xmax   <- min(max(ref_right$xmax), unified_xlim[2])
          if (is.null(ref_ribbon_data)) ref_ribbon_data <- list()
          ref_ribbon_data$right <- data.frame(
            x = c(is_xmax, focal_right_xmax, ref_right_xmax, ref_gap_xmax),
            y = c(Y_FOCAL, Y_FOCAL, Y_REF, Y_REF))
        }
      }
    }

    ref_df <- if (length(ref_rows)) do.call(rbind, ref_rows) else NULL
    gene_df <- rbind(focal_df, ref_df)

    # Flank identity labels
    contig_short <- gsub("^NZ_", "", el_contig)
    gap_label <- ""

    if (has_empty_site) {
      pident_left  <- round(rep_verify$left_pident, 1)
      pident_right <- round(rep_verify$right_pident, 1)
      ref_short <- gsub("^GCF_", "", gsub("\\.gbff$|\\.1$", "", ref_file))

      # TSD annotation
      tsd_text <- NULL
      if (!is.null(iv_updated) && nrow(iv_updated$validated_tsd)) {
        tsd_row <- iv_updated$validated_tsd %>%
          filter(element_id == rep_id) %>% slice(1)
        if (nrow(tsd_row) && !is.na(tsd_row$tsd_seq))
          tsd_text <- paste0("TSD: ", tsd_row$tsd_seq, " (", nchar(tsd_row$tsd_seq), " bp)")
      }

      gap_bp <- rep_verify$gap_bp
      gap_label <- if (!is.null(tsd_text)) {
        paste0(pident_left, "%  |  ", pident_right, "%    ", tsd_text)
      } else if (!is.na(gap_bp) && gap_bp == 0) {
        paste0(pident_left, "%  |  ", pident_right, "%    gap: 0 bp (clean insertion)")
      } else {
        paste0(pident_left, "%  |  ", pident_right, "%    gap: ", gap_bp, " bp")
      }
    }

    # Build synteny plot
    p_synteny <- ggplot()

    # Genome backbone lines (central axis through each row)
    x_full_range <- range(c(gene_df$xmin, gene_df$xmax))
    # Genome backbone lines — uniform linewidth
    x_full_range <- range(c(gene_df$xmin, gene_df$xmax))
    focal_start_bp <- round(x_full_range[1] + is_center)
    focal_end_bp <- round(x_full_range[2] + is_center)
    backbone_df <- data.frame(x = x_full_range[1], xend = x_full_range[2], y = Y_FOCAL, yend = Y_FOCAL)
    p_synteny <- p_synteny +
      geom_segment(data = backbone_df, aes(x = x, xend = xend, y = y, yend = yend),
                   color = "grey60", linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = x_full_range[1], y = Y_FOCAL - 0.15,
               label = format(focal_start_bp, big.mark = ","),
               size = 1.8, color = "grey45", hjust = 0) +
      annotate("text", x = x_full_range[2], y = Y_FOCAL - 0.15,
               label = format(focal_end_bp, big.mark = ","),
               size = 1.8, color = "grey45", hjust = 1)

    if (has_empty_site && !is.null(ref_df) && nrow(ref_df)) {
      ref_x_range <- range(ref_df$xmin, ref_df$xmax)
      # Reference actual positions (before shift)
      # Recover actual reference positions: ref gene was shifted by -ref_gap_center
      # Recover actual reference positions: subtract ref_offset then add ref_gap_center
      ref_start_bp <- round(ref_x_range[1] - ref_offset + ref_gap_center)
      ref_end_bp <- round(ref_x_range[2] - ref_offset + ref_gap_center)
      ref_bb_df <- data.frame(x = ref_x_range[1], xend = ref_x_range[2], y = Y_REF, yend = Y_REF)
      p_synteny <- p_synteny +
        geom_segment(data = ref_bb_df, aes(x = x, xend = xend, y = y, yend = yend),
                     color = "grey60", linewidth = 0.4, inherit.aes = FALSE) +
        annotate("text", x = ref_x_range[1], y = Y_REF + 0.15,
                 label = format(ref_start_bp, big.mark = ","),
                 size = 1.8, color = "grey45", hjust = 0) +
        annotate("text", x = ref_x_range[2], y = Y_REF + 0.15,
                 label = format(ref_end_bp, big.mark = ","),
                 size = 1.8, color = "grey45", hjust = 1)
    }

    # Ribbons
    if (!is.null(ref_ribbon_data)) {
      if (!is.null(ref_ribbon_data$left))
        p_synteny <- p_synteny + geom_polygon(data = ref_ribbon_data$left,
          aes(x = x, y = y), fill = "#42A5F5", alpha = 0.20, color = NA)
      if (!is.null(ref_ribbon_data$right))
        p_synteny <- p_synteny + geom_polygon(data = ref_ribbon_data$right,
          aes(x = x, y = y), fill = "#66BB6A", alpha = 0.20, color = NA)
    }

    # Gene arrows
    is_df <- gene_df %>% filter(type == "IS_element")
    non_is_df <- gene_df %>% filter(type != "IS_element")

    if (nrow(non_is_df))
      p_synteny <- p_synteny + gggenes::geom_gene_arrow(
        data = non_is_df,
        aes(xmin = xmin, xmax = xmax, y = y, fill = type, forward = forward),
        arrowhead_height = unit(4.5, "mm"), arrowhead_width = unit(2.5, "mm"),
        arrow_body_height = unit(4, "mm"), color = "grey50", size = 0.25)
    if (nrow(is_df))
      p_synteny <- p_synteny + gggenes::geom_gene_arrow(
        data = is_df,
        aes(xmin = xmin, xmax = xmax, y = y, fill = type, forward = forward),
        arrowhead_height = unit(4.5, "mm"), arrowhead_width = unit(2.5, "mm"),
        arrow_body_height = unit(4, "mm"), color = "grey50", size = 0.25)

    # Gene labels — alternating y-offset to prevent overlap
    focal_labels <- gene_df %>%
      filter(nzchar(gene_label) & type != "IS_element" & row == "focal") %>%
      arrange(xmin) %>%
      mutate(mid = (xmin + xmax) / 2,
             label_y = y + ifelse(seq_len(dplyr::n()) %% 2 == 1, 0.18, 0.28))
    ref_labels <- gene_df %>%
      filter(nzchar(gene_label) & type != "IS_element" & row == "reference") %>%
      arrange(xmin) %>%
      mutate(mid = (xmin + xmax) / 2,
             label_y = y - ifelse(seq_len(dplyr::n()) %% 2 == 1, 0.18, 0.28))
    all_labels <- rbind(focal_labels, ref_labels)

    if (nrow(all_labels))
      p_synteny <- p_synteny + geom_text(data = all_labels,
        aes(x = mid, y = label_y, label = gene_label),
        inherit.aes = FALSE, size = 2.0, color = "grey25", fontface = "italic")

    # IS label
    if (nrow(is_df)) {
      is_mid <- (is_df$xmin[1] + is_df$xmax[1]) / 2
      is_label_text <- paste0(fam_name, " (", formatC(el_len, format = "d", big.mark = ","), " bp)")
      p_synteny <- p_synteny + annotate("text", x = is_mid, y = Y_FOCAL + 0.22,
        label = is_label_text, size = 3, fontface = "bold", color = "#B71C1C")
    }

    # Gap/identity annotation between rows
    if (nzchar(gap_label)) {
      p_synteny <- p_synteny + annotate("text", x = 0, y = (Y_FOCAL + Y_REF) / 2,
        label = gap_label, size = 2.5, color = "grey30", fontface = "italic", hjust = 0.5)
    }

    # Row labels — 90° rotated, centered vertically on each row
    x_lbl <- -flank_bp - flank_bp * 0.06  # Position in left margin

    # Reference label: organism + accession, rotated 90°
    if (has_empty_site) {
      ref_info <- ref_organism_info[[ref_file]]
      ref_organism <- if (!is.null(ref_info)) ref_info$organism else gsub("\\.gbff$", "", ref_file)
      ref_acc <- if (!is.null(ref_info) && nzchar(ref_info$accession)) ref_info$accession else ""
      # Ref: Genus\nspecies\nstrain+rest\naccession
      ref_parts <- strsplit(ref_organism, "\\s+")[[1]]
      if (length(ref_parts) >= 3) {
        ref_label <- paste0(ref_parts[1], "\n", ref_parts[2], "\n",
                            paste(ref_parts[3:length(ref_parts)], collapse = " "), "\n", ref_acc)
      } else if (length(ref_parts) == 2) {
        ref_label <- paste0(ref_parts[1], "\n", ref_parts[2], "\n", ref_acc)
      } else {
        ref_label <- paste0(ref_organism, "\n", ref_acc)
      }
      p_synteny <- p_synteny +
        annotate("text", x = x_lbl, y = Y_REF, label = ref_label, size = 2.0,
                 fontface = "italic", color = "grey40", angle = 90, vjust = 0.5, hjust = 0.5,
                 lineheight = 0.8)
    }

    # Focal genome label — rotated 90° in left margin
    focal_source <- tryCatch({
      fl <- readLines(genbank_path, n = 15, warn = FALSE)
      si <- grep("^SOURCE\\s+", fl)
      if (length(si)) trimws(sub("^SOURCE\\s+", "", fl[si[1]])) else contig_short
    }, error = function(e) contig_short)
    # No focal genome label — only reference comparison genome is labeled

    p_synteny <- p_synteny +
      scale_fill_manual(values = type_cols, labels = type_labels, drop = FALSE) +
      scale_x_continuous(
        expand = c(0, 0),
        labels = function(x) paste0(ifelse(x >= 0, "+", ""), round(x / 1000, 1), " kb"),
        breaks = scales::pretty_breaks(n = 5)) +
      coord_cartesian(xlim = c(-flank_bp, flank_bp),
                      ylim = if (has_empty_site) c(-0.3, 0.95) else c(0.15, 0.95), clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 7, color = "grey40"),
        axis.ticks.x = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.2),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_blank(), legend.position = "none",
        plot.margin = margin(4, 8, 4, 85)
      )
  }

  # Fallback if no synteny could be built
  if (is.null(p_synteny)) {
    p_synteny <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No synteny data",
               size = 3, color = "grey50", fontface = "italic") +
      xlim(0, 1) + ylim(0, 1) + theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA))
  }

  # ── Combine 3 columns into one row ──
  row_plot <- cowplot::plot_grid(
    p_info, p_logo, p_synteny,
    ncol = 3, rel_widths = c(0.10, 0.12, 0.78),
    align = "h", axis = "tb"
  )

  # Add thin bottom border for row separation
  row_plot <- cowplot::ggdraw(row_plot) +
    cowplot::draw_line(x = c(0.01, 0.99), y = c(0, 0),
                       color = "grey80", size = 0.4)

  unified_rows[[fam_name]] <- row_plot
}

# ── Build panel title ──
n_empty <- if (!is.null(iv_for_panel) && nrow(iv_for_panel))
  sum(iv_for_panel$status == "empty_site", na.rm = TRUE) else 0
n_total_verified <- if (!is.null(iv_for_panel)) nrow(iv_for_panel) else 0
panel_subtitle <- paste0(
  length(show_families), " families  |  ",
  sum(census$element_count), " elements across ", nrow(census), " families  |  ",
  n_empty, "/", n_total_verified, " empty-site verified"
)

# B title aligned with B label (both at top of panel)
p_panel_title <- cowplot::ggdraw() +
  cowplot::draw_label("IS Element Family Summary",
                      fontface = "bold", size = 14, x = 0.025, y = 0.92, hjust = 0, vjust = 1) +
  cowplot::draw_label(panel_subtitle, size = 9, x = 0.025, y = 0.55,
                      hjust = 0, vjust = 1, color = "grey45")

# ── Column headers ──
p_col_headers <- cowplot::ggdraw() +
  cowplot::draw_label("Family", fontface = "bold", size = 9, x = 0.05, hjust = 0.5,
                      color = "grey30") +
  cowplot::draw_label("TSD SeqLogo", fontface = "bold", size = 9, x = 0.16, hjust = 0.5,
                      color = "grey30") +
  cowplot::draw_label("Insertion Context (Focal vs Reference)", fontface = "bold",
                      size = 9, x = 0.61, hjust = 0.5, color = "grey30")

# ── Shared legend ──
legend_df <- data.frame(
  xmin = 0, xmax = 1, y = 1, forward = TRUE,
  type = factor(names(type_cols)[names(type_cols) != "empty_site"],
                levels = names(type_cols)[names(type_cols) != "empty_site"]))
p_legend <- ggplot(legend_df, aes(xmin = xmin, xmax = xmax, y = y,
                                   fill = type, forward = forward)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(4, "mm"),
                            arrow_body_height = unit(3.5, "mm")) +
  scale_fill_manual(values = type_cols, labels = type_labels, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 8.5),
        legend.title = element_text(size = 9.5, face = "bold"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_blank())
shared_legend <- cowplot::get_legend(p_legend)

# Confidence legend
p_conf_legend <- ggplot(data.frame(lev = factor(c("High", "Medium", "Low"),
                                                 levels = c("Low", "Medium", "High")),
                                    cnt = 1),
                         aes(x = lev, y = cnt, fill = lev)) +
  geom_col() +
  scale_fill_manual(values = c(High = "#66BB6A", Medium = "#42A5F5", Low = "#FFC107"),
                    name = "Confidence") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))
conf_legend <- cowplot::get_legend(p_conf_legend)

combined_legends <- cowplot::plot_grid(conf_legend, shared_legend, nrow = 1,
                                       rel_widths = c(0.35, 0.65))

# ── Stack all rows ──
n_rows <- length(unified_rows)
if (n_rows > 0) {
  row_grid <- cowplot::plot_grid(
    plotlist = unified_rows, ncol = 1,
    align = "v", axis = "lr",
    rel_heights = rep(1, n_rows))

  gg_combined <- cowplot::plot_grid(
    p_panel_title,
    p_col_headers,
    row_grid,
    combined_legends,
    ncol = 1,
    rel_heights = c(0.05, 0.03, 0.87, 0.05),
    labels = c("B", "", "", ""),
    label_size = 22, label_fontface = "bold",
    label_x = 0, label_y = 1, hjust = 0, vjust = 1
  )
} else {
  gg_combined <- cowplot::ggdraw() +
    cowplot::draw_label("No IS family data available", size = 12, color = "grey50")
}

# ======================================================================
# Full-vector PDF: circos (Panel A) + ggplot (Panel B) via gridBase
# No raster intermediaries — everything stays as vector graphics.
# ======================================================================
bottom_height <- max(6, n_rows * 1.3 + 2)
top_height <- 14
total_height <- top_height + bottom_height

# Fraction of page for Panel A (top) and Panel B (bottom)
frac_top <- top_height / total_height
frac_bot <- bottom_height / total_height

cairo_pdf(out_file, width = 24, height = total_height, bg = "white")

# --- Panel A: Circos plot (base R graphics in top portion) ---
par(fig = c(0, 1, 1 - frac_top, 1), new = FALSE, mar = c(0.5, 0.2, 2.5, 0.2))

circos.clear()
circos.par(
  start.degree      = 90,
  gap.degree         = 3,
  cell.padding       = c(0, 0, 0, 0),
  track.margin       = c(0.015, 0.015),
  points.overflow.warning = FALSE
)

# Initialize — all contigs (chromosome + plasmid)
contig_lens <- genes %>% group_by(contig) %>% summarise(len = max(end, na.rm=TRUE), .groups="drop") %>% arrange(desc(len))
genome_bed <- data.frame(chr = contig_lens$contig, start = 0, end = contig_lens$len)
genome_len <- max(contig_lens$len)
circos.genomicInitialize(genome_bed, plotType = NULL,
                         major.by = 200000, axis.labels.cex = 0.9)

# Title — offset right to leave room for A label in cowplot
total_len_mb <- round(sum(contig_lens$len)/1e6, 2)
genome_accessions <- paste(contig_lens$contig, collapse = ", ")
mtext("IS Element Genome Distribution & Landing Pad Map", side = 3, line = 1.3,
      adj = 0.03, cex = 1.2, font = 2)
mtext(paste0(genome_accessions, "  |  ",
             nrow(contig_lens), " replicon(s), ", total_len_mb, " Mb  |  ",
             nrow(is_track), " IS elements, ",
             length(unique(is_track$fam)), " families"),
      side = 3, line = 0.2, adj = 0.03, cex = 0.85, col = "grey40")

# Track style constants
TH <- 0.07      # uniform track height
BD <- "grey50"   # border color
BW <- 0.5        # border width

# --- Track 0: Genomic coordinate axis (kb format) + sector accession label ---
circos.track(ylim = c(0, 1), track.height = 0.035, bg.border = NA,
             panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  sector_len <- contig_lens$len[contig_lens$contig == sector]
  # Accession label inside the outer rim
  circos.text(CELL_META$xcenter, 0.5, sector,
              facing = "bending.inside", niceFacing = TRUE,
              cex = 0.6, font = 2, col = "grey30")
  major <- if (sector_len > 500000) 200000 else 10000
  tick_at <- seq(0, sector_len, by = major)
  tick_labels <- paste0(format(tick_at / 1000, big.mark = ",", trim = TRUE), " kb")
  circos.genomicAxis(h = "top", major.at = tick_at, labels = tick_labels,
                     labels.cex = 0.85, labels.facing = "clockwise")
})

# --- Track 1: Gene essentiality heatmap ---
ess_col_fun <- colorRamp2(c(0, 0.15, 0.35, 0.55, 0.75, 1.0),
                           c("#E3F2FD","#90CAF9","#42A5F5",
                             "#1E88E5","#1565C0","#0D47A1"))
circos.track(ylim = c(0, 1), track.height = TH, bg.border = BD, bg.lwd = BW,
             panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  d <- ess_track[ess_track$contig == sector, , drop = FALSE]
  for (i in seq_len(nrow(d))) {
    circos.rect(d$start[i], 0, d$end[i], 1,
                col = ess_col_fun(d$essentiality_score[i]), border = NA)
  }
})

# --- Track 2: IS elements (family-colored) ---
circos.track(ylim = c(0, 1), track.height = TH, bg.border = BD, bg.lwd = BW,
             panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  d <- is_track[is_track$contig == sector, , drop = FALSE]
  for (i in seq_len(nrow(d))) {
    col <- fam_cols[d$fam[i]]
    if (is.na(col)) col <- "#BDBDBD"
    circos.rect(d$start[i], 0, d$end[i], 1,
                col = col, border = adjustcolor("grey30", 0.5), lwd = 0.25)
  }
})

# --- Track 3: Landing pad score (windowed max, score as height) ---
lp_col_fun <- colorRamp2(c(0, 0.20, 0.40, 0.60, 0.80, 1.0),
                          c("#E53935","#FF8F00","#FDD835","#66BB6A","#2E7D32","#1B5E20"))
circos.track(ylim = c(0, 1), track.height = TH, bg.border = BD, bg.lwd = BW,
             panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  d <- lp_track[lp_track$contig == sector, , drop = FALSE]
  sector_len <- contig_lens$len[contig_lens$contig == sector]
  win <- if (sector_len > 500000) 10000L else 2000L
  nw <- max(1L, ceiling(sector_len / win))
  for (w in seq_len(nw)) {
    ws <- (w-1)*win; we <- min(w*win, sector_len)
    lps_in <- d[d$region_start <= we & d$region_end >= ws, , drop = FALSE]
    if (nrow(lps_in)) {
      best_sc <- max(lps_in$landing_pad_score, na.rm = TRUE)
      # Height = score directly (0-1 range)
      h <- pmin(1, best_sc)
      circos.rect(ws, 0, we, h, col = lp_col_fun(best_sc), border = NA)
    }
  }
  # Reference line at "high" confidence threshold (0.60)
  circos.lines(CELL_META$cell.xlim, c(0.60, 0.60), lty = 3, col = "grey30", lwd = 0.4)
})

# --- Track 4: IS Target Site Density (z-score normalized, deviation from mean) ---
target_sites <- result$target_sites
if (!is.null(target_sites) && nrow(target_sites)) {
  circos.track(ylim = c(-1, 1), track.height = TH, bg.border = BD, bg.lwd = BW,
               panel.fun = function(x, y) {
    sector <- CELL_META$sector.index
    d <- target_sites[target_sites$contig == sector, , drop = FALSE]
    sector_len <- contig_lens$len[contig_lens$contig == sector]
    win <- if (sector_len > 500000) 50000L else 5000L
    nw <- max(1L, ceiling(sector_len / win))
    counts <- vapply(seq_len(nw), function(w) {
      ws <- (w-1)*win; we <- w*win
      sum(d$site_start <= we & d$site_end >= ws)
    }, integer(1))
    # Z-score normalize: highlight deviation from mean
    mu <- mean(counts); sigma <- max(sd(counts), 1)
    z <- pmin(1, pmax(-1, (counts - mu) / (2 * sigma)))
    mids <- seq(win/2, by = win, length.out = nw)
    # Above mean = red (high density = risky), below mean = blue (low density = safer)
    above <- pmax(0, z); below <- pmin(0, z)
    circos.lines(mids, above, col = adjustcolor("#E53935", 0.45), area = TRUE, border = NA, baseline = 0)
    circos.lines(mids, below, col = adjustcolor("#1565C0", 0.35), area = TRUE, border = NA, baseline = 0)
    circos.lines(CELL_META$cell.xlim, c(0, 0), col = "grey50", lwd = 0.3)
  })
}

# --- Track 5: Coding Density by strand (+strand above, -strand below) ---
circos.track(ylim = c(-1, 1), track.height = TH, bg.border = BD, bg.lwd = BW,
             panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  d <- ess_track[ess_track$contig == sector, , drop = FALSE]
  d_plus <- d[d$strand == "+", , drop = FALSE]
  d_minus <- d[d$strand == "-", , drop = FALSE]
  sector_len <- contig_lens$len[contig_lens$contig == sector]
  win <- if (sector_len > 500000) 20000L else 3000L
  nw <- max(1L, ceiling(sector_len / win))
  mids <- seq(win/2, by = win, length.out = nw)

  # + strand density (above baseline)
  plus_dens <- vapply(seq_len(nw), function(w) {
    ws <- (w-1)*win; we <- min(w*win, sector_len)
    sum(pmax(0, pmin(d_plus$end, we) - pmax(d_plus$start, ws))) / win
  }, numeric(1))
  # - strand density (below baseline)
  minus_dens <- vapply(seq_len(nw), function(w) {
    ws <- (w-1)*win; we <- min(w*win, sector_len)
    sum(pmax(0, pmin(d_minus$end, we) - pmax(d_minus$start, ws))) / win
  }, numeric(1))

  circos.lines(mids, plus_dens, col = adjustcolor("#43A047", 0.45), area = TRUE, border = NA, baseline = 0)
  circos.lines(mids, -minus_dens, col = adjustcolor("#FF8F00", 0.45), area = TRUE, border = NA, baseline = 0)
  circos.lines(CELL_META$cell.xlim, c(0, 0), col = "grey50", lwd = 0.3)
})

# ============================================================
# IS family links
# ============================================================
link_families <- is_track %>%
  count(fam, sort = TRUE) %>%
  filter(n >= 2) %>%
  pull(fam)

for (fam_name in link_families) {
  for (ctg_name in all_contigs) {
    fam_elements <- is_track %>% filter(fam == fam_name, contig == ctg_name) %>% arrange(start)
    if (nrow(fam_elements) < 2) next

    col <- fam_cols[fam_name]
    if (is.na(col)) col <- "#BDBDBD"
    link_col   <- adjustcolor(col, alpha.f = 0.18)
    border_col <- adjustcolor(col, alpha.f = 0.40)

    for (i in seq_len(nrow(fam_elements) - 1)) {
      s1 <- fam_elements$start[i]; e1 <- fam_elements$end[i]
      s2 <- fam_elements$start[i+1]; e2 <- fam_elements$end[i+1]
      dist <- s2 - e1

      if (dist > 30000) {
        h <- ifelse(dist > 500000, 0.5, ifelse(dist > 200000, 0.4, 0.3))
        circos.link(ctg_name, c(s1, e1), ctg_name, c(s2, e2),
                    col = link_col, border = border_col,
                    lwd = 0.4, h.ratio = h)
      }
    }
  }
}

# ============================================================
# Legends — drawn in the circos panel margins
# ============================================================
# Prepare legend objects
fam_counts <- sort(table(is_track$fam), decreasing = TRUE)
legend_fams <- names(fam_counts)
# Assign default colors to families not in fam_cols
extra_colors <- c("#795548","#607D8B","#FF5722","#9C27B0","#00BCD4","#CDDC39","#F44336","#3F51B5")
for (i in seq_along(legend_fams)) {
  if (is.na(fam_cols[legend_fams[i]])) {
    fam_cols[legend_fams[i]] <- extra_colors[((i - 1) %% length(extra_colors)) + 1]
  }
}
legend_labels <- paste0(legend_fams, " (", fam_counts[legend_fams], ")")

lgd_family <- Legend(
  labels = legend_labels,
  legend_gp = gpar(fill = fam_cols[legend_fams], col = NA),
  title = "IS element families",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(3.5, "mm"),
  grid_width = unit(3.5, "mm"),
  gap = unit(0.5, "mm"),
  ncol = 2
)

lgd_ess <- Legend(
  col_fun = ess_col_fun, at = c(0, 0.25, 0.5, 0.75, 1),
  title = "Gene essentiality",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  legend_width = unit(3.5, "cm"),
  direction = "horizontal"
)

lgd_lp <- Legend(
  col_fun = lp_col_fun, at = seq(0, 1.0, by = 0.25),
  title = "Landing pad score",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  legend_width = unit(3.5, "cm"),
  direction = "horizontal"
)

lgd_tracks <- Legend(
  labels = c("1. Gene essentiality",
             "2. IS elements (family)",
             "3. Landing pad score",
             "4. IS target density",
             "5. Coding (+/- strand)"),
  type = "lines",
  legend_gp = gpar(lwd = 3, col = c("#1565C0","#7E57C2","#43A047","#C62828","#FF8F00")),
  title = "Tracks (outer to inner)",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(3, "mm"),
  grid_width = unit(5, "mm")
)

# LEFT bottom of A panel: Tracks (bottom) + IS Element Families (above)
pushViewport(viewport(x = 0.02, y = 0.42, just = c("left", "bottom"),
                       width = 0.28, height = 0.22))
pushViewport(viewport(x = 0, y = 0, just = c("left", "bottom"),
                       width = 1, height = 0.30))
draw(lgd_tracks, just = c("left", "bottom"))
popViewport()
pushViewport(viewport(x = 0, y = 0.30, just = c("left", "bottom"),
                       width = 1, height = 0.70))
draw(lgd_family, just = c("left", "bottom"))
popViewport()
popViewport()

# RIGHT bottom of A panel: Landing Pad Score (bottom) + Gene Essentiality (above)
pushViewport(viewport(x = 0.90, y = 0.42, just = c("right", "bottom"),
                       width = 0.28, height = 0.14))
pushViewport(viewport(x = 0, y = 0, just = c("left", "bottom"),
                       width = 1, height = 0.48))
draw(lgd_lp, just = c("left", "bottom"))
popViewport()
pushViewport(viewport(x = 0, y = 0.48, just = c("left", "bottom"),
                       width = 1, height = 0.48))
draw(lgd_ess, just = c("left", "bottom"))
popViewport()
popViewport()

# (Footer note removed per user request)

# "A" label — top-left corner of Panel A
grid.text("A", x = unit(0.3, "lines"), y = unit(1, "npc") - unit(0.2, "lines"),
          just = c("left", "top"), gp = gpar(fontsize = 22, fontface = "bold"))

circos.clear()

# --- Panel B: ggplot (grid graphics in bottom portion) ---
# Use par(fig=...) to define a base-R plot region in the bottom portion,
# then use gridBase to switch to grid context and print the ggplot object.
par(fig = c(0, 1, 0, frac_bot), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()

# gridBase: capture the current base viewport and push into grid
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
print(gg_combined, vp = viewport(x = 0.5, y = 0.5, width = 1, height = 1))
popViewport(3)

dev.off()

cat("Saved:", out_file, "\n")
cat("Size:", round(file.size(out_file)/1024, 1), "KB\n")
