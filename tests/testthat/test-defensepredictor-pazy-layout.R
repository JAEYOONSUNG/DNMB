test_that("DefensePredictor visual tiers do not change the stringent call cutoff", {
  scores <- c(-1, 0, 1.99, 2, 3.99, 4, 6, 8, 10)
  bands <- as.character(DNMB:::.dnmb_defensepredictor_score_band(scores))

  expect_identical(
    bands,
    c(
      "DP_below0", "DP_0to2", "DP_0to2", "DP_2to4", "DP_2to4",
      "DP_4to6", "DP_6to8", "DP_8to10", "DP_10plus"
    )
  )
  expect_equal(DNMB:::.dnmb_defensepredictor_default_threshold(), 4)
})

test_that("DefensePredictor labels remain inside each replicon", {
  candidates <- data.frame(
    contig_id = rep("DNMB_CONTIG_002", 4),
    midpoint = c(20, 48, 51, 98),
    locus_tag = paste0("AC5Q2X_RS", c("15115", "15260", "15285", "15290")),
    gene = c("", "", "", ""),
    DefensePredictor_mean_log_odds = c(1.32, 0.51, 0.99, 0.48),
    stringsAsFactors = FALSE
  )
  contigs <- data.frame(
    contig_id = "DNMB_CONTIG_002",
    length_bp = 100,
    stringsAsFactors = FALSE
  )

  packed <- DNMB:::.dnmb_defensepredictor_pack_labels(
    candidates, contigs, panel_width_in = 7.8
  )

  expect_true(all(packed$label_x - packed$label_half_width >= 0))
  expect_true(all(packed$label_x + packed$label_half_width <= 100))
  expect_true(all(nzchar(packed$feature_label)))
})

test_that("DefensePredictor input avoids the upstream pseudo substring trap", {
  out_dir <- tempfile("defensepredictor-input-")
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  genes <- data.frame(
    locus_tag = "LOCUS123",
    protein_id = "WP_000001.1",
    translation = "MNNNNNNNNN",
    nt_seq = "ATGAACAACAACAACAACAACAACAACAAC",
    contig = "Genome with spaces",
    contig_number = 1,
    start = 1,
    end = 30,
    direction = "+",
    product = "pseudouridine synthase",
    gene = "truB",
    stringsAsFactors = FALSE
  )

  written <- DNMB:::.dnmb_write_defensepredictor_input(
    genes, output_dir = out_dir, assembly_id = "fixture"
  )
  feature_tbl <- utils::read.delim(
    written$feature_table, stringsAsFactors = FALSE, check.names = FALSE
  )

  expect_identical(feature_tbl$name, "pseudouridine synthase")
  expect_false(grepl("pseudo", feature_tbl$attributes, ignore.case = TRUE))
  expect_match(feature_tbl$attributes, "protein_id=WP_000001.1")
  expect_match(feature_tbl$attributes, "locus_tag=LOCUS123")
})

test_that("DefensePredictor cache metadata tracks input schema changes", {
  out_dir <- tempfile("defensepredictor-cache-")
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("score", file.path(out_dir, "defense_predictor_output.csv"))
  input_path <- file.path(out_dir, "fixture.gbff")
  writeLines("LOCUS fixture", input_path)
  signature <- DNMB:::.dnmb_file_signature(input_path)

  DNMB:::.dnmb_write_defensepredictor_signature(
    out_dir, genbank_signature = signature, version = "current"
  )
  current <- DNMB:::.dnmb_defensepredictor_reuse_status(
    out_dir, genbank_signature = signature, version = "current"
  )
  expect_true(current$reusable)

  legacy <- readRDS(DNMB:::.dnmb_defensepredictor_metadata_path(out_dir))
  legacy$input_schema_version <- NULL
  saveRDS(legacy, DNMB:::.dnmb_defensepredictor_metadata_path(out_dir))
  stale <- DNMB:::.dnmb_defensepredictor_reuse_status(
    out_dir, genbank_signature = signature, version = "current"
  )
  expect_false(stale$reusable)
  expect_match(stale$reason, "input_schema_changed")
})

test_that("PAZy hit annotations pack without truncation or edge clipping", {
  hits <- data.frame(
    start = c(2528575, 2589130, 2598147),
    end = c(2529453, 2590008, 2599088),
    pazy_label = c("214 (23%)", "Est8_89 (25%)", "PlaM9 (46%)"),
    product = c(
      "alpha/beta hydrolase family protein",
      "proline iminopeptidase-family hydrolase",
      "alpha/beta hydrolase"
    ),
    stringsAsFactors = FALSE
  )
  loc_start <- 2524000
  loc_end <- 2604000
  packed <- DNMB:::.dnmb_pazy_pack_hit_annotations(
    hits, loc_start = loc_start, loc_end = loc_end
  )

  expect_true(all(packed$annotation_x >= loc_start))
  expect_true(all(packed$annotation_x <= loc_end))
  expect_false(any(grepl("\\.\\.\\.", packed$product_label)))
  expect_true(all(grepl("<i>", packed$combined_label, fixed = TRUE)))
  expect_true(all(grepl("<b>", packed$combined_label, fixed = TRUE)))

  hit_chars <- DNMB:::.dnmb_pazy_max_line_chars(packed$pazy_label)
  product_chars <- DNMB:::.dnmb_pazy_max_line_chars(packed$product_label)
  width_in <- pmax(
    hit_chars * 6.55 * 0.52 / 72,
    product_chars * 5.12 * 0.52 / 72
  )
  span <- loc_end - loc_start
  half_width <- pmin(0.34, width_in / 3.65 / 2 + 0.012) * span
  expect_true(all(packed$annotation_x - half_width >= loc_start - 1e-7))
  expect_true(all(packed$annotation_x + half_width <= loc_end + 1e-7))

  for (tier in unique(packed$annotation_tier)) {
    ii <- which(packed$annotation_tier == tier)
    if (length(ii) > 1L) {
      ord <- ii[order(packed$annotation_x[ii])]
      expect_true(all(
        packed$annotation_x[ord[-1L]] - half_width[ord[-1L]] >=
          packed$annotation_x[ord[-length(ord)]] + half_width[ord[-length(ord)]]
      ))
    }
  }
})

test_that("PAZy rich labels escape markup-sensitive product text", {
  hits <- data.frame(
    start = 10,
    end = 20,
    pazy_label = "PAZY&1 (50%)",
    product = "A < B & C > D",
    stringsAsFactors = FALSE
  )
  packed <- DNMB:::.dnmb_pazy_pack_hit_annotations(hits, 0, 100)

  expect_match(packed$combined_label, "A &lt; B &amp; C &gt; D", fixed = TRUE)
  expect_match(packed$combined_label, "PAZY&amp;1", fixed = TRUE)
})

test_that("PAZy locus tags are packed below gene arrows", {
  genes <- data.frame(
    start = c(10, 11, 12),
    end = c(20, 21, 22),
    gene = c("a", "b", "c"),
    gene_label = c("LOCUS_0001", "LOCUS_0002", "LOCUS_0003"),
    is_pazy = FALSE,
    stringsAsFactors = FALSE
  )

  packed <- DNMB:::.dnmb_pazy_pack_context_labels(
    genes, loc_start = 0, loc_end = 100, max_tiers = 3L
  )

  expect_equal(nrow(packed), 3L)
  expect_true(all(packed$context_label_y < 1))
  expect_true(all(packed$context_label_y > 0.42))
  expect_equal(
    packed$context_label_y,
    0.83 - 0.20 * packed$context_tier,
    tolerance = 1e-8
  )
})

test_that("PAZy locus view includes complete edge genes and clamps to genome", {
  genes <- data.frame(
    start = c(80, 190),
    end = c(120, 240),
    stringsAsFactors = FALSE
  )
  bounds <- DNMB:::.dnmb_pazy_locus_view_bounds(
    genes, loc_start = 100, loc_end = 200,
    genome_len = 230, min_pad_bp = 10
  )

  expect_lte(bounds[["start"]], 70)
  expect_equal(bounds[["end"]], 230)
  expect_lte(bounds[["start"]], min(genes$start))
  expect_gte(bounds[["end"]], min(max(genes$end), 230))
})
