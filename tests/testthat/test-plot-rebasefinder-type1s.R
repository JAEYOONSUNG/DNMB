test_that("REBASE plot matching preserves numeric locus suffixes", {
  query <- c("ACHFCC_RS00020", "ACHFCC_RS05295")
  source <- c("ACHFCC_RS00020", "ACHFCC_RS05295")

  expect_identical(
    DNMB:::.dnmb_rebasefinder_match_plot_ids(query, source),
    c(1L, 2L)
  )

  # Cleaned keys are ambiguous here, so an inexact fallback must not attach
  # evidence to the wrong suffix sibling.
  expect_true(is.na(DNMB:::.dnmb_rebasefinder_match_plot_ids(
    "ACHFCC_RS99999", source
  )))
})

test_that("Type I-S plotting sidecar joins by exact locus without touching 3D fields", {
  root <- tempfile("dnmb-type1s-overlay-")
  module_dir <- file.path(root, "dnmb_module_rebasefinder")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  predictions <- data.frame(
    locus_tag = "ACHFCC_RS05295",
    type1s_prediction_status = "predicted_complete",
    type1s_predicted_recognition = "GGCANNNNNNTTC",
    type1s_overall_confidence = "low",
    stringsAsFactors = FALSE
  )
  write.table(
    predictions,
    file.path(module_dir, "DNMB_REBASEfinder_type1s_predictions.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  input <- data.frame(
    locus_tag = c("ACHFCC_RS00020", "ACHFCC_RS05295"),
    REBASEfinder_homology_geometry_status = c("geometry_a", "geometry_b"),
    stringsAsFactors = FALSE
  )
  observed <- DNMB:::.dnmb_rebasefinder_overlay_type1s_results(input, root)

  expect_true(is.na(observed$REBASEfinder_type1s_predicted_recognition[[1]]))
  expect_identical(
    observed$REBASEfinder_type1s_predicted_recognition[[2]],
    "GGCANNNNNNTTC"
  )
  expect_identical(
    observed$REBASEfinder_homology_geometry_status,
    c("geometry_a", "geometry_b")
  )
})

.type1s_plot_fixture <- function() {
  data.frame(
    locus_tag = "ACHFCC_RS05295",
    contig = "chromosome",
    start = 100,
    end = 1362,
    translation = paste(rep("A", 421), collapse = ""),
    REBASEfinder_family_id = "Type I",
    REBASEfinder_enzyme_role = "S",
    REBASEfinder_hit_label = "S.Saq4ORF5310P",
    REBASEfinder_blast_identity = 1,
    REBASEfinder_blast_length = 421,
    REBASEfinder_rec_seq = NA_character_,
    REBASEfinder_reference_rec_seq = "GGCANNNNNNTGG",
    REBASEfinder_typing_eligible = TRUE,
    REBASEfinder_type1s_prediction_status = "predicted_complete",
    REBASEfinder_type1s_predicted_recognition = "GGCANNNNNNTTC",
    REBASEfinder_type1s_spacer_length = 6L,
    REBASEfinder_type1s_overall_confidence = "low",
    REBASEfinder_type1s_spacer_method = "whole_hsds_scaffold_ensemble",
    REBASEfinder_type1s_spacer_model_agreement = TRUE,
    REBASEfinder_type1s_spacer_vote_support = 0.8527497,
    REBASEfinder_type1s_spacer_scaffold_vote_support = 0.5091657,
    REBASEfinder_type1s_spacer_tael_ruler_applied = FALSE,
    REBASEfinder_structure_coverage_status = "structure_missing",
    `Signature accession_Pfam` = "PF01420;PF01420",
    `Signature description_CDD` = "TRD1;TRD2",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("TRD prediction is a dedicated line and motif 3D stays separate", {
  input <- .type1s_plot_fixture()
  annotation <- DNMB:::.dnmb_rebasefinder_type1s_display_annotations(input)

  expect_identical(
    annotation$prediction_line,
    "TRD prediction: GGCA-N6-TTC [low]"
  )
  expect_identical(
    annotation$evidence_line,
    "TRD: 2×TRD · N6 match · E85/S51"
  )
  expect_identical(annotation$context_evidence_line, annotation$evidence_line)
  expect_false(grepl("3D", annotation$evidence_line, fixed = TRUE))

  display <- DNMB:::.dnmb_rebasefinder_display_labels(input)
  lines <- strsplit(as.character(display$display), "\n", fixed = TRUE)[[1]]
  expect_identical(
    lines[[1]],
    "ACHFCC_RS05295 | ref rec: GGCANNNNNNTGG"
  )
  expect_match(lines[[2]], "S.Saq4ORF5310P", fixed = TRUE)
  expect_identical(lines[[3]], annotation$prediction_line)
  expect_identical(lines[[4]], annotation$evidence_line)

  blast <- DNMB:::.dnmb_plot_rebasefinder_blast_quality(
    input,
    DNMB:::.dnmb_rebasefinder_role_palette(input),
    display,
    DNMB:::.dnmb_rebasefinder_palette(input$REBASEfinder_family_id)
  )
  typed_lines <- strsplit(
    as.character(blast$data$display_typed[[1]]), "<br>|\n", perl = TRUE
  )[[1]]
  expect_match(typed_lines[[1]], "^ACHFCC_RS05295 \\| ref rec:")
  expect_match(typed_lines[[2]], "[Type I]", fixed = TRUE)
})

test_that("REBASE overview expands for multiline Type I-S labels", {
  legacy <- DNMB:::.dnmb_rebasefinder_overview_widths(9L)
  expect_equal(unname(legacy[c("blast", "domain", "motif")]), c(6.9, 2.35, 4.87))

  display <- data.frame(
    display = paste(
      "ACHFCC_RS05295 | ref rec: GGCANNNNNNTGG",
      "A deliberately long REBASE label that needs additional horizontal room",
      "TRD prediction: GGCA-N6-TTC [low]",
      "TRD: 2×TRD · N6 match · E85/S51",
      sep = "\n"
    ),
    stringsAsFactors = FALSE
  )
  metrics <- DNMB:::.dnmb_rebasefinder_display_line_metrics(display)
  expanded <- DNMB:::.dnmb_rebasefinder_overview_widths(
    9L, metrics$max_line_width_in
  )

  expect_identical(metrics$max_lines, 4L)
  expect_identical(metrics$extra_lines, 2L)
  expect_gt(expanded[["blast"]], legacy[["blast"]])
  expect_lte(expanded[["blast"]], 9.5)
})
