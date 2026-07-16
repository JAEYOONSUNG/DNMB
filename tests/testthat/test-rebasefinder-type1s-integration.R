test_that("Type I-S sidecars atomically replace stale data with typed zero rows", {
  output <- tempfile("rebasefinder-type1s-sidecar-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE, force = TRUE), add = TRUE)

  stale_tsv <- file.path(output, "DNMB_REBASEfinder_type1s_predictions.tsv")
  stale_xlsx <- file.path(output, "DNMB_REBASEfinder_type1s_predictions.xlsx")
  writeLines("stale", stale_tsv)
  writeLines("stale", stale_xlsx)

  prediction <- DNMB:::.dnmb_type1s_empty_prediction("hsds_1")
  prediction$type1s_prediction_status <- "predicted_complete"
  prediction$type1s_predicted_recognition <- "CAANNNNNNCTC"
  prediction$type1s_overall_confidence <- "low"
  files <- DNMB:::.dnmb_rebasefinder_write_type1s_predictions(
    output,
    prediction
  )

  expect_true(all(file.exists(unlist(files))))
  observed <- read.delim(files$type1s_predictions_tsv, check.names = FALSE)
  expect_identical(observed$locus_tag, "hsds_1")
  expect_identical(
    observed$type1s_predicted_recognition,
    "CAANNNNNNCTC"
  )
  expect_identical(
    openxlsx::getSheetNames(files$type1s_predictions_xlsx),
    "Type_I_S_predictions"
  )

  empty <- DNMB:::.dnmb_type1s_empty_prediction(character())
  files <- DNMB:::.dnmb_rebasefinder_write_type1s_predictions(output, empty)
  cleared <- read.delim(files$type1s_predictions_tsv, check.names = FALSE)
  expect_equal(nrow(cleared), 0L)
  expect_identical(names(cleared), names(empty))
  expect_equal(
    nrow(openxlsx::read.xlsx(
      files$type1s_predictions_xlsx,
      sheet = "Type_I_S_predictions"
    )),
    0L
  )
  expect_false(any(grepl(
    "[.]DNMB_REBASEfinder_type1s_predictions-",
    list.files(output, all.files = TRUE)
  )))
})

test_that("Type I-S sidecar pair rolls back when the second install fails", {
  output <- tempfile("rebasefinder-type1s-sidecar-rollback-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE, force = TRUE), add = TRUE)

  staged <- file.path(output, c("new.tsv", "new.xlsx"))
  destination <- file.path(output, c("current.tsv", "current.xlsx"))
  writeLines(c("new tsv", "new row"), staged[[1L]])
  writeLines(c("new xlsx placeholder", "new row"), staged[[2L]])
  writeLines("old tsv", destination[[1L]])
  writeLines("old xlsx", destination[[2L]])

  rename_with_second_install_failure <- function(from, to) {
    if (identical(from, staged[[2L]]) && identical(to, destination[[2L]])) {
      return(FALSE)
    }
    file.rename(from, to)
  }

  expect_error(
    DNMB:::.dnmb_rebasefinder_atomic_replace_pair(
      staged,
      destination,
      rename_file = rename_with_second_install_failure
    ),
    "Could not install REBASEfinder sidecar.*current[.]xlsx"
  )
  expect_identical(readLines(destination[[1L]]), "old tsv")
  expect_identical(readLines(destination[[2L]]), "old xlsx")
  expect_false(any(grepl("-previous-", list.files(output, all.files = TRUE))))
})

test_that("Type I-S metadata attaches only to existing Type I HsdS hits", {
  hits <- data.frame(
    query = c("hsds_1", "hsds_2", "type2_s", "hsds_3", "hsds_high"),
    family_id = c("Type I", "Type I", "Type II", "Type I", "Type I"),
    enzyme_role = c("S", "R", "S", "S", "S"),
    evidence_mode = c("operon_context", "high_confidence", "direct", "review", "review"),
    substrate_label = c(NA, "existing-r", "existing-type2", NA, NA),
    rec_seq = c(NA, "AACNNNNNNGTT", "CCANNNNNNTGG", NA, NA),
    typing_eligible = c(FALSE, TRUE, TRUE, FALSE, FALSE),
    curation_tier = c("medium", "high", "high", "review", "review"),
    curation_keep = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  core_columns <- names(hits)
  core_before <- hits

  predictions <- DNMB:::.dnmb_type1s_empty_prediction(
    c("hsds_1", "hsds_2", "type2_s", "hsds_high")
  )
  predictions$type1s_prediction_status <- "predicted_complete"
  predictions$type1s_predicted_recognition <- c(
    "CAANNNNNNCTC", "GGCANNNNNNTTC", "ACGNNNNNNGTG", "GGANNNNNNTCC"
  )
  predictions$type1s_overall_confidence <- c("low", "high", "high", "high")
  predictions$type1s_prediction_eligible <- c(FALSE, TRUE, TRUE, TRUE)

  observed <- DNMB:::.dnmb_rebasefinder_attach_type1s_metadata(
    hits,
    predictions
  )

  expect_equal(nrow(observed), nrow(hits))
  expect_identical(observed[, core_columns, drop = FALSE], core_before)
  expect_identical(
    observed$type1s_predicted_recognition[[1]],
    "CAANNNNNNCTC"
  )
  expect_true(is.na(observed$type1s_predicted_recognition[[2]]))
  expect_true(is.na(observed$type1s_predicted_recognition[[3]]))
  # hsds_3 shares the generic cleaned key with two prediction IDs, so the
  # fallback is ambiguous and must abstain rather than leak a suffix sibling.
  expect_true(is.na(observed$type1s_predicted_recognition[[4]]))
  expect_identical(observed$type1s_predicted_recognition[[5]], "GGANNNNNNTCC")
  expect_false(observed$typing_eligible[[1]])
  # Even a prediction whose model-specific flag is high cannot bypass the
  # canonical review tier or promote the module's typing eligibility.
  expect_true(observed$type1s_prediction_eligible[[5]])
  expect_false(observed$typing_eligible[[5]])
  expect_false(observed$curation_keep[[5]])
  expect_identical(observed$curation_tier, core_before$curation_tier)
})

test_that("REBASEfinder Type I-S runner selects HsdS and reports prediction status", {
  genes <- data.frame(
    locus_tag = c("hsds_1", "hsdm_1"),
    product = c(
      "type I restriction-modification system specificity subunit",
      "type I restriction-modification system methylase subunit"
    ),
    translation = c(
      paste(rep("A", 421), collapse = ""),
      paste(rep("M", 500), collapse = "")
    ),
    stringsAsFactors = FALSE
  )
  rm <- data.frame(
    locus_tag = c("hsds_1", "hsdm_1"),
    rm_type = c("Type I", "Type I"),
    subunit = c("S", "M"),
    stringsAsFactors = FALSE
  )
  captured <- NULL
  local_mocked_bindings(
    dnmb_predict_type1s_recognition = function(candidates, ...) {
      captured <<- list(candidates = candidates, args = list(...))
      out <- DNMB:::.dnmb_type1s_empty_prediction(candidates$locus_tag)
      out$type1s_prediction_status <- "predicted_complete"
      out$type1s_predicted_recognition <- "CAANNNNNNCTC"
      out$type1s_overall_confidence <- "low"
      out
    },
    .package = "DNMB"
  )

  result <- DNMB:::.dnmb_rebasefinder_predict_type1s(
    genes,
    rm,
    cache_dir = tempfile("type1s-cache-"),
    cpu = 3L,
    verbose = FALSE
  )

  expect_identical(result$predictions$locus_tag, "hsds_1")
  expect_identical(result$status$component, "type1s_prediction")
  expect_identical(result$status$status, "ok")
  expect_match(result$status$detail, "1 HsdS candidates; 1 complete")
  expect_identical(captured$candidates$locus_tag, "hsds_1")
  expect_identical(captured$args$backend, "auto")
  expect_identical(captured$args$cpu, 3L)
})

test_that("Type I-S runner resolves private helpers from its enclosing namespace", {
  isolated_namespace <- new.env(parent = baseenv())
  runner <- DNMB:::.dnmb_rebasefinder_predict_type1s
  environment(runner) <- isolated_namespace
  assign(".dnmb_rebasefinder_predict_type1s", runner, isolated_namespace)
  assign(
    ".dnmb_rebasefinder_empty_type1s_predictions",
    function() {
      data.frame(
        locus_tag = character(),
        type1s_prediction_status = character(),
        type1s_prediction_eligible = logical(),
        stringsAsFactors = FALSE
      )
    },
    isolated_namespace
  )
  assign(
    ".dnmb_rebasefinder_status_row",
    function(component, status, detail) {
      data.frame(
        component = component,
        status = status,
        detail = detail,
        stringsAsFactors = FALSE
      )
    },
    isolated_namespace
  )
  assign(
    ".dnmb_type1s_candidate_ids",
    function(genes, rm_comprehensive) "hsds_1",
    isolated_namespace
  )
  assign(
    ".dnmb_type1s_candidate_rows",
    function(gene_ids, candidate_ids) gene_ids %in% candidate_ids,
    isolated_namespace
  )
  assign(
    ".dnmb_type1s_clean_id",
    function(x) trimws(as.character(x)),
    isolated_namespace
  )
  assign(
    "dnmb_predict_type1s_recognition",
    function(candidates, ...) {
      data.frame(
        locus_tag = candidates$locus_tag,
        type1s_prediction_status = "predicted_complete",
        type1s_prediction_eligible = FALSE,
        stringsAsFactors = FALSE
      )
    },
    isolated_namespace
  )

  genes <- data.frame(
    locus_tag = "hsds_1",
    translation = paste(rep("A", 421), collapse = ""),
    stringsAsFactors = FALSE
  )
  result <- runner(
    genes,
    data.frame(locus_tag = "hsds_1", stringsAsFactors = FALSE),
    cache_dir = tempfile("type1s-isolated-cache-"),
    verbose = FALSE
  )

  expect_identical(result$predictions$locus_tag, "hsds_1")
  expect_identical(result$status$status, "ok")
})

test_that("Type I-S runner distinguishes failed, partial, and successful searches", {
  genes <- data.frame(
    locus_tag = c("hsds_1", "hsds_2"),
    product = "type I restriction-modification system specificity subunit",
    translation = paste(rep("A", 421), collapse = ""),
    stringsAsFactors = FALSE
  )
  rm <- data.frame(
    locus_tag = genes$locus_tag,
    rm_type = "Type I",
    subunit = "S",
    stringsAsFactors = FALSE
  )
  prediction_status <- c("search_failed", "search_failed")
  local_mocked_bindings(
    dnmb_predict_type1s_recognition = function(candidates, ...) {
      out <- DNMB:::.dnmb_type1s_empty_prediction(candidates$locus_tag)
      out$type1s_prediction_status <- prediction_status
      out
    },
    .package = "DNMB"
  )
  run <- function() {
    DNMB:::.dnmb_rebasefinder_predict_type1s(
      genes,
      rm,
      cache_dir = tempfile("type1s-cache-"),
      verbose = FALSE
    )
  }

  expect_identical(run()$status$status, "failed")
  prediction_status <- c("predicted_complete", "search_failed")
  expect_identical(run()$status$status, "partial")
  prediction_status <- c("predicted_complete", "no_compatible_trd_pair")
  expect_identical(run()$status$status, "ok")
})

test_that("failed Type I-S searches return typed empty predictions", {
  genes <- data.frame(
    locus_tag = "hsds_1",
    product = "HsdS type I specificity subunit",
    translation = paste(rep("A", 421), collapse = ""),
    stringsAsFactors = FALSE
  )
  rm <- data.frame(
    locus_tag = "hsds_1",
    rm_type = "Type I",
    subunit = "S",
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    dnmb_predict_type1s_recognition = function(...) stop("forced search failure"),
    .package = "DNMB"
  )

  result <- DNMB:::.dnmb_rebasefinder_predict_type1s(
    genes,
    rm,
    cache_dir = tempfile("type1s-cache-"),
    verbose = FALSE
  )

  expect_equal(nrow(result$predictions), 0L)
  expect_identical(
    names(result$predictions),
    names(DNMB:::.dnmb_type1s_empty_prediction(character()))
  )
  expect_identical(result$status$status, "failed")
  expect_match(result$status$detail, "forced search failure", fixed = TRUE)
})
