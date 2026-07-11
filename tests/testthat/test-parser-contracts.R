test_that("PhiSpy headerless coordinates retain the first prophage", {
  coordinates <- tempfile(fileext = ".tsv")
  on.exit(unlink(coordinates, force = TRUE), add = TRUE)
  writeLines(c(
    "pp1\tcontig_A\t100\t500\tattachment evidence",
    "pp2\tcontig_A\t900\t1400\t"
  ), coordinates)

  parsed <- DNMB:::dnmb_prophage_parse_coordinates(coordinates)
  standardized <- DNMB:::.dnmb_prophage_standardize_coordinates(parsed)

  expect_equal(nrow(parsed), 2L)
  expect_identical(standardized$prophage_id, c("pp1", "pp2"))
  expect_identical(standardized$contig, c("contig_A", "contig_A"))
  expect_equal(standardized$prophage_start, c(100, 900))
  expect_equal(standardized$prophage_end, c(500, 1400))
})

test_that("PhiSpy coordinates still accept an explicit header", {
  coordinates <- tempfile(fileext = ".tsv")
  on.exit(unlink(coordinates, force = TRUE), add = TRUE)
  writeLines(c(
    "Prophage number\tContig\tStart\tStop",
    "pp1\tcontig_A\t100\t500"
  ), coordinates)

  parsed <- DNMB:::dnmb_prophage_parse_coordinates(coordinates)
  standardized <- DNMB:::.dnmb_prophage_standardize_coordinates(parsed)

  expect_equal(nrow(parsed), 1L)
  expect_identical(standardized$prophage_id, "pp1")
  expect_equal(standardized$prophage_start, 100)
  expect_equal(standardized$prophage_end, 500)
})

test_that("InterProScan 15-column TSV preserves GO and pathway annotations", {
  interpro_dir <- tempfile("interproscan-")
  dir.create(interpro_dir)
  on.exit(unlink(interpro_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit({
    if (exists("InterProScan_table", envir = .GlobalEnv, inherits = FALSE)) {
      rm("InterProScan_table", envir = .GlobalEnv)
    }
    if (exists("InterProScan_site", envir = .GlobalEnv, inherits = FALSE)) {
      rm("InterProScan_site", envir = .GlobalEnv)
    }
  }, add = TRUE)

  fields <- c(
    "protein_1", "0123456789abcdef", "300", "Pfam", "PF00145",
    "C-5 cytosine-specific DNA methylase", "20", "280", "51.2", "T",
    "11-07-2026", "IPR001525", "C-5 DNA methylase",
    "GO:0006304|GO:0008168", "MetaCyc:DNA-METHYLATION"
  )
  writeLines(paste(fields, collapse = "\t"), file.path(interpro_dir, "result.tsv"))

  expect_message(
    parsed <- InterProScan_annotations(InterProScan_dir = interpro_dir),
    "saved to the R environment"
  )

  expect_identical(parsed$query, "protein_1")
  expect_identical(parsed[["GO annotations_Pfam"]], "GO:0006304|GO:0008168")
  expect_identical(parsed[["Pathways annotations_Pfam"]], "MetaCyc:DNA-METHYLATION")
  expect_identical(
    get("InterProScan_table", envir = .GlobalEnv)[["GO annotations_Pfam"]],
    "GO:0006304|GO:0008168"
  )
})

test_that("InterProScan 14-column TSV preserves GO and supplies an empty pathway", {
  interpro_dir <- tempfile("interproscan-go-only-")
  dir.create(interpro_dir)
  on.exit(unlink(interpro_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit({
    if (exists("InterProScan_table", envir = .GlobalEnv, inherits = FALSE)) {
      rm("InterProScan_table", envir = .GlobalEnv)
    }
    if (exists("InterProScan_site", envir = .GlobalEnv, inherits = FALSE)) {
      rm("InterProScan_site", envir = .GlobalEnv)
    }
  }, add = TRUE)

  fields <- c(
    "protein_2", "0123456789abcdef", "280", "Pfam", "PF00145",
    "DNA methylase", "10", "260", "42.0", "T", "11-07-2026",
    "IPR001525", "DNA methylase", "GO:0008168"
  )
  writeLines(paste(fields, collapse = "\t"), file.path(interpro_dir, "result.tsv"))

  expect_message(
    parsed <- InterProScan_annotations(InterProScan_dir = interpro_dir),
    "saved to the R environment"
  )
  expect_identical(parsed[["GO annotations_Pfam"]], "GO:0008168")
  expect_identical(parsed[["Pathways annotations_Pfam"]], "-")
})
