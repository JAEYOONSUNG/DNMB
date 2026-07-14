module_plot_fixture <- function(root) {
  gbff <- c(
    "LOCUS       NZ_CHR                 1000 bp    DNA     circular BCT 01-JAN-2000",
    "DEFINITION  Fixture chromosome, complete genome.",
    "ACCESSION   NZ_CHR",
    "     source          1..1000",
    "//",
    "LOCUS       NZ_PLS                  200 bp    DNA     circular BCT 01-JAN-2000",
    "DEFINITION  Fixture strain plasmid pFIX, complete sequence.",
    "ACCESSION   NZ_PLS",
    "     source          1..200",
    "//"
  )
  writeLines(gbff, file.path(root, "fixture.gbff"))
  data.frame(
    contig_number = c(1L, 2L),
    contig = rep("Fixture strain", 2),
    start = c(100, 20),
    end = c(180, 80),
    locus_tag = c("FIX_RS00001", "FIX_RS00002"),
    gene = c("geneA", "geneB"),
    product = c("chromosome protein", "plasmid protein"),
    stringsAsFactors = FALSE
  )
}

test_that("module plot replicons do not collapse identical contig descriptions", {
  root <- tempfile("module-replicons-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  genes <- module_plot_fixture(root)

  out <- DNMB:::.dnmb_module_replicon_plot_data(genes, root)

  expect_equal(nrow(out$contigs), 2)
  expect_identical(out$contigs$replicon_id, c("DNMB_CONTIG_001", "DNMB_CONTIG_002"))
  expect_equal(out$contigs$length_bp, c(1000, 200))
  expect_identical(out$contigs$replicon_short, c("Chromosome", "Plasmid pFIX"))
})

test_that("PADLOC plot data preserves shared system memberships and empty plasmids", {
  root <- tempfile("padloc-plot-data-")
  module_dir <- file.path(root, "dnmb_module_padloc")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  genes <- module_plot_fixture(root)

  id_map <- data.frame(
    query_id = "WP_FIX.1",
    locus_tag = "FIX_RS00001",
    contig_id = "DNMB_CONTIG_001",
    stringsAsFactors = FALSE
  )
  utils::write.table(
    id_map, file.path(module_dir, "fixture_id_map.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  raw <- data.frame(
    system.number = c(7L, 8L),
    seqid = rep("DNMB_CONTIG_001", 2),
    system = c("RM_type_II", "VSPR"),
    target.name = rep("WP_FIX.1", 2),
    protein.name = rep("MTase_II", 2),
    domain.iE.value = rep(1e-40, 2),
    target.coverage = rep(0.90, 2),
    hmm.coverage = rep(0.95, 2),
    start = rep(100, 2),
    end = rep(180, 2),
    hmm.name = rep("MTase_fixture", 2),
    check.names = FALSE
  )
  utils::write.csv(
    raw, file.path(module_dir, "fixture_padloc.csv"),
    row.names = FALSE, quote = TRUE
  )

  out <- DNMB:::.dnmb_padloc_plot_data(genes, root)

  expect_equal(nrow(out$hits), 2)
  expect_equal(length(unique(out$hits$locus_tag)), 1)
  expect_setequal(out$hits$PADLOC_system, c("RM_type_II", "VSPR"))
  expect_equal(nrow(out$contigs), 2)
  expect_false("DNMB_CONTIG_002" %in% out$hits$replicon_id)
})

test_that("dbAPIS plot data retains every accepted family profile on its replicon", {
  root <- tempfile("dbapis-plot-data-")
  module_dir <- file.path(root, "dnmb_module_dbapis")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  genes <- module_plot_fixture(root)

  raw <- data.frame(
    query = c("FIX_RS00001", "FIX_RS00001", "FIX_RS00002"),
    family_id = c("APIS001", "APIS002", "APIS003"),
    hit_label = c("", "narp", "ardc"),
    defense_type = c("", "NARP", "RM"),
    clan_id = c("CLAN1", "CLAN1", "CLAN2"),
    clan_defense_type = c("NARP", "NARP", "RM"),
    description = "",
    i_evalue = c(1e-20, 1e-12, 1e-8),
    score = c(80, 55, 30),
    hmm_coverage = c(0.8, 0.6, 0.4),
    query_coverage = c(0.9, 0.7, 0.5),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    raw, file.path(module_dir, "dbapis_hits.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  out <- DNMB:::.dnmb_dbapis_plot_data(genes, root)

  expect_equal(nrow(out$hits), 3)
  expect_equal(length(unique(out$hits$module_family_id)), 3)
  expect_equal(table(out$hits$replicon_id)[["DNMB_CONTIG_001"]], 2)
  expect_equal(table(out$hits$replicon_id)[["DNMB_CONTIG_002"]], 1)
  expect_equal(nrow(out$contigs), 2)
})

test_that("shared replicon label packing stays within the sequence boundary", {
  labels <- data.frame(
    replicon_id = rep("DNMB_CONTIG_002", 3),
    start = c(1, 45, 95),
    end = c(4, 55, 99),
    midpoint = c(2.5, 50, 97),
    feature_label = c("left_edge", "middle", "right_edge"),
    priority = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  contigs <- data.frame(
    replicon_id = "DNMB_CONTIG_002",
    length_bp = 100,
    stringsAsFactors = FALSE
  )

  packed <- DNMB:::.dnmb_module_pack_replicon_labels(
    labels, contigs, priority_col = "priority"
  )

  expect_true(all(packed$label_x - packed$label_half_width >= 0))
  expect_true(all(packed$label_x + packed$label_half_width <= 100))
})
