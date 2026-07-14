test_that("overview headings keep panel labels and titles in one shared style", {
  expect_identical(
    DNMB:::.dnmb_overview_title("A", "Inventory"),
    "A   Inventory"
  )

  heading_theme <- DNMB:::.dnmb_overview_title_theme()
  expect_identical(heading_theme$plot.title$family, DNMB:::.dnmb_plot_font_family())
  expect_identical(heading_theme$plot.title$face, "bold")
  expect_equal(heading_theme$plot.title$size, 11)
  expect_equal(heading_theme$plot.title$hjust, 0)
  expect_identical(heading_theme$plot.title.position, "panel")

  tagged <- DNMB:::.dnmb_overview_tag_plot(
    ggplot2::ggplot() + ggplot2::labs(title = "Genome context"),
    "B"
  )
  expect_identical(tagged$labels$title, "B   Genome context")
  expect_equal(tagged$theme$plot.title$size, 11)
})

test_that("Prophage region heading renders its tag and title on one baseline", {
  header <- DNMB:::.dnmb_prophage_region_header_plot(
    "Prophage 1 (47.7 kb)",
    "HIGH | att: TCGA",
    label = "A"
  )
  meta <- attr(header$plot, "dnmb_overview_heading")
  expect_identical(meta$label, "A")
  expect_identical(meta$title, "Prophage 1 (47.7 kb)")
  expect_equal(meta$size, 11)

  built <- ggplot2::ggplot_build(header$plot)
  rendered_labels <- unlist(lapply(built$data, function(x) {
    if ("label" %in% names(x)) x$label else character()
  }))
  expect_true("A   Prophage 1 (47.7 kb)" %in% rendered_labels)
})

test_that("DNMB PDF device embeds editable Arial when the font is installed", {
  skip_if_not(capabilities("cairo"))
  pdffonts <- Sys.which("pdffonts")
  fc_match <- Sys.which("fc-match")
  skip_if(!nzchar(pdffonts) || !nzchar(fc_match))

  matched <- system2(fc_match, "Arial", stdout = TRUE, stderr = TRUE)
  skip_if(!length(matched) || !grepl("Arial", matched[[1]], fixed = TRUE))

  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  DNMB:::.dnmb_plot_pdf_device(path, width = 3, height = 2)
  grid::grid.newpage()
  grid::grid.text("Editable Arial", gp = grid::gpar(fontfamily = DNMB:::.dnmb_plot_font_family()))
  grDevices::dev.off()

  info <- system2(pdffonts, shQuote(path), stdout = TRUE, stderr = TRUE)
  info <- paste(info, collapse = "\n")
  expect_match(info, "Arial")
  expect_match(info, "yes\\s+yes\\s+yes")
})

test_that("mRNAcal component heatmap supports legacy summary-only workbooks", {
  tbl <- data.frame(
    locus_tag = c("gene_a", "gene_b"),
    mRNAcal_tir_score = c(60, 40),
    mRNAcal_codon_efficiency_score = c(55, 35),
    mRNAcal_cai = c(0.8, 0.4),
    mRNAcal_tai = c(0.7, 0.3),
    stringsAsFactors = FALSE
  )

  components <- DNMB:::.dnmb_mrnacal_component_table(tbl)

  expect_true(all(c("TIR", "CodonEff", "CAI", "tAI") %in% components$component))
  expect_equal(
    components$score[components$locus_tag == "gene_a" & components$component == "CAI"],
    80
  )
  expect_equal(
    components$score[components$locus_tag == "gene_a" & components$component == "tAI"],
    70
  )
  expect_true(all(components$score >= 0 & components$score <= 100))
})
