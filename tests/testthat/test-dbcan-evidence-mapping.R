make_dbcan_genes <- function() {
  data.frame(
    locus_tag = c("Athe_2769", "Athe_2770", NA_character_, "context_gene"),
    protein_id = c("P2769", "P2770", "CAD88572.1", NA_character_),
    translation = c("MAAAAA", "MBBBBB", "MCCCCC", NA_character_),
    contig = rep("same display definition", 4),
    contig_number = c(1L, 2L, 2L, 2L),
    start = c(1L, 1L, 100L, 220L),
    end = c(18L, 18L, 117L, 300L),
    direction = c("+", "-", "+", "+"),
    stringsAsFactors = FALSE
  )
}

write_dbcan_overview_fixture <- function(path) {
  fixture <- data.frame(
    "Gene ID" = c("DNMBCAZY_0000001", "DNMBCAZY_0000002", "DNMBCAZY_0000003", "legacy_0004"),
    "EC#" = c("3.2.1.1", "-", "-", "-"),
    dbCAN_hmm = c("CBM34(1-80)+GH13_39(81-300)", "-", "-", "GH5(1-100)"),
    dbCAN_sub = c("CBM34_e14|CBM34_e14|GH13_39", "GH20_e1", "-", "-"),
    DIAMOND = c("GH13_39", "GH20", "GT4", "-"),
    "#ofTools" = c(3L, 2L, 1L, 1L),
    "Recommend Results" = c("CBM34_e14|CBM34_e14|GH13_39", "GH20_e1", "-", "-"),
    Substrate = c("-;alpha-glucan", "hostglycan", "-", "cellulose"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  utils::write.table(fixture, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(path)
}

read_dbcan_pdf_geometry <- function(path) {
  pdfinfo <- Sys.which("pdfinfo")
  if (!nzchar(pdfinfo)) return(NULL)

  info <- system2(pdfinfo, shQuote(normalizePath(path)), stdout = TRUE, stderr = TRUE)
  page_line <- grep("^Pages:", info, value = TRUE)
  size_line <- grep("^Page size:", info, value = TRUE)
  if (length(page_line) != 1L || length(size_line) != 1L) return(NULL)

  size_match <- regmatches(
    size_line,
    regexec(
      "^Page size:[[:space:]]*([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts",
      size_line
    )
  )[[1]]
  if (length(size_match) != 3L) return(NULL)

  list(
    pages = as.integer(sub("^Pages:[[:space:]]*", "", page_line)),
    width_pt = as.numeric(size_match[[2]]),
    height_pt = as.numeric(size_match[[3]])
  )
}

test_that("dbCAN synthetic IDs preserve locus suffixes and separate replicons", {
  genes <- make_dbcan_genes()
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)

  expect_identical(id_map$mapping_id[1:3], c("Athe_2769", "Athe_2770", "CAD88572.1"))
  expect_length(unique(id_map$dbcan_query_id), nrow(genes))
  expect_identical(id_map$dbcan_contig_key[1:2], c("replicon_0001", "replicon_0002"))

  fasta <- tempfile(fileext = ".faa")
  gff <- tempfile(fileext = ".gff")
  on.exit(unlink(c(fasta, gff), force = TRUE), add = TRUE)
  fasta_info <- DNMB:::.dnmb_dbcan_write_mapped_query_fasta(id_map, fasta)
  gff_info <- DNMB:::.dnmb_dbcan_write_query_gff(genes, gff, id_map = id_map)

  expect_identical(fasta_info$n, 3L)
  expect_identical(DNMB:::.dnmb_fasta_headers(fasta), id_map$dbcan_query_id[1:3])
  expect_identical(nrow(gff_info$contig_map), 2L)
  gff_body <- readLines(gff)[-1]
  expect_true(any(startsWith(gff_body, "ctg1\t")))
  expect_true(any(startsWith(gff_body, "ctg2\t")))
  expect_true(any(grepl("DNMBCAZY_0000004", gff_body, fixed = TRUE)))
})

test_that("dbCAN overview retains consensus, multi-domain calls, and evidence tiers", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  write_dbcan_overview_fixture(path)
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(make_dbcan_genes())
  overview <- DNMB:::dnmb_dbcan_parse_overview(path, id_map = id_map)

  expect_identical(overview$query[1:3], c("Athe_2769", "Athe_2770", "CAD88572.1"))
  expect_identical(overview$dbcan_all_families[[1]], "CBM34; GH13_39")
  expect_identical(overview$dbcan_primary_family[[1]], "GH13_39")
  expect_identical(overview$dbcan_family_count[[1]], 2L)
  expect_identical(overview$dbcan_evidence_tier, c("very_high", "high", "audit", "medium"))
  expect_identical(overview$dbcan_overview_substrate[1:2], c("alpha-glucan", "host glycan"))
  expect_identical(
    DNMB:::.dnmb_dbcan_family_tokens(overview$dbcan_all_families[[1]]),
    c("CBM34", "GH13_39")
  )
})

test_that("CAZy 3-zone map retains every GH domain from multi-domain calls", {
  genes <- data.frame(
    locus_tag = c("gene_1", "gene_2", "gene_3", "gene_4"),
    gene = c("amyA", "celA", "nagA", "cbm_only"),
    dbCAN_dbcan_all_families = c("CBM34; GH13_39", NA_character_, "GT4; GH20_e1", "CBM50"),
    dbCAN_dbcan_hit = c("GH13_39", NA_character_, "GT4", "GH99"),
    dbCAN_family_id = c("GH99", "GH5", "GT4", "GH5"),
    stringsAsFactors = FALSE
  )

  gh <- DNMB:::.dnmb_cct_3zone_extract_gh(genes)

  expect_identical(gh$locus_tag, c("gene_1", "gene_2", "gene_3"))
  expect_identical(gh$gh_family, c("GH13_39", "GH5", "GH20"))
})

test_that("CAZy transporter packing separates dense glyph intervals", {
  packed <- DNMB:::.dnmb_cct_pack_transporters_lane(
    center_x = 10,
    half_spans = rep(0.14, 8),
    lane_ranks = 8:1,
    label_widths = rep(0.5, 8),
    label_dx = rep(0, 8),
    desired_x = rep(10, 8),
    y_memb = 8.5
  )

  expect_equal(nrow(packed), 8L)
  expect_true(all(is.finite(packed$tx)))
  expect_true(all(is.finite(packed$ty)))
  expect_gte(length(unique(packed$row_id)), 2L)
  expect_lte(max(packed$row_id), 3L)
  expect_lte(max(abs(packed$ty - 8.5)), 0.12)
  expect_gte(min(diff(sort(unique(packed$ty)))), 0.12 - 1e-10)

  for (row_id in unique(packed$row_id)) {
    row <- packed[packed$row_id == row_id, , drop = FALSE]
    row <- row[order(row$tx), , drop = FALSE]
    if (nrow(row) < 2L) next
    edge_gap <- diff(row$tx) - 0.14 - 0.14
    expect_true(all(edge_gap >= 0.08 - 1e-10))
  }
})

test_that("CAZy transporter buses consolidate substrate memberships deterministically", {
  memberships <- data.frame(
    cs_id = c("A", "A", "A", "B", "C", "D", "E", "F"),
    lane_x = c(1, 1, 1, 2, 3, 4, 5, 6),
    tx_draw = c(4, 5, 6, 4, 4, 5, 5, 6),
    ty_draw = c(8.5, 8.62, 8.5, 8.5, 8.5, 8.62, 8.62, 8.5),
    confidence = c("high", "high", "medium", "high", "high", "high", "medium", "high"),
    locus_tag = c("shared_1", "shared_2", "medium_only", "shared_1",
                  "shared_1", "shared_2", "shared_2", "direct_1"),
    color = c("#A00000", "#A00000", "#A00000", "#0060A0",
              "#008060", "#7040A0", "#A07000", "#303030"),
    stringsAsFactors = FALSE
  )

  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(memberships, y_memb = 8.5)
  shuffled <- DNMB:::.dnmb_cct_transporter_bus_layout(
    memberships[c(6, 2, 8, 1, 7, 4, 3, 5), ], y_memb = 8.5
  )

  expect_equal(
    layout$groups[order(layout$groups$group), ],
    shuffled$groups[order(shuffled$groups$group), ],
    ignore_attr = TRUE
  )
  expect_equal(
    layout$trunks[order(layout$trunks$group), ],
    shuffled$trunks[order(shuffled$trunks$group), ],
    ignore_attr = TRUE
  )
  expect_equal(sum(layout$memberships$draw), 7L)
  expect_false(layout$memberships$draw[layout$memberships$target_key == "medium_only"])
  expect_equal(sum(layout$groups$direct), 2L)
  expect_equal(nrow(layout$direct), 1L)

  bus_groups <- layout$groups[!layout$groups$direct, , drop = FALSE]
  expect_true(all(bus_groups$bus_y >= 8.5 + 0.70 - 1e-10))
  expect_true(all(bus_groups$bus_y <= 8.5 + 1.12 + 1e-10))
  for (tier in unique(bus_groups$tier)) {
    row <- bus_groups[bus_groups$tier == tier, , drop = FALSE]
    row <- row[order(row$x_left), , drop = FALSE]
    if (nrow(row) < 2L) next
    expect_true(all(head(row$x_right, -1L) + 0.05 <= tail(row$x_left, -1L) + 1e-10))
  }

  # Three substrate lanes share one physical target but produce one trunk.
  expect_equal(sum(layout$trunks$group == "target:shared_1"), 1L)
  expect_equal(sum(layout$trunks$group == "target:shared_2"), 1L)
  expect_true(all(pmin(layout$stems$y, layout$stems$yend) > 8.5 + 0.20))
  expect_true(all(layout$trunks$y > 8.5 + 0.20))
  expect_true(all(layout$trunks$yend >= 8.5 + 0.11 - 1e-10))
})

test_that("CAZy transporter bus tiers stay inside their reserved vertical band", {
  memberships <- data.frame(
    cs_id = sprintf("S%02d", seq_len(20)),
    lane_x = seq_len(20),
    tx_draw = rep(20.5, 20),
    ty_draw = rep(8.5, 20),
    confidence = rep("high", 20),
    locus_tag = rep("shared_target", 20),
    color = rep("#52616B", 20),
    stringsAsFactors = FALSE
  )

  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(
    memberships,
    y_memb = 8.5,
    bus_offset = 0.62,
    max_bus_offset = 0.86
  )
  buses <- layout$groups[!layout$groups$direct, , drop = FALSE]

  expect_gt(length(unique(buses$tier)), 10L)
  expect_gte(min(buses$bus_y), 8.5 + 0.62 - 1e-10)
  expect_lte(max(buses$bus_y), 8.5 + 0.86 + 1e-10)
  expect_equal(nrow(layout$trunks), 1L)
})

test_that("CAZy transporter labels use deterministic tracks above the membrane", {
  labels <- data.frame(
    x = c(1, 1.02, 1.04, 3, 4, 4.02),
    y = rep(8.5, 6),
    label = sprintf("T%02d", seq_len(6)),
    color = rep("#102A43", 6),
    size = rep(1.6, 6),
    fontface = rep("bold", 6),
    hjust = rep(0.5, 6),
    priority = c(20, 10, 20, 10, 20, 10),
    stringsAsFactors = FALSE
  )

  laid_out <- DNMB:::.dnmb_cct_layout_transporter_labels(
    labels, xlim = c(0, 5), y_memb = 8.5
  )
  shuffled <- DNMB:::.dnmb_cct_layout_transporter_labels(
    labels[c(4, 1, 6, 2, 5, 3), ], xlim = c(0, 5), y_memb = 8.5
  )

  expect_equal(
    laid_out[order(laid_out$label), c("x_lab", "y_lab", "track_level")],
    shuffled[order(shuffled$label), c("x_lab", "y_lab", "track_level")],
    ignore_attr = TRUE
  )
  expect_true(all(laid_out$track_level %in% 0:2))
  expect_equal(
    laid_out$y_lab,
    8.5 + 0.40 + laid_out$track_level * 0.16,
    tolerance = 1e-10
  )
  expect_true(all(laid_out$y_lab > 8.5))
  expect_true(all(laid_out$x_lab - laid_out$label_width / 2 >= 0))
  expect_true(all(laid_out$x_lab + laid_out$label_width / 2 <= 5))

  for (track in unique(laid_out$track_level)) {
    row <- laid_out[laid_out$track_level == track, , drop = FALSE]
    row <- row[order(row$x_lab), , drop = FALSE]
    if (nrow(row) < 2L) next
    required <- (head(row$label_width, -1) + tail(row$label_width, -1)) / 2 + 0.08
    expect_true(all(diff(row$x_lab) >= required - 1e-10))
  }
})

test_that("CAZy transporter labels use fixed text and short orthogonal leaders", {
  labels <- data.frame(
    x = c(1, 1.02, 2), y = rep(8.5, 3),
    label = c("T01", "T02", "T03"),
    color = c("#102A43", "#52616B", "#102A43"),
    size = rep(1.6, 3), fontface = rep("bold", 3),
    hjust = rep(0.5, 3), priority = c(20, 10, 20),
    stringsAsFactors = FALSE
  )
  laid_out <- DNMB:::.dnmb_cct_layout_transporter_labels(
    labels, xlim = c(0, 3), y_memb = 8.5
  )
  leaders <- DNMB:::.dnmb_cct_transporter_label_leaders(laid_out)

  expect_equal(nrow(leaders), 3L * sum(laid_out$draw_leader))
  expect_false(any(
    laid_out$track_level == 0L &
      abs(laid_out$x_lab - laid_out$anchor_x) <= 0.04 &
      laid_out$draw_leader
  ))
  for (segment in split(leaders, leaders$group)) {
    expect_equal(segment$x[1], segment$x[2], tolerance = 1e-10)
    expect_equal(segment$y[2], segment$y[3], tolerance = 1e-10)
    expect_true(segment$y[1] < segment$y[2])
  }

  layer <- DNMB:::.dnmb_cct_transporter_text_layer(laid_out)
  expect_s3_class(layer, "Layer")
  expect_true(inherits(layer$geom, "GeomText"))
  expect_false(inherits(layer$geom, "GeomTextRepel"))
  expect_no_error(ggplot2::ggplot_build(ggplot2::ggplot() + layer))
})

test_that("CAZy cleavage scissors have handles, crossed blades, and a pivot", {
  geom0 <- DNMB:::.dnmb_scissors_geometry_v2(2, 3, size = 0.1, angle = 0)
  geom90 <- DNMB:::.dnmb_scissors_geometry_v2(2, 3, size = 0.1, angle = 90)

  expect_equal(nrow(geom0$handles), 2L)
  expect_equal(nrow(geom0$shanks), 2L)
  expect_equal(nrow(geom0$blades), 2L)
  expect_equal(unname(unlist(geom0$pivot[1, c("x", "y")])), c(2, 3))
  expect_true(all(geom0$handles$x < geom0$pivot$x))
  expect_true(all(geom0$tips$x > geom0$pivot$x))

  distance_from_pivot <- function(points, pivot) {
    sqrt((points$x - pivot$x)^2 + (points$y - pivot$y)^2)
  }
  expect_equal(
    distance_from_pivot(geom0$handles, geom0$pivot),
    distance_from_pivot(geom90$handles, geom90$pivot),
    tolerance = 1e-10
  )
  layers <- DNMB:::.dnmb_scissors_grob_v2(2, 3, size = 0.1)
  expect_length(layers, 6L)
  expect_equal(sum(vapply(layers, function(layer) inherits(layer$geom, "GeomCustomAnn"), logical(1))), 3L)
  expect_false(any(vapply(layers, function(layer) inherits(layer$geom, "GeomPolygon"), logical(1))))
})

test_that("CAZy circular SNFG symbols use native circles rather than polygons", {
  circle_layers <- DNMB:::.dnmb_snfg_symbol_layers_v2(1, 1, "glucose", size = 0.1)
  square_layers <- DNMB:::.dnmb_snfg_symbol_layers_v2(1, 1, "glcnac", size = 0.1)

  expect_true(all(vapply(
    circle_layers,
    function(layer) inherits(layer$geom, "GeomCustomAnn"),
    logical(1)
  )))
  expect_false(any(vapply(
    circle_layers,
    function(layer) inherits(layer$geom, "GeomPolygon"),
    logical(1)
  )))
  expect_true(any(vapply(
    square_layers,
    function(layer) inherits(layer$geom, "GeomPolygon"),
    logical(1)
  )))
})

test_that("CAZy final annotation layer accepts multiline labels and node obstacles", {
  labels <- data.frame(
    x = c(1, 1.02), y = c(1, 1.01),
    label = c("transporter\ngene_1", "GH13_39\ngene_2"),
    color = c("#111111", "#C62828"), size = c(0.8, 1.0),
    fontface = "bold", hjust = 0.5, priority = c(10, 9),
    nudge_x = 0, nudge_y = c(0.1, -0.1),
    stringsAsFactors = FALSE
  )
  obstacles <- DNMB:::.dnmb_cct_obstacle_perimeter(1, 1, rx = 0.2, ry = 0.1)
  expect_equal(nrow(obstacles), 9L)
  expect_equal(range(obstacles$x), c(0.8, 1.2), tolerance = 1e-10)
  expect_equal(range(obstacles$y), c(0.9, 1.1), tolerance = 1e-10)

  layer <- DNMB:::.dnmb_cct_final_text_layer(
    labels,
    obstacles = obstacles,
    xlim = c(0, 2), ylim = c(0, 2)
  )

  expect_s3_class(layer, "Layer")
  expect_no_error(ggplot2::ggplot_build(ggplot2::ggplot() + layer))
})

test_that("CAZy transporter entities preserve supported shared memberships", {
  transporters <- data.frame(
    locus_tag = c(
      rep("AADSM_000052", 5),
      rep("AADSM_001386", 4),
      NA_character_
    ),
    pathway = c(
      "arabinose", "galactose", "xylose", "xylose", "sucrose",
      "arabinose", "galactose", "xylose", "maltose",
      "arabinose"
    ),
    step = c(
      "araE", "galP", "xylT", "xylT", "MFS-glucose",
      "araE", "galP", "xylT", "MFS-glucose",
      "araE"
    ),
    cs_id = c(
      "Arabinose", "Galactose", "Xylose", "Xylose", "Sucrose",
      "Arabinose", "Galactose", "Xylose", "Maltose",
      "Arabinose"
    ),
    confidence = c(
      "high", "high", "high", "medium", "none",
      "high", "high", "none", "low",
      "high"
    ),
    step_score = c(2, 2, 2, 1, 0, 2, 2, 2, 0, 2),
    stringsAsFactors = FALSE
  )

  result <- DNMB:::.dnmb_cct_transporter_entities(transporters)

  expect_setequal(result$entities$locus_tag, c("AADSM_000052", "AADSM_001386"))
  expect_equal(nrow(result$memberships), 6L)
  membership_counts <- table(result$memberships$locus_tag)
  expect_setequal(names(membership_counts), c("AADSM_000052", "AADSM_001386"))
  expect_identical(as.integer(membership_counts), c(3L, 3L))
  expect_setequal(result$memberships$step, c("araE", "galP", "xylT"))
  expect_false(any(result$memberships$step == "MFS-glucose"))
  expect_true(all(result$entities$shared))

  for (label in result$entities$label) {
    expect_match(label, "araE", fixed = TRUE)
    expect_match(label, "galP", fixed = TRUE)
    expect_match(label, "xylT", fixed = TRUE)
    expect_match(label, "Arabinose", fixed = TRUE)
    expect_match(label, "Galactose", fixed = TRUE)
    expect_match(label, "Xylose", fixed = TRUE)
    expect_match(label, "Shared: yes", fixed = TRUE)
    expect_false(grepl("+N genes", label, fixed = TRUE))
    expect_false(grepl("...", label, fixed = TRUE))
  }
})

test_that("CAZy transporter recognition includes xylT and SSS models", {
  expect_true(DNMB:::.dnmb_cct_is_transport_like_step("xylT"))
  expect_true(DNMB:::.dnmb_cct_is_transport_like_step("SSS-glucose"))
  expect_true(DNMB:::.dnmb_cct_is_transport_like_step(
    "BT0355", gene_name = "sodium:solute symporter family protein"
  ))
})

test_that("CAZy pathway presence separates active, partial, and reference metabolism", {
  steps <- data.frame(
    pathway_id = c("xylose", "xylose", "arabinose", "arabinose", "sucrose"),
    step_id = c("xylA", "xylB", "araE", "araA", "sucA"),
    confidence = c("high", "high", "high", "medium", "high"),
    locus_tag = paste0("LOC_", seq_len(5)),
    stringsAsFactors = FALSE
  )
  stats <- data.frame(
    pathway_id = c("xylose", "arabinose", "sucrose", "glucosamine"),
    fraction = c(1, 0.5, 1, 0),
    stringsAsFactors = FALSE
  )

  presence <- DNMB:::.dnmb_cct_pathway_presence(
    steps, stats,
    valid_pathways = c("xylose", "arabinose", "sucrose", "glucosamine")
  )
  state <- setNames(presence$cytoplasm_status, presence$pathway_id)

  expect_identical(unname(state["xylose"]), "active")
  expect_identical(unname(state["arabinose"]), "active")
  expect_identical(unname(state["sucrose"]), "partial")
  expect_identical(unname(state["glucosamine"]), "reference")
  expect_equal(presence$n_transport[presence$pathway_id == "arabinose"], 1L)
  expect_equal(presence$n_intracellular[presence$pathway_id == "xylose"], 2L)
})

test_that("CAZy route labels never collapse distinct genes at a shared target", {
  labels <- data.frame(
    x = c(1, 1.1, 2), y = c(1, 1.1, 2),
    label = c("stepA\nA", "stepD\nD", "stepA\nA"),
    base_label = c("stepA\nA", "stepD\nD", "stepA\nA"),
    step_id = c("stepA", "stepD", "stepA"),
    locus_tag = c("A", "D", "A"),
    member_loci = c("A;B;C", "D;E;F", "A"),
    member_keys = c("sA::A;sB::B;sC::C", "sD::D;sE::E;sF::F", "sA::A"),
    target_id = c("target_1", "target_1", "target_2"),
    priority = c(10, 9, 8),
    stringsAsFactors = FALSE
  )

  collapsed <- DNMB:::.dnmb_cct_collapse_route_labels(labels)
  target_1 <- collapsed[collapsed$target_id == "target_1", , drop = FALSE]
  target_2 <- collapsed[collapsed$target_id == "target_2", , drop = FALSE]

  expect_equal(nrow(collapsed), 3L)
  expect_setequal(target_1$label, c("stepA\nA", "stepD\nD"))
  expect_false(any(grepl("\\(\\+[0-9]+ genes?\\)", collapsed$label)))
  expect_setequal(target_1$locus_tag, c("A", "D"))
  expect_identical(target_2$label, "stepA\nA")
  expect_identical(target_2$member_loci, "A")
})

test_that("dbCAN output maps consensus-only genes by original gene row", {
  genes <- make_dbcan_genes()[1:3, ]
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)
  overview_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(overview_path, force = TRUE), add = TRUE)
  write_dbcan_overview_fixture(overview_path)
  overview <- DNMB:::dnmb_dbcan_parse_overview(overview_path, id_map = id_map)

  raw <- DNMB:::.dnmb_dbcan_empty_hits()
  raw <- rbind(
    raw,
    data.frame(
      query = "DNMBCAZY_0000001", profile_id = "GH13_39.hmm",
      profile_length = 300L, gene_length = 320L, evalue = 1e-40,
      profile_start = 1L, profile_end = 290L, gene_start = 10L,
      gene_end = 300L, coverage = 0.967, family_id = "GH13_39",
      hit_label = "GH13_39", substrate_label = NA_character_
    )
  )
  raw <- DNMB:::.dnmb_dbcan_restore_query_map(raw, id_map)
  hmm_hits <- DNMB:::dnmb_dbcan_normalize_hits(raw)
  hits <- DNMB:::.dnmb_dbcan_merge_overview_hits(hmm_hits, overview[1:3, ])
  out <- DNMB:::.dnmb_dbcan_output_table(
    genes = genes,
    hits = hits,
    overview = overview[1:3, ],
    id_map = id_map
  )

  expect_identical(out$locus_tag[1:2], c("Athe_2769", "Athe_2770"))
  expect_identical(out$dbcan_hit[1:3], c("CBM34; GH13_39", "GH20", "GT4"))
  expect_identical(out$dbcan_hmm_domain_count[[1]], 1L)
  expect_true(is.na(out$dbcan_hmm_domain_count[[2]]))
  expect_identical(sort(unique(hits$query)), sort(overview$query[1:3]))
  expect_false(hits$typing_eligible[hits$query == "CAD88572.1"][[1]])
})

test_that("dbCAN HMM fallback keeps duplicate locus tags separate by gene row", {
  genes <- data.frame(
    locus_tag = c("dup", "dup"),
    protein_id = c("p1", "p2"),
    translation = c("MAAAA", "MBBBB"),
    contig = c("same definition", "same definition"),
    contig_number = c(1L, 2L),
    start = c(1L, 1L),
    end = c(15L, 15L),
    direction = c("+", "+"),
    stringsAsFactors = FALSE
  )
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)
  hits <- tibble::tibble(
    query = c("dup", "dup"),
    source = "dbcan",
    family_system = "dbCAN",
    family_id = c("GH1", "GT2"),
    hit_label = c("GH1", "GT2"),
    enzyme_role = "CAZyme",
    evidence_mode = "hmm",
    substrate_label = NA_character_,
    support = "HMM",
    typing_eligible = TRUE,
    evalue = c(1e-20, 1e-30),
    coverage = c(0.8, 0.9),
    gene_start = c(1L, 1L),
    gene_end = c(50L, 60L),
    dbcan_gene_row = c(1L, 2L)
  )

  out <- DNMB:::.dnmb_dbcan_output_table(genes, hits, id_map = id_map)
  expect_identical(out$family_id, c("GH1", "GT2"))
  expect_identical(out$dbcan_hit, c("GH1", "GT2"))
})

test_that("dbCAN row mapping survives module merge and append", {
  genes <- data.frame(
    locus_tag = c("dup", "dup", NA_character_, NA_character_),
    protein_id = c("p1", "p2", "p3", "p4"),
    translation = c("MAAA", "MBBB", "MCCC", "MDDD"),
    contig = rep("same definition", 4),
    contig_number = c(1L, 2L, 2L, 3L),
    start = c(1L, 1L, 100L, 100L),
    end = c(12L, 12L, 111L, 111L),
    direction = "+",
    stringsAsFactors = FALSE
  )
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)
  hits <- tibble::tibble(
    query = c("dup", "dup", "p3", "p4"), source = "dbcan",
    family_system = "dbCAN", family_id = c("GH1", "GT2", "CE1", "CBM2"),
    hit_label = c("GH1", "GT2", "CE1", "CBM2"), enzyme_role = "CAZyme",
    evidence_mode = "hmm", substrate_label = NA_character_, support = "HMM",
    typing_eligible = TRUE, evalue = c(1e-20, 1e-30, 1e-25, 1e-22),
    coverage = 0.8, gene_start = 1L, gene_end = 50L,
    dbcan_gene_row = 1:4
  )
  dbcan_table <- DNMB:::.dnmb_dbcan_output_table(genes, hits, id_map = id_map)
  other_table <- genes[, intersect(DNMB:::dnmb_backbone_columns(), names(genes)), drop = FALSE]
  other_table$locus_tag <- c("legacy_dup", "legacy_dup", NA_character_, NA_character_)
  other_table$marker <- paste0("m", seq_len(nrow(other_table)))
  runs <- list(
    Other = structure(list(database = "Other", output_table = other_table), class = "dnmb_module_run"),
    dbCAN = structure(list(database = "dbCAN", output_table = dbcan_table), class = "dnmb_module_run")
  )

  appended <- DNMB:::append_module_results(genes, runs)
  expect_identical(appended$dbCAN_family_id, c("GH1", "GT2", "CE1", "CBM2"))
  expect_identical(appended$dbCAN_dbcan_gene_row, 1:4)

  shuffled <- genes[c(2, 1, 4, 3), , drop = FALSE]
  shuffled_appended <- DNMB:::append_module_results(shuffled, runs)
  expect_identical(shuffled_appended$protein_id, c("p2", "p1", "p4", "p3"))
  expect_identical(shuffled_appended$dbCAN_family_id, c("GT2", "GH1", "CBM2", "CE1"))
  expect_identical(shuffled_appended$dbCAN_dbcan_gene_row, c(2L, 1L, 4L, 3L))
})

test_that("dbCAN domtblout applies inclusive coverage before overlap filtering", {
  make_line <- function(query, profile, evalue, profile_from, profile_to, gene_from, gene_to) {
    fields <- rep("0", 19L)
    fields[c(1, 3, 4, 6, 13, 16, 17, 18, 19)] <- c(
      query, "200", profile, "100", format(evalue, scientific = TRUE),
      profile_from, profile_to, gene_from, gene_to
    )
    paste(fields, collapse = " ")
  }
  path <- tempfile(fileext = ".domtblout")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  writeLines(c(
    make_line("gene_0001", "GH1.hmm", 1e-40, 1, 20, 1, 80),
    make_line("gene_0001", "GH2.hmm", 1e-20, 1, 80, 10, 90),
    make_line("gene_0002", "GH3.hmm", 1e-25, 1, 35, 1, 35)
  ), path)

  hits <- DNMB:::dnmb_dbcan_parse_domtblout(path, evalue_threshold = 1e-15, coverage_threshold = 0.35)
  expect_identical(sort(hits$query), c("gene_0001", "gene_0002"))
  expect_identical(hits$family_id[hits$query == "gene_0001"], "GH2")
  expect_equal(hits$coverage[hits$query == "gene_0002"], 0.35)
})

test_that("dbCAN comparative tokens expand all domains once per gene", {
  calls <- c("CBM34_e14|CBM34_e14|GH13_39", "GH18(2-200)+CBM5_e25")
  expect_identical(
    DNMB:::.dnmb_comparative_dbcan_token(calls, level = "family"),
    c("CBM34", "GH13", "GH18", "CBM5")
  )
  expect_identical(
    DNMB:::.dnmb_comparative_dbcan_token(calls, level = "class"),
    c("CBM", "GH", "GH", "CBM")
  )
})

test_that("dbCAN family substrate prior is explicit and lower priority", {
  mapping <- tempfile(fileext = ".tsv")
  on.exit(unlink(mapping, force = TRUE), add = TRUE)
  writeLines(c(
    "Substrate_high_level\tSubstrate_curated\tFamily\tName\tEC_Number\tnew_Substrate_high_level\ttype",
    "xylan\txylan\tGH120\tenzyme\t3.2.1.-\t\t",
    "starch\tstarch\tGH13\tenzyme\t3.2.1.-\talpha-glucan\tCHANGED",
    "lignin\tlignin\tAA1\tenzyme\t1.1.1.-\t\tREMOVED",
    "\t\tGH3\t\"b-glucoside phosphorylase\t\"\t2.4.1.-\tbeta-glucan\tNEW"
  ), mapping)
  overview <- tibble::tibble(dbcan_all_families = c("GH120", "GH13_39", "AA1", "GH3"))
  annotated <- DNMB:::.dnmb_dbcan_add_family_substrate_prior(overview, mapping)
  expect_identical(
    annotated$dbcan_family_substrate_prior,
    c("xylan", "alpha-glucan", NA_character_, "beta-glucan")
  )
})

test_that("dbCAN consensus-only overview does not require HMM coordinates", {
  genes <- make_dbcan_genes()[2, , drop = FALSE]
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  fixture <- data.frame(
    "Gene ID" = id_map$dbcan_query_id,
    dbCAN_hmm = "-",
    dbCAN_sub = "GH20_e1",
    DIAMOND = "GH20",
    "#ofTools" = 2L,
    "Recommend Results" = "GH20_e1",
    check.names = FALSE
  )
  utils::write.table(fixture, path, sep = "\t", row.names = FALSE, quote = FALSE)
  overview <- DNMB:::dnmb_dbcan_parse_overview(path, id_map = id_map)
  hits <- DNMB:::.dnmb_dbcan_merge_overview_hits(
    DNMB:::.dnmb_module_empty_optional_long_table(),
    overview
  )

  out <- DNMB:::.dnmb_dbcan_output_table(genes, hits, overview = overview, id_map = id_map)
  expect_identical(out$dbcan_hit, "GH20")
  expect_false("dbcan_hmm_domain_count" %in% names(out))
})

test_that("dbCAN CGC overview keeps dense content on one landscape sheet", {
  n_clusters <- 25L
  rows <- lapply(seq_len(n_clusters), function(i) {
    data.frame(
      locus_tag = paste0("gene_", i, "_", 1:3),
      protein_id = paste0("protein_", i, "_", 1:3),
      contig = "same definition",
      contig_number = 1L,
      start = i * 10000L + c(1L, 1001L, 2001L),
      end = i * 10000L + c(900L, 1900L, 2900L),
      direction = c("+", "-", "+"),
      dbCAN_dbcan_all_families = c("GH13_39; CBM34", NA_character_, NA_character_),
      dbCAN_dbcan_evidence_tier = c("high", NA_character_, NA_character_),
      dbCAN_substrate_label = c("starch", "starch", "starch"),
      dbCAN_dbcan_contig_key = "replicon_0001",
      dbCAN_dbcan_cgc_id = paste0("replicon_0001|CGC", i),
      dbCAN_dbcan_cgc_gene_type = c("CAZyme", "prodoric", "Sulfatase"),
      dbCAN_dbcan_cgc_protein_family = c("CAZyme|GH13_39", "prodoric|TF", "Sulfatase|S1"),
      stringsAsFactors = FALSE
    )
  })
  genes <- dplyr::bind_rows(rows)
  inputs <- DNMB:::.dnmb_dbcan_read_plot_inputs(genes, tempdir())
  expect_true(all(c("CAZyme", "TF", "Sulfatase") %in% unique(inputs$cgc$gene_type)))
  pages <- DNMB:::.dnmb_dbcan_cgc_pages(inputs$cgc, clusters_per_page = 12L)
  expect_length(pages, 3L)
  sheet <- DNMB:::.dnmb_dbcan_one_page_layout(inputs)
  expect_s3_class(sheet$plot, "ggplot")
  expect_identical(sheet$n_cgc, n_clusters)
  expect_gte(sheet$width, 14)
  expect_gt(sheet$height, 9)

  output_dir <- tempfile("dbcan-plot-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  plot_dir <- file.path(output_dir, "visualizations")
  dir.create(plot_dir)
  cazy_path <- file.path(plot_dir, "CAZy_overview.pdf")
  writeLines("metabolic-map-sentinel", cazy_path, useBytes = TRUE)
  cazy_md5 <- unname(tools::md5sum(cazy_path))

  result <- DNMB:::.dnmb_plot_dbcan_module(genes, output_dir)

  expect_true(file.exists(result$pdf))
  expect_identical(basename(result$pdf), "dbcan_cgc_overview.pdf")
  expect_true(file.exists(cazy_path))
  expect_identical(unname(tools::md5sum(cazy_path)), cazy_md5)
  expect_gt(file.info(result$pdf)$size, 1000)
  expect_identical(result$pages, 1L)

  geometry <- read_dbcan_pdf_geometry(result$pdf)
  if (!is.null(geometry)) {
    expect_identical(geometry$pages, 1L)
    expect_equal(geometry$width_pt, sheet$width * 72, tolerance = 1)
    expect_equal(geometry$height_pt, sheet$height * 72, tolerance = 1)
    expect_gt(geometry$width_pt, geometry$height_pt)
  }

  pdftotext <- Sys.which("pdftotext")
  if (nzchar(pdftotext)) {
    page_text_path <- tempfile(fileext = ".txt")
    on.exit(unlink(page_text_path, force = TRUE), add = TRUE)
    status <- system2(
      pdftotext,
      c("-f", "1", "-l", "1", shQuote(normalizePath(result$pdf)), shQuote(page_text_path)),
      stdout = FALSE,
      stderr = FALSE
    )
    expect_identical(status, 0L)
    page_text <- paste(readLines(page_text_path, warn = FALSE), collapse = "\n")
    expect_match(page_text, "CGC1([^0-9]|$)")
    expect_match(page_text, paste0("CGC", n_clusters, "([^0-9]|$)"))
  }
})

test_that("module plotting keeps the CAZy map separate from dbCAN CGC", {
  cazy_stub <- NA_character_
  noop_plot <- function(...) NULL
  local_mocked_bindings(
    .dnmb_plot_gapmind_aa_pathway_map = noop_plot,
    .dnmb_plot_cazy_carbon_transport_map = function(genbank_table, output_dir, file_stub) {
      cazy_stub <<- file_stub
      list(pdf = file.path(output_dir, paste0(file_stub, ".pdf")))
    },
    .dnmb_plot_defensefinder_module = noop_plot,
    .dnmb_plot_dbapis_module = noop_plot,
    .dnmb_plot_acrfinder_module = noop_plot,
    .dnmb_plot_mrnacal_module = noop_plot,
    .dnmb_plot_padloc_module = noop_plot,
    .dnmb_plot_defensepredictor_module = noop_plot,
    .dnmb_plot_iselement_module = noop_plot,
    .dnmb_plot_dbcan_module = function(genbank_table, output_dir) {
      list(pdf = file.path(output_dir, "dbcan_cgc_overview.pdf"))
    },
    .dnmb_plot_merops_module = noop_plot,
    .dnmb_plot_pazy_module = noop_plot,
    .dnmb_plot_rebasefinder_module = noop_plot,
    .dnmb_plot_integrated_defense_module = noop_plot,
    .package = "DNMB"
  )

  plots <- DNMB:::dnmb_render_module_plots(
    data.frame(locus_tag = "gene_1", stringsAsFactors = FALSE),
    output_dir = tempdir()
  )

  expect_identical(cazy_stub, "CAZy_overview")
  expect_identical(basename(plots$CAZyTransport$pdf), "CAZy_overview.pdf")
  expect_identical(basename(plots$dbCAN$pdf), "dbcan_cgc_overview.pdf")
})

test_that("dbCAN summary renders CGCs when no CAZyme family is assigned", {
  genes <- make_dbcan_genes()[1:2, ]
  genes$dbCAN_dbcan_cgc_id <- c("replicon_0001|CGC1", "replicon_0002|CGC1")
  genes$dbCAN_dbcan_cgc_gene_type <- c("TC", "STP")
  genes$dbCAN_dbcan_cgc_protein_family <- c("TC|transporter", "STP|sensor")

  inputs <- DNMB:::.dnmb_dbcan_read_plot_inputs(genes, tempdir())
  expect_equal(nrow(inputs$cazy), 0L)
  expect_s3_class(DNMB:::.dnmb_dbcan_summary_page(inputs), "ggplot")
})

test_that("dbCAN genome map paginates hit-containing contigs only", {
  genes <- data.frame(
    locus_tag = paste0("gene_", 1:20),
    protein_id = paste0("protein_", 1:20),
    contig = paste0("contig ", 1:20),
    contig_number = 1:20,
    start = 1L,
    end = 900L,
    direction = "+",
    dbCAN_dbcan_all_families = c(rep("GH1", 17), rep(NA_character_, 3)),
    dbCAN_dbcan_evidence_tier = "high",
    stringsAsFactors = FALSE
  )
  inputs <- DNMB:::.dnmb_dbcan_read_plot_inputs(genes, tempdir())
  groups <- DNMB:::.dnmb_dbcan_hit_contig_groups(inputs, contigs_per_page = 8L)

  expect_length(groups, 3L)
  expect_lte(max(lengths(groups)), 8L)
  expect_setequal(unlist(groups), sprintf("replicon_%04d", 1:17))
})

test_that("dbCAN plot separates broad family substrate priors", {
  domains <- tibble::tibble(
    family = c("GH1", "GH2"),
    substrate = c("xylan", paste0("prior", 1:13, collapse = "; ")),
    substrate_source = c("overview", "family_prior")
  )
  supported <- DNMB:::.dnmb_dbcan_matrix_pages(
    domains,
    substrates_per_page = 6L,
    substrate_mode = "supported"
  )
  prior <- DNMB:::.dnmb_dbcan_matrix_pages(
    domains,
    substrates_per_page = 6L,
    substrate_mode = "family_prior"
  )

  expect_length(supported, 1L)
  expect_setequal(as.character(supported[[1]]$data$substrate), c("xylan", "unknown"))
  expect_length(prior, 3L)
  expect_false(any(vapply(prior, function(page) "unknown" %in% as.character(page$data$substrate), logical(1))))
})

test_that("dbCAN plot overview remaps duplicate locus tags by gene row", {
  genes <- make_dbcan_genes()[1:2, ]
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)
  genes$dbCAN_dbcan_gene_row <- 1:2
  genes$dbCAN_dbcan_all_families <- c("GH1", "GT2")
  genes$dbCAN_dbcan_evidence_tier <- "medium"
  output_dir <- tempfile("dbcan-plot-map-")
  run_dir <- file.path(output_dir, "dnmb_module_dbcan", "run_dbcan")
  dir.create(run_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  DNMB:::.dnmb_dbcan_write_gene_map(
    id_map,
    file.path(output_dir, "dnmb_module_dbcan", "dbcan_gene_id_map.tsv")
  )
  fixture <- data.frame(
    "Gene ID" = id_map$dbcan_query_id,
    dbCAN_hmm = c("GH3", "GT4"),
    dbCAN_sub = c("GH3_e1", "GT4_e1"),
    DIAMOND = c("GH3", "GT4"),
    "#ofTools" = 3L,
    "Recommend Results" = c("GH3_e1", "GT4_e1"),
    check.names = FALSE
  )
  utils::write.table(
    fixture,
    file.path(run_dir, "overview.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  inputs <- DNMB:::.dnmb_dbcan_read_plot_inputs(genes, output_dir)
  expect_identical(inputs$genes$dbcan_all_families, c("GH3", "GT4"))
  expect_identical(inputs$genes$dbcan_gene_row, 1:2)
})

test_that("dbCAN plot CGC fallback restores synthetic protein IDs", {
  genes <- make_dbcan_genes()[1:2, ]
  id_map <- DNMB:::.dnmb_dbcan_build_gene_map(genes)
  genes$dbCAN_dbcan_gene_row <- 1:2
  output_dir <- tempfile("dbcan-cgc-map-")
  run_dir <- file.path(output_dir, "dnmb_module_dbcan", "run_dbcan")
  dir.create(run_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  DNMB:::.dnmb_dbcan_write_gene_map(
    id_map,
    file.path(output_dir, "dnmb_module_dbcan", "dbcan_gene_id_map.tsv")
  )
  fixture <- data.frame(
    "CGC#" = c("CGC1", "CGC2"),
    "Gene Type" = c("TC", "STP"),
    "Contig ID" = c("ctg1", "ctg2"),
    "Protein ID" = id_map$dbcan_query_id,
    "Gene Start" = genes$start,
    "Gene Stop" = genes$end,
    "Gene Strand" = genes$direction,
    "Gene Annotation" = c("TC|transporter", "STP|sensor"),
    check.names = FALSE
  )
  utils::write.table(
    fixture,
    file.path(run_dir, "cgc_standard_out.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  inputs <- DNMB:::.dnmb_dbcan_read_plot_inputs(genes, output_dir)
  expect_identical(inputs$cgc$dbcan_gene_row, 1:2)
  expect_identical(inputs$cgc$dbcan_cgc_id, c("replicon_0001|CGC1", "replicon_0002|CGC2"))
  expect_identical(inputs$cgc$gene_type, c("TC", "STP"))
})

test_that("dbCAN cache signature includes the mapping contract", {
  expect_identical(DNMB:::.dnmb_dbcan_pipeline_contract_version(), 2L)
  body_text <- paste(deparse(body(DNMB:::dnmb_dbcan_run_standalone)), collapse = "\n")
  expect_match(body_text, "--threads", fixed = TRUE)
})

test_that("dbCAN cache signature tracks executable version and path", {
  tool_version <- "5.2.9"
  tool_path <- "/opt/dbcan-5.2.9/run_dbcan"
  local_mocked_bindings(
    .dnmb_db_manifest_identity = function(...) list(database = "fixture"),
    dnmb_detect_binary = function(...) list(found = TRUE, path = tool_path),
    .dnmb_dbcan_tool_version = function() tool_version,
    .package = "DNMB"
  )
  first <- DNMB:::.dnmb_collect_module_db_signatures("dbCAN")$dbCAN
  tool_version <- "5.3.0"
  tool_path <- "/opt/dbcan-5.3.0/run_dbcan"
  second <- DNMB:::.dnmb_collect_module_db_signatures("dbCAN")$dbCAN

  expect_identical(first$run_dbcan_version, "5.2.9")
  expect_identical(second$run_dbcan_version, "5.3.0")
  expect_false(identical(first, second))
})

test_that("dbCAN standalone readiness requires all v5 CGC assets", {
  root <- tempfile("dbcan-bundle-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  paths <- DNMB:::.dnmb_dbcan_supporting_paths(root)
  files <- c(
    paths$hmm_txt, paste0(paths$hmm_txt, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$fam_substrate_mapping, paths$pul_dmnd, paths$pul_excel,
    paths$tcdb_dmnd, paths$tf_hmm, paths$tf_dmnd, paths$stp_hmm,
    paste0(paths$stp_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paste0(paths$tf_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$cazy_dmnd, paths$peptidase_dmnd, paths$sulfatlas_dmnd,
    paths$dbcan_sub_hmm,
    paste0(paths$dbcan_sub_hmm, c(".h3f", ".h3i", ".h3m", ".h3p"))
  )
  file.create(files)
  dir.create(paths$pul_dir)
  expect_true(DNMB:::.dnmb_dbcan_standalone_ready(paths))
  unlink(paths$peptidase_dmnd)
  expect_false(DNMB:::.dnmb_dbcan_standalone_ready(paths))
  file.create(paths$peptidase_dmnd)
  unlink(paths$pul_dir, recursive = TRUE)
  expect_false(DNMB:::.dnmb_dbcan_standalone_ready(paths))
})

test_that("dbCAN recognizes the v5 synteny_pdf directory", {
  root <- tempfile("dbcan-synteny-")
  dir.create(file.path(root, "synteny_pdf"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  expect_identical(
    DNMB:::.dnmb_dbcan_synteny_output_dir(root),
    file.path(root, "synteny_pdf")
  )
})

test_that("dbCAN retries HMMER when standalone execution fails", {
  root <- tempfile("dbcan-fallback-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  fallback_called <- FALSE
  empty_result <- function(ok, component) list(
    ok = ok,
    status = DNMB:::.dnmb_dbcan_status_row(component, if (ok) "ok" else "failed", "fixture"),
    files = list(), manifest = NULL,
    raw_hits = DNMB:::.dnmb_dbcan_empty_hits(),
    hits = DNMB:::.dnmb_module_empty_optional_long_table(),
    overview = DNMB:::.dnmb_dbcan_empty_overview(),
    cgc_genes = tibble::tibble(), substrate = tibble::tibble()
  )
  local_mocked_bindings(
    dnmb_detect_binary = function(...) list(found = TRUE, message = "fixture"),
    dnmb_dbcan_run_standalone = function(...) empty_result(FALSE, "standalone_fixture"),
    dnmb_dbcan_run_hmmsearch = function(...) {
      fallback_called <<- TRUE
      empty_result(TRUE, "hmm_fixture")
    },
    .package = "DNMB"
  )
  genes <- make_dbcan_genes()[1, , drop = FALSE]
  result <- DNMB:::dnmb_run_dbcan_module(genes, root, install = FALSE)

  expect_true(result$ok)
  expect_true(fallback_called)
  expect_true(all(c("standalone_fixture", "run_dbcan_fallback", "hmm_fixture") %in% result$status$component))
})

test_that("comparative dbCAN collector accepts legacy overview headers", {
  root <- tempfile("dbcan-comparative-")
  genome <- file.path(root, "Genome_A")
  run_dir <- file.path(genome, "dnmb_module_dbcan", "run_dbcan")
  dir.create(run_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(c("LOCUS       FIXTURE", "//"), file.path(genome, "fixture.gbff"))
  fixture <- data.frame(
    "Protein ID" = c("p1", "p2", "p3"),
    HMMER = c("GH5(1-100)", "-", "-"),
    "dbCAN-sub" = c("-", "CBM2_e1", "-"),
    diamond = c("-", "-", "GT4"),
    "Recommended Results" = c("-", "-", "-"),
    check.names = FALSE
  )
  utils::write.table(
    fixture, file.path(run_dir, "overview.tsv"), sep = "\t",
    row.names = FALSE, quote = FALSE
  )
  collected <- DNMB:::.dnmb_comparative_collect_dbcan(
    root,
    marker_rel = file.path("dnmb_module_dbcan", "run_dbcan", "overview.tsv"),
    level = "family",
    verbose = FALSE
  )
  expect_setequal(collected$systems$subtype, c("GH5", "CBM2", "GT4"))
})

test_that("dbCAN zero-protein completion sentinel is cacheable", {
  root <- tempfile("dbcan-empty-cache-")
  module_dir <- file.path(root, "dnmb_module_dbcan")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  file.create(file.path(module_dir, "dbcan_query_proteins.faa"))
  writeLines("dbcan_gene_row\tdbcan_query_id", file.path(module_dir, "dbcan_gene_id_map.tsv"))
  DNMB:::.dnmb_module_write_completion_sentinel("dbCAN", wd = root, n_hits = 0L)

  expect_true(DNMB:::.dnmb_module_stage_artifacts_exist("dbCAN", wd = root))
})

test_that("dbCAN empty input removes stale standalone outputs", {
  root <- tempfile("dbcan-stale-")
  stale_dir <- file.path(root, "run_dbcan")
  dir.create(stale_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("stale", file.path(stale_dir, "overview.tsv"))
  genes <- make_dbcan_genes()[4, , drop = FALSE]

  result <- DNMB:::dnmb_run_dbcan_module(genes, root, install = FALSE)
  expect_true(result$ok)
  expect_false(dir.exists(stale_dir))
})

test_that("GapMind carbon step status preserves best-path primary and positive secondary loci", {
  root <- tempfile("gapmind-carbon-steps-")
  module_dir <- file.path(root, "dnmb_module_gapmindcarbon")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  steps <- data.frame(
    pathway = c("arabinose", "galactose", "xylose"),
    step = c("araE", "galP", "xylT"),
    onBestPath = c(1, 1, 0),
    score = c(0, 2, 2),
    locusId = c("LOC_PRIMARY_ZERO", "LOC_PRIMARY_HIGH", "LOC_OFF_PATH"),
    score2 = c(2, 0, 2),
    locusId2 = c("LOC_SECONDARY_HIGH", "LOC_SECONDARY_ZERO", "LOC_OFF_PATH_2"),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    steps,
    file.path(module_dir, "aa.sum.steps"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  status <- DNMB:::.dnmb_gapmind_carbon_step_status(data.frame(), output_dir = root)

  expect_identical(
    names(status),
    c("pathway_id", "step_id", "confidence", "locus_tag")
  )
  expect_setequal(
    paste(status$pathway_id, status$step_id, status$locus_tag, status$confidence, sep = "::"),
    c(
      "arabinose::araE::LOC_PRIMARY_ZERO::none",
      "arabinose::araE::LOC_SECONDARY_HIGH::high",
      "galactose::galP::LOC_PRIMARY_HIGH::high"
    )
  )
  expect_false(any(status$locus_tag %in% c("LOC_SECONDARY_ZERO", "LOC_OFF_PATH", "LOC_OFF_PATH_2")))
})

test_that("GapMind carbon workbook fallback keeps distinct loci per step", {
  workbook <- data.frame(
    GapMindCarbon_pathway_id = rep("arabinose", 4),
    GapMindCarbon_step_id = rep("araE", 4),
    GapMindCarbon_confidence = c("low", "high", "medium", "none"),
    GapMindCarbon_on_best_path = c(TRUE, TRUE, 1, FALSE),
    locus_tag = c("LOC_A", "LOC_A", "LOC_B", "LOC_IGNORED"),
    stringsAsFactors = FALSE
  )

  status <- DNMB:::.dnmb_gapmind_carbon_step_status(workbook)

  expect_identical(
    status[, c("locus_tag", "confidence")],
    data.frame(
      locus_tag = c("LOC_A", "LOC_B"),
      confidence = c("high", "medium"),
      stringsAsFactors = FALSE
    )
  )
})
