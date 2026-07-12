test_that("core metabolism requires corroborated annotations and rejects common confounders", {
  core <- data.frame(
    locus_tag = c(
      "core_cs", "core_acn", "core_icd", "core_odh_e1", "core_odh_e2",
      "core_suc_d", "core_suc_c", "core_sdh_a", "core_sdh_b", "core_sdh_c",
      "core_fum", "core_mdh", "core_ace_a", "core_ace_b",
      "noise_prpc", "noise_prpd", "noise_leub", "noise_prpb", "noise_bsha"
    ),
    gene = c(
      "gltA", "acnA", "icd", "sucA", "sucB",
      "sucD", "sucC", "sdhA", "sdhB", "sdhC",
      "fumB", "mdh", "aceA", "aceB",
      "prpC", "prpD", "leuB", "prpB", "bshA"
    ),
    product = c(
      "citrate synthase", "aconitate hydratase", "isocitrate dehydrogenase",
      "2-oxoglutarate dehydrogenase E1 component",
      "2-oxoglutarate dehydrogenase E2 dihydrolipoyl succinyltransferase",
      "succinate--CoA ligase subunit alpha", "succinate--CoA ligase subunit beta",
      "succinate dehydrogenase flavoprotein subunit",
      "succinate dehydrogenase iron-sulfur subunit",
      "succinate dehydrogenase cytochrome b subunit",
      "fumarate hydratase", "malate dehydrogenase", "isocitrate lyase",
      "malate synthase A", "citrate synthase", "2-methylcitrate dehydratase",
      "isocitrate/isopropylmalate family dehydrogenase",
      "methylisocitrate lyase", "N-acetylglucosaminyl L-malate synthase"
    ),
    EC_number = c(
      "2.3.3.1", "4.2.1.3", "1.1.1.42", "1.2.4.2", "2.3.1.61",
      "6.2.1.5", "6.2.1.5", "1.3.5.1", "1.3.5.1", "1.3.5.1",
      "4.2.1.2", "1.1.1.37", "4.1.3.1", "2.3.3.9",
      "2.3.3.1", "4.2.1.79", "1.1.1.41", "4.1.3.30", "2.4.1.251"
    ),
    EggNOG_Preferred_name = c(
      "gltA", "acnA", "icd", "sucA", "sucB", "sucD", "sucC",
      "sdhA", "sdhB", "sdhC", "fumB", "mdh", "aceA", "aceB",
      "prpC", "prpD", "leuB", "prpB", "bshA"
    ),
    EggNOG_EC = NA_character_,
    EggNOG_KEGG_ko = c(
      "ko:K01647", "ko:K01681", "ko:K00031", "ko:K00164", "ko:K00658",
      "ko:K01902", "ko:K01903", "ko:K00239", "ko:K00240", "ko:K00241",
      "ko:K01676", "ko:K00024", "ko:K01637", "ko:K01638",
      "ko:K01647,ko:K01659", "ko:K01720", "ko:K00030,ko:K00052",
      "ko:K03417", "ko:K00691"
    ),
    EggNOG_Description = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- DNMB:::.dnmb_cct_core_metabolism_evidence(core)

  core_steps <- result$steps[result$steps$pathway != "Glycolysis", , drop = FALSE]
  expect_true(all(core_steps$status == "active"))
  retained <- unlist(strsplit(core_steps$locus_tags, "; ", fixed = TRUE))
  expect_false(any(grepl("^noise_", retained)))
  expect_setequal(
    core_steps$reaction_id,
    sprintf("C%02d", 1:10)
  )

  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  labels <- DNMB:::.dnmb_cct_core_direct_labels(
    core_steps, result$display_hits, nodes
  )
  expect_setequal(labels$reaction_id, sprintf("C%02d", 1:10))
  expect_false(any(grepl("^C[0-9]{2}", labels$label)))
  expect_match(labels$label[labels$reaction_id == "C06"], "sdhA  core_sdh_a")
  expect_match(labels$label[labels$reaction_id == "C06"], "sdhC  core_sdh_c")
})

test_that("glycolysis reactions require specific corroborated annotations", {
  glycolysis <- data.frame(
    locus_tag = c(
      "gly_glk", "gly_pgi", "gly_pfk", "gly_fba", "gly_tpi",
      "gly_gap", "gly_pgk", "gly_gpm", "gly_eno", "gly_pyk",
      "noise_fruk", "noise_lacd", "noise_mtnw"
    ),
    gene = c(
      "glk", "pgi", "pfkA", "fbaB", "tpiA", "gapA", "pgk", "gpmI",
      "eno", "pyk", "fruK", "lacD", "mtnW"
    ),
    product = c(
      "glucokinase", "glucose-6-phosphate isomerase", "6-phosphofructokinase",
      "class I fructose-bisphosphate aldolase", "triose-phosphate isomerase",
      "glyceraldehyde-3-phosphate dehydrogenase", "phosphoglycerate kinase",
      "2,3-bisphosphoglycerate-independent phosphoglycerate mutase",
      "phosphopyruvate hydratase", "pyruvate kinase",
      "1-phosphofructokinase", "tagatose-bisphosphate aldolase",
      "2,3-diketo-5-methylthiopentyl-1-phosphate enolase"
    ),
    EC_number = c(
      "2.7.1.2", "5.3.1.9", "2.7.1.11", "4.1.2.13", "5.3.1.1",
      "1.2.1.12", "2.7.2.3", "5.4.2.12", "4.2.1.11", "2.7.1.40",
      "2.7.1.56", "4.1.2.40", "5.3.2.5"
    ),
    EggNOG_Preferred_name = c(
      "glk", "pgi", "pfkA", "fbaB", "tpiA", "gapA", "pgk", "gpmI",
      "eno", "pyk", "fruK", "lacD", "mtnW"
    ),
    EggNOG_EC = NA_character_,
    EggNOG_KEGG_ko = c(
      "ko:K00845", "ko:K01810", "ko:K00850", "ko:K11645", "ko:K01803",
      "ko:K00134", "ko:K00927", "ko:K15633", "ko:K01689", "ko:K00873",
      "ko:K00882", "ko:K01622", "ko:K01601"
    ),
    EggNOG_Description = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- DNMB:::.dnmb_cct_core_metabolism_evidence(glycolysis)
  steps <- result$steps[result$steps$pathway == "Glycolysis", , drop = FALSE]
  expect_setequal(steps$reaction_id, sprintf("C%02d", 11:20))
  expect_true(all(steps$status == "active"))
  retained <- unlist(strsplit(steps$locus_tags, "; ", fixed = TRUE))
  expect_false(any(grepl("^noise_", retained)))

  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  labels <- DNMB:::.dnmb_cct_core_direct_labels(
    steps, result$display_hits, nodes
  )
  expect_setequal(labels$reaction_id, sprintf("C%02d", 11:20))
  expect_true(all(grepl("gly_", labels$label, fixed = TRUE)))
  expect_false(any(grepl("^C[0-9]{2}", labels$label)))
})

test_that("the map omits the duplicate central-carbon summary shelf", {
  plot_body <- paste(
    deparse(body(DNMB:::.dnmb_plot_cazy_carbon_3zone_v2)),
    collapse = "\n"
  )

  expect_false(grepl("Central carbon metabolism", plot_body, fixed = TRUE))
  expect_false(grepl(
    "Conservative genome-annotation support", plot_body, fixed = TRUE
  ))
  expect_true(grepl("core_direct_label_df", plot_body, fixed = TRUE))
  expect_true(grepl("core_ledger_df", plot_body, fixed = TRUE))
  expect_true(grepl("Full functional evidence", plot_body, fixed = TRUE))
})

test_that("multi-subunit core reactions remain partial when a required component is absent", {
  partial <- data.frame(
    locus_tag = c("odh_e1", "sdh_a", "sdh_b", "frd_a"),
    gene = c("sucA", "sdhA", "sdhB", "frdA"),
    product = c(
      "2-oxoglutarate dehydrogenase E1 component",
      "succinate dehydrogenase flavoprotein subunit",
      "succinate dehydrogenase iron-sulfur subunit",
      "fumarate reductase flavoprotein subunit"
    ),
    EC_number = c("1.2.4.2", "1.3.5.1", "1.3.5.1", "1.3.5.4"),
    EggNOG_Preferred_name = c("sucA", "sdhA", "sdhB", "frdA"),
    EggNOG_EC = NA_character_,
    EggNOG_KEGG_ko = c("ko:K00164", "ko:K00239", "ko:K00240", "ko:K00244"),
    EggNOG_Description = NA_character_,
    stringsAsFactors = FALSE
  )

  steps <- DNMB:::.dnmb_cct_core_metabolism_evidence(partial)$steps
  expect_identical(steps$status[steps$reaction_id == "C04"], "partial")
  expect_identical(steps$status[steps$reaction_id == "C06"], "partial")
})

test_that("Ribose and Fucose routes have distinct biochemical semantics", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  nodes <- DNMB:::.dnmb_cct_merge_entry_into_mainstream(nodes)
  nodes <- DNMB:::.dnmb_cct_resolve_grid_overlaps(nodes)
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)

  has_edge <- function(from, to) any(edges$from == from & edges$to == to)
  expect_true(has_edge("Ribose", "R-5-P"))
  expect_true(has_edge("Fucose", "Fuculose"))
  expect_true(has_edge("Fuculose", "Fuculose-1-P"))
  expect_false(has_edge("Ribose", "Fuculose"))
})

test_that("deoxyribose has a dedicated directed mini-branch", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)

  expect_true("Deoxyribose" %in% nodes$id)
  expect_true("Deoxyribose-5-P" %in% nodes$id)
  expect_identical(
    nodes$sugar_type[nodes$id == "Deoxyribose"],
    "deoxyribose"
  )
  has_edge <- function(from, to) {
    any(edges$from == from & edges$to == to)
  }
  expect_true(has_edge("Deoxyribose", "Deoxyribose-5-P"))
  expect_true(has_edge("Deoxyribose-5-P", "GA3P"))
  expect_true(has_edge("Deoxyribose-5-P", "Acetaldehyde"))
  expect_true(has_edge("Acetaldehyde", "Acetate"))
  expect_true(has_edge("Acetate", "Acetyl-CoA"))
})

test_that("cytoplasmic disaccharide glyphs preserve component identity", {
  components <- function(source_id) {
    DNMB:::.dnmb_cct_carbon_source_render_spec(source_id)$sugar_type
  }

  expect_identical(components("Maltose"), c("glucose", "glucose"))
  expect_identical(components("Cellobiose"), c("glucose", "glucose"))
  expect_identical(components("Trehalose"), c("glucose", "glucose"))
  expect_identical(components("Lactose"), c("galactose", "glucose"))
  expect_identical(components("Sucrose"), c("glucose", "fructose"))

  expect_identical(
    unique(DNMB:::.dnmb_cct_carbon_source_render_spec("Cellobiose")$bond),
    "beta"
  )
  expect_identical(
    unique(DNMB:::.dnmb_cct_carbon_source_render_spec("Maltose")$bond),
    "alpha"
  )
})

test_that("every cytoplasmic carbon source has a registered glyph spec", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  source_ids <- nodes$id[nodes$type == "carbon_source"]
  specs <- lapply(source_ids, DNMB:::.dnmb_cct_carbon_source_render_spec)

  expect_true(length(source_ids) > 0L)
  expect_true(all(vapply(specs, function(x) all(x$registered), logical(1))))
  expect_true(all(vapply(specs, function(x) nrow(x) >= 1L, logical(1))))
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("NAG"), "GlcNAc")
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("deoxyribose"), "dRib")
})

test_that("exact deoxyribose mappings are case-insensitive and directed", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  evidence <- data.frame(
    pathway_id = c(
      "DeOxYrIbOsE", "DEOXYRIBOSE", "deoxyribose",
      "Deoxyribose", "DEOXYRIBOSE"
    ),
    step_id = c("DeOk", "DEOC", "adh", "ACS", "not-a-reaction"),
    confidence = c("high", "high", "medium", "high", "high"),
    rank = c(3L, 3L, 2L, 3L, 3L),
    locus_tag = c("deoK_1", "deoC_1", "adh_1", "acs_1", "noise_1"),
    stringsAsFactors = FALSE
  )

  matches <- DNMB:::.dnmb_cct_exact_step_edge_matches(evidence, edges)
  expect_setequal(
    unique(tolower(matches$mapped_step_id)),
    c("deok", "deoc", "adh", "acs")
  )
  expect_equal(sum(tolower(matches$mapped_step_id) == "deoc"), 2L)
  expect_false("noise_1" %in% matches$locus_tag)
  expect_true(all(
    paste(matches$mapped_from, matches$mapped_to, sep = "->") %in%
      paste(edges$from, edges$to, sep = "->")
  ))

  reversed <- edges
  idx <- reversed$from == "Acetate" & reversed$to == "Acetyl-CoA"
  reversed$from[idx] <- "Acetyl-CoA"
  reversed$to[idx] <- "Acetate"
  reversed_matches <- DNMB:::.dnmb_cct_exact_step_edge_matches(
    evidence[evidence$step_id == "ACS", , drop = FALSE], reversed
  )
  expect_equal(nrow(reversed_matches), 0L)
})

test_that("exact labels aggregate loci and never use arbitrary edge fallback", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  evidence <- data.frame(
    pathway_id = c("DEOXYRIBOSE", "deoxyribose", "deoxyribose", "ribose"),
    step_id = c("deoC", "DEOC", "adh", "deoC"),
    confidence = c("medium", "high", "medium", "high"),
    rank = c(2L, 3L, 2L, 3L),
    locus_tag = c("deoC_a", "deoC_b", "adh_a", "unmapped"),
    stringsAsFactors = FALSE
  )

  labels <- DNMB:::.dnmb_cct_exact_step_labels(evidence, edges)
  deoc <- labels[tolower(labels$step_id) == "deoc", , drop = FALSE]
  adh <- labels[tolower(labels$step_id) == "adh", , drop = FALSE]
  expect_equal(nrow(deoc), 1L)
  expect_match(deoc$label, "deoC_a", fixed = TRUE)
  expect_match(deoc$label, "deoC_b", fixed = TRUE)
  expect_identical(deoc$color, "#007C83")
  expect_identical(adh$color, "#B7791F")
  expect_false(any(grepl("unmapped", labels$label, fixed = TRUE)))

  unknown <- evidence[1, , drop = FALSE]
  unknown$step_id <- "unmapped-step"
  expect_equal(
    nrow(DNMB:::.dnmb_cct_exact_step_labels(unknown, edges)),
    0L
  )
})

test_that("hub routing avoids unrelated metabolite nodes", {
  node_x <- c(
    "R-5-P" = 19, "Fuculose" = 22, "Fuculose-1-P" = 22,
    "DHAP" = 9.5
  )
  node_y <- c(
    "R-5-P" = 5.5, "Fuculose" = 7, "Fuculose-1-P" = 6.5,
    "DHAP" = 5.5
  )
  obstacle_x <- c(node_x, Ribose = 21.75, Fucose = 22.25)
  obstacle_y <- c(node_y, Ribose = 8, Fucose = 8)

  ribose_path <- DNMB:::.dnmb_cct_points_from_start_to_nodes(
    21.75, 8, "R-5-P", node_x, node_y,
    obstacle_x = obstacle_x, obstacle_y = obstacle_y
  )
  fucose_path <- DNMB:::.dnmb_cct_points_from_start_to_nodes(
    22.25, 8, c("Fuculose", "Fuculose-1-P", "DHAP"), node_x, node_y,
    obstacle_x = obstacle_x, obstacle_y = obstacle_y
  )

  expect_lt(
    DNMB:::.dnmb_cct_path_obstacle_score(ribose_path, c(22, 22), c(7, 6.5)),
    1000
  )
  expect_lt(
    DNMB:::.dnmb_cct_path_obstacle_score(fucose_path, 19, 5.5),
    1000
  )
})

test_that("carbon-source lanes retain their optimized quarter-grid positions", {
  nodes <- data.frame(
    id = c("Ribose", "Fucose", "Fuculose"),
    x = c(21.75, 22.25, 22),
    y = c(8, 8, 7),
    label = c("Ribose", "Fucose", "Fuculose"),
    type = c("carbon_source", "carbon_source", "entry_intermediate"),
    sugar_type = c("ribose", "fucose", "fucose"),
    stringsAsFactors = FALSE
  )

  resolved <- DNMB:::.dnmb_cct_resolve_grid_overlaps(nodes, step = 0.5)
  expect_equal(resolved$x[resolved$id == "Ribose"], 21.75)
  expect_equal(resolved$x[resolved$id == "Fucose"], 22.25)
})

test_that("extracellular Rha uses its own rhamnose lane and evidence state", {
  lane_map <- DNMB:::.dnmb_cct_monomer_lane_map()

  expect_identical(unname(lane_map["fucose"]), "fucose")
  expect_identical(unname(lane_map["rhamnose"]), "rhamnose")
  expect_false(identical(unname(lane_map["fucose"]), unname(lane_map["rhamnose"])))
})

test_that("monomer aliases are separated on deterministic half-cell lanes", {
  preferred <- c(
    ribose = 24.75,
    fucose = 25.25,
    rhamnose = 25.75,
    glca = 26.50,
    gala = 26.50
  )

  placed <- DNMB:::.dnmb_cct_separate_lanes(
    preferred, min_gap = 0.5, snap_step = 0.25
  )
  repeated <- DNMB:::.dnmb_cct_separate_lanes(
    preferred, min_gap = 0.5, snap_step = 0.25
  )

  expect_identical(names(placed), names(preferred))
  expect_identical(placed, repeated)
  expect_equal(anyDuplicated(placed), 0L)
  expect_gte(min(diff(sort(placed))), 0.5)
  expect_gte(abs(placed["ribose"] - placed["fucose"]), 0.5)
  expect_gte(abs(placed["fucose"] - placed["rhamnose"]), 0.5)
  expect_gte(abs(placed["glca"] - placed["gala"]), 0.5)
})

test_that("short pentose routes stay inside endpoint bounds without loop-back", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  nodes <- DNMB:::.dnmb_cct_merge_entry_into_mainstream(nodes)
  nodes <- DNMB:::.dnmb_cct_resolve_grid_overlaps(nodes)
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  pentose_pairs <- list(
    c("Ru-5-P", "R-5-P"),
    c("Ru-5-P", "Xu-5-P"),
    c("Xu-5-P", "S-7-P"),
    c("R-5-P", "S-7-P"),
    c("S-7-P", "E-4-P")
  )

  for (pair in pentose_pairs) {
    edge <- edges[edges$from == pair[1] & edges$to == pair[2], , drop = FALSE]
    expect_equal(nrow(edge), 1L)
    path <- DNMB:::.dnmb_cct_edge_points_from_row(edge, grid_step = 0.5)
    expect_gte(min(path$x), min(edge$x, edge$xend) - 1e-8)
    expect_lte(max(path$x), max(edge$x, edge$xend) + 1e-8)
    expect_gte(min(path$y), min(edge$y, edge$yend) - 1e-8)
    expect_lte(max(path$y), max(edge$y, edge$yend) + 1e-8)

    path_length <- sum(sqrt(diff(path$x)^2 + diff(path$y)^2))
    manhattan <- abs(edge$xend - edge$x) + abs(edge$yend - edge$y)
    expect_lte(path_length, manhattan + 1e-8)
  }
})
