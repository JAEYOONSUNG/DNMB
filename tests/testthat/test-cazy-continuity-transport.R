test_that("continuity overlays retain their first explicit intracellular node", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  route_groups <- list(
    starch = c("maltose", "trehalose"),
    xylan = c("xylose", "arabinose"),
    deoxyribose = "deoxyribose",
    ppp = c("xylose", "arabinose", "ribose")
  )
  entry_map <- DNMB:::.dnmb_cct_auto_entry_map()

  expected_starts <- c(
    maltose = "Glucose",
    xylose = "Xylulose",
    deoxyribose = "Deoxyribose-5-P"
  )
  for (pathway_id in names(expected_starts)) {
    route <- DNMB:::.dnmb_cct_best_continuity_nodes(
      pathway_id = pathway_id,
      cyto_nodes = nodes,
      cyto_edges = edges,
      route_group_members = route_groups,
      entry_map = entry_map,
      matched_steps = data.frame(),
      transporters = NULL
    )
    expect_gte(length(route), 2L)
    expect_identical(unname(route[1L]), unname(expected_starts[pathway_id]))
  }
})

test_that("an aligned transporter trunk keeps a multi-target bus connected", {
  memberships <- data.frame(
    cs_id = c("Deoxyribose", "Deoxyribose"),
    lane_x = c(3, 3),
    tx_draw = c(3, 5),
    ty_draw = c(8.5, 8.5),
    confidence = c("high", "high"),
    locus_tag = c("aligned_target", "offset_target"),
    color = c("#8A5879", "#8A5879"),
    stringsAsFactors = FALSE
  )

  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(
    memberships, y_memb = 8.5, source_offset = 0.76
  )
  group_id <- layout$groups$group[1L]
  aligned <- layout$trunks[
    layout$trunks$group == "target:aligned_target", , drop = FALSE
  ]

  expect_false(layout$groups$direct[1L])
  expect_false(any(layout$stems$group == group_id))
  expect_equal(nrow(aligned), 1L)
  expect_equal(aligned$y, 8.5 + 0.76, tolerance = 1e-10)
  expect_lte(min(aligned$y, aligned$yend), layout$groups$bus_y[1L])
  expect_gte(max(aligned$y, aligned$yend), layout$groups$bus_y[1L])
})

test_that("deoP is retained as deoxyribose uptake with a source glyph lane", {
  expect_true(DNMB:::.dnmb_cct_is_transport_like_step("deoP"))

  genbank <- data.frame(
    locus_tag = "DEOP_0001",
    gene = "deoP",
    GapMindCarbon_pathway_id = "deoxyribose",
    GapMindCarbon_step_id = "deoP",
    GapMindCarbon_confidence = "high",
    GapMindCarbon_step_score = 2,
    stringsAsFactors = FALSE
  )
  transporters <- DNMB:::.dnmb_cct_3zone_extract_transporters(genbank)
  expect_equal(nrow(transporters), 1L)
  expect_identical(transporters$pathway, "deoxyribose")
  expect_identical(transporters$step, "deoP")

  expect_identical(
    unname(DNMB:::.dnmb_cct_monomer_lane_map()["deoxyribose"]),
    "deoxyribose"
  )
  expect_true(
    "deoxyribose" %in% DNMB:::.dnmb_cct_extracellular_monomer_ids()
  )
  source_nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  source_nodes <- source_nodes[source_nodes$type == "carbon_source", , drop = FALSE]
  expect_true("Deoxyribose" %in% source_nodes$id)
})
