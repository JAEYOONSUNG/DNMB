.cazy_route_directions <- function(points, tolerance = 1e-8) {
  if (is.null(points) || nrow(points) < 2L) return(character())
  dx <- diff(points$x)
  dy <- diff(points$y)
  out <- ifelse(
    abs(dx) <= tolerance & abs(dy) <= tolerance, "Z",
    ifelse(abs(dy) <= tolerance, "H", ifelse(abs(dx) <= tolerance, "V", "D"))
  )
  out[out != "Z"]
}

test_that("selected transporters resolve to visible extracellular anchors", {
  carbon_src <- data.frame(
    id = c("Maltose", "Galactose", "Sucrose", "Mannitol"),
    x = c(2, 4, 6, 8), y = rep(8, 4),
    sugar_type = c("glucose", "galactose", "glucose", "mannose")
  )
  anchors <- DNMB:::.dnmb_cct_transporter_source_anchors(
    cs_ids = carbon_src$id, carbon_src = carbon_src,
    mono_x_map = c(galactose = 4.5), y_mono = 9.5,
    extra_center_map = c(Maltose = 2.4),
    extra_y_map = c(Maltose = 10.5)
  )
  anchors <- anchors[match(carbon_src$id, anchors$cs_id), , drop = FALSE]

  expect_identical(anchors$source_kind,
                   c("complex_chain", "monosaccharide",
                     "explicit_source", "explicit_source"))
  expect_equal(anchors$source_x, c(2.4, 4.5, 6, 8))
  expect_equal(anchors$source_y, c(10.5, 9.5, 9.5, 9.5))
  expect_identical(anchors$create_glyph, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(anchors$cyto_x, carbon_src$x)
  expect_false("Ribose" %in% anchors$cs_id)
})

test_that("transporter pathways never substitute chemically distinct sugars", {
  pathway_map <- DNMB:::.dnmb_cct_transporter_pathway_map()
  pul_map <- DNMB:::.dnmb_cct_pul_substrate_to_cs()

  expect_identical(unname(pathway_map["glucose"]), "Glucose")
  expect_identical(unname(pathway_map["maltose"]), "Maltose")
  expect_true(is.na(unname(pathway_map["deoxyribonate"])))
  expect_identical(unname(pul_map["glucose"]), "Glucose")
  expect_true(is.na(unname(pul_map["pectin"])))
})

test_that("explicit uptake glyphs are packed away from visible monomers", {
  carbon_src <- data.frame(
    id = c("Sucrose", "Mannitol"),
    x = c(4, 4.1), y = c(8, 8),
    sugar_type = c("fructose", "mannose"),
    stringsAsFactors = FALSE
  )
  anchors <- DNMB:::.dnmb_cct_transporter_source_anchors(
    cs_ids = carbon_src$id,
    carbon_src = carbon_src,
    mono_x_map = c(glucose = 4, fructose = 5),
    y_mono = 9.5
  )

  expect_true(all(anchors$create_glyph))
  expect_gte(min(abs(outer(anchors$source_x, c(4, 5), "-"))), 0.5 - 1e-10)
  expect_gte(abs(diff(sort(anchors$source_x))), 0.5 - 1e-10)
})

test_that("exterior bus paths use visible source coordinates and one rounded elbow pass", {
  memberships <- data.frame(
    cs_id = c("Maltose", "Maltose"),
    lane_x = c(99, 99),
    source_x = c(2.4, 2.4), source_y = c(10.5, 10.5),
    tx_draw = c(4, 5), ty_draw = c(8.5, 8.5),
    confidence = c("high", "high"),
    locus_tag = c("T_A", "T_B"),
    color = c("#123456", "#123456"), stringsAsFactors = FALSE
  )
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(
    memberships, y_memb = 8.5, corner_radius = 0.08
  )

  expect_equal(length(layout$exterior_routes), 2L)
  expect_true(all(layout$memberships$lane_x == 2.4))
  for (route in layout$exterior_routes) {
    expect_equal(unlist(route$points[1, c("x", "y")]),
                 c(x = 2.4, y = 10.5))
    expect_equal(unlist(route$points[nrow(route$points), c("x", "y")]),
                 c(x = route$target_x, y = route$target_y))
    directions <- .cazy_route_directions(route$points)
    adjacent <- if (length(directions) > 1L) {
      paste(head(directions, -1L), tail(directions, -1L), sep = "")
    } else character()
    expect_false(any(adjacent %in% c("HV", "VH")))
    expect_true("D" %in% directions)
  }
})

test_that("aligned and empty transporter layouts never create orphan guides", {
  empty <- DNMB:::.dnmb_cct_transporter_bus_layout(NULL)
  expect_length(empty$exterior_routes, 0L)

  aligned <- data.frame(
    cs_id = "Ribose", lane_x = 3, source_x = 3, source_y = 9.5,
    tx_draw = 3, ty_draw = 8.5, confidence = "high",
    locus_tag = "RBS_T", color = "#008000", stringsAsFactors = FALSE
  )
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(aligned)
  expect_length(layout$exterior_routes, 1L)
  points <- layout$exterior_routes[[1L]]$points
  expect_equal(unlist(points[1L, c("x", "y")]), c(x = 3, y = 9.5))
  expect_equal(
    unlist(points[nrow(points), c("x", "y")]),
    c(x = 3, y = 8.61)
  )
  expect_true(all(abs(diff(points$x)) <= 1e-10))
  expect_false("D" %in% .cazy_route_directions(points))
})

test_that("selected transporter paths continue from membrane into the cytoplasmic source", {
  memberships <- data.frame(
    cs_id = c("Sucrose", "Mannitol"),
    locus_tag = c("SCR_T", "MTL_T"),
    tx_draw = c(5.4, 8.3), ty_draw = c(8.5, 8.62),
    cyto_x = c(6, 8), cyto_y = c(8, 8),
    confidence = c("high", "medium"),
    color = c("#00A651", "#3182BD"), stringsAsFactors = FALSE
  )
  routes <- DNMB:::.dnmb_cct_transporter_interior_routes(
    memberships, y_memb = 8.5, corner_radius = 0.08
  )

  expect_equal(length(routes), 2L)
  for (i in seq_along(routes)) {
    points <- routes[[i]]$points
    row <- memberships[memberships$locus_tag == routes[[i]]$target_key, ]
    expect_equal(unlist(points[1, c("x", "y")]),
                 c(x = row$tx_draw, y = row$ty_draw - 0.11))
    expect_equal(unlist(points[nrow(points), c("x", "y")]),
                 c(x = row$cyto_x, y = row$cyto_y + 0.11))
    directions <- .cazy_route_directions(points)
    adjacent <- if (length(directions) > 1L) {
      paste(head(directions, -1L), tail(directions, -1L), sep = "")
    } else character()
    expect_false(any(adjacent %in% c("HV", "VH")))
    expect_true("D" %in% directions)
  }

  expect_length(DNMB:::.dnmb_cct_transporter_interior_routes(NULL), 0L)
  expect_length(DNMB:::.dnmb_cct_transporter_interior_routes(
    memberships[0, , drop = FALSE]
  ), 0L)
})

test_that("every drawn membership has exterior and interior halves at the same T-ID", {
  memberships <- data.frame(
    cs_id = c("Glucose", "Glucose", "Maltose", "Mannitol"),
    locus_tag = c("T_SHARED", "T_MEDIUM", "T_SHARED", "T_MTL"),
    source_x = c(2, 2, 4, 6), source_y = c(9.5, 9.5, 10.5, 9.5),
    lane_x = c(2, 2, 4, 6),
    tx_draw = c(3, 3.5, 3, 6.5), ty_draw = rep(8.5, 4),
    cyto_x = c(2.5, 2.5, 4.5, 6), cyto_y = rep(8, 4),
    confidence = c("high", "medium", "high", "medium"),
    color = c("#004A7C", "#004A7C", "#004A7C", "#1B5A8A"),
    stringsAsFactors = FALSE
  )
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(
    memberships, y_memb = 8.5, suppress_redundant_medium = TRUE
  )
  drawn <- layout$memberships[layout$memberships$draw, , drop = FALSE]
  draw_keys <- paste(drawn$cs_id, drawn$target_key, sep = "\r")
  membership_keys <- paste(memberships$cs_id, memberships$locus_tag, sep = "\r")
  selected <- memberships[membership_keys %in% draw_keys, , drop = FALSE]
  interior <- DNMB:::.dnmb_cct_transporter_interior_routes(selected, y_memb = 8.5)

  exterior_keys <- vapply(layout$exterior_routes, function(route) {
    paste(route$cs_id, route$target_key, sep = "\r")
  }, character(1))
  interior_keys <- vapply(interior, function(route) {
    paste(route$cs_id, route$target_key, sep = "\r")
  }, character(1))
  expect_setequal(exterior_keys, interior_keys)
  expect_false(paste("Glucose", "T_MEDIUM", sep = "\r") %in% exterior_keys)
  expect_true(paste("Glucose", "T_SHARED", sep = "\r") %in% exterior_keys)
  expect_true(paste("Maltose", "T_SHARED", sep = "\r") %in% exterior_keys)
})

test_that("suppressed transporter entities remain ledger-only", {
  memberships <- data.frame(
    cs_id = c("Glucose", "Glucose"),
    lane_x = c(2, 2), source_x = c(2, 2), source_y = c(9.7, 9.7),
    tx_draw = c(3, 4), ty_draw = c(8.5, 8.5),
    confidence = c("high", "medium"),
    locus_tag = c("T_HIGH", "T_MEDIUM"),
    color = c("#004A7C", "#004A7C"), stringsAsFactors = FALSE
  )
  entities <- data.frame(
    locus_tag = c("T_HIGH", "T_MEDIUM"),
    anchor_id = c("T01", "T02"), stringsAsFactors = FALSE
  )
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(memberships)
  drawn <- DNMB:::.dnmb_cct_drawn_transporter_entities(
    entities, layout$memberships
  )

  expect_identical(drawn$locus_tag, "T_HIGH")
  expect_true("T_MEDIUM" %in% entities$locus_tag)
})
