test_that("deoxyribose acetaldehyde degradation alternatives have directed edges", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  mapping <- DNMB:::.dnmb_cct_exact_step_edge_map()

  expected <- data.frame(
    step = c("adh", "acka", "pta", "acs", "ald-dh-coa"),
    from = c("Acetaldehyde", "Acetate", "Acetyl-P", "Acetate", "Acetaldehyde"),
    to = c("Acetate", "Acetyl-P", "Acetyl-CoA", "Acetyl-CoA", "Acetyl-CoA"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(expected))) {
    edge_present <- edges$from == expected$from[i] & edges$to == expected$to[i]
    map_present <- tolower(mapping$pathway_id) == "deoxyribose" &
      tolower(mapping$step_id) == expected$step[i] &
      mapping$from == expected$from[i] & mapping$to == expected$to[i]
    expect_true(any(edge_present), info = expected$step[i])
    expect_true(any(map_present), info = expected$step[i])
  }
})

test_that("alternative deoxyribose steps match exact edges case-insensitively", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  evidence <- data.frame(
    pathway_id = rep(" DeOxYrIbOsE ", 5L),
    step_id = c("ADH", "AckA", "PTA", "AcS", "ALD-DH-COA"),
    confidence = "high",
    rank = 3L,
    locus_tag = paste0("locus_", seq_len(5L)),
    stringsAsFactors = FALSE
  )

  matches <- DNMB:::.dnmb_cct_exact_step_edge_matches(evidence, edges)
  expect_setequal(
    tolower(matches$mapped_step_id),
    c("adh", "acka", "pta", "acs", "ald-dh-coa")
  )
  expect_equal(nrow(matches), 5L)
})

test_that("exact reaction anchors lie on the rendered cumulative-midpoint path", {
  nodes <- DNMB:::.dnmb_cct_sugar_nodes()
  edges <- DNMB:::.dnmb_cct_3zone_cyto_edges(nodes)
  evidence <- data.frame(
    pathway_id = "deoxyribose",
    step_id = c("ackA", "pta", "acs", "ald-dh-CoA"),
    confidence = "high",
    rank = 3L,
    locus_tag = paste0("anchor_", seq_len(4L)),
    stringsAsFactors = FALSE
  )
  labels <- DNMB:::.dnmb_cct_exact_step_labels(evidence, edges)
  mapping <- DNMB:::.dnmb_cct_exact_step_edge_map()

  point_path_distance <- function(x, y, path) {
    distances <- vapply(seq_len(nrow(path) - 1L), function(i) {
      ax <- path$x[i]
      ay <- path$y[i]
      dx <- path$x[i + 1L] - ax
      dy <- path$y[i + 1L] - ay
      denom <- dx^2 + dy^2
      t <- if (denom > 0) ((x - ax) * dx + (y - ay) * dy) / denom else 0
      t <- min(1, max(0, t))
      sqrt((x - (ax + t * dx))^2 + (y - (ay + t * dy))^2)
    }, numeric(1))
    min(distances)
  }

  for (i in seq_len(nrow(labels))) {
    step <- tolower(labels$step_id[i])
    map_row <- mapping[
      tolower(mapping$pathway_id) == "deoxyribose" &
        tolower(mapping$step_id) == step & mapping$label_anchor,
      , drop = FALSE
    ][1, , drop = FALSE]
    edge <- edges[
      edges$from == map_row$from & edges$to == map_row$to,
      , drop = FALSE
    ][1, , drop = FALSE]
    path <- DNMB:::.dnmb_cct_edge_points_from_row(edge, grid_step = 0.5)
    expected_midpoint <- DNMB:::.dnmb_cct_route_label_point(path, frac = 0.5)

    expect_equal(labels$x[i], expected_midpoint$x, tolerance = 1e-10)
    expect_equal(labels$y[i], expected_midpoint$y, tolerance = 1e-10)
    expect_lt(point_path_distance(labels$x[i], labels$y[i], path), 1e-10)
  }

  acs_label <- labels[tolower(labels$step_id) == "acs", , drop = FALSE]
  acs_edge <- edges[edges$from == "Acetate" & edges$to == "Acetyl-CoA", , drop = FALSE]
  endpoint_mean <- c(mean(c(acs_edge$x, acs_edge$xend)), mean(c(acs_edge$y, acs_edge$yend)))
  expect_gt(sqrt((acs_label$x - endpoint_mean[1])^2 + (acs_label$y - endpoint_mean[2])^2), 0.1)
})

test_that("step target lookup is case-insensitive for mixed-case GapMind IDs", {
  expect_identical(
    DNMB:::.dnmb_cct_step_target_nodes(" NAgA ", "NAG"),
    c("GlcN-6-P")
  )
  expect_identical(DNMB:::.dnmb_cct_step_target_nodes("GalK", "galactose"), c("Gal-1-P"))
  expect_identical(DNMB:::.dnmb_cct_step_target_nodes("XYLa", "xylose"), c("Xylulose", "Xu-5-P"))
  expect_identical(DNMB:::.dnmb_cct_step_target_nodes("FucI", "fucose"), c("Fuculose", "Fuculose-1-P", "DHAP"))
  expect_identical(DNMB:::.dnmb_cct_step_target_nodes("LacZ", "lactose"), c("Gal-1-P", "Glc-6-P"))
  expect_identical(DNMB:::.dnmb_cct_step_target_nodes("ACKA", "DEOXYRIBOSE"), c("Acetate", "Acetyl-P"))
  expect_identical(DNMB:::.dnmb_cct_step_target_nodes("ald-DH-CoA", "deoxyribose"), c("Acetaldehyde", "Acetyl-CoA"))
})
