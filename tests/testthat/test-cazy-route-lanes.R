test_that("shared horizontal route segments receive deterministic bounded lanes", {
  routes <- list(
    alpha = data.frame(x = c(0, 0, 4, 4), y = c(2, 1, 1, 0)),
    beta = data.frame(x = c(1, 1, 5, 5), y = c(2, 1, 1, 0)),
    gamma = data.frame(x = c(2, 2, 6, 6), y = c(2, 1.02, 1.02, 0))
  )
  metadata <- data.frame(
    route_id = c("alpha", "beta", "gamma"),
    priority = c(3, 2, 1),
    group = c("hexose", "hexose", "pentose"),
    stringsAsFactors = FALSE
  )

  separated <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes, metadata,
    lane_step = 0.06, max_offset = 0.20,
    y_tolerance = 0.04, min_overlap = 0.20,
    connection_mode = "move_internal"
  )
  assignments <- attr(separated, "lane_assignments")
  lane_by_id <- stats::setNames(assignments$lane_y, assignments$route_id)

  shuffled <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes[c(3, 1, 2)], metadata[c(3, 1, 2), ],
    lane_step = 0.06, max_offset = 0.20,
    y_tolerance = 0.04, min_overlap = 0.20,
    connection_mode = "move_internal"
  )
  shuffled_assignments <- attr(shuffled, "lane_assignments")
  shuffled_by_id <- stats::setNames(
    shuffled_assignments$lane_y, shuffled_assignments$route_id
  )

  expect_equal(lane_by_id[sort(names(lane_by_id))],
               shuffled_by_id[sort(names(shuffled_by_id))])
  expect_equal(length(unique(round(assignments$lane_y, 8))), 3L)
  expect_true(all(abs(assignments$offset) <= 0.20 + 1e-10))
})

test_that("route lane transitions preserve endpoints and original vertices", {
  routes <- list(
    first = data.frame(x = c(0, 0, 4, 4), y = c(3, 2, 2, 0)),
    second = data.frame(x = c(1, 1, 5, 5), y = c(3, 2, 2, 0))
  )
  metadata <- data.frame(
    route_id = c("first", "second"), priority = c(2, 1), group = "entry"
  )
  separated <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes, metadata,
    lane_step = 0.06, max_offset = 0.12,
    connection_mode = "transition"
  )

  for (i in seq_along(routes)) {
    original <- routes[[i]]
    adjusted <- separated[[i]]
    expect_equal(unlist(adjusted[1, c("x", "y")]),
                 unlist(original[1, c("x", "y")]))
    expect_equal(
      unlist(adjusted[nrow(adjusted), c("x", "y")]),
      unlist(original[nrow(original), c("x", "y")])
    )
    original_keys <- paste(original$x, original$y, sep = ":")
    adjusted_keys <- paste(adjusted$x, adjusted$y, sep = ":")
    expect_true(all(original_keys %in% adjusted_keys))
  }
})

test_that("move-internal release lanes retain origin and target with rounded elbows", {
  routes <- list(
    glucose = data.frame(x = c(2, 2, 8, 8), y = c(7, 6.8, 6.8, 4)),
    galactose = data.frame(x = c(3, 3, 9, 9), y = c(7, 6.8, 6.8, 4))
  )
  metadata <- data.frame(
    route_id = c("release:glucose", "release:galactose"),
    priority = c(1, 1), group = c("release", "release")
  )
  separated <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes, metadata,
    lane_step = 0.07, max_offset = 0.14,
    y_tolerance = 0.04, min_overlap = 0.20,
    connection_mode = "move_internal", round_radius = 0.04
  )
  assignments <- attr(separated, "lane_assignments")

  expect_equal(length(unique(round(assignments$lane_y, 8))), 2L)
  for (i in seq_along(routes)) {
    expect_equal(unlist(separated[[i]][1, c("x", "y")]),
                 unlist(routes[[i]][1, c("x", "y")]))
    expect_equal(unlist(separated[[i]][nrow(separated[[i]]), c("x", "y")]),
                 unlist(routes[[i]][nrow(routes[[i]]), c("x", "y")]))
    expect_gt(nrow(separated[[i]]), nrow(routes[[i]]))
  }
})

test_that("continuity colour stops at the first central-carbon merge", {
  expect_identical(
    DNMB:::.dnmb_cct_route_to_first_core_merge(
      c("Glucose", "Glc-6-P", "Fru-6-P", "Pyruvate")
    ),
    c("Glucose", "Glc-6-P")
  )
  expect_identical(
    DNMB:::.dnmb_cct_route_to_first_core_merge(
      c("Fuculose", "Fuculose-1-P", "DHAP", "GA3P")
    ),
    c("Fuculose", "Fuculose-1-P", "DHAP")
  )
  expect_identical(
    DNMB:::.dnmb_cct_route_to_first_core_merge(
      c("Deoxyribose-5-P", "GA3P", "1,3-BPG")
    ),
    c("Deoxyribose-5-P", "GA3P")
  )
  expect_null(
    DNMB:::.dnmb_cct_route_to_first_core_merge(c("R-5-P", "S-7-P"))
  )
})

test_that("branched carbon sources retain every initial intracellular product", {
  node_ids <- DNMB:::.dnmb_cct_sugar_nodes()$id

  expect_setequal(
    DNMB:::.dnmb_cct_initial_entry_targets("Lactose", node_ids),
    c("Glucose", "Gal-1-P")
  )
  expect_setequal(
    DNMB:::.dnmb_cct_initial_entry_targets("Sucrose", node_ids),
    c("Glucose", "Fru-6-P")
  )
  expect_identical(
    DNMB:::.dnmb_cct_initial_entry_targets("Maltose", node_ids),
    "Glucose"
  )
})
