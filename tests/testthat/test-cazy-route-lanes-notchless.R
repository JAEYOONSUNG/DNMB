.cazy_notchless_direction_classes <- function(points, tolerance = 1e-8) {
  if (is.null(points) || nrow(points) < 2L) return(character())
  dx <- diff(points$x)
  dy <- diff(points$y)
  direction <- ifelse(
    abs(dx) <= tolerance & abs(dy) <= tolerance, "Z",
    ifelse(
      abs(dy) <= tolerance, "H",
      ifelse(abs(dx) <= tolerance, "V", "D")
    )
  )
  direction[direction != "Z"]
}

.cazy_notchless_has_horizontal_jog <- function(points, tolerance = 1e-8) {
  direction <- .cazy_notchless_direction_classes(points, tolerance = tolerance)
  if (length(direction) < 3L) return(FALSE)
  any(
    direction[-c(length(direction) - 1L, length(direction))] == "H" &
      direction[-c(1L, length(direction))] == "V" &
      direction[-c(1L, 2L)] == "H"
  )
}

.cazy_notchless_fixture <- function() {
  routes <- list(
    main = data.frame(
      x = c(0, 0, 3, 3, 6, 6),
      y = c(4, 3, 3, 3, 3, 2)
    ),
    left = data.frame(
      x = c(-0.5, -0.5, 2.5, 2.5),
      y = c(4, 3, 3, 2)
    ),
    right = data.frame(
      x = c(3.5, 3.5, 6.5, 6.5),
      y = c(4, 3, 3, 2)
    )
  )
  metadata <- data.frame(
    route_id = c("mid", "aaa", "zzz"),
    priority = c(1, 1, 1),
    group = "release",
    stringsAsFactors = FALSE
  )
  list(routes = routes, metadata = metadata)
}

.cazy_notchless_separate <- function(fixture, round_radius = 0) {
  DNMB:::.dnmb_cct_separate_route_lanes(
    fixture$routes,
    fixture$metadata,
    lane_step = 0.08,
    max_offset = 0.20,
    y_tolerance = 0.01,
    min_overlap = 0.20,
    connection_mode = "move_internal",
    round_radius = round_radius
  )
}

test_that("one route keeps one height across a continuous horizontal run", {
  fixture <- .cazy_notchless_fixture()
  separated <- .cazy_notchless_separate(fixture)
  assignments <- attr(separated, "lane_assignments")
  main_assignments <- assignments[assignments$route_id == "mid", , drop = FALSE]
  main_path <- separated[[which(fixture$metadata$route_id == "mid")]]

  expect_gt(nrow(main_assignments), 1L)
  expect_equal(length(unique(round(main_assignments$lane_y, 8))), 1L)
  expect_false(.cazy_notchless_has_horizontal_jog(main_path))

  dx <- diff(main_path$x)
  dy <- diff(main_path$y)
  horizontal_y <- main_path$y[-nrow(main_path)][
    abs(dx) > 1e-8 & abs(dy) <= 1e-8
  ]
  expect_equal(length(unique(round(horizontal_y, 8))), 1L)
  expect_equal(unlist(main_path[1L, c("x", "y")]), c(x = 0, y = 4))
  expect_equal(unlist(main_path[nrow(main_path), c("x", "y")]), c(x = 6, y = 2))
})

test_that("rounding is applied to the final offset lane without adding a notch", {
  fixture <- .cazy_notchless_fixture()
  unrounded <- .cazy_notchless_separate(fixture, round_radius = 0)
  rounded <- .cazy_notchless_separate(fixture, round_radius = 0.04)

  for (route_index in seq_along(unrounded)) {
    expected <- DNMB:::.dnmb_cct_round_orthogonal_route(
      unrounded[[route_index]], radius = 0.04
    )
    expect_equal(
      unname(as.matrix(rounded[[route_index]][, c("x", "y")])),
      unname(as.matrix(expected[, c("x", "y")]))
    )
  }

  main_path <- rounded[[which(fixture$metadata$route_id == "mid")]]
  expect_false(.cazy_notchless_has_horizontal_jog(main_path))
  expect_equal(unlist(main_path[1L, c("x", "y")]), c(x = 0, y = 4))
  expect_equal(unlist(main_path[nrow(main_path), c("x", "y")]), c(x = 6, y = 2))
})

test_that("rounded lane corners use symmetric incoming and outgoing tangents", {
  fixture <- .cazy_notchless_fixture()
  unrounded <- .cazy_notchless_separate(fixture, round_radius = 0)[[1L]]
  rounded <- .cazy_notchless_separate(fixture, round_radius = 0.04)[[1L]]
  radius <- 0.04

  checked_corners <- 0L
  for (corner_index in 2:(nrow(unrounded) - 1L)) {
    previous <- unlist(unrounded[corner_index - 1L, c("x", "y")])
    corner <- unlist(unrounded[corner_index, c("x", "y")])
    following <- unlist(unrounded[corner_index + 1L, c("x", "y")])
    incoming <- corner - previous
    outgoing <- following - corner
    incoming_length <- sqrt(sum(incoming^2))
    outgoing_length <- sqrt(sum(outgoing^2))
    incoming_unit <- incoming / incoming_length
    outgoing_unit <- outgoing / outgoing_length
    if (abs(sum(incoming_unit * outgoing_unit)) > 1e-8) next

    tangent_length <- min(radius, incoming_length * 0.30, outgoing_length * 0.30)
    expected_incoming <- corner - incoming_unit * tangent_length
    expected_outgoing <- corner + outgoing_unit * tangent_length
    incoming_distance <-
      (rounded$x - expected_incoming[[1L]])^2 +
      (rounded$y - expected_incoming[[2L]])^2
    outgoing_distance <-
      (rounded$x - expected_outgoing[[1L]])^2 +
      (rounded$y - expected_outgoing[[2L]])^2
    actual_incoming <- unlist(
      rounded[which.min(incoming_distance), c("x", "y")]
    )
    actual_outgoing <- unlist(
      rounded[which.min(outgoing_distance), c("x", "y")]
    )

    expect_equal(actual_incoming, expected_incoming, tolerance = 1e-10)
    expect_equal(actual_outgoing, expected_outgoing, tolerance = 1e-10)
    expect_equal(
      sqrt(sum((actual_incoming - corner)^2)),
      sqrt(sum((actual_outgoing - corner)^2)),
      tolerance = 1e-10
    )
    checked_corners <- checked_corners + 1L
  }
  expect_equal(checked_corners, 2L)
})

test_that("lane separation retains endpoints and unrelated route anchors", {
  routes <- list(
    anchored = data.frame(
      x = c(0, 0, 1, 1, 6, 6, 7, 7),
      y = c(5, 4, 4, 3, 3, 2, 2, 1)
    ),
    overlap = data.frame(
      x = c(2, 2, 5, 5),
      y = c(4, 3, 3, 2)
    )
  )
  metadata <- data.frame(
    route_id = c("anchored", "overlap"),
    priority = c(2, 1),
    group = "entry",
    stringsAsFactors = FALSE
  )
  separated <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes, metadata,
    lane_step = 0.08, max_offset = 0.16,
    connection_mode = "move_internal"
  )
  path <- separated[[1L]]
  preserved <- data.frame(
    x = c(0, 0, 1, 6, 7, 7),
    y = c(5, 4, 4, 2, 2, 1)
  )
  path_keys <- paste(path$x, path$y, sep = ":")
  anchor_keys <- paste(preserved$x, preserved$y, sep = ":")

  expect_true(all(anchor_keys %in% path_keys))
  expect_equal(unlist(path[1L, c("x", "y")]), c(x = 0, y = 5))
  expect_equal(unlist(path[nrow(path), c("x", "y")]), c(x = 7, y = 1))
})

test_that("short and same-axis routes never become diagonal", {
  routes <- list(
    point = data.frame(x = 1, y = 1),
    vertical = data.frame(x = c(2, 2, 2), y = c(4, 3, 2)),
    flat_a = data.frame(x = c(0, 2, 4), y = c(3, 3, 3)),
    flat_b = data.frame(x = c(0, 2, 4), y = c(3, 3, 3))
  )
  metadata <- data.frame(
    route_id = names(routes),
    priority = c(1, 1, 2, 1),
    group = "edge-case",
    stringsAsFactors = FALSE
  )
  separated <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes, metadata,
    lane_step = 0.08, max_offset = 0.16,
    connection_mode = "move_internal"
  )

  for (route_index in seq_along(routes)) {
    original <- routes[[route_index]]
    adjusted <- separated[[route_index]]
    expect_equal(
      unlist(adjusted[1L, c("x", "y")]),
      unlist(original[1L, c("x", "y")])
    )
    expect_equal(
      unlist(adjusted[nrow(adjusted), c("x", "y")]),
      unlist(original[nrow(original), c("x", "y")])
    )
    expect_false("D" %in% .cazy_notchless_direction_classes(adjusted))
    expect_false(.cazy_notchless_has_horizontal_jog(adjusted))
  }
})

test_that("notchless lane geometry is deterministic after route shuffling", {
  fixture <- .cazy_notchless_fixture()
  baseline <- .cazy_notchless_separate(fixture, round_radius = 0.04)
  baseline_by_id <- stats::setNames(baseline, fixture$metadata$route_id)

  order_shuffled <- c(3L, 1L, 2L)
  shuffled_fixture <- list(
    routes = fixture$routes[order_shuffled],
    metadata = fixture$metadata[order_shuffled, , drop = FALSE]
  )
  shuffled <- .cazy_notchless_separate(shuffled_fixture, round_radius = 0.04)
  shuffled_by_id <- stats::setNames(
    shuffled, shuffled_fixture$metadata$route_id
  )

  for (route_id in sort(fixture$metadata$route_id)) {
    expect_equal(
      unname(as.matrix(baseline_by_id[[route_id]][, c("x", "y")])),
      unname(as.matrix(shuffled_by_id[[route_id]][, c("x", "y")]))
    )
  }
})

test_that("hub-entry fallback stays orthogonal until final lane rounding", {
  raw_path <- DNMB:::.dnmb_cct_hub_entry_path(
    hub_x = 1, hub_y = 5,
    target_x = 5, target_y = 1,
    lane_rank = 2L, target_count = 3L,
    grid_step = 0.5,
    rounded = FALSE
  )
  direction <- .cazy_notchless_direction_classes(raw_path)

  expect_true(nrow(raw_path) >= 4L)
  expect_false("D" %in% direction)

  routes <- list(first = raw_path, second = raw_path)
  metadata <- data.frame(
    route_id = c("hub-fallback-a", "hub-fallback-b"),
    priority = c(2, 1), group = "hub",
    stringsAsFactors = FALSE
  )
  unrounded <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes = routes,
    metadata = metadata,
    lane_step = 0.08,
    max_offset = 0.16,
    round_radius = 0
  )
  rounded <- DNMB:::.dnmb_cct_separate_route_lanes(
    routes = routes,
    metadata = metadata,
    lane_step = 0.08,
    max_offset = 0.16,
    round_radius = 0.04
  )

  for (route_index in seq_along(routes)) {
    expected <- DNMB:::.dnmb_cct_round_orthogonal_route(
      unrounded[[route_index]], radius = 0.04
    )
    expect_equal(
      unname(as.matrix(rounded[[route_index]][, c("x", "y")])),
      unname(as.matrix(expected[, c("x", "y")]))
    )
  }

  separated <- rounded[[1L]]
  expect_equal(unlist(separated[1L, c("x", "y")]), c(x = 1, y = 5))
  expect_equal(
    unlist(separated[nrow(separated), c("x", "y")]),
    c(x = 5, y = 1)
  )
})
