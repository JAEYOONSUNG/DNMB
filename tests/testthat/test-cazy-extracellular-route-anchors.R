.cazy_extra_route_spec <- function(route_id, source_kind, source_id,
                                   target_kind, target_id, points) {
  list(
    route_id = route_id,
    source_kind = source_kind,
    source_id = source_id,
    target_kind = target_kind,
    target_id = target_id,
    points = points,
    group = "release",
    priority = 1,
    color = "#1565C0",
    linewidth = 0.3,
    alpha = 0.5,
    linetype = "solid",
    gradient = FALSE,
    trim_start = 0.06,
    trim_end = 0.06
  )
}

.cazy_extra_route_fixture <- function() {
  anchors <- data.frame(
    kind = c("scissor", "branch_sugar", "monosaccharide", "transporter"),
    id = c("cut:starch", "branch:pectin:1", "glucose", "T01"),
    x = c(1, 3, 7, 9),
    y = c(12, 11, 9.5, 8.5),
    stringsAsFactors = FALSE
  )
  specs <- list(
    valid = .cazy_extra_route_spec(
      "valid", "scissor", "cut:starch", "monosaccharide", "glucose",
      data.frame(x = c(1, 1, 7, 7), y = c(12, 11.4, 11.4, 9.5))
    ),
    near = .cazy_extra_route_spec(
      "near", "branch_sugar", "branch:pectin:1", "monosaccharide", "glucose",
      data.frame(
        x = c(3 + 5e-7, 3, 7, 7 - 5e-7),
        y = c(11 - 5e-7, 10.8, 10.8, 9.5 + 5e-7)
      )
    ),
    missing_source = .cazy_extra_route_spec(
      "missing-source", "scissor", "cut:absent", "monosaccharide", "glucose",
      data.frame(x = c(4, 4, 7, 7), y = c(12, 11.2, 11.2, 9.5))
    ),
    missing_target = .cazy_extra_route_spec(
      "missing-target", "scissor", "cut:starch", "transporter", "T99",
      data.frame(x = c(1, 1, 8, 8), y = c(12, 10.5, 10.5, 8.5))
    ),
    detached_source = .cazy_extra_route_spec(
      "detached-source", "scissor", "cut:starch", "monosaccharide", "glucose",
      data.frame(x = c(2, 2, 7, 7), y = c(12, 10.8, 10.8, 9.5))
    ),
    detached_target = .cazy_extra_route_spec(
      "detached-target", "scissor", "cut:starch", "monosaccharide", "glucose",
      data.frame(x = c(1, 1, 6, 6), y = c(12, 10.6, 10.6, 9.5))
    ),
    nonfinite = .cazy_extra_route_spec(
      "nonfinite", "scissor", "cut:starch", "monosaccharide", "glucose",
      data.frame(x = c(1, 1, NA_real_, 7), y = c(12, 10.4, 10.4, 9.5))
    )
  )
  list(anchors = anchors, specs = specs)
}

.cazy_extra_route_directions <- function(points, tolerance = 1e-8) {
  if (is.null(points) || nrow(points) < 2L) return(character())
  dx <- diff(points$x)
  dy <- diff(points$y)
  direction <- ifelse(
    abs(dx) <= tolerance & abs(dy) <= tolerance, "Z",
    ifelse(abs(dy) <= tolerance, "H", ifelse(abs(dx) <= tolerance, "V", "D"))
  )
  direction[direction != "Z"]
}

.cazy_extra_has_sharp_orthogonal_corner <- function(points) {
  direction <- .cazy_extra_route_directions(points)
  if (length(direction) < 2L) return(FALSE)
  adjacent <- paste(direction[-length(direction)], direction[-1L], sep = "")
  any(adjacent %in% c("HV", "VH"))
}

test_that("extracellular routes require visible source and target anchors", {
  fixture <- .cazy_extra_route_fixture()
  validated <- DNMB:::.dnmb_cct_validate_extracellular_routes(
    fixture$specs, fixture$anchors, tolerance = 1e-6
  )

  route_ids <- vapply(validated, `[[`, character(1), "route_id")
  expect_setequal(route_ids, c("valid", "near"))
  expect_false(any(c(
    "missing-source", "missing-target", "detached-source",
    "detached-target", "nonfinite"
  ) %in% route_ids))
})

test_that("validated extracellular endpoints snap to their drawn anchors", {
  fixture <- .cazy_extra_route_fixture()
  validated <- DNMB:::.dnmb_cct_validate_extracellular_routes(
    fixture$specs, fixture$anchors, tolerance = 1e-6
  )

  for (spec in validated) {
    source <- fixture$anchors[
      fixture$anchors$kind == spec$source_kind &
        fixture$anchors$id == spec$source_id,
      , drop = FALSE
    ]
    target <- fixture$anchors[
      fixture$anchors$kind == spec$target_kind &
        fixture$anchors$id == spec$target_id,
      , drop = FALSE
    ]
    expect_equal(nrow(source), 1L)
    expect_equal(nrow(target), 1L)
    expect_equal(
      unlist(spec$points[1L, c("x", "y")]),
      c(x = source$x, y = source$y),
      tolerance = 0
    )
    expect_equal(
      unlist(spec$points[nrow(spec$points), c("x", "y")]),
      c(x = target$x, y = target$y),
      tolerance = 0
    )
  }
})

test_that("extracellular lanes are offset before every corner is rounded", {
  anchors <- data.frame(
    kind = c("scissor", "scissor", "monosaccharide", "monosaccharide"),
    id = c("cut:a", "cut:b", "glucose", "galactose"),
    x = c(1, 2, 7, 8),
    y = c(12, 12, 9.5, 9.5),
    stringsAsFactors = FALSE
  )
  specs <- list(
    .cazy_extra_route_spec(
      "release-a", "scissor", "cut:a", "monosaccharide", "glucose",
      data.frame(x = c(1, 1, 7, 7), y = c(12, 11, 11, 9.5))
    ),
    .cazy_extra_route_spec(
      "release-b", "scissor", "cut:b", "monosaccharide", "galactose",
      data.frame(x = c(2, 2, 8, 8), y = c(12, 11, 11, 9.5))
    )
  )
  validated <- DNMB:::.dnmb_cct_validate_extracellular_routes(
    specs, anchors, tolerance = 1e-6
  )
  metadata <- data.frame(
    route_id = vapply(validated, `[[`, character(1), "route_id"),
    priority = vapply(validated, `[[`, numeric(1), "priority"),
    group = vapply(validated, `[[`, character(1), "group"),
    stringsAsFactors = FALSE
  )
  raw_lanes <- DNMB:::.dnmb_cct_separate_route_lanes(
    lapply(validated, `[[`, "points"), metadata,
    lane_step = 0.08, max_offset = 0.16,
    connection_mode = "move_internal", round_radius = 0
  )
  rounded_lanes <- DNMB:::.dnmb_cct_separate_route_lanes(
    lapply(validated, `[[`, "points"), metadata,
    lane_step = 0.08, max_offset = 0.16,
    connection_mode = "move_internal", round_radius = 0.08
  )

  for (route_index in seq_along(raw_lanes)) {
    expected <- DNMB:::.dnmb_cct_round_orthogonal_route(
      raw_lanes[[route_index]], radius = 0.08
    )
    expect_equal(
      unname(as.matrix(rounded_lanes[[route_index]][, c("x", "y")])),
      unname(as.matrix(expected[, c("x", "y")]))
    )
    expect_false(.cazy_extra_has_sharp_orthogonal_corner(rounded_lanes[[route_index]]))
    expect_true("D" %in% .cazy_extra_route_directions(rounded_lanes[[route_index]]))
    expect_equal(
      unlist(rounded_lanes[[route_index]][1L, c("x", "y")]),
      unlist(specs[[route_index]]$points[1L, c("x", "y")])
    )
    expect_equal(
      unlist(rounded_lanes[[route_index]][nrow(rounded_lanes[[route_index]]), c("x", "y")]),
      unlist(specs[[route_index]]$points[nrow(specs[[route_index]]$points), c("x", "y")])
    )
  }
})
