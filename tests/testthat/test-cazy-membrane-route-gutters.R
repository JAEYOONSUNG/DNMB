.cazy_membrane_route_pair_key <- function(route) {
  paste(route$cs_id, route$target_key, sep = "\r")
}

.cazy_membrane_horizontal_levels <- function(points, tolerance = 1e-8) {
  if (is.null(points) || nrow(points) < 2L) return(numeric())
  dx <- diff(points$x)
  dy <- diff(points$y)
  unique(points$y[-nrow(points)][
    abs(dx) > tolerance & abs(dy) <= tolerance
  ])
}

.cazy_membrane_route_fixture <- function() {
  data.frame(
    cs_id = c("Maltose", "Galactose", "Mannitol", "Sucrose"),
    locus_tag = c("T_SHARED", "T_SHARED", "T_MTL", "T_SCR"),
    source_x = c(1.4, 3.4, 6.2, 8.2),
    source_y = c(10.5, 9.5, 9.5, 9.5),
    lane_x = c(1.4, 3.4, 6.2, 8.2),
    tx_draw = c(4.4, 4.4, 6.8, 7.6),
    ty_draw = c(8.62, 8.62, 8.38, 8.50),
    cyto_x = c(2.0, 3.0, 6.0, 8.0),
    cyto_y = c(7.75, 7.65, 7.70, 7.80),
    confidence = c("high", "high", "high", "high"),
    color = c("#8C510A", "#01665E", "#5E3C99", "#C51B7D"),
    stringsAsFactors = FALSE
  )
}

.cazy_membrane_drawn_rows <- function(layout, fixture) {
  drawn <- layout$memberships[layout$memberships$draw, , drop = FALSE]
  drawn_keys <- paste(drawn$cs_id, drawn$target_key, sep = "\r")
  fixture_keys <- paste(fixture$cs_id, fixture$locus_tag, sep = "\r")
  fixture[fixture_keys %in% drawn_keys, , drop = FALSE]
}

.cazy_membrane_routes_by_pair <- function(routes) {
  keys <- vapply(routes, .cazy_membrane_route_pair_key, character(1))
  routes[order(keys)]
}

test_that("transporter route lanes leave visible gutters above and below the membrane", {
  y_memb <- 8.5
  membrane_half_height <- 0.20
  glyph_clearance <- 0.11
  route_gap <- 0.08
  fixture <- .cazy_membrane_route_fixture()

  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(
    fixture,
    y_memb = y_memb,
    glyph_clearance = glyph_clearance,
    corner_radius = 0.08
  )
  selected <- .cazy_membrane_drawn_rows(layout, fixture)
  interior <- DNMB:::.dnmb_cct_transporter_interior_routes(
    selected,
    y_memb = y_memb,
    glyph_clearance = glyph_clearance,
    source_clearance = glyph_clearance,
    corner_radius = 0.08
  )

  expect_gt(length(layout$exterior_routes), 0L)
  expect_gt(length(interior), 0L)

  for (route in layout$exterior_routes) {
    target <- selected[
      selected$cs_id == route$cs_id &
        selected$locus_tag == route$target_key,
      , drop = FALSE
    ]
    expect_equal(nrow(target), 1L)
    upper_envelope <- max(
      y_memb + membrane_half_height,
      target$ty_draw + glyph_clearance
    )
    horizontal_y <- .cazy_membrane_horizontal_levels(route$points)
    expect_gt(length(horizontal_y), 0L)
    expect_true(all(horizontal_y >= upper_envelope + route_gap - 1e-10))
  }

  for (route in interior) {
    target <- selected[
      selected$cs_id == route$cs_id &
        selected$locus_tag == route$target_key,
      , drop = FALSE
    ]
    expect_equal(nrow(target), 1L)
    lower_envelope <- min(
      y_memb - membrane_half_height,
      target$ty_draw - glyph_clearance
    )
    horizontal_y <- .cazy_membrane_horizontal_levels(route$points)
    expect_gt(length(horizontal_y), 0L)
    expect_true(all(horizontal_y <= lower_envelope - route_gap + 1e-10))
  }
})

test_that("exterior and interior route halves meet the same selected T-ID", {
  y_memb <- 8.5
  glyph_clearance <- 0.11
  fixture <- .cazy_membrane_route_fixture()
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(
    fixture, y_memb = y_memb, glyph_clearance = glyph_clearance
  )
  selected <- .cazy_membrane_drawn_rows(layout, fixture)
  interior <- DNMB:::.dnmb_cct_transporter_interior_routes(
    selected, y_memb = y_memb,
    glyph_clearance = glyph_clearance,
    source_clearance = glyph_clearance
  )

  exterior_keys <- vapply(
    layout$exterior_routes, .cazy_membrane_route_pair_key, character(1)
  )
  interior_keys <- vapply(
    interior, .cazy_membrane_route_pair_key, character(1)
  )
  expect_setequal(exterior_keys, interior_keys)

  exterior_ids <- stats::setNames(vapply(
    layout$exterior_routes, `[[`, character(1), "route_id"
  ), exterior_keys)
  interior_ids <- stats::setNames(vapply(
    interior, `[[`, character(1), "route_id"
  ), interior_keys)
  common_keys <- sort(intersect(names(exterior_ids), names(interior_ids)))
  expect_identical(
    unname(exterior_ids[common_keys]),
    unname(interior_ids[common_keys])
  )

  for (pair_key in common_keys) {
    exterior <- layout$exterior_routes[[match(pair_key, exterior_keys)]]
    inside <- interior[[match(pair_key, interior_keys)]]
    target <- selected[
      paste(selected$cs_id, selected$locus_tag, sep = "\r") == pair_key,
      , drop = FALSE
    ]
    expect_equal(nrow(target), 1L)

    exterior_end <- exterior$points[nrow(exterior$points), c("x", "y")]
    interior_start <- inside$points[1L, c("x", "y")]
    expect_equal(
      unlist(exterior_end),
      c(x = target$tx_draw, y = target$ty_draw + glyph_clearance)
    )
    expect_equal(
      unlist(interior_start),
      c(x = target$tx_draw, y = target$ty_draw - glyph_clearance)
    )

    if (nrow(exterior$points) > 1L) {
      expect_equal(
        exterior$points$x[nrow(exterior$points) - 1L],
        target$tx_draw
      )
    }
    if (nrow(inside$points) > 1L) {
      expect_equal(inside$points$x[2L], target$tx_draw)
    }
  }
})

test_that("shared-transporter route pairing and geometry are deterministic", {
  fixture <- .cazy_membrane_route_fixture()
  shuffled_fixture <- fixture[c(4, 2, 1, 3), , drop = FALSE]

  build_routes <- function(rows) {
    layout <- DNMB:::.dnmb_cct_transporter_bus_layout(rows, y_memb = 8.5)
    selected <- .cazy_membrane_drawn_rows(layout, rows)
    interior <- DNMB:::.dnmb_cct_transporter_interior_routes(
      selected, y_memb = 8.5
    )
    list(
      layout = layout,
      exterior = .cazy_membrane_routes_by_pair(layout$exterior_routes),
      render = .cazy_membrane_routes_by_pair(layout$exterior_render_routes),
      interior = .cazy_membrane_routes_by_pair(interior)
    )
  }

  original <- build_routes(fixture)
  shuffled <- build_routes(shuffled_fixture)

  compare_routes <- function(first, second) {
    first_keys <- vapply(first, .cazy_membrane_route_pair_key, character(1))
    second_keys <- vapply(second, .cazy_membrane_route_pair_key, character(1))
    expect_identical(first_keys, second_keys)
    expect_identical(
      vapply(first, `[[`, character(1), "route_id"),
      vapply(second, `[[`, character(1), "route_id")
    )
    for (i in seq_along(first)) {
      expect_equal(
        unname(as.matrix(first[[i]]$points[, c("x", "y")])),
        unname(as.matrix(second[[i]]$points[, c("x", "y")]))
      )
    }
  }

  compare_routes(original$exterior, shuffled$exterior)
  compare_routes(original$render, shuffled$render)
  compare_routes(original$interior, shuffled$interior)

  shared_exterior <- Filter(
    function(route) identical(route$target_key, "T_SHARED"),
    original$exterior
  )
  shared_interior <- Filter(
    function(route) identical(route$target_key, "T_SHARED"),
    original$interior
  )
  expect_equal(length(shared_exterior), 2L)
  expect_equal(length(shared_interior), 2L)
  expect_true(all(vapply(
    shared_exterior,
    function(route) isTRUE(all.equal(route$target_x, 4.4)),
    logical(1)
  )))
  expect_true(all(vapply(
    shared_interior,
    function(route) isTRUE(all.equal(route$points$x[1L], 4.4)),
    logical(1)
  )))
})

test_that("every exterior membership keeps its colour to the final transporter", {
  fixture <- .cazy_membrane_route_fixture()
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(fixture, y_memb = 8.5)
  drawn <- layout$memberships[layout$memberships$draw, , drop = FALSE]

  expected_ids <- vapply(seq_len(nrow(drawn)), function(i) {
    DNMB:::.dnmb_cct_transporter_route_id(
      drawn$cs_id[i], drawn$target_key[i]
    )
  }, character(1))
  rendered_ids <- vapply(
    layout$exterior_render_routes, `[[`, character(1), "route_id"
  )

  # A shared transporter must not replace membership paths with one neutral
  # transport-trunk route. Each coloured source path remains independently
  # traceable through its final segment and arrowhead.
  expect_setequal(rendered_ids, expected_ids)
  expect_identical(anyDuplicated(rendered_ids), 0L)
  expect_false(any(startsWith(rendered_ids, "transport-trunk:")))

  fixture_keys <- paste(fixture$cs_id, fixture$locus_tag, sep = "\r")
  for (i in seq_len(nrow(drawn))) {
    route <- layout$exterior_render_routes[[match(expected_ids[i], rendered_ids)]]
    fixture_row <- fixture[
      fixture_keys == paste(drawn$cs_id[i], drawn$target_key[i], sep = "\r"),
      , drop = FALSE
    ]
    expect_equal(nrow(fixture_row), 1L)
    expect_identical(route$color, fixture_row$color)
    expect_true(isTRUE(route$arrow_last))
    expect_equal(route$transporter_x, fixture_row$tx_draw)
    expect_lte(abs(route$landing_x - fixture_row$tx_draw), 0.07 + 1e-10)
    expect_equal(
      unlist(route$points[nrow(route$points), c("x", "y")]),
      c(x = fixture_row$tx_draw, y = fixture_row$ty_draw + 0.11)
    )
    expect_equal(c(route$target_x, route$target_y), unname(unlist(
      route$points[nrow(route$points), c("x", "y")]
    )))
  }

  shared_render <- Filter(
    function(route) identical(route$target_key, "T_SHARED"),
    layout$exterior_render_routes
  )
  expect_equal(length(shared_render), 2L)
  expect_equal(
    length(unique(vapply(shared_render, `[[`, numeric(1), "landing_x"))),
    2L
  )
  expect_true(all(vapply(shared_render, function(route) {
    if (abs(route$landing_x - route$target_x) <= 1e-8) return(TRUE)
    penultimate <- route$points[nrow(route$points) - 1L, c("x", "y")]
    isTRUE(all.equal(
      unname(unlist(penultimate)),
      c(route$landing_x, route$landing_y)
    ))
  }, logical(1))))
})

test_that("an aligned shared membership retains its coloured route", {
  memberships <- data.frame(
    cs_id = c("Galactose", "Xylose"),
    locus_tag = c("T_SHARED", "T_SHARED"),
    source_x = c(4.4, 2.0), source_y = c(9.5, 9.5),
    lane_x = c(4.4, 2.0),
    tx_draw = c(4.4, 4.4), ty_draw = c(8.5, 8.5),
    confidence = c("high", "high"),
    color = c("#CCA800", "#01665E"),
    stringsAsFactors = FALSE
  )
  layout <- DNMB:::.dnmb_cct_transporter_bus_layout(memberships, y_memb = 8.5)
  rendered_ids <- vapply(
    layout$exterior_render_routes, `[[`, character(1), "route_id"
  )
  aligned_id <- DNMB:::.dnmb_cct_transporter_route_id(
    "Galactose", "T_SHARED"
  )
  aligned <- layout$exterior_render_routes[[match(aligned_id, rendered_ids)]]

  expect_false(is.null(aligned))
  expect_identical(aligned$color, "#CCA800")
  expect_true(isTRUE(aligned$arrow_last))
  expect_equal(
    unlist(aligned$points[nrow(aligned$points), c("x", "y")]),
    c(x = 4.4, y = 8.61)
  )
})

test_that("distinct loci at one substrate anchor get separate route centers", {
  ids <- c("L03", "L01", "L02")
  desired <- rep(10, 3)
  spans <- c(0.30, 0.26, 0.30)
  ranks <- rep(2, 3)

  pack <- function(idx) {
    packed <- DNMB:::.dnmb_cct_pack_transporters_lane(
      center_x = mean(desired[idx]),
      half_spans = spans[idx],
      lane_ranks = ranks[idx],
      desired_x = desired[idx],
      stable_ids = ids[idx],
      y_memb = 8.5,
      xlim = c(9.0, 11.0)
    )
    packed$locus_tag <- ids[idx]
    packed[order(packed$locus_tag), , drop = FALSE]
  }

  original <- pack(seq_along(ids))
  shuffled <- pack(c(3, 1, 2))

  expect_equal(length(unique(original$tx)), length(ids))
  pairs <- utils::combn(seq_len(nrow(original)), 2L)
  for (j in seq_len(ncol(pairs))) {
    ii <- pairs[, j]
    pair_spans <- spans[match(original$locus_tag[ii], ids)]
    expect_gte(abs(diff(original$tx[ii])), max(pair_spans) + 0.08 - 1e-10)
  }
  ordered_spans <- spans[match(original$locus_tag, ids)]
  expect_true(all(original$tx - ordered_spans >= 9.0 - 1e-10))
  expect_true(all(original$tx + ordered_spans <= 11.0 + 1e-10))
  for (row in unique(original$row_id)) {
    ii <- which(original$row_id == row)
    if (length(ii) < 2L) next
    ii <- ii[order(original$tx[ii])]
    edge_gap <- diff(original$tx[ii]) -
      head(ordered_spans[ii], -1L) - tail(ordered_spans[ii], -1L)
    expect_true(all(edge_gap >= 0.08 - 1e-10))
  }
  expect_equal(original$tx, shuffled$tx)
  expect_equal(original$ty, shuffled$ty)
  expect_identical(original$row_id, shuffled$row_id)
})
