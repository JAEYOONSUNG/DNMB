.cazy_symbol_core_half_span <- function(x, cx = 0, cy = 0) {
  layers <- if (inherits(x, "ggplot")) x$layers else x
  if (length(layers) < 2L) {
    stop("SNFG symbols must contain a halo and a core layer")
  }

  # Every SNFG renderer puts its white clearance halo first. The remaining
  # layers describe the actual symbol, its divider, and optional text.
  core_layers <- layers[-1L]
  spans <- vapply(core_layers, function(layer) {
    if (inherits(layer$geom, "GeomCustomAnn")) {
      bounds <- unlist(layer$geom_params[c("xmin", "xmax", "ymin", "ymax")])
      if (length(bounds) == 4L && all(is.finite(bounds))) {
        return(max(abs(bounds - c(cx, cx, cy, cy))))
      }
    }

    dat <- layer$data
    if (!is.data.frame(dat) || !nrow(dat)) return(0)
    x_cols <- intersect(c("x", "xend"), names(dat))
    y_cols <- intersect(c("y", "yend"), names(dat))
    offsets <- c(
      unlist(lapply(x_cols, function(nm) abs(as.numeric(dat[[nm]]) - cx))),
      unlist(lapply(y_cols, function(nm) abs(as.numeric(dat[[nm]]) - cy)))
    )
    offsets <- offsets[is.finite(offsets)]
    if (length(offsets)) max(offsets) else 0
  }, numeric(1))

  max(spans)
}

test_that("single SNFG shapes share one nominal icon half-span", {
  requested_radius <- 0.10
  representative_sugars <- c(
    circle = "galactose",
    square = "glcnac",
    diamond = "glca",
    triangle = "fucose",
    star = "xylose",
    pentagon = "fructose"
  )

  spans <- vapply(representative_sugars, function(sugar) {
    layers <- DNMB:::.dnmb_snfg_symbol_layers_v2(
      2, 3, sugar, size = requested_radius
    )
    .cazy_symbol_core_half_span(layers, cx = 2, cy = 3)
  }, numeric(1))

  expect_equal(
    unname(spans),
    rep(requested_radius, length(representative_sugars)),
    tolerance = 1e-10
  )
})

test_that("SNFG polygon shapes use a comparable centered bounding box", {
  radius <- DNMB:::.dnmb_cct_sugar_icon_radius()
  shape_types <- c("circle", "square", "diamond", "triangle", "star", "pentagon")
  bounds <- vapply(shape_types, function(shape_type) {
    points <- DNMB:::.dnmb_snfg_polygon_coords_v2(
      shape_type, x = 0, y = 0, size = radius
    )
    c(
      width = diff(range(points$x)),
      height = diff(range(points$y)),
      center_x = mean(range(points$x)),
      center_y = mean(range(points$y))
    )
  }, numeric(4))

  expect_equal(
    unname(apply(bounds[c("width", "height"), , drop = FALSE], 2L, max)),
    rep(2 * radius, length(shape_types)),
    tolerance = 1e-10
  )
  expect_true(all(bounds[c("width", "height"), ] >= 1.70 * radius))
  expect_equal(
    max(abs(bounds[c("center_x", "center_y"), ])),
    0,
    tolerance = 1e-10
  )
})

test_that("single-sugar carbon sources preserve the requested icon radius", {
  requested_radius <- 0.10
  source_ids <- c("Galactose", "NAG", "Gluconate", "Fucose", "Xylose", "Fructose")

  spans <- vapply(source_ids, function(source_id) {
    plot <- DNMB:::.dnmb_cct_render_carbon_source_node(
      ggplot2::ggplot(), 2, 3, source_id, r = requested_radius
    )
    .cazy_symbol_core_half_span(plot, cx = 2, cy = 3)
  }, numeric(1))

  expect_equal(
    unname(spans),
    rep(requested_radius, length(source_ids)),
    tolerance = 1e-10
  )
})

test_that("direct and carbon-source monosaccharides use the same map footprint", {
  requested_radius <- 0.10
  direct <- DNMB:::.dnmb_snfg_symbol_layers_v2(
    2, 3, "galactose", size = requested_radius
  )
  source <- DNMB:::.dnmb_cct_render_carbon_source_node(
    ggplot2::ggplot(), 2, 3, "Galactose", r = requested_radius
  )
  default_source <- DNMB:::.dnmb_cct_render_carbon_source_node(
    ggplot2::ggplot(), 2, 3, "Galactose"
  )

  direct_span <- .cazy_symbol_core_half_span(direct, cx = 2, cy = 3)
  source_span <- .cazy_symbol_core_half_span(source, cx = 2, cy = 3)
  default_span <- .cazy_symbol_core_half_span(default_source, cx = 2, cy = 3)

  expect_equal(source_span, direct_span, tolerance = 1e-10)
  expect_equal(default_span, direct_span, tolerance = 1e-10)
})

test_that("composite sugars keep full-size non-overlapping residue icons", {
  requested_radius <- DNMB:::.dnmb_cct_sugar_icon_radius()
  plot <- DNMB:::.dnmb_cct_render_carbon_source_node(
    ggplot2::ggplot(), 2, 3, "Maltose", r = requested_radius
  )
  custom_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomCustomAnn"),
    plot$layers
  )
  bounds <- lapply(custom_layers, function(layer) layer$geom_params)
  radii <- vapply(bounds, function(item) {
    (as.numeric(item$xmax) - as.numeric(item$xmin)) / 2
  }, numeric(1))
  centers <- vapply(bounds, function(item) {
    (as.numeric(item$xmax) + as.numeric(item$xmin)) / 2
  }, numeric(1))
  core_centers <- sort(centers[abs(radii - requested_radius) <= 1e-10])

  expect_equal(requested_radius, 0.10)
  expect_equal(length(core_centers), 2L)
  expect_equal(
    unname(diff(core_centers)),
    2.18 * requested_radius,
    tolerance = 1e-10
  )
  expect_gt(diff(core_centers), 2 * requested_radius)
})
