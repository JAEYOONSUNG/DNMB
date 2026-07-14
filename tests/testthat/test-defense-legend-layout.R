test_that("integrated defense legend layout preserves natural legend height", {
  layout <- DNMB:::.dnmb_integrated_defense_legend_layout(
    meta_height_in = 0.32,
    fill_height_in = 0.74,
    main_rel_heights = c(2.05, 4.00, 1.50),
    main_height_in = 14.85,
    inter_legend_gap_in = 0.08,
    outer_padding_in = 0.06
  )

  expect_equal(layout$legend_height_in, 1.26, tolerance = 1e-12)
  expect_equal(layout$output_height_in, 16.11, tolerance = 1e-12)
  allocated_height <- layout$output_height_in * layout$legend_rel_height /
    sum(c(layout$main_rel_heights, layout$legend_rel_height))
  expect_equal(allocated_height, layout$legend_height_in, tolerance = 1e-12)

  taller <- DNMB:::.dnmb_integrated_defense_legend_layout(
    meta_height_in = 0.32,
    fill_height_in = 1.10
  )
  expect_gt(taller$legend_rel_height, layout$legend_rel_height)
  expect_gt(taller$output_height_in, layout$output_height_in)
})

test_that("integrated defense legend labels retain full subtype names", {
  palette <- c(
    "CAS_Class1-Subtyp..." = "#112233",
    "AbiE" = "#445566"
  )
  members <- data.frame(
    display_name = c("CAS_Class1-Subtyp...", "AbiE"),
    display_name_full = c("CAS_Class1-Subtype-IV-D", "AbiE"),
    stringsAsFactors = FALSE
  )

  labels <- DNMB:::.dnmb_integrated_defense_legend_labels(
    palette,
    members,
    wrap_width = 18L
  )

  expect_identical(names(labels), names(palette))
  expect_identical(gsub("\n", "", labels[[1]], fixed = TRUE), "CAS_Class1-Subtype-IV-D")
  expect_false(any(grepl("...", labels, fixed = TRUE)))
  expect_identical(labels[[2]], "AbiE")
})

test_that("DefenseFinder occupancy collector samples both ends and midpoint", {
  points <- DNMB:::.dnmb_defensefinder_legend_occupancy_points(
    data.frame(x = 1, y = 2, xend = 5, yend = 6),
    weight = 3,
    radius = 0.02
  )

  expect_equal(nrow(points), 3L)
  expect_equal(points$x, c(1, 5, 3))
  expect_equal(points$y, c(2, 6, 4))
  expect_true(all(points$weight == 3))
  expect_true(all(points$radius == 0.02))
})

test_that("DefenseFinder legend chooses the only unoccupied corner", {
  empty <- data.frame(x = numeric(), y = numeric())
  baseline <- DNMB:::.dnmb_defensefinder_auto_legend_corner(
    occupancy = empty,
    xlim = c(0, 1),
    ylim = c(0, 1),
    legend_labels = paste0("Type ", 1:4)
  )
  expect_identical(baseline$corner, "bottom-right")

  blocked_boxes <- baseline$candidates[
    baseline$candidates$corner != "top-left",
    ,
    drop = FALSE
  ]
  occupancy <- data.frame(
    x = (blocked_boxes$xmin + blocked_boxes$xmax) / 2,
    y = (blocked_boxes$ymin + blocked_boxes$ymax) / 2,
    weight = 1,
    radius = 0.01
  )
  selected <- DNMB:::.dnmb_defensefinder_auto_legend_corner(
    occupancy = occupancy,
    xlim = c(0, 1),
    ylim = c(0, 1),
    legend_labels = paste0("Type ", 1:4)
  )

  expect_identical(selected$corner, "top-left")
  expect_equal(
    selected$candidates$overlap_score[selected$candidates$corner == "top-left"],
    0
  )
})

test_that("DefenseFinder legend avoids a scale-bar-like occupied corner", {
  baseline <- DNMB:::.dnmb_defensefinder_auto_legend_corner(
    occupancy = data.frame(x = numeric(), y = numeric()),
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    legend_labels = c("RM Type I", "Wadjet II")
  )
  bottom_right <- baseline$candidates[
    baseline$candidates$corner == "bottom-right",
    ,
    drop = FALSE
  ]
  scale_bar <- data.frame(
    x = mean(c(bottom_right$xmin, bottom_right$xmax)) * 20 - 10,
    y = mean(c(bottom_right$ymin, bottom_right$ymax)) * 20 - 10,
    weight = 80,
    radius = 0.03
  )
  selected <- DNMB:::.dnmb_defensefinder_auto_legend_corner(
    occupancy = scale_bar,
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    legend_labels = c("RM Type I", "Wadjet II")
  )

  expect_false(identical(selected$corner, "bottom-right"))
  expect_gt(
    selected$candidates$score[selected$candidates$corner == "bottom-right"],
    selected$candidates$score[selected$candidates$corner == selected$corner]
  )
})

test_that("DefenseFinder legend corner validation rejects degenerate limits", {
  expect_error(
    DNMB:::.dnmb_defensefinder_auto_legend_corner(
      occupancy = data.frame(x = 0, y = 0),
      xlim = c(1, 1),
      ylim = c(0, 1),
      legend_labels = "Type I"
    ),
    "non-zero range"
  )
})

test_that("DefenseFinder legend corner padding supports zoomed panel insets", {
  selected <- DNMB:::.dnmb_defensefinder_auto_legend_corner(
    occupancy = data.frame(x = numeric(), y = numeric()),
    xlim = c(0, 1),
    ylim = c(0, 1),
    legend_labels = c("Type I", "Type II"),
    corner_padding = 0.10
  )

  bottom_left <- selected$candidates[
    selected$candidates$corner == "bottom-left",
    ,
    drop = FALSE
  ]
  expect_equal(bottom_left$position_x, 0.10)
  expect_equal(bottom_left$position_y, 0.10)
})
