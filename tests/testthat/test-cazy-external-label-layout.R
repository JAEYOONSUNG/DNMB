test_that("fixed extracellular labels retain regular chain-relative coordinates", {
  labels <- data.frame(
    x = c(1.0, 1.2, 1.4),
    y = c(3.27, 3.13, 2.90),
    label = c("Cellulose", "\u03b21-4", "\u03b11-3"),
    label_kind = c("title", "bond", "branch_bond"),
    anchor_x = c(1.0, 1.2, 1.4),
    anchor_y = c(3.0, 3.0, 3.0),
    priority = c(30, 8, 10),
    hjust = c(0.5, 0.5, 0),
    color = c("#333333", "#777777", "#795548"),
    stringsAsFactors = FALSE
  )

  placed <- DNMB:::.dnmb_cct_layout_external_labels(labels, xlim = c(0, 4))

  expect_equal(placed$x_lab, labels$x)
  expect_equal(placed$y_lab, labels$y)
  expect_identical(placed$row_level, rep(0L, nrow(labels)))
  expect_false(any(placed$draw_leader))
})

test_that("G-ID labels use one primary and at most two auxiliary rows", {
  labels <- data.frame(
    x = rep(2, 7),
    y = rep(3, 7),
    label = sprintf("G%02d", 1:7),
    label_kind = "gh",
    anchor_x = rep(2, 7),
    anchor_y = rep(3.18, 7),
    priority = 7:1,
    hjust = 0.5,
    color = "#C62828",
    stringsAsFactors = FALSE
  )

  placed <- DNMB:::.dnmb_cct_layout_external_labels(
    labels, xlim = c(0, 4), row_step = 0.15, max_aux_rows = 2L
  )

  expect_setequal(unique(placed$row_level), 0:2)
  expect_lte(max(placed$row_level), 2L)
  expect_true(any(placed$row_level > 0L))
  expect_true(all(placed$y_lab == 3 - placed$row_level * 0.15))
  expect_identical(placed$row_level[1:3], rep(0L, 3))
  expect_identical(placed$row_level[4:6], rep(1L, 3))
  expect_identical(placed$row_level[7], 2L)
  expect_equal(placed$x_lab[1:3], c(2.00, 2.28, 1.72), tolerance = 1e-8)

  for (yy in unique(placed$y_lab)) {
    row <- placed[placed$y_lab == yy, , drop = FALSE]
    if (nrow(row) < 2L) next
    row <- row[order(row$x_lab), , drop = FALSE]
    widths <- pmax(0.20, nchar(row$label) * 0.055)
    required <- (head(widths, -1L) + tail(widths, -1L)) / 2 + 0.08
    expect_true(all(diff(row$x_lab) >= required - 1e-8))
  }
})

test_that("spaced G-ID labels stay on the primary row without leaders", {
  labels <- data.frame(
    x = c(0.5, 1.5, 2.5),
    y = rep(3, 3),
    label = c("G01", "G02", "G03"),
    label_kind = "gh",
    anchor_x = c(0.5, 1.5, 2.5),
    anchor_y = rep(3.18, 3),
    priority = 1,
    hjust = 0.5,
    color = "#C62828",
    stringsAsFactors = FALSE
  )

  placed <- DNMB:::.dnmb_cct_layout_external_labels(labels, xlim = c(0, 3))

  expect_identical(placed$row_level, rep(0L, 3))
  expect_equal(placed$x_lab, labels$anchor_x)
  expect_false(any(placed$draw_leader))
  expect_equal(nrow(DNMB:::.dnmb_cct_external_label_leaders(placed)), 0L)
})

test_that("external G-ID layout is deterministic across input order", {
  labels <- data.frame(
    x = c(1, 1, 1.12, 1.15, 2),
    y = c(3, 3, 3, 3, 2.7),
    label = c("G05", "G01", "G03", "G02", "G04"),
    label_kind = "gh",
    anchor_x = c(1, 1, 1.12, 1.15, 2),
    anchor_y = c(3.18, 3.18, 3.18, 3.18, 2.88),
    priority = c(1, 5, 3, 4, 2),
    hjust = 0.5,
    color = "#C62828",
    stringsAsFactors = FALSE
  )

  a <- DNMB:::.dnmb_cct_layout_external_labels(labels, xlim = c(0, 3))
  b <- DNMB:::.dnmb_cct_layout_external_labels(labels[5:1, ], xlim = c(0, 3))
  a <- a[order(a$label), c("label", "x_lab", "y_lab", "row_level")]
  b <- b[order(b$label), c("label", "x_lab", "y_lab", "row_level")]
  rownames(a) <- NULL
  rownames(b) <- NULL

  expect_identical(a, b)
})
