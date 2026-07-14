.dnmb_plot_font_family <- function() {
  "Arial"
}

.dnmb_plot_pdf_device <- function(filename,
                                  width = 7,
                                  height = 7,
                                  pointsize = 12,
                                  onefile = TRUE,
                                  bg = "white",
                                  antialias = "default",
                                  fallback_resolution = 300,
                                  ...) {
  if (!isTRUE(base::capabilities("cairo"))) {
    base::stop(
      "DNMB PDF export requires Cairo so Arial text remains embedded and editable.",
      call. = FALSE
    )
  }

  dots <- list(...)
  reserved <- base::intersect(base::names(dots), c("family", "symbolfamily"))
  if (base::length(reserved)) {
    base::stop(
      "`.dnmb_plot_pdf_device()` controls `family` and `symbolfamily`; do not override them.",
      call. = FALSE
    )
  }

  base::do.call(
    grDevices::cairo_pdf,
    c(
      list(
        filename = filename,
        width = width,
        height = height,
        pointsize = pointsize,
        onefile = onefile,
        family = .dnmb_plot_font_family(),
        bg = bg,
        antialias = antialias,
        fallback_resolution = fallback_resolution,
        symbolfamily = .dnmb_plot_font_family()
      ),
      dots
    )
  )
}

.dnmb_with_plot_pdf_device <- function(expr, width = 7, height = 7) {
  path <- tempfile("dnmb-plot-layout-", fileext = ".pdf")
  cowplot_paths <- character()
  cowplot_device <- function(width = 6, height = 6) {
    cowplot_path <- tempfile("dnmb-cowplot-layout-", fileext = ".pdf")
    cowplot_paths <<- c(cowplot_paths, cowplot_path)
    .dnmb_plot_pdf_device(cowplot_path, width = width, height = height)
  }
  old_cowplot_device <- cowplot::set_null_device(cowplot_device)
  .dnmb_plot_pdf_device(path, width = width, height = height)
  device_id <- grDevices::dev.cur()
  on.exit({
    cowplot::set_null_device(old_cowplot_device)
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(which = device_id)
    }
    unlink(c(path, cowplot_paths), force = TRUE)
  }, add = TRUE)
  force(expr)
}
