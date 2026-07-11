test_that("eggNOG database readiness requires all v5 assets and records verified release", {
  data_dir <- tempfile("eggnog-data-")
  cache_root <- tempfile("eggnog-cache-")
  dir.create(data_dir)
  dir.create(cache_root)
  on.exit(unlink(c(data_dir, cache_root), recursive = TRUE, force = TRUE), add = TRUE)
  required <- c("eggnog_proteins.dmnd", "eggnog.db", "eggnog.taxa.db")
  invisible(lapply(file.path(data_dir, required), writeLines, text = "fixture"))

  local_mocked_bindings(
    dnmb_detect_emapper = function(required = FALSE) {
      list(found = TRUE, path = "/test/emapper.py", version = "2.1.15")
    },
    .dnmb_eggnog_database_release_from_db = function(data_dir) "5.0.2",
    .package = "DNMB"
  )

  complete <- DNMB:::dnmb_eggnog_ensure_database(
    data_dir = data_dir,
    install = FALSE,
    cache_root = cache_root
  )
  expect_true(complete$ok)
  expect_identical(complete$manifest$db_release, "5.0.2")
  expect_identical(complete$manifest$expected_db_release, "5.0.2")

  unlink(file.path(data_dir, "eggnog.taxa.db"))
  incomplete <- DNMB:::dnmb_eggnog_ensure_database(
    data_dir = data_dir,
    install = FALSE,
    cache_root = cache_root
  )
  expect_false(incomplete$ok)
})

test_that("eggNOG database release is read from its SQLite version table", {
  data_dir <- tempfile("eggnog-version-")
  dir.create(data_dir)
  on.exit(unlink(data_dir, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("fixture", file.path(data_dir, "eggnog.db"))

  local_mocked_bindings(
    dnmb_detect_binary = function(binary, required = FALSE) {
      list(found = identical(binary, "sqlite3"), path = binary)
    },
    dnmb_run_external = function(command, args = character(), ...) {
      list(ok = TRUE, stdout = "5.0.2", stderr = character())
    },
    .package = "DNMB"
  )

  expect_identical(DNMB:::.dnmb_eggnog_database_release_from_db(data_dir), "5.0.2")
})
