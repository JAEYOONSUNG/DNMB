test_that("release identities normalize tags without treating cache aliases as versions", {
  expect_true(DNMB:::.dnmb_db_release_ids_equal("2.0.1", "v2.0.1"))
  expect_true(DNMB:::.dnmb_db_release_ids_equal("5.73-104.0", "5.73-104.0"))
  expect_true(DNMB:::.dnmb_db_release_ids_equal(
    "abcdef123456",
    "abcdef1234567890abcdef1234567890abcdef12",
    kind = "commit"
  ))
  expect_true(DNMB:::.dnmb_db_is_symbolic_version("current"))
  expect_true(DNMB:::.dnmb_db_is_symbolic_version("LATEST"))
  expect_false(DNMB:::.dnmb_db_is_symbolic_version("v2.0.1"))
})

test_that("manifest identity uses real version fields before the cache key", {
  manifest <- list(
    version = "current",
    tool_version = "2.0.1",
    source_commit = "abcdef1234567890"
  )
  identity <- DNMB:::.dnmb_db_manifest_scalar(
    manifest,
    c("tool_version", "source_commit", "version"),
    skip_symbolic = TRUE
  )
  expect_identical(identity$field, "tool_version")
  expect_identical(identity$value, "2.0.1")

  unknown <- DNMB:::.dnmb_db_manifest_scalar(
    list(version = "current"),
    "version",
    skip_symbolic = TRUE
  )
  expect_true(is.na(unknown$value))
})

test_that("DefenseFinder compares tool and model releases, never written_at", {
  local_mocked_bindings(
    .dnmb_db_github_latest_release = function(repo) {
      if (identical(repo, "mdmparis/defense-finder")) {
        list(tag = "v3.0.0", name = "v3.0.0")
      } else {
        list(tag = "v3.1.0", name = "v3.1.0")
      }
    },
    .package = "DNMB"
  )

  current <- DNMB:::.dnmb_db_remote_check_defensefinder(list(
    version = "current",
    written_at = "1999-01-01T00:00:00Z",
    tool_version = "3.0.0",
    models_version = "3.1.0"
  ))
  expect_false(current$update_available)

  stale <- DNMB:::.dnmb_db_remote_check_defensefinder(list(
    version = "current",
    written_at = "2026-05-05T03:00:00Z",
    tool_version = "2.0.1",
    models_version = "3.1.0"
  ))
  expect_true(stale$update_available)

  legacy <- DNMB:::.dnmb_db_remote_check_defensefinder(list(
    version = "current",
    written_at = "3.0.0",
    models_version = "3.1.0"
  ))
  expect_true(is.na(legacy$update_available))
})

test_that("dbCAN uses tool and database releases rather than current", {
  local_mocked_bindings(
    .dnmb_db_github_latest_release = function(repo) {
      expect_identical(repo, "bcb-unl/run_dbcan")
      list(tag = "v5.2.9", name = "v5.2.9")
    },
    .dnmb_dbcan_release_info = function() list(
      version = "V14",
      source = "https://pro.unl.edu/dbCAN2/"
    ),
    .package = "DNMB"
  )

  current <- DNMB:::.dnmb_db_remote_check_dbcan(list(
    version = "current",
    tool_version = "5.2.9",
    resolved_release_version = "v14"
  ))
  expect_false(current$update_available)

  stale_tool <- DNMB:::.dnmb_db_remote_check_dbcan(list(
    version = "current",
    tool_version = "5.2.8",
    resolved_release_version = "V14"
  ))
  expect_true(stale_tool$update_available)

  unknown <- DNMB:::.dnmb_db_remote_check_dbcan(list(version = "current"))
  expect_true(is.na(unknown$update_available))
})

test_that("CLEAN compares repository commits and ignores split100", {
  remote_sha <- "abcdef1234567890abcdef1234567890abcdef12"
  local_mocked_bindings(
    .dnmb_db_github_commit = function(repo, ref = "HEAD") remote_sha,
    .dnmb_db_github_latest_release = function(repo) list(tag = "v1.0.0", name = "v1.0.0"),
    .package = "DNMB"
  )

  current <- DNMB:::.dnmb_db_remote_check_clean(list(
    version = "split100",
    repo_source = "https://github.com/tttianhao/CLEAN.git",
    source_commit = substr(remote_sha, 1L, 12L)
  ))
  expect_false(current$update_available)

  stale <- DNMB:::.dnmb_db_remote_check_clean(list(
    version = "split100",
    repo_source = "https://github.com/tttianhao/CLEAN.git",
    source_commit = "1234567890abcdef"
  ))
  expect_true(stale$update_available)

  unknown <- DNMB:::.dnmb_db_remote_check_clean(list(
    version = "split100",
    repo_source = "https://github.com/tttianhao/CLEAN.git"
  ))
  expect_true(is.na(unknown$update_available))
})

test_that("InterProScan uses the resolved tool release", {
  local_mocked_bindings(
    .dnmb_db_github_latest_release = function(repo) {
      expect_identical(repo, "ebi-pf-team/interproscan")
      list(tag = "5.73-104.0", name = "5.73-104.0")
    },
    .package = "DNMB"
  )

  current <- DNMB:::.dnmb_db_remote_check_interproscan(list(
    version = "current",
    tool_version = "5.73-104.0"
  ))
  expect_false(current$update_available)

  unknown <- DNMB:::.dnmb_db_remote_check_interproscan(list(version = "current"))
  expect_true(is.na(unknown$update_available))
})

test_that("unknown remote comparisons are not presented as latest", {
  cache_root <- tempfile("dnmb-freshness-")
  dir.create(cache_root)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)
  DNMB:::dnmb_db_write_manifest(
    "clean",
    "split100",
    manifest = list(install_ok = TRUE),
    cache_root = cache_root,
    overwrite = TRUE
  )
  local_mocked_bindings(
    .dnmb_db_check_remote_version = function(module, manifest) {
      list(remote_version = "commit:abcdef123456", update_available = NA)
    },
    .package = "DNMB"
  )

  expect_message(
    expect_true(DNMB:::dnmb_db_check_freshness(
      "clean", "split100", cache_root = cache_root, verbose = TRUE
    )),
    "local release identity unavailable"
  )
})

test_that("InterProScan cache manifests record the resolved release and hashes", {
  cache_root <- tempfile("dnmb-interproscan-manifest-")
  dir.create(cache_root)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)
  version <- "5.73-104.0"
  ipr_dir <- file.path(cache_root, "db_modules", "interproscan", version)
  dir.create(ipr_dir, recursive = TRUE)
  writeLines("#!/bin/sh", file.path(ipr_dir, "interproscan.sh"))

  DNMB:::.dnmb_interproscan_write_manifest(
    version = version,
    ipr_dir = ipr_dir,
    cache_root = cache_root,
    source_url = "https://example.test/interproscan.tar.gz",
    archive_md5 = "0123456789abcdef0123456789abcdef"
  )
  manifest <- DNMB:::dnmb_db_read_manifest(
    "interproscan", version, cache_root = cache_root, required = TRUE
  )

  expect_identical(manifest$tool_version, version)
  expect_identical(manifest$resolved_release_version, version)
  expect_match(manifest$launcher_md5, "^[0-9a-f]{32}$")
  expect_identical(manifest$archive_md5, "0123456789abcdef0123456789abcdef")
})

test_that("eggNOG freshness follows Bioconda and compares version order", {
  local_mocked_bindings(
    .dnmb_db_read_json_url = function(url) {
      if (grepl("api[.]anaconda[.]org", url)) {
        list(latest_version = "2.1.15")
      } else {
        list(info = list(version = "2.1.13"))
      }
    },
    .package = "DNMB"
  )

  current <- DNMB:::.dnmb_db_remote_check_eggnog(list(tool_version = "2.1.15"))
  expect_false(current$update_available)

  stale <- DNMB:::.dnmb_db_remote_check_eggnog(list(tool_version = "2.1.14"))
  expect_true(stale$update_available)

  newer_local <- DNMB:::.dnmb_db_remote_check_eggnog(list(tool_version = "2.1.16"))
  expect_false(newer_local$update_available)

  unknown <- DNMB:::.dnmb_db_remote_check_eggnog(list(version = "data"))
  expect_true(is.na(unknown$update_available))
})
