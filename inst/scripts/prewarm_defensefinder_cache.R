#!/usr/bin/env Rscript

cache_root <- Sys.getenv("DNMB_CACHE_ROOT", unset = "")
if (!nzchar(cache_root)) {
  cache_root <- NULL
}

repo_dir <- Sys.getenv("DNMB_DEFENSEFINDER_REPO_DIR", unset = "")
casfinder_dir <- Sys.getenv("DNMB_DEFENSEFINDER_CASFINDER_DIR", unset = "")

asset_urls <- list()
if (nzchar(casfinder_dir)) {
  asset_urls$casfinder_dir <- casfinder_dir
}

repo_source <- if (nzchar(repo_dir)) repo_dir else DNMB:::.dnmb_defensefinder_default_repo_url()

result <- DNMB:::dnmb_defensefinder_install_module(
  version = DNMB:::.dnmb_defensefinder_default_version(),
  cache_root = cache_root,
  install = TRUE,
  repo_url = repo_source,
  asset_urls = asset_urls,
  force = TRUE
)

print(result$status)

if (!isTRUE(result$ok)) {
  quit(status = 1L)
}

invisible(result)
