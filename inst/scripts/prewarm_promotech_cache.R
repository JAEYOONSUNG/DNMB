#!/usr/bin/env Rscript

cache_root <- Sys.getenv("DNMB_CACHE_ROOT", unset = "")
if (!nzchar(cache_root)) {
  cache_root <- NULL
}

repo_dir <- Sys.getenv("DNMB_PROMOTECH_REPO_DIR", unset = "")
model_dir <- Sys.getenv("DNMB_PROMOTECH_MODEL_DIR", unset = "")
model <- Sys.getenv("DNMB_PROMOTECH_MODEL", unset = "RF-HOT")
model_base_url <- Sys.getenv(
  "DNMB_PROMOTECH_MODEL_BASE_URL",
  unset = "https://www.cs.mun.ca/~lourdes/public/PromoTech_models"
)
download_model <- identical(tolower(Sys.getenv("DNMB_PROMOTECH_DOWNLOAD_MODEL", unset = "false")), "true")
force <- identical(tolower(Sys.getenv("DNMB_PROMOTECH_FORCE_INSTALL", unset = "false")), "true")

asset_urls <- list()
if (nzchar(model_dir)) {
  asset_urls$models_dir <- model_dir
}

repo_source <- if (nzchar(repo_dir)) repo_dir else DNMB:::.dnmb_promotech_default_repo_url()

result <- DNMB:::dnmb_promotech_install_module(
  version = DNMB:::.dnmb_promotech_default_version(),
  cache_root = cache_root,
  install = TRUE,
  repo_url = repo_source,
  asset_urls = asset_urls,
  model = model,
  download_model = download_model,
  model_base_url = model_base_url,
  force = force
)

print(result$status, n = Inf, width = Inf)

if (!isTRUE(result$ok)) {
  quit(status = 1L)
}

invisible(result)
