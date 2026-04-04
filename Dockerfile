## ─────────────────────────────────────────────────────────────
## DNMB – Full Dockerfile (All Modules + InterProScan)
## ─────────────────────────────────────────────────────────────
## Build:  docker build --build-arg DNMB_SOURCE_HASH="$(git rev-parse --short HEAD)" -t dnmb:latest .
## Run:    docker run --rm -v $(pwd):/data -v ~/.dnmb-cache:/opt/dnmb/cache dnmb:latest
## InterProScan: auto-downloaded to cache on first module_InterProScan=TRUE
## ─────────────────────────────────────────────────────────────

FROM --platform=linux/amd64 rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive
ENV DNMB_CACHE_ROOT=/opt/dnmb/cache
# InterProScan: auto-downloaded to cache on first use
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# ═══════════════════════════════════════════════════════════════
# 1. System packages
# ═══════════════════════════════════════════════════════════════
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
    libharfbuzz-dev libfribidi-dev \
    zlib1g-dev libbz2-dev liblzma-dev libpcre2-dev libicu-dev libgit2-dev \
    cmake pkg-config curl wget git procps ca-certificates locales \
    default-jdk \
    libwebp-dev \
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ═══════════════════════════════════════════════════════════════
# 2. Conda → /opt/biotools (isolated from system R)
# ═══════════════════════════════════════════════════════════════
RUN curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname -s)-$(uname -m).sh \
    -o /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /opt/miniforge \
    && rm /tmp/miniforge.sh

RUN /opt/miniforge/bin/conda create -y -p /opt/biotools \
    -c bioconda -c conda-forge \
    python=3.12 \
    hmmer blast prodigal diamond \
    eggnog-mapper \
    entrez-direct skani fastani \
    phispy \
    perl-dbi perl-lwp-simple perl-dbd-sqlite \
    && /opt/miniforge/bin/conda clean -afy

ENV PATH="/opt/biotools/bin:${PATH}"

# ═══════════════════════════════════════════════════════════════
# 3. CLEAN module env (PyTorch + ESM, pre-built)
# ═══════════════════════════════════════════════════════════════
RUN mkdir -p ${DNMB_CACHE_ROOT}/db_modules/clean/split100 \
    && /opt/biotools/bin/python -m venv ${DNMB_CACHE_ROOT}/db_modules/clean/split100/conda_env \
    && ${DNMB_CACHE_ROOT}/db_modules/clean/split100/conda_env/bin/pip install --no-cache-dir \
       torch --index-url https://download.pytorch.org/whl/cpu \
    && ${DNMB_CACHE_ROOT}/db_modules/clean/split100/conda_env/bin/pip install --no-cache-dir \
       "fair-esm==2.0.0" \
       "pandas>=1.4" \
       "scikit-learn>=1.2" \
       "scipy>=1.7" \
       "matplotlib>=3.7" \
       "tqdm>=4.64" \
       gdown

# VirSorter2 & DefenseFinder: auto-installed by DNMB on first run
RUN rm -rf /opt/miniforge

# ═══════════════════════════════════════════════════════════════
# 4. InterProScan (optional — x86_64 only)
# ═══════════════════════════════════════════════════════════════
# InterProScan: auto-downloaded to cache on first use (not baked into image)
# Keeps image small (~5GB) and cache reusable across rebuilds

# ═══════════════════════════════════════════════════════════════
# 5. R packages (CRAN + Bioconductor) — ALL dependencies
# ═══════════════════════════════════════════════════════════════
RUN R -e ' \
    options(Ncpus = parallel::detectCores()); \
    install.packages("BiocManager", repos = "https://cloud.r-project.org"); \
    BiocManager::install(c("Biostrings", "ComplexHeatmap"), ask = FALSE, update = FALSE); \
    install.packages(c( \
      "dplyr", "plyr", "tidyr", "data.table", "tibble", "reshape2", \
      "readr", "openxlsx", "seqinr", "stringr", "jsonlite", "gtools", \
      "ggplot2", "cowplot", "gggenes", "ggrepel", "ggtext", "ggseqlogo", \
      "ggforce", "gridExtra", "scales", "Peptides", "qdap", "circlize", \
      "devtools", "testthat", "tidyverse", \
      "ggnewscale", "patchwork", "gridBase", "gtable", "colorspace" \
    ), repos = "https://cloud.r-project.org"); \
'

# REBASEfinder dependency
RUN R -e 'devtools::install_github("JAEYOONSUNG/DefenseViz", quiet = TRUE)'

# ═══════════════════════════════════════════════════════════════
# 6. Install DNMB package (cache-bust: changes to R/ trigger rebuild)
# ═══════════════════════════════════════════════════════════════
ARG DNMB_SOURCE_HASH=dev
COPY R/ /tmp/DNMB/R/
COPY DESCRIPTION NAMESPACE /tmp/DNMB/
COPY man/ /tmp/DNMB/man/
COPY inst/ /tmp/DNMB/inst/
COPY docker/ /tmp/DNMB/docker/
RUN echo "$DNMB_SOURCE_HASH" > /tmp/DNMB/.source-hash \
    && R -e 'devtools::install("/tmp/DNMB", dependencies = FALSE)' \
    && rm -rf /tmp/DNMB

# ═══════════════════════════════════════════════════════════════
# 7. Runtime
# ═══════════════════════════════════════════════════════════════
RUN mkdir -p /data /results ${DNMB_CACHE_ROOT}

WORKDIR /data

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["R"]
