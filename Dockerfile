FROM ghcr.io/jaeyoonsung/dnmbsuite:latest AS dnmb_cache

FROM rocker/r-ver:4.4.1

ARG DNMB_PROMOTECH_DOWNLOAD_MODEL=false
ARG DNMB_PROMOTECH_MODEL=RF-HOT
ARG DNMB_PROMOTECH_MODEL_BASE_URL=https://www.cs.mun.ca/~lourdes/public/PromoTech_models
ARG VIENNARNA_VERSION=2.7.2

ENV DEBIAN_FRONTEND=noninteractive \
    DNMB_CACHE_ROOT=/opt/dnmb-cache \
    DNMB_DEFENSEFINDER_CASFINDER_DIR=/root/.macsyfinder/models/CasFinder \
    DNMB_DEFENSEFINDER_REPO_DIR=/opt/vendor/defense-finder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    diamond-aligner \
    git \
    libbio-perl-perl \
    libdbd-sqlite3-perl \
    libdbi-perl \
    libwww-perl \
    ncbi-blast+ \
    prodigal \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    libbz2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libicu-dev \
    libjpeg-dev \
    liblzma-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libuv1-dev \
    libxml2-dev \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

RUN curl -L http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz -o /tmp/hmmer-3.4.tar.gz && \
    tar -xzf /tmp/hmmer-3.4.tar.gz -C /tmp && \
    cd /tmp/hmmer-3.4 && \
    ./configure --prefix /usr/local && \
    make -j"$(nproc)" && \
    make install && \
    (cd easel && make install) && \
    rm -rf /tmp/hmmer-3.4 /tmp/hmmer-3.4.tar.gz

WORKDIR /opt/DNMB

RUN python3 -m pip install --no-cache-dir 'macsyfinder==2.1.4' && \
    macsydata install -u 'CasFinder==3.1.0'

RUN python3 -m pip install --no-cache-dir \
    'numpy<2' \
    'pandas<2' \
    joblib \
    biopython \
    progressbar2 \
    'scikit-learn<1.3'

RUN git clone --branch v2.0.1 --depth 1 https://github.com/mdmparis/defense-finder.git /opt/vendor/defense-finder

# Conda env at /opt/biotools holds the bioinformatics binaries that
# eggnog-mapper / padloc / phispy / GapMind expect on PATH. The cache lives
# under DNMB_CACHE_ROOT (host bind-mounted) so model/HMM data is reused
# across runs; only the executables need to be baked into the image. Once
# the env is set up, /opt/biotools/bin is prepended to PATH so module
# detect_runner functions find the tools.
RUN curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname -s)-$(uname -m).sh" -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/miniforge && \
    rm /tmp/miniforge.sh && \
    /opt/miniforge/bin/conda create -y -p /opt/biotools \
        -c bioconda -c conda-forge \
        python=3.10 \
        eggnog-mapper \
        phispy \
        entrez-direct && \
    /opt/miniforge/bin/conda create -y -p /opt/biotools/envs/padloc \
        -c bioconda -c conda-forge \
        padloc && \
    /opt/miniforge/bin/conda clean -afy && \
    rm -rf /opt/miniforge && \
    # PADLOC's bin/padloc bash script unconditionally mkdir's bin/../data
    # before honoring the --data flag, so non-root users (--user) hit
    # "Permission denied" on a fresh run. Pre-create the dir as world-
    # writable; the actual HMM data is served from the host-mounted
    # /opt/dnmb-cache/db_modules/padloc/current/data via the --data flag.
    mkdir -p /opt/biotools/envs/padloc/data && \
    chmod -R a+rwX /opt/biotools/envs/padloc/data

ENV PATH="/opt/biotools/bin:${PATH}"

RUN mkdir -p /tmp/dnmb-deps
COPY DESCRIPTION /tmp/dnmb-deps/DESCRIPTION

# cache-bust 2026-04-28: previous buildkit cache layers dropped the
# libssl-dev / libcurl4-openssl-dev provided headers, leaving the
# openssl R package compile broken. Tagged the RUN command with
# explicit verification so a regression is loud, not silent.
RUN R -q -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages(c('remotes', 'BiocManager')); BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask = FALSE, update = FALSE); remotes::install_deps('/tmp/dnmb-deps', dependencies = TRUE, upgrade = 'never'); stopifnot(requireNamespace('Biostrings', quietly = TRUE), requireNamespace('ComplexHeatmap', quietly = TRUE))"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      clustalw \
      emboss \
      libbio-perl-run-perl \
      libbio-tools-run-alignment-clustalw-perl \
      libdate-calc-perl \
      libgsl-dev \
      libjson-parse-perl \
      liblapack-dev \
      liblapacke-dev \
      libmpfr-dev \
      muscle \
      pkg-config \
      vmatch && \
    ln -sf /usr/bin/vmatch /usr/local/bin/vmatch2 && \
    ln -sf /usr/bin/mkvtree /usr/local/bin/mkvtree2 && \
    ln -sf /usr/bin/vsubseqselect /usr/local/bin/vsubseqselect2 && \
    rm -rf /var/lib/apt/lists/*

RUN curl -L "https://github.com/ViennaRNA/ViennaRNA/releases/download/v${VIENNARNA_VERSION}/ViennaRNA-${VIENNARNA_VERSION}.tar.gz" -o /tmp/ViennaRNA.tar.gz && \
    tar -xzf /tmp/ViennaRNA.tar.gz -C /tmp && \
    cd "/tmp/ViennaRNA-${VIENNARNA_VERSION}" && \
    ./configure --prefix=/usr/local --disable-lto --without-swig --without-perl --without-python --without-doc --without-forester --without-kinfold --without-rnalocmin && \
    make -j"$(nproc)" && \
    make install && \
    ldconfig && \
    RNAfold --version && \
    rm -rf "/tmp/ViennaRNA-${VIENNARNA_VERSION}" /tmp/ViennaRNA.tar.gz

RUN mkdir -p /opt/macsyfinder/models && \
    cp -a /root/.macsyfinder/models/CasFinder /opt/macsyfinder/models/CasFinder && \
    chmod -R a+rX /opt/macsyfinder

ENV DNMB_DEFENSEFINDER_CASFINDER_DIR=/opt/macsyfinder/models/CasFinder \
    HOME=/tmp

RUN git clone --depth 1 https://github.com/HaidYi/acrfinder.git /opt/vendor/acrfinder
RUN git clone --depth 1 --filter=blob:none --sparse https://github.com/BioinformaticsLabAtMUN/Promotech.git /opt/vendor/promotech && \
    git -C /opt/vendor/promotech sparse-checkout set promotech.py genome core models sequences

COPY . /opt/DNMB

RUN R -q -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); remotes::install_local('/opt/DNMB', dependencies = FALSE, upgrade = 'never'); stopifnot(requireNamespace('DNMB', quietly = TRUE))"

# Pre-install DefenseViz at image build time so REBASEfinder works under
# `--user $(id -u):$(id -g)` (non-root users cannot write to
# /usr/local/lib/R/site-library at runtime). Verbose install + explicit
# dependency resolution so any missing CRAN/Bioconductor dep surfaces in
# the build log instead of being swallowed by quiet=TRUE.
RUN R -q -e "options(repos = c(CRAN = 'https://cloud.r-project.org', BioCsoft = 'https://bioconductor.org/packages/release/bioc')); remotes::install_github('JAEYOONSUNG/DefenseViz', upgrade = 'never', quiet = FALSE); stopifnot(requireNamespace('DefenseViz', quietly = TRUE)); stopifnot(requireNamespace('DNMB', quietly = TRUE))"

COPY --from=dnmb_cache /opt/dnmb-cache/db_modules/defensefinder/current /opt/dnmb-cache/db_modules/defensefinder/current

RUN Rscript /opt/DNMB/inst/scripts/prewarm_defensefinder_cache.R
RUN Rscript -e 'cache_root <- Sys.getenv("DNMB_CACHE_ROOT"); result <- DNMB:::dnmb_acrfinder_install_module(cache_root = cache_root, install = TRUE, repo_url = "/opt/vendor/acrfinder", force = TRUE); print(result$status); stopifnot(isTRUE(result$ok))'
RUN DNMB_PROMOTECH_REPO_DIR=/opt/vendor/promotech \
    DNMB_PROMOTECH_DOWNLOAD_MODEL="${DNMB_PROMOTECH_DOWNLOAD_MODEL}" \
    DNMB_PROMOTECH_MODEL="${DNMB_PROMOTECH_MODEL}" \
    DNMB_PROMOTECH_MODEL_BASE_URL="${DNMB_PROMOTECH_MODEL_BASE_URL}" \
    Rscript /opt/DNMB/inst/scripts/prewarm_promotech_cache.R

WORKDIR /data

CMD ["Rscript", "-e", "suppressWarnings(suppressPackageStartupMessages(library(DNMB))); run_DNMB()"]
