FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive \
    DNMB_CACHE_ROOT=/opt/dnmb-cache \
    DNMB_DEFENSEFINDER_CASFINDER_DIR=/root/.macsyfinder/models/CasFinder \
    DNMB_DEFENSEFINDER_REPO_DIR=/opt/vendor/defense-finder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    git \
    python3 \
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
COPY . /opt/DNMB

RUN python3 -m pip install --no-cache-dir 'macsyfinder==2.1.4' && \
    macsydata install -u 'CasFinder==3.1.0'

RUN git clone --branch v2.0.1 --depth 1 https://github.com/mdmparis/defense-finder.git /opt/vendor/defense-finder

RUN R -q -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages(c('remotes', 'BiocManager')); BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask = FALSE, update = FALSE); remotes::install_local('/opt/DNMB', dependencies = TRUE, upgrade = 'never')"

RUN Rscript /opt/DNMB/inst/scripts/prewarm_defensefinder_cache.R

WORKDIR /work

CMD ["R"]
