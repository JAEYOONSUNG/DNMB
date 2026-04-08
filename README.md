# DNMB Package

## Installation

To install the DNMB package and its dependencies, please follow these steps:

### Check if devtools is installed, and install it if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

### Step 1: Install Bioconductor dependencies

First, install the necessary Bioconductor packages:

```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ComplexHeatmap", "Peptides"))


### Step 2: Install CRAN dependencies

Next, install the required CRAN packages. Some packages may not be installed automatically, so you can install them manually:

# sometimes it can not installed automatically followed
# Install CRAN packages
install.packages(c("qdap", "seqinr", "circlize", "splitstackshape"))

#remotes::install_github("trinker/qdap")

### Step 1: Install Java Development Kit (JDK)

The `rJava` package requires Java to be installed. Please install the Java Development Kit (JDK) from the following sources:

- **Windows/Mac**: [Oracle JDK Downloads](https://www.oracle.com/java/technologies/javase-jdk11-downloads.html)
- **Linux**: Use your package manager to install OpenJDK. For example:
  ```bash 
  sudo apt install openjdk-11-jdk

After installation, ensure that the JAVA_HOME environment variable is set correctly.

    Windows:

        1.    Go to System Properties -> Environment Variables -> System Variables -> New
        2.    Variable name: JAVA_HOME
        3.    Variable value: C:\Program Files\Java\jdk-11.0.1 (your JDK installation path)

    Mac/Linux:
        Add the following to your .bash_profile or .zshrc:
        ```bash
        export JAVA_HOME=$(/usr/libexec/java_home)

Step 2: Install rJava

Once Java is installed and configured, you can install the rJava package:
    Windows:
        ```bash
        Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-11.0.1")
        ```r
        install.packages("rJava")
    Mac/Linux:
        ```r
        install.packages("rJava")

### Step 3: Install DNMB package
devtools::install_github("JAEYOONSUNG/DNMB")

## MEROPS module

The DNMB `MEROPS` module can now run directly through the package module APIs.

### Requirements

- `makeblastdb` and `blastp` must be available in `PATH`
- the current implementation supports both legacy BLAST+ `2.2.31+` and newer BLAST+ releases
- MEROPS assets can be supplied as local files for offline or fixture-based runs

DNMB automatically adapts the `makeblastdb` arguments to the detected BLAST+ version. Older BLAST+ releases that do not support `-blastdb_version` are handled automatically.

When `module_version = "current"` and official MEROPS URLs are used, DNMB also records remote asset metadata from the `current_release` files and refreshes the cached MEROPS assets when the upstream files change.

### Local asset example

```r
merops_run <- DNMB::run_module_set(
  module_MEROPS = TRUE,
  genbank_table = genbank_table,
  module_install = TRUE,
  module_asset_urls = c(
    "pepunit.lib" = "/path/to/pepunit.lib",
    "dnld_list.txt" = "/path/to/dnld_list.txt"
  ),
  merge = TRUE
)
```

The module cache is stored under `DNMB::dnmb_db_home()` unless `module_cache_root` is supplied.
By default DNMB resolves the cache root in this order: `DNMB_CACHE_ROOT`, `~/.dnmb-cache`, then `tools::R_user_dir("DNMB", which = "cache")`.

### `run_module_set()` example

```r
merops_runs <- DNMB::run_module_set(
  module_MEROPS = TRUE,
  genbank_table = genbank_table,
  module_install = TRUE,
  module_cpu = 4
)

merged_merops <- DNMB::run_module_set(
  module_MEROPS = TRUE,
  genbank_table = genbank_table,
  module_install = TRUE,
  module_cpu = 4,
  merge = TRUE
)
```

### `run_DNMB()` example

```r
DNMB::run_DNMB(
  module_MEROPS = TRUE
)
```

When `module_MEROPS = TRUE`, the final DNMB table receives MEROPS-prefixed columns such as `MEROPS_family_id` and `MEROPS_hit_label`.

## dbCAN module

The DNMB `dbCAN` module now runs directly from the package by using the official dbCAN HMM database and `hmmsearch`.

### Requirements

- `hmmsearch`, `hmmpress`, and `hmmbuild` from HMMER must be available in `PATH`
- `run_dbcan`, `diamond`, and `makeblastdb` are used automatically when gene coordinates are available so DNMB can run `CGCFinder` and `dbCAN-PUL` substrate prediction
- when `module_version = "current"`, DNMB resolves the latest dbCAN HMMdb release from the official dbCAN site and downloads the current `dbCAN.txt` plus `fam-substrate-mapping.tsv`

The default dbCAN module path auto-installs its HMM database under `DNMB::dnmb_db_home()` if the cache is missing. The current installer uses the official `pro.unl.edu` dbCAN host.

When `run_dbcan` is available and the active gene table includes `contig`, `start`, `end`, and `direction`, DNMB automatically switches from the simple HMMER-only mode to the standalone `run_dbcan` workflow with:

- HMMER CAZyme annotation
- `CGCFinder`
- `dbCAN-PUL` substrate prediction

The generated dbCAN standalone outputs remain in the module output folder alongside the DNMB-wide table.

### `run_module_set()` example

```r
dbcan_runs <- DNMB::run_module_set(
  module_dbCAN = TRUE,
  genbank_table = genbank_table,
  module_install = TRUE,
  module_cpu = 4
)
```

### `run_DNMB()` example

```r
DNMB::run_DNMB(
  module_dbCAN = TRUE
)
```

When `module_dbCAN = TRUE`, the final DNMB table receives dbCAN-prefixed columns such as `dbCAN_family_id` and the legacy alias `dbCAN_dbcan_hit`.

## Mobileome / IS impact pipeline

The default mobileome pipeline now runs in `hybrid` mode:

- annotation-derived IS/mobile calls from the GenBank feature table
- ISEScan-like sequence evidence from `Prodigal + hmmsearch/phmmer + family-aware TIR/TSD heuristics`
- genome-native essentiality/redundancy scoring
- genome-native target-site modeling from observed TSD/context patterns
- optional ISfinder-style reference refinement when a curated FASTA or BLAST DB is supplied as a low-weight refinement layer
- optional comparative mode using related GenBank genomes for occupied/empty hotspot and chronology inference

### Current project status (2026-03-09)

Implemented and wired into the DNMB package:

- `run_DNMB_mobileome()` as the package entry point for single-GenBank mobileome analysis
- hybrid IS detection from annotation-derived calls plus ISEScan-like sequence evidence
- genome-native recognition modeling from observed TSD/context patterns
- essentiality-aware mutation point, SV scenario, and phenotype prioritization
- comparative mode from related GenBank genomes with occupied/empty hotspot and chronology inference
- integrated output tables including `mobileome_master_table.tsv` and `mobileome_compact_table.tsv`
- package plot helpers for integrated mobileome overview and comparative hotspot visualization

Current module layout:

- `R/Mobileome_pipeline.R`: end-to-end orchestration, exports, file writing
- `R/Mobileome_sequence_engine.R`: sequence-based IS candidate detection and merge logic
- `R/Mobileome_variant_engine.R`: target model, target site, variant catalog, evidence integration
- `R/Mobileome_comparative.R`: related-genome comparison, hotspot and chronology inference
- `R/Mobileome_plot.R`: integrated overview plot and comparative plotting helpers

Validated so far:

- example `data/examples/annot.gbk` runs end-to-end in `hybrid` mode
- `mobileome_compact_table.tsv` reports both detailed evidence labels (`ANN`, `ANN+SEQ`, `ANN+SEQ+NAT`, `SEQ`, `SEQ+NAT`) and simplified support groups (`NATIVE`, `OVERLAP`, `ANN`, `SEQ`)
- integrated overview plot now combines:
  - observed merged IS calls
  - targetable-region underlay
  - recognition panel
  - per-family evidence count panel

Current known gaps:

- some IS families still fall back to context-only or weak motif models
- RNA-guided/hairpin-driven families such as `IS110` and `IS200/IS605` still need a stronger native recognition layer
- `ggseqlogo` emits deprecation warnings from the upstream package, although rendering is currently successful
- comparative mode has been validated on local/example runs, but accuracy depends strongly on the quality and proximity of related input genomes

For a single GenBank-driven mobileome analysis, use:

```r
DNMB::run_DNMB_mobileome(
  genbank = "path/to/input.gbk",
  output_dir = "path/to/output",
  detection_mode = "hybrid",
  site_scan_strategy = "balanced",
  related_genbanks = c("path/to/related1.gbk", "path/to/related2.gbk")
)
```

From the source tree you can also run:

```bash
Rscript inst/scripts/run_mobileome_pipeline.R \
  --genbank path/to/input.gbk \
  --output path/to/output
```

Comparative mode from the CLI:

```bash
Rscript inst/scripts/run_mobileome_pipeline.R \
  --genbank path/to/input.gbk \
  --output path/to/output \
  --related-genbanks path/to/related1.gbk,path/to/related2.gbk,path/to/related3.gbk \
  --related-metadata path/to/related_metadata.tsv \
  --min-related-completeness 95 \
  --min-related-ani 95 \
  --related-detection-mode annotation
```

Comparative outputs are written under `output/comparative/` and include:

- `comparative_loci.tsv`
- `occupied_empty_matrix.tsv`
- `comparative_hotspots.tsv`
- `event_chronology.tsv`
- `family_locus_master.tsv`

`related_metadata.tsv` can follow the BPGA/NCBI-style frame used in your other
scripts. DNMB will directly recognize columns such as:

- `RefSeq assembly accession`
- `CheckM completeness`
- `ANI`
- `SOURCE`
- `DEFINITION`

The integrated focal result directory now also includes:

- `mobileome_master_table.tsv`
- `mobileome_compact_table.tsv` with detailed integrated evidence labels (`ANN`, `ANN+SEQ`, `ANN+SEQ+NAT`, `SEQ`, `SEQ+NAT`) and simplified support groups (`NATIVE`, `OVERLAP`, `ANN`, `SEQ`)
- `is_target_models.tsv`
- `is_targetable_regions.tsv`
- `is_targetable_regions.bed`

Plot helpers:

```r
DNMB::plot_DNMB_mobileome_distribution("path/to/output")
DNMB::plot_DNMB_mobileome_comparative(result$comparative, "path/to/output")
```

If you have a curated IS reference FASTA or BLAST DB compatible with ISfinder-style family refinement:

```bash
Rscript inst/scripts/run_mobileome_pipeline.R \
  --genbank path/to/input.gbk \
  --output path/to/output \
  --isfinder-fasta path/to/is_reference.fna
```
