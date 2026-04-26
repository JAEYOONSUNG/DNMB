# DNMB
[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Fig  1](https://github.com/user-attachments/assets/ef35c8cc-5e75-40ea-a13e-42022ccbc0ca)

The domestication of non-model bacteria (DNMB) presents a unique set of challenges and opportunities. To address these, we introduce a novel approach—the DNMB pipeline—a comprehensive genomics analysis tool tailored specifically for non-model bacterial species. Unlike traditional model organisms, non-model bacteria often lack well-defined genetic resources and established analytical pipelines. Consequently, researchers face hurdles in elucidating the genetic basis underlying phenotypic traits relevant to domestication efforts. Our pipeline integrates multiple scripts and computational tools to streamline various stages of genomic analysis, from making comprehensive functional annotation in user-friendly table format to genomic features including codon usage and ribosomal binding site preference and distance. 
Herein, we provide a detailed overview of the DNMB pipeline, highlighting its key components and functionalities. Additionally, we demonstrate its utility through a case study involving the domestication of a non-model bacterial strain. The DNMB pipeline not only accelerates genomic analysis but also enhances our understanding of non-model bacterial physiology, thereby facilitating the exploitation of microbial diversity for biotechnological applications.

## Key Features
Diverse Scripts: Our collection includes a range of scripts, each developed to address unique challenges in systems biology research.

Curated Database: Access to a comprehensive database, meticulously compiled to support and enhance your research projects. 
                  We provide a package to facilitate tabulating data from various databases including [REBASE](http://rebase.neb.com), [MEROPS](https://www.ebi.ac.uk/merops/download_list.shtml), and [CAZy_dbCAN3](). The tables, which can be easily converted into                   FASTA format, allow for seamless integration with various sequence analysis tools, 
                  providing flexibility and ease of use for researchers. enabling users to extract desired information using various sequence analysis tools, including BLAST.

User-Friendly Documentation: Detailed documentation is available to guide you through the installation, setup, and utilization of both the scripts and the database.

### Algorithms for analysis

#### Ribosomal binding site: 
The RBS algorithm begins by extracting the last 9 nucleotides of the full-length 16S ribosomal RNA. It then identifies all instances of the reverse complement sequence within the genome, allowing up to 2 mismatches. From these, it selects sequences based on their proximity to an annotated start codon—typically within a range of 1 to 10 nucleotides (this range is the default setting). Finally, it statistically calculates the preference for selected RBS sequences and the distance (spacer) between the start codon and the RBS sequence.

#### Codon usage:
Codon usage analysis quantifies the total count of amino acids and nucleotides (codons) across the entire coding sequence. This analysis is performed using the uco function from the seqinr package, calculating three key indices: eff for codon counts, freq for relative frequencies of codons, and rscu for the Relative Synonymous Codon Usage index. The terms "eff", "freq", and "rscu" are equivalent to "R0", "R1", and "R3", respectively, as defined in Suzuki et al. (2005) under the section "2.2 Normalization of codon usage data". Furthermore, "eff" and "rscu" correspond to "AF" and "RSCU", respectively, in Suzuki et al. (2008) "2.2. Definitions of codon usage data".


## Getting Started:
To begin using our resources, please follow the steps outlined in our documentation. 
Whether you're looking to integrate our scripts into your existing projects or explore our database for new insights, we've provided all the necessary instructions to get you started.

## Installation
### Requirements

The DNMB is supported for macOS, Linux and Windows machines, which can provide an environment for using R.
It requires R version >=4.2.1 for release, and R version >=4.3 for devel.

One of the third-party functionalities is not available for Windows and MacOS machines (InterProScan).

The [EggNOG-mapper webserver](http://eggnog-mapper.embl.de), allows users to input sequences in FASTA format based on locus_tag identifiers and receive results in either XLSX or CSV format. Additionally, the standalone version available on GitHub is compatible with DNMB.

InterProScan requires a Linux operating system. Without access to Linux, you can proceed with the analysis up to Eggnog-mapper in the annotation stage, but you won't be able to obtain information about motif analysis.


To download and install R, see the [R-project website](https://www.r-project.org/).

To download and install InterProScan, see the [InterProScan github](https://github.com/ebi-pf-team/interproscan).

To download and install EggNOG-mapper, see the [EggNOG-mapper github](https://github.com/eggnogdb/eggnog-mapper).


#### Warning
The basic file for genomic analysis, known as a GenBank file, requires both sequence and annotation in full-format files such as gbff, gb, or gbk. Additionally, GenBank prefers a format based on the GeneMarkS2+ pipeline, and using a different annotation pipeline to obtain GenBank files may lead to errors.



## Anaylsis flow
![Fig S1](https://github.com/user-attachments/assets/7fb20e1c-1d72-463b-9ff0-caeb82b7cf2f)

## Prerequisites
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ComplexHeatmap", "Peptides"))

install.packages(c("qdap", "seqinr", "stringr", "stringi", "splitstackshape", "gtools", "ggplot2", "ggseqlogo", "circlize", "grid", "gridExtra","plyr", "dplyr", "tidyr", "readr", "reshape2", "data.table", "tibble", "qdap", "openxlsx"))
```
- **Note:** If you encounter issues installing the qdap package, try installing it with the following command:
```r
install.packages("qdap", INSTALL_opts = "--no-multiarch")
```

   
- **Note:** Java and the rJava package must be installed and configured to enable .xlsx output using this package.
1. Install Java Development Kit (JDK):
        Download and install the appropriate JDK for your operating system from the Oracle website or OpenJDK.
2. Install rJava Package in R:
```r
install.packages("rJava")
library(rJava)
```

3. Set $JAVA_HOME Path:

You need to set the environment variable JAVA_HOME to point to the location of your JDK installation.
	
 •	On Windows:

1.	Install Java jdk (https://www.oracle.com/kr/java/technologies/downloads/)
2.	Check the “System Variables,” :
	•	Variable name: JAVA_HOME
	•	Variable value: The path to your JDK installation 

```bash
# Print the current value of the JAVA_HOME environment variable.
echo %JAVA_HOME%  #(e.g., C:\Program Files\Java\jdk-18)

# Set the JAVA_HOME environment variable to point to the Java Development Kit (JDK) installation.
# Replace [version] with your installed JDK version (e.g., jdk-18).
setx JAVA_HOME "C:\Program Files\Java\jdk[version]"

# Update the system PATH to include the bin directory of the JDK.
setx PATH "%JAVA_HOME%\bin;%path%"

# Check the installed Java version to confirm that the correct version is being used.
JAVA -version
```
3.	Restart R or RStudio.
	
 •	On macOS/Linux:
-If Xcode is not installed, you may encounter compiler issues during package installation. To resolve this, install Xcode from the App Store.
Add this line to your .bash_profile or .bashrc (depending on the shell):
```bash
# Navigate to your Java installation directory to check available Java versions
/Library/Java/JavaVirtualMachines/[my_java_folder]/Contents/Home # check my java list

# Open your .bash_profile (or .bashrc) file for editing
vi ~/.bash_profile # edit bash profile

# Press 'i' to enter insert mode in the vi editor
i # insert mode

# Add or update the JAVA_HOME environment variable with the path to your Java installation
export JAVA_HOME=/Library/Java/JavaVirtualMachines/[my_java_folder]/Contents/Home
# Add Java's bin directory to the system PATH variable so that Java commands can be run from the terminal
export PATH=${PATH}:$JAVA_HOME/bin

# Save the changes and exit the vi editor. ":wq!" means "write" (save) and "quit" (exit) forcefully
: # activate command line
wq! # save

# Apply the changes made to the .bash_profile or .bashrc immediately (without needing to restart the terminal)
source ~/.bash_profile  ## or ~/.bashrc #apply changes

# Verify that JAVA_HOME is set correctly by printing its value
echo $JAVA_HOME # validation
```

## Install DNMB R package
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("JAEYOONSUNG/DNMB")
```

# Quick start
## Run DNMB analysis
```r
setwd([GenBank directory]) # Set the working directory to the location where your GenBank files are stored.
library(DNMB)
run_DNMB()
```

### GPU-gated defaults (CLEAN and PIDE)

`run_DNMB()` probes for an NVIDIA GPU via `nvidia-smi -L` at call time.
`module_CLEAN` and `module_PIDE` default to `TRUE` only when a CUDA
device is detected; they default to `FALSE` otherwise. This keeps
CPU-only runs fast (both modules use large neural models — CLEAN's
LayerNormNet and PIDE's ESM-650M — and are ~50–100× slower without a
GPU).

- Force enable: `run_DNMB(module_CLEAN = TRUE, module_PIDE = TRUE)` or
  set `Sys.setenv(DNMB_CUDA = "1")` before the call.
- Force disable: pass `module_CLEAN = FALSE` / `module_PIDE = FALSE`
  explicitly, or set `DNMB_CUDA=0`.
- In the DNMBsuite Docker wrapper, the same probe runs on the host and
  automatically attaches `--gpus all` when CUDA is present.

### With comparative heatmaps across sibling genomes
Passing `comparative = TRUE` runs the single-genome pipeline as usual and,
at the end, renders the full suite of comparative heatmaps across every
sibling folder that holds a GenBank file. By default the parent directory
of `getwd()` is scanned; pass `comparative_data_root` to point elsewhere.

```r
setwd([per-genome directory])
library(DNMB)
run_DNMB(comparative = TRUE)
# or: run_DNMB(comparative = TRUE, comparative_data_root = "/path/to/parent")
```

Comparative stage renders 14 heatmaps into `<data_root>/comparative/`:

| Module | Plotter | Subtype axis |
|---|---|---|
| DefenseFinder | `dnmb_plot_comparative_defensefinder` | system |
| PADLOC | `dnmb_plot_comparative_padloc` | system |
| DefensePredictor | `dnmb_plot_comparative_defensepredictor` | system |
| REBASEfinder | `dnmb_plot_comparative_rebasefinder` | enzyme type |
| MEROPS | `dnmb_plot_comparative_merops` | family (C26, S8, …) |
| MEROPS | `dnmb_plot_comparative_merops_catalytic` | catalytic type (Cysteine, Serine, …) |
| dbCAN | `dnmb_plot_comparative_dbcan` | class (GH, GT, PL, …) |
| dbCAN | `dnmb_plot_comparative_dbcan_family` | family (GH13, GT2, …) |
| CGC | `dnmb_plot_comparative_cgc` | signature mix (CAZyme+TC+TF, …) |
| CGC | `dnmb_plot_comparative_cgc_substrate` | predicted substrate |
| PAZy | `dnmb_plot_comparative_pazy` | family |
| PhiSpy | `dnmb_plot_comparative_phispy` | region size bucket |
| VirSorter2 | `dnmb_plot_comparative_virsorter2` | max_score_group (dsDNAphage, ssDNA, …) |
| PIDE | `dnmb_plot_comparative_pide` | region size bucket |

Each plotter auto-runs its own module on any genome that has not been
analyzed yet, so a fresh sibling folder just needs a GenBank file.

---
**DefenseFinderViz** (Optional)

![DefenseFinder_Heatmap](https://github.com/user-attachments/assets/68e953c1-d568-4b56-89f2-1107a60a6f6c)
```r
DefenseFinder_Heatmap()
```
- **Note:** [Strain of interest].defense_finder_systems.tsv output are used for merging data. GenBank’s SOURCE field is used for extracting names.
- **Note:** protien coding sequence (.faa) output was used for defense-finder analysis (https://github.com/mdmparis/defense-finder)

---
**Comparative per-module heatmaps across genomes** (Optional)

The easiest entry point is `run_DNMB(comparative = TRUE)` from any
per-genome folder — the per-genome analysis runs as usual and the full
comparative suite renders against the parent directory at the end. The
individual plotters below are useful when you want to render only a
subset, override colors, or point at a non-sibling parent directory.

Point each plotter at a parent directory containing one subfolder per genome.
Every subfolder that holds a GenBank file (`*.gbff` / `*.gbk` / `*.gb`) is
treated as a genome. Genomes missing the relevant module output are
analyzed on the fly (`auto_run_missing = TRUE`, default) — each plotter
triggers only its own module via `run_module_set(db = ...)`, not the
full DNMB pipeline. Genomes whose module has already run are read from
disk; genomes that truly have no hits still render as empty rows so
"analyzed, empty" is distinguishable from "not yet analyzed".

```r
library(DNMB)

data_root <- "/path/to/parent-dir-of-genome-folders"

# Defense-system heatmaps (purple palette)
dnmb_plot_comparative_defensefinder(data_root)    # DefenseFinder
dnmb_plot_comparative_padloc(data_root)           # PADLOC
dnmb_plot_comparative_defensepredictor(data_root) # DefensePredictor
dnmb_plot_comparative_rebasefinder(data_root)     # REBASEfinder

# Enzyme / CAZyme heatmaps (module-specific palettes)
dnmb_plot_comparative_merops(data_root)             # MEROPS family (C26, S8, …)
dnmb_plot_comparative_merops_catalytic(data_root)   # MEROPS catalytic type (Cysteine, Serine, …)
dnmb_plot_comparative_dbcan(data_root)              # dbCAN class (GH, GT, PL, …)
dnmb_plot_comparative_dbcan_family(data_root)       # dbCAN family (GH13, GT2, …)
dnmb_plot_comparative_cgc(data_root)                # CGC signature mix (CAZyme+TC+TF, …)
dnmb_plot_comparative_cgc_substrate(data_root)      # CGC substrate (starch, melibiose, …)
dnmb_plot_comparative_pazy(data_root)               # PAZy families

# Prophage heatmaps (purple palette)
dnmb_plot_comparative_phispy(data_root)     # PhiSpy regions bucketed by size
dnmb_plot_comparative_virsorter2(data_root) # VirSorter2 max_score_group
dnmb_plot_comparative_pide(data_root)       # PIDE regions bucketed by size
```

Outputs are written under `<data_root>/comparative/` as
`Comparative_<Module>_Heatmap.pdf` alongside the underlying count matrix.

Pass `auto_run_missing = FALSE` to skip the on-the-fly analysis and only
render what already exists.

**EggNOG-mapper** (Optional)

```python
emapper.py --cpu 20 --mp_start_method forkserver --data_dir [eggnog_data directory] -o out --output_dir [eggnog_output] --temp_dir [eggnog_output] --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i [fasta] --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel

```

- **Note:** [Strain of interest].emapper.annotations.xlsx or [Strain of interest]emapper.annotations.csv output are used for merging data.

**InterProScan** (Optional)

```python
./interproscan.sh -i [input_file] -f tsv -iprlookup -etra -goterms -pa -cpu 20
```

- **Note:** Files with [Strain of interest].tsv and [Strain of interest].tsv.sites extensions are used for merging data.


**Promotech** (Optional)

https://github.com/BioinformaticsLabAtMUN/PromoTech

Promotech can be appended as a DNMB module. It is disabled by default because
the upstream RF-HOT/RF-TETRA models are large and are not bundled with DNMB.
DNMB caches Promotech runtime files under the module cache
(`DNMB_CACHE_ROOT` or `~/.dnmb-cache`). If live prediction is requested without
a precomputed predictions file, the selected model is downloaded into that same
cache unless `promotech_download_model = FALSE`.
Each run writes `dnmb_module_promotech/promotech_promoter_feature_for_gb` for
copying into a GenBank FEATURES block and, when a GenBank input is available,
`dnmb_module_promotech/promotech_promoters_annotated.gbk` for SnapGene import.
Promoter feature labels include both the stable Promotech id and score, for
example `Promotech_000001 (score=0.62883)`.

```r
# Import a precomputed Promotech genome_predictions.csv/TSV file
run_DNMB(
  module_Promotech = TRUE,
  promotech_predictions = "genome_predictions.csv",
  promotech_threshold = 0.6
)

# Or run only the Promotech module against the active genbank_table
run_module_set(
  db = "Promotech",
  promotech_predictions = "genome_predictions.csv",
  merge = TRUE
)

# Live prediction: caches the Promotech repo and selected model first
run_DNMB(
  module_Promotech = TRUE,
  promotech_model = "RF-HOT",
  promotech_threshold = 0.6
)
```
- **Note:** Live Promotech execution requires the upstream Python dependencies
  and enough RAM. Precomputed `genome_predictions.csv` import works without
  downloading the model or running Promotech itself.
- **Runtime:** RF-HOT live prediction is intentionally heavy because upstream
  Promotech scans every 40-nt window on both strands. A 2.36 Mb bacterial
  genome took about 26 minutes and roughly 12 GB RAM in Docker during local
  validation; DNMB promoter-to-gene mapping and GenBank/SnapGene artifact
  generation then completed in about 15 seconds.


## Contributing
We welcome contributions from the community! If you have suggestions for improvements, additional scripts, or updates to the database, please see our contributing guidelines for more information on how to get involved.



## License
This project is released under MIT licence, which allows for both personal and commercial use, modification, and distribution of our work, provided that proper credit is given.

We hope our resources will prove invaluable to your research in systems biology. For any questions or feedback, please don't hesitate to reach out through our GitHub issues or contact section.

## Citation
If you use this piepline, please cite:
```
[DNMB] DNMB: Programmable domestication of thermophilic bacteria through removal of non-canonical defense systems.
			 Sung, J.Y., Lee, M.H., Park, J.S., Kim, H.B., Ganbat, D., Kim, D.G., Cho, H.W., Suh, M.K., Lee, J.S., Lee, S.J., Kim, S.B.*, and Lee, D.W.*.
			 *bioRxiv* 2026.03.21.173436. (2026)  
```
Please, cite also the underlying algorithm/database if it was used for the search step of DNMB:
```
  [EggNOG-mapper v2]    eggNOG-mapper v2: Functional annotation, orthology assignments, and domain prediction at
                        the metagenomic scale. Carlos P. Cantalapiedra, Ana Hernandez-Plaza,
                        Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021. Molecular Biology and Evolution
                        38(12):5825-5829. https://doi.org/10.1093/molbev/msab293

  [CLEAN]               Enzyme function prediction using contrastive learning. Tianhao Yu, Haiyang Cui, Jianan Canal Li,
                        Yunan Luo, Guangde Jiang, Huimin Zhao. 2023. Science 379(6639):1358-1363. 
                        https://doi.org/10.1126/science.adf2465

  [InterProScan]        InterProScan 5: genome-scale protein function classification.
                        Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla,
                        Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn,
                        Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, Sarah Hunter.
                        2014. Bioinformatics 30(9):1236-1240. https://doi.org/10.1093/bioinformatics/btu031

  [DefenseFinder]       DefenseFinder: Systematic and quantitative view of the antiviral arsenal of prokaryotes.
                        Florian Tesson, Alexandre Herve, Ernest Mordret, Marie Touchon, Camille d'Humieres, Jean Cury,
                        Aude Bernheim. 2022. Nature Communications 13:2561. https://doi.org/10.1038/s41467-022-30269-9

  [REBASE]              REBASE-a database for DNA restriction and modification: enzymes, genes and genomes.
                        Richard J. Roberts, Tamas Vincze, Janos Posfai, Dana Macelis. 2010. Nucleic Acids Research
                        38(Database issue):D234-D236. https://doi.org/10.1093/nar/gkp874

  [GapMindAA]           GapMind: Automated annotation of amino acid biosynthesis.
                        Morgan N. Price, Adam M. Deutschbauer, Adam P. Arkin. 2020. mSystems 5(3):e00291-20.
                        https://doi.org/10.1128/mSystems.00291-20

  [GapMindCarbon]       Filling gaps in bacterial catabolic pathways with computation and high-throughput genetics.
                        Morgan N. Price, Adam M. Deutschbauer, Adam P. Arkin. 2022. PLoS Genetics 18(4):e1010156.
                        https://doi.org/10.1371/journal.pgen.1010156

  [MEROPS]              The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a
                        comparison with peptidases in the PANTHER database. Neil D. Rawlings, Alan J. Barrett,
                        Paul D. Thomas, Xiaosong Huang, Alex Bateman, Robert D. Finn. 2018. Nucleic Acids Research
                        46(D1):D624-D632. https://doi.org/10.1093/nar/gkx1134

  [dbCAN]               dbCAN3: automated carbohydrate-active enzyme and substrate annotation.
                        Jinfang Zheng, Qiwei Ge, Yuchen Yan, Xinpeng Zhang, Le Huang, Yanbin Yin. 2023.
                        Nucleic Acids Research 51(W1):W115-W121. https://doi.org/10.1093/nar/gkad328

  [PAZy]                Plastics degradation by hydrolytic enzymes: The plastics-active enzymes database-PAZy.
                        Patrick C. F. Buchholz, Golo Feuerriegel, Hongli Zhang, Pablo Perez-Garcia,
                        Lena-Luisa Nover, Jennifer Chow, Wolfgang R. Streit, Jurgen Pleiss. 2022.
                        Proteins 90(7):1443-1456. https://doi.org/10.1002/prot.26325

  [ISelement]           ISEScan: automated identification of insertion sequence elements in prokaryotic genomes.
                        Zhiqun Xie, Haixu Tang. 2017. Bioinformatics 33(21):3340-3347.
                        https://doi.org/10.1093/bioinformatics/btx433

						ISfinder: the reference centre for bacterial insertion sequences. 
                        Philippe Siguier, Jerome Perochon, Lucie Lestrade, Jacques Mahillon,
                        Michael Chandler. 2006. Nucleic Acids Research 34(Database issue):D32-D36.
                        https://doi.org/10.1093/nar/gkj014

  [PhiSpy]              PhiSpy: a novel algorithm for finding prophages in bacterial genomes that combines similarity-
                        and composition-based strategies. Sajia Akhter, Rashedul Aziz,
                        Robert A. Edwards. 2012. Nucleic Acids Research 40(16):e126. https://doi.org/10.1093/nar/gks406

  [VirSorter2]          VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses.
                        Jiarong Guo, Ben Bolduc, Ahmed A Zayed, Arvind Varsani, Guillermo Dominguez-Huerta, Tom O Delmont,
                        Akbar Adjie Pratama, M Consuelo Gazitúa, Dean Vik, Matthew B Sullivan, Simon Roux. 2021. Microbiome 9:37.
                        https://doi.org/10.1186/s40168-020-00990-y

  [PIDE]                PIDE: a deep learning-based tool for prophage identification using genome-wide features.
                        https://github.com/BackofenLab/PIDE
```
