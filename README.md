# Welcome to the DNMB, a combined scripts-based pipeline for comprehensive genomics analysis of non-model bacteria.
![Fig  1](https://github.com/user-attachments/assets/ef35c8cc-5e75-40ea-a13e-42022ccbc0ca)

In the realm of microbiology, the domestication of non-model bacteria (DNMB) presents a unique set of challenges and opportunities. To address these, we introduce a novel approach—the DNMB pipeline—a comprehensive genomics analysis tool tailored specifically for non-model bacterial species. Unlike traditional model organisms, non-model bacteria often lack well-defined genetic resources and established analytical pipelines. Consequently, researchers face hurdles in elucidating the genetic basis underlying phenotypic traits relevant to domestication efforts. Our pipeline integrates multiple scripts and computational tools to streamline various stages of genomic analysis, from making comprehensive functional annotation in user-friendly table format to genomic features including codon usage and ribosomal binding site preference and distance. 
Herein, we provide a detailed overview of the DNMB pipeline, highlighting its key components and functionalities. Additionally, we demonstrate its utility through a case study involving the domestication of a non-model bacterial strain. The DNMB pipeline not only accelerates genomic analysis but also enhances our understanding of non-model bacterial physiology, thereby facilitating the exploitation of microbial diversity for biotechnological applications.


## Project Introduction

Welcome to our GitHub repository, where we're excited to share a series of workflows designed to streamline processes in systems biology. This repository is composed of various scripts, each tailored to specific tasks within our broader research framework. Additionally, we're providing access to a curated database to enhance your research capabilities.
- Principal investigator: Dong-Woo Lee
- Project lead: Jae-Yoon Sung
- Maintainers: Jae-Yoon Sung
- Contributors: Jae-Yoon Sung, Seong Do Kim

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
---
**DefenseFinder** (Optional)

![DefenseFinder_Heatmap](https://github.com/user-attachments/assets/68e953c1-d568-4b56-89f2-1107a60a6f6c)
```r
DefenseFinder_Heatmap()
```
- **Note:** [Strain of interest].defense_finder_systems.tsv output are used for merging data. GenBank’s SOURCE field is used for extracting names.

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
```python
python promotech.py -pg -m RF-HOT -f examples/genome/[my_fasta].fna -g -o results 

python promotech.py -g -m RF-HOT -t 0.6 -i results -o results
```
- **Note:** fasta must have only capital letters


## Contributing
We welcome contributions from the community! If you have suggestions for improvements, additional scripts, or updates to the database, please see our contributing guidelines for more information on how to get involved.



## License
This project is released under MIT licence, which allows for both personal and commercial use, modification, and distribution of our work, provided that proper credit is given.

We hope our resources will prove invaluable to your research in systems biology. For any questions or feedback, please don't hesitate to reach out through our GitHub issues or contact section.

## Citation
If you use this piepline, please cite:
```
[DNMB] DNMB: A Strategic Blueprint for the Domestication of Geobacillus stearothermophilus as a Thermophilic Platform using the DNMB Suite.
             Jae-Yoon Sung, Mun Hoe Lee, Hyungbin Kim, Dariimaa Ganbat, Hyun-Woo Cho, Sang Jae Lee, Seong Bo Kim, and Dong-Woo Lee. 2024.
             XXX, XXX, https://doi.org/XXX
```
Please, cite also the underlying algorithm if it was used for the search step of DNMB:
```
[eggNOG-mapper v2] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
                   prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
                   Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
                   Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293

[InterProScan] InterProScan 5: genome-scale protein function classification Philip Jones, David Binns,
			       Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam,
			       John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn,
			       Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez,
			       Sarah Hunter Bioinformatics (2014), PMID: 24451626

[Promotech] Promotech: A general tool for bacterial promoter recognition. Ruben Chevez-Guardado and Lourdes Peña-Castillo.
		       Genome Biol 22(1):318 (2021). PMID: 34789306. (https://doi.org/10.1186/s13059-021-02514-9)

```
Please, cite also the underlying Database if it was used for the search step of DNMB:
```
[REBASE]
[MEROPS]
[CAZy]
[LED]
```

