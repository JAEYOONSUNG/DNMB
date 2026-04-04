\name{run_DNMB_mobileome}
\alias{run_DNMB_mobileome}
\title{Run the DNMB mobileome impact pipeline}
\usage{
run_DNMB_mobileome(
  genbank,
  output_dir = NULL,
  related_genbanks = NULL,
  related_detection_mode = c("annotation", "hybrid", "sequence"),
  detection_mode = c("hybrid", "sequence", "annotation"),
  prodigal_mode = c("single", "meta"),
  isescan_dir = NULL,
  isfinder_db = NULL,
  isfinder_fasta = NULL,
  flank_bp = 10000,
  pair_max_dist_bp = 200000,
  export_excel = TRUE,
  save_rds = TRUE,
  verbose = TRUE,
  max_target_sites_per_family = 2000,
  max_variant_targets_per_element = 50,
  site_scan_strategy = c("balanced", "full")
)
}
\arguments{
\item{genbank}{Path to a single \code{.gb}, \code{.gbk}, or \code{.gbff} file.}

\item{output_dir}{Directory where result tables will be written. If \code{NULL},
a directory named \code{<genbank_basename>_dnmb_mobileome} is created beside the
input file.}

\item{related_genbanks}{Optional vector of related-strain GenBank paths used
for comparative hotspot and chronology inference.}

\item{related_detection_mode}{Detection mode used for related genomes in the
comparative layer. Defaults to \code{annotation} for speed.}

\item{detection_mode}{Detection mode. \code{hybrid} combines annotation and
sequence-driven calls, \code{sequence} uses only the sequence engine, and
\code{annotation} uses only annotation-derived calls.}

\item{prodigal_mode}{Prodigal mode used by the sequence engine, either
\code{single} or \code{meta}.}

\item{isescan_dir}{Optional path to an ISEScan installation directory. If
\code{NULL}, the directory is inferred from \code{isescan.py} in \code{PATH}.}

\item{isfinder_db}{Optional BLAST database prefix built from curated ISfinder
references or a compatible custom IS reference set.}

\item{isfinder_fasta}{Optional nucleotide FASTA used to create a BLAST
database for reference-based family refinement.}

\item{flank_bp}{Number of bases on each side of an element used when calling
nearby mutation-prone genes.}

\item{pair_max_dist_bp}{Maximum distance between two candidate elements on
the same contig when generating pair-mediated structural-variant scenarios.}

\item{export_excel}{Logical; if \code{TRUE}, also write an Excel workbook.}

\item{save_rds}{Logical; if \code{TRUE}, save the full result object as an RDS.}

\item{verbose}{Logical; if \code{TRUE}, print progress messages.}

\item{max_target_sites_per_family}{Maximum number of predicted target sites
retained per IS family after recognition-sequence and viability scoring.}

\item{max_variant_targets_per_element}{Maximum number of target-site variants
retained per source element when building the variant catalog.}

\item{site_scan_strategy}{Target-site scan strategy. \code{balanced} keeps
genome-native, higher-confidence models for speed; \code{full} scans all
available models.}
}
\value{
An invisible list containing parsed features, genes, IS/mobile elements,
mutation points, structural-variant scenarios, phenotype predictions, and
output file paths.
}
\description{
This pipeline accepts a single GenBank file and performs a hybrid mobile
element analysis focused on IS-like elements and their possible local genome
consequences. The workflow parses annotated features, optionally adds
ISEScan-inspired sequence evidence through ORF calling plus pHMM/singleton
searches and family-aware TIR/TSD heuristics, merges annotation and sequence
evidence, scores essentiality and redundancy from genome-native annotations,
builds genome-native target-site models from observed TSD/context information,
predicts viable target sites, infers
deletion/inversion/transposition scenarios that could be mediated by element
pairs, assembles a genome-variant catalog, and summarizes plausible phenotype
changes from affected genes. When related GenBank genomes are supplied, the
pipeline also infers comparative loci, occupied/empty hotspot matrices, and
proposed single-step versus multi-step event chronology.
}
