\name{plot_DNMB_mobileome_distribution}
\alias{plot_DNMB_mobileome_distribution}
\title{Plot DNMB mobileome distribution}
\usage{
plot_DNMB_mobileome_distribution(
  output_dir,
  genbank = NULL,
  prefix = "mobileome_distribution",
  width = 15,
  height = 10,
  min_contig_fraction = 0.12
)
}
\arguments{
\item{output_dir}{DNMB mobileome result directory.}

\item{genbank}{Optional GenBank path. If \code{NULL}, the function tries to
infer it from \code{mobileome_summary.tsv}.}

\item{prefix}{Prefix used for generated plot filenames.}

\item{width}{Plot width in inches.}

\item{height}{Plot height in inches.}

\item{min_contig_fraction}{Minimum displayed contig width as a fraction of the
largest contig, used to keep short plasmids visible.}
}
\value{
An invisible list containing plot objects and saved file paths.
}
\description{
Create integrated ggplot summaries of merged IS distributions,
sequence-supported calls, targetable regions, and family-wise recognition and
count summaries.
}
