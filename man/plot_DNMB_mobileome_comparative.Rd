\name{plot_DNMB_mobileome_comparative}
\alias{plot_DNMB_mobileome_comparative}
\title{Plot comparative mobileome loci}
\usage{
plot_DNMB_mobileome_comparative(comparative, output_dir, top_n = 6)
}
\arguments{
\item{comparative}{A comparative result object returned by the DNMB
comparative layer.}

\item{output_dir}{Directory where plots should be written.}

\item{top_n}{Number of loci to plot.}
}
\value{
An invisible list with plot file paths.
}
\description{
Visualize focal and related genomes around comparative mobileome loci using
gggenes, highlighting insertion, deletion, and transposition-related locus
differences.
}
