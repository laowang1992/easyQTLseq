% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterDP.R
\name{filterDP}
\alias{filterDP}
\alias{filterDP.BothParent}
\alias{filterDP.HighParent}
\alias{filterDP.LowParent}
\alias{filterDP.WithoutParent}
\title{Filter SNP according depth}
\usage{
filterDP(x, ...)

\method{filterDP}{BothParent}(
  x,
  minHPdp = 6,
  minLPdp = 6,
  minHBdp = 6,
  minLBdp = 6,
  maxHPdp,
  maxLPdp,
  maxHBdp,
  maxLBdp
)

\method{filterDP}{HighParent}(x, minHPdp = 6, minHBdp = 6, minLBdp = 6, maxHPdp, maxHBdp, maxLBdp)

\method{filterDP}{LowParent}(x, minLPdp = 6, minHBdp = 6, minLBdp = 6, maxLPdp, maxHBdp, maxLBdp)

\method{filterDP}{WithoutParent}(x, minHBdp = 6, minLBdp = 6, maxHBdp, maxLBdp)
}
\arguments{
\item{x}{A QTLseq S3 object.}

\item{minHPdp}{Minimum coverage depth of high phenotype parent.}

\item{minLPdp}{Minimum coverage depth of low phenotype parent.}

\item{minHBdp}{Minimum coverage depth of high phenotype bulk.}

\item{minLBdp}{Minimum coverage depth of low phenotype bulk.}

\item{maxHPdp}{Maximum coverage depth of high phenotype parent, default if average+3*sd.}

\item{maxLPdp}{Maximum coverage depth of low phenotype parent, default if average+3*sd.}

\item{maxHBdp}{Maximum coverage depth of high phenotype bulk, default if average+3*sd.}

\item{maxLBdp}{Maximum coverage depth of low phenotype bulk, default if average+3*sd.}
}
\value{
filtered QTLseq S3 object.
}
\description{
For low coverage depth SNP may have low reliability and accuracy, and extremely high coverage depth may be derived from repetitive sequence.
This function filter low coverage depth and extremely high coverage depth SNP.
}
\examples{
library(easyQTLseq)
# Example with sample data from a GATK table.
file_path <- system.file("extdata", "subset.table.gz", package = "easyQTLseq")
# readr::read_tsv() has a faster speed than read.table() when reading a file.
data <- readr::read_tsv(file = file_path)
x <- select_sample_and_SNP(data = data, highP = "qY", lowP = "R3", highB = "Y", lowB = "R", popType = "F2", bulkSize = c(30, 30))
x_filter <- filterDP(x = x)
}
