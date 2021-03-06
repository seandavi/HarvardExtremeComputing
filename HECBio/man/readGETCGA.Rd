% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/TCGA.R
\name{readGETCGA}
\alias{readGETCGA}
\title{Read Level 3 TCGA gene expression data}
\usage{
readGETCGA(uri)
}
\arguments{
\item{uri}{a url or filename of a TCGA data matrix}
}
\value{
A \code{SummarizedExperiment} object with
values in the \code{assay(object, 'exprs')},
sample names as colnames, and identifiers as
rownames.
}
\description{
This function reads a file of Level 3 gene expression
data and returns a basic \code{SummarizedExperiment}
object with the data in \code{assay(object, 'exprs')}.
}

