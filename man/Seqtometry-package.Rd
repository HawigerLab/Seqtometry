\name{Seqtometry-package}
\alias{Seqtometry-package}
\alias{Seqtometry}
\docType{package}
\title{
  Signature scoring for single cell analysis
}
\description{
  This package provides functions used in Seqtometry (Kousnetsov et al. 2024), a method for analyzing single cell (scRNA-seq or scATAC-seq) data via signature (gene set) enrichment scores. The Seqtometry scores may be useful for annotating or characterizing cells, either in a flow cytometry like workflow (where scores are standalone features used for progressive partitoning as demonstrated in the Seqtometry publication) or in a cluster-based workflow (as features of clusters). The exported impute function (a port of Python's MAGIC-impute, van Dijk et al. 2018), may also be useful for single cell analysis on its own.
}
\author{
Robert Kousnetsov <robert.kousnetsov@health.slu.edu>

Maintainer: Robert Kousnetsov <robert.kousnetsov@health.slu.edu>
}
