% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggbipca.R
\name{prepare_pca_data}
\alias{prepare_pca_data}
\title{Collects and computes needed metrics for PCA biplot.}
\usage{
prepare_pca_data(prcomp.res, dt.annot, PCs, scale)
}
\arguments{
\item{prcomp.res}{A PCA result of classes \code{prcomp} or
\code{irlba_prcomp} resulting from stats::prcomp() or
irlba::prcomp_irlba().}

\item{dt.annot}{A \code{data.table} containing annotations matching data in
'prcomp.res'.}

\item{PCs}{An \code{integer} vector matching principal components to
be used to generate the PCA biplot.}

\item{scale}{A \code{double} scaling parameter, disabled by 0.}
}
\value{
A \code{list} containing PC IDs, scaled values of the PCA, and PCs
        labels to be displayed on the margins of each axis.
}
\description{
Collects and computes needed metrics for PCA biplot.
}
\author{
Yoann Pageaud.
}
