% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_ggplot_fun.R
\name{build.col_table}
\alias{build.col_table}
\title{Builds color tables for legends.}
\usage{
build.col_table(annot.grps, annot.pal, facet = NULL)
}
\arguments{
\item{annot.grps}{A \code{list} of characters or numerics specifying the
categories for legends, in their order of input.}

\item{annot.pal}{A \code{list} of character vectors, or a \code{vector} of
characters matching valid codes for R colors.}

\item{facet}{A \code{character} matching an annotation name in
'annot.grps' to be used to split heatmap in separate panels
following the annotation.}
}
\value{
A \code{data.table} list of categories matching their colors for each
        legend.
}
\description{
Builds color tables for legends.
}
\author{
Yoann Pageaud.
}
\keyword{internal}
