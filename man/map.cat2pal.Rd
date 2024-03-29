% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_ggplot_fun.R
\name{map.cat2pal}
\alias{map.cat2pal}
\title{Maps annotations categories to palettes.}
\usage{
map.cat2pal(origin.grps, groups, annot.pal)
}
\arguments{
\item{origin.grps}{A \code{list} of characters or numerics specifying the
categories for legends, in their order of input.}

\item{groups}{A \code{list} of characters or numerics specifying the
categories for legends, in their order of appearance on
the annotation sidebar of the plot.}

\item{annot.pal}{A \code{list} of character vectors matching valid codes
for R colors.}
}
\value{
A \code{list} of data.tables where categories of each annotation are
        matched with their color palette.
}
\description{
Maps annotations categories to palettes.
}
\author{
Yoann Pageaud.
}
\keyword{internal}
