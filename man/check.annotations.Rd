% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_ggplot_fun.R
\name{check.annotations}
\alias{check.annotations}
\title{Checks matching between annotation groups and annotation palettes.}
\usage{
check.annotations(data, annot.grps, annot.pal, verbose = FALSE)
}
\arguments{
\item{data}{A \code{matrix} or a \code{data.frame} with column names.}

\item{annot.grps}{A \code{list} of vectors of groups to which variables
belongs for the annotation sidebars. Vectors' lengths have
to match the number of variables.}

\item{annot.pal}{A \code{vector} or a list of vectors containing colors as
characters for the annotation sidebars. The length of
vectors has to match the number of levels of vectors listed
in 'annot.grps'. If a list is provided, its length must
match the length of the list provided to 'annot.grps'.}

\item{verbose}{A \code{logical} to specify wether the function should be
run on verbose mode (verbose = TRUE) or not
(Default: verbose = FALSE).}
}
\value{
An error message if something goes wrong during annotations checks.
}
\description{
Checks matching between annotation groups and annotation palettes.
}
\author{
Yoann Pageaud.
}
\keyword{internal}
