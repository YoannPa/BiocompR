% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_stat_tests.R
\name{prepare_annot_asso}
\alias{prepare_annot_asso}
\title{Prepares annotations to be tested for associations.}
\usage{
prepare_annot_asso(annot.table, verbose = FALSE)
}
\arguments{
\item{annot.table}{A \code{data.frame} containing all annotations, 1
annotation per column.}

\item{verbose}{A \code{logical} to display information about the
step-by-step processing of the data if TRUE
(Default: verbose = FALSE).}
}
\value{
A \code{list} containing updated annotations, the number of
        annotations available, and the annotation table itself.
}
\description{
Prepares annotations to be tested for associations.
}
\examples{
# Create mtcars annotation table
mtcars_annot <- data.table::as.data.table(
    mtcars[, c("cyl", "vs", "am", "gear", "carb")], keep.rownames = "cars")
# Prepare the annotations for testing
ls_annot <- prepare_annot_asso(annot.table = mtcars_annot, verbose = TRUE)
}
\author{
Yoann Pageaud.
}
