% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggasso.R
\name{ggasso.annot_pc}
\alias{ggasso.annot_pc}
\title{Plots association tests' results between some annotations and some PCs.}
\usage{
ggasso.annot_pc(
  annot.table,
  prcomp.res,
  perm.count = 10000,
  max.PCs = 8,
  p.treshold = 0.05,
  dataset.name = "dataset",
  PC.origin = NULL,
  verbose = FALSE
)
}
\arguments{
\item{annot.table}{A \code{data.frame} containing annotations, 1 annotation
per column.}

\item{prcomp.res}{A PCA result of classes \code{prcomp} or
\code{irlba_prcomp} resulting from stats::prcomp() or
irlba::prcomp_irlba().}

\item{perm.count}{An \code{integer} specifying the number of permutations to
realize on a vector, for the permutations matrix
initialization, to be used for calculating the significance
of a correlation test (Default: perm.count = 10000).}

\item{max.PCs}{An \code{integer} specifying the maximum number of
principal components to consider for association tests with
annotations (Default: max.PCs = 8).}

\item{p.treshold}{A \code{numeric} specifying a value for the P-value
significance threshold (Default: p.treshold = 0.05).}

\item{dataset.name}{A \code{character} string specifying the name you want to
give to your dataset, that will be included in the plot
title (Default: dataset.name = "dataset").}

\item{PC.origin}{A \code{character} string specifying what kind of data the
prcomp object has been computed from (e.g. 'mtcars',
'iris', 'mpg', ...) to be appended to the plot title.}

\item{verbose}{A \code{logical} to display information about the
step-by-step processing of the data if TRUE
(Default: verbose = FALSE).}
}
\value{
A \code{type} object returned description.
}
\description{
Plots association tests' results between some annotations and some PCs.
}
\examples{
# Create mtcars annotation table
mtcars_annot <- data.table::as.data.table(
    mtcars[, c("cyl", "vs", "am", "gear", "carb")], keep.rownames = "cars")
# Create mtcars matrix with scaled and transposed data
mat_mtcars <- scale(as.matrix(
    mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")]))
# Compute PCA on mtcars matrix
pcs_mtcars <- prcomp(x = mat_mtcars)
# Plot association tests' results between annotations and PCs
ggasso.annot_pc(
    annot.table = mtcars_annot, prcomp.res = pcs_mtcars,
    dataset.name = "mtcars", PC.origin = "mpg, disp, hp, drat, wt & qsec",
    verbose = TRUE)
}
\author{
Yoann Pageaud.
}
