% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_stat_tests.R
\name{check_fun}
\alias{check_fun}
\title{Checks if a function exists and package of origin.}
\usage{
check_fun(fun, param.name = "fun", ncores = 1)
}
\arguments{
\item{fun}{A \code{character} string matching a function you are
looking for. The syntax must either be 'function' or
'package::function'.}

\item{param.name}{A \code{character} specifying the name of the parameter to
which the function is related. This name is used in error
messages.}

\item{ncores}{An \code{integer} to specify the number of cores/threads to
be used to parallel-compute distances for dendrograms.}
}
\value{
A \code{type} object returned description.
}
\description{
Checks if a function exists and package of origin.
}
\examples{
# Looking for the function "abbreviate"
check_fun(fun = "abbreviate") #Result should be: "base::abbreviate"
}
\author{
Yoann Pageaud.
}
