% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poc_functions.R
\name{poc1}
\alias{poc1}
\title{Function poc1}
\usage{
poc1(n, p = 0.5, structure = "any")
}
\arguments{
\item{n}{The number of nodes in the graph}

\item{p}{The probability of including each element in the subsets}

\item{structure}{The desired graph structure: 'any', 'tree', or 'chain'}
}
\value{
A list containing the perfect ordering of cliques and the separators
}
\description{
Draws (randomly, non-uniformly) a perfect ordering of cliques of n nodes (and thus a decomposable graph). Its output CLI is the input to the main procedure.
}
\examples{
poc1(10, structure = 'tree')
}
\keyword{graph,}
\keyword{ordering}
\keyword{perfect}
\keyword{random}
