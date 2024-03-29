% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poc_functions.R
\name{get.MH.poc}
\alias{get.MH.poc}
\title{get.MH.poc}
\usage{
get.MH.poc(CLIQUES, mh_burn = 10)
}
\arguments{
\item{CLIQUES}{A list of cliques.}

\item{mh_burn}{Integer indicating the number of Metropolis-Hastings burn-in iterations to be used.}
}
\value{
Returns the generated POC.
}
\description{
The function get.MH.poc generates a perfect ordering of cliques (POC) using the Metropolis-Hastings method.
}
\examples{
get.MH.poc(CLIQUES,1000)
}
