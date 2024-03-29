% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callECNV.R
\name{performFitPoisson}
\alias{performFitPoisson}
\title{Fit data}
\usage{
performFitPoisson(mergedCovFile, parameterFile, path = NULL, lowdepth = 10)
}
\arguments{
\item{mergedCovFile}{Path of merged coverage file.}

\item{parameterFile}{Path of metrics file.}

\item{path}{Path to write to.}

\item{lowdepth}{A numeric value, regions that avaerage depth less than this
value will replaced by NA.}
}
\description{
Fit data
}
\author{
Zhan-Ni Chen
}
