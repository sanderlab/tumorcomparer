% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_weighted_corr.R
\name{calc_weighted_corr}
\alias{calc_weighted_corr}
\title{Calculate weighted correlation}
\usage{
calc_weighted_corr(a, b, w = rep(1, nrow(a))/nrow(a))
}
\arguments{
\item{a}{numeric matrix}

\item{b}{numeric matrix}

\item{w}{alteration weights}
}
\value{
a weighted correlation (similarity) matrix
}
\description{
Calculate weighted correlation
}
\author{
Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
}
\concept{
tumorcomparer
}
