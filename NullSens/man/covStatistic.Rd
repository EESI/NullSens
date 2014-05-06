\name{covStatistic}
\alias{covStatistic}
\title{Test for significant positive/negative community covariation.}
\description{Covarition test statistic (species pair-wise) for significant community positive/negative covariance hypothesis testing}
\usage{covStatistic(matrix, sites_sel, mutual_reject)}
\arguments{
  \item{matrix}{Data matrix for computing test statistic}
  \item{sites_sel}{Nested list of booleans for including/excluding sites per species}
  \item{mutual_reject}{Threshold for excluding species pairs from test statistic}
}
\details{
mutual_reject (default 8): Any species pair having less than mutual_reject sites in common will be excluded from the test statistic calculation.
}
\value{
\item{index}{Covariation type test statistic computed on input matrix}
}
\author{Steve Essinger}