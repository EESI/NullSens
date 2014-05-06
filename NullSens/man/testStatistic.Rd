\name{testStatistic}
\alias{testStatistic}
\title{Test Statistic for Significant Covariation Test}
\description{Test Statistic employed to generate null distribution}
\usage{testStatistic(matrix, sites_sel, test_stat, mutual_reject)}
\arguments{
  \item{matrix}{Data matrix for computing test statistic}
  \item{sites_sel}{Nested list of booleans for including/excluding sites per species}
  \item{test_stat}{Selection of test statistic employed for computing indices for null distribution}
  \item{mutual_reject}{Threshold for excluding species pairs from test statistic}
}
\details{
test_stat (default c(1,1))

test_stat[1] = 1: absolute-covariance = 2: squared-covariance = 3: absolute-correlation = 4: squared-correlation

test_stat[2] = 1: absolute-covariance = 2: absolute-correlation = 3: sum of sites included per species pair = 4: 1 (constant)

mutual_reject (default 8): Any species pair having less than mutual_reject sites in common will be excluded from the test statistic calculation.
}
\value{
\item{index}{Covariation test statistic computed on input matrix}
\item{CR}{Pair-wise species correlation matrix}
\item{CV}{Pair-wise species covariation matrix}
}
\author{Steve Essinger}