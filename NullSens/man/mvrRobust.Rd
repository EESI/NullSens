\name{mvrRobust}
\alias{mvrRobust}
\title{Robust Regression using bisquare estimator objective function}
\description{Multivariate regression using iteratively reweighted least-squares bisquare estimator}
\note{Requires MASS package} 
\usage{mvrRobust(CDM, X, sites_sel)}
\arguments{
  \item{CDM}{Community Data Matrix (sites x species)}
  \item{X}{Abiotic Factors (sites x factors) + vector of ones for intercept}
  \item{sites_sel}{Nested list of booleans for including/excluding sites per species}
}
\value{
	\item{Yhat}{Predicted (fitted) responses}
 	\item{Yres}{Residual responses}
 	\item{B_est}{Model coefficient estimates}
}
\references{Robust Statistics, Peter. J. Huber, Wiley, 1981 (republished in paperback, 2004)}
\author{Steve Essinger}
