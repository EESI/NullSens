\name{mvrStandard}
\alias{mvrStandard}
\title{OLS Linear Regression}
\description{Ordinary Least Squares Linear Regression}
\usage{mvrStandard(CDM, X, sites_sel)}
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
\references{Rao, C.R. (1973). Linear statistical inference and its applications (2nd ed.). New York: John Wiley & Sons.}
\author{Steve Essinger}
