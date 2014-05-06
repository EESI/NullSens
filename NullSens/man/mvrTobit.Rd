\name{mvrTobit}
\alias{mvrTobit}
\title{Censored Regression (Tobit) Model}
\description{Fitting a linear model with a censored dependent variable}
\usage{mvrTobit(CDM, X, sites_sel)}
\note{Requires censReg package}
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
\references{
Greene, W.H. (2008): Econometric Analysis, Sixth Edition, Prentice Hall, p. 871-875.

Kleiber, C. and Zeileis, A. (2008): Applied Econometrics with R, Springer, p. 141-143.

Tobin, J. (1958): Estimation of Relationships for Limited Dependent Variables. Econometrica 26,
p. 24-36.
}
\author{Steve Essinger}

