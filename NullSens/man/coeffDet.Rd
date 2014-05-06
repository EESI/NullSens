\name{coeffDet}
\alias{coeffDet}
\title{Compute the coefficient of determination}
\description{NullSens internal function: Computes the coefficient of determination for each species and the overall community.}
\usage{coeffDet(CDM, X, Yhat)}
\arguments{
  \item{CDM}{Community Data Matrix (sites x species)}
  \item{X}{Abiotic Factors (sites x # of abiotic factors + 1 (intercept))}
  \item{Yhat}{Predicted (fitted) responses}
}
\details{Avg_Adj_R2 is the environmentally explained variation}
\value{
\item{R2}{Coefficient of Multiple Determination, per species}
\item{Adj_R2}{Adjusted R2, per species}
\item{Avg_R2}{Community R2}
\item{Avg_Adj_R2}{Community Adjusted R2}
}
\author{Steve Essinger}