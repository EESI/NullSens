\name{sitesSelect}
\alias{sitesSelect}
\title{Site Selection Procedure}
\description{Procedure used to select sites for inclusion in analysis, per species.}
\usage{sitesSelect(CDM, X)}
\arguments{
  \item{CDM}{Community Data Matrix (sites x species)}
  \item{X}{Abiotic Factors (sites x factors) + vector of ones for intercept}
}
\value{
\item{sites_sel}{Nested list of booleans for including/excluding sites per species}
}
\author{Steve Essinger}
