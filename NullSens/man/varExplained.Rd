\name{varExplained}
\alias{varExplained}
\title{Partitioning Abiotic and Biotic Variation}
\description{Partitions community variation attributed to abiotic and biotic factors}
\usage{varExplained(CDM,p_value,sites_sel,CR, Avg_Adj_R2,Adj_R2,alpha)}
\arguments{
  \item{CDM}{Community Data Matrix (sites x species)}
  \item{p_value}{p-value from species interaction significance test}
  \item{sites_sel}{Nested list of booleans for including/excluding sites per species}
  \item{CR}{Pair-wise species correlation matrix}
  \item{Avg_Adj_R2}{Community Adjusted R2}
  \item{Adj_R2}{Adjusted R2, per species}
  \item{alpha}{Significance level for hypothesis test}
}
\value{
	\item{IND_summary}{Abiotic, Biotic and Unexplained Variation (per species)}
	\item{COM_summary}{Abiotic, Biotic and Unexplained Variation (community)}
}
\author{Steve Essinger}
