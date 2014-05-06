\name{SpeciesResponseModel}
\alias{SpeciesResponseModel}
\title{Species Response Model - for use with simulations}
\description{This function generates a community data matrix (CDM) with n sites and p species to q environmental factors (X), using a linear model. The number (CD) of species pairs will covary within the community. The type of covariation for each pair depends on Type (0 negative, 1 positive, 2 mixed). The magnitude of covariation depends on Magnitude(M), while the noise dependson Noise(N). Both the noise and the covariation magnitude are added to the environmental species responses.
}
\usage{SpeciesResponseModel(n, p, q, N, M, CD, Type, Noise, Magnitude)}
\arguments{
  \item{n}{Number of Sites to Include in Simulation}
  \item{p}{Number of Species to Include in Simulation}
  \item{q}{Number of Abiotic Factors to Include in Simulation}
  \item{N}{Noise Parameter Selection (integer 1 thru 12)}
  \item{M}{Covariation Magnitude Selection (integer 1 thru 5)}
  \item{CD}{Number of Covarying Species Pairs to Include in Simulation}
  \item{Type}{Type of Species Covariation in Simulation (see details)}
  \item{Noise}{Vector of Noise Parameters}
  \item{Magnitude}{Vector of Covariation Magnitude Parameters}
}
\value{
  \item{CDM}{Community Data Matrix (Sites x Species)}
  \item{X}{Environmental Data (Explanatory Variables) (Includes 1's intercept)}
  \item{Y}{Environmental Species Responses}
  \item{YN}{Noise (Normal and Additive)}
  \item{Mag}{Magnitude of Species Covariation}
  \item{CDMB}{Community Data Matrix prior to Zero Censoring}
  \item{B}{Species Response Parameters}
  \item{cvid}{Type of species covariation for each pair (mixed dynamics only)}
}
\author{Steve Essinger}
