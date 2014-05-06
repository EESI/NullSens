\name{nullModel}
\alias{nullModel}
\title{Null model}
\description{Randomize matrix using null model}
\usage{nullModel(matrix, sites_sel)}
\details{This null model randomly permutes species abundance using the respective sites selected in sites_sel. Each species is randomized independently of the other species.}
\arguments{
  \item{matrix}{Data matrix for randomization}
  \item{sites_sel}{Nested list of booleans for including/excluding sites per species}
}
\value{
\item{rand_matrix}{Return matrix randomized by null model}
}
\author{Steve Essinger}
