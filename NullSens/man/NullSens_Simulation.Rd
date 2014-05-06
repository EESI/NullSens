\name{NullSens_Simulation}
\alias{NullSens_Simulation}
\title{Evaluating NullSens on Simulated Community Data}
\description{This function is useful for performing sensitivity analyses on NullSens using simulated communities (see SpeciesResponseModel.R)}
\usage{
	NullSens_Simulation(n, p, q, N, M, CD, Type, num_exp = 100, select = TRUE, 
			reg_method = "robust", null_reps = 200, test_stat = c(1, 1), mutual_reject = 7, alpha = 0.05)
}
\arguments{
  \item{n}{Number of Sites to Include in Simulation}
  \item{p}{Number of Species to Include in Simulation}
  \item{q}{Number of Abiotic Factors to Include in Simulation}
  \item{N}{Noise Parameter Selection (integer 1 thru 12)}
  \item{M}{Covariation Magnitude Selection (integer 1 thru 5)}
  \item{CD}{Number of Covarying Species Pairs to Include in Simulation}
  \item{Type}{Type of Species Covariation in Simulation (see details)}
  \item{num_exp}{Number of Experiments to Include in Simulation}
  \item{select}{Enable site selection procedure -- TRUE (default)}
  \item{reg_method}{Regression method -- "robust" (defaut), "tobit", "standard"}
  \item{null_reps}{Number of random matrices to generate for null distribution -- 200 (default)}
  \item{test_stat}{Test statistic employed for computing indices}
  \item{mutual_reject}{Ignore species pairs with less than m mutual sites in testStatistic -- 8 (default)}
  \item{alpha}{Significance level -- 0.05 (default)}
}
\details{
Type:
	
0 -> All negatively covarying pairs

1 -> All positively covarying pairs

2 -> Mixed type covarying pairs (randomly selected)
}

\value{
  \item{Parameters}{List of parameters used in function call}
  \item{COM_summary}{Community variation paritioning, for each experiment}
  \item{p_values}{Detection p-value for each experiment}
  \item{avg_det_rate}{Number of communities showing significant iteraction over total experiments}
  \item{avg_det_type}{Number of communities with significant [positive,negative] interactions}
  \item{run_time}{Run-time of entire simulation (in seconds)}
}
\references{S.D. Essinger, "Partitioning Abiotic and Biotic Contributions to Community Variation," Ph.D. dissertation, Dept. Elect. Eng., Drexel Univ., Philadelphia, PA, 2013}
\author{Steve Essinger}