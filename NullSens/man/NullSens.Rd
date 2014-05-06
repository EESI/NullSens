\name{NullSens}
\alias{NullSens}
\title{Partitioning Abiotic and Biotic Contributions to Community Variation}
\description{
NullSens is a method that employs a nullmodel based on existing species covariation models to assess community-wide biotic covariation after removing individual species responses to environmental gradients. Analytical output from NullSens provides the percentages of the variation that can be explained due to abiotic, biotic and unexplained factors. Within the biotic variation results, the algorithm also outputs if it detects positive, negative, or mixed community covariation. NullSens can be used to verify hypotheses on species interactions when environmental gradients are present.}
\usage{NullSens(CDM, X, select = TRUE, reg_method = "robust", null_reps = 200,
			test_stat = c(1, 1), mutual_reject = 8, alpha = 0.05)}
\arguments{
  \item{CDM}{Community Data Matrix (sites x species) *Must be a matrix}
  \item{X}{Abiotic Factors (sites x factors) *Must be a matrix (i.e. use as.matrix())}
  \item{select}{Enable site selection procedure -- TRUE (default)}
  \item{reg_method}{Regression method -- "robust" (defaut), "tobit", "standard"}
  \item{null_reps}{Number of random matrices to generate for null distribution -- 200 (default)}
  \item{test_stat}{Test statistic employed for computing indices}
  \item{mutual_reject}{Ignore species pairs with less than m mutual sites in testStatistic -- 8 (default)}
  \item{alpha}{Significance level -- 0.05 (default)}
}
\details{
Test Statistic Options: test_stat[1] -> Statistic, test_stat[2] -> Weight

test_stat[1] = 1 -> abs(cvar), = 2 -> cvar^2, = 3 -> abs(ccorr), = 4 -> ccorr^2

test_stat[2] = 1 -> abs(cvar), = 2 -> abs(ccorr), = 3 -> sum(mutsel), = 4 -> 1}

\value{
\item{CDM}{Community Data Matrix}
\item{X}{Abiotic Factors, vector of ones appended if missing}
\item{Yhat}{Predicted (fitted) responses}
\item{Yres}{Residual Responses}
\item{B_est}{Estimated Regression Parameters}
\item{sites_sel}{List Sites Selected for Analysis, per Species}
\item{p_value}{P-Value of Covariation Significance Test}
\item{index}{List of index values for random and test matrices (test is last element)}
\item{CR}{Pairwise Residual (Yres) Correlation Matrix from testStatistic}
\item{CV}{Pairwise Residual (Yres) Covaration Matrix from testStatistic}
\item{R2}{Coefficient of Multiple Determination, per species}
\item{Adj_R2}{Adjusted R2, per species}
\item{Avg_R2}{Community Averaged R2}
\item{Avg_Adj_R2}{Community Averaged Adjusted R2}
\item{IND_summary}{abiotic, biotic, unexplained variation, per species}
\item{COM_summary}{abiotic, biotic, unexplained variation, community}
\item{COM_var_type}{[avg_covar, p-value_pos, p_value_neg]}
}
\references{S.D. Essinger, "Partitioning Abiotic and Biotic Contributions to Community Variation," Ph.D. dissertation, Dept. Elect. Eng., Drexel Univ., Philadelphia, PA, 2013}
\author{Steve Essinger}