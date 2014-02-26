NullSensWrap <- function(CDM, X, select = TRUE, reg_method="robust", null_reps=200, test_stat_type = 1, mutual_reject=7, alpha=0.05) {

######################################################################################
# 02/19/14
# Steven D. Essinger

# INPUT DATA:
# CDM -- Community Data Matrix (sites x species) *Must be a matrix 
# X -- Abiotic Factors (sites x factors) *Must be a matrix (i.e. use as.matrix())

# INPUT PARAMETERS:
# select = TRUE (default) -- enable site selection procedure
# reg_method = "robust" (defaut), "tobit", "standard" -- regression method
# null_reps = 1000 (default) -- number of random matrices to generate for null distribution
# test_stat_type = 1 (default) -- test statistic employed for computing indices
# mutual_reject = 7 (default) -- ignore pairs with less than 7 mutual sites in testStatistic
# alpha = 0.05 (default) -- significance level

# OUTPUT DATA:
# CDM -- Community Data Matrix
# X -- Abiotic Factors, vector of ones appended if missing
# Yhat -- Predicted (fitted) responses
# Yres -- Residual Responses
# B_est -- Estimated Regression Parameters
# sites_sel -- List Sites Selected for Analysis, per Species
# p_value -- P-Value of Covariation Significance Test
# index -- List of index values for random and test matrices (test is last element)
# CR -- Pairwise Residual (Yres) Correlation Matrix  from testStatistic
# CV -- Pairwise Residual (Yres) Covaration Matrix from testStatistic
# R2 -- Coefficient of Multiple Determination, per species
# Adj_R2 -- Adjusted R2, per species
# Avg_R2 -- Community Averaged R2
# Avg_Adj_R2 -- Community Averaged Adjusted R2
# summary -- abiotic, biotic, unexplained variation, per species
# AVGsummary -- summary averaged over all species

######################################################################################
# Source the support functions
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/sitesSelect.R")
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/mvrRobust.R")
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/coeffDet.R")
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/nullModel.R")
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/testStatistic.R")

#####################################################################################
# Check input data
# Insert column of ones in X if missing
if (identical(X[,1],matrix(1,nrow(X),1))==FALSE){
	print('X missing column of ones. Automatically included.')
	X <- cbind(matrix(1,nrow(X),1),X)
}

# Check CDM for negative abundances
if (length(which(CDM < 0)) > 0){
	print('CDM contains negative abundance(s). Aborted.')
	return()
}
#####################################################################################

n <- nrow(CDM) # Number of Sites
p <- ncol(CDM) # Number of Species
q <- ncol(X) # Number of Abiotic Factors

# CHOOSE SITES TO INCLUDE IN ANALYSIS
if (select) sites_sel <- sitesSelect(CDM,X,n,p,q) # Select sites
else sites_sel <- list(rep(list(rep(TRUE,n)),p)) # Include all sites

# CHOOSE MULTIVARIATE REGRESSION METHOD
if (reg_method == "standard") mvr_out <- mvrStandard(CDM,X,sites_sel,n,p,q)
else if (reg_method == "robust") mvr_out <- mvrRobust(CDM,X,sites_sel,n,p,q)
else if (reg_method == "tobit") mvr_out <- mvrTobit(CDM,X,sites_sel,n,p,q)
else print ('Invalid regression method.')

# ENVIORNMENTAL VARIATION EXPLAINED - COFF. DET.
coeff_out <- coeffDet(CDM,X,mvr_out$Yhat)

# GENERATE NULL DISTRIBUTION
index <- c(rep(0,null_reps))
for (i in 1:null_reps-1){
	rand_matrix <- nullModel(mvr_out$Yres,sites_sel)
	tS_out <- testStatistic(rand_matrix,sites_sel,n,p,q,test_stat_type,mutual_reject)
	index[i] = tS_out$index
}
# Index Computed on Residuals of CDM under Test
tS_out = testStatistic(mvr_out$Yres,sites_sel,n,p,q,test_stat_type,mutual_reject)
index[i+1] = tS_out$index
CR = tS_out$CR # Correlation matrix for test residuals
CV = tS_out$CV # Covariation matrix for test residuals

# SIGNIFICANCE TEST
# Determine if residuals of CDM covary significantly via p-value
count_for_p <- which(index >= index[null_reps]) # Number of random matrices with index greater than the test matrix
p_value <- length(count_for_p)/null_reps # p-value of covaration significance test

# TYPE OF COMMUNITY COVARIATION
# Determine if community exhibits significant positive or negative covariation

# VARIATION PARTITIONING


#####################################################################################

result = list('CDM'=CDM,'X'=X,'Yhat'=mvr_out$Yhat,'Yres'=mvr_out$Yres,'B_est'=mvr_out$B_est,'sites_sel'=sites_sel,'p_value'=p_value,'test_indices'=index,'CR'=CR,'CR'=CV,'R2'=coeff_out$R2,'Adj_R2'=coeff_out$Adj_R2, 'Avg_R2'=coeff_out$Avg_R2, 'Avg_Adj_R2'=coeff_out$Avg_Adj_R2)
}