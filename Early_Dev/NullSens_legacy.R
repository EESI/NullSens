NullSens <- function(CDM, X, reps=200,  index_rej=7){

######################################################################################
# Version 1.0
# 10/10/13
# Steve Essinger

# INPUT DATA:
# CDM -- Community Data Matrix (sites x species) *Must be a matrix 
# X -- Abiotic Factors (sites x factors) *Must be a matrix (i.e. use as.matrix())

# INPUT PARAMETERS:
# reps = 200 (default) -- Number of random matrices to generate for null distribution
# index_rej = 7 (default) -- ignore covarying pairs with less than 7 mutual sites

# OUTPUT RETURNED:
# CDM -- Community Data Matrix
# X -- Abiotic Factors, vector of ones appended if missing
# Yhat -- Predicted (fitted) responses
# Yres -- Residual Responses
# BP -- Estimated Regression Parameters
# SitesSel -- Sites Selected for Analysis, per Species
# NN -- Number of Sites Selected per Species
# PRSqYres -- P-Value of Covariation Significance Test
# AvgRSqYres -- Index value of random and test matrices
# CR -- Pairwise Residual (Yres) Correlation Matrix (using SitesSel)
# CV -- Pairwise Residual (Yres) Covaration Matrix (using SitesSel)
# R2 -- Coefficient of Multiple Determination, per species
# Avg_R2 -- Community Averaged R2
# Adj_R2 -- Adjusted R2, per species
# Avg_Adj_R2 -- Community Averaged Adjusted R2
# R2_RDA -- Coefficient of Multiple Determination, RDA Method
# Adj_R2_RDA -- Adjusted Coefficient of Multiple Determination, RDA Method
# summary -- abiotic, biotic, unexplained variation, per species
# AVGsummary -- summary averaged over all species
######################################################################################

# Insert column of ones in X if missing
if (identical(X[,1],matrix(1,nrow(X),1))==FALSE){
	print('X missing column of ones. Automatically included.')
	X <- cbind(matrix(1,nrow(X),1),X)
}

if (length(which(CDM < 0)) > 0){
	print('CDM contains negative abundances. Aborted.')
	return()
}

n <- nrow(CDM) # Number of Sites
p <- ncol(CDM) # Number of Species
q <- ncol(X) # Number of Abiotic Factors

SitesSel <- vector('list',p) # Stores Sites Selected for each Species
BP <- matrix(0,q,p)  # Regression Parameters
Yhat <- matrix(0,n,p) # Predicted (Fitted) Species Responses
Yres <- matrix(0,n,p) # Residual Species Responses
R2 <- c(rep(0,p)) # Coefficient of multiple determination
NN <- c(rep(0,p)) # Number of Sites included in R2 calculation (Number of SitesSel)

# Site Selection and SitesSelion
for (i in 1:p){ # For each species
	posAb <- which(CDM[,i] > 0) # Choose sites with positive abundance for theshold
	remove <- NA # initialize sites for removal
	centroid <- c(rep(0,q)) 
	thresh <- c(rep(0,q))
	
	for (grad in 2:q){ # For each abiotic factor
		centroid[grad] <- quantile(X[posAb,grad],0.5) # Centroid of positive abundance
		distance_pos <- as.matrix(dist(c(centroid[grad],X[posAb,grad]))) # distance of positive abundances to centroid
		thresh[grad] <- mean(distance_pos[-1,1]) + sd(distance_pos[-1,1]) # threshold for site selection
	
		distance_all <- as.matrix(dist(c(centroid[grad],X[,grad]))) # Distance of all sites (positive and zero) to centroid
		remove_temp <- which(distance_all[-1,1] > thresh[grad]) # remove sites beyond threshold from further analysis
		remove <- c(remove,remove_temp) # keep a list of all sites removed for each species
	}
	
	keep = setdiff(1:n,remove) # Keep sites not marked for removal
	SitesSel[[i]] <- c(rep(FALSE,n))
	if (length(keep) > 3){ # remove species with < 4 selected sites from analysis
		SitesSel[[i]][keep] = TRUE	
	}
	
	# Regression Parameter Estimation -- Normal Equation
	tryCatch({
		BP[,i] <- solve(t(X[SitesSel[[i]],])%*%X[SitesSel[[i]],])%*%t(X[SitesSel[[i]],])%*%CDM[SitesSel[[i]],i]	
	}, error = function(ex) {
		paste('Inverse computation -- singular. Cannot compute regression parameters for Species', i) # solve fails due to singular (X'*X) given SitesSel[[i]]
	})
	
	rmnan <- which(BP[,i]==NA)
	BP[rmnan,i] <- 0
	
	Yhat[,i] <- X %*% BP[,i] # Compute fitted values
	rmneg <- which(Yhat[,i]<=0) # Remove negative fitted abundances
	Yhat[rmneg,i] = 0
	Yres[,i] <- CDM[,i]-Yhat[,i] # Compute Residuals

	# Compute the squared correlation for R2 Value
	if (sum(SitesSel[[i]]) > q+1){
		R2[i] <- cor(Yhat[SitesSel[[i]],i],CDM[SitesSel[[i]],i])^2
		NN[i] <- sum(SitesSel[[i]]) # Number of sites included in R2 calculation
	} else {
		R2[i] <- 0
		NN[i] <- 0
	}
}

# Weighted Community Averaged R2
Avg_R2 <- sum(R2*colSums(CDM))/sum(CDM)

# Weighted Community Averaged Adjusted R2
Adj_R2 <- 1-((1-R2)*(NN-1)/(NN-q-1))
setzero <- which(R2 == 0)
Adj_R2[setzero] <- 0
Avg_Adj_R2 <- sum(Adj_R2*colSums(CDM))/sum(CDM)

# Compute Adjusted RDA R2 -- using approach from RDA in vegan package
CDMM <- (diag(n)-1/n*matrix(1,n,n))%*%CDM # Center CDM
M <- colMeans(X)
XM <- matrix(0,n,q-1)
for (i in 1:(q-1)){
	XM[,i] <- X[,i+1]-M[i+1] # Remove means from abiotic factors
}
Yhat_RDA <- XM %*% solve(t(XM) %*% XM) %*% t(XM) %*% CDMM
R2_RDA <- sum(diag(t(Yhat_RDA) %*% Yhat_RDA))/sum(diag(t(CDMM) %*% CDMM))
Adj_R2_RDA <- 1-(1-R2_RDA)*(n-1)/(n-q-1); 

# Nullmodel and Index Computation
# reps: Number of random matrices generated in nullmodel
AvgRSqYres <- c(rep(0,reps))
for (i in 1:(reps-1)){
	# Nullmodel - generates random matrices
	RANDYres <- Yres # For each randomization always start with original Yres
	for (species in 1:p){
		permNon <- which(SitesSel[[species]]==TRUE)
		neworderNon <- sample(permNon, length(permNon), replace = FALSE, prob = NULL)
		RANDYres[permNon,species] <- Yres[neworderNon,species]
	}
	# Index calculation
	w <- 0
	AvgRSq <- 0
	for (species1 in 1:(p-1)){
		for (species2 in (species1+1):p){
			  temp <- SitesSel[[species1]] & SitesSel[[species2]] # Choose mutually selected sites of both species in pair
			  if (sum(temp) <= n*0.05 | sum(temp) <= index_rej){
			  	cvar <- 0
			  }	else {
			  	COV = cov(RANDYres[temp,c(species1,species2)]) # Covariance
			  	if (is.na(COV[1,2]) == TRUE | COV[1,2] == 0){
	           	cvar <- 0
	           } else {
	            	cvar <- COV[1,2];
	           } 
	         }
	         AvgRSq <- AvgRSq + cvar^2; # Index Calculation
	         w <- w + abs(cvar); # Weights for averaging the index
		}
	}
	if (w == 0){
		AvgRSqYres[i] <- 0
	} else {
		AvgRSqYres[i] <- AvgRSq/w  # Index of random matrix i of reps
	}
}

# Computation of index on CDM under test
w <- 0
AvgRSq <- 0
CR <- matrix(0,p-1,p)
CV <- matrix(0,p-1,p)
for (species1 in 1:(p-1)){
	for (species2 in (species1+1):p){
		  temp <- SitesSel[[species1]] & SitesSel[[species2]]
		  if (sum(temp) <= n*0.05 | sum(temp) <= index_rej){
		  	cvar <- 0
		  }	else {
		  	options(warn=-1) # Turn warnings off: OK for cor to return NA
		  	COV = cov(Yres[temp,c(species1,species2)]) # Covariation
		  	CORR = cor(Yres[temp,c(species1,species2)]) # Correlation
		  	if (is.na(COV[1,2]) == TRUE | COV[1,2] == 0){
           		cvar <- 0
           		CV[species1,species2] <- cvar
           		CR[species1,species2] <- 0
           } else {
            	cvar <- COV[1,2];
            	CV[species1,species2] <- cvar # Keep track of pairwise covariance
           		CR[species1,species2] <- CORR[1,2] # Tracking of pairwise correlation
           } 
           options(warn=0) # Warnings back on
         }
         AvgRSq <- AvgRSq + cvar^2; # Index Calculation
         w <- w + abs(cvar); # Weights for averaging the index
	}
}
if (w == 0){
	AvgRSqYres[reps] <- 0
} else {
	AvgRSqYres[reps] <- AvgRSq/w  # Index of test matrix
}

# Get P-Values to determine if residuals of CDM covary significantly
indRYres <- which(AvgRSqYres >= AvgRSqYres[reps]) # Number of random matrices with index greater than the test matrix
PRSqYres <- length(indRYres)/reps # P-value of covaration significance test

# Return Abiotic, Biotic, Unexplained Variances
absum <- c(rep(0,p)) # Compute raw abundance totals for each species
for (i in 1:p){
	absum[i] <- sum(CDM[SitesSel[[i]],i])
}
if (PRSqYres > 0.05){ # If no significant covariation found
	summary <- cbind(Adj_R2,matrix(0,p,1),1-Adj_R2) # Abiotic, Biotic, Unexplained
	AVGsummary <- cbind(sum(Avg_Adj_R2*absum)/sum(absum),0,sum((1-Avg_Adj_R2)*absum)/sum(absum)) # Summary Averaged
} else{ # If significant covariation found
	sumENV <- Adj_R2
	unexp <- 1-sumENV # Abiotic variation explained
	cr <- rbind(CR,matrix(0,1,p))
	BR <- cr + t(cr)
	BR2 <- BR^2 # Squared correlation matrix of pairwise species
	pair <- apply(BR2,2,max) # Max squared correlation value for each species
	ind <- apply(BR2,2,which.max) # index of species pair with max correlation
	sumBIO <- pair*unexp # Biotic variation computation
	sumUNX <- matrix(1,p,1)-sumENV-sumBIO # Remainder variation is unexplained
	summary <- cbind(sumENV,sumBIO,sumUNX) # Abiotic, Biotic, Unexplained Variation
	AVGsummary = cbind(sum(sumENV*absum)/sum(absum), sum(sumBIO*absum)/sum(absum),sum(sumUNX*absum)/sum(absum)) # Summary Averaged, weighted by abundance
}

result = list('CDM'=CDM,'X'=X,'Yhat'=Yhat,'Yres'=Yres,'BP'=BP,'SitesSel'=SitesSel,'NN'=NN,'PRSqYres'=PRSqYres,'AvgRSqYres'=AvgRSqYres,'CR'=CR,'CV'=CV,'R2'=R2,'Avg_R2'=Avg_R2,'Adj_R2'=Adj_R2,'Avg_Adj_R2'=Avg_Adj_R2,'R2_RDA'=R2_RDA,'Adj_R2_RDA'=Adj_R2_RDA,'summary'=summary,'AVGsummary'=AVGsummary)

print('PRSqYres =')
print(PRSqYres)
colnames(AVGsummary) <- c('Abiotic','Biotic','Unexplained')
print('Community Variation Partitioned')
print(AVGsummary)

return(result)
}