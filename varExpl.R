varExpl <- function(CDM,PRSqYres,n,p,q,sites_sel,CR,Avg_Adj_R2,Adj_R2) {

# Return Abiotic, Biotic, Unexplained Variances
absum <- c(rep(0,p)) # Compute raw abundance totals for each species
for (i in 1:p){
	absum[i] <- sum(CDM[sites_sel[[i]],i])
}

if (PRSqYres > 0.05){ # If no significant covariation found
	INDsummary <- cbind(Adj_R2,matrix(0,p,1),1-Adj_R2) # Abiotic, Biotic, Unexplained, per species
	COMsummary <- cbind(Avg_Adj_R2,0,1-Avg_Adj_R2) # Abiotic, Biotic, Unexplained, community

} else{ # If significant covariation found
	# Individual species variation partitioning
	sumENV <- Adj_R2
	unexp <- 1-sumENV # Abiotic variation explained
	cr <- rbind(CR,matrix(0,1,p))
	BR <- cr + t(cr)
	BR2 <- BR^2 # Squared correlation matrix of pairwise species
	pair <- apply(BR2,2,max) # Max squared correlation value for each species
	ind <- apply(BR2,2,which.max) # index of species pair with max correlation
	sumBIO <- pair*unexp # Biotic variation computation
	sumUNX <- matrix(1,p,1)-sumENV-sumBIO # Remainder variation is unexplained
	INDsummary <- cbind(sumENV,sumBIO,sumUNX) # Abiotic, Biotic, Unexplained Variation
		
	# Community variation partitioning
	ENV <- Avg_Adj_R2 # Environmental variation explained
	unexp <- 1-ENV # Abiotic variation explained
	cr <- rbind(CR,matrix(0,1,p))
	BR <- cr + t(cr)
	BR2 <- BR^2 # Squared correlation matrix of pairwise species
	pair <- apply(BR2,2,max) # Max squared correlation value for each species
	ind <- apply(BR2,2,which.max) # index of species pair with max correlation
	BIO <- sum(pair*unexp*absum)/sum(absum) # Biotic variation computation
	UNX <- 1-ENV-BIO # Remainder variation is unexplained
	COMsummary = cbind(ENV, BIO, UNX)
}
return(list(INDsummary,COMsummary))
}
