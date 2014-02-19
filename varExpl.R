varExpl <- function() {

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

}