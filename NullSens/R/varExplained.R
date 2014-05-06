varExplained <-
function(CDM,p_value,sites_sel,CR,Avg_Adj_R2,Adj_R2,alpha) {

p <- ncol(CDM) # Number of Species

# Return Abiotic, Biotic, Unexplained Variances
absum <- c(rep(0,p)) # Compute raw abundance totals for each species
for (i in 1:p){
	absum[i] <- sum(CDM[sites_sel[[i]],i])
}

if (p_value > alpha){ # If no significant covariation found
	env = Adj_R2; bio = matrix(0,p,1); unx = 1-Adj_R2
	ENV = Avg_Adj_R2; BIO = 0; UNX = 1-Avg_Adj_R2
	IND_summary <- cbind(Adj_R2,matrix(0,p,1),1-Adj_R2) # Abiotic, Biotic, Unexplained, per species
	COM_summary <- cbind(Avg_Adj_R2,0,1-Avg_Adj_R2) # Abiotic, Biotic, Unexplained, community

} else { # If significant covariation found
	# Individual species variation partitioning
	env <- Adj_R2
	unx <- 1-env # Abiotic variation explained
	cr <- rbind(CR,matrix(0,1,p))
	BR <- cr + t(cr)
	BR2 <- BR^2 # Squared correlation matrix of pairwise species
	pair <- apply(BR2,2,max) # Max squared correlation value for each species
	ind <- apply(BR2,2,which.max) # index of species pair with max correlation
	bio <- pair*unx # Biotic variation computation
	unx <- matrix(1,p,1)-env-bio # Remainder variation is unexplained
	IND_summary <- cbind(env,bio,unx) # Abiotic, Biotic, Unexplained Variation
		
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
	COM_summary <- cbind(ENV, BIO, UNX)
}
return(list('IND_summary' = IND_summary,'COM_summary' = COM_summary))
}
