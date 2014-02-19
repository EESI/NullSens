mvrRobust <- function() {	
	
	
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

}