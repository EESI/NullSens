mvrStandard <- function(CDM,X,sites_sel,n,p,q) {	

B_est <- matrix(0,q,p)  # Regression Parameters
Yhat <- matrix(0,n,p) # Predicted (Fitted) Species Responses
Yres <- matrix(0,n,p) # Residual Species Responses
	
for (i in 1:p) { # For each species
	# Regression Parameter Estimation -- Normal Equation
	tryCatch({
		B_est[,i] <- solve(t(X[sites_sel[[i]],])%*%X[sites_sel[[i]],])%*%t(X[sites_sel[[i]],])%*%CDM[sites_sel[[i]],i]	
	}, error = function(ex) {
		paste('Inverse computation -- singular. Cannot compute regression parameters for Species', i) # solve fails due to singular (X'*X) given sites_sel[[i]]
	})
	
	rmnan <- which(B_est[,i]==NA)
	B_est[rmnan,i] <- 0
	
	Yhat[,i] <- X %*% B_est[,i] # Compute fitted values
	rmneg <- which(Yhat[,i]<=0) # Remove negative fitted abundances
	Yhat[rmneg,i] = 0
	Yres[,i] <- CDM[,i]-Yhat[,i] # Compute Residuals

}
return(list("Yhat"=Yhat,"Yres"=Yres,"B_est"=B_est))
}