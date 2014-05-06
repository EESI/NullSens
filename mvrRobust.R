mvrRobust <- function(CDM,X,sites_sel,n,p,q) {	

B_est <- matrix(0,q,p)  # Regression Parameters
Yhat <- matrix(0,n,p) # Predicted (Fitted) Species Responses
Yres <- matrix(0,n,p) # Residual Species Responses
	
# get MASS rlm for robust regression

for (i in 1:p) { # For each species
	# Regression Parameter Estimation -- Robust Regression using BiSquare Function
	result <- rlm(X[sites_sel[[i]],],CDM[sites_sel[[i]],i],,psi.bisquare,method="M",maxit=50)
	B_est[,i] <- as.numeric(result$coefficients)

	rmnan <- which(B_est[,i]==NA)
	B_est[rmnan,i] <- 0

	Yhat[sites_sel[[i]],i] <- X[sites_sel[[i]],] %*% B_est[,i] # Compute fitted values
	rmneg <- which(Yhat[,i]<=0) # Remove negative fitted abundances
	Yhat[rmneg,i] = 0
	Yres[,i] <- CDM[,i]-Yhat[,i] # Compute Residuals

}
return(list("Yhat"=Yhat,"Yres"=Yres,"B_est"=B_est))
}