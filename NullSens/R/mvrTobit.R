mvrTobit <-
function(CDM,X,sites_sel) {	

n <- nrow(matrix) # Number of sites
p <- ncol(matrix) # Number of species
q <- ncol(X) # Number of abotic factors + intercept

B_est <- matrix(0,q,p)  # Regression Parameters
Yhat <- matrix(0,n,p) # Predicted (Fitted) Species Responses
Yres <- matrix(0,n,p) # Residual Species Responses
	
# get censReg for tobit regression

for (i in 1:p) { # For each species
	# Regression Parameter Estimation -- Tobit (Censored) Regression
	# Tobit automatically computes intercept so exclude column of 1s
	data = data.frame(X[sites_sel[[i]],2:ncol(X)], cdm=CDM[sites_sel[[i]],i])
	result <- censReg(as.formula("cdm ~ ."), left=0, right=Inf, data)
	B_est[,i] <- as.matrix(as.numeric(coef(result))[1:ncol(X)]) # Omit last column (logSigma)
	
	rmnan <- which(B_est[,i]==NA)
	B_est[rmnan,i] <- 0

	Yhat[,i] <- X %*% B_est[,i] # Compute fitted values
	rmneg <- which(Yhat[,i]<=0) # Remove negative fitted abundances
	Yhat[rmneg,i] = 0
	Yres[,i] <- CDM[,i]-Yhat[,i] # Compute Residuals

}
return(list("Yhat"=Yhat,"Yres"=Yres,"B_est"=B_est))
}
