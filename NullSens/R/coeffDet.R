coeffDet <-
function(CDM,X,Yhat) {

# OUTPUTS
# R2 -- Coefficient of Multiple Determination, per species
# Adj_R2 -- Adjusted R2, per species
# Avg_R2 -- Community Averaged R2
# Avg_Adj_R2 -- Community Averaged Adjusted R2

n <- nrow(CDM) # Number of Sites
p <- ncol(CDM) # Number of Species
q <- ncol(X)-1 # Number of Abiotic Factors (minus intercept column)

R2 <- c(rep(0,p))
Adj_R2 <- c(rep(0,p))
for (i in 1:p) { # For each species
	# R2, per species - Same as weighted community. No trace explicitly needed since only vectors.
	R2[i] = (Yhat[,i] %*% Yhat[,i])/(CDM[,i] %*% CDM[,i])

	# Adjusted R2, per species
	Adj_R2[i] = 1-(1-R2[i])*(n-1)/(n-q-1)
}

# Weighted (by Abundance) Community Averaged R2
Avg_R2 <- sum(diag(t(Yhat) %*% Yhat))/sum(diag(t(CDM) %*% CDM))

# Weighted (by Abundance) Community Averaged Adjusted R2
Avg_Adj_R2 <- 1-(1-Avg_R2)*(n-1)/(n-q-1)

return(list("R2"=R2, "Adj_R2"=Adj_R2, "Avg_R2"=Avg_R2, "Avg_Adj_R2"=Avg_Adj_R2))
}
