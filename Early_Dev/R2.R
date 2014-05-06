R2 <- function() {}

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

}