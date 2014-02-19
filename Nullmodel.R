Nullmodel <- function() {

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

}