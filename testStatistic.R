testStatistic <- function(matrix,sites_sel,n,p,q,test_stat_type,mutual_reject) {

w <- 0
AvgRSq <- 0
CR <- matrix(0,p-1,p)
CV <- matrix(0,p-1,p)
for (species1 in 1:(p-1)){
	for (species2 in (species1+1):p){
		  temp <- sites_sel[[species1]] & sites_sel[[species2]] # Choose mutually selected sites of both species in pair
		  if (sum(temp) <= n*0.05 | sum(temp) <= mutual_reject){
		  	cvar <- 0
		  }	else {
		  	options(warn=-1) # Turn warnings off: OK for cor to return NA
		  	COV = cov(matrix[temp,c(species1,species2)]) # Covariation
		  	CORR = cor(matrix[temp,c(species1,species2)]) # Correlation
		  	if (is.na(COV[1,2]) == TRUE | COV[1,2] == 0){
           		cvar <- 0
           		CV[species1,species2] <- cvar
           		CR[species1,species2] <- 0
           } else {
            	cvar <- COV[1,2];
            	CV[species1,species2] <- cvar # Keep track of pairwise covariance
           		CR[species1,species2] <- CORR[1,2] # Tracking of pairwise correlation
           } 
           options(warn=0) # Warnings back on
         }
         if (test_stat_type == 1) {
         AvgRSq <- AvgRSq + cvar^2; # Index Calculation
         w <- w + abs(cvar); # Weights for averaging the index
         }
         else print('Invalid Test Statistic Selection') 
         
	}
}
if (w == 0){
	index <- 0
} else {
	index <- AvgRSq/w  # index of test matrix
}

return(list("index"=index,"CR"=CR,"CV"=CV))
}