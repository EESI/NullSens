testStatistic <- function(matrix,sites_sel,n,p,q,test_stat,mutual_reject) {

w <- 0
AvgRSq <- 0
CR <- matrix(0,p-1,p)
CV <- matrix(0,p-1,p)
for (species1 in 1:(p-1)){
	for (species2 in (species1+1):p){
		  mutsel <- sites_sel[[species1]] & sites_sel[[species2]] # Choose mutually selected sites of both species in pair
		  if (sum(mutsel) <= n*0.05 | sum(mutsel) <= mutual_reject){
		  	cvar <- 0
		  }	else {
		  	options(warn=-1) # Turn warnings off: OK for cor to return NA
		  	COV = cov(matrix[mutsel,c(species1,species2)]) # Covariation
		  	CORR = cor(matrix[mutsel,c(species1,species2)]) # Correlation
		  	if (is.na(COV[1,2]) == TRUE | COV[1,2] == 0){
           		cvar <- 0
           		CV[species1,species2] <- cvar
           		CR[species1,species2] <- 0
           } else {
            	cvar <- COV[1,2];
            	ccorr <- CORR[1,2]
            	CV[species1,species2] <- cvar # Keep track of pairwise covariance
           		CR[species1,species2] <- ccorr # Tracking of pairwise correlation
           } 
           options(warn=0) # Warnings back on
         }
         
         if (test_stat[1] == 1) {
        	STAT = abs(cvar);
         } else if (test_stat[1] == 2){
         	STAT = cvar^2
         } else if (test_stat[1] == 3){
         	STAT = abs(ccorr)
         } else if (test_stat[1] == 4){
         	STAT = ccorr^2
         } else stop('Invalid Test Statistic Selection') 
         
         if (test_stat[2] == 1){
			WEIGHT = abs(cvar)
         } else if (test_stat[2] == 2) {
         	WEIGHT = abs(ccorr)
         } else if (test_stat[2] == 3){
         	WEIGHT = sum(mutsel)
         } else if (test_stat[2] == 4){
         	WEIGHT = 1
         } else stop('Invalid Test Statistic Weight') 
         
         AvgRSq <- AvgRSq + STAT*WEIGHT; # Index Calculation
         w <- w + WEIGHT; # Weights for averaging the index
	}
}
if (w == 0){
	index <- 0
} else {
	index <- AvgRSq/w  # index of test matrix
}

return(list("index"=index,"CR"=CR,"CV"=CV))
}