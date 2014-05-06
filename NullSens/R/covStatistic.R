covStatistic <-
function(matrix,sites_sel,mutual_reject) {

n <- nrow(matrix) # Number of Sites
p <- ncol(matrix) # Number of Species

w <- 0
AvgRSq <- 0
for (species1 in 1:(p-1)){
	for (species2 in (species1+1):p){
		  mutsel <- sites_sel[[species1]] & sites_sel[[species2]] # Choose mutually selected sites of both species in pair
		  if (sum(mutsel) <= n*0.05 | sum(mutsel) <= mutual_reject){
		  	cvar <- 0
		  }	else {
		  	options(warn=-1) # Turn warnings off: OK for cor to return NA
		  	COV = cov(matrix[mutsel,c(species1,species2)]) # Covariation
		  	if (is.na(COV[1,2]) == TRUE | COV[1,2] == 0){
           		cvar <- 0
           } else {
            	cvar <- COV[1,2]; # Covariation
           } 
           options(warn=0) # Warnings back on
         }
         
         AvgRSq <- AvgRSq + cvar; # Index Calculation
         w <- w + 1; # Weights for averaging the index
	}
}
if (w == 0){
	index <- 0
} else {
	index <- AvgRSq/w  # index of test matrix
}

return("index"=index)
}
