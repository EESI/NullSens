nullModel <-
function(matrix,sites_sel) {

p = ncol(matrix) # Number of Species

# Null model - generates random matrices
rand_matrix <- matrix # For each randomization always start with original Yres
for (species in 1:p){
	permNon <- which(sites_sel[[species]]==TRUE)
	neworderNon <- sample(permNon, length(permNon), replace = FALSE, prob = NULL)
	rand_matrix[permNon,species] <- matrix[neworderNon,species]
}

return(rand_matrix)
}
