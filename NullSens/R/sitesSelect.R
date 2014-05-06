sitesSelect <-
function(CDM,X) {

n <- nrow(matrix) # Number of sites
p <- ncol(matrix) # Number of species
q <- ncol(X) # Number of abotic factors + intercept

sites_sel <- vector('list',p) # Stores Sites Selected for each Species

# Site Selection
for (i in 1:p){ # For each species
	posAb <- which(CDM[,i] > 0) # Choose sites with positive abundance for theshold
	remove <- NA # initialize sites for removal
	centroid <- c(rep(0,q)) # Median of each gradient
	thresh <- c(rep(0,q)) # Threshold for site removal
	
	for (grad in 2:q){ # For each abiotic factor
		centroid[grad] <- quantile(X[posAb,grad],0.5) # Centroid (Median) of positive abundance
		distance_pos <- as.matrix(dist(c(centroid[grad],X[posAb,grad]))) # distances of positive abundances to centroid
		thresh[grad] <- mean(distance_pos[-1,1]) + sd(distance_pos[-1,1]) # threshold for site selection
	
		distance_all <- as.matrix(dist(c(centroid[grad],X[,grad]))) # distances of all sites (positive and zero) to centroid
		remove_temp <- which(distance_all[-1,1] > thresh[grad]) # remove sites beyond threshold from further analysis
		remove <- c(remove,remove_temp) # keep a list of all sites removed for each species
	}
	
	keep = setdiff(1:n,remove) # Keep sites not marked for removal
	sites_sel[[i]] <- c(rep(FALSE,n))
	if (length(keep) > 3){ # remove species with < 4 selected sites from analysis
		sites_sel[[i]][keep] = TRUE	
	}
}
return(sites_sel)
}