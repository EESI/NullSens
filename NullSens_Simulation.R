NullSens_Simulation <- function(n,p,q,N,M,CD,Type,num_exp=100,select = TRUE,reg_method="robust",null_reps=200,test_stat = c(1,1),mutual_reject=7,alpha=0.05) {

source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/SpeciesResponseModel.R")
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/NullSens.R")

Noise = c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2) # Noise Parameters
Magnitude = Noise[N]*c(1,3,5,10,20) # Covariation Magnitude Parameters

# Preallocation
COM_summary <- matrix(0,num_exp,3)
pos_cov_p_values <- c(rep(0,num_exp))
neg_cov_p_values <- c(rep(0,num_exp))
p_values <- c(rep(0,num_exp))

tic() # Start the stopwatch
for (exp in 1:num_exp) {
	# Generate Simulations using Species Response Model
	SRM <- SpeciesResponseModel(100,20,2,2,3,0,2,Noise,Magnitude)
	# Run NullSens on the Simulated Data
	results <- NullSens(SRM$CDM,SRM$X)
	# Store the results
	COM_summary[exp,] <- results$COM_summary
	pos_cov_p_values[exp] <- results$COM_variation_type$p_value_pos
	neg_cov_p_values[exp] <- results$COM_variation_type$p_value_neg
	p_values[exp] <- results$p_value
}

# Number of communities with significant iteraction over total experiments
avg_detection_rate <- length(which(p_values <= alpha))/num_exp
# Number of communities with significant pos/neg covariation (considering only those where a significant interaction was found)
if (avg_detection_rate > 0) {
	num_sig_pos <- length(which(pos_cov_p_values[which(p_values <= alpha)] <= alpha))/length(which(p_values <= alpha))
	num_sig_neg <- length(which(neg_cov_p_values[which(p_values <= alpha)] <= alpha))/length(which(p_values <= alpha))
} else {
	num_sig_pos = 0; num_sig_neg = 0
}
avg_detection_type <- list(num_sig_pos,num_sig_neg)

'run_time' = toc() # Total Elapsed Time

# Summarize all input parameters
Parameters <- list('n'=n,'p'=p,'q'=q,'N'=N,'M'=M,'CD'=CD,'Type'=Type,'num_exp'=num_exp,'select'=select,'reg_method'=reg_method,'null_reps'=null_reps,'test_stat'=test_stat,'mutual_reject'=mutual_reject,'alpha'=alpha)

# Return results from all experiemnts
return(list('Paramters' = Parameters, 'COM_summary' = COM_summary, 'p_values' = p_values, 'avg_detection_rate'= avg_detection_rate, 'avg_detection_type' = avg_detection_type, 'run_time'=run_time))
}

# Helper Functions
# Stopwatch function in R
# http://stackoverflow.com/questions/1716012/stopwatch-function-in-r
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   #print(toc - tic)
   invisible(toc)
}