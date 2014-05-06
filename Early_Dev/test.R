#############################################################################
# Package it up
#############################################################################

# Need to include MASS and censReg packages
rm(list = ls()) # Clear the workspace
library(MASS)
library(censReg)

source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/SpeciesResponseModel.R")
Noise = c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2) # Noise Parameters
Magnitude = Noise[2]*c(1,3,5,10,20) # Covariation Magnitude Parameters
SRM <- SpeciesResponseModel(100,20,2,2,3,2,0,Noise,Magnitude)
CDM = SRM$CDM; X = SRM$X

#set.seed(100)
source("/Users/Dizzy/Desktop/NullSens_R/NullSens-R/NullSens.R")
result <- NullSens(CDM,X)


#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/sitesSelect.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/mvrRobust.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/mvrTobit.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/mvrStandard.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/coeffDet.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/nullModel.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/testStatistic.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/covariationType.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/varExpl.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/NullSens_Simulation.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/SpeciesResponseModel.R")
#source("/Users/Dizzy/Desktop/NullSens_R/NullSens/NullSens.R")


#############################################################################
# Simulation Sweep

#rm(list = ls()) # Clear the workspace


#n = 100; p = 20; q=2; N=6; M=3; CD=2; Type=2; num_exp = 10;
#results <- NullSens_Simulation(n,p,q,N,M,CD,Type,num_exp)