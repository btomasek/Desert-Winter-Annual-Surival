library(Rcpp)
library(bindrcpp)
library(sandwich)
library(gmm)
library(Matrix)

source('Functions/survivalLikelihood.R')
source('Functions/simSurvivalFunction.R')
sourceCpp('Functions/cppWeightedRowSums.cpp')

#As an example, use the coefficients estimated for EVMU
BetaSim<- c(-4.98528049,  0.13378620, -0.15117765,  0.20951936 , 0.18524818, -0.02342517 )
AlphaSim<- c(0.8666072, 0.9996096, 0.9530009 ) # This species has 3 strong-chronic effects

simInterval<- 14 #biweekly collection
n.per.group<- 100 ## number of individuals per year observed
Xnames.All<- c('int','Tmax','SM','TxSM','hab', 'age')
chronic.X<- c('Tmax','SM','TxSM')
acute.X<- c('int','hab','int2')

# Load in the rds file that has the simulation predictors in it
simList<- readRDS('data/simulationData.rds')
# This data has been saved in a list of 2 objects: 1), the 3-dimensional array that corresponds to the predictors
# 1) The first dimension of the array corresponds to predictor, the second corresponds to individual, the third correponds to day of year
# 2) The maximum time of censusing corresponding to each row. For example, row 3 has a maximum census time of 226

dim(simList$simX) # 6x46x224. For simulation, only one row is needed per year. In practice, you will need a row for each individual
simList$maxT[3] # 226 is the maximum time for row 3

set.seed(101)
sim.census=TRUE ### do you want to simulate with interval censoring, if set to false will return ONLY daily survival matrix

survivalSimulation<- simulateSurvival(x=simList$simX,tmax=simList$maxT, n.per.group = n.per.group, 
                           census=sim.census, alpha.type='Different', beta=BetaSim, 
                           alpha=AlphaSim, X.chronic=c(2,3,4), X.acute=c(1,5,6), 
                           average.census.interval = simInterval)

#survSim is a list with 2 objects. The first is the interval censored simulated data. the second object is the underlying daily survival that generated the interval data

survTrue<- survivalSimulation[[2]]
survSim<- survivalSimulation[[1]]

#### now do a fit on the simulated data
N<- nrow(survSim)

# Below, wrapper functions provided for fitting the model. These functions are available in Functions/survivalLikelihood.R
# It can be helpful to take the code out of the functions for debugging purposes.

initial.fit <- initializeFit(survivalMatrix = survSim, X = survivalSimulation$sim.X, iterations = 10000) # Recommend doing at least 10k iterations on real data, where the initialization of Beta is unknown

# Good practice to save out this initial (no chronic parameter). This initial fit is used to initialize values on the chronic-effect fit
tmpFit<- list(Beta=initial.fit$Beta, Alpha=initial.fit$Alpha, covarianceBeta = initial.fit$covarianceBeta, likVector=initial.fit$likVector, ng = nrow(initial.fit$Beta))
save(tmpFit, file=paste0('Data/Sim_InitialFit.rdata'))

# This wrapper function for fitting model with chronic parameters requires an initial fit formatted as above to run. 
#Recommend >20-100k iterations on real data, as convergence can take a while. It can be helpful to do more than one run to assess convergence
chronicFit<- chronicSurvivalFit(initialFit=tmpFit, 
                                survivalMatrix=survSim, 
                                X=simList$simX, 
                                iterations=20000) 

finalFit<- list(Beta=chronicFit$Beta, Alpha=chronicFit$Alpha, 
                covarianceBeta=chronicFit$covarianceBeta, covarianceAlpha=chronicFit$covarianceAlpha, # Proposition covariances
                BetaSim=BetaSim, AlphaSim=AlphaSim, X=simList$simX, # true values used in simulation
                likVector=chronicFit$likVector, ng=chronicFit$ng,
                survMat=survSim, survTrue=survTrue) # These two matrices have the interval censored and true daily survival simulation

save(finalFit, file=paste0('Sim_FinalFit.rdata'))

# To continue a chain and assure convergence, we can just use the previous chronicFit as the starting point
continueFit<- chronicSurvivalFit(initialFit=finalFit, 
                                survivalMatrix=survSim, 
                                X=simList$simX, 
                                iterations=20000, Xnames=colnames(chronicFit$Beta)) 

save(continueFit, file=paste0('Sim_FinalFit_Continued.rdata'))