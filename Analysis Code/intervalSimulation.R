library(Rcpp)
library(bindrcpp)
library(sandwich)
library(gmm)
source('Functions/survivalLikelihood.R')
sourceCpp('Functions/cppWeightedRowSums.cpp')
library(Matrix)
source('Functions/simSurvivalFunction.R')

#As an example, use the coefficients estimated for EVMU
BetaSim<- c(-4.98528049,  0.13378620, -0.15117765,  0.20951936 , 0.18524818, -0.02342517 )
AlphaSim<- c(0.8666072, 0.9996096, 0.9530009 )

simInterval<- 14 #biweekly collection
n.per.group<- 100 ## number of n per year observed
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

survSim<- simulateSurvival(x=simList$simX,tmax=simList$maxT, n.per.group = n.per.group, 
                           census=sim.census, alpha.type='Different', beta=BetaSim, 
                           alpha=AlphaSim, X.chronic=c(2,3,4), X.acute=c(1,5,6), 
                           average.census.interval = simInterval)

#survSim is a list with 2 objects. The first is the interval censored simulated data. the second object is the underlying daily survival that generated the interval data

survTrue<- survSim[[2]]
survSim<- survSim[[1]]

#### now do a fit on the simulated data
N<- nrow(survSim)

now = Sys.time()
initial.fit <- initializeFit(survSim, simList$simX, 2000) # Approximately 1/2 second per iteration. Recommend doing at least 5k iterations on real data, where the initialization of Beta is unknown
Sys.time() - now

tmpFit<- list(Beta=initial.fit$Beta, alpha=initial.fit$Alpha, covarianceBeta = initial.fit$covarianceBeta, likVector=initial.fit$likVector, ng = nrow(initial.fit$Beta))
save(tmpFit, file=paste0('Sim_InitialFit.rdata'))

chronicFit<- chronicSurvivalFit(tmpFit, survSim, simList$simX, 10000) #Approximately 1 second per iteration. Recommend >20-100k iterations on real data, as convergence can take a while

finalFit<- list(Beta=chronicFit$Beta, alpha=chronicFit$Alpha, covb=chronicFit$covb, cova=chronicFit$cova, 
                betasim=BetaSim,alphasim=chronicFit$Alpha, X=simList$simX, 
                likVector=chronicFit$likVector, survMat=survSim, survTrue=survTrue)

save(finalFit, file=paste0('Sim_FinalFit.rdata'))

