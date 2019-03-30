sourceCpp('Functions/cppWeightedRowSums.cpp')

spec.run<- 'pehe'

tmp.fit<- readRDS(paste0("Data/", spec.run,'_fit.rds'))

# Below function will if given a object like tmp.fit will return a plot of the predictions 
# Function will also return 95% credible intervals for cumulative survival probability across all individuals
species_prediction<- plotSurvivalFit(tmp.fit, nsim=10) # These are posterior simulated predictors, larger nsim is more accurate but can be time-consuming

# Calculate the % chronic vs. acute. Here, we stratify by the habitat type

X<- tmp.fit$X
Beta<- tmp.fit$Beta
alpha = tmp.fit$alpha
maxCensus<- tmp.fit$maxCensusTime
ng<- nrow(Beta)
hab = tmp.fit$X[5,,1] # The 5th index is habitat (binary factor), we can figure out which individuals are in shrub vs. light by taking first column

# need to specify the indices that match the chronic versus acute effects
tmp.hab.Abs<-  partitionSurvivalEffectsP.Abs(X = X[,hab==1,], Beta, alpha, maxRecensus  = maxCensus[hab==1], npred=2,up.to.day=100,
                                             X.chronic=c(2,3,4), X.acute=c(1,5,6))
tmp.nohab.Abs<- partitionSurvivalEffectsP.Abs(X = X[,hab==0,], Beta, alpha, maxRecensus = maxCensus[hab==0], npred=2,up.to.day=100,
                                              X.chronic=c(2,3,4), X.acute=c(1,5,6))

# You can also calculate the cumulative survival probabilities with and without the chronic effects
