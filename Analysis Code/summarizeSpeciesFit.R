sourceCpp('Functions/cppWeightedRowSums.cpp')

spec.run<- 'scba'

tmp.fit<- readRDS(paste0("Data/", spec.run,'_fit.rds'))

# summarize the fit from this species
# Coefficient Values
summary(tmp.fit$Beta)
# For 95% credible intervals...
t(apply(tmp.fit$Beta,2,quantile, c(0.025,0.5,0.975)))

# Chronic effects
summary(tmp.fit$alpha)
#95% credible intervals
t(apply(tmp.fit$alpha,2,quantile, c(0.025,0.5,0.975))) #only one strong chronic predictor


# Below function will if given a object like tmp.fit will return a plot of the predictions 
# Function will also return 95% credible intervals for cumulative survival probability across all individuals
species_prediction<- plotSurvivalFit(tmp.fit, nsim=10)

# Calculate the % chronic vs. acute. Here, we stratify by the habitat type

X<- tmp.fit$X
Beta<- tmp.fit$Beta
alpha = tmp.fit$alpha
maxCensus<- tmp.fit$maxCensusTime
ng<- nrow(Beta)
hab = tmp.fit$X[5,,1] # The 5th index is habitat (binary factor), we can figure out which individuals are in shrub vs. light by taking first column

# need to specify the indices that match the chronic versus acute effects,
# Note that unobserved weather time series can be accomodated by changing X and maxRecensus
# Because algorithm relies on draws from posterior distribution, may have slightly different results from run to run
lightHabitat-  partitionSurvivalEffectsP.Abs(X = X[,hab==1,], Beta, alpha, 
                                                maxRecensus  = maxCensus[hab==1], 
                                                npred=100, up.to.day=100,
                                             X.chronic=c(2,3,4), X.acute=c(1,5,6))

shrubHabitat<- partitionSurvivalEffectsP.Abs(X = X[,hab==0,], Beta, alpha, 
                                               maxRecensus = maxCensus[hab==0], 
                                               npred=100, up.to.day=100,
                                              X.chronic=c(2,3,4), X.acute=c(1,5,6))

# Get summary values for the overall effects with or without chronic effects
mean(lightHabitat$absolutelyMortalityDiffMean)
lightHabitat$absolutelyMortalityDiffSD

mean(shrubHabitat$absolutelyMortalityDiffMean)
shrubHabitat$absolutelyMortalityDiffSD

# Plot differences, one individual. Note that all individuals within a year will have the same curves
plot.individual<- 100
par(bty='n')
plot(lightHabitat$cumulativeToday[plot.individual,],type='l',col='blue',lwd=4,ylim=c(0,1),xlim=c(0,sum(!is.na(lightPartition$cumulativeToday[plot.individual,]))), ylab='Cumulative survival probability',xlab='Days since germination')
lines(lightHabitat$cumulativeAll[plot.individual,],type='l',col='red',lwd=4)
legend('bottomleft', c('Acute only', 'Acute & chronic'), bty='n', col=c('blue','red'),lwd=4,cex=2)
