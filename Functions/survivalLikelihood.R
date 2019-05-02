library(mvtnorm)
library(truncnorm)

extractProbability<- function(survData, X,Beta, i){
  N = length(survData)
  pmat <- rep(0, N) #holds daily probabilities of death
  p <- dim(X[[1]])[2]
  XB<- matrix(unlist(X[[i]]), ncol=p)%*%t(Beta)
  p.tmp <- 1-exp(-exp(XB))
  return(p.tmp)
}

proposeBetaAccumulation<- function(beta, mu0, S0, cov=1E-8*diag(length(beta)),maxBeta=Inf){  ### need to truncate so the logit doesn't completely blow up
  betanew<-beta + rmvnorm(n=1,sigma=cov, method='svd')   # Draw from symmetric MV t-distv 
  # Log Prior
  lpold<-dmvnorm(t(beta),mu0,S0,log=T)    # Old Log Prior
  lpnew<-dmvnorm(betanew,mu0,S0,log=T)    # New Log Prior
  return(list('betaold'=beta, 'betanew'=betanew, 'likOld'=lpold, 'likNew'=lpnew)) 
}

proposeBeta<- function(beta, mu0, S0, cov=1E-8*diag(length(beta))){ 
  betanew<-beta + rmvnorm(n=1,sigma=cov, method='chol')   # Draw from symmetric MV t-distv 
  # Log Prior
  lpold<-dmvnorm(t(beta),mu0,S0,log=T)    # Old Log Prior
  lpnew<-dmvnorm(betanew,mu0,S0,log=T)    # New Log Prior
  return(list('betaold'=beta, 'betanew'=betanew, 'likOld'=lpold, 'likNew'=lpnew)) 
}

adaptiveCovariance<- function(B.bar, Beta.t, Ct,t, t0, epsilon=1*10^(-6), C0){
  if(t<=t0){
    return(C0)
  }else{
    s = (2.4)^2 / ncol(Beta)
    X1 = B.bar[1,]
    X2 = B.bar[2,]
    Ct.1 = (t-1)/t * Ct + s/t * (t*(X1%*%t(X1)) - (t+1)*(X2 %*%t(X2)) + Beta.t%*%t(Beta.t) + diag(epsilon,ncol(Beta))) ### should be diag epsilon    
    return(Ct.1)
  }
}

updateBeta<- function(beta, covb, r1, r2, mu0, S0){
  betanew<- as.matrix(t(beta) + rmvnorm(n=1,sigma=covb), method="chol")  # Draw from MVN
  lpold<-dmvnorm(t(beta),mu0,S0,log=T)  	# Old Log Prior
  lpnew<-dmvnorm(betanew,mu0,S0,log=T)		# New Log Prior  
  return(list('betaOld' = beta, 'betaNew' = betanew, 'oldLik' = lpold, 'lpnew' = lpnew))
}


ilogit <- function(X){
  exp(X)/(1+exp(X))
}

prodInv<- function(X){
  prod(1 - X) 
}

library(dplyr)
survLikelihood<- function(Y, p, ID){
  df.tmp = data.frame('P.tmp'=p, "Y" = Y, "ID"=ID)
  grp<- group_by(df.tmp, ID, Y) # set up the grouping
  
  use.dplyr<-summarise(grp,  "term"= prod(P.tmp)) #set up aggregation by groups  
  use.dplyr<-collect(use.dplyr)
  use.dplyr$term[use.dplyr$Y==1] = 1 - use.dplyr$term[use.dplyr$Y==1]
  
  int = aggregate(use.dplyr$term, by=list(use.dplyr$ID), prod)
  loglik<-  sum(log(int[,2]))
  return(loglik)
}

blank.names <- function(dat){
  for(i in 1:ncol(dat)){
    names(dat)[i] <- paste(rep(" ",i),collapse="")
  }
  return(dat) 
}


library(Matrix)
updatePropAdaptive<- function(covBeta, acceptRate, iters,mult=1, adjust=(-0.3), adaptiveAdjust=T,nearestPD=T){
  tiny = 1E-8
  p= ncol(covBeta)
  covBeta<- covBeta #+ diag(nrow(covBeta)) * tiny ### adding the tiny in will make the matrix positive definite if chain is not moving much, set tiny to 0 if not wanted 
  if(length(p) == 0){
    p<- length(covBeta)
  }
  if(adaptiveAdjust){ ### attempts to smooth out adjustment
    if(acceptRate <= 0.35){
      adjust = -0.2 + (  acceptRate) * -(0.8 / 0.35)
    }else{
      adjust = -0.8/0.65 + (  acceptRate) * (0.8 / 0.65) - 0.2
    }
    
  }
  if(acceptRate<0.35){
    covBeta<- mult  * (1 - iters^(adjust)) * (2.4^2 / p) * covBeta 
  }
  if(acceptRate>=0.35){
    covBeta<- mult  * (1 + iters^(adjust)) * (2.4^2 / p) * covBeta 
  }
  cb<- covBeta + diag(nrow(covBeta)) * tiny ### will make matrix positive definite usually
  if(nearestPD){
    cb<- nearPD(cb)$mat
  }
  return(cb)
}
library(tmvtnorm)

#### prediction 
plotSurvivalFit<- function(speciesFitList, nsim=100){
  print('Plotting posterior distribution simulated predictions, this may take a few minutes...')
  tmp.fit<- speciesFitList
  ng<- nrow(tmp.fit$Beta)
  N<- nrow(tmp.fit$survMat)
  
  individualMaxCensus<- tmp.fit$maxCensusTime
  
  npred<- nsim
  holdPred<- matrix(NA, nrow=N, ncol=npred)
  par(mfrow=c(1,1))
  for(S in 1:npred){  
    beta.samp<- tmp.fit$Beta[sample(round((3*ng/4)):ng,1),]
    alpha.samp<- tmp.fit$alpha[sample(round((3*ng/4)):ng,1),]
    dL.predict<- predict.Survival(as.vector(beta.samp), c(alpha.samp), bigX.pred=tmp.fit$X,tRow=tmp.fit$maxCensusTime) # returns a matrix of daily probabilites
    
    tmp.p = t(apply(dL.predict,1,cumprod)) #Translate survival to cumulative survival
    
    rand = runif(N)
    for(t in 1:length(rand)){
      tmpD= which(tmp.p[t,1:individualMaxCensus[t]] <= rand[t])[1]
      if(is.na(tmpD))tmpD=NA
      holdPred[t,S] = tmpD
    }
  }
  actT<- calcActualSurvival(tmp.fit$survMat)
  
  t1<- actT[,2]
  t2<- actT[,3]
  maxTime<- tmp.fit$maxCensusTime[which(is.na(t1))]
  d<- is.na(t1) # did not die, so !d is death
  t1[is.na(t1)]<- maxTime
  t2[is.na(t2)]<- maxTime
  
  test1<- my.KM(t1, !d)
  test2<- my.KM(t2, !d)
  
  holdKM<- matrix(NA, nrow=ncol(holdPred), ncol=max(tmp.fit$maxCensusTime))
  for(j in 1:ncol(holdPred)){
    pp<- holdPred[,j]
    d.tmp<- is.na(pp)
    pp[is.na(pp)]<- tmp.fit$maxCensusTime[is.na(pp)]
    tmp.km<- my.KM(pp, !d.tmp)
    holdKM[j,1:length(tmp.km)]<- tmp.km
    if(length(tmp.km)<max(tmp.fit$maxCensusTime)){
      tmp.max<- length(tmp.km)
      holdKM[j,(tmp.max+1):max(tmp.fit$maxCensusTime)]<- holdKM[j,tmp.max]
    }
  }
  quant.KM<- apply(holdKM,2,quantile,c(0.025,0.5,0.975))
  
  
  plot(1:length(test1), test1,type='l', ylim=range(quant.KM),xlab='Days since germination',ylab='Proportion surviving individuals',cex.axis=1.4,cex.lab=1.3)
  lines(1:length(test2), test2)
  
  plotConfidenceRegion(1:ncol(holdKM), quant.KM[1,], quant.KM[3,],col.border='grey38')
  lines(quant.KM[2,],col='grey45',lty=1,lwd=4)
  lines(1:length(test1), test1,lwd=3)
  lines(1:length(test2), test2,lwd=3)
  legend('bottomleft', c('"Observed" cumulative survival', 'Median predicted', '95% Prediction Interval'), col=c(1,'grey45','grey'), lwd=3)
  abline(v= which(colSums(survMat != (0))  > round(0.90 * nrow(survMat)) )[1], col='red',lty=3,lwd=2)
  
  return(quant.KM)
}

predict.Survival<- function(betaSamp, alphaSamp, bigX.pred, tRow=newtRow, X.chronic=c(2,3,4), X.acute=c(1,5,6)){ #### alpha.type controls which type to do, {'Different', 'Same', 'One', 'Zero'}
  
  nt<- dim(bigX.pred)[3]
  ni<- dim(bigX.pred)[2]
  dL<- matrix(0, nrow=ni, ncol=nt)
  dL.static<- matrix(0, nrow=ni, ncol=nt)
  dL.acc<- array(0, dim=c(ni, nt, length(X.chronic)))
  
  if(alpha.type=='Different'){ 
    for(j in 1:nt){
      dL.static[,j]<-  t(bigX.pred[X.acute,,j])%*%c(betaSamp[X.acute])
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- bigX.pred[X.chronic[h],,j] * betaSamp[X.chronic][h]
      }
    }
    dL<- dL.static
    for(h in 1:length(X.chronic)){
      dL<- dL + cppWeightedRowSums(dL.acc[,,h], tRow, alphaSamp[h])
    }
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
  }
  
  if(alpha.type=='Zero'){ 
    for(j in 1:nt){
      dL.static[,j] <-  t(bigX.pred[,,j])%*%betaSamp
    }
    dL2 <-   1 - ( 1-exp(-exp(dL.static)))  # survival prob
  }
  
  return(dL2)
}

plotConfidenceRegion<- function(x, ylo, yhi,col.reg='grey',col.border='black',alpha=1){
  xx<- c(x,rev(x))
  yy<- c(ylo, rev(yhi))
  col.plot<- t(col2rgb(col.reg) )
  col.plot<- rgb(col.plot,alpha=alpha*255 ,maxColorValue=255)
  polygon(xx, yy, col=col.plot, border=col.border)
}


# Given a survival matrix, return the KM estimates
calcActualSurvival<- function(allIndividuals){
  N<- nrow(allIndividuals)
  individualList<- list()
  for( i in 1:N){
    individualList[[i]]<- allIndividuals[i,][allIndividuals[i,]!=(-1)] 
  }
  act.surv<- matrix(NA, nrow=N,ncol=3)
  for(i in 1:N){
    tt<- 1:length(individualList[[i]])
    act.surv[i,1]<- median(tt[individualList[[i]] == 1],na.rm=T)
    act.surv[i,2]<- min(tt[individualList[[i]] == 1],na.rm=T)
    act.surv[i,3]<- max(tt[individualList[[i]] == 1],na.rm=T)
  }
  act.surv[abs(act.surv)==Inf]<- NA
  return(act.surv)
}

my.KM<- function(x,d){
  tmax<- max(x,na.rm=T)
  surv.tmp<- rep(1, tmax)
  for(i in 2:tmax){
    ni<- length(x[x>=i])
    di<- sum(d[x==i])
    surv.tmp[i]<- (ni - di)/ni
  }
  return(cumprod(surv.tmp))
}

####
#### calculate end of season survival from the vector of predictions 
calcAcceptance<- function(betaVec){
  mean(1 - (diff(betaVec)==0))
}

conditionalMVN <- function(xx, mu, sigma, cindex){  
  
  # x and mu are vectors, cindex is vector index for conditional
  tiny <- min(diag(sigma))*.0001
  nm   <- length(mu)
  if(length(xx) != nm)stop('x and mu different length in conditionalMVN')
  
  xx <- matrix(xx,nrow=1)
  mu <- matrix(mu,nrow=1)
  
  testv <- try(chol(sigma[-cindex,-cindex]),T)
  if(inherits(testv,'try-error')){
    return( list(mu = numeric(0), vr = numeric(0)) )
  }
  
  sin <- chol2inv(testv)
  p1  <- sigma[cindex,-cindex]%*%sin
  
  mu1 <- mu[cindex] +  p1%*%(xx[-cindex] - mu[-cindex]) 
  vr1 <- sigma[cindex,cindex] - p1%*%sigma[-cindex,cindex]
  
  list(mu = mu1, vr = vr1)
}

tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}

.tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  
  tiny <- 10e-6
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] - tiny
  z
}

tnorm.mvt <- function(avec,muvec,smat,lo=rep(-Inf,length(avec)),
                      hi=rep(Inf,length(avec)),
                      whichSample=c(1:length(avec)),times=1){   
  
  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  # whichSample indicates which variables to sample
  
  if(length(lo) == 1)lo <- rep(lo,length(avec))
  if(length(hi) == 1)hi <- rep(hi,length(avec))
  
  for(j in 1:times){
    for(k in whichSample){
      
      tmp <- conditionalMVN(avec,muvec,smat,k)
      muk <- tmp$mu
      sgk <- tmp$vr
      
      if(length(muk) == 0)next
      
      avec[k] <- .tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
    }
  }
  avec
}

conditionalMVN <- function(xx, mu, sigma, cindex){  
  
  # x and mu are vectors, cindex is vector index for conditional
  
  tiny <- min(diag(sigma))*.0001
  nm   <- length(mu)
  if(length(xx) != nm)stop('x and mu different length in conditionalMVN')
  
  xx <- matrix(xx,nrow=1)
  mu <- matrix(mu,nrow=1)
  
  testv <- try(chol(sigma[-cindex,-cindex]),T)
  if(inherits(testv,'try-error')){
    return( list(mu = numeric(0), vr = numeric(0)) )
  }
  
  sin <- chol2inv(testv)
  p1  <- sigma[cindex,-cindex]%*%sin
  
  mu1 <- mu[cindex] +  p1%*%(xx[-cindex] - mu[-cindex]) 
  vr1 <- sigma[cindex,cindex] - p1%*%sigma[-cindex,cindex]
  
  list(mu = mu1, vr = vr1)
}

# Function to simulate survival given a covariate array x and a number of parameters
simulateSurvival<- function(x=NULL,tmax=150,n.per.group = 50, beta=NULL,alpha=NULL,alpha.type='Different',X.chronic=c(2,3,4),X.acute=c(1,5), census=TRUE,average.census.interval=7){
  yrs<- dim(x)[2]
  groups=yrs
  censusTimes<- list()
  if(census==T){
    for(y in 1:yrs){
      days.between<- rep(average.census.interval, 100)
      census.Days<- cumsum(days.between)
      census.Days<- census.Days[census.Days < tmax[y]]
      census.Days<- unique(c(1,census.Days, tmax[y]))
      censusTimes[[y]]<- census.Days
    }
  }else{
    for(y in 1:yrs){
      censusTimes[[y]]<- 1:tmax[y]
    }
  }
  if(alpha.type=='Different'){
    if(is.null(beta)){
      beta<- c(-5.5,0.25, -1, -3, 0.5)
    }
    
    if(is.null(alpha)){
      alpha<- c(.8, .4, .1)
    }
    dL<- dL.static <- matrix(0, nrow=dim(x)[2], ncol=max(tmax))
    dL.acc2<-dL.acc <- array(0, dim=c(nrow(dL), ncol(dL), length(X.chronic)))
    
    #List Xstatic, List Xacc, NumericVector Bstatic, NumericVector Bacc, int N, int NT, double alpha
    for(j in 1:ncol(dL)){
      dL.static[,j]<-  t(x[X.acute,,j])%*%beta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(x[X.chronic[h],,j]) * beta[X.chronic][h]
      }
    }
    dL<- dL.static
    for(h in 1:length(X.chronic)){
      dL<- dL + cppWeightedRowSums(dL.acc[,,h], max(tmax), alpha[h])
    }
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
    dL2[,1]<- 1
    holdRisk<- dL2
    ### now make sure dL2 is right and interval censored properly    
    sim.Individuals<- matrix(-1, nrow=n.per.group * yrs, ncol=ncol(dL))
    sim.Daily<- matrix(-1, nrow=n.per.group * yrs, ncol=ncol(dL))
    for(j in 1:nrow(dL2)){
      tmp.cens<- censusTimes[[j]]
      rand<- runif(n.per.group)
      tmp.tmax<- tmax[j]
      for(t in 1:length(rand)){
        tt<- which(cumprod(dL2[j,1:tmp.tmax]) <= rand[t])[1]
        if(census == TRUE){
          if(!is.na(tt)){
            sim.Daily[n.per.group*(j-1)+t,1:(tt-1)]<- 0
            sim.Daily[n.per.group*(j-1)+t,tt]<- 1
          }else{
            sim.Daily[n.per.group*(j-1)+t,1:tmp.tmax]<- 0
          }
          if(!is.na(tt)){
            if(tt == max(tmp.cens)){
              ll<- length(tmp.cens)
              start.int<- ll-1
              end.int<- ll
              sim.Individuals[n.per.group*(j-1)+t,1:tmp.cens[start.int]]<- 0 
              sim.Individuals[n.per.group*(j-1)+t,(1+tmp.cens[start.int]):tmp.cens[end.int]]<- 1 
            }else{
              end.int<- which(tmp.cens >= tt)
              end.int<- end.int[1]
              start.int<- which(tmp.cens < tt)
              start.int<- start.int[length(start.int)]
              sim.Individuals[n.per.group*(j-1)+t,1:tmp.cens[start.int]]<- 0 
              sim.Individuals[n.per.group*(j-1)+t,(1+tmp.cens[start.int]):tmp.cens[end.int]]<- 1 
            }
          }else{
            sim.Individuals[n.per.group*(j-1)+t,1:tmp.tmax]<- 0
          }       
        }else{
          if(!is.na(tt)){
            sim.Individuals[n.per.group*(j-1)+t,1:(tt-1)]<- 0
            sim.Individuals[n.per.group*(j-1)+t,tt]<- 1
          }else{
            sim.Individuals[n.per.group*(j-1)+t,1:tmp.tmax]<- 0
          }
        }
      }
    }   
  }
  
  sim.Individuals[,1]<- 0 
  allSim<- list(sim.Censored=sim.Individuals, sim.Daily=sim.Daily)
  return(allSim) 
}

# Unnecessary function
covariateArray<- function(path_to_covariate_data, covar_list=c('int','tmax','sm','txsm','light','age')){ #given the path to covariate data for simulation, will construct the array for modeling
  holdArray<- list()
  for(j in 1:length(covar_list)){
    holdArray[[j]]<- as.matrix(read.csv(paste0(path_to_covariate_data, covar_list[j], '_simulation.csv')))
  }
  
  out_array <- array(
    data = do.call(rbind, lapply(holdArray, as.vector)), 
    dim = c(length(holdArray), dim(holdArray[[1]])))
  
  return(out_array)
}

# Function for setting up fitting structures

# Set up for fitting data, and running through without alpha parameter
initializeFit<- function(survivalMatrix, X, iterations=1000, 
                         progressBar = T, X.chronic = c(2,3,4), X.acute=c(1,5,6), 
                         beta_starting=NULL, #can start chain at specific beta values
                         alpha_starting=NULL, # can start chain at specific alpha values, won't be used here
                         Xnames=seq(1, dim(X)[1])){
  run.chronic<- Xnames[X.chronic]
  ng<- iterations
  updatecov<- ng/100
  adapt.mult<-0.1
  adapt.adjust<- (-0.7)
  
  ww <- which(!(survivalMatrix == (-1)), arr.ind=T)
  dL <- matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))
  dL.static <- dL.acc = matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))
  
  IDall <- NULL
  for(i in 1:N){
    IDall <- c(IDall, rep(i, sum(ww[,1]==i)))
  }
  
  survivalMatrix2 <- as.vector(t(survivalMatrix))
  survivalMatrix2 <- survivalMatrix2[survivalMatrix2!=(-1)]
  
  YY<- survivalMatrix2
  
  dL.acc2 <- array(0, dim=c(nrow(survivalMatrix), ncol(survivalMatrix), length(X.chronic)))
  
  dL <- dL.static <- dL.acc <- matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))

  Beta = matrix(0, nrow=ng,ncol=dim(X)[1])
  Beta[1,1]<- -4 # starting value of beta with high survival probability
  if(!is.null(beta_starting))Beta[1, ]<- beta_starting
  colnames(Beta)<- Xnames
  p <- ncol(Beta)
  mu0 <- rep(0, p)
  S0 <- 100 * diag(p)
  nT <- ncol(survivalMatrix)
  likVector <- rep(-Inf, ng)
  rmsurvSim <- which(!(survivalMatrix==(-1)),arr.ind=T)
  reorder <- order(rmsurvSim[,1] + (rmsurvSim[,2] /(max(rmsurvSim[,2])+1)))
  rmsurvmat<- rmsurvSim <- rmsurvSim[reorder,]
  covb =  diag(ncol(Beta))*1E-5

  alphastep <- rep(.001,length(run.chronic))
  alpha <- matrix(0.01,nrow=ng,ncol=length(X.chronic))
  colnames(alpha)<- run.chronic
  cova<- diag(ncol(alpha)) * 0.0001
  dL.acc<- dL.acc2
  
  if(!is.null(alphasim))alpha[1,]<- alphasim
  aliveT<- function(X)return(length(X[X!=(-1)]))
  tRow <- apply(survivalMatrix, 1, aliveT)
  likVector <- rep(-Inf, ng)
  nt<- ncol(survivalMatrix)
  
  now <- Sys.time()
  acceptance <- 0
  if(progressBar)pb=txtProgressBar(min = 0, max = 1, initial = 0, char = "=", width = NA, title, label, style = 3, file = "")
  
  updatecov <- round(seq(ng/100, ng, length.out=100))
  ### first fit, no alpha
  for(n in 1:(ng-1)){
    if(n %in% updatecov){ # Updating the proposition covariance for metropolis can be tricky, below tackles that in multiple ways
      if(n > (ng/10)){
        covb <- as.matrix(updatePropAdaptive(cov(Beta[round(n/2):n,]), acceptRate = length(unique(Beta[round(n/2):n,1])) / n,adjust=(-0.7), n ,adaptiveAdjust=T))
      }else{
        covb <- (2.38^2)/p*signif(cov(Beta[round(n/2):n,]),5)
        if(sum(diag(covb)==0)>0){
          covb<- diag(p)*1E-8
        }
      }
    }
    
    PB<- Beta[n,] + tnorm.mvt(avec=rep(0,length(Beta[n,])),rep(0,length(Beta[n,])),as.matrix(covb)) 
    
    propBeta<- PB
    propAlpha<- alpha[1,]
    
    lpold<-dmvnorm(Beta[n,],mu0,S0,log=T)    # Old Log Prior
    lpnew<-dmvnorm(propBeta,mu0,S0,log=T)    # New Log Prior
    
    for(j in 1:nT){
      dL.static[,j]<-  t(X[X.acute,,j])%*%propBeta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(X[X.chronic[h],,j]) * propBeta[X.chronic][h]
      }
    }
    
    dL<- dL.static
    
    for(h in 1:length(X.chronic)){
      dL<- dL + dL.acc[,,h]	
    }
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
    dL2[,1]<- 1 #### always survival on the first day of observation
    dL2 <-  dL2[rmsurvmat]
    dL2[dL2==0]<- 1E-16 # Needed for algorithm stability in early iterations when still converging
    survLik <-  survLikelihood(Y = YY, p = dL2, ID = IDall)
    newLik <-  lpnew + survLik 
    
    # Accept or reject based on likelihood ratio (log scale)
    if((newLik - likVector[n]) > log(runif(1))){
      Beta[n+1,] <- propBeta
      likVector[n+1] <- newLik
      alpha[n+1,]<- alpha[n,]
    }else{ # Reject
      Beta[n+1,] <- Beta[n,]
      likVector[n+1] <- likVector[n]
      alpha[n+1,]<- alpha[n,]
    }
    setTxtProgressBar(pb, n/ng)
  }
  
  return(list('Beta'=Beta, 'Alpha'=alpha, 'likVector'=likVector, 'covarianceBeta'=covb, 'ng'=ng))
}

chronicSurvivalFit<- function(initialFit, survivalMatrix, X, iterations=1000, progressBar = T, X.chronic = c(2,3,4), X.acute=c(1,5,6), betasim=NULL, alphasim=NULL, Xnames=seq(1, dim(X)[1])){
  run.chronic<- Xnames[X.chronic]
  ng= iterations
  updatecov = ng/100
  adapt.mult=0.1
  adapt.adjust= (-0.7)
  
  ww = which(!(survivalMatrix == (-1)), arr.ind=T)
  dL = matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))
  dL.static = dL.acc = matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))
  
  IDall = NULL
  for(i in 1:N){
    IDall = c(IDall, rep(i, sum(ww[,1]==i)))
  }
  
  survivalMatrix2 = as.vector(t(survivalMatrix))
  survivalMatrix2 = survivalMatrix2[survivalMatrix2!=(-1)]
  
  YY<- survivalMatrix2
  
  dL.acc2 = array(0, dim=c(nrow(survivalMatrix), ncol(survivalMatrix), length(X.chronic)))
  
  dL = dL.static = dL.acc = matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))
  
  Beta = matrix(0, nrow=ng,ncol=dim(X)[1])
  Beta[1,1]<- -4
  if(!is.null(betasim))Beta[1, ]<- betasim
  colnames(Beta)<- Xnames
  p = ncol(Beta)
  mu0 = rep(0, p)
  S0 = 100 * diag(p)
  nT = ncol(survivalMatrix)
  likVector = rep(-Inf, ng)
  rmsurvSim = which(!(survivalMatrix==(-1)),arr.ind=T)
  reorder = order(rmsurvSim[,1] + (rmsurvSim[,2] /(max(rmsurvSim[,2])+1)))
  rmsurvmat<- rmsurvSim <- rmsurvSim[reorder,]
  covb =  diag(ncol(Beta))*1E-5
  
  alphastep = rep(.001,length(run.chronic))
  alpha = matrix(0.01,nrow=ng,ncol=length(X.chronic))
  colnames(alpha)<- run.chronic
  cova<- diag(ncol(alpha)) * 0.0001
  dL.acc<- dL.acc2
  
  if(!is.null(alphasim))alpha[1,]<- alphasim
  aliveT<- function(X)return(length(X[X!=(-1)]))
  tRow = apply(survivalMatrix, 1, aliveT)
  likVector = rep(-Inf, ng)
  nt<- ncol(survivalMatrix)
  
  ng.old<- initialFit$ng
  BetaBurn<- initialFit$Beta[ng.old,]
  likVector.old<- initialFit$likVector
  Beta<- matrix(NA, nrow=ng, ncol=ncol(initialFit$Beta))
  Beta[1,]<- BetaBurn
  colnames(Beta)<- colnames(initialFit$Beta)
  covb<- initialFit$covarianceBeta
  alpha<- matrix(0.0001,nrow=ng,ncol=ncol(initialFit$Alpha))
  cova<- diag(ncol(alpha)) * 0.0001
  colnames(alpha)<- run.chronic
  updatecov = round(seq(ng/1000, ng, length.out=1000))
  likVector<- rep(likVector.old[ng.old], ng)
  
  
  now = Sys.time()
  acceptance = 0
  if(progressBar)pb=txtProgressBar(min = 0, max = 1, initial = 0, char = "=", width = NA, title, label, style = 3, file = "")
  
  updatecov = round(seq(ng/100, ng, length.out=100))
  
  for(n in 1:(ng-1)){
    if(n %in% updatecov){
      if(n > (ng/20)){
        covb = as.matrix(updatePropAdaptive(cov(Beta[round(n/2):n,]), acceptRate = length(unique(Beta[round(n/2):n,1])) / n,adjust=(-0.7), n ,adaptiveAdjust=T))
      }    
    }
    
    PB<- Beta[n,] + tnorm.mvt(avec=rep(0,length(Beta[n,])),rep(0,length(Beta[n,])),as.matrix(covb)) 
    propBeta<- PB
    propAlpha<- alpha[n,]
    
    lpold<-dmvnorm(Beta[n,],mu0,S0,log=T)    # Old Log Prior
    lpnew<-dmvnorm(propBeta,mu0,S0,log=T)    # New Log Prior
    
    for(j in 1:nT){
      dL.static[,j]<-  t(X[X.acute,,j])%*%propBeta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(X[X.chronic[h],,j]) * propBeta[X.chronic][h]
      }
    }
    
    dL<- dL.static
    for(h in 1:length(X.chronic)){
      dL<- dL + cppWeightedRowSums(dL.acc[,,h], tRow, propAlpha[h])
    }
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
    dL2[,1]<- 1 #### always survival on the first day of observation
    dL2 <-  dL2[rmsurvmat]
    dL2[dL2==0]<- 1E-16
    survLik <-  survLikelihood(Y = YY, p = dL2, ID = IDall)
    newLik <-  lpnew + survLik 
    
    if((newLik - likVector[n]) > log(runif(1))){
      Beta[n+1,] = propBeta
      alpha[n+1,] = alpha[n,]
      likVector[n+1] = newLik
    }else{
      Beta[n+1,] = Beta[n,]
      alpha[n+1,] = alpha[n,]
      likVector[n+1] = likVector[n]
    }
    
    #sample and R/A alpha
    if(n %in% updatecov){
      if(n > (ng/20)){
        if(ncol(alpha)==1){
          cova<- as.matrix(updatePropAdaptive(as.matrix(sd(alpha[(round(n/2):n),1])), acceptRate = length(unique(alpha[(round(n/2):n),1])) / n, n,adjust=(-0.7)))
        }else{
          cova<- updatePropAdaptive(cov(alpha[(round(n/2):n),]), acceptRate = length(unique(alpha[(round(n/2):n),1]))/ n, n, adaptiveAdjust=T,adjust=(-0.7))
        }
      }
    }
    if(ncol(alpha)==1){
      PB<- alpha[n,] + tnorm(1,0,as.matrix(cova), lo=-alpha[n,], hi=1-alpha[n,]) 
    }
    if(ncol(alpha)>1){
      PB<- alpha[n,] + tnorm.mvt(avec=rep(0,length(alpha[n,])),rep(0,length(alpha[n,])),as.matrix(cova), lo=-alpha[n,], hi=1-alpha[n,]) 
    }
    propBeta<- Beta[n+1,]
    propAlpha<- PB
    
    lpnew<-lpold<-dmvnorm(Beta[n+1,],mu0,S0,log=T)    # Old Log Prior
    
    for(j in 1:nT){
      dL.static[,j]<-  t(X[X.acute,,j])%*%propBeta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(X[X.chronic[h],,j]) * propBeta[X.chronic][h]
      }
    }
    
    dL<- dL.static
    for(h in 1:length(X.chronic)){
      dL<- dL + cppWeightedRowSums(dL.acc[,,h], tRow, propAlpha[h])
    }
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
    dL2[,1]<- 1 # always survival on the first day of observationl
    dL2[dL2==0]<- 1E-16
    dL2 <-  dL2[rmsurvmat]
    survLik <-  survLikelihood(Y = YY, p = dL2, ID = IDall)
    newLik <-  lpnew + survLik 
    
    if((newLik - likVector[n]) > log(runif(1))){
      alpha[n+1,] = propAlpha
      likVector[n+1] = newLik
    }else{
      alpha[n+1,] = alpha[n,]
      likVector[n+1] = likVector[n]
    }
    setTxtProgressBar(pb, n/ng)
  }
  
  #Calculate DIC value
  meanBeta<- apply(Beta[round(n/2):n,],2,mean)
  meanAlpha<- apply(alpha[round(n/2):n,],2,mean)
  
  for(j in 1:nt){
    dL.static[,j] = t(X[X.acute,,j])%*%meanBeta[X.acute]
    for(h in 1:length(X.chronic)){
      dL.acc[,j,h]<- t(X[X.chronic[h],,j]) * meanBeta[X.chronic][h]
    }
  }
  dL<- dL.static
  for(h in 1:length(X.chronic)){
    dL<- dL + cppWeightedRowSums(dL.acc2[,,h], tRow, meanAlpha[h])
  }
  dL2 =  1 - ( 1-exp(-exp(dL)))  # survival prob
  dL2 = dL2[rmsurvmat]
  
  dic<- mean(-2*likVector[round(n/2):n])
  dBar = -2*survLikelihood(Y = YY, p = dL2, ID = IDall)
  DIC<- dBar  + 2*(dic - dBar)
  pd<- 1/2 * var(-2*likVector[round(n/2):n])
  DIC2<- pd + dBar
  
  tmpFit<- list(Beta=Beta, Alpha=alpha, ng=ng, likVector=likVector, covb=covb, cova=cova, dic1=DIC, dic2=DIC2)
  return(tmpFit)
}

partitionSurvivalEffectsP.Abs<- function(X, Beta, alpha, maxRecensus, npred = 1000, up.to.day = NULL,
                                         X.acute=c(1,5,6), X.chronic=c(2,3,4)){ #based on mortality-weighted effects instead of predicted
  N<- length(maxRecensus)
  nt<- max(maxRecensus)
  chronic.Contribution<- NULL
  chronic.odds<- NULL
  holdDiffs<- NULL
  dL.static<- matrix(0, nrow= dim(X)[2], ncol=dim(X)[3])
  dL.acc<- array(0, dim = c(dim(X)[2], dim(X)[3],ncol(alpha)))
  holdOdds<- rep(NA, npred)
  mysamp<- sample(round(nrow(Beta) / 2) : nrow(Beta), npred)
  for(n in 1:npred){
    betaSamp<- Beta[mysamp[n],]
    alphaSamp<- alpha[mysamp[n],]
    for(j in 1:nt){
      dL.static[,j]<-  t(X[X.acute,,j])%*%c(betaSamp[X.acute])
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- X[X.chronic[h],,j] * betaSamp[X.chronic][h]
      }
    }
    dL<- dL.static
    
    dL.Today<- dL.static
    for(h in 1:length(X.chronic)){
      dL.Today<- dL.Today + dL.acc[,,h] 
    }
    
    
    for(h in 1:length(X.chronic)){
      dL<- dL + cppWeightedRowSums(dL.acc[,,h], maxRecensus, alphaSamp[h])
    }
    dL.All<- dL
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
    dL2.static<- 1 - (1 - exp(dL.static))
    dL2.Today<- 1 - ( 1-exp(-exp(dL.Today)))
    
    p.All<- 1 - dL2
    p.Today<- 1 - dL2.Today
    up.to.day<- min(c(up.to.day, max(maxRecensus) -1))
    if(!is.null(up.to.day)){
      p.All[,(up.to.day+1):ncol(p.All)]<- NA
      p.Today[,(up.to.day+1):ncol(p.Today)]<- NA
    }
    for(i in 1:nrow(p.All)){
      if(maxRecensus[i]<ncol(p.All)){
        p.All[i,(maxRecensus[i]+1):ncol(p.All)]<- NA
        p.Today[i,(maxRecensus[i]+1):ncol(p.Today)]<- NA
      }
    }
    chronic.odds<- (p.All / (1 - p.All)) / (p.Today / (1 - p.Today))
    
    hio<- NULL
    for(i in 1:nrow(chronic.odds)){
      hio[i]<- mean(abs(log(chronic.odds[i,1:maxRecensus[i]])))
      if(maxRecensus[i]<ncol(chronic.odds)){
        chronic.odds[i,(maxRecensus[i]+1):ncol(chronic.odds)]<- NA
      }
    }
    holdOdds[n]<- mean(hio,na.rm=T)
    
    allDiffs<- apply(abs(p.Today - p.All), 1, mean,na.rm=T)
    
    holdDiffs<- rbind(holdDiffs, cbind(n,1:nrow(p.Today),allDiffs) )
  }
  
  holdDiffMeans<- aggregate(holdDiffs[,3],by=list(holdDiffs[,1]),mean)
  
  tmp<- list('hio'=hio,'pAll'=p.All, 'pToday'= p.Today,'averageOddsIncrease' = holdOdds,'medianOddsIncrease' = median(abs(1-log(chronic.odds))),'tmax'=maxRecensus, 'oddsTimeSeries'=abs(log(chronic.odds)), 'oddsTimeSeriesLog'=(log(chronic.odds)), "holdDiffs" = holdDiffs, 'diffMeans' = holdDiffMeans)
  
  hab.prodT<- t(apply((1-tmp$pToday),1,cumprod) )
  hab.prodA<- t(apply(1-tmp$pAll,1, cumprod) )
  
  partition<- list(allDiffs = abs(tmp$pToday - tmp$pAll)[!is.na(abs(tmp$pToday - tmp$pAll))], survivalDiff=apply(abs(diff(hab.prodA - hab.prodT)), 1, max,na.rm=T), 
                   absolutelyMortalityDiffMean=apply(abs(tmp$pToday - tmp$pAll), 1, mean,na.rm=T), absolutelyMortalityDiffSD=sd(tmp$diffMeans[,2]),
                   cumulativeToday = hab.prodT, cumulativeAll = hab.prodA)
  
  
  return(partition)
}

formatData_Survival<- function(survData, censusDates, dt = 'daily', censusGroups = 'year', dateFormat='%m/%d/%Y'){
  survGroup <- which(colnames(survData)==censusGroups)
  censGroup <- which(colnames(censusDates)==censusGroups)
  
  if(length(survGroup) == 0 | length(censGroup) == 0){
    return('Error: Grouping variable not in both dataframes')
  }
  
  #groups.surv = unique(survData[,survGroup]) #extract unique groups in survival data
  groups.surv = survData[,survGroup]
  Ns = length(unique(groups.surv)) #number of unique groups in the survival data
  groups.cens = censusDates[,censGroup] #extract grouping column from the census data
  unique.groups <- unique(groups.cens) #the unique groupings in census data
  unique.groups <- unique.groups[unique.groups %in% groups.surv]
  Nc = length(unique.groups) # number of unique groups in the census data
  same.length = (sum(unique.groups %in% unique(groups.surv)) & sum(unique(groups.surv) %in% unique.groups))  ### test if census data and survival data have same unique groups
  
  index = 1:nrow(survData)
  plotData<- NULL ### will hold the location in the list and original row number
  
  newSurv<- list()
  for(i in 1:Nc){
    surv.tmp <- survData[groups.surv==unique.groups[i],] ### extract the survival data from the group(year) i
    censdates.tmp <-  censusDates[groups.cens == unique.groups[i],] #extract census data for group i
    index.tmp <- index[groups.surv==unique.groups[i]]      
    
    jd.census<- unique(julian(as.Date(censdates.tmp[,'date'], format=dateFormat))) ### extract julian days
    jd.census<- jd.census[order(jd.census)]
    jd.census<- jd.census[!is.na(jd.census)]
    jd.germ <- julian(as.Date(surv.tmp[,'germ_date'], format=dateFormat))
    jd.death <- julian(as.Date(surv.tmp[,'death_date'], format=dateFormat))
    
    jd.deathcensored <- as.vector(surv.tmp[,'death_date']) == "-" ### individuals we don't observe die
    
    rm.ind <- !(jd.germ == jd.death) ### individuals that were observed to die on the same day as germinate are removed
    rm.ind[is.na(rm.ind)] = TRUE
    jd.germ=jd.germ[rm.ind]
    jd.death=jd.death[rm.ind]      
    index.tmp = index.tmp[rm.ind]
    
    for(j in 1:length(jd.germ)){      
      if(jd.deathcensored[j]==F){ ### if censored this code is different
        t = jd.germ[j]:jd.death[j]
        newSurv = append(newSurv, list())
        survMat<- matrix(NA, nrow=length(t), ncol=2)
        survMat[,1] =  t
        survMat[,2] = rep(0, length(t))
        t.death = which(jd.census == jd.death[j])
        
        start.death<- which(survMat[,1] == jd.census[t.death-1]) + 1 #### we observe death one day after the previous census
        end.death <- which(survMat[,1] == jd.census[t.death]) #### census death was observed
        
        survMat[start.death:end.death,2] = 1 ### code death as 1 for an interval
        t.ind <- nrow(survMat)
        newSurv = append(newSurv, list(survMat))
        plotData = rbind(plotData, cbind(index.tmp[j], min(survMat[,1]), max(survMat[,1]), t.ind))
      }else{
        t = jd.germ[j]:max(jd.census)
        survMat = matrix(NA, nrow=length(t), ncol=2)
        survMat[,1] = t ### this individual was observed to survive the whole period
        survMat[,2] = 0
        t.ind <- nrow(survMat)
        newSurv = append(newSurv, list(survMat))
        plotData = rbind(plotData, cbind(index.tmp[j], min(survMat[,1]), max(survMat[,1]), t.ind))
      }
    }
  }
  return(list(plotData, newSurv))
}
