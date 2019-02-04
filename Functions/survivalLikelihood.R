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
  #use.dplyr<-arrange(use.dplyr, ID,Y) # order the data
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
predict.Survival<- function(betaSamp, alphaSamp, bigX.pred, tRow=newtRow){ #### alpha.type controls which type to do, {'Different', 'Same', 'One', 'Zero'}
  if(alpha.type=='Different'){ 
    for(j in 1:nt){
      dL.static[,j]<-  t(bigX.pred[X.static,,j])%*%c(betaSamp[X.static])
      for(h in 1:length(X.acc)){
        dL.acc[,j,h]<- bigX.pred[X.acc[h],,j] * betaSamp[X.acc][h]
      }
    }
    dL<- dL.static
    for(h in 1:length(X.acc)){
      dL<- dL + cppRowCumsumNew(dL.acc[,,h], tRow, alphaSamp[h])
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
simulateSurvival<- function(x=NULL,tmax=150,n.per.group = 50, beta=NULL,alpha=NULL,alpha.type='Different',X.acc=c(2,3,4),X.static=c(1,5), census=TRUE,average.census.interval=7){
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
    dL.acc2<-dL.acc <- array(0, dim=c(nrow(dL), ncol(dL), length(X.acc)))
    
    #List Xstatic, List Xacc, NumericVector Bstatic, NumericVector Bacc, int N, int NT, double alpha
    for(j in 1:ncol(dL)){
      dL.static[,j]<-  t(x[X.static,,j])%*%beta[X.static]
      for(h in 1:length(X.acc)){
        dL.acc[,j,h]<- t(x[X.acc[h],,j]) * beta[X.acc][h]
      }
    }
    dL<- dL.static
    for(h in 1:length(X.acc)){
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
initializeFit<- function(survivalMatrix, X, iterations=1000, progressBar = T, X.chronic = c(2,3,4), X.acute=c(1,5,6), betasim=NULL, alphasim=NULL, Xnames=seq(1, dim(X)[1])){
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
  bigX<- newbigX
  
  dL = dL.static = dL.acc = matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))

  Beta = matrix(0, nrow=ng,ncol=dim(bigX)[1])
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
  
  now = Sys.time()
  acceptance = 0
  if(progressBar)pb=txtProgressBar(min = 0, max = 1, initial = 0, char = "=", width = NA, title, label, style = 3, file = "")
  
  ### first fit, no alpha
  updatecov = round(seq(ng/100, ng, length.out=100))
  ### first fit, no alpha
  for(n in 1:(ng-1)){
    if(n %in% updatecov){
      if(n > (ng/10)){
        covb = as.matrix(updatePropAdaptive(cov(Beta[round(n/2):n,]), acceptRate = length(unique(Beta[round(n/2):n,1])) / n,adjust=(-0.7), n ,adaptiveAdjust=T))
      }else{
        covb = (2.38^2)/p*signif(cov(Beta[round(n/2):n,]),5)
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
      dL.static[,j]<-  t(bigX[X.acute,,j])%*%propBeta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(bigX[X.chronic[h],,j]) * propBeta[X.chronic][h]
      }
    }
    
    dL<- dL.static
    
    for(h in 1:length(X.chronic)){
      dL<- dL + dL.acc[,,h]	
    }
    dL2 <-   1 - ( 1-exp(-exp(dL)))  # survival prob
    dL2[,1]<- 1 #### always survival on the first day of observation
    dL2 <-  dL2[rmsurvmat]
    dL2[dL2==0]<- 1E-16
    survLik <-  survLikelihood(Y = YY, p = dL2, ID = IDall)
    newLik <-  lpnew + survLik 
    
    
    
    if((newLik - likVector[n]) > log(runif(1))){
      Beta[n+1,] = propBeta
      likVector[n+1] = newLik
      alpha[n+1,]<- alpha[n,]
    }else{
      Beta[n+1,] = Beta[n,]
      likVector[n+1] = likVector[n]
      alpha[n+1,]<- alpha[n,]
    }
    setTxtProgressBar(pb, n/ng)
  }
  
  return(list('Beta'=Beta, 'Alpha'=alpha, 'likVector'=likVector, 'covarianceBeta'=covb))
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
  bigX<- newbigX
  
  dL = dL.static = dL.acc = matrix(0, nrow=nrow(survivalMatrix), ncol=ncol(survivalMatrix))
  
  Beta = matrix(0, nrow=ng,ncol=dim(bigX)[1])
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
  
  BetaBurn<- initialFit$Beta[ng.old,]
  likVector.old<- initialFit$likVector
  Beta<- matrix(NA, nrow=ng, ncol=ncol(initialFit$Beta))
  Beta[1,]<- BetaBurn#apply(BetaBurn[floor(3*ng.old/4):ng.old,],2,median)
  colnames(Beta)<- colnames(initialFit$Beta)
  covb<- initialFit$covarianceBeta
  alpha<- matrix(0.0001,nrow=ng,ncol=ncol(initialFit$alpha))
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
      dL.static[,j]<-  t(bigX[X.acute,,j])%*%propBeta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(bigX[X.chronic[h],,j]) * propBeta[X.chronic][h]
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
      dL.static[,j]<-  t(bigX[X.acute,,j])%*%propBeta[X.acute]
      for(h in 1:length(X.chronic)){
        dL.acc[,j,h]<- t(bigX[X.chronic[h],,j]) * propBeta[X.chronic][h]
      }
    }
    
    dL<- dL.static
    for(h in 1:length(X.chronic)){
      dL<- dL + cppRowCumsumNew(dL.acc[,,h], tRow, propAlpha[h])
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
    dL.static[,j] = t(bigX[X.acute,,j])%*%meanBeta[X.acute]
    for(h in 1:length(X.chronic)){
      dL.acc[,j,h]<- t(bigX[X.chronic[h],,j]) * meanBeta[X.chronic][h]
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
  
}
