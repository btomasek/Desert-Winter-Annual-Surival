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
    # X.static<- c(1,5)
    #  X.acc<- c(2,3,4)
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
