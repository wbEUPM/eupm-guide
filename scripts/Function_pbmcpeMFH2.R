###############################################################################
###
### Function pbmcpeMFH2() obtains estimators of the mean cross product errors of 
### the EBLUPs obtained at different time instants under a multivariate 
### Fay-Herriot (MFH) model with homoscedastic area-time effects (MFH2 model), 
### by a parametric bootstrap method. 
###
### Developed for Contract 7216304 between Complutense University of Madrid and 
### The World Bank.
###
### Authors: Isabel Molina and Eva Romero 
### File name: pbmcpeMFH2.R
### Updated: April 25, 2025
###
###############################################################################

pbmcpeMFH2 <-function(formula,vardir,nB=100,data)
{
  nD<-nrow(data)
  nT<-length(formula)
  M <- nD*nT
  
  Xdt1<-model.matrix(formula[[1]],data)
  p<-ncol(Xdt1)
  
  Xdt <- array(0,c(nD,p,nT))
  Xdt[,,1]<-Xdt1
  for (t in 2:nT){
    Xdt[,,t] <- model.matrix(formula[[t]],data)
  }

  # Fit the STFH model and obtain model parameter estimates
  
  result <- eblupMFH2(formula=formula,vardir=vardir,data=data)
  
  beta<-matrix(result$fit$estcoef[,1],nr=p,nc=nT)
  varu2<-result$fit$refvar
  rho<-result$fit$rho[,1]
  
  # Construct covariance matrices and generate bootstrap area-time effects,
  # bootstrap errors and bootstrap vectors of true indicators
  
  Unomenrho2_05 <- (1-rho^2)^(-0.5)
  
  sigmaedts<-data[,vardir]
  
  Sigmaed <- array(0, c(nT, nT, nD))
  for (d in 1:nD) {
    idx <- nT + 1
    for (t1 in 1:(nT - 1)) {
      Sigmaed[t1, t1, d] <- sigmaedts[d, t1]
      for (t2 in (t1 + 1):nT) {
        Sigmaed[t1, t2, d] <- sigmaedts[d, idx]
        Sigmaed[t2, t1, d] <- sigmaedts[d, idx]
        idx <- idx + 1
      }
    }
    Sigmaed[nT, nT, d] <- sigmaedts[d, nT]
  }
  Sigmaed[,,1]
  # Initialize MCPE matrix
  
  mcpedt<-matrix(0,nr=nD,nc=nT+nT*(nT-1)/2)
  countfail<-0
  
  b<-1
  while(b<=nB){ ## bootstrap cycle
    
    print(b)
    
    adt_b<-rnorm(M,mean=0,sd=sqrt(varu2))
    udt_b<-rep(0,M)
    edt_b<-rep(0,M)
    meandt_b<-rep(0,M)
    i <- 1
    for (d in 1:nD){
      udt_b[i] <- Unomenrho2_05*adt_b[i]
      
      edt_b[i:(i+nT-1)] <- mvrnorm(n=1,mu=rep(0,nT),Sigma=Sigmaed[,,d])
      for (t in (1:T)){
        meandt_b[i+t-1] <- Xdt[d,,t]%*%beta[,t]}
      
      for (t in 2:nT){
        i <- i+1
        udt_b[i] <- rho*udt_b[i-1]+adt_b[i]
      }
      i<-i+1
    }
    
    # Calculate bootstrap vectors of direct estimates
    
    mudt_b <- meandt_b+udt_b
    ydt_b  <- mudt_b + edt_b
    
    # Obtain the bootstrap EBLUPs
    
    ydt.mat<-matrix(0,nr=nD,nc=nT)
    mudt.mat<-matrix(0,nr=nD,nc=nT)
    Xdt.mat2<-NULL
    formula.b<-vector(mode = "list", length = nT)
    
    for (t in 1:nT){
      ydt.mat[,t]<-ydt_b[seq(from=t,to=M,by=nT)]
      mudt.mat[,t]<-mudt_b[seq(from=t,to=M,by=nT)]
      Xdt.mat2<-cbind(Xdt.mat2,Xdt[,-1,t])
      formula.b[[t]]<-as.formula(paste("ydt.mat[,",t,"]~Xdt[,-1,",t,"]",sep=""))
    }
    
    ydt.df<-as.data.frame(ydt.mat)
    Xdt.df<-as.data.frame(Xdt.mat2)
    
    data.b<-cbind(ydt.df,Xdt.df,sigmaedts)
    
    an.error.occured <- FALSE
    tryCatch( { result.b <- eblupMFH2(formula=formula.b,vardir=vardir,data=data.b)}
              , error = function(e) {an.error.occured <<- TRUE})
    
    if (an.error.occured){
      countfail<-countfail+1
      next
    } 
    
    if (result.b$fit$convergence==FALSE){
      countfail<-countfail+1
      next
    } 
    
    # 6: Parametric bootstrap MSEs and MCPEs
    
    # Prediction errors
    
    dif <- result.b$eblup - mudt.mat
    
    # Matrix with squared and crossed product prediction errors
    
    dif.b<-matrix(0,nr=nD,nc=nT+nT*(nT-1)/2)
    
    pos<-nT
    for (t1 in 1:(nT-1)){
      dif.b[,t1]<-dif[,t1]^2
      for (t2 in (t1+1):nT){
        pos<-pos+1
        dif.b[,pos]<-dif[,t1]*dif[,t2]
      }
    }
    dif.b[,nT]<-dif[,nT]^2
    
    # Matrix with mean crossed product errors 
    
    mcpedt <- mcpedt + dif.b
    
    b <- b+1     
  } # End of bootstrap cycle
  
  mcpedt.b <- mcpedt/nB
  
  mcpe=mcpedt.b[,(nT+1):(nT+nT*(nT-1)/2)]
  
  combinaciones <- combn(nT, 2)
  nombres <- apply(combinaciones, 2, function(par) {
    paste0("(", par[1], ",", par[2], ")")
  })
  if (nT > 2) {colnames(mcpe)<- nombres}   
  
  # Matrix of EBLUPs for each time instant in columns
  
  eblup<-result$eblup
  
  return (list(eblup=eblup,mse=mcpedt.b[,1:nT],mcpe=mcpe,fails=countfail))
}