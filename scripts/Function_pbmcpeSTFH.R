###############################################################################
###
### Function pbmcpeSTFH() obtains estimators of the mean cross product errors of 
### the EBLUPs obtained at different time instants under a spatio-temporal 
### Fay-Herriot (STFH) model, by a parametric bootstrap method. 
###
### Developed for Contract 7216304 between Complutense University of Madrid and 
### The World Bank.
###
### Authors: Isabel Molina and Eva Romero 
### File name: pbmcpeSTFH.R
### Updated: April 25, 2025
###
###############################################################################

pbmcpeSTFH <-function(formula,nD,nT,vardir,proxmat,nB=100,model="ST",
                      sigma21_start = 0.5 * median(vardir), rho1_start = 0.5, 
                      sigma22_start = 0.5 * median(vardir), rho2_start = 0.5,
                      data)
{
  if (!missing(data)) {
    Xdt <- model.matrix(formula, data)
    names(data)[names(data) == deparse(substitute(vardir))] <- "vardir"
    vardir <- data[,"vardir"]
  } else    Xdt <- as.matrix(model.matrix(formula))
  
  p<-ncol(Xdt)
  
  if (model!="S" && model!="ST")
  {
    print("Error: Model must be S or ST.")
    return (-1)    
  }
  
  # 1: Fit the STFH model and obtain model parameter estimates
  result <- eblupSTFH(formula=formula,D=nD,T=nT,vardir = vardir,proxmat=proxmat,
                      model=model,sigma21_start=sigma21_start,
                      rho1_start=rho1_start,sigma22_start=sigma22_start,
                      rho2_start=rho2_start,data=data)
  
  beta<-matrix(result$fit$estcoef$beta,nr=p,nc=1)
  
  sigma21<-result$fit$estvarcomp$estimate[1]
  rho1<-result$fit$estvarcomp$estimate[2]
  sigma22<-result$fit$estvarcomp$estimate[3]
  
  # 2: Construct covariance matrices and generate bootstrap area and area-time effects
  
  Id <- diag(1,nrow=nD, ncol=nD)
  Omega1rho1 <- solve(crossprod(Id-rho1*proxmat))  
  sigma21Omega1rho1 <- sigma21*Omega1rho1
  M <- nD*nT
  
  if (model=="ST")
  {
    rho2<-result$fit$estvarcomp$estimate[4]
    
    Unomenrho22_05 <- (1-rho2^2)^(-0.5)
    u2dt_b <- matrix(0, nrow=M, ncol=1)
  }
  
  # Initialize MCPE matrix
  
  mcpedt<-matrix(0,nr=nD,nc=nT+nT*(nT-1)/2)
  countfail<-0
  
  b<-1
  while(b<=nB){ ## bootstrap cycle
    
    u1_b   <- matrix(data=mvrnorm(n=1, mu=rep(0,nD),
                                  Sigma=sigma21Omega1rho1), nrow=nD, ncol=1)
    u1dt_b <- matrix(data=rep(u1_b, each=nT),nrow=M,ncol=1)   
    
    if(model=="S"){       
      u2dt_b <- matrix(data=rnorm(M,mean=0,sd=sqrt(sigma22)),nrow=M,ncol=1) 
    }
    if(model=="ST")
    {       
      epsilondt_b<-matrix(data=rnorm(M,mean=0,sd=sqrt(sigma22)),nrow=M, ncol=1)
      i <- 1
      for (d in 1:nD)
      {
        u2dt_b[i] <- Unomenrho22_05*epsilondt_b[i]
        for (t in 2:nT)
        {
          i <- i+1
          u2dt_b[i] <- rho2*u2dt_b[i-1]+epsilondt_b[i]
        }
        i<-i+1
      }
    } # End of if
    
    # 3: Generate bootstrap model errors
    
    edt_b<-rnorm(M,mean=0,sd=sqrt(vardir))
    
    # 4: Calculate bootstrap vectors of true indicators and bootstrap vectors of direct estimates
    
    mudt_b <- as.vector(Xdt%*%beta + u1dt_b + u2dt_b)
    ydt_b  <- mudt_b + edt_b
    
    # 5: Obtain the bootstrap EBLUPs
    
    result.b <- eblupSTFH(ydt_b~Xdt-1,D=nD,T=nT,vardir,proxmat=proxmat,
                          model=model,sigma21_start=sigma21_start,
                          rho1_start=rho1_start,sigma22_start=sigma22_start,
                          rho2_start=rho2_start,data=data)
    
    if (result.b$fit$convergence==FALSE)
    {
      countfail<-countfail+1
      next
    } 
    
    print(b)

    # 6: Parametric bootstrap MSEs and MCPEs
    
    # Prediction errors
    
    dif <- result.b$eblup - mudt_b
    
    # Matrix with squared and crossed product prediction errors
    
    dif.b<-matrix(0,nr=nD,nc=nT+nT*(nT-1)/2)
    
    pos<-nT
    for (t1 in 1:(nT-1)){
      dif.b[,t1]<-dif[seq(from=t1,to=M,by=nT)]^2
      for (t2 in (t1+1):nT){
        pos<-pos+1
        dif.b[,pos]<-dif[seq(from=t1,to=M,by=nT)]*dif[seq(from=t2,to=M,by=nT)]
      }
    }
    
    dif.b[,nT]<-dif[seq(from=nT,to=M,by=nT)]^2
    
    # Matrix with mean crossed product errors 
    
    mcpedt <- mcpedt + dif.b
    
    b <- b+1     
  } # End of bootstrap cycle
  
  mcpedt.b <- mcpedt/nB
  
  # Matrix of EBLUPs for each time instant in columns
  
  eblup<-result$eblup
  eblup.time<-matrix(0,nr=nD,nc=nT)
  
  for (t1 in 1:nT){
    eblup.time[,t1]<-eblup[seq(from=t1,to=M,by=nT)]
  }
  
  mcpe=mcpedt.b[,(nT+1):(nT+nT*(nT-1)/2)]
  
  combinaciones <- combn(nT, 2)
  nombres <- apply(combinaciones, 2, function(par) {
    paste0("(", par[1], ",", par[2], ")")
  })
  if (nT > 2) {colnames(mcpe)<- nombres}     
  
  return (list(eblup=eblup.time,mse=mcpedt.b[,1:nT],mcpe=mcpe,fails=countfail))
}
