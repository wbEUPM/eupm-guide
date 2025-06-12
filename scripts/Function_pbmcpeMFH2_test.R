pbmcpeMFH2 <- function(formula, vardir, nB = 100, data) {
  library(Matrix)
  library(MASS)
  library(matrixcalc)
  
  nD <- nrow(data)
  nT <- length(formula)
  M <- nD * nT
  
  Xdt1 <- model.matrix(formula[[1]], data)
  p <- ncol(Xdt1)
  
  Xdt <- array(0, c(nD, p, nT))
  Xdt[,,1] <- Xdt1
  for (t in 2:nT) {
    Xdt[,,t] <- model.matrix(formula[[t]], data)
  }
  
  result <- eblupMFH2(formula = formula, vardir = vardir, data = data)
  
  beta <- matrix(result$fit$estcoef[,1], nr = p, nc = nT)
  varu2 <- result$fit$refvar
  rho <- result$fit$rho[,1]
  Unomenrho2_05 <- (1 - rho^2)^(-0.5)
  
  sigmaedts <- data[, vardir]
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
    
    # 🔧 Ensure positive definiteness
    Sigma_tmp <- Sigmaed[,,d]
    if (!is.positive.definite(Sigma_tmp)) {
      Sigmaed[,,d] <- as.matrix(nearPD(Sigma_tmp, keepDiag = TRUE)$mat)
    }
  }
  
  mcpedt <- matrix(0, nr = nD, nc = nT + nT * (nT - 1) / 2)
  countfail <- 0
  b <- 1
  
  while (b <= nB) {
    message(b)
    
    adt_b <- rnorm(M, mean = 0, sd = sqrt(varu2))
    udt_b <- rep(0, M)
    edt_b <- rep(0, M)
    meandt_b <- rep(0, M)
    i <- 1
    for (d in 1:nD) {
      udt_b[i] <- Unomenrho2_05 * adt_b[i]
      edt_b[i:(i + nT - 1)] <- mvrnorm(n = 1, mu = rep(0, nT), Sigma = Sigmaed[,,d])
      for (t in 1:nT) {
        meandt_b[i + t - 1] <- Xdt[d,,t] %*% beta[,t]
      }
      for (t in 2:nT) {
        i <- i + 1
        udt_b[i] <- rho * udt_b[i - 1] + adt_b[i]
      }
      i <- i + 1
    }
    
    mudt_b <- meandt_b + udt_b
    ydt_b <- mudt_b + edt_b
    
    ydt.mat <- matrix(0, nr = nD, nc = nT)
    mudt.mat <- matrix(0, nr = nD, nc = nT)
    Xdt.mat2 <- NULL
    formula.b <- vector(mode = "list", length = nT)
    
    for (t in 1:nT) {
      ydt.mat[, t] <- ydt_b[seq(from = t, to = M, by = nT)]
      mudt.mat[, t] <- mudt_b[seq(from = t, to = M, by = nT)]
      Xdt.mat2 <- cbind(Xdt.mat2, Xdt[, -1, t])
      formula.b[[t]] <- as.formula(paste("ydt.mat[,", t, "] ~ Xdt[, -1, ", t, "]", sep = ""))
    }
    
    ydt.df <- as.data.frame(ydt.mat)
    Xdt.df <- as.data.frame(Xdt.mat2)
    data.b <- cbind(ydt.df, Xdt.df, sigmaedts)
    
    an.error.occured <- FALSE
    tryCatch({
      result.b <- eblupMFH2(formula = formula.b, vardir = vardir, data = data.b)
    }, error = function(e) {
      an.error.occured <<- TRUE
    })
    
    if (an.error.occured || !result.b$fit$convergence) {
      countfail <- countfail + 1
      next
    }
    
    dif <- result.b$eblup - mudt.mat
    dif.b <- matrix(0, nr = nD, nc = nT + nT * (nT - 1) / 2)
    pos <- nT
    for (t1 in 1:(nT - 1)) {
      dif.b[, t1] <- dif[, t1]^2
      for (t2 in (t1 + 1):nT) {
        pos <- pos + 1
        dif.b[, pos] <- dif[, t1] * dif[, t2]
      }
    }
    dif.b[, nT] <- dif[, nT]^2
    
    mcpedt <- mcpedt + dif.b
    b <- b + 1
  }
  
  mcpedt.b <- mcpedt / nB
  mcpe <- mcpedt.b[, (nT + 1):(nT + nT * (nT - 1) / 2)]
  
  combinaciones <- combn(nT, 2)
  nombres <- apply(combinaciones, 2, function(par) {
    paste0("(", par[1], ",", par[2], ")")
  })
  if (nT > 2) colnames(mcpe) <- nombres
  
  eblup <- result$eblup
  
  return(list(eblup = eblup, mse = mcpedt.b[, 1:nT], mcpe = mcpe, fails = countfail))
}
