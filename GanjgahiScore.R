# Whiten function give the SY (whitened data), SX (transformed covariates) as
# Y is the traits, imaging features N times V
# X is the covariates
# vector is the eigenvectors of K design matrix (S matrix in the paper)
Whitened <- function(Y, X, vector) {
  if (is.null(X)) {
    sy <- vector%*%Y
    sx <- X
  } else {
    sy <- vector%*%Y
    sx <- vector%*%X
  }
  return (list(sy=sy, sx=sx))
}


#Score function return the score statistic
#res - residual after OLS
#Z   - design matrix for the auxiliary model [ones(nS,1) diag(value)]
#      cbind(rep(1,nS),eigenvalue(K))
#b   - (Z^TZ)^{-1} (save the time to compute it)
Score <- function(res, Z, b) {
  nS <- nrow(res)
  Ff <- res^2
  SigmaP <- colMeans(Ff);
  F <- Ff / matrix(rep(SigmaP, nS),nrow=nS,byrow = T) - 1
  

  a <- t(F)%*%Z
  
  c <- a%*%b
  
  tt <- c[,1]*a[,1]
  tt2 <- c[,2]*a[,2]
  
  # calculate T_S statistic image
  Ts <- (tt+tt2)/2
  Score <-1/2*((t(Z[,2])%*%F)/SigmaP)
  Ts[which(Score<0)]=0
  
  return (Ts)
}

# Scoretest return the original and permuted score test staistics
# Y  - N*V matrix contains V features for N images
# X  - design matrix
# K  - the kinship matrix
# nP - number of permutations 
Scoretest <- function(Y, X, K, nP) {
  
  # obtain the design matrix for the auxiliary model
  eK <- eigen(K)
  vec <- eK$vectors
  nS <- nrow(Y)
  Z <- cbind(rep(1,nS),eK$values)
  
  # whiten the data
  whiteres <- Whitened(Y,X,vec)
  sy <- whiteres$sy
  sx <- whiteres$sx
  
  if (is.null(sx)) {
    res <- sy
    cov <- 0
    return ("no X")
  } else {
    beta <- solve(t(sx)%*%sx)%*%t(sx)%*%sy
    hat <- diag(nS)-sx%*%solve(t(sx)%*%sx)%*%t(sx)
    cov <- sx%*%beta
    res <- hat%*%sy
  }
  b <- solve(t(Z)%*%Z)
  TS <- Score(res,Z,b) #orginal test statistics
  maxT <- rep(0, nP)
  IdS <- rep(0,length(TS))
  
  ### store the permuted test statistics
  permStat <- matrix(0, nrow=nP, ncol=length(TS))
  for (p in c(1:nP)) {
    ## We used the Null model residual permutation (the second permutation method P2)
    syP <- res[sample(nS),]+cov
    resP <- hat%*%syP
    TSp <- Score(resP, Z, b)
    IdS <- IdS + (TSp >= TS)
    maxT[p] <- max(TSp)
    permStat[p,] <- TSp
  }
  
  UncorPval <- (IdS + 1)/(nP+1)

  return (list(Upval= UncorPval, TS=TS, maxT=maxT,permStat=permStat))
}








