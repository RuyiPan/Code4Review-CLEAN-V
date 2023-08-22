# ymat        : a V times N matrix of imaging data (V: # of vertices, N: # of images)
# mod         : a N times p covariate matrix (p: # of covariates)
#             : It can be generated easily by using the model.matrix() function.
# K           : a N times N matrix specifying between-image dependencies.
# distmat     : a V times V distance matrix
# sacf        : spatial autocorrelation function
#             : The exponential function is assumed as a default.
#             : Other choices include "gau" (Gaussian) 
#             : and "mix" (mixture of exponential and Gaussian)
# max.radius  : The maximum radius for cluster enhancement. 20 is assumed as a default.
# spatial     : Indicate modelling spatial autocorrelation function or not
# nperm       : number of permutations to be used. At least 5000 permutation is recommended.
# alpha       : A desired FWER. alpha=0.05 is assumed as a default.
# seed        : A random seed to be used for permutation.
#             : It is important to use the same seed for integrating results from both hemispheres.
# cores       : The number of cores when parallel computing is executed.
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("/Users/ruyipan/Desktop/Neuroimage/CLEAN-V/Review/Codes4Review/Codes4Post/Clean_support.cpp")

#CLEAN-V: set max.radius to a postive integer, set spatial =T
#CLEAN-V without spatial autocorrelation: set max.radius to a postive intege, set spatial=F
#CLEAN-V without clusterwise enhancement: set max.radius=0, set spatial=T
#Massive Univariate: set max.radius=0, set spatial=F

CleanV=function(ymat,
                distmat, 
                cortex = NULL,
                mod = NULL,
                K,
                sacf = "exp",
                max.radius = 20,
                nperm = 5000, 
                alpha = 0.05,
                seed = NULL, 
                nngp = T,
                nngp.J = 50,
                spatial = F,
                # npartition = NULL, 
                # parallel = F, 
                ncores = 1){
  
  if ( alpha < 0 |alpha > 1){
    stop("[CLEAN] alpha should range between 0 and 1.")
  }
  if ( nrow(distmat) != nrow(ymat) ){
    stop("[CLEAN] The number of rows of distmat and the number of rows of ymat needs to be the same (# vertices).")
  }
  if (is.null(seed)){ 
    seed = sample(1e6, 1) 
  }
  
  V = nrow(ymat)
  if (!is.null(cortex)){
    ymat = ymat[cortex, ]
    distmat = distmat[cortex, cortex]
  } else{
    cortex= 1:V
  }
  
  if (is.null(mod)){
    mod=rep(1, ncol(ymat))
  }
  if (spatial) {
    ymat.leverage = spLeverage(data=ymat, distmat=distmat, mod0=mod, sacf=sacf, nngp=nngp, nngp.J=nngp.J)$out
  } else {
    ymat.leverage = ymat
  }
  
    # s
  NNmatrix = buildNNmatrixDist(distmat, max.radius = max.radius)
  
  K=Matrix(K, sparse=T)

  out = CleanVarC(ymat.leverage, NNmatrix, K, nperm, seed)
  out$seed = seed
  out$nlocations = ncol(NNmatrix)
  out$alternative = "greater"

  set.seed(NULL)
  out$Tstat = apply(matrix(out$Tstat, out$nlocations), 1, max)
  pU <- lapply(c(1:nperm), function(b) apply(matrix(out$permU[,b], out$nlocations), 1, max))
  out$permU <- do.call(cbind, pU)
  return(out)
}


buildNNGPmat=function(distMat, NNGPinfo, params, kernel = "exp"){
  m=nrow(distMat)
  A=matrix(0,m,m)
  D=matrix(0,m,m)
  phi=params$phi
  sigma2=params$sigma2
  tau2=params$tau2
  k=params$K
  
  f.exp=function(phi, d){ exp(-phi*d) }
  f.gau=function(phi, d){ exp(-phi*d^2/2)}
  
  for (i in 1:(nrow(NNGPinfo$NN)-1)){
    nn=na.omit(NNGPinfo$NN[i+1,])
    lnn=length(nn)
    coordip1=nn[lnn]
    
    if (kernel=="exp"){
      K=sigma2*f.exp(phi,distMat[nn,nn,drop=F])+tau2*diag(lnn)
    } else if (kernel=="gau"){
      K=sigma2*f.gau(phi,distMat[nn,nn,drop=F])+tau2*diag(lnn)
    } else if (kernel=="mix"){
      K=sigma2[1]*f.exp(phi[1],distMat[nn,nn,drop=F])+
        sigma2[2]*f.gau(phi[2],distMat[nn,nn,drop=F])+
        tau2*diag(lnn)
    }
    
    A[coordip1,nn[-lnn]]=solve(K[-lnn,-lnn], K[lnn,-lnn])
    D[coordip1,coordip1]=K[lnn,lnn]-sum(K[lnn,-lnn]*solve(K[-lnn,-lnn], K[lnn,-lnn]))
  }
  D[NNGPinfo$NN[1,1],NNGPinfo$NN[1,1]]=sum(sigma2)+tau2
  
  IA=Matrix(diag(ncol(A))-A,sparse=T)
  sD=Matrix(diag(1/diag(D)),sparse=T)
  NNGPprec=IA%*%sD%*%Matrix::t(IA)
  return(list(A=IA, D=sD, NNGPprec=NNGPprec))
}


buildNNmatrixDist=function(distMat, max.radius=20){
  p=nrow(distMat)
  dist_range = sort(c(ceiling(min(distMat, na.rm=T)),ceiling(min(c(max(distMat,na.rm=T),max.radius),na.rm=T))))
  dist_range_sequence = seq(dist_range[1], dist_range[2], by=1) # different radii to consider
  
  out=foreach(i=1:p,.combine="rbind")%do%{
    dt=unname(distMat[i,])
    ind = which(dt<= dist_range[2]) 
    cbind(i,ind, dt[ind])
  }
  
  NNmatrix=foreach(r=1:length(dist_range_sequence), .combine="rbind")%do%{
    dist=dist_range_sequence[r]
    ind2=which(out[,3]<=dist)
    out2=out[ind2,]
    sp=sparseMatrix(i=out2[,1], j=out2[,2], x=1, dims=c(p,p))
    sp
  }
  
  return(NNmatrix)
}

constructNNGPinfo=function(distMat, NN=50, start.vertex=NULL){
  orders=NULL
  m=nrow(distMat)
  out=matrix(NA, m,NN)
  if (is.null(start.vertex)){
    sum.dist=apply(distMat,1,sum)
    prev.vertex=which.min(sum.dist)
  }
  else{prev.vertex=start.vertex}
  
  current.indices=1:m
  out[1,1]=prev.vertex
  orders=c(orders, prev.vertex)
  
  for (j in 2:m){
    current.indices=setdiff(current.indices, prev.vertex)
    vertex=current.indices[which.min(distMat[prev.vertex, current.indices])]
    if (j<=NN){
      out[j,1:j]=c(na.omit(out[j-1,]),vertex) 
    } else{
      distvec=distMat[orders, vertex]
      elements=sort(distvec)[1:(NN-1)]
      ind.elements=(orders[which(distvec%in%elements)])[1:(NN-1)]
      out[j,]=c(ind.elements, vertex)
      # 
      # out[j,]=c(out[j-1, 2:min(j-1,NN)], vertex)
    }
    orders=c(orders, vertex)
    prev.vertex=vertex
  }
  
  return(list(NN=out, orders=orders))
}


CovReg=function(epsilon,  distmat, kernel="exp", n.covariates=NULL, sparse=T, qtl=0.5, maxdist=NULL){
  if (is.null(n.covariates)){    
    n.covariates=0
  }
  if (!is.null(qtl) & !is.null(maxdist)){ stop("Only one of qtl or maxdist should be specified.")}
  if (sparse & is.null(qtl) & is.null (maxdist)){ stop("One of qtl or maxdist should be specified to enable the sparse option.") }
  
  q=NULL
  if (!is.null(qtl)){  q=quantile(distmat, qtl) }
  if (!is.null(maxdist)){ q=maxdist }
  
  n=ncol(epsilon); p=nrow(epsilon)
  if (kernel%in%c("exp", "gau")){
    if (kernel=="exp"){corMat.base=exp(-distmat)}
    else if (kernel=="gau"){corMat.base=exp(-distmat^2/2)}
    
    if (!is.null(q)){ corMat.base=ifelse(distmat<q, corMat.base, 0) }
    
    if (sparse){
      corMat.base=Matrix(corMat.base,sparse=T)
      phi.hat=optimize(CovRegOptim,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarComps(phi.hat, epsilon, corMat.base)
    } else{
      phi.hat=optimize(CovRegOptimC,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarCompsC(phi.hat, epsilon, corMat.base)
    }
  } else if (kernel=="mix"){
    corMat.base1=exp(-distmat)
    corMat.base2=exp(-distmat^2/2)
    
    if (!is.null(q)){
      corMat.base1=ifelse(distmat<q, corMat.base1, 0)
      corMat.base2=ifelse(distmat<q, corMat.base2, 0)
    }
    
    corMat.base1=Matrix(corMat.base1,sparse=T)
    corMat.base2=Matrix(corMat.base2,sparse=T)
    phi.hat=optim(c(0.01, 0.01), CovRegOptim2,
                  epsilon=epsilon, 
                  corMat_base1=corMat.base1,
                  corMat_base2=corMat.base2)$par
    varcomps=ObtainVarComps2(phi.hat, epsilon, corMat.base1, corMat.base2)
  }
  return(list(sigma2=varcomps$sigma2*(n/(n-n.covariates)), 
              tau2=varcomps$tau2*(n/(n-n.covariates)), 
              phi=phi.hat,
              kernel=kernel))
}

spLeverage=function(data, distmat=NULL, mod0=NULL, sacf="exp", nngp=T, nngp.J=50){
  if (!is.null(mod0)){
    data=t(lm(t(data)~mod0)$residuals)
    q=ncol(mod0)
  } else{
    q=0
  }
  
  covreg.fit=CovReg(data, distmat, n.covariates=q, kernel=sacf)
  
  if (nngp & nngp.J<nrow(distmat)){
    NNGPinfo=constructNNGPinfo(distmat, NN=min(nngp.J))
    NNGPprec=buildNNGPmat(distmat, NNGPinfo, covreg.fit, sacf)$NNGPprec
    data.new=as.matrix(NNGPprec%*%data)
  } else{
    
    phi=covreg.fit$phi
    sigma2=covreg.fit$sigma2
    tau2=covreg.fit$tau2
    
    if (sacf=="exp"){
      f.exp=exp(-phi*distmat)
      Sigma=sigma2*f.exp+tau2*diag(nrow(distmat))
    } else if (sacf=="gau"){
      f.gau=exp(-phi*distmat^2/2)
      Sigma=sigma2*f.gau+tau2*diag(nrow(distmat))
    } else if (sacf=="mix"){
      f.exp=exp(-phi*distmat)
      f.gau=exp(-phi*distmat^2/2)
      
      Sigma=sigma2[1]*f.exp+
        sigma2[2]*f.gau+
        tau2*diag(nrow(distmat))
    }
    
    data.new=solve(Sigma, data)
  }
  
  return(list(out = data.new, 
              params.covariance = covreg.fit))
}


CovRegOptim=function(phi, epsilon, corMat_base){
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat=corMat_base^phi
  corMat_norm=sum(corMat^2)
  if (corMat_norm>p+1e-10){
    y1=sum(epsilon*(corMat%*%epsilon))/n
    y2=sum(epsilon^2)/n
    sigma2= (y1-y2)/(corMat_norm-p)
    tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p^2)
    ss= -2*(sigma2*y1*n+tau2*y2*n)+p*tau2*tau2+corMat_norm*sigma2*sigma2+2*sigma2*tau2*p;
  } else{
    sigma2=-1;
    tau2=-1;
    ss=10^10;
  }
  return(ss)
}
ObtainVarComps=function(phi, epsilon, corMat_base){
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat=corMat_base^phi
  corMat_norm=sum(corMat^2)
  
  y1=sum(epsilon*(corMat%*%epsilon))/n
  y2=sum(epsilon^2)/n
  sigma2= (y1-y2)/(corMat_norm-p)
  tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p^2)
  
  return(list(sigma2=sigma2, tau2=tau2))
}