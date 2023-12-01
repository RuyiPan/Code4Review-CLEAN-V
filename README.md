# Supplementary Codes 

## [CELANV.R](https://github.com/RuyiPan/Code4Review-CLEAN-V/blob/main/CLEANV.R)
This file includes both **CLEAN-V** and the three simplied versions:
1. **CLEAN-V without spatial correlation**
2. **CLEAN-V without cluster enhancement**
3. **Massive univariate**

The way to call them has been documented in the file (line 21-24).
To use the functions, please make sure source the file [Clean_support.cpp](https://github.com/RuyiPan/Code4Review-CLEAN-V/blob/main/Clean_support.cpp).

## [Codes4Simulations.R](https://github.com/RuyiPan/Code4Review-CLEAN-V/blob/main/Code4Simulations.R)
This file contains the way to simulate data for power analysis and FWER checking.

## [GanjgahiScore.R](https://github.com/RuyiPan/Code4Review-CLEAN-V/blob/main/GanScore.R)
This file contains the Ganjgahi et al Scoretest method.  We used the Score test that Ganjgahi showed in the their paper [Ganjgahi et al, 2015](https://doi.org/10.1016/j.neuroimage.2015.03.005). The p-value calculation is based on their second permutation method (Null model residual permutation (P2)).


An example of running the scoretest 
```R
N=100
V=1000
Y=matrix(rnorm(N*V,0,1),N,V)
X=cbind(rep(1,nrow(Y)))
K=kronecker(diag(N/2),matrix(c(1,1,1,1),nrow=2,ncol=2))
nP=2000

res <- Scoretest(Y,X,K,nP)


# maxT procedure can be achieved by  the following.
# check if the maxminum of the test statistics is larger than the 0.95th quantile of the permuted maximum test statistics.
max(res$TS) >= quantile(res$maxT,0.95) 
```