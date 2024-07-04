# RobustDC
This R package implements a Huber-type robust distance correlation measure (called "Robust DC" or "RDC") for sure screening in ultrahigh-dimensional heavy-tailed data.

# Installation
To install the `RobustDC` package, you will first need to install `devtools` package and then execute the following code: 
```
devtools::install_github('kehongjie/RobustDC')
```

# Main Functions
Two main functions in this package are `RDC` and `RDC.lasso`. `RDC` is for screening with robust distance correlation, and `RDC.lasso` is for 
fitting the final LASSO regularized regression model after RDC. You can always use the following command to see more details:
```
library(RobustDC)
?RDC
?RDC.lasso
```

# Application of the Real Data
We apply our RDC function to the real genomic dataset from PAAD to get the estimated distance correlation for each gene. It can be executed using the following codes.

```
load("PAAD_gene_prot_matched.RData")
x <- t(paad_tpm_v2)
x <- scale(x,scale=F)
y <- unlist(paad_prot_v2["MAPK1",])
n <- nrow(x)
p <- ncol(x)

#
trunc <- function(Z,t,n){
  f <- function(tau){
    mean(  ( (Z^4 + tau^4 ) - abs( tau^4 - Z^4 )  ) /2  /(tau^4)  ) -  t/n
  }
  tau = uniroot(f, interval=c(1e-20,1e20),extendInt="yes")$root
  return(tau)
}

Zy=as.numeric(dist(y))
tau_y=trunc(Z=Zy,t=3*log(p),n=n)

y2 <- pmin(abs(y),tau_y)*sign(y)
y2=as.matrix(y2)

x2=matrix(nrow=n,ncol=p)
for (i in 1:p) {
  Zx=as.numeric(dist(x[,i]))
  tau_x=trunc(Z=Zx,t=3*log(p),n=n)
  x2[,i] <- pmin(abs(x[,i]),tau_x)*sign(x[,i])
  print(i)
}

rdc<-RDC(x2,y2,c4=0.2,c2=0.2,auto=FALSE)

```