##' Robust Distance Correlation for Sure Screening.
##'
##' This function implements a Huber-type robust distance correlation measure
##' (called "Robust DC") for sure screening in ultrahigh-dimensional heavy-tailed
##' data.
##'
##' @param x Predictors, a matrix of dimensions n * p.
##' @param y Response variables, a vector of length n.
##' @param c2 Truncation parameter, which will further be used to solve for the
##' robustification parameter \eqn{\tau \prime}. See the paper for details. Default is 2.
##' @param c4 Truncation parameter, which will further be used to solve for the
##' robustification parameter \eqn{\tau}. See the paper for details. Default is 2.
##' @param auto If TRUE, then the parameter c2 would automatically be tuned
##' based on c4. Default is FALSE.
##'
##' @return A vector of screening statistics for each predictor.
##' @export
##' @examples
##' \dontrun{
##'
##' set.seed(1234)
##' ## parameter setting
##' n <- 50
##' p <- 500
##' s <- 2
##' x <- matrix(NA,n,p)
##' true.ind <- round(seq(1,p,length.out=s))
##' noise.ind <- (1:p)[-true.ind]
##' beta <- rep(0,p)
##' beta[true.ind] <- c(1,-1)
##' ## simulate data
##' for(j in 1:p){
##'   x[,j] <- rlnorm(n)
##' }
##' y <- x%*%beta + rnorm(n)
##' ## Robust DC
##' rdc.stat <- RDC(x,y)
##'
##' }



RDC <- function(x, y, c4=2, c2=2, auto=FALSE) {
  x=as.matrix(x)
  y=as.matrix(y)
  n=dim(x)[1]
  p=dim(x)[2]
  t4=c4*log(p)

  trunc <- function(Z,a,t,n){
    f <- function(tau){
      mean(  ( (Z^a + tau^a ) - abs( tau^a - Z^a )  ) /2  /(tau^a)  ) -  t/n
    }
    tau = uniroot(f, interval=c(1e-20,1e20),extendInt="yes")$root
    return(tau)
  }

  ##basic ingredients for y
  Zy=as.numeric(dist(y))
  tau_y4=trunc(Z=Zy,a=4,t=t4,n=n)

  if (auto) {
    r_y=(mean(Zy^2))^0.5/(mean(Zy^4))^0.25

    cq=seq(from=c4+0.2,by=0.2,length.out=10)

    c2=c4
    tau_y2=trunc(Z=Zy,a=2,t=c2*log(p),n=n)
    k=1

    while (r_y*(t4/n)^0.25*tau_y4<(c2*log(p)/n)^0.5*tau_y2 && k<=10){
      c2=cq[k]
      tau_y2=trunc(Z=Zy,a=2,t=c2*log(p),n=n)
      k=k+1
    }

    tau_y2=max(tau_y2,tau_y4)
  } else {
    tau_y2=trunc(Z=Zy,a=2,t=c2*log(p),n=n)
  }


  y_outer=as.matrix(dist(y))
  y_outer2=pmin(y_outer,tau_y2)
  y_outer4=pmin(y_outer,tau_y4)
  y_outer4_rowsum=apply(y_outer4,1,sum)
  ##

  rdc=c()
  for (i in 1:p){
    ##basic ingredients for x
    Zx=as.numeric(dist(x[,i]))
    tau_x4=trunc(Z=Zx,a=4,t=t4,n=n)

    if (auto) {
      r_x=(mean(Zx^2))^0.5/(mean(Zx^4))^0.25

      cq=seq(from=c4+0.2,by=0.2,length.out=10)

      c2=c4
      tau_x2=trunc(Z=Zx,a=2,t=c2*log(p),n=n)
      k=1

      while (r_x*(t4/n)^0.25*tau_x4<(c2*log(p)/n)^0.5*tau_x2 && k<=10){
        c2=cq[k]
        tau_x2=trunc(Z=Zx,a=2,t=c2*log(p),n=n)
        k=k+1
      }

      tau_x2=max(tau_x2,tau_x4)
    } else {
      tau_x2=trunc(Z=Zx,a=2,t=c2*log(p),n=n)
    }

    x_outer=as.matrix(dist(x[,i]))
    x_outer2=pmin(x_outer,tau_x2)
    x_outer4=pmin(x_outer,tau_x4)
    x_outer4_rowsum=apply(x_outer4,1,sum)
    ##

    dxy=1/(n*(n-1))*sum(x_outer4*y_outer4)+1/(n*(n-1))*sum(x_outer2)*1/(n*(n-1))*sum(y_outer2)-2*1/(n*(n-1)*(n-2))*(sum(x_outer4_rowsum*y_outer4_rowsum)-sum(x_outer4*y_outer4))
    dxx=1/(n*(n-1))*sum(x_outer4*x_outer4)+1/(n*(n-1))*sum(x_outer2)*1/(n*(n-1))*sum(x_outer2)-2*1/(n*(n-1)*(n-2))*(sum(x_outer4_rowsum*x_outer4_rowsum)-sum(x_outer4*x_outer4))
    dyy=1/(n*(n-1))*sum(y_outer4*y_outer4)+1/(n*(n-1))*sum(y_outer2)*1/(n*(n-1))*sum(y_outer2)-2*1/(n*(n-1)*(n-2))*(sum(y_outer4_rowsum*y_outer4_rowsum)-sum(y_outer4*y_outer4))

    rdc[i]=dxy/sqrt(abs(dxx*dyy))

  }

  return(rdc)
}


