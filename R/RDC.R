#' @useDynLib RobustDC
#' @importFrom Rcpp sourceCpp
NULL

##' Robust Distance Correlation for Sure Screening.
##'
##' This function implements a Huber-type robust distance correlation measure
##' (called "Robust DC") for sure screening in ultrahigh-dimensional heavy-tailed
##' data.
##'
##' @param x Predictors, a matrix of dimensions n * p.
##' @param y Response variables, a vector of length n.
##' @param c1 Truncation parameter, which will further be used to sovle for the
##' robustification parameter \code{gamma}. See the paper for details. The larger
##' c1, the more truncation. Default is 0.5.
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
##' rdc.stat <- RDC(x,y,c1=0.5)
##'
##' }



RDC <- function(x, y, c1=0.5) {
  n <- nrow(x)
  p <- ncol(x)
  t <- c1*log(p)
  y <- as.numeric(y)
  s1xy <- s2xy <- s3xy <- s1xx <- s2xx <- s3xx <- rep(NA,p)

  index.comb3 <- t(combn(1:n,  3))
  index.comb2 <- t(combn(1:n,  2))

  y_diff <- as.numeric(dist(y))

  y_trun <- trunc(Z=y_diff,Y=y,t,n)
  y_mat_trun <- abs(outer(y_trun,y_trun,'-'))
  y_diff_trun <- as.numeric(y_mat_trun)

  s1yy <- mean(y_diff_trun*y_diff_trun)
  s2yy <- mean(y_diff_trun)*mean(y_diff_trun)
  s3yy <- (sum(S3yycpp(index.comb3, y_mat_trun)) +
             sum(S3yycpp_v2(index.comb2, y_mat_trun))) /n^3

  rdvary <- s1yy + s2yy - 2*s3yy

  rdcov <- rdvarx <- rdc <- rep(NA,p)

  #start <- proc.time()

  for(j in 1:p) {
    # if(j%%100==0) print(j)

    xj = x[,j]

    xj_diff <- as.numeric(dist(xj))

    xj_trun <- trunc(Z=xj_diff,Y=xj,t,n)
    xj_mat_trun <- abs(outer(xj_trun,xj_trun,'-'))
    xj_diff_trun <- as.numeric(xj_mat_trun)

    s1xy[j] = mean(xj_diff_trun*y_diff_trun)
    s2xy[j] = mean(xj_diff_trun)*mean(y_diff_trun)
    s3xy[j] <- (sum(S3xycpp(index.comb3, xj_mat_trun,y_mat_trun)) +
                  sum(S3xycpp_v2(index.comb2, xj_mat_trun, y_mat_trun )))/n^3

    s1xx[j] = mean(xj_diff_trun*xj_diff_trun)
    s2xx[j] = mean(xj_diff_trun)*mean(xj_diff_trun)
    s3xx[j] = (sum(S3xxcpp(index.comb3, xj_mat_trun)) +
                 sum(S3xxcpp_v2(index.comb2, xj_mat_trun)))/n^3

    rdvarx[j] <- s1xx[j] + s2xx[j] - 2*s3xx[j]
    rdcov[j] <- s1xy[j] + s2xy[j] - 2*s3xy[j]
    rdc[j] <- ifelse(rdcov[j]<0,0,sqrt(rdcov[j]/sqrt(rdvarx[j]*rdvary)))
  }

  #end <- proc.time()
  #print(end - start)

  return(rdc)

}


trunc <- function(Z,Y,t,n){
  f <- function(tau){
    mean(  ( (Z^4 + tau^4 ) - abs( tau^4 - Z^4 )  ) /2  /(tau^4)  ) -  t/n
  }
  tau = uniroot(f, interval=c(1e-20,1e20),extendInt="yes")$root
  Y_trunc <- pmin(abs(Y),tau)*sign(Y)
  return(Y_trunc)
}
