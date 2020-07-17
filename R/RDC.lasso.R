##' Fit the Final Regression Model by Regularized Method after RDC Screening.
##'
##' This function takes the output of the robust distance correlation (RDC)
##' sure screening, and then fits the final regression model using the R
##' packages \pkg{glmnet} for the LASSO regularized loglikelihood for the
##' variables selected by RDC.
##'
##' @param x Predictors, a matrix of dimensions n * p.
##' @param y Response variables, a vector of length n. It should be
##' quantitative.
##' @param rdc.stat Robust distance correlation for each predictors, usually
##' the output of \code{RDC} function.
##' @param nrdc Number of pedictors recuited by RDC.
##' @param lambda Optional user-supplied lambda sequence, as in package
##' \pkg{glmnet}.The default is \code{NULL}, and \code{glmnet} chooses its
##' own sequence.
##' @param nfolds Number of folds for cross-validation, default is 10.
##' @param type.measure Loss to use for cross-validation, as in package
##' \pkg{glmnet}. It can be "deviance", "mse" or "mae" for the qutitative
##' \code{y} here. The default is \code{type.measure="deviance"}.
##' @param which.lambda Which \code{lambda} value from the cross-validation
##' to be used in the final model. The default is "lambda.min", which is the
##' value of \code{lambda} that gives minimum mean cross-validated error.
##' The other choice is "lambda.1se", which gives the most regularized model
##' such that error is within one standard error of the minimum.
##'
##'
##' @return Returns an list with \item{rdc.idx}{The vector of indices
##' selected by only RDC screening step. It is of length \code{nrdc}.}
##' \item{lasso.idx}{The vector of indices selected by the LASSO regularization
##' step after the RDC screening step.}
##' \item{coef.est}{The vector of coefficients of the final model with intercept
##' selected by LASSO. The vector is of length \code{nrdc+1}, with the intercept
##' term being the first one.}
##' \item{fit}{A fitted object of type \code{cv.glmnet} for the final model
##' selected by the RDC+LASSO.}
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
##' ## LASSO
##' lasso.fit <- RDC.lasso(x,y,rdc.stat,nrdc=100)
##' print(lasso.fit$rdc.idx)
##' print(lasso.fit$lasso.idx)
##' }



RDC.lasso <- function(x, y, rdc.stat, nrdc, lambda=NULL, nfolds=10,
                    type.measure="deviance", which.lambda="lambda.min") {
  n <- nrow(x)
  p <- ncol(x)
  y <- as.numeric(y)

  rdc.idx <- sort(order(rdc.stat, decreasing=T)[1:nrdc])
  newx <- x[,rdc.idx]
  fit.cv <- glmnet::cv.glmnet(newx, y, lambda=lambda, nfolds=nfolds,
                             type.measure=type.measure)

  coef.cv <- as.vector(coef(fit.cv, s = which.lambda))
  lasso.idx <- rdc.idx[which(coef.cv[-1] != 0)]
  return(list(rdc.idx=rdc.idx, lasso.idx=lasso.idx, coef.est=coef.cv,
              fit=fit.cv))

}



