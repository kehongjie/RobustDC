// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double S1cpp(NumericVector x, NumericVector y){
  int n = x.size();
  double s1 =0 ;
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      s1 += abs(x(i) - x(j))*abs(y(i) - y(j));
    }
  }
  return(s1 / (n*n));
}

// [[Rcpp::export]]
double S2cpp(NumericVector x, NumericVector y){
  int n = x.size();
  double s2x =0 ;
  double s2y =0 ;
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      s2x += abs(x(i) - x(j));
      s2y += abs(y(i) - y(j));
    }
  }
  return((s2x / (n*n))*(s2y / (n*n)));
}

// [[Rcpp::export]]
double S3cpp(NumericVector x, NumericVector y){
  int n = x.size();
  double s3 = 0 ;
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      for(int k = 0; k < n; ++k){
        s3 += abs(x(i) - x(j))*abs(y(i) - y(k));
      }
    }
  }
  return(s3 / (n*n*n));
}

// [[Rcpp::export]]
arma::cube S3ijk(NumericVector x, NumericVector y){
  int n = x.size();
  arma::cube s3(n,n,n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      for(int k = 0; k < n; ++k){
        s3(i,j,k) = abs(x(i) - x(j))*abs(y(i) - y(k));
      }
    }
  }
  return(s3);
}


// [[Rcpp::export]]
arma::vec U3xycpp(arma::mat comb_mat){
  int n = comb_mat.n_rows;
  arma::vec outvec(n);

  for(int i=0; i<n; i++){
    double x1=comb_mat(i,0);
    double x2=comb_mat(i,1);
    double x3=comb_mat(i,2);
    double y1=comb_mat(i,3);
    double y2=comb_mat(i,4);
    double y3=comb_mat(i,5);
    outvec(i) = (abs(x1-x2)*abs(y1-y3) + abs(x1-x3)*abs(y1-y2) + abs(x2-x1)*abs(y2-y3) +
      abs(x2-x3)*abs(y2-y1) + abs(x3-x1)*abs(y3-y2) + abs(x3-x2)*abs(y3-y1)  ) /6;
  }
  return(outvec);
}

// [[Rcpp::export]]
arma::vec U3xcpp(arma::mat comb_mat){
  int n = comb_mat.n_rows;
  arma::vec outvec(n);

  for(int i=0; i<n; i++){
    double x1=comb_mat(i,0);
    double x2=comb_mat(i,1);
    double x3=comb_mat(i,2);
    outvec(i) = (2*abs(x1-x2)*abs(x1-x3) + 2*abs(x2-x1)*abs(x2-x3) + 2*abs(x3-x1)*abs(x3-x2) ) /6;
  }
  return(outvec);
}

// [[Rcpp::export]]
arma::vec U3ycpp(arma::mat comb_mat){
  int n = comb_mat.n_rows;
  arma::vec outvec(n);

  for(int i=0; i<n; i++){
    double y1=comb_mat(i,0);
    double y2=comb_mat(i,1);
    double y3=comb_mat(i,2);
    outvec(i) = (2*abs(y1-y2)*abs(y1-y3) + 2*abs(y2-y1)*abs(y2-y3) + 2*abs(y3-y1)*abs(y3-y2) ) /6;
  }
  return(outvec);
}

// [[Rcpp::export]]
arma::vec S3xycpp(arma::mat index_mat,
                  arma::mat x_mat,
                  arma::mat y_mat) {
  int N = index_mat.n_rows;
  arma::vec S3xy(N);

  for(int i=0; i<N; i++){
    arma::vec ind_vec = vectorise(index_mat.rows(i,i));
    double ind1 = ind_vec[0]-1;
    double ind2 = ind_vec[1]-1;
    double ind3 = ind_vec[2]-1;
    double x12 = x_mat(ind1,ind2);
    double x13 = x_mat(ind1,ind3);
    double x23 = x_mat(ind2,ind3);
    double y12 = y_mat(ind1,ind2);
    double y13 = y_mat(ind1,ind3);
    double y23 = y_mat(ind2,ind3);

    S3xy[i] = x12*y13 + x12*y23 + x23*y12 + x23*y13 +
      x13*y12 + x13*y23;
  }
  return(S3xy);
}

// [[Rcpp::export]]
arma::vec S3xycpp_v2(arma::mat index_mat,
                     arma::mat x_mat,
                     arma::mat y_mat) {
  int N = index_mat.n_rows;
  arma::vec S3xy(N);

  for(int i=0; i<N; i++){
    arma::vec ind_vec = vectorise(index_mat.rows(i,i));
    double ind1 = ind_vec[0]-1;
    double ind2 = ind_vec[1]-1;
    double x12 = x_mat(ind1,ind2);
    double y12 = y_mat(ind1,ind2);
    S3xy[i] = 2*x12*y12;
  }
  return(S3xy);
}


// [[Rcpp::export]]
arma::vec S3xxcpp(arma::mat index_mat,
                  arma::mat x_mat) {
  int N = index_mat.n_rows;
  arma::vec S3xx(N);

  for(int i=0; i<N; i++){
    arma::vec ind_vec = vectorise(index_mat.rows(i,i));
    double ind1 = ind_vec[0]-1;
    double ind2 = ind_vec[1]-1;
    double ind3 = ind_vec[2]-1;
    double x12 = x_mat(ind1,ind2);
    double x13 = x_mat(ind1,ind3);
    double x23 = x_mat(ind2,ind3);

    S3xx[i] = x12*x13 + x12*x23 + x23*x12 + x23*x13 +
      x13*x12 + x13*x23  ;
  }
  return(S3xx);
}

// [[Rcpp::export]]
arma::vec S3xxcpp_v2(arma::mat index_mat,
                     arma::mat x_mat) {
  int N = index_mat.n_rows;
  arma::vec S3xx(N);

  for(int i=0; i<N; i++){
    arma::vec ind_vec = vectorise(index_mat.rows(i,i));
    double ind1 = ind_vec[0]-1;
    double ind2 = ind_vec[1]-1;
    double x12 = x_mat(ind1,ind2);

    S3xx[i] = 2*x12*x12  ;
  }
  return(S3xx);
}


// [[Rcpp::export]]
arma::vec S3yycpp(arma::mat index_mat,
                  arma::mat y_mat) {
  int N = index_mat.n_rows;
  arma::vec S3yy(N);

  for(int i=0; i<N; i++){
    arma::vec ind_vec = vectorise(index_mat.rows(i,i));
    double ind1 = ind_vec[0]-1;
    double ind2 = ind_vec[1]-1;
    double ind3 = ind_vec[2]-1;
    double y12 = y_mat(ind1,ind2);
    double y13 = y_mat(ind1,ind3);
    double y23 = y_mat(ind2,ind3);

    S3yy[i] = y12*y13 + y12*y23 + y23*y12 + y23*y13 +
      y13*y12 + y13*y23 ;
  }
  return(S3yy);
}


// [[Rcpp::export]]
arma::vec S3yycpp_v2(arma::mat index_mat,
                     arma::mat y_mat) {
  int N = index_mat.n_rows;
  arma::vec S3yy(N);

  for(int i=0; i<N; i++){
    arma::vec ind_vec = vectorise(index_mat.rows(i,i));
    double ind1 = ind_vec[0]-1;
    double ind2 = ind_vec[1]-1;
    double y12 = y_mat(ind1,ind2);
    S3yy[i] = 2*y12*y12 ;
  }
  return(S3yy);
}


