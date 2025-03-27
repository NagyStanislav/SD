// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


#include <Rcpp.h>
using namespace Rcpp;
using namespace std; /* min */

#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]    
 
// [[Rcpp::export]]   
double cl2(arma::vec x, int k){

    int n = x.size();
    arma::mat A(k+1,2);

    A(0,0) = 1;
    A(0,1) = 1;
    for(int i = 0; i < n; i++){
      if(x[i]<0){
        for(int j = 1; j < k+1; j++){
          A(j,0) = A(j,0) + A(j-1,1);
        }
        } else {
        for(int j = 1; j < k+1; j++){
          A(j,1) = A(j,1) + A(j-1,0);
        }
        }
    }
    // Rcout << A << std::endl;
    double res = 0;
    res = A(k,0) + A(k,1);
    return(res);
    } 
 
// [[Rcpp::export]]  
arma::mat SD2C(arma::mat data, int n, int k){
// data (n*2)

  arma::vec sgns(n-1,arma::fill::zeros);
  arma::vec res(n,arma::fill::zeros);
  for(int i = 0; i<n; i++){
  arma::vec datax(n,arma::fill::value(arma::datum::inf));
    for(int j = 0; j<n; j++){
      if(j!=i){
        datax(j) = (data(j,1)-data(i,1))/(data(j,0)-data(i,0));
        }
    }
    // if(i==0) Rcout << datax << std::endl;
    arma::uvec ordr = arma::sort_index(datax, "ascend");
    for(int j = 0; j<n-1; j++){
      sgns(j) = data(ordr(j),0) - data(i,0);
  }
  // if(i==0) Rcout << sgns << std::endl;
  res(i) = cl2(sgns, k);
  }
  return(res);
}

 
 
// [[Rcpp::export]]   
double cl(NumericVector x){

    int k = 3;
    int n = x.size();
    NumericMatrix A(k+1,2);

    A(0,0) = 1;
    A(0,1) = 1;
    for(int i = 0; i < n; i++){
      if(x[i]<0){
        for(int j = 1; j < k+1; j++){
          A(j,0) = A(j,0) + A(j-1,1);
        }
        } else {
        for(int j = 1; j < k+1; j++){
          A(j,1) = A(j,1) + A(j-1,0);
        }
        }
    }
    // Rcout << A << std::endl;
    double res = 0;
    res = A(k,0) + A(k,1);
    return(res);
    }

// [[Rcpp::export]]
NumericVector csample( int m,
                       int size,
                       bool replace,
                       NumericVector prob = NumericVector::create()
                       ) {
  NumericVector v(m);
  for(int i = 0; i<m; i++) v(i) = i;
  NumericVector ret = Rcpp::RcppArmadillo::sample(v, size, replace, prob);
  return ret;
}

// [[Rcpp::export]]
NumericVector fsample( NumericVector v,
                       int size,
                       bool replace,
                       NumericVector prob = NumericVector::create()
                       ) {
  NumericVector ret = Rcpp::RcppArmadillo::sample(v, size, replace, prob);
  return ret;
}

// [[Rcpp::export]]
arma::mat mSD(NumericVector Y, int n, int mu, int ml, int S){
  arma::mat res(S,2);
  NumericVector indsu(mu), indsl(ml);
  for(int s = 0; s < S; s++){
    indsu = csample(n, mu, false);
    std::sort(indsu.begin(), indsu.end());
    res(s,0) = cl(Y[indsu]);
    //
    indsl = fsample(indsu,ml,false);
    std::sort(indsl.begin(), indsl.end());
    res(s,1) = cl(Y[indsl]);        
  }  
  return(res);
}

// [[Rcpp::export]]
arma::mat Amat(double rho){
  arma::mat A(2,2);
  double p = sqrt(1+rho)/2, m = sqrt(1-rho)/2;
  A(0,0) = m+p;
  A(0,1) = p-m;
  A(1,0) = p-m;
  A(1,1) = p+m;
  return(A);
}

// [[Rcpp::export]]
arma::mat gen(int n, double rho){
  arma::mat X(2,n, arma::fill::randn);  
  arma::mat A = Amat(rho);
  return(A*X);
}

// [[Rcpp::export]]
double cmb(int n){
  double res = n*(n-1)*(n-2)/6;
  return(res);
}

// [[Rcpp::export]]
NumericVector stl_nth_element(NumericVector x, int n) {
   NumericVector y = clone(x);
   std::nth_element(y.begin(), y.begin()+n-1, y.end());
   return y;
}

// [[Rcpp::export]]
arma::mat gammasim(int n, arma::mat x, int mu, int ml, int S, int B, double rho){
  int xl = x.n_rows;
  arma::mat gamma(B,xl);
  arma::mat X(2,n);
  arma::uvec beta(n);
  NumericVector Y(n);
  arma::mat tempres(S,2);
  arma::mat res(B,xl);
  double th, tl, tu;
  NumericVector tmp(S);
  int Sh = std::floor(S/2);
  
  double cn = cmb(n), cml = cmb(ml), cmu = cmb(mu);
  
  for(int b = 0; b < B; b++){
    X = gen(n,rho);
    for(int xi = 0; xi < xl; xi++){
      // Rcout << "X: " << X << std::endl;
      // Rcout << "rtio: " << (X.row(1)-x(xi,1))/(X.row(0)-x(xi,0)) << std::endl;
      beta = sort_index((X.row(1)-x(xi,1))/(X.row(0)-x(xi,0)));
      // Rcout << beta << std::endl;
      // Rcout << X.row(0)-x(xi,0) << std::endl;
      for(int ii = 0; ii < n; ii++) Y(ii) = X(0,beta(ii))-x(xi,0);
      // Rcout << Y << std::endl; 
      th = cl(Y);
      tempres = mSD(Y, n, mu, ml, S);
      // medians of absolute deviations
      tmp = abs(tempres.col(0)/cmu - th/cn);
      tu = stl_nth_element(tmp, Sh)[Sh];
      tmp = abs(tempres.col(1)/cml - th/cn);
      tl = stl_nth_element(tmp, Sh)[Sh]; 
      // final estimator of gamma     
      res(b,xi) = (log(tl)-log(tu))/(log(mu)-log(ml));          
     }  
  } 
  return(res);
}