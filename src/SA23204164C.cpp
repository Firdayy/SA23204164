#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//' @title BH in linear regression models using Rcpp
//' @description BH implementation of variable selection in linear regression models using Rcpp
//' @param X the observed value of the predictor (numeric)
//' @param y the observed value of the response (numeric)
//' @param alpha the level of significance of variable selection (hypothesis testing)
//' @return the subscript vector of the selected predictor
//' @examples
//' \dontrun{
//' X = matrix(rnorm(100*50),100)
//' nonzero = sample(50, 10)
//' beta = 0.5*(1:50 %in% nonzero)
//' y = X %*% beta + rnorm(100)
//' s<-BH_lmC(X,y,0.1)
//' }
//' @export
// [[Rcpp::export]]
arma::uvec BH_lmC(
    arma::mat& X,
    arma::colvec& y,
    double alpha) {
  
  // Dimension information
  int n = X.n_rows, p = X.n_cols;
  
  // Fit model y ~ X
  arma::colvec coef = arma::solve(X, y);
  
  // Compute the residuals
  arma::colvec res = y - X*coef;
  
  // Estimated variance of the random error
  double s2 = std::inner_product(
    res.begin(), res.end(),
    res.begin(), 0.0) / (n - p);
  
  // Standard error matrix of coefficients
  arma::colvec std_err = arma::sqrt(
    s2 * arma::diagvec(arma::pinv(X.t()*X)));
  
  // compute the t-statistics
  arma::colvec t = coef/std_err;
  
  // compute p-value
  arma::vec pvalue(p);
  for ( int j=0; j<p;j++ ){
    if (t(j)<0){t(j)=-t(j);}
    pvalue(j)=2*R::pt(t(j),n-p,0,0);
  }
  // adjust p-value
  
  arma::uvec ind = arma::sort_index(pvalue);
  arma::vec dp = arma::sort(pvalue);
  
  int k = p-1;
  // Returns the set of rejected subscripts
  for ( ;k>-1;k-- ){
    if(dp(k)<=(k+1)*alpha/p) break;
  };
  
  if(k==-1){
    arma::uvec s(1);
    s(0)=0;
    return(s);
  }
  else{
    arma::uvec selected = ind.head(k+1)+1;
    return(selected);
  }
  
}

