# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title BH in linear regression models using Rcpp
#' @description BH implementation of variable selection in linear regression models using Rcpp
#' @param X the observed value of the predictor (numeric)
#' @param y the observed value of the response (numeric)
#' @param alpha the level of significance of variable selection (hypothesis testing)
#' @return the subscript vector of the selected predictor
#' @examples
#' \dontrun{
#' X = matrix(rnorm(100*50),100)
#' nonzero = sample(50, 10)
#' beta = 0.5*(1:50 %in% nonzero)
#' y = X %*% beta + rnorm(100)
#' s<-BH_lmC(X,y,0.1)
#' }
#' @export
BH_lmC <- function(X, y, alpha) {
    .Call('_SA23204164_BH_lmC', PACKAGE = 'SA23204164', X, y, alpha)
}

