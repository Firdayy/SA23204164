#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of R functions (\code{BH_lmR}) and Cpp functions (\code{BH_lmC}).
#' @examples
#' \dontrun{
#' X = matrix(rnorm(100*50),100)
#' nonzero = sample(50, 10)
#' beta =0.5 * (1:50 %in% nonzero)
#' y = X %*% beta + rnorm(100)
#' tm <- microbenchmark(
#'   tR = BH_lmR(X,y,0.1),
#'   tC = BH_lmC(X,y,0.1)
#' )
#' print(summary(tm)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @useDynLib SA23204164
NULL

#' @title R packages used in my homework.
#' @name homework
#' @description R packages used in my homework.
#' @import DAAG
#' @import boot
#' @import bootstrap
#' @import coda
#' @import dslabs
#' @import ggplot2
#' @import stats
#' @importFrom DT datatable
#' @importFrom htmltools withTags
#' @useDynLib SA23204164
NULL


#' @title FPR and TPR computation using R.
#' @description FPR and TPR computation using R.
#' @param beta the true parameters of the predictor (numeric)
#' @param selected the subscript vector of the selected predictor (numeric)
#' @return a list of false positive rate (FPR) and true positive rate (TPR)
#' @examples
#' \dontrun{
#' X = matrix(rnorm(100*50),100)
#' nonzero = sample(50, 10)
#' beta = 1*(1:50 %in% nonzero)
#' y = X %*% beta + rnorm(100)
#' s<-BH_lmR(X,y,0.1)
#' rR<-unlist(FTPR(beta,sR))
#' }
#' @export
FTPR <- function(beta,selected){
  p=length(beta)
  k=p-sum(beta== 0)
  # If the reject set is empty
  if (length(selected)==1 && selected==0){
    FPR=0
    if(k==0){
      TPR=1
    }else{
      TPR=0
    }
  }else{
    FPR<-sum(beta[selected] == 0) / max(1, length(selected))
    if(k==0){
      TPR=1-length(selected) / p
    }else{
      TPR=sum(beta[selected] != 0) / k
    }
    return(list(FPR=FPR,TPR=TPR))
  } 
}

#' @title BH in linear regression models using R
#' @description BH implementation of variable selection in linear regression models using R
#' @param X the observed value of the predictor (numeric)
#' @param y the observed value of the response (numeric)
#' @param alpha the level of significance of variable selection (hypothesis testing)
#' @return the subscript vector of the selected predictor
#' @examples
#' \dontrun{
#' X = matrix(rnorm(100*50),100)
#' nonzero = sample(50, 10)
#' beta = 1*(1:50 %in% nonzero)
#' y = X %*% beta + rnorm(100)
#' s<-BH_lmR(X,y,0.1)
#' }
#' @export
BH_lmR <- function(X,y,alpha){
  # X: the predictor 
  # y: the response
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # estimate the parameter: beta
  XX <- t(X) %*% X 
  XX_inv <- solve(XX)
  b_hat <- XX_inv %*% t(X) %*% y
  
  # compute the t-statistics
  px <- X %*% XX_inv %*% t(X)
  se <- t(y) %*% (diag(n) - px ) %*% y
  sigma <- as.vector(sqrt(se/(n-p)))*sqrt(diag(XX_inv))
  t <- b_hat/sigma
  
  # compute p-value
  p_value=2*pt(abs(t),df=n-p,lower.tail = FALSE)
  
  #adjust p-value
  dp=p.adjust(p_value,"BH")
  #Returns the set of rejected subscripts
  selected=which(dp <= alpha)
  if (length(selected) ==0) selected=0 
  return(selected)
}
