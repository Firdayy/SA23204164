## ----eval=FALSE---------------------------------------------------------------
#  BH_lmR <- function(X,y,alpha){
#    # X: the predictor
#    # y: the response
#  
#    n <- dim(X)[1]
#    p <- dim(X)[2]
#  
#    # estimate the parameter: beta
#    XX <- t(X) %*% X
#    XX_inv <- solve(XX)
#    b_hat <- XX_inv %*% t(X) %*% y
#  
#    # compute the t-statistics
#    px <- X %*% XX_inv %*% t(X)
#    se <- t(y) %*% (diag(n) - px ) %*% y
#    sigma <- as.vector(sqrt(se/(n-p)))*sqrt(diag(XX_inv))
#    t <- b_hat/sigma
#  
#    # compute p-value
#    p_value=2*pt(abs(t),df=n-p,lower.tail = FALSE)
#  
#    #adjust p-value
#    dp=p.adjust(p_value,"BH")
#    #Returns the set of rejected subscripts
#    selected=which(dp <= alpha)
#    if (length(selected) ==0) selected=0
#    return(selected)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  arma::uvec BH_lmC(
#      arma::mat& X,
#      arma::colvec& y,
#      double alpha) {
#  
#    // Dimension information
#    int n = X.n_rows, p = X.n_cols;
#  
#    // Fit model y ~ X
#    arma::colvec coef = arma::solve(X, y);
#  
#    // Compute the residuals
#    arma::colvec res = y - X*coef;
#  
#    // Estimated variance of the random error
#    double s2 = std::inner_product(
#      res.begin(), res.end(),
#      res.begin(), 0.0) / (n - p);
#  
#    // Standard error matrix of coefficients
#    arma::colvec std_err = arma::sqrt(
#      s2 * arma::diagvec(arma::pinv(X.t()*X)));
#  
#    // compute the t-statistics
#    arma::colvec t = coef/std_err;
#  
#    // compute p-value
#    arma::vec pvalue(p);
#    for ( int j=0; j<p;j++ ){
#      if (t(j)<0){t(j)=-t(j);}
#      pvalue(j)=2*R::pt(t(j),n-p,0,0);
#    }
#    // adjust p-value
#  
#    arma::uvec ind = arma::sort_index(pvalue);
#    arma::vec dp = arma::sort(pvalue);
#  
#    int k = p-1;
#    // Returns the set of rejected subscripts
#    for ( ;k>-1;k-- ){
#      if(dp(k)<=(k+1)*alpha/p) break;
#    };
#  
#    if(k==-1){
#      arma::uvec s(1);
#      s(0)=0;
#      return(s);
#    }
#    else{
#      arma::uvec selected = ind.head(k+1)+1;
#      return(selected);
#    }
#  
#  }

## ----eval=FALSE---------------------------------------------------------------
#  FTPR <- function(beta,selected){
#    # the number of all candidate predictors
#    p=length(beta)
#    # the number of nonzero parameters
#    k=p-sum(beta== 0)
#    # If the reject set is empty and k>0,FPR=0 and TPR=0;
#    # If the reject set is empty and k=0,FPR=0 but TPR=1;
#    if (length(selected)==1 && selected==0){
#      FPR=0
#      if(k==0){
#        TPR=1
#      }else{
#        TPR=0
#      }
#    }else{
#      # If the reject set is nonempty
#      FPR<-sum(beta[selected] == 0) / max(1, length(selected))
#      if(k==0){
#        TPR=1-length(selected) / p
#      }else{
#        TPR=sum(beta[selected] != 0) / k
#      }
#      return(list(FPR=FPR,TPR=TPR))
#    }
#  }

## ----eval=TRUE----------------------------------------------------------------
library(SA23204164)
n = 100; # number of observation
p = 50; # number of candidate predictor
k = 10; # number of true predictors
A =0.5; # amplitude
rho = 0.2 # Correlations between predictors
Sigma = toeplitz(rho^(0:(p-1)))

M=1e3 # Number of replicates
rR<-matrix(NA,nrow=2,ncol=M)
rC<-matrix(NA,nrow=2,ncol=M)
sR<-sC<-NULL
for (i in 1:M){
  #Generate simulation samples
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = sample(p, k)
  beta = A * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n)
  
  sR<-BH_lmR(X,y,0.1)
  rR[,i]<-unlist(FTPR(beta,sR))
  
  sC<-BH_lmC(X,y,0.1)
  rC[,i]<-unlist(FTPR(beta,sC))

}

R<-rowMeans(rR)
C<-rowMeans(rC)
m<-cbind(R,C)
rownames(m)<-c("FDR","TDR")
knitr::kable(m)

## ----eval=TRUE----------------------------------------------------------------
library(SA23204164)
library(microbenchmark)

n = 100 # number of observation
p = 50 # number of candidate predictor
k = 10  # number of true predictors 
A =0.5  # amplitude
X = matrix(rnorm(n*p),n) 
nonzero = sample(p, k)
beta = A * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)

tm <- microbenchmark(
  tR = BH_lmR(X,y,0.1),
  tC = BH_lmC(X,y,0.1)
)
knitr::kable(summary(tm)[,c(1,3,5,6)])

