## ----heights------------------------------------------------------------------
library(dslabs)
attach(heights)
library(ggplot2)
## Violin plot + box plot
ggplot(heights,aes(x=sex,y=height,fill=sex))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white")+
  theme_classic()

## ----blood pressure, echo=FALSE-----------------------------------------------
blood=data.frame(
  weight=c(76.0,91.5,85.5,82.5,79.0,80.5,74.5,
           79.0,85.0,76.5,82.0,95.0,92.5),
  age=c(50,20,20,30,30,50,60,50,40,55,40,40,20),
  BloodPressure=c(120,141,124,126,117,125,123,125,
             132,123,132,155,147)
)
head(blood,10)
plot(blood)


## ----blood pressure 2, echo=FALSE---------------------------------------------

lm.blood=lm(BloodPressure~weight+age,data=blood)

summary(lm.blood)$coef

## ----plot, echo=FALSE---------------------------------------------------------

plot(lm.blood)

## ----sample-------------------------------------------------------------------
# "x" represents the population 
# "size" represents the sample size. The default value is the population size.
# "prob"represents the distribution of the population. "prob=NULL" means that each case has the same probability 1/length(x)

my.sample<-function(x,size=length(x),prob=NULL){
  n<- size
  u<- runif(n)
  if (is.null(prob)){
    prob=rep(1/length(x),length(x))
  }
  cp <- cumsum(prob)
  r <- x[findInterval(u,cp)+1]
  return(r)
}
#test 1
p=c(.2, .3, .5)
x1 <- my.sample(1:3, size = 1000,prob = p )
ct <- as.vector(table(x1))
ct/sum(ct)/p

#test 2:"prob=null" 
x2 <- my.sample(1:3, size = 1000)
ct <- as.vector(table(x2))
3*ct/sum(ct)

#test 3:size=length(x),prob=NULL
my.sample(letters)


## ----Inverse transform--------------------------------------------------------
n<- 1000
u<- runif(n)
a<-which(u<0.5)
x<-c(log(2*u[a]),-log(2-2*u[-a]))
hist(x, prob = TRUE, main = expression(f(x)==1/2*exp(abs(-x))))
y <-seq(-5, 5, .01)
lines(y,0.5*exp(-abs(y)))


## ----acceptance-rejection-----------------------------------------------------
my.beta<-function(a,b,size){
n <- size;y <- numeric(n)
j<-k<-0;
while (k < n) {
u <- runif(1)
j <- j + 1
x <- runif(1) #random variate from g(.)
if (x^(a-1) * (1-x)^(b-1) > u) {
#we accept x
k <- k + 1
y[k] <- x
}
}
return (y)
}

y <- my.beta(a = 3, b = 2,size=1000)
hist(y, prob = TRUE, main = "Beta(3,2)")
x <- seq(0, 1, .01)
fx <- 12 * x^2 * (1 - x)
lines(x, fx)

## ----rescaled Epanechnikov kernel---------------------------------------------
Ek<-function(size){
n <- size;y <- numeric(n)
u1 <- runif(n,-1,1)
u2 <- runif(n,-1,1)
u3 <- runif(n,-1,1)

for (i in 1:n) {
  if (abs(u3[i])>=abs(u2[i]) && abs(u3[i])>=abs(u1[i])) {
  y[i] <- u2[i]
  }else{
  y[i] <- u3[i]
  }
}

return (y)
}

y <- Ek(size=10000)
hist(y, prob = TRUE, main = "rescaled Epanechnikov kernel")
x <- seq(-1, 1, .01)
fx <- 3/4* (1 - x^2)
lines(x, fx)

## ----Buffon’s niddle experiment-----------------------------------------------
K=100
pihat1=vector(length=K)
pihat2=vector(length=K)
pihat3=vector(length=K)
for (i in 1:K){
  rho=c(0.5,0.8,1)
  d <- 2
  l <- d*rho
  n <- 1e6
  X <- runif(n,0,d/2)
  Y <- runif(n,0,pi/2)
  pihat1[i] <- 2*l[1]/d/mean(l[1]/2*sin(Y)>X)
  pihat2[i] <- 2*l[2]/d/mean(l[2]/2*sin(Y)>X)
  pihat3[i] <- 2*l[3]/d/mean(l[3]/2*sin(Y)>X)
}

sample_variance=c(var(pihat1),var(pihat2),var(pihat3))
sample_variance

## ----antithetic variate approach----------------------------------------------
K=100
thetahat1=vector(length=K)
thetahat2=vector(length=K)
for (i in 1:K){
  m=1e5
  X <- runif(2*m,0,1)
  #simple Monte Carlo method
  thetahat1[i]=mean(exp(X))
  #antithetic variate approach
  thetahat2[i]=(sum(exp(X[1:m]))+sum(exp(rep(1,m)-X[1:m])))/2/m
}

1-var(thetahat2)/var(thetahat1)

## ----importance functions-----------------------------------------------------
x <- seq(1.01, 10, .01)
g <- function(x) {
exp(-x^2/2) * (x^2/sqrt(2*pi))*(x > 1)
}
y <- g(x)

par(mfrow=c(1,2))

plot(x, y,xlim = c(1,10), ylim = c(0, 1),type="l",lwd=2)
lines(x, 2*dnorm(x,1), lty = 2,lwd=3,col="red")
lines(x, dgamma(x-1, 2, 2), lty = 3,lwd=3,col="blue")
legend("topright", inset = 0.02, legend = c("g(x)", "f1",
 "f2"),col=c("black","red","blue"), lty = 1:3)

plot(x, y/(2*dnorm(x,1)),type="l",lty = 2,lwd=2, ylab = "",col="red")
lines(x, y/(dgamma(x-1, 2, 2)), lty = 3,lwd=2,col="blue")
legend("topright", inset = 0.02, legend = c("g/f1", "g/f2"),col=c("red","blue"),lty = 2:3)

## ----importance sampling------------------------------------------------------
m <- 100000
n <-100
theta.hat1 <- se1 <- numeric(n)
theta.hat2 <- se2 <- numeric(n)
g <- function(x) {
exp(-x^2/2) * (x^2/sqrt(2*pi))*(x > 1)
}

for (i in 1:n) {
  x <- abs(rnorm(m)) + 1 #using f1
  fg <- g(x)/(2*dnorm(x, 1))
  theta.hat1[i] <- mean(fg)
  
  x <- rgamma(m, 2,2) #using f2
  fg <- g(x+1) /dgamma(x, 2, 2)
  theta.hat2[i] <- mean(fg)
}

c(mean(theta.hat1),mean(theta.hat2))
c(var(theta.hat1),var(theta.hat2))

## ----stratified importance sampling-------------------------------------------
M <- 10000; 
k <- 5 
r <- M/k #replicates per stratum
T2<-var2<- numeric(k)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)

u <- runif(r) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
est1 <- mean(fg)#importance sampling
var1<- var(fg)

for(j in 1:k){
    u <- runif(r,(j-1)/k,j/k)
    x <- - log(1 - u * (1 - exp(-1)))
    fg <- g(x) / (k*exp(-x) / (1 - exp(-1)))
    T2[j]<-mean(fg)
    var2[j]<-var(fg)
}
est2<-sum(T2)#stratified importance sampling
var22<-mean(var2)

rbind(c(est1,est2),c(var1,var22))

## -----------------------------------------------------------------------------
m <- 1e5; n <- 20
c1 <- c2 <-UCL<- numeric(m)
alpha <- .05
t <- qt(alpha/2 , df = n-1)
for(i in 1:m){
  x <- rchisq(n,df=2)
  s <- sd(x)
  c1[i] <- mean(x)+t*s/sqrt(n)
  c2[i] <- mean(x)-t*s/sqrt(n)
  UCL[i] <- (n-1) * var(x) / qchisq(alpha, df=n-1)
}
mean(c1 <= 2 & c2 >= 2)
mean(UCL>4)

## -----------------------------------------------------------------------------
mu<- 1 # null hypothesis
m <- 1e4; n <- 10; set.seed(123)
t1 <- t2<-t3 <- numeric(m)
alpha <- 0.05
t <- qt(alpha/2 , df = n-1)
for(i in 1:m){
  x1 <- rchisq(n,df=1)
  x2 <- runif(n,0,2)
  x3 <- rexp(n,rate=1)
  xbar1=mean(x1);xbar2=mean(x2);xbar3=mean(x3);
  s1=sd(x1);s2=sd(x2);s3=sd(x3);
  t1[i]=(xbar1-mu)*sqrt(n-1)/s1
  t2[i]=(xbar2-mu)*sqrt(n-1)/s2
  t3[i]=(xbar3-mu)*sqrt(n-1)/s3
}
q1=2*(1-pt(abs(t1),n-1)) #real p values
q2=2*(1-pt(abs(t2),n-1))
q3=2*(1-pt(abs(t3),n-1))

print(c(mean(q1<=alpha),mean(q2<=alpha),mean(q3<=alpha)))

## ----Mulitple test------------------------------------------------------------
m<- 1e3    # the total number of multiple tests
pho<- 0.95 # the percentage of null hypothesis
m0<- m*pho
alpha<- 0.1
M<- 1e3 ## replicate M times
p<- FWER1<- FDR1<-TPR1<-numeric(m)
FWER2<- FDR2<-TPR2<-numeric(m)
for (i in 1:M){
  p[1:m0]<-runif(m0)
  p[(m0+1):m]<- rbeta(m-m0,0.1,1)
  pbonf<-p.adjust(p,method ='bonferroni')
  pfdr<-p.adjust(p,method ='BH')
  reject1<- c(which(pbonf<=alpha))
  FWER1[i]<- (reject1[1]<= m0)
  FDR1[i]<- sum(reject1 <= m0)/max(1,length(reject1))
  TPR1[i]<- sum(reject1>m0)/(m-m0) 
  
  reject2<- c(which(pfdr<=alpha))
  FWER2[i]<- (reject2[1]<= m0)
  FDR2[i]<- sum(reject2 <= m0)/max(1,length(reject2))
  TPR2[i]<- sum(reject2>m0)/(m-m0) 
}
FWER<-round(c(mean(FWER1),mean(FWER2)),3)
FDR<-round(c(mean(FDR1),mean(FDR2)),3)
TPR<-round(c(mean(TPR1),mean(TPR2)),3)

res<- data.frame(FWER,FDR,TPR,row.names=c('Bonferroni','BH'))
knitr::kable(res)


## ----bootstrap in exponential distribution------------------------------------
lambda<-2
n<- c(5,10,20)
B=1000 #bootstrap repeat times
lambdastar<-numeric(B)
bias1<-se1<-numeric(m)
bias<-se<-matrix(0,length(n),2)
m=1000 #simulation repeat times
for(j in 1:length(n)){
  for(i in 1:m){
    x=rexp(n[j],rate=lambda)
    lambdahat=1/mean(x)
    for(b in 1:B){
      xstar=sample(x,replace=TRUE)
      lambdastar[b]=1/mean(xstar)
    }
    bias1[i]=mean(lambdastar)-lambdahat
    se1[i]=sd(lambdastar)
  }
  bias[j,2]=lambda/(n[j]-1)
  bias[j,1]=mean(bias1)
  se[j,2]=lambda*n[j]/(n[j]-1)/sqrt(n[j]-2)
  se[j,1]=mean(se1)
}
res<- data.frame(n,bias,se)
sketch = htmltools::withTags(table(
  class='display',
  thead(
    tr(
      th(rowspan = 2, 'n'),
      th(rowspan = 1, colspan = 2, 'bias'),
      th(rowspan = 1, colspan = 2, 'standard error')
    ),
    tr(
      th(rowspan = 1, colspan = 1, 'bootstrap'),
      th(rowspan = 1, colspan = 1, 'theoretical'),
      th(rowspan = 1, colspan = 1, 'bootstrap'),
      th(rowspan = 1, colspan = 1, 'theoretical')
    )
    )
  )
)
DT::datatable(
  res,
  rownames = FALSE,
  escape = FALSE, 
  container = sketch, 
  options = list(scrollY = FALSE,
                 autoWidth = TRUE)
)


## ----bootstrap in CI----------------------------------------------------------
library(bootstrap) #for the law data
#set up the bootstrap
B <- 200 #number of replicates
n <- nrow(law) #sample size
R<-R1<-R2<-numeric(B) #storage for replicates
alpha=0.05
Rhat<-cor(law$LSAT, law$GPA)
#bootstrap estimate of standard error of R
set.seed(12345)
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  LSAT <- law$LSAT[i] #i is a vector of indices
  GPA <- law$GPA[i]
  R[b] <- cor(LSAT, GPA)
  
  # bootstrap estimate of the standard error of R_b
  for(j in 1:B){
    k <- sample(i, size = n, replace = TRUE)
    LSAT <- law$LSAT[k]
    GPA <- law$GPA[k]
    R2[j] <- cor(LSAT, GPA)
  }
  #t distribution sample
  R1[b]=(R[b]-Rhat)/sd(R2)
}
se<- sd(R)
t=unname(quantile(R1,c(1-alpha/2,alpha/2)))
Rhat-t*se


## ----bootstrap CI-------------------------------------------------------------
library(boot)
failuret<-aircondit[1]
set.seed(12345)
boot.mean <- function(x,i) return(mean(as.matrix(x[i,])))
de <- boot(failuret,statistic=boot.mean, R = 999)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
de
ci
plot(de)

## ----jackknife----------------------------------------------------------------
library(bootstrap)
attach(scor)
x<-as.matrix(scor)
n<-nrow(x)
Sigma<-cov(x)
lambda<-eigen(Sigma)$values
theta.hat <- lambda[1]/sum(lambda)
theta.jack <- numeric(n)
for(i in 1:n){
  sigma<-cov(x[-i,])
  lambda<-eigen(sigma)$values
  theta.jack[i] <- lambda[1]/sum(lambda)
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),3)

## ----leave-two-out------------------------------------------------------------
library(DAAG); attach(ironslag)
a <- seq(10, 40, .1)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)
# for leave-two-out cross validation
# fit models on leave-two-out samples
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    k<-c(i,j)
    y <- magnetic[-k]
    x <- chemical[-k]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] *chemical[k]
    e1[i,j] <- sum((magnetic[k] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2]*chemical[k] +J2$coef[3] * chemical[k]^2
    e2[i,j] <- sum((magnetic[k] - yhat2)^2)
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    e3[i,j] <- sum((magnetic[k] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4 <- exp(logyhat4)
    e4[i,j] <- sum((magnetic[k] - yhat4)^2)
  }
}
c(sum(sum(e1)),sum(sum(e2)),sum(sum(e3)),sum(sum(e4)))/(n*(n-1)/2)


## ----Cramer-von Mises---------------------------------------------------------
library(boot)
attach(chickwts)
x1 <- sort(as.vector(weight[feed == "soybean"]))
x2 <- sort(as.vector(weight[feed == "linseed"]))
x3<-sort(as.vector(weight[feed == "sunflower"]))
detach(chickwts)
CM<-function(x,y,R){
  z <- c(x, y) #pooled sample
  n<-length(x)
  m<-length(y)
  Mn <- function(z, sizes) {
    n <- sizes[1]; m <- sizes[2]
    x=z[1:n]
    y=z[(n+1):(n+m)]
    Fn<-numeric(n+m)
    Gm<-numeric(n+m)
    for (i in 1:(m+n)){
      Fn[i]=sum(x<=z[i])/n
      Gm[i]=sum(y<=z[i])/m
    }
    m*n/(m+n)^2*sum((Fn-Gm)^2)
  }
  N<-c(n,m)
  t<-numeric(R)
  for (j in 1:R){
  t0=Mn(z,N)
  k=sample(1:(n+m))
  t[j]=Mn(z[k],N)
  }
  t=c(t0,t)
  p.value <- mean(t>=t0)
  return(c(p.value,t0))
}
set.seed(12345)
t1=CM(x1,x2,R=999)#linseed vs. soybean
t2=CM(x1,x3,R=999)#soybean vs. sunflower
t3=CM(x2,x3,R=999)#linseed vs. sunflower

res<- data.frame(t1,t2,t3,row.names=c('p-value','statistic'))
colnames(res)<-c("linseed vs. soybean","soybean vs. sunflower","linseed vs. sunflower")
knitr::kable(res)

## -----------------------------------------------------------------------------
maxi<- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}

n1 <- 20 
n2 <- 40
n=n1+n2
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1 #NULL
sigma3<-5#alternative
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

#NULL
z=c(x,y)
d0=maxi(x, y)
set.seed(1123)
R<-999
d=numeric(R)
for(i in 1:R){
  k=sample(1:n)
  k1=k[1:n1]
  k2=k[(n1+1):n]
  d[i]=maxi(z[k1], z[k2])
}
d=c(d0,d)
mean(d>=d0)

#alternative
x3 <- rnorm(n1, mu1, sigma3)
d10=maxi(x3,y)
d1=numeric(R)
z1=c(x3,y)
for(i in 1:R){
  k=sample(1:n)
  k1=k[1:n1]
  k2=k[(n1+1):n]
  d1[i]=maxi(z1[k1], z1[k2])
}
d1=c(d10,d1)
mean(d1>=d10)

## -----------------------------------------------------------------------------
solve_a<-function(N,b1,b2,b3,f0){
  x1<-rpois(N,1)
  x2<-rexp(N,1)
  x3<-sample(0:1,size=N,replace=TRUE,prob=c(0.5,0.5))
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2); p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-20,20))
  return(round(unlist(solution),5)[1])
}

f0<-c( 0.1, 0.01, 0.001, 0.0001)
N=1e6;b1 = 0; b2 = 1; b3 =-1
set.seed(12345)
a=numeric(length(f0))
for (i in 1:length(f0)) {
  a[i]=solve_a(N,b1,b2,b3,f0[i])
}
a
lf<- -log(f0)
plot(lf,a,type = "o",xlab = "-log(f_0)")

## ----rwMc---------------------------------------------------------------------
f<-function(x){exp(-abs(x))/2} # aim pdf
set.seed(12345)
rw.Metropolis <- function(n, sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0 #initial value
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma) # proposal distribution
    if (u[i] <= (f(y) / f(x[i-1])))
      x[i] <- y 
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
return(list(x=x, k=k))
}

N <- 2000
sigma <- c(.05, .5, 2.5, 16)
x0 <- 25
rw1 <- rw.Metropolis(n, sigma[1], x0, N)
rw2 <- rw.Metropolis(n, sigma[2], x0, N)
rw3 <- rw.Metropolis(n, sigma[3], x0, N)
rw4 <- rw.Metropolis(n, sigma[4], x0, N)

#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k)/N)


# trace plot

qe<-qexp(0.9875)
# 0.025 and 0.975 quantile of laplace distribution
refline<-c(-qe,qe)
rws<-cbind(rw1$x,rw2$x,rw3$x,rw4$x)
for (j in 1:4) {
  plot(rws[,j],type="l",xlab=bquote(sigma==.(round(sigma[j],3))),
       ylab="X",ylim=range(rws[,j]))
  abline(h=refline)
}


# histogram 

b <- seq(.025,.975,.01)
y <- qexp(b, 1)
x <- c(-rev(y), y)
fx<-f(x)
for (j in 1:4) {
  hist(rws[501:N,j],breaks="scott",prob=TRUE,ylab ="density",xlab=bquote(sigma==.(round(sigma[j],3))),main = "Histogram",ylim = c(0,0.8))
  lines(x,fx)
}

#QQ plot

for (j in 1:4) {
  qqplot(x,quantile(rws[501:N,j],b),xlab = " laplace quantiles",ylab =bquote(sigma==.(round(sigma[j],3))),main="sample quantiles v.s. laplace quantiles",cex=0.4)
   abline(0, 1)
}


## -----------------------------------------------------------------------------
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- .9 #correlation
mu<- 0;sigma <- 1
s <- sqrt(1-rho^2)*sigma

###### generate the chain #####
X[1, ] <- c(mu, mu) #initialize
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu + rho * (x2 - mu) 
  X[i, 1] <- rnorm(1, m1, s)
  x1 <- X[i, 1]
  m2 <- mu + rho * (x1 - mu) 
  X[i, 2] <- rnorm(1, m2, s)
}

b <- burn + 1
x <- X[b:N, ]

# compare sample statistics to parameters

colMeans(x)
cov(x)
cor(x)

plot(x, main="", cex=.5, xlab=bquote(X[1]),
ylab=bquote(X[2]), ylim=range(x[,2]))


# linear regression model Y = β_0 + β_1X
X<-x[,1]
Y<-x[,2]
lm<-lm(Y~X)
summary(lm)
res<-lm$residuals
mean(res)
# qqplot of residuals: normality
qqnorm(res)
qqline(res)
ks.test(res,"pnorm")
#residuals: constant variance
plot(lm$fitted.values,res)
abline(h=0)


## -----------------------------------------------------------------------------

Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

# Rayleigh(σ) density
f <- function(x, sigma) {
  if (any(x < 0)) return (0)
    stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

rl.chain <- function(sigma, N, x0) {
#generates a Metropolis chain for Rayleigh(σ)
#with chi_square(df=X[t]) proposal distribution
# initial value x0
  x <- numeric(N)
  x[1] <-x0
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i]<=num/den)
      x[i] <- y
    else x[i] <- xt
  }
  return(x)
}

sigma<- 4
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 2000 #burn-in length

#choose overdispersed initial values
x0 <- c(5, 10, 20, 30)

set.seed(1234)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k) X[i,]<- rl.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains

for (i in 1:k)
plot(psi[i, (b+1):n], type="l",
xlab=i, ylab=bquote(psi))


#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)


# coda[212]
library(coda)
x1 <- as.mcmc(X[1, ])
x2 <- as.mcmc(X[2, ])
x3 <- as.mcmc(X[3, ])
x4 <- as.mcmc(X[4, ])
y <- mcmc.list(x1, x2, x3, x4)
print(gelman.diag(y))
gelman.plot(y, col = c(1, 1))

## ----EM-----------------------------------------------------------------------
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
n<-length(u)
 # the first order derivative of log-likelihood
dl<-function(lambda){
    d<-0
    for (i in 1:n) {
      t1<-exp(lambda*u[i])
      t2<-exp(lambda*v[i])
      d=d+(v[i]*t1-u[i]*t2)/(t2-t1)
    }
    return(d)
}
 # the second order derivative of log-likelihood
ddl<-function(lambda){
    d<-0
    for (i in 1:n) {
      t1<-exp(lambda*u[i])
      t2<-exp(lambda*v[i])
      d=d-t1*t2*(u[i]-v[i])^2/(t1-t2)^2
    }
    return(d)
}


#Newton-Raphsen algorithm

NR<-function(lambda_0,epsilon=1e-5,maxite=1e5){
  lambda<-numeric(2)
  lambda[1]<-lambda_0
  k<-0
  while(isTRUE(abs(lambda[2]-lambda[1]) >= epsilon) && k<=maxite){
    lambda[2]<-lambda[1]-dl(lambda[1])/ddl(lambda[1])
    k<-k+1
    lambda=rev(lambda)
  }
  return(list(MLE=lambda[1],ite=k,eps=abs(lambda[2]-lambda[1])))
}

#EM
EM<-function(lambda_0,epsilon=1e-5,maxite=1e5){
  lambda<-numeric(2)
  lambda[1]<-lambda_0
  k<-0
  while (abs(lambda[2]-lambda[1])>=epsilon && k<=maxite ) {
    lambda[2]<-(n*lambda[1])/(n-lambda[1]*dl(lambda[1]))
    k<-k+1
    lambda=rev(lambda)
  }
  return(list(MLE=lambda[1],ite=k,eps=abs(lambda[2]-lambda[1])))
}

s1<-NR(0.1)
s2<-EM(0.1)
s<-cbind(as.matrix(s1),as.matrix(s2))
colnames(s)=c('Netwon-Raphsen','EM')
rownames(s)<-c("MLE","iterations times","epsilon")
knitr::kable(s)


## ----Morra game---------------------------------------------------------------
#enter the payoff matrix
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)

B<- A+2

solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}

library(boot) #needed for simplex function

sA <- solve.game(A) # solution of the original game A
sB <- solve.game(B) # solution of the game B

cbind(sA$v,sB$v) #value of game A and B

round(cbind(sA$x, sA$y), 7) #strategy of game A

round(cbind(sB$x, sB$y), 7) #strategy of game B

## -----------------------------------------------------------------------------
# atomic
x <- 1:3
y <- c("x","y")

# list
l <- list(x,y)
is.vector(l)
typeof(as.vector(l))
typeof(unlist(l))
unlist(l)

## -----------------------------------------------------------------------------
#vector
x<-1:6
dim(x)
#Use dim() to convert a vector to a matrix
dim(x)<-c(2,3)
x

## -----------------------------------------------------------------------------
#x is a matrix
x<-matrix(1:6,nrow = 2)
is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
d<-data.frame(a=1:3,b=c("x","y","z"),c=list(1:3))
as.matrix(d)

d2<-data.frame(a=1:3,b=c("x","y","z"),c=I(list(1:2,1:3,1:4)))
as.matrix(d2)

## -----------------------------------------------------------------------------
x<-y<-NULL
d<-data.frame(x,y)
d

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
d<-data.frame(a=1:3,b=c(7,5,9))
outd<-sapply(d,function(x) scale01(x))
outd
# we can convert the output to a data frame
outd<-data.frame(outd)
outd

## -----------------------------------------------------------------------------
d<-data.frame(a=1:3,b=c(TRUE,TRUE,FALSE),c=c("x","y","z"))
outd<-sapply(d,function(x) if (is.numeric(x)) scale01(x) else x)
outd
# we can convert the output to a data frame
outd<-data.frame(outd)
outd

## -----------------------------------------------------------------------------
d<-data.frame(a=1:3,b=c(7,9,18))
outd<-vapply(d,function(x) sd(x),FUN.VALUE=0)
outd

## -----------------------------------------------------------------------------
d<-data.frame(a=1:3,b=c(7,9,18),c=c("x","y","z"))
d1<-d[vapply(d, is.numeric, FUN.VALUE = TRUE)]
outd<-vapply(d1,function(x) sd(x),FUN.VALUE=0)
outd

## -----------------------------------------------------------------------------
#R function
f<-function(n,a,b,N){
  #n,a,b are parameters of the target distribution
  #N is the number of samples that the Gibbs sampler will generate
  X <- matrix(0, N, 2) #the chain, a bivariate sample

  ###### generate the chain #####
  X[1,1] <- rbinom(1,n,0.5)#initialize
  X[1,2] <- rbeta(1,X[1,1]+a,n-X[1,1]+b)
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    X[i, 1] <- rbinom(1,n,x2)
    x1 <- X[i, 1]
    X[i, 2] <- rbeta(1,x1+a,n-x1+b)
  }
  return(X)
}

library(Rcpp) # Attach R package "Rcpp"
# Define function "cf"
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cf( double n, double a, double b, int N) {
  NumericMatrix X(N,2);
  
  X(0,0)= R::rbinom(n,0.5);
  double x1 = X(0,0);
  X(0,1)= R::rbeta(x1+a,n-x1+b);
  
  for(int t=1; t < N; t++){
    X(t,0) = R::rbinom(n,X(t-1,1));
    X(t,1) = R::rbeta(X(t,0)+a,n-X(t,0)+b);
  }
  return X;
}')

n<-10
a<-b<-1
N<-1e5

library(stats)
set.seed(1234)
## R function v.s. Rcpp function
C = cf(n,a,b,N)
R = f(n,a,b,N)
R1 = mahalanobis(R,colMeans(R),cov(R))
C1 = mahalanobis(C,colMeans(C),cov(C))
qqplot(R1,C1,main ="R function v.s. Rcpp function" )
abline(0, 1, col = 'black')


library(microbenchmark)
# Compare the computation time of the two functions with the function “microbenchmark”.
ts = microbenchmark(R = f(n,a,b,N), C = cf(n,a,b,N))
summary(ts)[,c(1,3,5,6)]

