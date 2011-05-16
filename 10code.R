##The functions for the project

sim.X <- function(n,m,rho=0.5,sd=1,burn.in=300,...) {
    e <- rnorm(n*m+burn.in,sd=sd)
    res <- filter(e,filter=rho,method="recursive")
    res[-burn.in:-1]
}

sim.Y.finite <- function(n,m,k,sde=1,tfun,xfun,...) {
   x <- xfun(n+k+1,m,...)
   theta <- tfun(1:(m*(k+1))-1,...)
   X <- x[sort(length(x)-(1:(n*m)-1))]
   start <- length(x)-n*m
   midas <- foreach(i=1:n,.combine=c) %do% {
        sum(x[start:1+i*m]*theta)
    }
   Y <- midas+rnorm(n,sd=sde)
   res <- list(X=X,Y=Y,theta=theta,midas=midas,x=x,n=n,m=m)
   class(res) <- "midas_sim"
   res
}

sim.Y.finite.list <- function(n,m,k,M,sde=1,tfun,xfun,...) {
    res <- foreach(i=1:M,.combine=c) %do% {
        list(sim.Y.finite(n,m,k,sde=sde,tfun=tfun,xfun=xfun,...))
    }
    class(res) <- "midas_sim_list"
    res
}

reflow <- function(object,...) UseMethod("reflow")


reflow.midas_sim <- function(object,k) {
    m <- length(object$X)%/%length(object$Y)
    n <- length(object$Y)
    Y <- object$Y[(k+1):n]
    X <- rcppreflow(object$X,c(n,m,k))
    data.frame(Y,X)
}


theta.r214 <- function(index,lambda,beta,...) {
    pol <- poly((index+1),2,raw=TRUE) %*%lambda
    epol <- exp(pol)
    beta*epol/sum(epol)
}

theta.rg214 <- function(index,lambda0,lambda1,beta,...) {
    i <- (index+1)/100
    pl <- poly(i,2,raw=TRUE)
    pol <- pl %*%c(lambda0,lambda1)
    epol <- exp(pol)[,,drop=TRUE]
    sepol <- sum(epol)
    res <- beta*epol/sepol
    
    ple <- pl*epol
    sple <- colSums(ple)
    Z <- cbind(beta*sweep(pl*sepol,2,sple,"-")*epol/sepol^2,epol/sepol)
    dimnames(Z) <- list(NULL,c("lambda0","lambda1","beta"))
    attr(res,"gradient") <- Z
    res
}

theta.2r214 <- function(index,lambda0,lambda1,beta,...) {
    pol <- poly((index+1)/100,2,raw=TRUE) %*%c(lambda0,lambda1)
    epol <- exp(pol)
    beta*epol/sum(epol)
}
  
theta.u214 <- function(index,lambda,beta,...) {
    pol <- poly(index+1,2,raw=TRUE) %*%lambda
    beta*exp(pol)
}

theta.uab214 <- function(index,lambda,alpha,beta) {
    i <- (index+1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    (alpha+beta*i)*exp(pol)
}
theta.rs214 <- function(index,lambda,beta,...) {
    pol <- poly((index+1),2,raw=TRUE) %*%lambda
    epol <- exp(pol)*((-1)^index)
    beta*epol/sum(epol)
}

theta.us214 <- function(index,lambda,beta,...) {
    pol <- poly(index+1,2,raw=TRUE) %*%lambda
    beta*exp(pol)*((-1)^index)
}


prep.nls.finite.old <- function(object,k) {
    mod <- lm(Y~.-1,data=reflowold.midas_sim(object,k=k))
    that <- coef(mod)
    lst <- list(that=that,index=1:(object$m*(k+1))-1)
    lst
}

prep.nls.finite <- function(object,k) {
    dt <- reflow(object,k=k)
    that <- coef(lsfit(dt[,-1],dt[,1],intercept=FALSE))
    list(that=that,index=1:(object$m*(k+1))-1)
}
  
gen.IC <- function(object,kmax) {
    res <- foreach(k=0:kmax,.combine=rbind) %do% {
        mod <- lm(Y~.-1,data=reflow(object,k))
        c(k,AIC(mod),BIC(mod),KZIC(mod))
    }
    colnames(res) <- c("Lag","AIC","BIC","KZIC")
    res
}

KZIC <- function(object) {
    xx <- 1/min(svd(object$model[,-1])$d)
    n <- nrow(object$model)
    k <- ncol(object$model)-1
    sigma2 <- sum(object$residuals^2)/(n-k)
    (log(n)+log(xx)+log(sigma2))/2
    
}

fit.lambda<- function(object,formula,k,...)UseMethod("fit.lambda")

fit.lambda.midas_sim<- function(object,formula,k,...) {
##formula must be of form that~function(index,...)
    bb <- prep.nls.finite(object,k=k)
    res <- try(nls(formula,data=bb,...))
    res
}

fit.lambda.midas_sim_list <- function(object,formula,k,...) {
##formula must be of form that~function(index,...)
    if(length(k)==1) res <- lapply(object,fit.lambda,formula=formula,k=k,...)
    else res <- foreach(l=object,ki=k,.combine=c) %do% list(fit.lambda(l,formula,ki,...))
    res
}

compare.lambda <- function(x,y) {
    if(length(x)!=length(y))stop("Incomparable objects")
    foreach(xx=x,yy=x,.combine=c) %do% {
        list(cbind(predict(xx),predict(yy)))
    }
}
