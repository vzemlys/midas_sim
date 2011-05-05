##The functions for the project

sim.X <- function(n,m,rho=0.5,sd=1,burn.in=300,...) {
    e <- rnorm(n*m+burn.in,sd=sd)
    res <- filter(e,filter=rho,method="recursive")
    res[-burn.in:-1]
}

sim.Y.infty<- function(n,m,infty=10,sde=1,tfun,xfun,...) {
    x <- xfun((infty+1)*n,m,...)
    theta <- tfun(expand.grid(1:m,1:(infty*n)-1),...)
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

reflow <- function(object,...) UseMethod("reflow")

reflowold.midas_sim <- function(object,k) {
    m <- length(object$X)%/%length(object$Y)
    n <- length(object$Y)
    Y <- object$Y[(k+1):n]
    X <- foreach(i=1:(n-k),.combine=rbind) %do% {
        object$X[(m*(k+1)):1+(i-1)*m]
    }
    data.frame(Y,X)
}

reflow.midas_sim <- function(object,k) {
    m <- length(object$X)%/%length(object$Y)
    n <- length(object$Y)
    Y <- object$Y[(k+1):n]
    X <- rcppreflow(object$X,c(n,m,k))
    data.frame(Y,X)
}

theta.216 <- function(index,lambda,...) {
    poly(1/(index+1),length(lambda)-1,raw=TRUE) %*%lambda[-1]+lambda[1]
}

theta.214 <- function(index,lambda,...) {
    pol <- poly(1/(index+1),2,raw=TRUE) %*%lambda
    epol <- exp(pol)
    epol/sum(epol)
}

theta.r214 <- function(index,lambda,beta,...) {
    pol <- poly((index+1),2,raw=TRUE) %*%lambda
    epol <- exp(pol)
    beta*epol/sum(epol)
}

theta.u214 <- function(index,lambda,beta,...) {
    pol <- poly(index+1,2,raw=TRUE) %*%lambda
    beta*exp(pol)
}

theta.ar <- function(grid,theta=0.3,...) {    
    theta^(grid[,2])
}

theta.mid <- function(grid,theta=0.7,...) {
    m <- max(grid[,1])
    h <- grid[,2]*m+grid[,1]
    1000*theta^h
}

theta.mid2 <- function(grid,theta=0.7,C=1000,...) {
    m <- max(grid[,1])
    h <- grid[,2]*m+grid[,1]
    C*theta^(h)
}

theta.sin <- function(grid,theta=0.7,...) {
    m <- max(grid[,1])
    h <- grid[,2]*m+grid[,1]
    1000*theta^(h)*sin(pi*(2*h+3)/6)
}

prep.nls.infty <- function(object,k) {
    mod <- lm(Y~.-1,data=reflow(object,k=k))
    that <- coef(mod)
    lst <- list(that=that,grid=expand.grid(1:object$m,0:k))
    lst
}
# mod2 <- nls(that~theta.mid(grid,theta),data=lst,start(theta=0.4))
 #list(ols=mod,nls=mod2)

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
  
gen.k <- function(object,kmax) {
    res <- foreach(k=0:kmax,.combine=rbind) %do% {
        mod <- lm(Y~.-1,data=reflow(object,k))
        c(k,AIC(mod),BIC(mod),KZ.k(mod))
    }
    colnames(res) <- c("Lag","AIC","BIC","KZIC")
    res
}

KZ.k <- function(object) {
    xx <- 1/min(svd(object$model[,-1])$d)
    n <- nrow(object$model)
    k <- ncol(object$model)-1
    sigma2 <- sum(object$residuals^2)/(n-k)
    (log(n)+log(xx)+log(sigma2))/2
    
}

#fit.lambda <- 

calc.with.best <- function(sim,olsmod,method=KZ.k ) {
##Method should be a function which gets lm object and returns
##the information criteria such as AIC, BIC. 

    M <- length(sim)
    if(length(models)!=M)stop("Number of simulations do not coincide")

    res <- foreach(i=1:M,.combine=c) %do% {

        k <- which.min(sapply(olsmod[[i]],method))-1
        bb <- prep.nls.finite(sim[[i]],k)
        
        sn4u <- try(nls(that~theta.u214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=FALSE))

        sn4r <- try(nls(that~theta.r214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=FALSE))
    
        list(list(k=k,bb=bb,u=sn4u,r=sn4r))
    }
    class(res) <- "midas_sim_modtheta"
    res
}

calc.with.fixed <- function(sim,olsmod,k=5) {
    M <- length(sim)
    if(length(models)!=M)stop("Number of simulations do not coincide")

    res <- foreach(i=1:M,.combine=c) %do% {

        bb <- prep.nls.finite(sim[[i]],k=k)
        
        sn4u <- try(nls(that~theta.u214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=FALSE))

        sn4r <- try(nls(that~theta.r214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=FALSE))
    
        list(list(k=k,bb=bb,u=sn4u,r=sn4r))
    }
    class(res) <- "midas_sim_modtheta"
    res
}

