##The functions for the project

sim.X <- function(n,m,rho=0.5,sd=1,burn.in=300,...) {
    e <- rnorm(n*m+burn.in,sd=sd)
    res <- filter(e,filter=rho,method="recursive")
    res[-burn.in:-1]
}

sim.Y <- function(n,m,infty=10,sde=1,tfun,xfun,...) {
    x <- xfun((infty+1)*n,m,...)
    theta <- tfun(expand.grid(1:m,1:(infty*n)-1),...)
    X <- x[sort(length(x)-(1:(n*m)-1))]
    start <- length(x)-n*m
    midas <- foreach(i=1:n,.combine=c) %do% {
        sum(x[start:1+i*m]*theta)
    }
    Y <- 1+midas+rnorm(n,sd=sde)
    browser()
    res <- list(X=X,Y=Y,theta=theta,midas=midas,x=x)
    class(res) <- "midas_sim"
    res
}

reflow <- function(object,...) UseMethod("reflow")

reflow.midas_sim <- function(object,k) {
    m <- length(object$X)%/%length(object$Y)
    Y <- object$Y[(k+1):n]
    X <- foreach(i=1:(n-k),.combine=rbind) %do% {
        object$X[(m*(k+1)):1+(i-1)*m]
    }
    data.frame(Y,X)
}
  

theta.ar <- function(grid,theta=0.3,...) {    
    theta^(grid[,2])
}
