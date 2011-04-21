##The functions for the project

sim.X <- function(n,m,rho=0.5,sd=1,burn.in=300,...) {
    e <- rnorm(n*m+burn.in,sd=sd)
    res <- filter(e,filter=rho,method="recursive")
    res[-burn.in:-1]
}

sim.Y <- function(n,m,infty=10,sde=1,tfun,xfun,...) {
    x <- xfun((infty+1)*n,m,...)
    theta <- tfun(expand.grid(1:m,1:(infty*n)-1),...)
    X <- x[length(x)-(1:(n*m)-1)]
    start <- length(x)-n*m
    midas <- foreach(i=1:n,.combine=c) %do% {
        sum(x[1:start+i*m]*theta)
    }
    Y <- midas+rnorm(n,sd=sde)
    res <- list(X=X,Y=Y)
    class(res) <- "midas_sim"
    res
}

reflow <- function(object,...) UseMethod("reflow")

reflow.midas_sim <- function(object,k) {
    m <- length(object$X)%/%length(object$Y)
    Y <- object$Y[(k+1):n]
    X <- foreach(i=(k+1):n,.combine=rbind) %do% {
        object$X[1:(m*k)+i*m]
    }
    data.frame(Y,X)
}
  

theta.ar <- function(grid,theta=0.3,...) {    
    theta^(grid[,2])
}
