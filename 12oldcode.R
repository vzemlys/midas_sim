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

reflowold.midas_sim <- function(object,k) {
    m <- length(object$X)%/%length(object$Y)
    n <- length(object$Y)
    Y <- object$Y[(k+1):n]
    X <- foreach(i=1:(n-k),.combine=rbind) %do% {
        object$X[(m*(k+1)):1+(i-1)*m]
    }
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
