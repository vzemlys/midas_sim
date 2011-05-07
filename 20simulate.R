source("10code.R")
source("11code.R")

library(foreach)
library(stats4)
library(ggplot2)

sim5 <- sim.Y.finite(100,3,1000,tfun=theta.us214,xfun=sim.X,lambda=c(-0.01,-0.001),beta=1)


IC5 <- gen.IC(sim5,kmax=15)

sn5u <- fit.lambda(sim5,that~theta.us214(index,lambda,beta),k=12,start=list(lambda=c(-0.01,-0.001),beta=5),trace=TRUE)

sn5r <- fit.lambda(sim5,that~theta.rs214(index,lambda,beta),k=12,start=list(lambda=c(-0.01,-0.001),beta=5),trace=TRUE)


