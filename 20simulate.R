source("10code.R")
source("11code.R")

library(foreach)
library(stats4)
library(ggplot2)

sim5 <- sim.Y.finite(100,4,1000,tfun=theta.us214,xfun=sim.X,lambda=c(-0.01,-0.001),beta=1)


IC5 <- gen.IC(sim5,kmax=15)

sn5u <- fit.lambda(sim5,that~theta.us214(index,lambda,beta),k=12,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=TRUE)

sn5r <- fit.lambda(sim5,that~theta.r214(index,lambda,beta),k=8,start=list(lambda=c(-2,-0.1),beta=-.5),trace=TRUE)

sim5.100 <- sim.Y.finite.list(100,3,1000,100,tfun=theta.us214,xfun=sim.X,lambda=c(-0.01,-0.001),beta=1)


sn5u.100 <- fit.lambda.midas_sim_list(sim5.100,that~theta.us214(index,lambda,beta),k=12,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=TRUE)

sn5r.100 <- fit.lambda.midas_sim_list(sim5.100,that~theta.r214(index,lambda,beta),k=12,start=list(lambda=c(0.2,-0.001),beta=-0.7),trace=TRUE)
