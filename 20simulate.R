sim1 <-  sim.Y(100,3,sde=1,tfun=theta.ar,xfun=sim.X)

summary(lm(Y~.,data=reflow(sim1,k=3)))


source("10code.R")
library(foreach)

sim3 <- sim.Y.finite(100,3,5,tfun=theta.216,xfun=sim.X,lambda=c(1,2,-3))


modl <- gen.k(sim3,17)
k <- which.min(sapply(modl,KZ.k))-1
bb <- prep.nls.finite(sim3,5)



summary(nls(that~theta.216(index,lambda),data=bb,start=list(lambda=c(1,1,1))))

summary(nls(that~theta.214(index,lambda),data=bb,start=list(lambda=c(1,1))))


##do100

sim3.100 <- foreach(i=1:100,.combine=c) %do% {
    list(sim.Y.finite(100,3,5,tfun=theta.216,xfun=sim.X,lambda=c(1,2,-3)))
}

pnls3.100.1 <- lapply(sim3.100,prep.nls.finite,k=1)

nls3.100.216 <- lapply(pnls3.100.1,function(l){
    nls(that~theta.216(index,lambda),data=l,start=list(lambda=c(1,1,1)))
})

nls3.100.214 <- lapply(pnls3.100.1,function(l){
    nls(that~theta.214(index,lambda),data=l,start=list(lambda=c(1,1)))
})

modl3.100 <- lapply(sim3.100,gen.k,kmax=17)

k3.100 <- lapply(modl3.100,function(l)which.min(sapply(l,KZ.k))-1)
