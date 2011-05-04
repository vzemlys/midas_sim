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


###The work with 2.14

sim4 <- sim.Y.finite(100,3,1000,tfun=theta.u214,xfun=sim.X,lambda=c(-0.01,-0.001),beta=1)

modl <- gen.k(sim4,17)

k <- which.min(sapply(modl,KZ.k))-1

bb <- prep.nls.finite(sim4,k)

sn4u <- nls(that~theta.u214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=TRUE)

sn4r <- nls(that~theta.r214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=TRUE)


sim4.100 <- foreach(i=1:100,.combine=c) %do% {
    list(sim.Y.finite(100,3,1000,tfun=theta.u214,xfun=sim.X,lambda=c(-0.01,-0.001),beta=1))
}

modl4.100 <- lapply(sim4.100,gen.k,kmax=17)

sim4ur <- foreach(i=1:100,.combine=c) %do% {

    k <- which.min(sapply(modl4.100[[i]],KZ.k))-1
    bb <- prep.nls.finite(sim4.100[[i]],k)
    
    sn4u <- try(nls(that~theta.u214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=FALSE))

    sn4r <- try(nls(that~theta.r214(index,lambda,beta),data=bb,start=list(lambda=c(-0.01,-0.001),beta=0.5),trace=FALSE))
    
    list(list(k=k,bb=bb,u=sn4u,r=sn4r))
}

sim4dif <- foreach(i=1:100,.combine=c) %do% {
    bb <- sim4ur[[i]]$bb
    sn4r <- sim4ur[[i]]$r
    sn4u <- sim4ur[[i]]$u
    if(class(sn4r)=="try-error" | class(sn4u)=="try-error") NULL
    else
    list(cbind(theta.r214(bb$index,coef(sn4r)[1:2],coef(sn4r)[3]),  theta.u214(bb$index,coef(sn4u)[1:2],coef(sn4u)[3])))
}

