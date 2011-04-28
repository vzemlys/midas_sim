sim1 <-  sim.Y(100,3,sde=1,tfun=theta.ar,xfun=sim.X)

summary(lm(Y~.,data=reflow(sim1,k=3)))


sim3 <- sim.Y.finite(100,3,5,tfun=theta.216,xfun=sim.X,lambda=c(1,2,-3))


modl <- gen.k(sim3,17)
k <- which.min(sapply(modl,KZ.k))-1
bb <- prep.nls.finite(sim3,5)


summary(nls(that~theta.216(index,lambda),data=bb,start=list(lambda=c(1,1,1))))
