sim1 <-  sim.Y(100,3,sde=1,tfun=theta.ar,xfun=sim.X)

summary(lm(Y~.,data=reflow(sim1,k=3)))
