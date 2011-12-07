###TASK 2. Check whether generated process (AR) is really AR(1).
###On the other hand check both.

n <- steb[n.s]
kmax <- trunc(n^laipsn)
if(laipsn==0){kmax <- trunc(12*(n/100)^0.25)}
n.x <- n*m+3000

type <- "ar"


yx <- dgp(n,n.x,max(kmax,n.xout/m),theta)
