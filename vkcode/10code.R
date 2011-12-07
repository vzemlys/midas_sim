#################
### Funkcijos ###
#################
theta.h0 <- function(index,alpha,beta,lambda.1,lambda.2) {
    i <- (index-1)/100
    pol <- lambda.1*i+lambda.2*i^2
    (alpha+beta*i)*exp(pol)
}
theta.h1 <- function(index,alpha,beta,lambda.1,lambda.2) {
    i <- (index-1)/100
    pol <- lambda.1*i+lambda.2*i^2
    epol <- exp(pol)
    alpha+beta*epol/sum(epol)
}

##TASK 1: Here be gradients, check them!
make.g.h0<-function(object,...){
    alpha<-coefficients(object)[1]
    beta<-coefficients(object)[2]
    lambda<-c(coefficients(object)[3],coefficients(object)[4])
    index <- c(1:length(resid(object)))
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    a<-(alpha+beta*i)*exp(pol)
    cbind(a,a*i,a*i*(alpha+beta*i),a*i^2*(alpha+beta*i))
}

make.g.h1<-function(object,...){
    alpha<-coefficients(object)[1]
    beta<-coefficients(object)[2]
    lambda<-c(coefficients(object)[3],coefficients(object)[4])
    index <- c(1:length(resid(object)))
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    epol<-exp(pol)
    b0 <- epol/sum(epol)
    b1 <- sum(epol*i)/sum(epol)
    b2 <- sum(epol*i^2)/sum(epol)
    cbind(index/index,b0,beta*b0*(i-b1),beta*b0*(i^2-b2))
}

###TASK 2. Check whether generated process (AR) is really AR(1).
###On the other hand check both.
dgp.x <- function(n.x,rho,sd,burn.in=300,...) {
    e <- rnorm(n.x+burn.in,sd=sd.x)
	if(type=="ar"){res <- filter(e,filter=rho,method="recursive",sides=1)}
	if(type=="ma"){res<-c(0,e[2:length(e)]+rho*e[1:(length(e)-1)])}
    res[-burn.in:-1]
}

##TASK 3. Check whether optimisation is needed. Lm is usualy slow.
gen.IC <- function(yx,kmax,theta,n){
    res <- foreach(k=0:(kmax-1),.combine=rbind) %do% {
        c.1 <- k+1
        mod <- lm(yx[c.1:n,1]~yx[c.1:n,2:(c.1*m+1)]-1)
        nobs <- length(resid(mod))
        c(k,AIC(mod,k=2),AIC(mod,k=log(nobs)),AIC(mod,k=2*log(log(nobs))),IC(mod,k,nobs))
    }
    colnames(res) <- c("Lag","AIC","BIC","HQ","KZIC")
    res
}


IC <- function(object,k,n) {
    xx <- 1/min(eigen(t(object$model[,-1])%*%object$model[,-1],symmetric=TRUE)$values)
    sigma2 <- sum(object$residuals^2)/(length(resid(object))-length(coef(object)))
    log(n-k)+log(xx)+log(sigma2)
}

##TASK 4. Is this really MIDAS?
dgp <- function(n,n.x,km,theta){
                theta.d <- cbind(theta,1:n.x)[order(-(1:n.x)),1]
                x <- dgp.x(n.x,rho.v,sd.x)
                idx <- m*c((n.x/m-n+1):(n.x/m))
                X <- foreach(h.x=0:((km+1)*m-1), .combine='cbind')%do%{
                         x[idx-h.x]
                     }
                y <-foreach(h.y= 1:n,.combine='c')%do%{
                             t(x[1:(n.x-(n-h.y)*m)])%*%theta.d[((n-h.y)*m+1):n.x]

                       }+rnorm(n,sd=sd.y)
                yx <- cbind(y,X)
            }


powerf<-function(obj,krit,reiksm) {
    c((sum(obj[,2]>krit,na.rm=T)+sum(is.na(obj[,2])))/length(obj[,2]),(sum(obj[,2]>qchisq(1-reiksm,df=obj[,10]),na.rm=T)+sum(is.na(obj[,2])))/length(cbind(obj[,2],obj[,10])[,1]))
}


powerf.NA<-function(obj,krit,reiksm) {
    c(sum(obj[,2]>krit,na.rm=T)/length(na.omit(cbind(obj[,2],krit))[,1]),sum(obj[,2]>qchisq(1-reiksm,df=obj[,10]),na.rm=T)/length(cbind(obj[,2],obj[,10])[,1]))
}
