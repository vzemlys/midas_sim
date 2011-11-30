##############
### PIRMAS ###
##############
rm(list=ls())
iter<-10000
steb<- c(250,500,1000,2000)# - jei m-12, tai nuo 200
type<-"ar" ## "ma" arba "ar"
rho.v<-0.75
m <- 4
l.0<-c(-10,-10)      # c(-10,-10)
a.0<--0.2 #  -0.2
b.0<-10 # -10
sd.x <- 1
sd.y <- 1
laipsn <- 0 ## 0 - kmax= Scwert12, kitaip kmax=n^laipsn

ptm <- proc.time()
###############
###############
#library("nls2",lib="/home/scratch/R/site-library")
library("multicore",lib="/home/scratch/R/site-library")
library("iterators",lib="/home/scratch/R/site-library")
library("foreach",lib="/home/scratch/R/site-library")
library("doMC",lib="/home/scratch/R/site-library")
#library("stats",lib="/home/scratch/R/site-library")

#library("Rcpp",lib="/home/scratch/lib64/R/library")
#library("inline",lib="/home/scratch/lib64/R/library")
#library("MASS",lib="/home/scratch/lib64/R/library")

registerDoMC(16)
#set.seed(1234)
#################
### Funkcijos ###
#################
theta.h0 <- function(index,alpha,beta,lambda) {
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    (alpha+beta*i)*exp(pol)
}
theta.h1 <- function(index,alpha,beta,lambda) {
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    epol <- exp(pol)
    alpha+beta*epol/sum(epol)
}
dgp.x <- function(n.x,rho,sd,burn.in=300,...) {
    e <- rnorm(n.x+burn.in,sd=sd.x)
	if(type=="ar"){res <- filter(e,filter=rho,method="recursive",sides=1)}
	if(type=="ma"){res<-c(0,e[2:length(e)]+rho*e[1:(length(e)-1)])}
    res[-burn.in:-1]
}
gen.IC <- function(yx,kmax,theta,n){
    res <- foreach(k=0:(kmax-1),.combine=rbind) %do% {
        c.1 <- k+1
        mod <- lm(yx[c.1:n,1]~yx[c.1:n,2:(c.1*m+1)]-1)
        nobs <- length(resid(mod))
        c(k,AIC(mod,k=2),AIC(mod,k=log(nobs)),AIC(mod,k=2*log(log(nobs))),IC(mod,k,nobs),APE(mod),PPE(mod,k,theta),MSE(mod,k,theta))
    }
    colnames(res) <- c("Lag","AIC","BIC","HQ","KZIC","APE","PPE","MSE")
    res
}
IC <- function(object,k,n) {
    xx <- 1/min(eigen(t(object$model[,-1])%*%object$model[,-1],symmetric=TRUE)$values)
    sigma2 <- sum(object$residuals^2)/(length(resid(object))-length(coef(object)))
    log(n-k)+log(xx)+log(sigma2)
}
APE <- function(object){
	sum(object$residuals^2)/(length(resid(object))-length(coef(object)))
}
PPE <- function(object,k,theta) {
	t(as.matrix(coef(object))-theta[1:length(coef(object))])%*%t(as.matrix(object$model[2:length(object$model)])) %*% as.matrix(object$model[2:length(object$model)])%*%(as.matrix(coef(object))-theta[1:length(coef(object))])
}
MSE<-function(object,k,theta) {
	 t(as.matrix(coef(object))-theta[1:length(coef(object))])%*%(as.matrix(coef(object))-theta[1:length(coef(object))]) /length(as.matrix(coef(object)))

}
dgp <- function(n,n.x,kmax,theta){
                theta.d <- cbind(theta,1:n.x)[order(-(1:n.x)),1]
                x <- dgp.x(n.x,rho.v,sd.x)
                idx <- m*c((n.x/m-n+1):(n.x/m))
                X <- foreach(h.x=0:((kmax+1)*m-1), .combine='cbind')%do%{
                         x[idx-h.x]
                     }
                y <-foreach(h.y= 1:n,.combine='c')%do%{
                             t(x[1:(n.x-(n-h.y)*m)])%*%theta.d[((n-h.y)*m+1):n.x]

                       }+rnorm(n,sd=sd.y)
                yx <- cbind(y,X)
            }
##########
##########
##########
rez <- foreach(n.s=1:length(steb),.combine='cbind')%:%
          foreach(i=1:iter,.combine='rbind') %dopar% {
                n <- steb[n.s]
#                kmax <- max(trunc(n^laipsn),trunc(12*(n/100)^0.25))
                kmax <- trunc(n^laipsn)
                if(laipsn==0){kmax <- trunc(12*(n/100)^0.25)}
                n.x <- n*m+3000
                theta <- theta.h0(1:n.x,alpha=a.0,beta=b.0,lambda=l.0)
                IC6<-gen.IC(dgp(n,n.x,kmax,theta),kmax,theta,n)
		ic<-c(apply(IC6[,2:5],2,which.min),trunc(4*(n/100)^0.25),trunc(12*(n/100)^0.25),trunc(n^0.33),trunc(n^0.3))
		names(ic) <- c("AIC","BIC","HQ","X","Schw4","Schw12","n0.33","n0.3")
		prec<-c(IC6[ic[1],6],IC6[ic[2],6],IC6[ic[3],6],IC6[ic[4],6],IC6[ic[5],6],IC6[ic[6],6],IC6[ic[7],6],IC6[ic[8],6],IC6[ic[1],7],IC6[ic[2],7],IC6[ic[3],7],IC6[ic[4],7],IC6[ic[5],7],IC6[ic[6],7],IC6[ic[7],7],IC6[ic[8],7],IC6[ic[1],8],IC6[ic[2],8],IC6[ic[3],8],IC6[ic[4],8],IC6[ic[5],8],IC6[ic[6],8],IC6[ic[7],8],IC6[ic[8],8])
		names(prec) <- c("APE_aic","APE_bic","APE_hq","APE_x","APE_schw4","APE_schw12","APE_n0.33","APE_n0.3","PPE_aic","PPE_bic","PPE_hq","PPE_x","PPE_schw4","PPE_schw12","PPE_n0.33","PPE_n0.3","MSE_aic","MSE_bic","MSE_hq","MSE_x","MSE_schw4","MSE_schw12","MSE_n0.33","MSE_n0.3")
                names(n) <- c("n")
		c(ic,prec,n,1e+10)
                }
means<-apply(rez,2,mean)
quart<-apply(rez,2,quantile, probs = c(0,0.25,0.5,0.75,1),na.rm=T)
nas<-colSums(is.na(rez))
non.nas<-colSums(!is.na(rez))
#quart
#means

kiek <- length(rez[1,])/length(steb)
rez.1 <- foreach(n.s=1:length(steb))%do%{
    n <- steb[n.s]
    quart.n<-sprintf("quart.%1.f",n)
    assign(quart.n,quart[,((n.s-1)*kiek+1):(n.s*kiek)]) # atvirkscias yra get() !!!
    means.n<-sprintf("means.%1.f",n)
    assign(means.n,means[((n.s-1)*kiek+1):(n.s*kiek)]) # atvirkscias yra get() !!!
}
#####################
#####################
#####################
s.file<-sprintf("/home/scratch/Virmantas/dat_midas/0_I_schwert/IC_iter.%1.f_type.%s%.2f_m.%1.f.Rdata",iter,type,rho.v,m)
save(iter,steb,type,rho.v,m,a.0,b.0,l.0,sd.x,sd.y,rez,rez.1,means,quart,nas,non.nas,file=s.file)

trukme <-(proc.time()- ptm)/60
trukme
