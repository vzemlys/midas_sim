##############
### ANTRAS ###
##############
rm(list=ls())
iter<-2000
steb<-c(125,250,500,1000,2000)#c(250,500,1000,2000)
type<-"ar" ## "ma" arba "ar"
rho.v<-0.75
m <- 4
info.type <- 4 # AIC - 1, BIC - 2, HQ - 3, KZIC - 4
laipsn <- 0 ## 0 - kmax= Scwert12, kitaip kmax=n^laipsn
k0 <- 4 ## 0 - L=all pagal IC, >0 - nurodytas, tik perstumtas per viena, t.y. k0=L-1

no.of.coef <- m*(k0+1)

l.0<-c(-20,-10)
a.0<- -0.1
b.0<-10
sd.x <- 1
sd.y <- 1

dal.interv <- 6

init.0 <- c(alpha=a.0,beta=b.0,lambda.1=l.0[1],lambda.2=l.0[2])
init<- expand.grid(alpha=seq(-1, 1, len = dal.interv),beta=seq(-20, 30, len = dal.interv), lambda.1=seq(-30, 20, len = dal.interv),lambda.2=seq(-30, 20, len = dal.interv))

min.iter.sk <- 1 # daugiau kokio min realizaciju skaiciaus skaiciuojama galia ir reiksmingumas (galima paskui suskaiciuoti rez.1-ame, kai rez.0 jau yra)

ptm <- proc.time()
###############
###############
###############
library("nls2",lib="/home/scratch/R/site-library")
library("multicore",lib="/home/scratch/R/site-library")
library("iterators",lib="/home/scratch/R/site-library")
library("foreach",lib="/home/scratch/R/site-library")
library("doMC",lib="/home/scratch/R/site-library")
#library("stats",lib="/home/scratch/R/site-library")
#library("Rcpp",lib="/home/scratch/lib64/R/library")
#library("inline",lib="/home/scratch/lib64/R/library")
library("MASS",lib="/home/scratch/lib64/R/library")

registerDoMC(16)
#set.seed(1234)
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
powerf<-function(obj,krit,reiksm) {
    c((sum(obj[,2]>krit,na.rm=T)+sum(is.na(obj[,2])))/length(obj[,2]),(sum(obj[,2]>qchisq(1-reiksm,df=obj[,10]),na.rm=T)+sum(is.na(obj[,2])))/length(cbind(obj[,2],obj[,10])[,1]))
}
powerf.NA<-function(obj,krit,reiksm) {
    c(sum(obj[,2]>krit,na.rm=T)/length(na.omit(cbind(obj[,2],krit))[,1]),sum(obj[,2]>qchisq(1-reiksm,df=obj[,10]),na.rm=T)/length(cbind(obj[,2],obj[,10])[,1]))
}
set.seed(1)
#################
#################
#################
rez.0 <- foreach(n.s=1:length(steb),.combine='rbind')%:%
          foreach(i=1:iter,.combine='rbind') %dopar% {
                n <- steb[n.s]
                kmax <- trunc(n^laipsn)
                if(laipsn==0){kmax <- trunc(12*(n/100)^0.25)}
                n.x <- n*m+3000
                theta <- theta.h0(1:n.x,alpha=a.0,beta=b.0,lambda.1=l.0[1],lambda.2=l.0[2])
                yx <- dgp(n,n.x,kmax,theta)
                IC6<-gen.IC(yx,kmax,theta,n)
                k.x<-apply(IC6[,2:5],2,which.min)[info.type]
                c.1 <- k.x+1
                mod <- lm(yx[k.x:n,1]~yx[k.x:n,2:(c.1*m+1)]-1)

                n.o <- length(coef(mod))
                if(no.of.coef>0){
                    n.o <- no.of.coef
#                    c.1 <- k0+1
                }
                if(n.o>length(coef(mod))){
                    c.1 <- n.o/m
                    mod <- lm(yx[k.x:n,1]~yx[k.x:n,2:(n.o+1)]-1)
                }

#                if(n.o<=length(coef(mod))){
                    dat <- data.frame(that=coef(mod)[1:n.o],index=1:n.o)
                    xtx <- t(yx[k.x:n,2:(c.1*m+1)])%*%yx[k.x:n,2:(c.1*m+1)]
                    xtx <- xtx[1:n.o,1:n.o] # ar nereiktu daugiklio pakoreguojancio del n ir n.o skirtumo ???
                    Axtx<-solve(xtx)
                    Ch <- chol(Axtx)
                    W <- solve(t(Ch))
 #               }
#                W <- t(chol(xtx))

                sol.meth.h0 <- 0
                h0.try <- try({
                    sol.meth.h0 <- 1
                    lamb <- c(a.0,b.0,l.0[1],l.0[2])
                    names(lamb) <- c("alpha","beta","lambda.1","lambda.2")
                    g <- function(lamb) {coef(mod)[1:n.o]-theta.h0(1:n.o,lamb[1],lamb[2],lamb[3],lamb[4])}
                    fn0 <- function(lamb) {t(g(lamb))%*%xtx%*%g(lamb)}
                    e0 <- optim(lamb,fn0)[[1]]
                    eq.u <- list(nobs=length(dat[,1]),residuals=(dat[,1]-theta.h0(dat[,2],alpha=e0[1],beta=e0[2],lambda.1=e0[3],lambda.2=e0[4])),coefficients=e0)
                    gr.h0<-make.g.h0(eq.u)
                   h0 <- t(W%*%resid(eq.u))%*%ginv(sd(resid(mod))^2*(diag(length(resid(eq.u)))-W%*%gr.h0%*%ginv(t(W%*%gr.h0)%*%W%*%gr.h0)%*%t(W%*%gr.h0))%*%t(diag(length(resid(eq.u)))-W%*%gr.h0%*%ginv(t(W%*%gr.h0)%*%W%*%gr.h0)%*%t(W%*%gr.h0)))%*%W%*%resid(eq.u)
                })
                if(class(h0.try)=="try-error") h0 <- NA

                sol.meth.h1 <- 0
                h1.try <- try({
                    sol.meth.h1 <- 1
                    g <- function(lamb) {coef(mod)[1:n.o]-theta.h1(1:n.o,lamb[1],lamb[2],lamb[3],lamb[4])}
                    fn1 <- function(lamb) {t(g(lamb))%*%xtx%*%g(lamb)}
                    e1 <- optim(lamb,fn1)[[1]]
                    eq.r <- list(nobs=length(dat[,1]),residuals=(dat[,1]-theta.h1(dat[,2],alpha=e1[1],beta=e1[2],lambda.1=e1[3],lambda.2=e1[4])),coefficients=e1)
                    gr.h1<-make.g.h1(eq.r)
                   h1 <- t(W%*%resid(eq.r))%*%ginv(sd(resid(mod))^2*(diag(length(resid(eq.r)))-W%*%gr.h1%*%ginv(t(W%*%gr.h1)%*%W%*%gr.h1)%*%t(W%*%gr.h1))%*%t(diag(length(resid(eq.r)))-W%*%gr.h1%*%ginv(t(W%*%gr.h1)%*%W%*%gr.h1)%*%t(W%*%gr.h1)))%*%W%*%resid(eq.r)
                })
                if(class(h1.try)=="try-error") h1 <- NA
                ifelse(!is.na(h0),
#                    {vec <- c(h0,h1,k.x,kmax,e1,1111,1112,0,sol.meth.h0,sol.meth.h1,n,1e+10)

                    {vec <- c(h0,h1,k.x,kmax,coefficients(eq.u)[1:4],length(resid(eq.u)),length(resid(eq.u))-length(coef(eq.u)),1113,sol.meth.h0,sol.meth.h1,n,1e+10)
                    names(vec) <- c("H0","H1","IC","kmax","lambda.1","lambda.2","alpha","beta","nobs_eq.u","df_eq.u","omit.last","sol.meth.H0","sol.meth.H1","n","1e+10")
                    },
                    {vec <- c(h0,h1,k.x,kmax,c(1:6)*NA,0,sol.meth.h0,sol.meth.h1,n,1e+10)
                    })
                    names(vec) <- c("H0","H1","IC","kmax","lambda.1","lambda.2","alpha","beta","nobs_eq.u","df_eq.u","omit.last","sol.meth.H0","sol.meth.H1","n","1e+10")
                    vec
       }
apply(rez.0,2,mean)
summary(rez.0[,1])
summary(rez.0[,2])
#rez.0
#################
#################
#################
dfmm <- foreach(n.s=1:length(steb),.combine='rbind')%do%{
   a <- c(steb[n.s],no.of.coef/m-1,no.of.coef/m-1)
   if(no.of.coef==0){a <- c(steb[n.s],min(rez.0[rez.0[,14]==steb[n.s],3]),max(rez.0[rez.0[,14]==steb[n.s],3]))}
   a
}
rez.1 <- foreach(n.s=1:length(steb),.combine='rbind')%:%
          foreach(z.s=(m*c(dfmm[dfmm[,1]==steb[n.s],2]:dfmm[dfmm[,1]==steb[n.s],3])),.combine='rbind')%dopar%{
              n <- steb[n.s]
              rez <- rez.0[rez.0[,14]==n&rez.0[,10]==z.s,]
              ifelse(length(rez)>15*min.iter.sk,{
              nas <- apply(is.na(rez),2,sum)
              non.nas <- apply(!is.na(rez),2,sum)
              sol.meth.h0 <- c(sum(na.omit(rez[,12]==0)),sum(na.omit(rez[,12]==1)),sum(na.omit(rez[,12]==2)),sum(na.omit(rez[,12]==3)))/length(rez[,1])*100
              sol.meth.h1 <- c(sum(na.omit(rez[,13]==0)),sum(na.omit(rez[,13]==1)),sum(na.omit(rez[,13]==2)),sum(na.omit(rez[,13]==3)))/length(rez[,1])*100
              n.z <- length(rez[,1])
              quart<-apply(rez,2,quantile, probs = c(0.99,0.95,0.9),na.rm=T)
              ### sizes ###
              size.1proc <- (sum(rez[,1]>qchisq(0.99,df=z.s),na.rm=T)+sum(is.na(rez[,1])))/length(rez[,1])
              size.1proc.NA<- sum(rez[,1]>qchisq(0.99,df=z.s),na.rm=T)/sum(!is.na(rez[,1]))
              size.5proc <- (sum(rez[,1]>qchisq(0.95,df=z.s),na.rm=T)+sum(is.na(rez[,1])))/length(rez[,1])
              size.5proc.NA<- sum(rez[,1]>qchisq(0.95,df=z.s),na.rm=T)/sum(!is.na(rez[,1]))
              size.10proc <- (sum(rez[,1]>qchisq(0.90,df=z.s),na.rm=T)+sum(is.na(rez[,1])))/length(rez[,1])
              size.10proc.NA<- sum(rez[,1]>qchisq(0.90,df=z.s),na.rm=T)/sum(!is.na(rez[,1]))
              size <-t(c(size.1proc, size.5proc,size.10proc))
              size.NA <- c(size.1proc.NA, size.5proc.NA,size.10proc.NA)
              ### power ###
              power.1proc <-  powerf(obj=rez,krit=quart[1,1],reiksm=0.01)
              power.1proc.NA <-  powerf.NA(obj=rez,krit=quart[1,1],reiksm=0.01)
              power.5proc <-  powerf(obj=rez,krit=quart[2,1],reiksm=0.05)
              power.5proc.NA <-  powerf.NA(obj=rez,krit=quart[2,1],reiksm=0.05)
              power.10proc <-  powerf(obj=rez,krit=quart[3,1],reiksm=0.1)
              power.10proc.NA <-  powerf.NA(obj=rez,krit=quart[3,1],reiksm=0.1)
              power <- c(power.1proc, power.5proc,power.10proc)
              power.NA <- c(power.1proc.NA, power.5proc.NA,power.10proc.NA)
              vect <- c(n,z.s,n.z,t(size),t(size.NA), t(power), t(power.NA),t(nas[1:2]),t(sol.meth.h0),t(sol.meth.h1))
              },
              vect <- c(n,z.s,t(c(1:29)*NA)))
              names(vect) <- c("n","df","no of df cases","size_1","size_5","size_10","size.beNA_1","size.beNA_5","size.beNA_10","power.adj_1","power.nom_1","power.adj_5","power.nom_5","power.adj_10","power.nom_10","power.beNA.adj_1","power.beNA.nom_1","power.beNA.adj_5","power.beNA.nom_5","power.beNA.adj_10","power.beNA.nom_10","nas.H0","nas.H1","sol.meth.H0=0","sol.meth.H0=1","sol.meth.H0=2","sol.meth.H0=3","sol.meth.H1=0","sol.meth.H1=1","sol.meth.H1=2","sol.meth.H1=3")
              vect
}
rez.1
rez.2 <- foreach(n.s=1:length(steb),.combine='rbind')%dopar%{
         n <- steb[n.s]
         a <- t(rez.1[rez.1[,1]==n,])
         if(length(rez.1[rez.1[,1]==n,1])>1){a <- apply(rez.1[rez.1[,1]==n,],2,mean,na.rm=T)}
         a
}
rez.2

trukme <- (proc.time()-ptm)/60
trukme
#####################
#####################
#####################
s.file<-sprintf("/home/scratch/Virmantas/dat_midas/_II_k0%1.f_beW/H0_iter.%1.f_steb.%1.f_type.%s.%.2f_m.%1.f_ic.%1.f_omit.%1.f_power.%.2f_dal.interv.%1.f.Rdata",k0,iter,steb,type,rho.v,m,info.type,0,laipsn,dal.interv)[1]

#save(iter,steb,type,rho.v,m,l.0,a.0,b.0,info.type,0,sd.x,sd.y,rez.0,rez.1,rez.2,file=s.file)
