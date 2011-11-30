##################
### PIRMAS-rez ###
##################
rm(list=ls())
iter<-1000
steb<- c(125,250,500,1000,2000)# - jei m-12, tai nuo 200
type<-"ar" ## "ma" arba "ar"
rho.v<-0.75
m <- 4
l.0<-c(-10,-10)      # c(-10,-10)
a.0<--0.2 #  -0.2
b.0<-10 # -10
sd.x <- 1
sd.y <- 1
###############
###############
###############
library("foreach",lib="/home/scratch/R/site-library")
s.file<-sprintf("/home/scratch/Virmantas/dat_midas/0_I/IC_iter.%1.f_type.%s%.2f_m.%1.f.Rdata",iter,type,rho.v,m)
load(file=s.file)
###############
kiek <- length(rez[1,])/length(steb)
rez.1 <- foreach(n.s=1:length(steb))%do%{
    n <- steb[n.s]
    quart.n<-sprintf("quart.%1.f",n)
    assign(quart.n,quart[,((n.s-1)*kiek+1):(n.s*kiek)])
    means.n<-sprintf("means.%1.f",n)
    assign(means.n,means[((n.s-1)*kiek+1):(n.s*kiek)])
}
###############
###############
###############
#n.s <- 1000
#cbind(t(get(sprintf("quart.%1.f",n.s)))[1:4,],cbind(get(sprintf("means.%1.f",n.s))[5:8],get(sprintf("means.%1.f",n.s))[9:12],get(sprintf("means.%1.f",n.s))[13:16]))

#t(get(sprintf("quart.%1.f",n.s)))[1:4,]
#cbind(get(sprintf("means.%1.f",n.s))[5:8],get(sprintf("means.%1.f",n.s))[9:12],get(sprintf("means.%1.f",n.s))[13:16])
###############
###############
###############

kint <- 8

a <- rbind(cbind(t(get(sprintf("quart.%1.f",250)))[1:kint,],cbind(get(sprintf("means.%1.f",250))[(kint+1):(2*kint)],get(sprintf("means.%1.f",250))[(2*kint+1):(3*kint)],get(sprintf("means.%1.f",250))[(3*kint+1):(4*kint)])),cbind(t(get(sprintf("quart.%1.f",500)))[1:kint,],cbind(get(sprintf("means.%1.f",500))[(kint+1):(2*kint)],get(sprintf("means.%1.f",500))[(2*kint+1):(3*kint)],get(sprintf("means.%1.f",500))[(3*kint+1):(4*kint)])),cbind(t(get(sprintf("quart.%1.f",1000)))[1:kint,],cbind(get(sprintf("means.%1.f",1000))[(kint+1):(2*kint)],get(sprintf("means.%1.f",1000))[(2*kint+1):(3*kint)],get(sprintf("means.%1.f",1000))[(3*kint+1):(4*kint)])),cbind(t(get(sprintf("quart.%1.f",2000)))[1:kint,],cbind(get(sprintf("means.%1.f",2000))[(kint+1):(2*kint)],get(sprintf("means.%1.f",2000))[(2*kint+1):(3*kint)],get(sprintf("means.%1.f",2000))[(3*kint+1):(4*kint)])))
a<- cbind(a[,1:6],a[,8])
colnames(a) <- c("0%","25%","50%","75%","100%","PMSE","MSE")
a
