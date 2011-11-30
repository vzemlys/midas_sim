#########################
### ANTRAS-NLS2-k-rez ###
#########################
rm(list=ls())
iter<-1000
steb <- c(250,500,1000,2000)
type<-"ar" ## "ma" arba "ar"
m <- 4
info.type <- 3 # AIC - 1, BIC - 2, HQ - 3, KZIC - 4
omit <- 0 # number of last estimated coefficients to omit
laipsn <- 0 #5/17
L <- 6 ## 0 - L=all, >0 - L=nurodytas

dal.interv <- 6

#l.0<-c(-10,-10)
#a.0<- -0.2
#b.0<-10
#sd.x <- 1
#sd.y <- 1


library("multicore",lib="/home/scratch/R/site-library")
library("iterators",lib="/home/scratch/R/site-library")
library("foreach",lib="/home/scratch/R/site-library")
library("doMC",lib="/home/scratch/R/site-library")
registerDoMC(16)
###############
###############
###############
rho <- c(0.75,0.50,0.25)
stebej <- c(250,500,1000,2000)
a <- foreach(rho.s=1:length(rho),.combine='rbind')%do%{
      rho.v <-rho[rho.s]
      steb <- steb[1]
s.file<-sprintf("/home/scratch/Virmantas/dat_midas/_II_L%1.f/H0_iter.%1.f_steb.%1.f_type.%s.%.2f_m.%1.f_ic.%1.f_omit.%1.f_power.%.2f_dal.interv.%1.f.Rdata",L,iter,steb,type,rho.v,m,info.type,omit,laipsn,dal.interv)[1]
      load(file=s.file)

      cbind(numeric(1:length(rez.1[,1]))+rho[rho.s],rez.1[,1:6])
      cbind(rez.1[,10],rez.1[,12],rez.1[,14],(rez.1[,25]+rez.1[,26]),(rez.1[,29]+rez.1[,30]))

}
a
#rez.2
