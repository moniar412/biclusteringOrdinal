st0=1+5*(0)
source("randmixBicl.R")
source("margX.R")
source("idx.R")
source("EMbicl.R")
source("parl2vr.R")
source("parv2lr.R")
source("M.step1.R")
source("int.R")
source("E.step1.R")
source("pcllikf1.R")
source("pcllikf1g.R")
source("A_update.R")
source("stPoint.R")
source("qmix.R")
source("t3postProb1.R")
source("t3int_post.R")
source("functionsStPoint.R")
source("simH4.R")
library(mvtnorm)
library(mclust)
#library(psych)
library("optimParallel")
library(Matrix)
library(polycor)

#in locale

#cores=detectCores()
#n_cor=round((cores-2))
#cl <- makeCluster(n_cor,outfile = "")
#setDefaultCluster(cl=cl) # set 'cl' as default cluster

#terastat
cores=detectCores()
n_cor=round((cores-2))
message("detectCores: ", cores)
message("round((cores-2)): ", n_cor)
n_cores <- as.integer(Sys.getenv("SLURM_NTASKS", unset = "1"))
# Mostra il numero dei core
print(paste("Numero di core allocati:", n_cores))
cl <- makeCluster(n_cores,outfile = "")
setDefaultCluster(cl=cl) # set 'cl' as default cluster


N=c(1000, 2500)
P=8; Q=c(2, 3, 3); G=3;
no.TH=rep(5,P);nk=100


out=simH4(N,P,Q,G,no.TH,nk,st0,cl)
