rm(list=ls())
library('mvtnorm')
set.seed(1)

#get useful functions
setwd('U:\\GIT_models\\NNGP')
source('gibbs sampler functions.R')

#get data
setwd('U:\\GIT_models\\NNGP\\simulat')
dat=read.csv('1 sim data.csv')
nobs=nrow(dat)
xmat=cbind(1,dat$cov1,dat$cov2)
ncov=ncol(xmat)

#get precalc stuff
tx_invM_x=data.matrix(read.csv('3 precal tx_invM_x.csv'))
phi.pot=read.csv('3 precal phi_pot.csv')$x
nphi.pot=length(phi.pot)
log.detD=read.csv('3 precal log_detD.csv')
tx.invM.y=data.matrix(read.csv('3 precal tx_invM_y.csv'))
ty.invM.y=data.matrix(read.csv('3 precal ty_invM_y.csv'))
tmp=read.csv('3 precal tx_invM_x.csv')
tx.invM.x=array(unlist(tmp),dim=c(ncov,ncov,nphi.pot))

#potential sig2
var.resp=var(dat$resp)*2
nsig2.pot=30
sig2.pot=seq(from=var.resp/1000,to=var.resp,length.out=nsig2.pot)

#gibbs stuff
ngibbs=10000
tmp=floor(nsig2.pot/2)
param=list(betas=matrix(0,ncov,1),
           ind.sig2=tmp,
           sig2=sig2.pot[tmp],
           ind.phi=floor(nphi.pot/2))
store.betas=matrix(NA,ngibbs,ncov)
store.sig2.phi=matrix(NA,ngibbs,2)
for (i in 1:ngibbs){
  print(i)
  param$betas=sample.betas(param=param,
                           ncov=ncov,
                           tx.invM.x=tx.invM.x,
                           tx.invM.y=tx.invM.y)
  tmp=sample.phi.sig2(param=param,
                      nphi.pot=nphi.pot,
                      nsig2.pot=nsig2.pot,
                      nobs=nobs,
                      sig2.pot=sig2.pot,
                      log.detD=log.detD,
                      ty.invM.y=ty.invM.y,
                      tx.invM.y=tx.invM.y,
                      tx.invM.x=tx.invM.x)
  param$ind.phi=tmp$ind.phi
  param$ind.sig2=tmp$ind.sig2      
  param$sig2=sig2.pot[param$ind.sig2]
  
  #store results
  store.betas[i,]=param$betas
  store.sig2.phi[i,]=c(param$sig2,
                       phi.pot[param$ind.phi])
}

#examine results
burn.in=ngibbs/2
seq1=burn.in:ngibbs
plot(store.sig2.phi[seq1,2],type='l')
abline(h=3,col='red')

z=table(store.sig2.phi[seq1,2])
plot(z,type='h')
abline(v=3,col='red')

plot(store.sig2.phi[seq1,1],type='l')
abline(h=0.3,col='red')

cor(store.sig2.phi)

betas.true=c(0.5,0.3,0.1)
for (i in 1:ncov){
  plot(store.betas[seq1,i],type='l')
  abline(h=betas.true[i],col='red')
}