rm(list=ls())
library('Matrix')
set.seed(1)

#get useful functions
setwd('U:\\GIT_models\\NNGP')
source('aux functions.R')

#get data
setwd('U:\\GIT_models\\NNGP\\simulat')
dat=read.csv('1 sim data.csv')
nobs=nrow(dat)

#get design matrix
xmat=cbind(1,dat$cov1,dat$cov2)
ncov=ncol(xmat)
t.xmat=t(xmat)

#get neighbors
neighb=read.csv('2 neighbors.csv')

#calculate distance
coord=dat[,c('x','y')]
dist1=as.matrix(dist(coord))

#potential values of phi
max1=max(dist1)
med=max1/2
cor1=seq(from=0.01,to=0.5,length.out=30)
phi.pot=-log(cor1)/med
nphi.pot=length(phi.pot)
seq1=seq(from=0,to=max1,length.out=100)
par(ask=F)
plot(NA,NA,ylim=c(0,1),xlim=c(0,max1))
for (i in 1:nphi.pot){
  lines(seq1,exp(-phi.pot[i]*seq1),col=i)
}

#get A and D based on different values of phi
log.detD=ty.invM.y=rep(NA,nphi.pot)
tx.invM.y=matrix(NA,ncov,nphi.pot)
tx.invM.x=array(NA,dim=c(ncov,ncov,nphi.pot))
for (i in 1:nphi.pot){
  print(i)
  tmp=get.DandA(n=nobs,
                N=neighb,
                coord=coord,
                theta=phi.pot[i])
  
  #get log.detD
  diag.elem=tmp$D
  log.detD[i]=sum(log(diag.elem))
  
  #get invM
  invD=Matrix(diag(1/diag.elem),sparse=T)
  i.minus.a=diag(1,nobs)-tmp$A
  invM=t(i.minus.a)%*%invD%*%i.minus.a
  
  #get other stuff
  tx.invM.y[,i]=as.numeric(t.xmat%*%invM%*%dat$resp)
  ty.invM.y[i]=as.numeric(t(dat$resp)%*%invM%*%dat$resp)
  tx.invM.x[,,i]=as.matrix(t.xmat%*%invM%*%xmat)
}

#export calculations
setwd('U:\\GIT_models\\NNGP\\simulat')
write.csv(log.detD,'3 precal log_detD.csv',row.names=F)
write.csv(tx.invM.y,'3 precal tx_invM_y.csv',row.names=F)
write.csv(ty.invM.y,'3 precal ty_invM_y.csv',row.names=F)
tmp=matrix(tx.invM.x,ncov*ncov*nphi.pot,1)
write.csv(tmp,'3 precal tx_invM_x.csv',row.names=F)
write.csv(phi.pot,'3 precal phi_pot.csv',row.names=F)