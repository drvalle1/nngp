rm(list=ls())
set.seed(1)

#get useful functions
setwd('U:\\GIT_models\\NNGP')
source('aux functions.R')

#get data
setwd('U:\\GIT_models\\NNGP\\simulat')
dat=read.csv('1 sim data.csv')

#if I have K(theta), can I get matrices A and D? Yes
n=nrow(dat)
coord=dat[,c('x','y')]

#get neighbors
neighb.max=5
N=get.neighb(x=coord$x,
             y=coord$y,
             neighb.max=neighb.max)

#visualize pattern of neighbors
xord=sort(coord$x)
par(ask=T)
for (i in 2:n){
  plot(y~x,data=coord)
  
  #get focus observation
  ind.focus=which(xord[i]==coord$x)
  points(y~x,data=coord[ind.focus,],col='red',pch=19)
  abline(v=coord$x[ind.focus],lty=3,col='red')
  
  #get neighbors
  ind.neighb=N[ind.focus,]
  tmp=ind.neighb[!is.na(ind.neighb)]
  points(y~x,data=coord[tmp,],col='red')
}

#export neighbors
setwd('U:\\GIT_models\\NNGP\\simulat')
write.csv(N,'2 neighbors.csv',row.names=F)