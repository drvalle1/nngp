rm(list=ls())
set.seed(1)

setwd('/Volumes/Users/drvalle/independent studies/NNGP')
source('aux functions.R')

#if I have K(theta), can I get matrices A and D? Yes
n=20
coord=data.frame(x=runif(n,min=-10,max=10),
                 y=runif(n,min=-10,max=10))
plot(y~x,data=coord)
dist1=as.matrix(dist(coord))
theta=0.5
K=exp(-theta*dist1)

#look at the implied exponential decay
seq1=seq(from=0,to=max(dist1),length.out=100)
plot(seq1,exp(-theta*seq1),type='l')

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
