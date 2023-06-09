rm(list=ls())
set.seed(1)

setwd('/Volumes/Users/drvalle/independent studies/NNGP')
source('aux functions.R')

#if I have K(theta), can I get matrices A and D? Yes
n=20
coord=data.frame(x=runif(n,min=-10,max=10),
                 y=runif(n,min=-10,max=10))
coord=coord[order(coord$x),]
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

#notice how we fill A row by row
#we fill D element by element afterwards
D=A=matrix(0,n,n)
D[1,1]=K[1,1]
for (i in 1:(n-1)){
  tmp=N[i+1,]
  Pa=tmp[!is.na(tmp)]
  A[i+1,Pa]=solve(K[Pa,Pa],K[i+1,Pa])
  D[i+1,i+1]=K[i+1,i+1]-K[i+1,Pa]%*%A[i+1,Pa]
}

tmp=solve(diag(1,n)-A)
Kestim=tmp%*%D%*%t(tmp)

rango=range(c(K,Kestim))
plot(K,Kestim,ylim=rango,xlim=rango)
lines(rango,rango,col='red')

hist(K-Kestim)

image(A!=0) #A is sparse (in each row, at most 5 non-zero numbers)