rm(list=ls())
set.seed(1)

#if I have K(theta), can I get matrices A and D? Yes
n=5
tmp=1:n
coord=data.frame(x=tmp,y=tmp)
dist1=as.matrix(dist(coord))
theta=0.5
K=exp(-theta*dist1)

#look at the implied exponential decay
seq1=seq(from=0,to=max(dist1),length.out=100)
plot(seq1,exp(-theta*seq1),type='l')

#notice how we fill A row by row
#we fill D element by element afterwards
D=A=matrix(0,n,n)
D[1,1]=K[1,1]
for (i in 1:(n-1)){
  A[i+1,1:i]=solve(K[1:i,1:i],K[1:i,i+1])
  D[i+1,i+1]=K[i+1,i+1]-K[i+1,1:i]%*%A[i+1,1:i]
}

tmp=solve(diag(1,n)-A)
Kestim=tmp%*%D%*%t(tmp)

rango=range(c(K,Kestim))
plot(K,Kestim,ylim=rango,xlim=rango)
lines(rango,rango,col='red')
