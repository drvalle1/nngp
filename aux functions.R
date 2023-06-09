get.neighb=function(x,y,neighb.max){
  xord=sort(x)
  dist1=data.matrix(dist(cbind(x,y)))
  n=length(x)
  neighb=matrix(NA,n,neighb.max)
  for (i in 2:n){
    ind.focus=which(xord[i]==x)
    ind.pot=which(x<xord[i]) #potential neighbors
    n.pot=length(ind.pot)    #number of potential neighbors
    if (n.pot <= neighb.max) sel2=ind.pot
    if (n.pot >  neighb.max){
      n.pot=neighb.max
      tmp=dist1[ind.focus,ind.pot]
      sel1=order(tmp,decreasing = F) #order index of closest neighbors 
      sel2=ind.pot[sel1[1:neighb.max]]
    }
    neighb[ind.focus,1:n.pot]=sel2
  }
  neighb
}
#------------------------------
get.DandA=function(n,N,coord,theta){
  A=Matrix(0,nrow=n,ncol=n,sparse=T)
  D=rep(NA,n)
  D[1]=1
  for (i in 1:(n-1)){
    tmp=N[i+1,]
    Pa=tmp[!is.na(tmp)]
    ind=c(i+1,Pa)
    dist1=data.matrix(dist(coord[ind,]))
    K=exp(-theta*dist1)
    A[i+1,Pa]=solve(K[-1,-1],K[1,-1])
    D[i+1]=K[1,1]-K[1,-1]%*%A[i+1,Pa]
  }
  list(A=A,D=D)
}
