sample.betas=function(param,ncov,tx.invM.x,tx.invM.y){
  #get covariance matrix
  invV.beta=diag(1,ncov)
  prec=tx.invM.x[,,param$ind.phi]+param$sig2*invV.beta
  var1=solve(prec)
  
  #get part of mean
  pmedia=tx.invM.y[,param$ind.phi]
  t(rmvnorm(1,var1%*%pmedia,var1))
}
sample.phi.sig2=function(param,nphi.pot,nsig2.pot,nobs,
                         sig2.pot,
                         log.detD,ty.invM.y,tx.invM.y,
                         tx.invM.x){
  
  #calculate stuff involving phi
  p1=-(1/2)*log.detD$x
  p2=ty.invM.y-t(2*t(param$betas)%*%tx.invM.y)
  for (i in 1:nphi.pot){
    p2[i]=p2[i]+t(param$betas)%*%tx.invM.x[,,i]%*%param$betas
  }
  p2=as.numeric(p2)
  
  #calculate stuff involving sig2
  p1.sig2=-(nobs/2)*log(sig2.pot)
  p2.sig2=-(1/(2*sig2.pot))
  
  combo=expand.grid(ind.phi=1:nphi.pot,
                    ind.sig2=1:nsig2.pot)
  combo$fim=p1.sig2[combo$ind.sig2]+p1[combo$ind.phi]+
            p2.sig2[combo$ind.sig2]*p2[combo$ind.phi]  

  #normalize
  fim=combo$fim-max(combo$fim)
  fim1=exp(fim)
  fim1=fim1/sum(fim1)
  # plot(1:(nphi.pot*nsig2.pot),unlist(fim1))
  
  #sample
  tmp=rmultinom(1,size=1,prob=unlist(fim1))
  ind=which(tmp==1)
  combo[ind,]
}