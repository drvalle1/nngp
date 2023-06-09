sample.betas=function(param,ncov,tx.invM.x,tx.invM.y){
  #get covariance matrix
  invV.beta=diag(1,ncov)
  prec=tx.invM.x[,,param$ind.phi]+param$sig2*invV.beta
  var1=solve(prec)
  
  #get part of mean
  pmedia=tx.invM.y[,param$ind.phi]
  t(rmvnorm(1,var1%*%pmedia,var1))
}
sample.sig2=function(param,a.prior,b.prior,nobs,
                     ty.invM.y,tx.invM.x){
  a1=(nobs/2)+a.prior
  quad=ty.invM.y[param$ind.phi]-
    2*t(param$betas)%*%tx.invM.y[,param$ind.phi]+
    t(param$betas)%*%tx.invM.x[,,param$ind.phi]%*%param$betas
  b1=(1/2)*quad+(b.prior)
  1/rgamma(1,a1,b1)
}
sample.phi=function(param,nphi.pot,log.detD,
                    ty.invM.y,tx.invM.y,
                    tx.invM.x){
  #calculate fcd
  p1=-(1/2)*log.detD
  p2=ty.invM.y-t(2*t(param$betas)%*%tx.invM.y)
  for (i in 1:nphi.pot){
    p2[i]=p2[i]+t(param$betas)%*%tx.invM.x[,,i]%*%param$betas
  }
  fim=p1-(1/(2*param$sig2))*p2
  
  #normalize
  fim=fim-max(fim)
  fim1=exp(fim)
  fim1=fim1/sum(fim1)
  # plot(1:nphi.pot,unlist(fim1))
  
  #sample
  tmp=rmultinom(1,size=1,prob=unlist(fim1))
  ind=which(tmp==1)
  ind
}