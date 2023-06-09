rm(list=ls())
library('ggplot2')
library('mvtnorm')
set.seed(3)

#settings
nloc=2000
coord=data.frame(x=runif(nloc),
                 y=runif(nloc))
dist1=as.matrix(dist(coord))

#make spatial pattern for covariates
x1=coord$x*coord$y; x1=(x1-mean(x1))/sd(x1)
x2=coord$x/coord$y; x2=(x2-mean(x2))/sd(x2)
xmat=cbind(1,x1,x2)
cor(x1,x2)

# tmp=data.frame(x=coord$x,y=coord$y,res=x2)
# tmp.loess <- loess (res ~ x * y, data=tmp,
#                     degree = 2, span = 0.2)
# x <- seq (min (tmp$x), max (tmp$x), .05)
# y <- seq (min (tmp$y), max (tmp$y), .05)
# interpolated <- predict (tmp.loess, expand.grid (x = x, y = y))
# image (x= x, y= y, z = interpolated, asp = 1)

#parameters
phi=3
seq1=seq(from=0,to=max(dist1),length.out=100); plot(exp(-phi*seq1),type='l')
sig2=0.3
betas=c(0.5,0.3,0.1)

#simulate data
Mtilde=exp(-phi*dist1)
media=xmat%*%betas
resp=rmvnorm(1,media,sig2*Mtilde)

hist(resp)

#look at spatial distribution
err=t(resp)-media
tmp=data.frame(x=coord[,1],
               y=coord[,2],
               err=err,
               resp=t(resp))
range(tmp$err)
# plot(y~x,data=tmp)

#look at errors
tmp.loess <- loess (err ~ x * y, data=tmp, 
                     degree = 2, span = 0.2)
x <- seq (min (tmp$x), max (tmp$x), .05)
y <- seq (min (tmp$y), max (tmp$y), .05)
interpolated <- predict (tmp.loess, expand.grid (x = x, y = y))
image (x= x, y= y, z = interpolated, asp = 1)

#look at response
tmp.loess <- loess (resp ~ x * y, data=tmp, 
                    degree = 2, span = 0.2)
x <- seq (min (tmp$x), max (tmp$x), .05)
y <- seq (min (tmp$y), max (tmp$y), .05)
interpolated <- predict (tmp.loess, expand.grid (x = x, y = y))
image (x= x, y= y, z = interpolated, asp = 1)

#export info
fim=data.frame(x=coord$x,
               y=coord$y,
               cov1=xmat[,'x1'],
               cov2=xmat[,'x2'],
               resp=t(resp))
setwd('U:\\GIT_models\\NNGP\\simulat')
write.csv(fim[order(fim$x),],'1 sim data.csv',row.names=F)