library(mixtools)
library(mpcmp)
source("/project/dyo227_uksr/JSTP/Dongying_functions.R")
source("/project/dyo227_uksr/JSTP/Derek_Functions.R")
sids2 <- read.delim("sids.txt")
setwd("/scratch/dyo227/")


x <- sids2$BIR74
y <- sids2$SID74


# 2 component
set.seed(1)
out.mpcmp.sids2 <- mix.mpcmp2(x=x,y=y,k=2,eps=1e-4)
out.poi.sids2 <- poisregmixEM(y=y, x=x, epsilon = 1e-4, k=2, lambda=out.mpcmp.sids2$lambda, beta = out.mpcmp.sids2$betas)
out.nb.sids2 <- nb.mixEMReg(y=y, x=x, eps=1e-4, k=2, theta.t=out.mpcmp.sids2[1:3])


beta.t <- matrix(tail(out.nb.sids2,1)[5:8],ncol=2)
pi.t <- tail(out.nb.sids2,1)[3:4]
phi.t <- tail(out.nb.sids2,1)[9:10]
mu.t <- exp(cbind(1,x) %*% beta.t)
y.k <- matrix(NA,nrow=length(x), ncol=2)
for (i in 1:length(x)){
  for (j in 1:2) {
    y.k[i,j] <- dnbinom(y[i],mu=mu.t[i,j],size=phi.t[j])
  }
}
nb.post <- t(t(y.k)*pi.t) / rowSums(t(t(y.k)*pi.t))




set.seed(1)
B <- 1000
poi.bs <- nb.bs <- cmp.bs <- NULL
n <- length(x)
for(i in 1:B){
  if(i%in%seq(50,B,by=50)) print(i)
  
  y1 <- y2 <- y3 <- rep(NA,n)
  
  z.poi <- rbinom(n,1,round(out.poi.sids2$posterior[,1],4))
  z.nb <- rbinom(n,1,round(nb.post[,1],4))
  z.cmp <- rbinom(n,1,round(out.mpcmp.sids2$posterior[,1]))
  
  y1[z.poi==1] <- rpois(sum(z.poi==1),exp(out.poi.sids2$beta[1,1]+out.poi.sids2$beta[2,1]*x[z.poi==1]))
  y1[z.poi==0] <- rpois(sum(z.poi==0),exp(out.poi.sids2$beta[1,2]+out.poi.sids2$beta[2,2]*x[z.poi==0]))
  y2[z.nb==1] <- rnegbin(sum(z.nb==1),exp(tail(out.nb.sids2,1)[5]+tail(out.nb.sids2,1)[6]*x[z.nb==1]),theta=tail(out.nb.sids2,1)[9])
  y2[z.nb==0] <- rnegbin(sum(z.nb==0),exp(tail(out.nb.sids2,1)[7]+tail(out.nb.sids2,1)[8]*x[z.nb==0]),theta=tail(out.nb.sids2,1)[10])
  y3[z.cmp==1] <- rcomp(sum(z.cmp==1),exp(out.mpcmp.sids2$beta[1,1]+out.mpcmp.sids2$beta[2,1]*x[z.cmp==1]),nu=out.mpcmp.sids2$nus[1])
  y3[z.cmp==0] <- rcomp(sum(z.cmp==0),exp(out.mpcmp.sids2$beta[1,2]+out.mpcmp.sids2$beta[2,2]*x[z.cmp==0]),nu=out.mpcmp.sids2$nus[2])
  
  poi.bs <- cbind(poi.bs,y1)
  nb.bs <- cbind(nb.bs,y2)
  cmp.bs <- cbind(cmp.bs,y3)
  
}

poi.bs.fits <- nb.bs.fits <- cmp.bs.fits <- vector("list",B)
set.seed(1)  
n <- length(x)
for(i in 1:B){
  print(paste("Bootstrap",i))
  poi.bs.fits[[i]] <- poisregmixEM(x=x,y=poi.bs[,i], lambda=out.poi.sids2$lambda, beta=out.poi.sids2$beta,  epsilon = 1e-4, k=2)
  temp.nb <- try(tail(nb.mixEMReg(x=x,y=nb.bs[,i], theta.t=list(tail(out.nb.sids2,1)[9:10],matrix(tail(out.nb.sids2,1)[5:8],ncol=2),tail(out.nb.sids2,1)[9:10]),eps = 1e-4, k=2),1),silent=TRUE)
  while(class(temp.nb)[1]=="try-error"){
    #    nb.bs.fits[[i]] <- tail(nb.mixEMReg(x=x,y=nb.bs[,i],eps = 1e-4, k=2),1)
    y2 <- rep(NA,n)
    z.nb <- rbinom(n,1,tail(out.nb.sids2,1)[3])
    y2[z.nb==1] <- rnegbin(sum(z.nb==1),exp(tail(out.nb.sids2,1)[5]+tail(out.nb.sids2,1)[6]*x[z.nb==1]),theta=tail(out.nb.sids2,1)[9])
    y2[z.nb==0] <- rnegbin(sum(z.nb==0),exp(tail(out.nb.sids2,1)[7]+tail(out.nb.sids2,1)[8]*x[z.nb==0]),theta=tail(out.nb.sids2,1)[10])
    nb.bs[,i] <- y2
    temp.nb <- try(tail(nb.mixEMReg(x=x,y=nb.bs[,i], theta.t=list(tail(out.nb.sids2,1)[9:10],matrix(tail(out.nb.sids2,1)[5:8],ncol=2),tail(out.nb.sids2,1)[9:10]),eps = 1e-4, k=2),1),silent=TRUE)
  }
  nb.bs.fits[[i]] <- temp.nb
  temp.cmp <- try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],start=list(nus=out.mpcmp.sids2$nu,lambda=out.mpcmp.sids2$lambda,beta=out.mpcmp.sids2$beta),k=2,eps=1e-4),timeout=30,onTimeout = "silent"),silent=TRUE)
  while(class(temp.cmp)[1]=="try-error"|is.null(temp.cmp)){
    #    try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],k=2,eps=1e-4),timeout=40,onTimeout = "silent"),silent=TRUE)
    y3 <- rep(NA,n)
    z.cmp <- rbinom(n,1,out.mpcmp.sids2$lambda[1])
    y3[z.cmp==1] <- rcomp(sum(z.cmp==1),exp(out.mpcmp.sids2$beta[1,1]+out.mpcmp.sids2$beta[2,1]*x[z.cmp==1]),nu=out.mpcmp.sids2$nus[1])
    y3[z.cmp==0] <- rcomp(sum(z.cmp==0),exp(out.mpcmp.sids2$beta[1,2]+out.mpcmp.sids2$beta[2,2]*x[z.cmp==0]),nu=out.mpcmp.sids2$nus[2])
    cmp.bs[,i] <- y3
    temp.cmp <- try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],start=list(nus=out.mpcmp.sids2$nu,lambda=out.mpcmp.sids2$lambda,beta=out.mpcmp.sids2$beta),k=2,eps=1e-4),timeout=30,onTimeout = "silent"),silent=TRUE)
  }
  cmp.bs.fits[[i]] <- temp.cmp
}






set.seed(1)
B <- 1000
n <- length(x)
poi.bs <- nb.bs <- cmp.bs <- matrix(sample(1:n,n*B,replace=TRUE),nrow=n)

poi.bs.fits <- nb.bs.fits <- cmp.bs.fits <- vector("list",B)
set.seed(SS)  
n <- length(x)
for(i in 1:B){
  print(paste("Bootstrap",i))
  poi.bs.fits[[i]] <- poisregmixEM(x=x[poi.bs[,i]],y=y[poi.bs[,i]], lambda=out.poi.sids2$lambda, beta=out.poi.sids2$beta,  epsilon = 1e-3, k=2)
  temp.nb <- try(tail(nb.mixEMReg(x=x[nb.bs[,i]],y=y[nb.bs[,i]], theta.t=list(tail(out.nb.sids2,1)[9:10],matrix(tail(out.nb.sids2,1)[5:8],ncol=2),tail(out.nb.sids2,1)[9:10]),eps = 1e-3, k=2),1),silent=TRUE)
  while(class(temp.nb)[1]=="try-error"){
    #    nb.bs.fits[[i]] <- tail(nb.mixEMReg(x=x,y=nb.bs[,i],eps = 1e-4, k=2),1)
    nb.bs[,i] <- sample(1:n,n,replace=TRUE)
    temp.nb <- try(tail(nb.mixEMReg(x=x[nb.bs[,i]],y=y[nb.bs[,i]], theta.t=list(tail(out.nb.sids2,1)[9:10],matrix(tail(out.nb.sids2,1)[5:8],ncol=2),tail(out.nb.sids2,1)[9:10]),eps = 1e-3, k=2),1),silent=TRUE)
  }
  nb.bs.fits[[i]] <- temp.nb
  temp.cmp <- try(withTimeout(mix.mpcmp2(x=x[cmp.bs[,i]],y=y[cmp.bs[,i]],start=list(nus=out.mpcmp.sids2$nu,lambda=out.mpcmp.sids2$lambda,beta=out.mpcmp.sids2$beta),k=2,eps=1e-3),timeout=60,onTimeout = "silent"),silent=TRUE)
  while(class(temp.cmp)[1]=="try-error"|is.null(temp.cmp)){
    #    try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],k=2,eps=1e-4),timeout=40,onTimeout = "silent"),silent=TRUE)
    cmp.bs[,i] <- sample(1:n,n,replace=TRUE)
    temp.cmp <- try(withTimeout(mix.mpcmp2(x=x[cmp.bs[,i]],y=y[cmp.bs[,i]],start=list(nus=out.mpcmp.sids2$nu,lambda=out.mpcmp.sids2$lambda,beta=out.mpcmp.sids2$beta),k=2,eps=1e-3),timeout=60,onTimeout = "silent"),silent=TRUE)
  }
  cmp.bs.fits[[i]] <- temp.cmp
  if(i %in% seq(10,B,by=10)) save.image("sids_results.RData")
}



