library(mixreg)
library(mixtools)
library(mpcmp)
setwd("/Users/derekyoung/Documents/CMP Research/JSTP Revision/Revision Code (Submitted)/4 - Aphids Data Analysis/")
source("Dongying_functions.R")
source("Derek_Functions.R")

x=aphids$aphRel
y=aphids$plntsInf

# 1 component
out.mpcmp.aphids1 <- glm.cmp(plntsInf~aphRel,data=aphids)
out.poi.aphids1 <- glm(plntsInf~aphRel,data=aphids,family="poisson")
out.nb.aphids1 <- glm.nb(plntsInf~aphRel,data=aphids)

#MPCMP BIC
BIC(out.mpcmp.aphids1) #283.5356

#Poisson BIC
BIC(out.poi.aphids1) #403.3398

#NB BIC
BIC(out.nb.aphids1) #283.9228



# 2 components
set.seed(1)
out.mpcmp.aphids2 <- mix.mpcmp2(x=aphids$aphRel,y=aphids$plntsInf,start=list(nus=c(1,1),lambda=c(.5,.5)),k=2,eps=1e-5)
out.poi.aphids2 <- poisregmixEM(x=aphids$aphRel,y=aphids$plntsInf, epsilon = 1e-5, k=2)
out.nb.aphids2 <- nb.mixEMReg(x=aphids$aphRel,y=aphids$plntsInf, eps=1e-5, k=2)

#MPCMP BIC
-2*tail(out.mpcmp.aphids2$all.ll,1)+7*log(51) #282.0920

#Poisson BIC
-2*out.poi.aphids2$loglik+5*log(51) #275.0845

#NB BIC
-2*tail(out.nb.aphids2[,2],1)+7*log(51) #282.4683


# 3 components
set.seed(10)
out.mpcmp.aphids3 <- mix.mpcmp2(x=aphids$aphRel,y=aphids$plntsInf,maxit=14,k=3,eps=1e-4,start=list(betas=out.poi.aphids3$beta,lambda=out.poi.aphids3$lambda,nus=c(1,1,1)))
out.poi.aphids3 <- poisregmixEM(x=aphids$aphRel,y=aphids$plntsInf, epsilon = 1e-4, k=3)
out.nb.aphids3 <- nb.mixEMReg(x=aphids$aphRel,y=aphids$plntsInf, eps=1e-4, k=3)

#MPCMP BIC
-2*tail(out.mpcmp.aphids3$all.ll,1)+11*log(51) #294.1541

#Poisson BIC
-2*out.poi.aphids3$loglik+8*log(51) #285.0741

#NB BIC
-2*tail(out.nb.aphids3[,2],1)+11*log(51) #298.1957



#Replace with the ggplot code used for Figure 3.
plot(x=aphids$aphRel,y=aphids$plntsInf,pch=19,main="Aphids Data",xlab="# of Aphids Release",ylab="# of Infected Plants")
points(x=aphids$aphRel,y=aphids$plntsInf,col=apply(out.mpcmp.aphids2$posteriors,1,which.max)+1,pch=19)
X <- 0:350
lines(X,exp(out.mpcmp.aphids2$betas[1,1]+out.mpcmp.aphids2$betas[2,1]*X),col=2)
lines(X,exp(out.mpcmp.aphids2$betas[1,2]+out.mpcmp.aphids2$betas[2,2]*X),col=3)


set.seed(1)
B <- 1000
poi.bs <- nb.bs <- cmp.bs <- NULL
n <- length(x)
for(i in 1:B){
  if(i%in%seq(50,B,by=50)) print(i)
  
  y1 <- y2 <- y3 <- rep(NA,n)

  z.poi <- rbinom(n,1,out.poi.aphids2$lambda[1])
  z.nb <- rbinom(n,1,tail(out.nb.aphids2,1)[3])
  z.cmp <- rbinom(n,1,out.mpcmp.aphids2$lambda[1])
  
  y1[z.poi==1] <- rpois(sum(z.poi==1),exp(out.poi.aphids2$beta[1,1]+out.poi.aphids2$beta[2,1]*x[z.poi==1]))
  y1[z.poi==0] <- rpois(sum(z.poi==0),exp(out.poi.aphids2$beta[1,2]+out.poi.aphids2$beta[2,2]*x[z.poi==0]))
  y2[z.nb==1] <- rnegbin(sum(z.nb==1),exp(tail(out.nb.aphids2,1)[5]+tail(out.nb.aphids2,1)[6]*x[z.nb==1]),theta=tail(out.nb.aphids2,1)[9])
  y2[z.nb==0] <- rnegbin(sum(z.nb==0),exp(tail(out.nb.aphids2,1)[7]+tail(out.nb.aphids2,1)[8]*x[z.nb==0]),theta=tail(out.nb.aphids2,1)[10])
  y3[z.cmp==1] <- rcomp(sum(z.cmp==1),exp(out.mpcmp.aphids2$beta[1,1]+out.mpcmp.aphids2$beta[2,1]*x[z.cmp==1]),nu=out.mpcmp.aphids2$nus[1])
  y3[z.cmp==0] <- rcomp(sum(z.cmp==0),exp(out.mpcmp.aphids2$beta[1,2]+out.mpcmp.aphids2$beta[2,2]*x[z.cmp==0]),nu=out.mpcmp.aphids2$nus[2])
  
  poi.bs <- cbind(poi.bs,y1)
  nb.bs <- cbind(nb.bs,y2)
  cmp.bs <- cbind(cmp.bs,y3)
  
}

poi.bs.fits <- nb.bs.fits <- cmp.bs.fits <- vector("list",B)
set.seed(1)  
n <- length(x)
for(i in 1:B){
  print(paste("Bootstrap",i))
  poi.bs.fits[[i]] <- poisregmixEM(x=x,y=poi.bs[,i], lambda=out.poi.aphids2$lambda, beta=out.poi.aphids2$beta,  epsilon = 1e-4, k=2)
  temp.nb <- try(tail(nb.mixEMReg(x=x,y=nb.bs[,i], theta.t=list(tail(out.nb.aphids2,1)[9:10],matrix(tail(out.nb.aphids2,1)[5:8],ncol=2),tail(out.nb.aphids2,1)[9:10]),eps = 1e-4, k=2),1),silent=TRUE)
  while(class(temp.nb)[1]=="try-error"){
#    nb.bs.fits[[i]] <- tail(nb.mixEMReg(x=x,y=nb.bs[,i],eps = 1e-4, k=2),1)
    y2 <- rep(NA,n)
    z.nb <- rbinom(n,1,tail(out.nb.aphids2,1)[3])
    y2[z.nb==1] <- rnegbin(sum(z.nb==1),exp(tail(out.nb.aphids2,1)[5]+tail(out.nb.aphids2,1)[6]*x[z.nb==1]),theta=tail(out.nb.aphids2,1)[9])
    y2[z.nb==0] <- rnegbin(sum(z.nb==0),exp(tail(out.nb.aphids2,1)[7]+tail(out.nb.aphids2,1)[8]*x[z.nb==0]),theta=tail(out.nb.aphids2,1)[10])
    nb.bs[,i] <- y2
    temp.nb <- try(tail(nb.mixEMReg(x=x,y=nb.bs[,i], theta.t=list(tail(out.nb.aphids2,1)[9:10],matrix(tail(out.nb.aphids2,1)[5:8],ncol=2),tail(out.nb.aphids2,1)[9:10]),eps = 1e-4, k=2),1),silent=TRUE)
  }
  nb.bs.fits[[i]] <- temp.nb
  temp.cmp <- try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],start=list(nus=out.mpcmp.aphids2$nu,lambda=out.mpcmp.aphids2$lambda,beta=out.mpcmp.aphids2$beta),k=2,eps=1e-4),timeout=30,onTimeout = "silent"),silent=TRUE)
  while(class(temp.cmp)[1]=="try-error"|is.null(temp.cmp)){
#    try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],k=2,eps=1e-4),timeout=40,onTimeout = "silent"),silent=TRUE)
    y3 <- rep(NA,n)
    z.cmp <- rbinom(n,1,out.mpcmp.aphids2$lambda[1])
    y3[z.cmp==1] <- rcomp(sum(z.cmp==1),exp(out.mpcmp.aphids2$beta[1,1]+out.mpcmp.aphids2$beta[2,1]*x[z.cmp==1]),nu=out.mpcmp.aphids2$nus[1])
    y3[z.cmp==0] <- rcomp(sum(z.cmp==0),exp(out.mpcmp.aphids2$beta[1,2]+out.mpcmp.aphids2$beta[2,2]*x[z.cmp==0]),nu=out.mpcmp.aphids2$nus[2])
    cmp.bs[,i] <- y3
    temp.cmp <- try(withTimeout(mix.mpcmp2(x=x,y=cmp.bs[,i],start=list(nus=out.mpcmp.aphids2$nu,lambda=out.mpcmp.aphids2$lambda,beta=out.mpcmp.aphids2$beta),k=2,eps=1e-4),timeout=30,onTimeout = "silent"),silent=TRUE)
  }
  cmp.bs.fits[[i]] <- temp.cmp
}

#poi.bs.fits <- lapply(1:B, function(i) poisregmixEM(x=x,y=poi.bs[,i], lambda=out.poi.aphids2$lambda, beta=out.poi.aphids2$beta,  epsilon = 1e-4, k=2))
#nb.bs.fits <- lapply(1:B, function(i) nb.mixEMReg(x=x,y=nb.bs[,i], theta.t=list(tail(out.nb.aphids2,1)[9:10],matrix(tail(out.nb.aphids2,1)[5:8],ncol=2),tail(out.nb.aphids2,1)[9:10]),eps = 1e-4, k=2))
#cmp.bs.fits <- lapply(1:B, function(i) nb.mixEMReg(x=x,y=nb.bs[,i], theta.t=list(tail(out.nb.aphids2,1)[9:10],matrix(tail(out.nb.aphids2,1)[5:8],ncol=2),tail(out.nb.aphids2,1)[9:10]),eps = 1e-4, k=2))

#cmp.bs.fits <- lapply(1:B, function(i) mix.mpcmp2(x=x,y=nb.bs[,i],start=list(nus=out.mpcmp.aphids2$nu,lambda=out.mpcmp.aphids2$lambda,beta=out.mpcmp.aphids2$beta),k=2,eps=1e-4))

save.image("aphids_results2.RData")

min.max <- function(bs.beta,est.beta){
  ind.1 <- sum((bs.beta-est.beta)^2)
  ind.2 <- sum((bs.beta[,2:1]-est.beta)^2)
  if(ind.1<ind.2) 1:2 else 2:1
}

poi.ident <- t(sapply(1:B, function(i) min.max(poi.bs.fits[[i]]$beta,est.beta=out.poi.aphids2$beta) ) )
nb.ident <- t(sapply(1:B, function(i) min.max(matrix(nb.bs.fits[[i]][5:8],nrow=2),est.beta=matrix(tail(out.nb.aphids2,1)[5:8],nrow=2) ) ) )
cmp.ident <- t(sapply(1:B, function(i) min.max(cmp.bs.fits[[i]]$beta,est.beta=out.mpcmp.aphids2$beta) ) )

all.poi.par <- all.nb.par <- all.cmp.par <- NULL
for(i in 1:B){
  all.poi.par <- rbind(all.poi.par,c(poi.bs.fits[[i]]$lambda[poi.ident[i,]],c(poi.bs.fits[[i]]$beta[,poi.ident[i,]])) )
  all.nb.par <- rbind(all.nb.par, c(nb.bs.fits[[i]][3:4][nb.ident[i,]],c(matrix(nb.bs.fits[[i]][5:8],nrow=2)[,nb.ident[i,]]),nb.bs.fits[[i]][9:10][nb.ident[i,]]   ))
  all.cmp.par <- rbind(all.cmp.par,c(cmp.bs.fits[[i]]$lambda[cmp.ident[i,]],c(cmp.bs.fits[[i]]$beta[,cmp.ident[i,]]),cmp.bs.fits[[i]]$nu[cmp.ident[i,]]) )
}

apply(all.poi.par,2,sd)
apply(all.nb.par,2,sd)
apply(all.cmp.par,2,sd)

ttt=c(which(all.cmp.par[,7]>quantile(all.cmp.par[,7],.98)),which(all.cmp.par[,8]>quantile(all.cmp.par[,8],.98)))
apply(all.cmp.par[-ttt,],2,sd)



