library(mixtools)
library(readxl)
library(mpcmp)
setwd("/Users/derekyoung/Documents/CMP Research/JSTP Revision/Revision Code (Submitted)/5 - SIDS Data Analysis/")
source("Dongying_functions.R")
source("Derek_Functions.R")
sids2 <- read_excel("sids2.xlsx")


x <- sids2$BIR74
y <- sids2$SID74

# 1 component
out.mpcmp.sids1 <- glm.cmp(y~x)
out.poi.sids1 <- glm(y~x,family="poisson")
out.nb.sids1 <- glm.nb(y~x)

#MPCMP BIC
BIC(out.mpcmp.sids1) #542.2042

#Poisson BIC
BIC(out.poi.sids1) #637.1011

#NB BIC
BIC(out.nb.sids1) #542.9037





# 2 component
set.seed(1)
out.mpcmp.sids2 <- mix.mpcmp2(x=x,y=y,k=2,eps=1e-4)
#out.poi.sids2 <- poisregmixEM(y=y, x=x, epsilon = 1e-4)
#out.nb.sids2 <- nb.mixEMReg(y=y, x=x, eps=1e-4, k=2)
out.poi.sids2 <- poisregmixEM(y=y, x=x, epsilon = 1e-4, k=2, lambda=out.mpcmp.sids2$lambda, beta = out.mpcmp.sids2$betas)
out.nb.sids2 <- nb.mixEMReg(y=y, x=x, eps=1e-4, k=2, theta.t=out.mpcmp.sids2[1:3])
#out.mpcmp2.sids2 <- mix.mpcmp2(x=x,y=y,k=2,eps=1e-4,start=list(betas=out.poi2.sids2$beta,lambda=out.poi2.sids2$lambda,nus=c(1,1)))


#MPCMP BIC
-2*tail(out.mpcmp.sids2$all.ll,1)+7*log(100) #536.3535

#Poisson BIC
-2*out.poi.sids2$loglik+5*log(100) #537.0027

#NB BIC
-2*tail(out.nb.sids2[,2],1)+7*log(100) #538.0002






# 3 component
set.seed(1)
#out.mpcmp.sids3 <- mix.mpcmp2(x=x,y=y,k=3,eps=1e-4)
#out.poi.sids3 <- poisregmixEM(y=y, x=x, epsilon = 1e-4, k=3, lambda=out.mpcmp.sids3$lambda, beta = out.mpcmp.sids3$betas)
#out.nb.sids3 <- nb.mixEMReg(y=y, x=x, eps=1e-4, k=3, theta.t=out.mpcmp.sids3[1:3])
out.poi.sids3 <- poisregmixEM(y=y, x=x, epsilon = 1e-4,k=3)
out.nb.sids3 <- nb.mixEMReg(y=y, x=x, eps=1e-4, k=3)
out.mpcmp.sids3 <- mix.mpcmp2(x=x,y=y,k=2,eps=1e-4,start=list(betas=out.poi.sids3$beta,lambda=out.poi.sids3$lambda,nus=c(1,1,1)))


#MPCMP BIC
-2*tail(out.mpcmp.sids3$all.ll,1)+11*log(100) #554.7741

#Poisson BIC
-2*out.poi.sids3$loglik+8*log(100) #550.6334

#NB BIC
-2*tail(out.nb.sids3[,2],1)+11*log(100) #556.4209

#Replace with the ggplot code like what was used for Figure 3.
X <- 0:25000
plot(x,y,pch=19,main="SIDS Data in North Carolina (1974-78)",xlab="# of Live Births",ylab="# of SIDS Deaths")
points(x,y,col=apply(out.mpcmp.sids2$posteriors,1,which.max)+1,pch=19)
lines(X,exp(out.mpcmp.sids2$betas[1,1]+out.mpcmp.sids2$betas[2,1]*X),col=2)
lines(X,exp(out.mpcmp.sids2$betas[1,2]+out.mpcmp.sids2$betas[2,2]*X),col=3)




set.seed(1)
B <- 1000
poi.bs <- nb.bs <- cmp.bs <- NULL
n <- length(x)
for(i in 1:B){
  if(i%in%seq(50,B,by=50)) print(i)
  
  y1 <- y2 <- y3 <- rep(NA,n)
  
  z.poi <- rbinom(n,1,out.poi.sids2$lambda[1])
  z.nb <- rbinom(n,1,tail(out.nb.sids2,1)[3])
  z.cmp <- rbinom(n,1,out.mpcmp.sids2$lambda[1])
  
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

set.seed(10)
poi.bs <- nb.bs <- cmp.bs <- matrix(sample(1:n,n*B,replace=TRUE),nrow=n)

poi.bs.fits <- nb.bs.fits <- cmp.bs.fits <- vector("list",B)
set.seed(10)  
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
  save.image("sids_out.RData")
}




