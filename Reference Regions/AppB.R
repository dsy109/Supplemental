###############################################################
###############################################################
### Code for "Appendix B: Coverage Study for Rectangular 
### Central Tolerance Regions Using the Parametric Bootstrap in
### Algorithm 5" in "A Review of Statistical Reference Regions 
### in Laboratory Medicine: Theory and Computation" by Thomas 
### Mathew and Derek S. Young
###############################################################
###############################################################

library(mvtnorm)
library(nlme)
library(ggplot2)

alpha <- 0.05
P <- 0.90

BS.fn <- function(x,B=1000,P,alpha){
  x.bar <- apply(x,2,mean)
  p <- length(x.bar)
  n <- nrow(x)
  S <- cov(x)
  s.ii <- diag(S)
  R <- cov2cor(S)
  f <- function(c.rho,x.bar,S,P) pmvnorm(lower=x.bar-c.rho*sqrt(diag(S)),upper=x.bar+c.rho*sqrt(diag(S)),mean=x.bar,sigma=S)-P 
  c.rho <- uniroot(f,interval=c(0,100),x.bar=x.bar,S=S,P=P)$root
  x.bar.bs <- rmvnorm(B,sigma=R/n)
  S.bar.bs <- rWishart(B,n-1,R)/(n-1)
  W.b <- sapply(1:B,function(i) max((abs(x.bar.bs[i,])+c.rho)/sqrt(diag(S.bar.bs[,,i]))) )
  c1.rho <- quantile(W.b,1-alpha)
  c1.rho
}

f <- function(c.rho,x.bar,S,P) pmvnorm(lower=x.bar-c.rho*sqrt(diag(S)),upper=x.bar+c.rho*sqrt(diag(S)),mean=x.bar,sigma=S)-P 

M <- 5000
n <- c(15,50,100,500,1500)


#### d=2

rho <- seq(-.9,.9,by=.1)
mu <- c(0,0)

all.CPs_2 <- matrix(NA,ncol=length(n),nrow=length(rho))

set.seed(1)
for(j in 1:length(n)){
  for(k in 1:length(rho)){
    Sigma <- rbind(c(1,rho[k]),c(rho[k],1))
    true.c.rho <- uniroot(f,interval=c(0,20),x.bar=mu,S=Sigma,P=P)$root
    true.region <- cbind(mu-true.c.rho*sqrt(diag(Sigma)),mu+true.c.rho*sqrt(diag(Sigma)))
    MC.regions <- vector("list",M)
    for(i in 1:M){
      x <- rmvnorm(n[j],mean=mu,sigma=Sigma)
      x.bar <- apply(x,2,mean)
      S <- cov(x)
      c1.rho <- BS.fn(x,B=1000,P=P,alpha=alpha)
      MC.regions[[i]] <- cbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
    }
    
    all.CPs_2[k,j] <- mean(sapply(1:M,function(i) all(c(MC.regions[[i]][,1]<=true.region[,1],MC.regions[[i]][,2]>=true.region[,2]))))
    print(c(j,k))
  }
  print(all.CPs_2)
}
all.CPs_2a <- data.frame(rho=rep(rho,5),CP=c(all.CPs_2),n=as.factor(rep(n,each=19)))

ggplot(all.CPs_2a, aes(rho, CP)) + geom_point(aes(rho, CP,col=n),size=3)  + geom_hline(yintercept=0.95) + 
  labs(x = expression(rho), y = "Coverage Probability") +
  geom_line(aes(rho, CP,col=n),lty=1) + theme(text = element_text(size = 20)) + ggtitle("Coverage Probabilities for Rectangular Central Tolerance Regions \n(Bivariate Normal)")


#### d=3

mu <- c(0,0,0)
Sigma <- rbind(c(1,-.5,.1),c(-.5,1,-.8),c(.1,-.8,1))

all.CPs_3 <- rep(NA,length(n))

set.seed(1)
for(j in 1:length(n)){
  true.c.rho <- uniroot(f,interval=c(0,20),x.bar=mu,S=Sigma,P=P)$root
  true.region <- cbind(mu-true.c.rho*sqrt(diag(Sigma)),mu+true.c.rho*sqrt(diag(Sigma)))
  MC.regions <- vector("list",M)
  for(i in 1:M){
    x <- rmvnorm(n[j],mean=mu,sigma=Sigma)
    x.bar <- apply(x,2,mean)
    S <- cov(x)
    c1.rho <- BS.fn(x,B=1000,P=P,alpha=alpha)
    MC.regions[[i]] <- cbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
  }
  
  all.CPs_3[j] <- mean(sapply(1:M,function(i) all(c(MC.regions[[i]][,1]<=true.region[,1],MC.regions[[i]][,2]>=true.region[,2]))))
  print(all.CPs_3)
}



#### d=6

mu <- rep(0,6)
Sigma <- toeplitz(c(1,-.4,.2,0,0,0))

all.CPs_6 <- rep(NA,length(n))

set.seed(1)
for(j in 1:length(n)){
  true.c.rho <- uniroot(f,interval=c(0,20),x.bar=mu,S=Sigma,P=P)$root
  true.region <- cbind(mu-true.c.rho*sqrt(diag(Sigma)),mu+true.c.rho*sqrt(diag(Sigma)))
  MC.regions <- vector("list",M)
  for(i in 1:M){
    x <- rmvnorm(n[j],mean=mu,sigma=Sigma)
    x.bar <- apply(x,2,mean)
    S <- cov(x)
    c1.rho <- BS.fn(x,B=1000,P=P,alpha=alpha)
    MC.regions[[i]] <- cbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
  }
  
  all.CPs_6[j] <- mean(sapply(1:M,function(i) all(c(MC.regions[[i]][,1]<=true.region[,1],MC.regions[[i]][,2]>=true.region[,2]))))
  print(all.CPs_6)
}


#### d=10

mu <- rep(0,10)
cor1 <- corAR1(0.5,form=~1)
cor1. <- Initialize(cor1,data=data.frame(lag=1:10))
Sigma <- corMatrix(cor1.)

all.CPs_10 <- rep(NA,length(n))

set.seed(1)
for(j in 1:length(n)){
  true.c.rho <- uniroot(f,interval=c(0,20),x.bar=mu,S=Sigma,P=P)$root
  true.region <- cbind(mu-true.c.rho*sqrt(diag(Sigma)),mu+true.c.rho*sqrt(diag(Sigma)))
  MC.regions <- vector("list",M)
  set.seed(1)
  for(i in 1:M){
    x <- rmvnorm(n[j],mean=mu,sigma=Sigma)
    x.bar <- apply(x,2,mean)
    S <- cov(x)
    c1.rho <- BS.fn(x,B=1000,P=P,alpha=alpha)
    MC.regions[[i]] <- cbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
  }
  
  all.CPs_10[j] <- mean(sapply(1:M,function(i) all(c(MC.regions[[i]][,1]<=true.region[,1],MC.regions[[i]][,2]>=true.region[,2]))))
  print(all.CPs_10)
}

save(all.CPs_2a,all.CPs_3,all.CPs_6,all.CPs_10,file="CP_results.RData")




