#setwd("")
source("Prelim_Functions.R")

B <- 1000
n <- 100
beta0 <- -0.5
beta1 <- 1
p <- 0.15
X <- seq(0,1,length=n)
lambda.X <- exp(beta0+beta1*X)

set.seed(2)
ZIP.data <- sapply(1:B,function(i) rzipois(n,lambda=lambda.X,pstr0=p))
ZIB.data <- sapply(1:B,function(i) rzibinom(n,size=1,prob=inv.logit(beta0+beta1*X)))
ZINB.data <- sapply(1:B,function(i) rzinegbin(n,size=5,munb=lambda.X,pstr0=p))
ZIG.data <- sapply(1:B,function(i) rzigeom(n,prob=1/(1+lambda.X),pstr0=p))
ZIGP1.data <- sapply(1:B,function(i) sapply(1:n,function(j) rzigp(1,mu=lambda.X[j],phi=1-lambda.X[j]/(4*sqrt(exp(1))),omega=p))) #NB, phi=1-0.5*mu/mu=0.5
ZIGP2.data <- sapply(1:B,function(i) sapply(1:n,function(j) rzigp(1,mu=lambda.X[j],phi=1+2*lambda.X[j],omega=p)))
ZICMP1.data <- sapply(1:B,function(i) rzicmp(n,lambda.X,nu=5,p=p))
ZICMP2.data <- sapply(1:B,function(i) rzicmp(n,lambda.X,nu=0.3,p=p))
ZIsCMP1.data <- sapply(1:B,function(i) sapply(1:n,function(j) r.zi.scompoisson(1,NuOfVar=2,lambda=lambda.X[j],nu=5,p=p)))
ZIsCMP2.data <- sapply(1:B,function(i) sapply(1:n,function(j) r.zi.scompoisson(1,NuOfVar=2,lambda=lambda.X[j],nu=0.3,p=p)))
ZIsCMP3.data <- sapply(1:B,function(i) sapply(1:n,function(j) r.zi.scompoisson(1,NuOfVar=3,lambda=lambda.X[j],nu=5,p=p)))
ZIsCMP4.data <- sapply(1:B,function(i) sapply(1:n,function(j) r.zi.scompoisson(1,NuOfVar=3,lambda=lambda.X[j],nu=0.3,p=p)))

save.image("/Users/derekyoung/Documents/Derek's MacBook/Documents/ZIsCOM Poisson Regression/Manuscript/codes/Second Sim_Results/n100.RData")
