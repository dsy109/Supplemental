#setwd("")
source("Prelim_Functions.R")
shark <- read.table("SharksGBRMP.txt",header=T)

#1:23
#shark.type <- apply(shark[,17:23],1,sum)

zip.rqres <- function(y, x, z, eps = 1e-6, b0, b1, a0, a1)
{
  n <- length(y)
  FL <- pzipois(y - eps, lambda=exp(b0+b1*x),pstr0=1/(1+exp(-(a0+z*a1))))
  FU <- pzipois(y, lambda=exp(b0+b1*x),pstr0=1/(1+exp(-(a0+z*a1))))
  u <- runif(n, min = FL, max = FU)
  qres <- qnorm(u)
  return(qres)
}

nb.rqres <- function(y, x, eps = 1e-6, b0, b1, theta)
{
  n <- length(y)
  FL <- pnbinom(y - eps, mu=exp(b0+b1*x),size=theta)
  FU <- pnbinom(y, mu=exp(b0+b1*x),size=theta)
  u <- runif(n, min = FL, max = FU)
  qres <- qnorm(u)
  return(qres)
}

zig.rqres <- function(y, x, z, eps = 1e-6, b0, b1, a0, a1)
{
  n <- length(y)
  FL <- pzinegbin(y - eps, munb=exp(b0+b1*x),pstr0=1/(1+exp(-(a0+z*a1))),size=1)
  FU <- pzinegbin(y, munb=exp(b0+b1*x),pstr0=1/(1+exp(-(a0+z*a1))),size=1)
  u <- runif(n, min = FL, max = FU)
  qres <- qnorm(u)
  return(qres)
}


###################################################################
###################################################################
### SNS Analysis
###################################################################
###################################################################

shark.type <- shark[,15]
n <- length(shark.type)
table(shark.type)

P.shark <- glm(shark.type~shark$SoakTime,family=poisson)
NB.shark <- glm.nb(shark.type~shark$SoakTime)

ZIP.shark <- zeroinfl(shark.type~shark$SoakTime,dist="poisson")
ZIG.shark <- zeroinfl(shark.type~shark$SoakTime,dist="geometric")
ZINB.shark <- zeroinfl(shark.type~shark$SoakTime,dist="negbin")

ZIGP.shark <- est.zigp2(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)

ZICMP.shark <- fit.zicmp.reg2(y=shark.type,X=cbind(1,shark$SoakTime),S=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta.init=c(-4.680621250,0.004074857),gamma.init=1.1,zeta.init=c(-0.914944825,0.003777524),max=250,optim.control=list(maxit=2000))

ZIsCMP2.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=c(-10.1901,0.0435748),gamma0=10.9283,xi0=c(2.68428,0.00172647),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

ZIsCMP3.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZIsCMP2.shark$Parameters[1:2,1],gamma0=log(ZIsCMP2.shark$Parameters[3,1]),xi0=ZIsCMP2.shark$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

P.ll <- logLik(P.shark)
NB.ll <- logLik(NB.shark)
ZIP.ll <- logLik(ZIP.shark)
ZIG.ll <- logLik(ZIG.shark)
ZINB.ll <- logLik(ZINB.shark)
ZIGP.ll <- tail(ZIGP.shark,1)
ZICMP.ll <- ZICMP.shark$loglik
ZIsCMP2.ll <- ZIsCMP2.shark[[2]]
ZIsCMP3.ll <- ZIsCMP3.shark[[2]]

P.pars <- length(coef(P.shark))
NB.pars <- P.pars + 1
ZIP.pars <- length(coef(ZIP.shark))
ZIG.pars <- ZIP.pars
ZINB.pars <- ZIP.pars+1
ZIGP.pars <- length(ZIGP.shark)-1
ZICMP.pars <- ZIGP.pars
ZIsCMP2.pars <- ZIGP.pars
ZIsCMP3.pars <- ZIGP.pars

P.AIC <- -2*P.ll+2*P.pars
NB.AIC <- -2*NB.ll+2*NB.pars
ZIP.AIC <- -2*ZIP.ll+2*ZIP.pars
ZIG.AIC <- -2*ZIG.ll+2*ZIG.pars
ZINB.AIC <- -2*ZINB.ll+2*ZINB.pars
ZIGP.AIC <- -2*ZIGP.ll+2*ZIGP.pars
ZICMP.AIC <- -2*ZICMP.ll+2*ZICMP.pars
ZIsCMP2.AIC <- -2*ZIsCMP2.ll+2*ZIsCMP2.pars
ZIsCMP3.AIC <- -2*ZIsCMP3.ll+2*ZIsCMP3.pars

P.BIC <- -2*P.ll+log(n)*P.pars
NB.BIC <- -2*NB.ll+log(n)*NB.pars
ZIP.BIC <- -2*ZIP.ll+log(n)*ZIP.pars
ZIG.BIC <- -2*ZIG.ll+log(n)*ZIG.pars
ZINB.BIC <- -2*ZINB.ll+log(n)*ZINB.pars
ZIGP.BIC <- -2*ZIGP.ll+log(n)*ZIGP.pars
ZICMP.BIC <- -2*ZICMP.ll+log(n)*ZICMP.pars
ZIsCMP2.BIC <- -2*ZIsCMP2.ll+log(n)*ZIsCMP2.pars
ZIsCMP3.BIC <- -2*ZIsCMP3.ll+log(n)*ZIsCMP3.pars

R2.P <- 1-(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(P.ll)+P.pars+.5)/(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(logLik(glm(shark.type~1,family=poisson))))
R2.NB <- 1-(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(NB.ll)+NB.pars+.5)/(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(logLik(glm.nb(shark.type~1))))

d.zip <- function(x, lambda , pi) {
  ifelse(x == 0, pi + (1-pi)*dpois(0,lambda = lambda), (1-pi)*dpois(x,lambda = lambda))
}
ll.zip.sat <- sum(log(d.zip(shark.type,shark.type,as.numeric(shark.type == 0))))

d.zinb <- function(x,size,mu,pi) {
  ifelse(x == 0, pi + (1-pi)*dnbinom(0,size = size, mu = mu), (1-pi)*dnbinom(x,size = size, mu = mu))
}
ll.zinb.sat <- sum(log(d.zinb(shark.type,size = NB.shark$theta , mu = shark.type, pi = as.numeric(shark.type == 0))))

d.geom <- function(x,prob,pi) {
  ifelse(x == 0, pi + (1-pi)*dgeom(0,prob = prob), (1-pi)*dgeom(x,prob = prob))
}
ll.zig.sat <- sum(log(d.geom(shark.type,prob=1/(1+mean(shark.type)), pi = as.numeric(shark.type == 0))))

d.genpoi <- function(x,mu,phi,pi) {
  ifelse(x == 0, pi + (1-pi)*dzigp(0,mu=mu,phi=phi,omega=0), (1-pi)*dzigp(x,mu=mu,phi=phi,omega=0))
}
ll.zigp.sat <- sum(sapply(1:n,function(i) log(d.genpoi(shark.type[i],mu=shark.type[i], phi=1+exp(ZIGP.shark[3]), pi = as.numeric(shark.type[i] == 0)))))

d.cmp <- function(x,lambda,nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dcom(x,lambda=lambda,nu=nu))
}
ll.zicmp.sat <- sum(sapply(1:n,function(i) log(d.cmp(shark.type[i],lambda=shark.type[i], nu=exp(ZICMP.shark$theta.hat$gamma), pi = as.numeric(shark.type[i] == 0)))))

d.scmp <- function(x,NuOfVar,Lambda,Nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dSCMP(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu))
}
ll.ziscmp2.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=2,Lambda=shark.type[i], Nu=ZIsCMP2.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))
ll.ziscmp3.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=3,Lambda=shark.type[i], Nu=ZIsCMP3.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))

R2.ZIP <- 1-((ll.zip.sat-as.numeric(ZIP.ll)+ZIP.pars+.5)/(ll.zip.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="poisson")))))
R2.ZINB <- 1-((ll.zinb.sat-as.numeric(ZINB.ll)+ZINB.pars+.5)/(ll.zinb.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="negbin")))))
R2.ZIG <- 1-((ll.zig.sat-as.numeric(ZIG.ll)+ZIG.pars+.5)/(ll.zig.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="geometric")))))
R2.ZIGP <- 1-((ll.zigp.sat-as.numeric(ZIGP.ll)+ZIGP.pars+.5)/(ll.zigp.sat-tail(est.zigp2(Yin=shark.type,fm.X=~1,fm.W=~1,fm.Z=~1,init=T,tex=FALSE),1)))
R2.ZICMP <- 1-((ll.zicmp.sat-as.numeric(ZICMP.ll)+ZICMP.pars+.5)/(ll.zicmp.sat-fit.zicmp.reg2(y=shark.type,X=cbind(rep(1,n)),S=cbind(rep(1,n)),W=cbind(rep(1,n)),beta.init=1,gamma.init=1,zeta.init=1,max=250,optim.control=list(maxit=2000))$loglik))
R2.ZIsCMP2 <- 1-((ll.ziscmp2.sat-as.numeric(ZIsCMP2.ll)+ZIsCMP2.pars+.5)/(ll.ziscmp2.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))
R2.ZIsCMP3 <- 1-((ll.ziscmp3.sat-as.numeric(ZIsCMP3.ll)+ZIsCMP3.pars+.5)/(ll.ziscmp3.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))

sum.out <- data.frame(ll=c(P.ll,NB.ll,ZIP.ll,ZIG.ll,ZINB.ll,ZIGP.ll,ZICMP.ll,ZIsCMP2.ll,ZIsCMP3.ll),
                      AIC=c(P.AIC,NB.AIC,ZIP.AIC,ZIG.AIC,ZINB.AIC,ZIGP.AIC,ZICMP.AIC,ZIsCMP2.AIC,ZIsCMP3.AIC),
                      BIC=c(P.BIC,NB.BIC,ZIP.BIC,ZIG.BIC,ZINB.BIC,ZIGP.BIC,ZICMP.BIC,ZIsCMP2.BIC,ZIsCMP3.BIC),
                      R2.adj=c(R2.P,R2.NB,R2.ZIP,R2.ZIG,R2.ZINB,R2.ZIGP,R2.ZICMP,R2.ZIsCMP2,R2.ZIsCMP3))
rownames(sum.out) <- c("Poi","NegBin","ZIP","ZIG","ZINB","ZIGP","ZICMP","ZIsCMP2","ZIsCMP3")

sum.out

save.image("SNS.RData")

ZICMP.JK <- vector("list",n)
ZIsCMP2.JK <- vector("list",n)
ZIsCMP3.JK <- vector("list",n)
for(i in 1:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-4.680621250,0.004074857),gamma.init=1.1,zeta.init=c(-0.914944825,0.003777524),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-10.1901,0.0435748),gamma0=10.9283,xi0=c(2.68428,0.00172647),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}


#i=453
ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-10.1901,0.0435748),gamma0=10.9283,xi0=c(2.68428,0.00172647),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-5, rel.tol = 1e-5, xf.tol=1e-5, step.min=.001, step.max=1))
ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-5, rel.tol = 1e-5, xf.tol=1e-5, step.min=.001, step.max=1))


for(i in 454:n){
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-10.1901,0.0435748),gamma0=10.9283,xi0=c(2.68428,0.00172647),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}

NB.JK <- rep(NA,n)
for(i in 1:n){
  NB.JK[i] <- log(glm.nb(shark.type[-i]~shark$SoakTime[-i])$theta)
  print(i)
}

cbind(coef(P.shark),sqrt(diag(summary(P.shark)$cov.scaled)))
#cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),summary(NB.shark)$SE.theta^2/summary(NB.shark)$theta^2))
cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),sqrt(var(NB.JK)*((n-1)^2/n^2))))
cbind(unlist(coef(ZIP.shark)),sqrt(diag(ZIP.shark$vcov)))
cbind(unlist(coef(ZIG.shark)),sqrt(diag(ZIG.shark$vcov)))
cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),c(sqrt(diag(ZINB.shark$vcov)),ZINB.shark$SE.logtheta))
est.zigp(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)
ZICMP_1 <- sapply(1:n,function(i) unlist(ZICMP.JK[[i]]$theta));cbind(apply(ZICMP_1,1,mean),sqrt(apply(ZICMP_1,1,var)*((n-1)^2/n^2)))
ZIsCMP2_1 <- sapply(1:n,function(i) unlist(ZIsCMP2.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP2_1,1,mean),sqrt(apply(ZIsCMP2_1,1,var)*((n-1)^2/n^2)))
ZIsCMP3_1 <- sapply(1:n,function(i) unlist(ZIsCMP3.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP3_1,1,mean),sqrt(apply(ZIsCMP3_1,1,var)*((n-1)^2/n^2)))

save.image("SNS_SE.RData")

set.seed(100)

rq.1 <- zip.rqres(y=shark.type,x=shark$SoakTime,z=shark$SoakTime,
                  b0=coef(ZIP.shark)[1],b1=coef(ZIP.shark)[2],
                  a0=coef(ZIP.shark)[3],a1=coef(ZIP.shark)[4])

rq.2 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 2)

rq.3 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 3)

df=data.frame(rq.1,rq.2,rq.3)

ggplot(df, aes(sample = rq.1))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SNS Randomized Quantile Residuals (ZIP Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.2))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SNS Randomized Quantile Residuals (ZISCMP(m=2) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.3))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SNS Randomized Quantile Residuals (ZISCMP(m=3) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

###################################################################
###################################################################
### SBS Analysis
###################################################################
###################################################################

shark.type <- shark[,8]
n <- length(shark.type)
table(shark.type)

P.shark <- glm(shark.type~shark$SoakTime,family=poisson)
NB.shark <- glm.nb(shark.type~shark$SoakTime)

ZIP.shark <- zeroinfl(shark.type~shark$SoakTime,dist="poisson")
ZIG.shark <- zeroinfl(shark.type~shark$SoakTime,dist="geometric")
ZINB.shark <- zeroinfl(shark.type~shark$SoakTime,dist="negbin")

ZIGP.shark <- est.zigp2(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)

ZICMP.shark <- fit.zicmp.reg2(y=shark.type,X=cbind(1,shark$SoakTime),S=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta.init=ZIGP.shark[1:2],gamma.init=ZIGP.shark[3],zeta.init=ZIGP.shark[4:5],max=250,optim.control=list(maxit=2000))

ZIsCMP2.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

ZIsCMP3.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZIsCMP2.shark$Parameters[1:2,1],gamma0=log(ZIsCMP2.shark$Parameters[3,1]),xi0=ZIsCMP2.shark$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

P.ll <- logLik(P.shark)
NB.ll <- logLik(NB.shark)
ZIP.ll <- logLik(ZIP.shark)
ZIG.ll <- logLik(ZIG.shark)
ZINB.ll <- logLik(ZINB.shark)
ZIGP.ll <- tail(ZIGP.shark,1)
ZICMP.ll <- ZICMP.shark$loglik
ZIsCMP2.ll <- ZIsCMP2.shark[[2]]
ZIsCMP3.ll <- ZIsCMP3.shark[[2]]

P.pars <- length(coef(P.shark))
NB.pars <- P.pars + 1
ZIP.pars <- length(coef(ZIP.shark))
ZIG.pars <- ZIP.pars
ZINB.pars <- ZIP.pars+1
ZIGP.pars <- length(ZIGP.shark)-1
ZICMP.pars <- ZIGP.pars
ZIsCMP2.pars <- ZIGP.pars
ZIsCMP3.pars <- ZIGP.pars

P.AIC <- -2*P.ll+2*P.pars
NB.AIC <- -2*NB.ll+2*NB.pars
ZIP.AIC <- -2*ZIP.ll+2*ZIP.pars
ZIG.AIC <- -2*ZIG.ll+2*ZIG.pars
ZINB.AIC <- -2*ZINB.ll+2*ZINB.pars
ZIGP.AIC <- -2*ZIGP.ll+2*ZIGP.pars
ZICMP.AIC <- -2*ZICMP.ll+2*ZICMP.pars
ZIsCMP2.AIC <- -2*ZIsCMP2.ll+2*ZIsCMP2.pars
ZIsCMP3.AIC <- -2*ZIsCMP3.ll+2*ZIsCMP3.pars

P.BIC <- -2*P.ll+log(n)*P.pars
NB.BIC <- -2*NB.ll+log(n)*NB.pars
ZIP.BIC <- -2*ZIP.ll+log(n)*ZIP.pars
ZIG.BIC <- -2*ZIG.ll+log(n)*ZIG.pars
ZINB.BIC <- -2*ZINB.ll+log(n)*ZINB.pars
ZIGP.BIC <- -2*ZIGP.ll+log(n)*ZIGP.pars
ZICMP.BIC <- -2*ZICMP.ll+log(n)*ZICMP.pars
ZIsCMP2.BIC <- -2*ZIsCMP2.ll+log(n)*ZIsCMP2.pars
ZIsCMP3.BIC <- -2*ZIsCMP3.ll+log(n)*ZIsCMP3.pars

R2.P <- 1-(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(P.ll)+P.pars+.5)/(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(logLik(glm(shark.type~1,family=poisson))))
R2.NB <- 1-(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(NB.ll)+NB.pars+.5)/(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(logLik(glm.nb(shark.type~1))))

d.zip <- function(x, lambda , pi) {
  ifelse(x == 0, pi + (1-pi)*dpois(0,lambda = lambda), (1-pi)*dpois(x,lambda = lambda))
}
ll.zip.sat <- sum(log(d.zip(shark.type,shark.type,as.numeric(shark.type == 0))))

d.zinb <- function(x,size,mu,pi) {
  ifelse(x == 0, pi + (1-pi)*dnbinom(0,size = size, mu = mu), (1-pi)*dnbinom(x,size = size, mu = mu))
}
ll.zinb.sat <- sum(log(d.zinb(shark.type,size = NB.shark$theta , mu = shark.type, pi = as.numeric(shark.type == 0))))

d.geom <- function(x,prob,pi) {
  ifelse(x == 0, pi + (1-pi)*dgeom(0,prob = prob), (1-pi)*dgeom(x,prob = prob))
}
ll.zig.sat <- sum(log(d.geom(shark.type,prob=1/(1+mean(shark.type)), pi = as.numeric(shark.type == 0))))

d.genpoi <- function(x,mu,phi,pi) {
  ifelse(x == 0, pi + (1-pi)*dzigp(0,mu=mu,phi=phi,omega=0), (1-pi)*dzigp(x,mu=mu,phi=phi,omega=0))
}
ll.zigp.sat <- sum(sapply(1:n,function(i) log(d.genpoi(shark.type[i],mu=shark.type[i], phi=1+exp(ZIGP.shark[3]), pi = as.numeric(shark.type[i] == 0)))))

d.cmp <- function(x,lambda,nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dcom(x,lambda=lambda,nu=nu))
}
ll.zicmp.sat <- sum(sapply(1:n,function(i) log(d.cmp(shark.type[i],lambda=shark.type[i], nu=exp(ZICMP.shark$theta.hat$gamma), pi = as.numeric(shark.type[i] == 0)))))

d.scmp <- function(x,NuOfVar,Lambda,Nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dSCMP(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu))
}
ll.ziscmp2.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=2,Lambda=shark.type[i], Nu=ZIsCMP2.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))
ll.ziscmp3.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=3,Lambda=shark.type[i], Nu=ZIsCMP3.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))

R2.ZIP <- 1-((ll.zip.sat-as.numeric(ZIP.ll)+ZIP.pars+.5)/(ll.zip.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="poisson")))))
R2.ZINB <- 1-((ll.zinb.sat-as.numeric(ZINB.ll)+ZINB.pars+.5)/(ll.zinb.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="negbin")))))
R2.ZIG <- 1-((ll.zig.sat-as.numeric(ZIG.ll)+ZIG.pars+.5)/(ll.zig.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="geometric")))))
R2.ZIGP <- 1-((ll.zigp.sat-as.numeric(ZIGP.ll)+ZIGP.pars+.5)/(ll.zigp.sat-tail(est.zigp2(Yin=shark.type,fm.X=~1,fm.W=~1,fm.Z=~1,init=T,tex=FALSE),1)))
R2.ZICMP <- 1-((ll.zicmp.sat-as.numeric(ZICMP.ll)+ZICMP.pars+.5)/(ll.zicmp.sat-fit.zicmp.reg2(y=shark.type,X=cbind(rep(1,n)),S=cbind(rep(1,n)),W=cbind(rep(1,n)),beta.init=1,gamma.init=1,zeta.init=1,max=250,optim.control=list(maxit=2000))$loglik))
R2.ZIsCMP2 <- 1-((ll.ziscmp2.sat-as.numeric(ZIsCMP2.ll)+ZIsCMP2.pars+.5)/(ll.ziscmp2.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))
R2.ZIsCMP3 <- 1-((ll.ziscmp3.sat-as.numeric(ZIsCMP3.ll)+ZIsCMP3.pars+.5)/(ll.ziscmp3.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))

sum.out <- data.frame(ll=c(P.ll,NB.ll,ZIP.ll,ZIG.ll,ZINB.ll,ZIGP.ll,ZICMP.ll,ZIsCMP2.ll,ZIsCMP3.ll),
                      AIC=c(P.AIC,NB.AIC,ZIP.AIC,ZIG.AIC,ZINB.AIC,ZIGP.AIC,ZICMP.AIC,ZIsCMP2.AIC,ZIsCMP3.AIC),
                      BIC=c(P.BIC,NB.BIC,ZIP.BIC,ZIG.BIC,ZINB.BIC,ZIGP.BIC,ZICMP.BIC,ZIsCMP2.BIC,ZIsCMP3.BIC),
                      R2.adj=c(R2.P,R2.NB,R2.ZIP,R2.ZIG,R2.ZINB,R2.ZIGP,R2.ZICMP,R2.ZIsCMP2,R2.ZIsCMP3))
rownames(sum.out) <- c("Poi","NegBin","ZIP","ZIG","ZINB","ZIGP","ZICMP","ZIsCMP2","ZIsCMP3")

sum.out

save.image("SBS.RData")

ZICMP.JK <- vector("list",n)
ZIsCMP2.JK <- vector("list",n)
ZIsCMP3.JK <- vector("list",n)
for(i in 1:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-4.680621250,0.004074857),gamma.init=1.1,zeta.init=c(-0.914944825,0.003777524),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-10.1901,0.0435748),gamma0=10.9283,xi0=c(2.68428,0.00172647),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}

ZIG.JK <- vector("list",n)
ZINB.JK <- vector("list",n)
for(i in 1:n){
  ZIG.JK[[i]] <- zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="geometric")
  ZINB.JK[[i]] <- zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="negbin")
  print(i)
}

NB.JK <- rep(NA,n)
ZIG.JK <- NULL
ZINB.JK <- NULL
for(i in 1:n){
  NB.JK[i] <- log(glm.nb(shark.type[-i]~shark$SoakTime[-i])$theta)
  ZIG.JK <- cbind(ZIG.JK,unlist(coef(zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="geometric"))))
  tmp <- zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="negbin")
  ZINB.JK <- cbind(ZINB.JK,c(unlist(tmp$coefficients),log(tmp$theta)))
  print(i)
}

cbind(coef(P.shark),sqrt(diag(summary(P.shark)$cov.scaled)))
#cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),summary(NB.shark)$SE.theta^2/summary(NB.shark)$theta^2))
cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),sqrt(var(NB.JK)*((n-1)^2/n^2))))
cbind(unlist(coef(ZIP.shark)),sqrt(diag(ZIP.shark$vcov)))
#cbind(unlist(coef(ZIG.shark)),sqrt(diag(ZIG.shark$vcov)))
cbind(unlist(coef(ZIG.shark)),sqrt(apply(ZIG.JK,1,var)*((n-1)^2/n^2)))
#ZIG_1 <- sapply(1:n,function(i) unlist(ZIG.JK[[i]]$coefficients));cbind(unlist(coef(ZIG.shark)),sqrt(apply(ZIG_1,1,var)*((n-1)^2/n^2)))
#cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),c(sqrt(diag(ZINB.shark$vcov)),ZINB.shark$SE.logtheta))
cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),sqrt(apply(ZINB.JK,1,var)*((n-1)^2/n^2)))
est.zigp(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)
ZICMP_1 <- sapply(1:n,function(i) unlist(ZICMP.JK[[i]]$theta));cbind(apply(ZICMP_1,1,mean),sqrt(apply(ZICMP_1,1,var)*((n-1)^2/n^2)))
ZIsCMP2_1 <- sapply(1:n,function(i) unlist(ZIsCMP2.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP2_1,1,mean),sqrt(apply(ZIsCMP2_1,1,var)*((n-1)^2/n^2)))
ZIsCMP3_1 <- sapply(1:n,function(i) unlist(ZIsCMP3.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP3_1,1,mean),sqrt(apply(ZIsCMP3_1,1,var)*((n-1)^2/n^2)))


save.image("SBS_SE.RData")

set.seed(100)

rq.1 <-  nb.rqres(y=shark.type,x=shark$SoakTime,
                  b0=coef(NB.shark)[1],b1=coef(NB.shark)[2],theta=NB.shark$theta)

rq.2 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 2)

rq.3 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 3)

df=data.frame(rq.1,rq.2,rq.3)

ggplot(df, aes(sample = rq.1))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SBS Randomized Quantile Residuals (NB Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.2))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SBS Randomized Quantile Residuals (ZISCMP(m=2) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.3))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SBS Randomized Quantile Residuals (ZISCMP(m=3) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)




###################################################################
###################################################################
### WSH Analysis
###################################################################
###################################################################

shark.type <- shark[,9]
n <- length(shark.type)
table(shark.type)

P.shark <- glm(shark.type~shark$SoakTime,family=poisson)
NB.shark <- glm.nb(shark.type~shark$SoakTime)

ZIP.shark <- zeroinfl(shark.type~shark$SoakTime,dist="poisson")
ZIG.shark <- zeroinfl(shark.type~shark$SoakTime,dist="geometric")
ZINB.shark <- zeroinfl(shark.type~shark$SoakTime,dist="negbin")

ZIGP.shark <- est.zigp2(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)

ZICMP.shark <- fit.zicmp.reg2(y=shark.type,X=cbind(1,shark$SoakTime),S=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta.init=c(-2.0440418307,0.0006571716),gamma.init=-8.67944,zeta.init=c(1.02962938,0.00026627),max=250,optim.control=list(maxit=2000))

ZIsCMP2.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=c(-2.0440418307,0.0006571716),gamma0=exp(-8.67944),xi0=c(1.02962938,0.00026627),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

ZIsCMP3.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZIsCMP2.shark$Parameters[1:2,1],gamma0=log(ZIsCMP2.shark$Parameters[3,1]),xi0=ZIsCMP2.shark$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

P.ll <- logLik(P.shark)
NB.ll <- logLik(NB.shark)
ZIP.ll <- logLik(ZIP.shark)
ZIG.ll <- logLik(ZIG.shark)
ZINB.ll <- logLik(ZINB.shark)
ZIGP.ll <- tail(ZIGP.shark,1)
ZICMP.ll <- ZICMP.shark$loglik
ZIsCMP2.ll <- ZIsCMP2.shark[[2]]
ZIsCMP3.ll <- ZIsCMP3.shark[[2]]

P.pars <- length(coef(P.shark))
NB.pars <- P.pars + 1
ZIP.pars <- length(coef(ZIP.shark))
ZIG.pars <- ZIP.pars
ZINB.pars <- ZIP.pars+1
ZIGP.pars <- length(ZIGP.shark)-1
ZICMP.pars <- ZIGP.pars
ZIsCMP2.pars <- ZIGP.pars
ZIsCMP3.pars <- ZIGP.pars

P.AIC <- -2*P.ll+2*P.pars
NB.AIC <- -2*NB.ll+2*NB.pars
ZIP.AIC <- -2*ZIP.ll+2*ZIP.pars
ZIG.AIC <- -2*ZIG.ll+2*ZIG.pars
ZINB.AIC <- -2*ZINB.ll+2*ZINB.pars
ZIGP.AIC <- -2*ZIGP.ll+2*ZIGP.pars
ZICMP.AIC <- -2*ZICMP.ll+2*ZICMP.pars
ZIsCMP2.AIC <- -2*ZIsCMP2.ll+2*ZIsCMP2.pars
ZIsCMP3.AIC <- -2*ZIsCMP3.ll+2*ZIsCMP3.pars

P.BIC <- -2*P.ll+log(n)*P.pars
NB.BIC <- -2*NB.ll+log(n)*NB.pars
ZIP.BIC <- -2*ZIP.ll+log(n)*ZIP.pars
ZIG.BIC <- -2*ZIG.ll+log(n)*ZIG.pars
ZINB.BIC <- -2*ZINB.ll+log(n)*ZINB.pars
ZIGP.BIC <- -2*ZIGP.ll+log(n)*ZIGP.pars
ZICMP.BIC <- -2*ZICMP.ll+log(n)*ZICMP.pars
ZIsCMP2.BIC <- -2*ZIsCMP2.ll+log(n)*ZIsCMP2.pars
ZIsCMP3.BIC <- -2*ZIsCMP3.ll+log(n)*ZIsCMP3.pars

R2.P <- 1-(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(P.ll)+P.pars+.5)/(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(logLik(glm(shark.type~1,family=poisson))))
R2.NB <- 1-(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(NB.ll)+NB.pars+.5)/(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(logLik(glm.nb(shark.type~1))))

d.zip <- function(x, lambda , pi) {
  ifelse(x == 0, pi + (1-pi)*dpois(0,lambda = lambda), (1-pi)*dpois(x,lambda = lambda))
}
ll.zip.sat <- sum(log(d.zip(shark.type,shark.type,as.numeric(shark.type == 0))))

d.zinb <- function(x,size,mu,pi) {
  ifelse(x == 0, pi + (1-pi)*dnbinom(0,size = size, mu = mu), (1-pi)*dnbinom(x,size = size, mu = mu))
}
ll.zinb.sat <- sum(log(d.zinb(shark.type,size = NB.shark$theta , mu = shark.type, pi = as.numeric(shark.type == 0))))

d.geom <- function(x,prob,pi) {
  ifelse(x == 0, pi + (1-pi)*dgeom(0,prob = prob), (1-pi)*dgeom(x,prob = prob))
}
ll.zig.sat <- sum(log(d.geom(shark.type,prob=1/(1+mean(shark.type)), pi = as.numeric(shark.type == 0))))

d.genpoi <- function(x,mu,phi,pi) {
  ifelse(x == 0, pi + (1-pi)*dzigp(0,mu=mu,phi=phi,omega=0), (1-pi)*dzigp(x,mu=mu,phi=phi,omega=0))
}
ll.zigp.sat <- sum(sapply(1:n,function(i) log(d.genpoi(shark.type[i],mu=shark.type[i], phi=1+exp(ZIGP.shark[3]), pi = as.numeric(shark.type[i] == 0)))))

d.cmp <- function(x,lambda,nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dcom(x,lambda=lambda,nu=nu))
}
ll.zicmp.sat <- sum(sapply(1:n,function(i) log(d.cmp(shark.type[i],lambda=shark.type[i], nu=exp(ZICMP.shark$theta.hat$gamma), pi = as.numeric(shark.type[i] == 0)))))

d.scmp <- function(x,NuOfVar,Lambda,Nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dSCMP(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu))
}
ll.ziscmp2.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=2,Lambda=shark.type[i], Nu=ZIsCMP2.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))
ll.ziscmp3.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=3,Lambda=shark.type[i], Nu=ZIsCMP3.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))

R2.ZIP <- 1-((ll.zip.sat-as.numeric(ZIP.ll)+ZIP.pars+.5)/(ll.zip.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="poisson")))))
R2.ZINB <- 1-((ll.zinb.sat-as.numeric(ZINB.ll)+ZINB.pars+.5)/(ll.zinb.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="negbin")))))
R2.ZIG <- 1-((ll.zig.sat-as.numeric(ZIG.ll)+ZIG.pars+.5)/(ll.zig.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="geometric")))))
R2.ZIGP <- 1-((ll.zigp.sat-as.numeric(ZIGP.ll)+ZIGP.pars+.5)/(ll.zigp.sat-tail(est.zigp2(Yin=shark.type,fm.X=~1,fm.W=~1,fm.Z=~1,init=T,tex=FALSE),1)))
R2.ZICMP <- 1-((ll.zicmp.sat-as.numeric(ZICMP.ll)+ZICMP.pars+.5)/(ll.zicmp.sat-fit.zicmp.reg2(y=shark.type,X=cbind(rep(1,n)),S=cbind(rep(1,n)),W=cbind(rep(1,n)),beta.init=1,gamma.init=1,zeta.init=1,max=250,optim.control=list(maxit=2000))$loglik))
R2.ZIsCMP2 <- 1-((ll.ziscmp2.sat-as.numeric(ZIsCMP2.ll)+ZIsCMP2.pars+.5)/(ll.ziscmp2.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))
R2.ZIsCMP3 <- 1-((ll.ziscmp3.sat-as.numeric(ZIsCMP3.ll)+ZIsCMP3.pars+.5)/(ll.ziscmp3.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))

sum.out <- data.frame(ll=c(P.ll,NB.ll,ZIP.ll,ZIG.ll,ZINB.ll,ZIGP.ll,ZICMP.ll,ZIsCMP2.ll,ZIsCMP3.ll),
                      AIC=c(P.AIC,NB.AIC,ZIP.AIC,ZIG.AIC,ZINB.AIC,ZIGP.AIC,ZICMP.AIC,ZIsCMP2.AIC,ZIsCMP3.AIC),
                      BIC=c(P.BIC,NB.BIC,ZIP.BIC,ZIG.BIC,ZINB.BIC,ZIGP.BIC,ZICMP.BIC,ZIsCMP2.BIC,ZIsCMP3.BIC),
                      R2.adj=c(R2.P,R2.NB,R2.ZIP,R2.ZIG,R2.ZINB,R2.ZIGP,R2.ZICMP,R2.ZIsCMP2,R2.ZIsCMP3))
rownames(sum.out) <- c("Poi","NegBin","ZIP","ZIG","ZINB","ZIGP","ZICMP","ZIsCMP2","ZIsCMP3")

sum.out

save.image("WSH.RData")

ZICMP.JK <- vector("list",n)
ZIsCMP2.JK <- vector("list",n)
ZIsCMP3.JK <- vector("list",n)
for(i in 1:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-2.0440418307,0.0006571716),gamma.init=-8.67944,zeta.init=c(1.02962938,0.00026627),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-2.0440418307,0.0006571716),gamma0=exp(-8.67944),xi0=c(1.02962938,0.00026627),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}


#i=264
ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-2.0440418307,0.0006571716),gamma0=exp(-8.67944),xi0=c(1.02962938,0.00026627),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-5, rel.tol = 1e-5, xf.tol=1e-5, step.min=.001, step.max=1))
ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-5, rel.tol = 1e-5, xf.tol=1e-5, step.min=.001, step.max=1))

for(i in 265:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-2.0440418307,0.0006571716),gamma.init=-8.67944,zeta.init=c(1.02962938,0.00026627),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=c(-2.0440418307,0.0006571716),gamma0=exp(-8.67944),xi0=c(1.02962938,0.00026627),NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}


NB.JK <- rep(NA,n)
for(i in 1:n){
  NB.JK[i] <- log(glm.nb(shark.type[-i]~shark$SoakTime[-i])$theta)
  print(i)
}

cbind(coef(P.shark),sqrt(diag(summary(P.shark)$cov.scaled)))
#cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),summary(NB.shark)$SE.theta^2/summary(NB.shark)$theta^2))
cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),sqrt(var(NB.JK)*((n-1)^2/n^2))))
cbind(unlist(coef(ZIP.shark)),sqrt(diag(ZIP.shark$vcov)))
cbind(unlist(coef(ZIG.shark)),sqrt(diag(ZIG.shark$vcov)))
cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),c(sqrt(diag(ZINB.shark$vcov)),ZINB.shark$SE.logtheta))
est.zigp(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)
ZICMP_1 <- sapply(1:n,function(i) unlist(ZICMP.JK[[i]]$theta));cbind(apply(ZICMP_1,1,mean),sqrt(apply(ZICMP_1,1,var)*((n-1)^2/n^2)))
ZIsCMP2_1 <- sapply(1:n,function(i) unlist(ZIsCMP2.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP2_1,1,mean),sqrt(apply(ZIsCMP2_1,1,var)*((n-1)^2/n^2)))
ZIsCMP3_1 <- sapply(1:n,function(i) unlist(ZIsCMP3.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP3_1,1,mean),sqrt(apply(ZIsCMP3_1,1,var)*((n-1)^2/n^2)))


save.image("WSH_SE.RData")

set.seed(100)

rq.1 <- zip.rqres(y=shark.type,x=shark$SoakTime,z=shark$SoakTime,
                  b0=coef(ZIP.shark)[1],b1=coef(ZIP.shark)[2],
                  a0=coef(ZIP.shark)[3],a1=coef(ZIP.shark)[4])

rq.2 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 2)

rq.3 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 3)

df=data.frame(rq.1,rq.2,rq.3)

ggplot(df, aes(sample = rq.1))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("WSH Randomized Quantile Residuals (ZIP Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.2))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("WSH Randomized Quantile Residuals (ZISCMP(m=2) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.3))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("WSH Randomized Quantile Residuals (ZISCMP(m=3) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)


###################################################################
###################################################################
### UBT Analysis
###################################################################
###################################################################

shark.type <- shark[,10]
n <- length(shark.type)
table(shark.type)

P.shark <- glm(shark.type~shark$SoakTime,family=poisson)
NB.shark <- glm.nb(shark.type~shark$SoakTime)

ZIP.shark <- zeroinfl(shark.type~shark$SoakTime,dist="poisson")
ZIG.shark <- zeroinfl(shark.type~shark$SoakTime,dist="geometric")
ZINB.shark <- zeroinfl(shark.type~shark$SoakTime,dist="negbin")

ZIGP.shark <- est.zigp2(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)

ZICMP.shark <- fit.zicmp.reg2(y=shark.type,X=cbind(1,shark$SoakTime),S=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta.init=c(-0.075125957,0.001710425),gamma.init=0.1934497,zeta.init=c(2.083981045,0.004844035),max=250,optim.control=list(maxit=2000))

ZIsCMP2.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

ZIsCMP3.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZIsCMP2.shark$Parameters[1:2,1],gamma0=log(ZIsCMP2.shark$Parameters[3,1]),xi0=ZIsCMP2.shark$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

P.ll <- logLik(P.shark)
NB.ll <- logLik(NB.shark)
ZIP.ll <- logLik(ZIP.shark)
ZIG.ll <- logLik(ZIG.shark)
ZINB.ll <- logLik(ZINB.shark)
ZIGP.ll <- tail(ZIGP.shark,1)
ZICMP.ll <- ZICMP.shark$loglik
ZIsCMP2.ll <- ZIsCMP2.shark[[2]]
ZIsCMP3.ll <- ZIsCMP3.shark[[2]]

P.pars <- length(coef(P.shark))
NB.pars <- P.pars + 1
ZIP.pars <- length(coef(ZIP.shark))
ZIG.pars <- ZIP.pars
ZINB.pars <- ZIP.pars+1
ZIGP.pars <- length(ZIGP.shark)-1
ZICMP.pars <- ZIGP.pars
ZIsCMP2.pars <- ZIGP.pars
ZIsCMP3.pars <- ZIGP.pars

P.AIC <- -2*P.ll+2*P.pars
NB.AIC <- -2*NB.ll+2*NB.pars
ZIP.AIC <- -2*ZIP.ll+2*ZIP.pars
ZIG.AIC <- -2*ZIG.ll+2*ZIG.pars
ZINB.AIC <- -2*ZINB.ll+2*ZINB.pars
ZIGP.AIC <- -2*ZIGP.ll+2*ZIGP.pars
ZICMP.AIC <- -2*ZICMP.ll+2*ZICMP.pars
ZIsCMP2.AIC <- -2*ZIsCMP2.ll+2*ZIsCMP2.pars
ZIsCMP3.AIC <- -2*ZIsCMP3.ll+2*ZIsCMP3.pars

P.BIC <- -2*P.ll+log(n)*P.pars
NB.BIC <- -2*NB.ll+log(n)*NB.pars
ZIP.BIC <- -2*ZIP.ll+log(n)*ZIP.pars
ZIG.BIC <- -2*ZIG.ll+log(n)*ZIG.pars
ZINB.BIC <- -2*ZINB.ll+log(n)*ZINB.pars
ZIGP.BIC <- -2*ZIGP.ll+log(n)*ZIGP.pars
ZICMP.BIC <- -2*ZICMP.ll+log(n)*ZICMP.pars
ZIsCMP2.BIC <- -2*ZIsCMP2.ll+log(n)*ZIsCMP2.pars
ZIsCMP3.BIC <- -2*ZIsCMP3.ll+log(n)*ZIsCMP3.pars

R2.P <- 1-(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(P.ll)+P.pars+.5)/(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(logLik(glm(shark.type~1,family=poisson))))
R2.NB <- 1-(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(NB.ll)+NB.pars+.5)/(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(logLik(glm.nb(shark.type~1))))

d.zip <- function(x, lambda , pi) {
  ifelse(x == 0, pi + (1-pi)*dpois(0,lambda = lambda), (1-pi)*dpois(x,lambda = lambda))
}
ll.zip.sat <- sum(log(d.zip(shark.type,shark.type,as.numeric(shark.type == 0))))

d.zinb <- function(x,size,mu,pi) {
  ifelse(x == 0, pi + (1-pi)*dnbinom(0,size = size, mu = mu), (1-pi)*dnbinom(x,size = size, mu = mu))
}
ll.zinb.sat <- sum(log(d.zinb(shark.type,size = NB.shark$theta , mu = shark.type, pi = as.numeric(shark.type == 0))))

d.geom <- function(x,prob,pi) {
  ifelse(x == 0, pi + (1-pi)*dgeom(0,prob = prob), (1-pi)*dgeom(x,prob = prob))
}
ll.zig.sat <- sum(log(d.geom(shark.type,prob=1/(1+mean(shark.type)), pi = as.numeric(shark.type == 0))))

d.genpoi <- function(x,mu,phi,pi) {
  ifelse(x == 0, pi + (1-pi)*dzigp(0,mu=mu,phi=phi,omega=0), (1-pi)*dzigp(x,mu=mu,phi=phi,omega=0))
}
ll.zigp.sat <- sum(sapply(1:n,function(i) log(d.genpoi(shark.type[i],mu=shark.type[i], phi=1+exp(ZIGP.shark[3]), pi = as.numeric(shark.type[i] == 0)))))

d.cmp <- function(x,lambda,nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dcom(x,lambda=lambda,nu=nu))
}
ll.zicmp.sat <- sum(sapply(1:n,function(i) log(d.cmp(shark.type[i],lambda=shark.type[i], nu=exp(ZICMP.shark$theta.hat$gamma), pi = as.numeric(shark.type[i] == 0)))))

d.scmp <- function(x,NuOfVar,Lambda,Nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dSCMP(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu))
}
ll.ziscmp2.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=2,Lambda=shark.type[i], Nu=ZIsCMP2.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))
ll.ziscmp3.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=3,Lambda=shark.type[i], Nu=ZIsCMP3.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))

R2.ZIP <- 1-((ll.zip.sat-as.numeric(ZIP.ll)+ZIP.pars+.5)/(ll.zip.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="poisson")))))
R2.ZINB <- 1-((ll.zinb.sat-as.numeric(ZINB.ll)+ZINB.pars+.5)/(ll.zinb.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="negbin")))))
R2.ZIG <- 1-((ll.zig.sat-as.numeric(ZIG.ll)+ZIG.pars+.5)/(ll.zig.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="geometric")))))
R2.ZIGP <- 1-((ll.zigp.sat-as.numeric(ZIGP.ll)+ZIGP.pars+.5)/(ll.zigp.sat-tail(est.zigp2(Yin=shark.type,fm.X=~1,fm.W=~1,fm.Z=~1,init=T,tex=FALSE),1)))
R2.ZICMP <- 1-((ll.zicmp.sat-as.numeric(ZICMP.ll)+ZICMP.pars+.5)/(ll.zicmp.sat-fit.zicmp.reg2(y=shark.type,X=cbind(rep(1,n)),S=cbind(rep(1,n)),W=cbind(rep(1,n)),beta.init=1,gamma.init=1,zeta.init=1,max=250,optim.control=list(maxit=2000))$loglik))
R2.ZIsCMP2 <- 1-((ll.ziscmp2.sat-as.numeric(ZIsCMP2.ll)+ZIsCMP2.pars+.5)/(ll.ziscmp2.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))
R2.ZIsCMP3 <- 1-((ll.ziscmp3.sat-as.numeric(ZIsCMP3.ll)+ZIsCMP3.pars+.5)/(ll.ziscmp3.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))

sum.out <- data.frame(ll=c(P.ll,NB.ll,ZIP.ll,ZIG.ll,ZINB.ll,ZIGP.ll,ZICMP.ll,ZIsCMP2.ll,ZIsCMP3.ll),
                      AIC=c(P.AIC,NB.AIC,ZIP.AIC,ZIG.AIC,ZINB.AIC,ZIGP.AIC,ZICMP.AIC,ZIsCMP2.AIC,ZIsCMP3.AIC),
                      BIC=c(P.BIC,NB.BIC,ZIP.BIC,ZIG.BIC,ZINB.BIC,ZIGP.BIC,ZICMP.BIC,ZIsCMP2.BIC,ZIsCMP3.BIC),
                      R2.adj=c(R2.P,R2.NB,R2.ZIP,R2.ZIG,R2.ZINB,R2.ZIGP,R2.ZICMP,R2.ZIsCMP2,R2.ZIsCMP3))
rownames(sum.out) <- c("Poi","NegBin","ZIP","ZIG","ZINB","ZIGP","ZICMP","ZIsCMP2","ZIsCMP3")

sum.out

save.image("UBT.RData")

ZICMP.JK <- vector("list",n)
ZIsCMP2.JK <- vector("list",n)
ZIsCMP3.JK <- vector("list",n)
for(i in 1:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-0.075125957,0.001710425),gamma.init=0.1934497,zeta.init=c(2.083981045,0.004844035),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}


#i=59,166,175,249,251,309,445,468
ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-4, rel.tol = 1e-4, xf.tol=1e-4, step.min=.001, step.max=1))
ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-4, rel.tol = 1e-4, xf.tol=1e-4, step.min=.001, step.max=1))

for(i in 469:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-0.075125957,0.001710425),gamma.init=0.1934497,zeta.init=c(2.083981045,0.004844035),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.JK[[i]]$Parameters[1:2,1],gamma0=log(ZIsCMP2.JK[[i]]$Parameters[3,1]),xi0=ZIsCMP2.JK[[i]]$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}

NB.JK <- rep(NA,n)
for(i in 1:n){
  NB.JK[i] <- log(glm.nb(shark.type[-i]~shark$SoakTime[-i])$theta)
  print(i)
}

cbind(coef(P.shark),sqrt(diag(summary(P.shark)$cov.scaled)))
#cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),summary(NB.shark)$SE.theta^2/summary(NB.shark)$theta^2))
cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),sqrt(var(NB.JK)*((n-1)^2/n^2))))
cbind(unlist(coef(ZIP.shark)),sqrt(diag(ZIP.shark$vcov)))
cbind(unlist(coef(ZIG.shark)),sqrt(diag(ZIG.shark$vcov)))
cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),c(sqrt(diag(ZINB.shark$vcov)),ZINB.shark$SE.logtheta))
est.zigp(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)
ZICMP_1 <- sapply(1:n,function(i) unlist(ZICMP.JK[[i]]$theta));cbind(apply(ZICMP_1,1,mean),sqrt(apply(ZICMP_1,1,var)*((n-1)^2/n^2)))
ZIsCMP2_1 <- sapply(1:n,function(i) unlist(ZIsCMP2.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP2_1,1,mean),sqrt(apply(ZIsCMP2_1,1,var)*((n-1)^2/n^2)))
ZIsCMP3_1 <- sapply(1:n,function(i) unlist(ZIsCMP3.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP3_1,1,mean),sqrt(apply(ZIsCMP3_1,1,var)*((n-1)^2/n^2)))


save.image("UBT_SE.RData")

set.seed(100)

rq.1 <-  nb.rqres(y=shark.type,x=shark$SoakTime,
                  b0=coef(NB.shark)[1],b1=coef(NB.shark)[2],theta=NB.shark$theta)

rq.2 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 2)

rq.3 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 3)

df=data.frame(rq.1,rq.2,rq.3)

ggplot(df, aes(sample = rq.1))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("UBT Randomized Quantile Residuals (NB Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.2))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("UBT Randomized Quantile Residuals (ZISCMP(m=2) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.3))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("UBT Randomized Quantile Residuals (ZISCMP(m=3) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)



###################################################################
###################################################################
### GHH Analysis
###################################################################
###################################################################

shark.type <- shark[,7]
n <- length(shark.type)
table(shark.type)

P.shark <- glm(shark.type~shark$SoakTime,family=poisson)
NB.shark <- glm.nb(shark.type~shark$SoakTime)

ZIP.shark <- zeroinfl(shark.type~shark$SoakTime,dist="poisson")
ZIG.shark <- zeroinfl(shark.type~shark$SoakTime,dist="geometric")
ZINB.shark <- zeroinfl(shark.type~shark$SoakTime,dist="negbin")

ZIGP.shark <- est.zigp2(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)

ZICMP.shark <- fit.zicmp.reg2(y=shark.type,X=cbind(1,shark$SoakTime),S=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta.init=ZIGP.shark[1:2],gamma.init=ZIGP.shark[3],zeta.init=ZIGP.shark[4:5],max=250,optim.control=list(maxit=2000))

ZIsCMP2.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

ZIsCMP3.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZIsCMP2.shark$Parameters[1:2,1]+c(-.5,.0001),gamma0=log(ZIsCMP2.shark$Parameters[3,1])+2,xi0=ZIsCMP2.shark$Parameters[4:5,1]+c(.5,-.0001),NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

P.ll <- logLik(P.shark)
NB.ll <- logLik(NB.shark)
ZIP.ll <- logLik(ZIP.shark)
ZIG.ll <- logLik(ZIG.shark)
ZINB.ll <- logLik(ZINB.shark)
ZIGP.ll <- tail(ZIGP.shark,1)
ZICMP.ll <- ZICMP.shark$loglik
ZIsCMP2.ll <- ZIsCMP2.shark[[2]]
ZIsCMP3.ll <- ZIsCMP3.shark[[2]]

P.pars <- length(coef(P.shark))
NB.pars <- P.pars + 1
ZIP.pars <- length(coef(ZIP.shark))
ZIG.pars <- ZIP.pars
ZINB.pars <- ZIP.pars+1
ZIGP.pars <- length(ZIGP.shark)-1
ZICMP.pars <- ZIGP.pars
ZIsCMP2.pars <- ZIGP.pars
ZIsCMP3.pars <- ZIGP.pars

P.AIC <- -2*P.ll+2*P.pars
NB.AIC <- -2*NB.ll+2*NB.pars
ZIP.AIC <- -2*ZIP.ll+2*ZIP.pars
ZIG.AIC <- -2*ZIG.ll+2*ZIG.pars
ZINB.AIC <- -2*ZINB.ll+2*ZINB.pars
ZIGP.AIC <- -2*ZIGP.ll+2*ZIGP.pars
ZICMP.AIC <- -2*ZICMP.ll+2*ZICMP.pars
ZIsCMP2.AIC <- -2*ZIsCMP2.ll+2*ZIsCMP2.pars
ZIsCMP3.AIC <- -2*ZIsCMP3.ll+2*ZIsCMP3.pars

P.BIC <- -2*P.ll+log(n)*P.pars
NB.BIC <- -2*NB.ll+log(n)*NB.pars
ZIP.BIC <- -2*ZIP.ll+log(n)*ZIP.pars
ZIG.BIC <- -2*ZIG.ll+log(n)*ZIG.pars
ZINB.BIC <- -2*ZINB.ll+log(n)*ZINB.pars
ZIGP.BIC <- -2*ZIGP.ll+log(n)*ZIGP.pars
ZICMP.BIC <- -2*ZICMP.ll+log(n)*ZICMP.pars
ZIsCMP2.BIC <- -2*ZIsCMP2.ll+log(n)*ZIsCMP2.pars
ZIsCMP3.BIC <- -2*ZIsCMP3.ll+log(n)*ZIsCMP3.pars

R2.P <- 1-(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(P.ll)+P.pars+.5)/(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(logLik(glm(shark.type~1,family=poisson))))
R2.NB <- 1-(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(NB.ll)+NB.pars+.5)/(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(logLik(glm.nb(shark.type~1))))

d.zip <- function(x, lambda , pi) {
  ifelse(x == 0, pi + (1-pi)*dpois(0,lambda = lambda), (1-pi)*dpois(x,lambda = lambda))
}
ll.zip.sat <- sum(log(d.zip(shark.type,shark.type,as.numeric(shark.type == 0))))

d.zinb <- function(x,size,mu,pi) {
  ifelse(x == 0, pi + (1-pi)*dnbinom(0,size = size, mu = mu), (1-pi)*dnbinom(x,size = size, mu = mu))
}
ll.zinb.sat <- sum(log(d.zinb(shark.type,size = NB.shark$theta , mu = shark.type, pi = as.numeric(shark.type == 0))))

d.geom <- function(x,prob,pi) {
  ifelse(x == 0, pi + (1-pi)*dgeom(0,prob = prob), (1-pi)*dgeom(x,prob = prob))
}
ll.zig.sat <- sum(log(d.geom(shark.type,prob=1/(1+mean(shark.type)), pi = as.numeric(shark.type == 0))))

d.genpoi <- function(x,mu,phi,pi) {
  ifelse(x == 0, pi + (1-pi)*dzigp(0,mu=mu,phi=phi,omega=0), (1-pi)*dzigp(x,mu=mu,phi=phi,omega=0))
}
ll.zigp.sat <- sum(sapply(1:n,function(i) log(d.genpoi(shark.type[i],mu=shark.type[i], phi=1+exp(ZIGP.shark[3]), pi = as.numeric(shark.type[i] == 0)))))

d.cmp <- function(x,lambda,nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dcom(x,lambda=lambda,nu=nu))
}
ll.zicmp.sat <- sum(sapply(1:n,function(i) log(d.cmp(shark.type[i],lambda=shark.type[i], nu=exp(ZICMP.shark$theta.hat$gamma), pi = as.numeric(shark.type[i] == 0)))))

d.scmp <- function(x,NuOfVar,Lambda,Nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dSCMP(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu))
}
ll.ziscmp2.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=2,Lambda=shark.type[i], Nu=ZIsCMP2.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))
ll.ziscmp3.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=3,Lambda=shark.type[i], Nu=ZIsCMP3.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))

R2.ZIP <- 1-((ll.zip.sat-as.numeric(ZIP.ll)+ZIP.pars+.5)/(ll.zip.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="poisson")))))
R2.ZINB <- 1-((ll.zinb.sat-as.numeric(ZINB.ll)+ZINB.pars+.5)/(ll.zinb.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="negbin")))))
R2.ZIG <- 1-((ll.zig.sat-as.numeric(ZIG.ll)+ZIG.pars+.5)/(ll.zig.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="geometric")))))
R2.ZIGP <- 1-((ll.zigp.sat-as.numeric(ZIGP.ll)+ZIGP.pars+.5)/(ll.zigp.sat-tail(est.zigp2(Yin=shark.type,fm.X=~1,fm.W=~1,fm.Z=~1,init=T,tex=FALSE),1)))
R2.ZICMP <- 1-((ll.zicmp.sat-as.numeric(ZICMP.ll)+ZICMP.pars+.5)/(ll.zicmp.sat-fit.zicmp.reg2(y=shark.type,X=cbind(rep(1,n)),S=cbind(rep(1,n)),W=cbind(rep(1,n)),beta.init=1,gamma.init=1,zeta.init=1,max=250,optim.control=list(maxit=2000))$loglik))
R2.ZIsCMP2 <- 1-((ll.ziscmp2.sat-as.numeric(ZIsCMP2.ll)+ZIsCMP2.pars+.5)/(ll.ziscmp2.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))
R2.ZIsCMP3 <- 1-((ll.ziscmp3.sat-as.numeric(ZIsCMP3.ll)+ZIsCMP3.pars+.5)/(ll.ziscmp3.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))

sum.out <- data.frame(ll=c(P.ll,NB.ll,ZIP.ll,ZIG.ll,ZINB.ll,ZIGP.ll,ZICMP.ll,ZIsCMP2.ll,ZIsCMP3.ll),
                      AIC=c(P.AIC,NB.AIC,ZIP.AIC,ZIG.AIC,ZINB.AIC,ZIGP.AIC,ZICMP.AIC,ZIsCMP2.AIC,ZIsCMP3.AIC),
                      BIC=c(P.BIC,NB.BIC,ZIP.BIC,ZIG.BIC,ZINB.BIC,ZIGP.BIC,ZICMP.BIC,ZIsCMP2.BIC,ZIsCMP3.BIC),
                      R2.adj=c(R2.P,R2.NB,R2.ZIP,R2.ZIG,R2.ZINB,R2.ZIGP,R2.ZICMP,R2.ZIsCMP2,R2.ZIsCMP3))
rownames(sum.out) <- c("Poi","NegBin","ZIP","ZIG","ZINB","ZIGP","ZICMP","ZIsCMP2","ZIsCMP3")

sum.out

save.image("GHH.RData")


ZICMP.JK <- vector("list",n)
ZIsCMP2.JK <- vector("list",n)
ZIsCMP3.JK <- vector("list",n)
for(i in 1:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=ZIGP.shark[1:2],gamma.init=ZIGP.shark[3],zeta.init=ZIGP.shark[4:5],max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.shark$Parameters[1:2,1]+c(-.5,.0001),gamma0=log(ZIsCMP2.shark$Parameters[3,1])+2,xi0=ZIsCMP2.shark$Parameters[4:5,1]+c(.5,-.0001),NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}


NB.JK <- rep(NA,n)
ZIP.JK <- NULL
ZINB.JK <- NULL
for(i in 1:n){
  NB.JK[i] <- log(glm.nb(shark.type[-i]~shark$SoakTime[-i])$theta)
  ZIP.JK <- cbind(ZIP.JK,unlist(coef(zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="poisson"))))
  tmp <- zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="negbin")
  ZINB.JK <- cbind(ZINB.JK,c(unlist(tmp$coefficients),log(tmp$theta)))
  print(i)
}


cbind(coef(P.shark),sqrt(diag(summary(P.shark)$cov.scaled)))
#cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),summary(NB.shark)$SE.theta^2/summary(NB.shark)$theta^2))
cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),sqrt(var(NB.JK)*((n-1)^2/n^2))))
#cbind(unlist(coef(ZIP.shark)),sqrt(diag(ZIP.shark$vcov)))
cbind(unlist(coef(ZIP.shark)),sqrt(apply(ZIP.JK,1,var)*((n-1)^2/n^2)))
cbind(unlist(coef(ZIG.shark)),sqrt(diag(ZIG.shark$vcov)))
#cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),c(sqrt(diag(ZINB.shark$vcov)),ZINB.shark$SE.logtheta))
cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),sqrt(apply(ZINB.JK,1,var)*((n-1)^2/n^2)))
est.zigp(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)
ZICMP_1 <- sapply(1:n,function(i) unlist(ZICMP.JK[[i]]$theta));cbind(apply(ZICMP_1,1,mean),sqrt(apply(ZICMP_1,1,var)*((n-1)^2/n^2)))
ZIsCMP2_1 <- sapply(1:n,function(i) unlist(ZIsCMP2.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP2_1,1,mean),sqrt(apply(ZIsCMP2_1,1,var)*((n-1)^2/n^2)))
ZIsCMP3_1 <- sapply(1:n,function(i) unlist(ZIsCMP3.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP3_1,1,mean),sqrt(apply(ZIsCMP3_1,1,var)*((n-1)^2/n^2)))

save.image("GHH_SE.RData")

set.seed(100)

rq.1 <-  zig.rqres(y=shark.type,x=shark$SoakTime,z=shark$SoakTime,
                   b0=coef(ZIG.shark)[1],b1=coef(ZIG.shark)[2],
                   a0=coef(ZIG.shark)[3],a1=coef(ZIG.shark)[4])

rq.2 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 2)

rq.3 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 3)

df=data.frame(rq.1,rq.2,rq.3)

ggplot(df, aes(sample = rq.1))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("GHH Randomized Quantile Residuals (ZIG Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.2))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("GHH Randomized Quantile Residuals (ZISCMP(m=2) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.3))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("GHH Randomized Quantile Residuals (ZISCMP(m=3) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)



###################################################################
###################################################################
### SLE Analysis
###################################################################
###################################################################

shark.type <- shark[,12]
n <- length(shark.type)
table(shark.type)

P.shark <- glm(shark.type~shark$SoakTime,family=poisson)
NB.shark <- glm.nb(shark.type~shark$SoakTime)

ZIP.shark <- zeroinfl(shark.type~shark$SoakTime,dist="poisson")
ZIG.shark <- zeroinfl(shark.type~shark$SoakTime,dist="geometric")
ZINB.shark <- zeroinfl(shark.type~shark$SoakTime,dist="negbin")

ZIGP.shark <- est.zigp2(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)

ZICMP.shark <- fit.zicmp.reg2(y=shark.type,X=cbind(1,shark$SoakTime),S=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta.init=c(-4.296546728,0.002949773),gamma.init=0.5816472,zeta.init=c(-1.482435057,0.001814978),max=250,optim.control=list(maxit=2000))

ZIsCMP2.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

ZIsCMP3.shark <- ziscompreg.ll(y=shark.type,X=cbind(1,shark$SoakTime),G=cbind(rep(1,length(shark$SoakTime))),W=cbind(1,shark$SoakTime),beta0=ZIsCMP2.shark$Parameters[1:2,1],gamma0=log(ZIsCMP2.shark$Parameters[3,1]),xi0=ZIsCMP2.shark$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))

P.ll <- logLik(P.shark)
NB.ll <- logLik(NB.shark)
ZIP.ll <- logLik(ZIP.shark)
ZIG.ll <- logLik(ZIG.shark)
ZINB.ll <- logLik(ZINB.shark)
ZIGP.ll <- tail(ZIGP.shark,1)
ZICMP.ll <- ZICMP.shark$loglik
ZIsCMP2.ll <- ZIsCMP2.shark[[2]]
ZIsCMP3.ll <- ZIsCMP3.shark[[2]]

P.pars <- length(coef(P.shark))
NB.pars <- P.pars + 1
ZIP.pars <- length(coef(ZIP.shark))
ZIG.pars <- ZIP.pars
ZINB.pars <- ZIP.pars+1
ZIGP.pars <- length(ZIGP.shark)-1
ZICMP.pars <- ZIGP.pars
ZIsCMP2.pars <- ZIGP.pars
ZIsCMP3.pars <- ZIGP.pars

P.AIC <- -2*P.ll+2*P.pars
NB.AIC <- -2*NB.ll+2*NB.pars
ZIP.AIC <- -2*ZIP.ll+2*ZIP.pars
ZIG.AIC <- -2*ZIG.ll+2*ZIG.pars
ZINB.AIC <- -2*ZINB.ll+2*ZINB.pars
ZIGP.AIC <- -2*ZIGP.ll+2*ZIGP.pars
ZICMP.AIC <- -2*ZICMP.ll+2*ZICMP.pars
ZIsCMP2.AIC <- -2*ZIsCMP2.ll+2*ZIsCMP2.pars
ZIsCMP3.AIC <- -2*ZIsCMP3.ll+2*ZIsCMP3.pars

P.BIC <- -2*P.ll+log(n)*P.pars
NB.BIC <- -2*NB.ll+log(n)*NB.pars
ZIP.BIC <- -2*ZIP.ll+log(n)*ZIP.pars
ZIG.BIC <- -2*ZIG.ll+log(n)*ZIG.pars
ZINB.BIC <- -2*ZINB.ll+log(n)*ZINB.pars
ZIGP.BIC <- -2*ZIGP.ll+log(n)*ZIGP.pars
ZICMP.BIC <- -2*ZICMP.ll+log(n)*ZICMP.pars
ZIsCMP2.BIC <- -2*ZIsCMP2.ll+log(n)*ZIsCMP2.pars
ZIsCMP3.BIC <- -2*ZIsCMP3.ll+log(n)*ZIsCMP3.pars

R2.P <- 1-(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(P.ll)+P.pars+.5)/(sum(dpois(shark.type,shark.type,log=TRUE))-as.numeric(logLik(glm(shark.type~1,family=poisson))))
R2.NB <- 1-(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(NB.ll)+NB.pars+.5)/(sum(dnbinom(shark.type,size=NB.shark$theta,mu=shark.type,log=TRUE))-as.numeric(logLik(glm.nb(shark.type~1))))

d.zip <- function(x, lambda , pi) {
  ifelse(x == 0, pi + (1-pi)*dpois(0,lambda = lambda), (1-pi)*dpois(x,lambda = lambda))
}
ll.zip.sat <- sum(log(d.zip(shark.type,shark.type,as.numeric(shark.type == 0))))

d.zinb <- function(x,size,mu,pi) {
  ifelse(x == 0, pi + (1-pi)*dnbinom(0,size = size, mu = mu), (1-pi)*dnbinom(x,size = size, mu = mu))
}
ll.zinb.sat <- sum(log(d.zinb(shark.type,size = NB.shark$theta , mu = shark.type, pi = as.numeric(shark.type == 0))))

d.geom <- function(x,prob,pi) {
  ifelse(x == 0, pi + (1-pi)*dgeom(0,prob = prob), (1-pi)*dgeom(x,prob = prob))
}
ll.zig.sat <- sum(log(d.geom(shark.type,prob=1/(1+mean(shark.type)), pi = as.numeric(shark.type == 0))))

d.genpoi <- function(x,mu,phi,pi) {
  ifelse(x == 0, pi + (1-pi)*dzigp(0,mu=mu,phi=phi,omega=0), (1-pi)*dzigp(x,mu=mu,phi=phi,omega=0))
}
ll.zigp.sat <- sum(sapply(1:n,function(i) log(d.genpoi(shark.type[i],mu=shark.type[i], phi=1+exp(ZIGP.shark[3]), pi = as.numeric(shark.type[i] == 0)))))

d.cmp <- function(x,lambda,nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dcom(x,lambda=lambda,nu=nu))
}
ll.zicmp.sat <- sum(sapply(1:n,function(i) log(d.cmp(shark.type[i],lambda=shark.type[i], nu=exp(ZICMP.shark$theta.hat$gamma), pi = as.numeric(shark.type[i] == 0)))))

d.scmp <- function(x,NuOfVar,Lambda,Nu,pi) {
  ifelse(x == 0, pi + (1-pi)*0, (1-pi)*dSCMP(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu))
}
ll.ziscmp2.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=2,Lambda=shark.type[i], Nu=ZIsCMP2.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))
ll.ziscmp3.sat <- sum(sapply(1:n,function(i) log(d.scmp(shark.type[i],NuOfVar=3,Lambda=shark.type[i], Nu=ZIsCMP3.shark$Parameters[3,1], pi = as.numeric(shark.type[i] == 0)))))

R2.ZIP <- 1-((ll.zip.sat-as.numeric(ZIP.ll)+ZIP.pars+.5)/(ll.zip.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="poisson")))))
R2.ZINB <- 1-((ll.zinb.sat-as.numeric(ZINB.ll)+ZINB.pars+.5)/(ll.zinb.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="negbin")))))
R2.ZIG <- 1-((ll.zig.sat-as.numeric(ZIG.ll)+ZIG.pars+.5)/(ll.zig.sat-as.numeric(logLik(zeroinfl(shark.type~1,dist="geometric")))))
R2.ZIGP <- 1-((ll.zigp.sat-as.numeric(ZIGP.ll)+ZIGP.pars+.5)/(ll.zigp.sat-tail(est.zigp2(Yin=shark.type,fm.X=~1,fm.W=~1,fm.Z=~1,init=T,tex=FALSE),1)))
R2.ZICMP <- 1-((ll.zicmp.sat-as.numeric(ZICMP.ll)+ZICMP.pars+.5)/(ll.zicmp.sat-fit.zicmp.reg2(y=shark.type,X=cbind(rep(1,n)),S=cbind(rep(1,n)),W=cbind(rep(1,n)),beta.init=1,gamma.init=1,zeta.init=1,max=250,optim.control=list(maxit=2000))$loglik))
R2.ZIsCMP2 <- 1-((ll.ziscmp2.sat-as.numeric(ZIsCMP2.ll)+ZIsCMP2.pars+.5)/(ll.ziscmp2.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))
R2.ZIsCMP3 <- 1-((ll.ziscmp3.sat-as.numeric(ZIsCMP3.ll)+ZIsCMP3.pars+.5)/(ll.ziscmp3.sat-ziscompreg.ll(y=shark.type,X=cbind(rep(1,n)),G=cbind(rep(1,n)),W=cbind(rep(1,n)),beta0=1,gamma0=1,xi0=1,NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))[[2]]))

sum.out <- data.frame(ll=c(P.ll,NB.ll,ZIP.ll,ZIG.ll,ZINB.ll,ZIGP.ll,ZICMP.ll,ZIsCMP2.ll,ZIsCMP3.ll),
                      AIC=c(P.AIC,NB.AIC,ZIP.AIC,ZIG.AIC,ZINB.AIC,ZIGP.AIC,ZICMP.AIC,ZIsCMP2.AIC,ZIsCMP3.AIC),
                      BIC=c(P.BIC,NB.BIC,ZIP.BIC,ZIG.BIC,ZINB.BIC,ZIGP.BIC,ZICMP.BIC,ZIsCMP2.BIC,ZIsCMP3.BIC),
                      R2.adj=c(R2.P,R2.NB,R2.ZIP,R2.ZIG,R2.ZINB,R2.ZIGP,R2.ZICMP,R2.ZIsCMP2,R2.ZIsCMP3))
rownames(sum.out) <- c("Poi","NegBin","ZIP","ZIG","ZINB","ZIGP","ZICMP","ZIsCMP2","ZIsCMP3")

sum.out

save.image("SLE.RData")

ZICMP.JK <- vector("list",n)
ZIsCMP2.JK <- vector("list",n)
ZIsCMP3.JK <- vector("list",n)
for(i in 1:n){
  ZICMP.JK[[i]] <- fit.zicmp.reg2(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),S=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta.init=c(-4.296546728,0.002949773),gamma.init=0.5816472,zeta.init=c(-1.482435057,0.001814978),max=250,optim.control=list(maxit=2000))
  ZIsCMP2.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZICMP.shark$theta.hat$beta,gamma0=exp(ZICMP.shark$theta.hat$gamma),xi0=ZICMP.shark$theta.hat$zeta,NuOfVar=2, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1.e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  ZIsCMP3.JK[[i]] <- ziscompreg.ll(y=shark.type[-i],X=cbind(1,shark$SoakTime[-i]),G=cbind(rep(1,length(shark$SoakTime[-i]))),W=cbind(1,shark$SoakTime[-i]),beta0=ZIsCMP2.shark$Parameters[1:2,1],gamma0=log(ZIsCMP2.shark$Parameters[3,1]),xi0=ZIsCMP2.shark$Parameters[4:5,1],NuOfVar=3, method="nlminb",hessian=FALSE,alg.options=list(trace = 1, abs.tol = 1e-6, rel.tol = 1e-6, xf.tol=1e-6, step.min=.0001, step.max=1))
  print(i)
}


NB.JK <- rep(NA,n)
ZIG.JK <- NULL
for(i in 1:n){
  NB.JK[i] <- log(glm.nb(shark.type[-i]~shark$SoakTime[-i])$theta)
  ZIG.JK <- cbind(ZIG.JK,unlist(coef(zeroinfl(shark.type[-i]~shark$SoakTime[-i],dist="geometric"))))
  print(i)
}


cbind(coef(P.shark),sqrt(diag(summary(P.shark)$cov.scaled)))
#cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),summary(NB.shark)$SE.theta^2/summary(NB.shark)$theta^2))
cbind(c(coef(NB.shark),log(NB.shark$theta)),c(sqrt(diag(summary(NB.shark)$cov.scaled)),sqrt(var(NB.JK)*((n-1)^2/n^2))))
cbind(unlist(coef(ZIP.shark)),sqrt(diag(ZIP.shark$vcov)))
#cbind(unlist(coef(ZIG.shark)),sqrt(diag(ZIG.shark$vcov)))
cbind(unlist(coef(ZIG.shark)),sqrt(apply(ZIG.JK,1,var)*((n-1)^2/n^2)))
cbind(c(unlist(coef(ZINB.shark)),log(ZINB.shark$theta)),c(sqrt(diag(ZINB.shark$vcov)),ZINB.shark$SE.logtheta))
est.zigp(Yin=shark.type,fm.X=~shark$SoakTime,fm.W=~1,fm.Z=~shark$SoakTime,init=T,tex=FALSE)
ZICMP_1 <- sapply(1:n,function(i) unlist(ZICMP.JK[[i]]$theta));cbind(apply(ZICMP_1,1,mean),sqrt(apply(ZICMP_1,1,var)*((n-1)^2/n^2)))
ZIsCMP2_1 <- sapply(1:n,function(i) unlist(ZIsCMP2.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP2_1,1,mean),sqrt(apply(ZIsCMP2_1,1,var)*((n-1)^2/n^2)))
ZIsCMP3_1 <- sapply(1:n,function(i) unlist(ZIsCMP3.JK[[i]]$Parameters[,1]));cbind(apply(ZIsCMP3_1,1,mean),sqrt(apply(ZIsCMP3_1,1,var)*((n-1)^2/n^2)))


save.image("SLE_SE.RData")

set.seed(100)

rq.1 <-  zig.rqres(y=shark.type,x=shark$SoakTime,z=shark$SoakTime,
                   b0=coef(ZIG.shark)[1],b1=coef(ZIG.shark)[2],
                   a0=coef(ZIG.shark)[3],a1=coef(ZIG.shark)[4])

rq.2 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 2)

rq.3 <- rqres.ziscmp.reg(y=shark.type,X=cbind(1,shark$SoakTime),W=cbind(1,shark$SoakTime),
                         beta0=ZIsCMP3.shark$Parameters[1:2,1],
                         gamma0=log(ZIsCMP3.shark$Parameters[3,1]),
                         xi0=ZIsCMP3.shark$Parameters[4:5,1],NuOfVar = 3)

df=data.frame(rq.1,rq.2,rq.3)

ggplot(df, aes(sample = rq.1))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SLE Randomized Quantile Residuals (ZIG Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.2))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SLE Randomized Quantile Residuals (ZISCMP(m=2) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

ggplot(df, aes(sample = rq.3))+stat_qq()+
  geom_abline(intercept = 0, slope = 1,col=2)+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  ggtitle("SLE Randomized Quantile Residuals (ZISCMP(m=3) Fit)")+
  theme(text = element_text(size = 13))+xlim(-4,4)+ylim(-4,4)

