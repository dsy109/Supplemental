library(mixtools)
library(mixR)
library(ggplot2)

GMM.fid <- function(x,w=NULL,k,M){
  n <- length(x)
  if(is.null(w)) w <- sample(1:k,n,rep=T)
  R.sig2.j <- function(x,w) sapply(1:max(w),function(i) (sum(w==i)-1)*var(x[w==i])/rchisq(1,sum(w==i)-1) ) 
  R.mu.j <- function(x,w,R.sig2.j) sapply(1:max(w),function(i) mean(x[w==i])-rnorm(1)*sqrt(R.sig2.j[i]/sum(w==i))) 
  all.R.pi <- all.R.sig2 <- all.R.mu <- NULL
  for(h in 1:M){
    i.new <- sample(1:n,1)
    w.star <- w
    w.star[i.new] <- sample(c(1:k)[-w[i.new]],1)
    a <- w[i.new]
    b <- w.star[i.new]
    
    A <- x[w==a]; n.a <- length(A); s2.a <- var(A)
    A.star <- x[w.star==a]; n.a.star <- length(A.star); s2.a.star <- var(A.star)
    B <- x[w==b]; n.b <- length(B); s2.b <- var(B)
    B.star <- x[w.star==b]; n.b.star <- length(B.star); s2.b.star <- var(B.star)
    
    r.num <- log(n.b-1)+lgamma((n.a.star-1)/2)+lgamma((n.b.star-1)/2)+((n.a-2)/2)*(log(n.a-1)+log(s2.a))+((n.b-2)/2)*(log(n.b-1)+log(s2.b))
    r.den <- log(n.b.star-1)+lgamma((n.a-1)/2)+lgamma((n.b-1)/2)+((n.a.star-2)/2)*(log(n.a.star-1)+log(s2.a.star))+((n.b.star-2)/2)*(log(n.b.star-1)+log(s2.b.star))
    r <- exp(r.num-r.den)
    if(is.na(r)) r <- 1
    accept <- rbinom(1,1,min(r,1))
    if(accept & min(table(w.star))>1) w <- w.star
    
    r.j <- cumsum(sapply(1:max(w),function(i) sum(w==i)))
    U <- c(0,sort(runif(n)),1)
    D.j <- runif(k)
    R.pi.j <- U[r.j[1]+1] + D.j[1]*(U[r.j[1]+2]-U[r.j[1]+1])
    if(k > 2) for(j in 2:(k-1)) R.pi.j <- c(R.pi.j,U[r.j[j]+1] + D.j[j]*(U[r.j[j]+2]-U[r.j[j]+1])-tail(R.pi.j,1))
    R.pi.j <- c(R.pi.j,1-sum(R.pi.j))
    
    Rsig2 <- R.sig2.j(x,w)
    Rmu <- R.mu.j(x,w,Rsig2)
    
    if(any(is.na(Rsig2))) stop()
    
    all.R.pi <- rbind(all.R.pi,R.pi.j)
    all.R.mu <- rbind(all.R.mu,Rmu)
    all.R.sig2 <- rbind(all.R.sig2,Rsig2)
    print(h)
  }
  ind <- t(sapply(1:M,function(i) order(all.R.mu[i,])))
  all.R.mu <- t(sapply(1:M,function(i) all.R.mu[i,][ind[i,]]))
  all.R.sig2 <- t(sapply(1:M,function(i) all.R.sig2[i,][ind[i,]]))
  all.R.pi <- t(sapply(1:M,function(i) all.R.pi[i,][ind[i,]]))
  out <- list(R.pi=all.R.pi,R.mu=all.R.mu,R.sig2=all.R.sig2)
  out
}

#Simulated Data
burn.in <- 5000
set.seed(1)
x <- c(rnorm(50),rnorm(25,6),rnorm(25,12))[sample(1:100,size=100)]
out.fid <- GMM.fid(x,k=3,M=10000)
out.EM <- normalmixEM(x,k=3)
sapply(1:3,function(i) apply(out.fid[[i]][-c(1:burn.in),],2,mean))
out.EM[2:4]

#Hidalgo Stamp Data
data(Stamp)
burn.in <- 50000
set.seed(10)
out.fid.stamp <- GMM.fid(Stamp,k=4,M=100000)
set.seed(10)
out.EM.stamp <- normalmixEM(Stamp,k=4,eps=1e-10,maxit=10000)
sapply(1:3,function(i) apply(out.fid.stamp[[i]][-c(1:burn.in),],2,mean))
out.EM.stamp[2:4]









burn.in <- 50000
ind <- t(sapply(1:M,function(i) order(all.R.mu[i,])))
all.R.mu.red <- t(sapply((burn.in+1):M,function(i) all.R.mu[i,][ind[i,]]))
all.R.sig2.red <- t(sapply((burn.in+1):M,function(i) all.R.sig2[i,][ind[i,]]))
all.R.pi.red <- t(sapply((burn.in+1):M,function(i) all.R.pi[i,][ind[i,]]))

apply(all.R.mu.red,2,mean)
apply(all.R.sig2.red,2,mean)
apply(all.R.pi.red,2,mean)

ii <- 1
plot(all.R.mu.red[,ii],type="l")
plot(all.R.sig2.red[,ii],type="l")
plot(all.R.pi.red[,ii],type="l")





