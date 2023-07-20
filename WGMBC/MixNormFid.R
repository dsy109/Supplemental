library(mixtools)
library(mixR)
library(ggplot2)

GMM.fid <- function(x,w=NULL,k,M,min.nj=2){
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
    if(accept & min(table(w.star))>=min.nj) w <- w.star
    
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
out.fid2 <- sapply(1:3,function(i) apply(out.fid[[i]][-c(1:burn.in),],2,mean))
out.fid2
out.EM[2:4]

ggplot(data.frame(i=5001:10000,out.fid[[1]][-c(1:burn.in),]), aes(x=i,y=X1)) + geom_line(col="#56B4E9") + labs(x="Post-Burn-In Iteration",y=expression(R[pi[1]])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", pi[1]))) +   theme(text = element_text(size = 25))
ggplot(data.frame(i=5001:10000,out.fid[[1]][-c(1:burn.in),]), aes(x=i,y=X2)) + geom_line(col="#56B4E9") + labs(x="Post-Burn-In Iteration",y=expression(R[pi[2]])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", pi[2]))) +   theme(text = element_text(size = 25))
ggplot(data.frame(i=5001:10000,out.fid[[1]][-c(1:burn.in),]), aes(x=i,y=X3)) + geom_line(col="#56B4E9") + labs(x="Post-Burn-In Iteration",y=expression(R[pi[3]])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", pi[3]))) +   theme(text = element_text(size = 25))

ggplot(data.frame(i=5001:10000,out.fid[[2]][-c(1:burn.in),]), aes(x=i,y=X1)) + geom_line(col="#56E98B") + labs(x="Post-Burn-In Iteration",y=expression(R[mu[1]])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", mu[1]))) +   theme(text = element_text(size = 25))
ggplot(data.frame(i=5001:10000,out.fid[[2]][-c(1:burn.in),]), aes(x=i,y=X2)) + geom_line(col="#56E98B") + labs(x="Post-Burn-In Iteration",y=expression(R[mu[2]])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", mu[2]))) +   theme(text = element_text(size = 25))
ggplot(data.frame(i=5001:10000,out.fid[[2]][-c(1:burn.in),]), aes(x=i,y=X3)) + geom_line(col="#56E98B") + labs(x="Post-Burn-In Iteration",y=expression(R[mu[3]])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", mu[3]))) +   theme(text = element_text(size = 25))

ggplot(data.frame(i=5001:10000,out.fid[[3]][-c(1:burn.in),]), aes(x=i,y=X1)) + geom_line(col="#E956B4") + labs(x="Post-Burn-In Iteration",y=expression(R[sigma[1]^2])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", sigma[1]^2))) +   theme(text = element_text(size = 25))
ggplot(data.frame(i=5001:10000,out.fid[[3]][-c(1:burn.in),]), aes(x=i,y=X2)) + geom_line(col="#E956B4") + labs(x="Post-Burn-In Iteration",y=expression(R[sigma[2]^2])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", sigma[2]^2))) +   theme(text = element_text(size = 25))
ggplot(data.frame(i=5001:10000,out.fid[[3]][-c(1:burn.in),]), aes(x=i,y=X3)) + geom_line(col="#E956B4") + labs(x="Post-Burn-In Iteration",y=expression(R[sigma[3]^2])) +
  ggtitle(expression(paste("Trace Plot of Generalized Fiducial Distribution for ", sigma[3]^2))) +   theme(text = element_text(size = 25))


#Hidalgo Stamp Data
data(Stamp)
burn.in <- 50000
set.seed(10)
out.fid.stamp <- GMM.fid(Stamp,k=4,M=100000)
set.seed(10)
out.EM.stamp <- normalmixEM(Stamp,k=4,eps=1e-10,maxit=10000)
out.fid.stamp2 <- sapply(1:3,function(i) apply(out.fid.stamp[[i]][-c(1:burn.in),],2,mean))
out.fid.stamp2
out.EM.stamp[2:4]

Stamps <- data.frame(Stamp)
ggplot(Stamps, aes(x=Stamp)) + geom_histogram(color="#FFFFFF") + ggtitle("Hidalgo Stamp Data") + theme(text = element_text(size = 15)) + xlim(c(0.055,0.135))
X <- seq(0.055,0.135,length=200)
DF <- data.frame(x=X,y1=out.fid.stamp2[1,1]*dnorm(X,out.fid.stamp2[1,2],sqrt(out.fid.stamp2[1,3])),y2=out.fid.stamp2[2,1]*dnorm(X,out.fid.stamp2[2,2],sqrt(out.fid.stamp2[2,3])),
                 y3=out.fid.stamp2[3,1]*dnorm(X,out.fid.stamp2[3,2],sqrt(out.fid.stamp2[3,3])),y4=out.fid.stamp2[4,1]*dnorm(X,out.fid.stamp2[4,2],sqrt(out.fid.stamp2[4,3])))
ggplot(Stamps, aes(x=Stamp)) + geom_histogram(aes(y=..density..),color="#FFFFFF") + 
  layer(geom="area",data=DF,stat="identity",position="identity",mapping=aes(x=x,y=y1),params = list(fill = "#56B4E9", alpha = 0.5)) + 
  layer(geom="area",data=DF,stat="identity",position="identity",mapping=aes(x=x,y=y2),params = list(fill = "#56E98B", alpha = 0.5)) + 
  layer(geom="area",data=DF,stat="identity",position="identity",mapping=aes(x=x,y=y3),params = list(fill = "#E956B4", alpha = 0.5)) + 
  layer(geom="area",data=DF,stat="identity",position="identity",mapping=aes(x=x,y=y4),params = list(fill = "#E98B56", alpha = 0.5)) + 
  geom_line(data=DF,aes(x=x,y=y1+y2+y3+y4),color="black",lwd=.18,alpha=.8) + ggtitle("Fiducial Solution") +
  theme(text = element_text(size = 15))
DF.EM <- data.frame(x=X,y1=out.EM.stamp$lambda[1]*dnorm(X,out.EM.stamp$mu[1],out.EM.stamp$sigma[1]),y2=out.EM.stamp$lambda[2]*dnorm(X,out.EM.stamp$mu[2],out.EM.stamp$sigma[2]),
                 y3=out.EM.stamp$lambda[3]*dnorm(X,out.EM.stamp$mu[3],out.EM.stamp$sigma[3]),y4=out.EM.stamp$lambda[4]*dnorm(X,out.EM.stamp$mu[4],out.EM.stamp$sigma[4]))
ggplot(Stamps, aes(x=Stamp)) + geom_histogram(aes(y=..density..),color="#FFFFFF") + 
  layer(geom="area",data=DF.EM,stat="identity",position="identity",mapping=aes(x=x,y=y1),params = list(fill = "#56B4E9", alpha = 0.5)) + 
  layer(geom="area",data=DF.EM,stat="identity",position="identity",mapping=aes(x=x,y=y2),params = list(fill = "#56E98B", alpha = 0.5)) + 
  layer(geom="area",data=DF.EM,stat="identity",position="identity",mapping=aes(x=x,y=y3),params = list(fill = "#E956B4", alpha = 0.5)) + 
  layer(geom="area",data=DF.EM,stat="identity",position="identity",mapping=aes(x=x,y=y4),params = list(fill = "#E98B56", alpha = 0.5)) + 
  geom_line(data=DF.EM,aes(x=x,y=y1+y2+y3+y4),color="black",lwd=.18,alpha=.8) + ggtitle("EM Solution") +
  theme(text = element_text(size = 15))






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





