#################################################################################
###############################   Function List   ###############################
#################################################################################

## nb.mixEM()
## pois.mix()
## cmpmixEM() 
## --------------------------------------------------------
## Functions used in cmpmixEM()
# Z(lambda,nu)
# mu.fun(lambda,nu)
# mean.fun(lambda,nu)
# lambda.fun(mu,nu)
# pmf(mu,nu)
# var.fun(lambda,nu)
# mean_logfactorialy.fun(lambda,nu)
# mean_ylogfactorialy.fun(lambda,nu)
# nu.fun(x,mu)
# --------------------------------------------------------

#################################################################################
# nb.mixEM() 
# Fit mixture of Negative Binomials
#################################################################################

nb.mixEM <- function(y,k=k,theta.t=NULL,eps=1e-6,maxit=1000){
  
  library(MASS)
  
  ## initial values
  n <- length(y)
  if(is.null(theta.t)){
    
    pi.t <- rep(1,k)/k
    
    y.label <- kmeans(y,k)$cluster
    yi <- lapply(1:k, function(i) assign(paste0("y",i),y[which(y.label==i)]) ) 
    glm.out <- lapply(1:k, function(i) glm.nb(yi[[i]]~1))
    
    theta.t <- list(pi.t, sapply(1:k,function(i) exp(glm.out[[i]]$coefficients)),
                    sapply(1:k,function(i) glm.out[[i]]$theta))
  }
  pi.t <- theta.t[[1]]
  mu.t <- theta.t[[2]]
  phi.t <- theta.t[[3]]
  
  ## observed loglikelihood 
  x.k <- matrix(NA,nrow=n, ncol=k)
  for (i in 1:n){
    for (j in 1:k) {
      x.k[i,j] <- dnbinom(y[i],mu=mu.t[j],size=phi.t[j])
    }
  }
  ll <- sum(log(rowSums(t(t(x.k)*pi.t))))
  
  
  ## update  
  iter <- 0
  diff <- 1
  
  ## output
  out <- c(pi.t,mu.t,phi.t,ll,iter) 
  
  while((iter <= maxit) & (diff > eps)){
    
    print(iter)
    iter <- iter+1
    
    ## glm.nb fit with weights
    z <- sapply(1:k, function(i) pi.t[i]*dnbinom(y,size=phi.t[i],mu=mu.t[i]))
    z <- z/apply(z,1,sum)
    glm.out <- lapply(1:k,function(i) glm.nb(y~1,weights=z[,i]))
    theta.t[[1]] <- rbind(theta.t[[1]],apply(z,2,mean))
    theta.t[[2]] <- rbind(theta.t[[2]],sapply(1:k,function(i) exp(glm.out[[i]]$coefficients)))
    theta.t[[3]] <- rbind(theta.t[[3]],sapply(1:k,function(i) glm.out[[i]]$theta))
    
    ## update parameters
    pi.t <- tail(theta.t[[1]],1)
    mu.t <- tail(theta.t[[2]],1)
    phi.t <- tail(theta.t[[3]],1)
    
    ## calculate loglikelihood
    x.k <- matrix(NA,nrow=n, ncol=k)
    for (i in 1:n){
      for (j in 1:k) {
        x.k[i,j] <- dnbinom(y[i],mu=mu.t[j],size=phi.t[j])
      }
    }
    ll.new <- sum(log(rowSums(t(t(x.k)*colMeans(z)))))
    
    diff <- abs(ll-ll.new) ; print(diff)
    ll <- ll.new
    
    ## output
    out <- rbind(out,c(pi.t,mu.t,phi.t,ll.new,iter))
    
  }
  
  colnames(out) <- c(paste("pi",1:k,sep=""), 
                     paste("mu",1:k,sep=""),
                     paste("phi",1:k,sep=""),"ll", "iter")
  
  return(out)
}

#################################################################################



#################################################################################
# pois.mix() 
# Fit mixture of Poissons 
#################################################################################

pois.mix <- function(y,k){
  
  library(flexmix)
  
  df <- data.frame(x)
  
  pois.mod <- flexmix(x~1, data=df, k=k, model=FLXMRglm(family="poisson"))
  
  ## output
  out <- c(summary(pois.mod)@comptab$prior,
           exp(as.vector(parameters(pois.mod))),
           as.numeric(summary(pois.mod)@logLik))
  
  names(out) <- c(paste("pi",1:k,sep=""), 
                  paste("mu",1:k,sep=""), "ll")
  
  return(out)
  
}

#################################################################################



#################################################################################
# cmpmixEM() 
# Fit mixture of Mean-Parameterized CMPs
#################################################################################

cmpmixEM <- function(x,k=k, 
                     mu=NULL,nu=NULL,pi=NULL,
                     nu.star=seq(0.1,10,0.1),
                     eps=1e-6,maxit=1000){
  ## initial data
  x <- sort(x)
  n <- length(x)
  
  ## initial parameters
  if (is.null(mu)) mu <- sort(as.vector(kmeans(x,k)$centers)) # k means
  if (is.null(nu)) nu <- rep(1,k)
  if (is.null(pi)) pi <- rep(1/k,k)
  
  
  ## initial observations in columns
  x.k <- matrix(nrow=n, ncol=k)
  for (i in 1:k) {
    x.k[,i] <- pmf(x,mu[i],nu[i])
  }
  
  ## observed loglikelihood  
  obs.ll <- sum(log(rowSums(t(t(x.k)*pi))))
  
  ## iteration to update the parameters
  iter <- 0
  dif <- 1
  
  # output
  out <- c(pi,mu,nu,obs.ll,iter)
  
  while(iter < maxit && dif > eps){
    
    ## update parameters
    z.t <- t(t(x.k)*pi) / rowSums(t(t(x.k)*pi))
    pi <- colMeans(z.t);pi
    
    for (i in 1:k) {
      mu[i] <- weighted.mean(x,z.t[,i])
      nu[i] <- nu.fun(x,z.t[,i],mu[i],nu=nu.star)
    }
    
    ## update observations
    for (i in 1:k) {
      x.k[,i] <- pmf(x,mu[i],nu[i])
    }
    
    new.obs.ll <- sum(log(rowSums(t(t(x.k)*pi))))
    dif <- abs(new.obs.ll-obs.ll)
    print(dif)
    
    obs.ll <- new.obs.ll
    print(iter)
    iter <- iter+1
    
    # output dataframe
    out <- rbind(out,c(pi,mu,nu,new.obs.ll,iter))
  }
  
  colnames(out) <- c(paste("pi",1:k,sep=""), 
                     paste("mu",1:k,sep=""),
                     paste("nu",1:k,sep=""),"ll","iter")
  return(out)
}


#################################################################################
# Z(lambda,nu)
# function to approximate the normalizing constant Z
#################################################################################

Z <- function(lambda, nu, maxj=150, eps=1e-6){
  j <- 0
  Z.value <- term <- exp(j*log(lambda)-nu*lgamma(j+1))
  while(j < maxj && term > eps ){
    j <- j+1
    term <- exp(j*log(lambda)-nu*lgamma(j+1))
    Z.value <- Z.value + term
  }
  
  return(Z.value)
}


#################################################################################
# mu.fun(lambda,nu)
# Solve for mu for given lambda and nu
#################################################################################
library(optimr)

mu.fun <- function(lambda,nu,ylim=150){
  y <- 0:ylim
  fun.s <- function(lambda,mu,nu) abs(sum((y-mu)*exp(y*log(lambda)-nu*lgamma(y+1))))
  mu <- hjn(0,fun.s,lower=0,upper=100,lambda=lambda,nu=nu)$par
  return(mu)
}

#################################################################################
# mean.fun(lambda,nu)
# solve for the mu i.e. the expectation for given lambda and nu
#################################################################################

mean.fun <- function(lambda,nu,maxy=100,eps=1e-6){
  y <- 1
  sum.value <- term <- y*exp(y*log(lambda)-nu*lgamma(y+1))
  while(y < maxy && term > eps){
    y <- y+1
    term <- y*exp(y*log(lambda)-nu*lgamma(y+1))
    sum.value <- sum.value + term
  }
  return(sum.value/Z(lambda,nu))
}

#################################################################################
# lambda.fun(mu,nu)
# Solve for the lambda for given mu and nu
#################################################################################

library(optimr)

lambda.fun <- function(mu,nu,ylim=150){
  y <- seq(0,ylim,by=1)
  fun.s <- function(lambda,mu,nu) abs(mu.fun(lambda,nu)-mu)
  lambda <- hjn(5,fun.s,lower=0.1,upper=1000,mu=mu,nu=nu)$par
  return(lambda)
}


#################################################################################
# pmf(mu,nu)
# Evaluate pmf given mu and nu
#################################################################################

pmf <- function(y,mu,nu){
  lambda <- lambda.fun(mu,nu)
  return(exp(y*log(lambda)-nu*lgamma(y+1)) / Z(lambda,nu))
}


#################################################################################
# The following three functions are used to solve nu 
#################################################################################

#################################################################################
# var.fun(lambda,nu)
# variance given lambda and nu
#################################################################################

var.fun <- function(lambda,nu,maxy=150,eps=1e-6){
  y <- 1
  sum.value <- term <- (y^2)*exp(y*log(lambda)-nu*lgamma(y+1))
  while(y < maxy && term > eps){
    y <- y+1
    term <- (y^2)*exp(y*log(lambda)-nu*lgamma(y+1))
    sum.value <- sum.value + term
  }
  mu <- mean.fun(lambda,nu)
  return(sum.value/Z(lambda,nu) - mu^2)
}

#################################################################################
# mean_logfactorialy.fun(lambda,nu)
# Find E[log(Y!)] given lambda and nu
#################################################################################

mean_logfactorialy.fun <- function(lambda,nu,maxy=150,eps=1e-6){
  y <- 2
  sum.value <- term <- lgamma(y+1)*exp(y*log(lambda)-nu*lgamma(y+1))
  while(y < maxy && term > eps){
    y <- y+1
    term <- lgamma(y+1)*exp(y*log(lambda)-nu*lgamma(y+1))
    sum.value <- sum.value + term
  }
  return(sum.value/Z(lambda,nu))
}

#################################################################################
# mean_ylogfactorialy.fun(lambda,nu)
# Find E[Ylog(Y!)] given lambda and nu
#################################################################################

mean_ylogfactorialy.fun <- function(lambda,nu,maxy=150,eps=1e-6){
  y <- 2
  sum.value <- term <- y*lgamma(y+1)*exp(y*log(lambda)-nu*lgamma(y+1))
  while(y < maxy && term > eps){
    y <- y+1
    term <- y*lgamma(y+1)*exp(y*log(lambda)-nu*lgamma(y+1))
    sum.value <- sum.value + term
  }
  return(sum.value/Z(lambda,nu))
}

#################################################################################



#################################################################################
# nu.fun(x,mu) 
# Solve the MLE of nu by searching  
#################################################################################

nu.fun <- function(y,z,mu,nu=seq(0.1,10,0.1)){
  lambda <- c()
  f <- c()
  for (i in 1:length(nu)) {
    lambda[i] <- lambda.fun(mu,nu[i])
    f[i] <- ((mean_ylogfactorialy.fun(lambda[i],nu[i])-mu*mean_logfactorialy.fun(lambda[i],nu[i])) * ((t(z) %*% y)-mu*sum(z)) / var.fun(lambda[i],nu[i]) - t(z) %*% (lgamma(y+1)) + sum(z)* mean_logfactorialy.fun(lambda[i],nu[i]))
  }
  return(nu[which.min(abs(f))])
}

#################################################################################





