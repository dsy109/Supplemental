#################################################################################
###############################   Functions List   ##############################
#################################################################################

## ------- Regression -------------------------------------

# cmp.mixEMReg() 
# pois.mixReg.bestfit()
# pois.mixReg()
# nb.mixEMReg() 

## ------- Uni-variate -----------------------------------

# nb.mixEM()
# pois.mix()
# cmp.mixEM() 

## ------- Functions used in cmpmixEM() ------------------
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
#####   Packages   ##############################################################
#################################################################################

library(optimr)  # function hjn()

library(MASS)    # for NB mixtures
library(flexmix) # for Poisson mixtures

library(optimr)  # for mean-parameterized CMP pmf

# for CMPs regression
library(Matrix)  # to make diagonal matrices
library(nloptr)  # to minimize -Q 


#################################################################################
# cmpmixEMReg() 
# Regress mixture of CMPs
#################################################################################

cmp.mixEMReg <- function(y,x,k=k,
                        beta=NULL,nu=NULL,Pi=NULL,
                        eps=1e-3,maxit=1000){
  
  
  ## create dataframe for Poisson mixtures
  df <- data.frame(y,x)
  
  ## covariates
  X <- cbind(1,x)
  
  ## number of observations
  n <- length(y)
  
  ## number of components k
  
  ## column number of beta's 
  q <- ncol(X)
  
  ## initial beta's from Poisson mixtures
  if (is.null(beta)) {
    out.pois <- pois.mixReg.bestfit(y,x,k,maxr=100)
    beta <- matrix(out.pois[c("beta01","beta11","beta02","beta12")],nrow=q,ncol=k)
  }
  
  ## initial Pi's and nu's
  if (is.null(nu)) nu <- rep(1,k)
  if (is.null(Pi)) Pi <- rep(1/k,k)
  
  ## initial observations
  x.beta <- X %*% beta
  mu0 <- exp(X %*% beta)
  lambda0 <- matrix(nrow=n, ncol=k)
  y.k <- matrix(nrow=n, ncol=k)
  for (i in 1:n) {
    for (j in 1:k) {
      lambda0[i,j] <- lambda.fun(mu0[i,j],nu=nu[j],ylim=150)
      y.k[i,j] <- pmf(y[i], mu=mu0[i,j], nu=nu[j])
    }
  }
  
  ## initial observed loglikelihood
  obs.ll <- sum(log(rowSums(t(t(y.k)*Pi)))) 
  
  ## define posterior probabilities
  z.t <- t(t(y.k)*Pi) / rowSums(t(t(y.k)*Pi)) 
  
  ## initial estimates summary 
  param <- c(c(beta),c(lambda0),nu)
  
  
  ######################################################  
  ## Define the functions for using in EM algorithm ##
  ######################################################
  
  #################################################
  ## objective function to maximize Q (i.e. minimize -Q)
  Q.f <- function(param, k=k, y,X,z.t){
    
    q <- ncol(X)
    n <- length(y)
    
    ## parameters
    beta <- param[1:(q*k)]
    beta <- matrix(beta, nrow=q, ncol=k)
    
    lambda <- param[(q*k+1):(q*k+n*k)]
    lambda <- matrix(lambda, nrow=n, ncol=k)
    
    nu <- param[(q*k+n*k+1):(q*k+n*k+k)]
    
    ## values in Q
    nu.lfactorial <- matrix(nrow=n, ncol=k)
    ZZ <- matrix(nrow=n, ncol=k)
    for (i in 1:n) {
      for (j in 1:k) {
        nu.lfactorial[i,j] <- nu[j] * lgamma(y[i]+1) 
        ZZ[i,j] <- Z(lambda[i,j],nu[j])
      }
    }
    
    Q <- sum(t(z.t) * log(Pi)) + sum(z.t * ( y * log(lambda)  - nu.lfactorial - log(ZZ))) 
    
    return(-Q) 
  }
  
  #################################################
  
  #################################################
  ## gradients of the objective function
  Q.g <- function(param, k=k, y,X,z.t){
    
    q <- ncol(X)
    n <- length(y)
    
    ## parameters
    beta <- param[1:(q*k)]
    beta <- matrix(beta, nrow=q, ncol=k)
    
    lambda <- param[(q*k+1):(q*k+n*k)]
    lambda <- matrix(lambda, nrow=n, ncol=k)
    
    nu <- param[(q*k+n*k+1):(q*k+n*k+k)]
    mu <- exp(X %*% beta)
    
    ## values in gradients of Q
    mean_logfacy <- matrix(nrow=n, ncol=k)
    for (i in 1:n) {
      for (j in 1:k) {
        mean_logfacy[i,j] <- mean_logfactorialy.fun(lambda[i,j],nu[j])
      }
    }
    
    beta.grad <- rep(0,q*k)
    lambda.grad <- (z.t * y - z.t * mu)/lambda
    nu.grad <- colSums(z.t * (-lgamma(y+1)) + z.t *  mean_logfacy)
    
    return(-c(beta.grad,c(lambda.grad),nu.grad)) 
  }
  
  #################################################
  
  #################################################
  ## equality constraint function
  mu.con <- function(param, k=k, y,X,z.t){
    
    q <- ncol(X)
    n <- length(y)
    
    ## parameters
    beta <- param[1:(q*k)]
    beta <- matrix(beta, nrow=q, ncol=k)
    
    lambda <- param[(q*k+1):(q*k+n*k)]
    lambda <- matrix(lambda, nrow=n, ncol=k)
    
    nu <- param[(q*k+n*k+1):(q*k+n*k+k)]
    
    mu <- matrix(nrow=n, ncol=k)
    for (i in 1:n) {
      for (j in 1:k) {
        mu[i,j] <- mean.fun(lambda[i,j],nu[j],maxy=100,eps=1e-6)
      }
    }
    
    x.beta <- X %*% beta
    
    return(c(exp(x.beta)-mu))
    
  }
  
  #################################################
  
  #################################################
  ## gradients of equality constraint function
  mu.con.g <- function(param, k=k, y,X,z.t){
    
    q <- ncol(X)
    n <- length(y)
    
    ## parameters
    beta <- param[1:(q*k)]
    beta <- matrix(beta, nrow=q, ncol=k)
    
    lambda <- param[(q*k+1):(q*k+n*k)]
    lambda <- matrix(lambda, nrow=n, ncol=k)
    
    nu <- param[(q*k+n*k+1):(q*k+n*k+k)]
    
    
    ## beta gradients of constraint
    mu <- exp(X %*% beta)
    
    grad.beta <- t(X) %*% diag(mu[,1])
    if (k > 1) {
      for (i in 2:k) {
        grad.beta <- bdiag(grad.beta, t(X) %*% diag(mu[,i]) )
      }
    }
    
    
    ## lambda gradients of constraint
    V <- matrix(nrow=n, ncol=k)
    for (i in 1:n) {
      for (j in 1:k) {
        V[i,j] <- var.fun(lambda[i,j],nu[j],maxy=150,eps=1e-6)
      }
    } 
    
    grad.lambda <- diag(V[,1]/lambda[,1])
    if (k > 1) {
      for (i in 2:k) {
        grad.lambda <- bdiag(grad.lambda, diag(V[,i]/lambda[,i]) )
      }
    }
    
    
    ## nu gradients of constraint
    grad.nu.m <- matrix(nrow=n, ncol=k)
    for (i in 1:n) {
      for (j in 1:k) {
        grad.nu.m[i,j] <- mean_ylogfactorialy.fun(lambda[i,j],nu[j],maxy=150,eps=1e-6) - mu[i,j] * mean_logfactorialy.fun(lambda[i,j],nu[j],maxy=150,eps=1e-6)
      }
    } 
    
    grad.nu <- grad.nu.m[,1]
    if (k > 1){
      for (i in 2:k) {
        grad.nu <- bdiag(grad.nu, grad.nu.m[,i] )
      }
    }
    
    
    return(as.matrix(cbind(t(grad.beta), grad.lambda, grad.nu)))
    
  }
  
  #################################################
  
  #################################################
  # nloptr package to solve for the estimates by minimizing -Q
  
  # algorithms applicable for case with equality constraint
  # agrm <- "NLOPT_GN_ISRES"
  # agrm <- "NLOPT_LD_AUGLAG"
  
  # agrm <- "NLOPT_LD_SLSQP"
  
  # my_opts <- list("algorithm"= agrm,"maxeval" = 1000,
  #                 "local_opts" = list("algorithm" = agrm, "xtol_rel" = 0.01))
  #################################################
  
  
  #################################################
  # define the upper and lower bounds for the parameters
  
  ## lower bounds
  beta_l <- rep(-Inf, q*k)
  lambda_l <- rep(0.1, n*k)
  nu_l <- rep(0.5, k)
  
  param_l <- c(beta_l, lambda_l, nu_l)
  
  ## upper bounds
  beta_u <- rep(Inf, q*k)
  lambda_u <- rep(200, n*k)
  nu_u <- rep(10, k)
  
  param_u <- c(beta_u, lambda_u, nu_u)
  
  #################################################
  
  
  
  ######################################################  
  ## EM algorithm ######################################
  ######################################################
  
  ## iteration starts
  iter <- 0
  dif <- 1
  
  ## output summary
  out <- c(iter, obs.ll, Pi, c(beta), nu)
  
  ## iteration for EM algorithm
  while(iter < maxit && dif > eps){
    
    ## nloptr to minimize -Q
    # library(nloptr)
    fit <- nloptr::nloptr(x0 = param,
                          eval_f = Q.f,
                          eval_grad_f = Q.g,
                          lb = param_l,
                          ub = param_u,
                          eval_g_eq = mu.con,
                          eval_jac_g_eq = mu.con.g,
                          opts = list("algorithm"= "NLOPT_LD_SLSQP","maxeval" = 1000,
                                      "local_opts" = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 0.001)),
                          k=k,y=y,X=X,z.t=z.t)
    
    ## update estimates
    param <- fit$solution
    beta <- param[1:(q*k)]
    beta <- matrix(beta, nrow=q, ncol=k)
    nu <- param[(q*k+n*k+1):(q*k+n*k+k)]
    
    ## update observations
    x.beta <- X %*% beta
    mu0 <- exp(X %*% beta)
    y.k <- matrix(nrow=n, ncol=k)
    for (i in 1:n) {
      for (j in 1:k) {
        y.k[i,j] <- pmf(y[i], mu=mu0[i,j], nu=nu[j])
      }
    }
    
    ## update loglikelihood  
    new.obs.ll <- sum(log(rowSums(t(t(y.k)*Pi)))) 
    
    ## update posterior probabilites
    z.t <- t(t(y.k)*Pi) / rowSums(t(t(y.k)*Pi)) 
    
    ## update mixing proportions
    Pi <- colMeans(z.t)
    
    ## print
    iter <- iter+1
    print(iter)
    dif <- abs(new.obs.ll-obs.ll)
    print(dif)
    
    obs.ll <- new.obs.ll
    
    ## iteration ends 
    
    ## output dataframe
    out <- rbind(out,c(iter, obs.ll, Pi, c(beta), nu))
  }
  
  colnames(out) <- c("iter", "ll", 
                     paste("Pi",1:k,sep=""), 
                     paste(rep(paste("beta",0:(q-1),sep=""),k),rep(1:k,each=q),sep=""),
                     paste("nu",1:k,sep=""))
  
  return(out)
  
}

#################################################################################

#################################################################################
# pois.mixReg.bestfit() 
# Find the best fit of Poisson mixture regression
#################################################################################

pois.mixReg.bestfit <- function(y,x,k,maxr=100){
  
  # library(flexmix)
  df <- data.frame(y,x)
  
  # list to store Poisson regression
  pois.mod <- list(maxr)
  
  # assign a random fit as the best fit
  pois.mod.bestfit <- flexmix(y~x, data=df, k=k, model=FLXMRglm(family="poisson"))
  
  # compare pois.mod.bestfit with each of the other pois.mod
  # return the one with maximum likelihood as the best fit
  r <- 1 
  while(r < maxr){
    
    print(paste("r =",r))
    
    pois.mod[[r]] <- flexmix(y~x, data=df, k=k, model=FLXMRglm(family="poisson"))
    
    if (pois.mod.bestfit@logLik < pois.mod[[r]]@logLik) {
        pois.mod.bestfit <- pois.mod[[r]]}
    
    r <- r+1
    
  }
  
  
  ## output
  out <- c(as.numeric(BIC(pois.mod.bestfit)),
           as.numeric(summary(pois.mod.bestfit)@logLik),
           summary(pois.mod.bestfit)@comptab$prior,
           as.vector(parameters(pois.mod.bestfit)))
  
  # some empty component may be dropped, new k retrieved
  k <- length(summary(pois.mod.bestfit)@comptab$prior)
  
  # q calculated
  X <- cbind(1,x); q <- ncol(X)
  
  names(out) <- c("BIC","ll", 
                  paste0("Pi",1:k,sep=""), 
                  paste(rep(paste("beta",0:(q-1),sep=""),k),rep(1:k,each=q),sep=""))
  
  return(out)
  
  # use sort.R script to sort the output
  
}



#################################################################################
# pois.mixReg() 
# Regress mixture of Poissons 
#################################################################################

pois.mixReg <- function(y,x,k){
  
  # library(flexmix)
  df <- data.frame(y,x)
  
  pois.mod <- flexmix(y~x, data=df, k=k, model=FLXMRglm(family="poisson"))
  
  X <- cbind(1,x)
  q <- ncol(X)
  
  ## output
  out <- c(as.numeric(BIC(pois.mod)),
           as.numeric(summary(pois.mod)@logLik),
           summary(pois.mod)@comptab$prior,
           as.vector(parameters(pois.mod)))
  
  names(out) <- c("BIC","ll", 
                  paste0("Pi",1:k,sep=""), 
                  paste(rep(paste("beta",0:(q-1),sep=""),k),rep(1:k,each=q),sep=""))
  
  return(out)
  
}

#################################################################################





#################################################################################
# nb.mixEMReg() 
# Regress mixture of Negative Binomials
#################################################################################


nb.mixEMReg <- function(y,x,k=k,theta.t=NULL,eps=1e-6,maxit=1000){
  
  # library(MASS)
  
  # initial useful variables
  n <- length(y)
  X <- cbind(1,x)
  q <- ncol(X)
  
  ## initial values if is.null
  if(is.null(theta.t)){
    
    out.pois <- pois.mixReg.bestfit(y,x,k,maxr=100)
    
    ## initial beta's
    beta.t <- matrix(out.pois[c("beta01","beta11","beta02","beta12")],nrow=q,ncol=k)
    
    ## initial mixing proportions
    pi.t <- rep(1,k)/k
    
    ## initial dispersions
    phi.t <- rep(1,k)
    
    theta.t <- list(pi.t, beta.t, phi.t)
  }
  
  ## initial values
  pi.t <- theta.t[[1]]
  beta.t <- theta.t[[2]]
  phi.t <- theta.t[[3]]
  
  mu.t <- exp(X %*% beta.t)
  
  ## observed loglikelihood 
  y.k <- matrix(NA,nrow=n, ncol=k)
  for (i in 1:n){
    for (j in 1:k) {
      y.k[i,j] <- dnbinom(y[i],mu=mu.t[i,j],size=phi.t[j])
    }
  }
  ll <- sum(log(rowSums(t(t(y.k)*pi.t))))
  
  ## update  
  iter <- 0
  diff <- 1
  
  ## output
  out <- c(iter, ll, pi.t, c(beta.t), phi.t) 
  
  while((iter <= maxit) & (diff > eps)){
    
    print(iter)
    iter <- iter+1
    
    ## glm.nb fit with weights
    z.t <- t(t(y.k)*pi.t) / rowSums(t(t(y.k)*pi.t)) 
    glm.out <- lapply(1:k,function(i) glm.nb(y~x,weights=z.t[,i]))
    theta.t[[1]] <- rbind(theta.t[[1]],apply(z.t,2,mean))
    theta.t[[2]] <- rbind(theta.t[[2]],sapply(1:k,function(i) glm.out[[i]]$coefficients))
    theta.t[[3]] <- rbind(theta.t[[3]],sapply(1:k,function(i) glm.out[[i]]$theta))
    
    ## update parameters
    pi.t <- as.vector(tail(theta.t[[1]],1))
    mu.t <- exp(X %*% tail(theta.t[[2]],2))
    phi.t <- tail(theta.t[[3]],1)
    
    ## calculate loglikelihood
    y.k <- matrix(NA,nrow=n, ncol=k)
    for (i in 1:n){
      for (j in 1:k) {
        y.k[i,j] <- dnbinom(y[i],mu=mu.t[i,j],size=phi.t[j])
      }
    }
    ll.new <- sum(log(rowSums(t(t(y.k)*colMeans(z.t)))))
    
    diff <- abs(ll-ll.new) ; print(diff)
    ll <- ll.new
    
    ## output
    out <- rbind(out,c(iter, ll.new, pi.t, c(tail(theta.t[[2]],2)), phi.t))
    
  }
  
  colnames(out) <- c("iter", "ll", paste("pi",1:k,sep=""), 
                                   paste(rep(paste("beta",0:(q-1),sep=""),k),rep(1:k,each=q),sep=""),
                                   paste("phi",1:k,sep="") )
  
  return(out)
}

#################################################################################





#################################################################################
# nb.mixEM() 
# Fit mixture of Negative Binomials
#################################################################################

nb.mixEM <- function(y,k=k,theta.t=NULL,eps=1e-6,maxit=1000){
  
  # library(MASS)
  
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

pois.mix <- function(x,k){
  
  # library(flexmix)
  
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
# cmp.mixEM() 
# Fit mixture of Mean-Parameterized CMPs
#################################################################################

cmp.mixEM <- function(x,k=k, 
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


###############################################################
# Z(lambda,nu)
# function to approximate the normalizing constant Z
###############################################################

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


###############################################################
# mu.fun(lambda,nu)
# Solve for mu for given lambda and nu
###############################################################

# library(optimr)

mu.fun <- function(lambda,nu,ylim=150){
  y <- 0:ylim
  fun.s <- function(lambda,mu,nu) abs(sum((y-mu)*exp(y*log(lambda)-nu*lgamma(y+1))))
  mu <- hjn(0,fun.s,lower=0,upper=100,lambda=lambda,nu=nu)$par
  return(mu)
}


###############################################################
# mean.fun(lambda,nu)
# solve for the mu i.e. the expectation for given lambda and nu
###############################################################

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


###############################################################
# lambda.fun(mu,nu)
# Solve for the lambda for given mu and nu
###############################################################

#library(optimr)

lambda.fun <- function(mu,nu,ylim=150){
  y <- seq(0,ylim,by=1)
  fun.s <- function(lambda,mu,nu) abs(mu.fun(lambda,nu)-mu)
  lambda <- hjn(5,fun.s,lower=0.1,upper=Inf,mu=mu,nu=nu)$par
  return(lambda)
}


###############################################################
# pmf(mu,nu)
# Evaluate pmf given mu and nu
###############################################################

pmf <- function(y,mu,nu){
  lambda <- lambda.fun(mu,nu)
  return(exp(y*log(lambda)-nu*lgamma(y+1)) / Z(lambda,nu))
}


###############################################################
# The following three functions are used to solve nu 
###############################################################

###############################################################
# var.fun(lambda,nu)
# variance given lambda and nu
###############################################################

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


###############################################################
# mean_logfactorialy.fun(lambda,nu)
# Find E[log(Y!)] given lambda and nu
###############################################################

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


###############################################################
# mean_ylogfactorialy.fun(lambda,nu)
# Find E[Ylog(Y!)] given lambda and nu
###############################################################

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


###############################################################
# nu.fun(x,mu) 
# Solve the MLE of nu by searching  
###############################################################

nu.fun <- function(y,z,mu,nu=seq(0.1,10,0.1)){
  lambda <- c()
  f <- c()
  for (i in 1:length(nu)) {
    lambda[i] <- lambda.fun(mu,nu[i])
    f[i] <- ((mean_ylogfactorialy.fun(lambda[i],nu[i])-mu*mean_logfactorialy.fun(lambda[i],nu[i])) * ((t(z) %*% y)-mu*sum(z)) / var.fun(lambda[i],nu[i]) - t(z) %*% (lgamma(y+1)) + sum(z)* mean_logfactorialy.fun(lambda[i],nu[i]))
  }
  return(nu[which.min(abs(f))])
}


###############################################################





