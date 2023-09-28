suppressPackageStartupMessages(library(R.utils))
library(mixtools)
library(glmmTMB)
library(mpcmp)
library(nnet)

inv.logit.mix <- function(w,alpha){
  lp <- matrix(w%*%alpha,nrow=nrow(w))
  out <- exp(lp)
  out <- out/(1+apply(out,1,sum))
  #  out <- cbind(out,1-apply(out,1,sum)) #kth component is baseline
  out <- cbind(1-apply(out,1,sum),out) #1st component is baseline
  out
}

mix.mpcmp1 <- function(y,x=NULL,w=NULL,z=NULL,k=2,start=list(betas=NULL,mus=NULL,alpha=NULL,lambda=NULL,omegas=NULL,nus=NULL),eps=1e-5,maxit=10000){
  diff <- Inf
  old.ll <- -Inf
  n <- length(y)
  X <- cbind(1,x)
  W <- cbind(1,w)
  Z <- cbind(1,z)
  p <- ncol(X)
  q <- ncol(W)
  r <- ncol(Z)
  if(is.null(x)){
    if(!all(is.null(w),is.null(z),is.null(start$betas),is.null(start$alpha),is.null(start$omegas))) warning(paste("Only a mixture of univariate MCMP1s is being fit, so supplied covariates are ignored!", "\n"),immediate.=TRUE,call.=FALSE)
    w <- z <- start$betas <- start$alpha <- start$omegas <- NULL
    betas <- start$betas
    alpha <- start$alpha
    omegas <- start$omegas
    out.cmp <- normalmixEM(pmax(0,y+runif(n,-.1,.1)),k=k,lambda=start$lambda,mu=start$mu,eps=1e-5)
    z.post <- out.cmp$posterior
    lambda <- out.cmp$lambda
    mus <- out.cmp$mu
    if(is.null(start$nu)) start$nu=rep(exp(1),k)
    all.ll <- NULL
    iter <- 0
    while(abs(diff)>eps && iter < maxit){
      iter <- iter+1
      if(iter==1){
        out.cmp <- lapply(1:k,function(i) glmmTMB(y~1,disp=~1,weights=z.post[,i],data=data.frame(y),family=compois,
                                                  start=list(beta=log(mus[i]),betad=log(start$nu[i]))))
      } else{
        out.cmp <- lapply(1:k,function(i) glmmTMB(y~1,disp=~1,weights=z.post[,i],data=data.frame(y),family=compois,
                                                  start=list(beta=fixef(out.cmp[[i]])$cond,betad=fixef(out.cmp[[i]])$disp)))
      }
      nus <- 1/sapply(1:k,function(i) exp(fixef(out.cmp[[i]])$disp))
      mus <- sapply(1:k,function(i) exp(fixef(out.cmp[[i]])$cond))
      comp <- sapply(1:k,function(i) lambda[i]*dcomp(y, mus[i],nu=nus[i]))
      new.ll <- sum(log(apply(comp,1,sum)))
      diff <- new.ll-old.ll
      old.ll <- new.ll
      all.ll <- c(all.ll,old.ll)
      z.post <- comp/apply(comp,1,sum)
      lambda <- apply(z.post,2,mean)
      print(diff)
    }
  } else{
    if(!is.null(start$alpha)){
      pp <- exp(W%*%start$alpha)
      pp <- pp/(1+apply(pp,1,sum))
      pp <- cbind(1-apply(pp,1,sum),pp)
      start$lambda <- apply(pp,2,mean)
    }
    out.cmp <- poisregmixEM(y=y,x=x,k=k,beta=start$betas,lambda=start$lambda)
    z.post <- out.cmp$posterior
    if(!is.null(start$omegas)&is.null(start$nus)){
      start$nus <- apply(sapply(1:k,function(i) exp(Z%*%start$omegas)),2,mean)
    } 
    if(!is.null(z)){
      if(is.null(start$omegas)){
        out.nb <- suppressWarnings(sapply(1:k,function(i) fixef(glmmTMB(y~x,dispformula = ~z,weights=z.post[,i],family=compois, start=list(beta=out.cmp$beta[,i])))$disp))
      } else out.nb <- suppressWarnings(sapply(1:k,function(i) fixef(glmmTMB(y~x,dispformula = ~z,weights=z.post[,i],family=compois, start=list(beta=out.cmp$beta[,i],betad=start$omegas[,i])))$disp))
    } else{
      if(is.null(start$nus)){
        out.nb <- suppressWarnings(sapply(1:k,function(i) fixef(glmmTMB(y~x,weights=z.post[,i],family=compois, start=list(beta=out.cmp$beta[,i])))$disp))
      } else out.nb <- suppressWarnings(sapply(1:k,function(i) fixef(glmmTMB(y~x,dispformula = ~1,weights=z.post[,i],family=compois, start=list(beta=out.cmp$beta[,i],betad=log(start$nus[i]))))$disp))
    }
    out.nb <- matrix(out.nb,ncol=k)
    #  out.nb <- sapply(1:k,function(i) glmmTMB(y~x,dispformula = ~Z,weights=z.post[,i],family=nbinom2)
    all.ll <- NULL
    alpha <- start$alpha
    omegas <- start$omegas
    iter <- 0
    while(abs(diff)>eps && iter < maxit){
      iter <- iter+1
      if(is.null(w)){
        lambda <- apply(z.post,2,mean)
      } else{
        z.cat <- multinom(z.post~w)
        alpha <- coef(z.cat)
        lambda <- fitted(z.cat)
      }
      if(iter==1){
        if(is.null(z)){
          out.cmp <- lapply(1:k,function(i) glmmTMB(y~x,disp=~1,weights=z.post[,i],data=data.frame(y,x),family=compois,
                                                    start=list(beta=out.cmp$beta[,i],betad=out.nb[,i])))
        } else{
          out.cmp <- lapply(1:k,function(i) glmmTMB(y~x,disp=~z,weights=z.post[,i],data=data.frame(y,x,z),family=compois,
                                                    start=list(beta=out.cmp$beta[,i],betad=out.nb[,i])))
        }
      } else{
        if(is.null(z)){
          out.cmp <- lapply(1:k,function(i) glmmTMB(y~x,disp=~1,weights=z.post[,i],data=data.frame(y,x),family=compois,
                                                    start=list(beta=fixef(out.cmp[[i]])$cond,betad=fixef(out.cmp[[i]])$disp)))
        } else{
          out.cmp <- lapply(1:k,function(i) glmmTMB(y~x,disp=~z,weights=z.post[,i],data=data.frame(y,x,z),family=compois,
                                                    start=list(beta=fixef(out.cmp[[i]])$cond,betad=fixef(out.cmp[[i]])$disp)))
        }
      }
      betas <- sapply(1:k,function(i) fixef(out.cmp[[i]])$cond)
      if(is.null(z)){
        nus <- 1/sapply(1:k,function(i) exp(fixef(out.cmp[[i]])$disp))
      } else{
        omegas <- sapply(1:k,function(i) fixef(out.cmp[[i]])$disp)
        nus <- 1/sapply(1:k,function(i) exp(Z%*%cbind(omegas[,i])))
      }
      xbeta <- X %*% betas
      #    if(is.null(w)){
      #      comp <- sapply(1:k,function(i) lambda[i]*dcomp(y, exp(xbeta[,i]),nu=nus[i]))
      #    } else comp <- sapply(1:k,function(i) lambda[,i]*dcomp(y, exp(xbeta[,i]),nu=nus[i]))
      if(is.null(w)){
        comp <- sapply(1:k,function(i) lambda[i]*dcomp(y, exp(xbeta[,i]),nu=ifelse(is.null(z),nus[i],nus[,i])))
      } else comp <- sapply(1:k,function(i) lambda[,i]*dcomp(y, exp(xbeta[,i]),nu=ifelse(is.null(z),nus[i],nus[,i])))
      new.ll <- sum(log(apply(comp,1,sum)))
      diff <- new.ll-old.ll
      old.ll <- new.ll
      all.ll <- c(all.ll,old.ll)
      z.post <- comp/apply(comp,1,sum)
      print(diff)
    }
    mus <- xbeta
  }
  if(!is.null(betas)){
    rownames(betas) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
    colnames(betas) <- c(paste("comp", ".", 1:k, sep = ""))
  }
  if(is.matrix(mus)){
    colnames(mus) <- c(paste("comp", ".", 1:k, sep = ""))
  } else names(mus) <- c(paste("comp", ".", 1:k, sep = ""))
  if(!is.null(alpha)){
    alpha <- t(alpha)
    rownames(alpha) <- c(paste("alpha", ".", 0:(q - 1), sep = ""))
    colnames(alpha) <- c(paste("comp", ".", 2:k, sep = ""))
  }
  if(is.matrix(lambda)){
    colnames(lambda) <- c(paste("comp", ".", 1:k, sep = ""))
  } else names(lambda) <- c(paste("comp", ".", 1:k, sep = ""))
  if(!is.null(omegas)){
    rownames(omegas) <- c(paste("omega", ".", 0:(r - 1), sep = ""))
    colnames(omegas) <- c(paste("comp", ".", 1:k, sep = ""))
  }
  if(is.matrix(nus)){
    colnames(nus) <- c(paste("comp", ".", 1:k, sep = ""))
  } else names(nus) <- c(paste("comp", ".", 1:k, sep = ""))
#  names(nus) <- c(paste("comp", ".", 1:k, sep = ""))
  colnames(z.post) <- c(paste("comp", ".", 1:k, sep = ""))
  
  vc <- lapply(1:k,function(i) vcov(out.cmp[[i]],full=TRUE))
#  se.nus <- sqrt(sapply(1:k,function(i) (-exp(-fixef(out.cmp[[i]])$disp))^2*tail(diag(vc[[i]]),1)))
#  names(se.nus) <- c(paste("comp", ".", 1:k, sep = ""))
  if(is.null(betas)){
    se.betas <- NULL
    se.mus <- sqrt(sapply(1:k,function(i) (exp(fixef(out.cmp[[i]])$cond))^2*head(diag(vc[[i]]),1)))
    names(se.mus) <- c(paste("comp", ".", 1:k, sep = ""))
  } else{
    se.mus <- NULL
    se.betas <- sapply(1:k,function(i) sqrt(diag(vc[[i]]))[1:p])
    rownames(se.betas) <- c(paste("se.beta", ".", 0:(p - 1), sep = ""))
    colnames(se.betas) <- c(paste("comp", ".", 1:k, sep = "")) 
  }
  if(is.null(alpha)){
    se.alpha <- NULL
    se.lambda <- sqrt(lambda*(1-lambda)/n)
    names(se.lambda) <- c(paste("comp", ".", 1:k, sep = ""))
  } else{
    se.lambda <- NULL
    se.alpha <- t(summary(z.cat)$standard.errors)
    rownames(se.alpha) <- c(paste("se.alpha", ".", 0:(q - 1), sep = ""))
    colnames(se.alpha) <- c(paste("comp", ".", 2:k, sep = ""))
  }
  if(is.null(omegas)){
    se.omegas <- NULL
    se.nus <- sqrt(sapply(1:k,function(i) (-exp(-fixef(out.cmp[[i]])$disp))^2*tail(diag(vc[[i]]),1)))
    names(se.nus) <- c(paste("comp", ".", 1:k, sep = ""))
    if(any(is.nan(se.nus))) se.nus[which(is.nan(se.nus))] <- 0
  } else{
    se.nus <- NULL
    se.omegas <- sapply(1:k,function(i) sqrt(diag(summary(out.cmp[[i]])$vcov$disp)))
    rownames(se.omegas) <- c(paste("se.alpha", ".", 0:(r - 1), sep = ""))
    colnames(se.omegas) <- c(paste("comp", ".", 1:k, sep = ""))
  }  
  out <- list(alpha=alpha,lambda=lambda,betas=betas,mus=mus,omegas=omegas,nus=nus,posteriors=z.post,all.ll=all.ll,
              se.alpha=se.alpha,se.lambda=se.lambda,se.betas=se.betas,se.mus=se.mus,se.omegas=se.omegas,se.nus=se.nus)
  out
}


mix.mpcmp2 <- function(y,x=NULL,k=2,start=list(betas=NULL,lambda=NULL,nus=NULL),eps=1e-5,maxit=10000){
  diff <- Inf
  old.ll <- -Inf
  n <- length(y)
  X <- cbind(1,x)
  p <- ncol(X)
  out <- poisregmixEM(y=y,x=x,k=k,beta=start$betas,lambda=start$lambda,epsilon=eps)
  z.post <- out$posterior
  betas <- start$betas
  lambda <- start$lambda
  nus <- start$nus
  if(is.null(betas)) betas <- out$beta
  if(is.null(lambda)) lambda <- out$lambda
  if(is.null(nus)) nus <- 1/exp(suppressWarnings(sapply(1:k,function(i) fixef(glmmTMB(y~x,weights=z.post[,i],family=compois, start=list(beta=betas[,i])))$disp)))
  if(is.null(nus)) nus <- rep(exp(1),k)
  iter <- 0
  all.ll <- NULL
  out.cmp <- lapply(1:k,function(i) glmmTMB(y~x,disp=~1,weights=z.post[,i],data=data.frame(y,x),family=compois,
                                            start=list(beta=betas[,i],betad=log(nus[i])),
                                            control=glmmTMBControl(optCtrl=list(iter.max=100, eval.max=150),conv_check = "skip")))
#  print(out.cmp)
  while(abs(diff)>eps && iter < maxit){
    iter <- iter+1
    out.cmp <- lapply(1:k,function(i) glmmTMB(y~x,disp=~1,weights=z.post[,i],data=data.frame(y,x),family=compois,
                                              start=list(beta=fixef(out.cmp[[i]])$cond,betad=fixef(out.cmp[[i]])$disp),
                                              control=glmmTMBControl(optCtrl=list(iter.max=100, eval.max=150),conv_check = "skip")))
    nus <- 1/sapply(1:k,function(i) exp(fixef(out.cmp[[i]])$disp))
    betas <- sapply(1:k,function(i) fixef(out.cmp[[i]])$cond)
    xbeta <- X %*% betas
    comp <- sapply(1:k,function(i) lambda[i]*dcomp(y, exp(xbeta[,i]),nu=nus[i]))
    new.ll <- sum(log(apply(comp,1,sum)))
    diff <- new.ll-old.ll
    old.ll <- new.ll
    all.ll <- c(all.ll,old.ll)
    z.post <- comp/apply(comp,1,sum)
    lambda <- apply(z.post,2,mean)
    print(diff)    
    
  }
  rownames(betas) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
  colnames(betas) <- c(paste("comp", ".", 1:k, sep = ""))
  names(lambda) <- c(paste("comp", ".", 1:k, sep = ""))
  names(nus) <- c(paste("comp", ".", 1:k, sep = ""))
  out <- list(lambda=lambda,betas=betas,nus=nus,posteriors=z.post,all.ll=all.ll)
  out
}  


model.sel <- function(out){  
  n <- nrow(out$posteriors)
  k <- ncol(out$posteriors)
  df <- ifelse(is.null(out$alpha),k-1,length(out$alpha))+length(out$betas)+length(out$nus)
  finalloglik <- tail(out$all.ll,1)
  z <- out$posteriors
  IC <- list()
  IC$AIC <- 2 * finalloglik - df * 2
  IC$BIC <- 2 * finalloglik - df * log(n)
  IC$AIC3 <- 2 * finalloglik - df * 3
  IC$CAIC <- 2 * finalloglik - df * (1 + log(n))
  IC$AWE <- 2 * finalloglik - 2 * df * (3/2 + log(n))
  if (n - df - 1 > 0) {
    IC$AICc <- IC$AIC - (2 * df * (df + 1))/(n - 
                                               df - 1)
    IC$AICu <- IC$AICc - n * log(n/(n - df - 1))
  }
  else {
    IC$AICc <- IC$AICu <- -Inf
  }
  z.const <- (z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z
  hard.z <- (matrix(rep(apply(z, 1, max), k), n, k, 
                    byrow = F) == z) * 1
  ECM <- sum(hard.z * log(z.const))
  IC$ICL <- IC$BIC + ECM
  out <- unlist(IC)
  out
}

rmix.mpcmp1 <- function(n,alpha.true=NULL,lambda.true=NULL,
                       beta.true=NULL,mu.true=NULL,
                       omega.true=NULL,nu.true=NULL){
  if(!is.null(mu.true)){
    k <- length(mu.true)
    if(is.null(lambda.true)|is.null(nu.true)) stop(paste("Must specify both of lambda.true and nu.true!", "\n"))
    Z <- t(sapply(1:n, function(i) rmultinom(1,size=1,prob=lambda.true)))
    y <- vector("numeric",n)
    ind <- sapply(1:n,function(i) which(Z[i,]==1))
    for(j in 1:k){
      ind.j <- which(ind==j)
      y[ind.j] <- rcomp(length(ind.j),mu.true[j],nu=nu.true[j])
    }
    x.mat <- NULL
  } else{
    k <- ncol(beta.true)
    p <- nrow(beta.true)
    if(p<2|p>3) stop(paste("Must have only 2 or 3 covariates!", "\n"))
    if(is.null(alpha.true)&is.null(lambda.true)) stop(paste("Must specify one of alpha.true or lambda.true!", "\n"))
    if(is.null(omega.true)&is.null(nu.true)) stop(paste("Must specify one of omega.true or nu.true!", "\n"))
    if(p==2){
      x.mat <- seq(0,10,length=n)
    } else{
      x1 <- seq(0,10,length=n/10)
      x2 <- seq(0,1,length=10)
      x.mat <- as.matrix(expand.grid(x1,x2))
    }
    mu <- sapply(1:k,function(i) exp(cbind(1,x.mat)%*%beta.true[,i]))
    if(is.null(alpha.true)){
      probs <- matrix(rep(lambda.true,n),nrow = n,byrow=T)
    } else{
      w <- cbind(1,x.mat)
      probs <- inv.logit.mix(w,alpha.true)
    }
    if(is.null(omega.true)){
      disp <- matrix(rep(nu.true,n),nrow = n,byrow=T)
    } else{
      z <- cbind(1,x.mat)
      disp <- exp(z%*%omega.true)
    }
    Z <- t(sapply(1:n, function(i) rmultinom(1,size=1,prob=probs[i,])))
    y <- vector("numeric",n)
    ind <- sapply(1:n,function(i) which(Z[i,]==1))
    for(j in 1:k){
      ind.j <- which(ind==j)
      y[ind.j] <- rcomp(length(ind.j),mu[ind.j,j],nu=disp[ind.j,j])
    }
  }
  out <- list(y=y,x=x.mat,z=ind)
  out
}
  
  
  