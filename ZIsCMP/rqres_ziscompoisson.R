rqres <- function(y, F, eps = 1e-6)
{
  n <- length(y)
  FL <- F(y - eps)
  FU <- F(y)
  u <- runif(n, min = FL, max = FU)
  qres <- qnorm(u)
  return(qres)
}

rqres.ziscmp <- function(y, lambda, nu, p, NuOfVar)
{
	n <- length(y)
	if (length(lambda) == 1) lambda <- rep(lambda, n)
	if (length(p) == 1) p <- rep(p, n)
	if (length(nu) == 1) lambda <- rep(nu, n)
	
	F <- function(y) {
		ret <- numeric(n)
		for (i in 1:n) {
			ret[i] <- p.zi.scompoisson(y[i], lambda=lambda[i], nu=nu[i], p=p[i], NuOfVar=NuOfVar)
		}
		return(ret)
	}
	rqres(y, F)
}

rqres.ziscmp.reg <- function(y, X=NULL, G=NULL, W=NULL, beta0, gamma0, xi0, NuOfVar)
{
	n <- length(y)
	if(is.null(X)){
	  X <- cbind(rep(1,n))
	  if(length(beta0)>1) beta0 <- beta0[1]
	}
	if(is.null(G)){
	  G <- cbind(rep(1,n))
	  if(length(gamma0)>1) gamma0 <- gamma0[1]
	}
	if(is.null(W)){
	  W <- cbind(rep(1,n))
	  if(length(xi0)>1) xi0 <- xi0[1]
	}
	lambda <- exp(X %*% beta0)
	nu <- exp(G %*% gamma0)
	p <- plogis(W %*% xi0)

	F <- function(y) {
		ret <- numeric(n)
		for (i in 1:n) {
			ret[i] <- p.zi.scompoisson(y[i], lambda=lambda[i], nu=nu[i], p=p[i], NuOfVar=NuOfVar)
		}
		return(ret)
	}

	rqres(y, F)
}

