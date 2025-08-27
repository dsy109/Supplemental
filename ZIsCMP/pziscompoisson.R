p.zi.scompoisson <- function (q,NuOfVar,lambda,nu,p) {
  if (length(length(lambda)) > 1) 
    stop("Currently, lambda must be a single number")
  if (length(length(nu)) > 1) 
    stop("Currently, nu must be a single number")
  fn <- Vectorize(pSCMP,vectorize.args="count")
  fx <- suppressWarnings(p * (q >= 0) + (1 - p) * fn(count=q,NuOfVar=NuOfVar,lambda=lambda,nu=nu))
  fx[which(is.nan(fx))] <- 1
  return(fx)
}