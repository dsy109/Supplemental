qSCMP <- function(NuOfVar,lambda,nu, pr)
# This code determines the quantile q such that "P(X <= q) = pr" for the sCOM-Poisson(NuOfVar,lambda, nu) distribution.
# Code written by Li Zhu and Kimberly Sellers (2014)

{
  if(NuOfVar == 1)
  {
    return("NuOfVar must be > 1 to use this function. For NuOfVar=1, use the compoisson package.")
  }


  i <- 0
  sum1 <- 0
  if(pr == 0)
  {
    return(0)
  }
  if(pr > 1){return("pr must be between 0 and 1.")}  

  while(sum1 < pr)
  {
    sum1 <- pSCMP(NuOfVar,lambda,nu,i)
    i <- i + 1
  }
  return((i - 1))
}