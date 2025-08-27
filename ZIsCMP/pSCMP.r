pSCMP <- function(NuOfVar,lambda,nu,count)
# Code written by Li Zhu and Kimberly Sellers (2014)
# This code computes P(X <= count) for a sCOM-Poisson(lambda, Nu, NuOfVar) distribution

{
  if(NuOfVar == 1)
  {
    print("NuOfVar must be > 1. If interested in computing a cumulative probability where NuOfVar = 1, run dcom in compoisson package.")
    return("TRY AGAIN!")
  }


  i <- 0
  sum1 <- 0
  while(i <= count)
  {
    sum1 <- sum1 + dSCMP(NuOfVar,lambda,nu,i)
    i <- i + 1
  }
  return(sum1)
}