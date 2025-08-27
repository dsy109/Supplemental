rSCMP <- function(NuOfVar,lambda,nu,size)
# Code written by Li Zhu and Kimberly Sellers (2014)
# This code produces a random number generator of length "size" via the sCOM-Poisson(NuOfVar,lambda, nu) distribution

{
  if(NuOfVar == 1)
  {
    return("NuOfVar must be > 1 to use this function. For NuOfVar=1, use rcom in the compoisson package.")
  }

  i <- 1
  keep <- numeric(0)
  while(i <= size)
  {
    w <- 0
    sum1 <- 0
    pr <- runif(1)
    while(sum1 <= pr)
    {
      sum1 <- sum1 + dSCMP(NuOfVar,lambda,nu,w)
      w <- w + 1
    }
    keep <- c(keep, (w - 1))
    i <- i + 1
  }
  return(keep)
}

