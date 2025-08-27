dSCMP <- function(NuOfVar,Lambda,Nu,count)
# This code computes P(X=count) for a sCOM-Poisson(lambda, Nu, NuOfVar) distribution.
# Code written by Li Zhu and Kimberly Sellers (2014)

{
  lambda <- rep.int(Lambda,NuOfVar)
  nu <- rep.int(Nu,NuOfVar)
  
  if(length(lambda) == 1)
  {
    print("NuOfVar must be > 1. If interested in computing probability where NuOfVar = 1, run dcom in compoisson package.")
    return("TRY AGAIN!")
  }
  
  m <- 1
  while(m <= length(lambda))
  {                                  
    if((lambda[m] <= 0)|(nu[m] <= 0))
    {
      print("lambda and nu have to be greater than 0.")
      return("TRY AGAIN!")
    }
    m <- m + 1
  }
  
  a <- numeric(0)
  keep <- numeric(0)
  total <- 0
    
  if (length(lambda) == 2)
  {
    a[1] <- 0
    a[2] <- count
    while(a[1] <= a[2])
    {
      total <- ((lambda[1])^(a[1])*(lambda[2])^(a[2]-a[1])/(((factorial(a[1]))^(nu[1])*(factorial(a[2]-a[1]))^(nu[2]))))
      keep <- c(keep, total)
      a[1] <- a[1] + 1
    }
    result <- 1/(com.compute.z(lambda[1],nu[1],exp(-30))*(com.compute.z(lambda[2],nu[2],exp(-30))))* sum(keep)
    return(result)
  }
  
  recur <- function(number,iterator)
  {
    a <- numeric(0)
    t <- number
    a[length(lambda)] <- count 
    storage <- vector(mode = "list", length = t)
    storage[[1]] <- numeric(0)
    sum1 <- 0
    sum2 <- 0
    
    if (number == 2)
    {
      a[1] <- 0
      a[2] <- iterator
      storage[[1]] <- numeric(0)
      while(a[1] <= a[2])
      {
        sum1 <- ((lambda[1])^(a[1])*(lambda[2])^(a[2]-a[1])/(((factorial(a[1]))^(nu[1])*(factorial(a[2]-a[1]))^(nu[2]))))
        storage[[1]] <- c(storage[[1]], sum1)
        a[1] <- a[1] + 1
      }
      result <- 1/(com.compute.z(lambda[1],nu[1],exp(-10))*(com.compute.z(lambda[2],nu[2],exp(-10))))* sum(storage[[1]])
      return (sum(storage[[1]]))
    }
    else
    {
      a[t-1] <- 0
      a[t] <- iterator 
      while(a[t-1] <= a[t])
      {
        sum2 <- (recur(number - 1,a[t-1]) * (lambda[t])^(a[t]-a[t-1])/(factorial(a[t]-a[t-1]))^(nu[t]))
        storage[[t]] <- c(storage[[t]], sum2)
        a[t-1] <- a[t-1] + 1
      }
      if (t < length(lambda))
      {
        return(sum(storage[[t]]))
      }
    } 
    
    w <- 1
    result2 <- 1
    while(w <= length(lambda))
    {
      result2 <- result2 /(com.compute.z(lambda[w],nu[w],exp(-10)))
      w <- w + 1
    }
    result2 <- result2 * (sum(storage[[t]]))  
    #cat("lambdas:",lambda,"\n")
    #cat("nus:    ",nu,"\n")
    #cat("Prob of",count,"counts: ", result2,"\n")
    return(result2)  
  }
  

   
  recur(length(lambda),count)
}
