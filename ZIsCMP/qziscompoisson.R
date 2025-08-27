q.zi.scompoisson <- function(pr,NuOfVar,lambda,nu,p)
  {
    if(NuOfVar == 1)
    {
      return("NuOfVar must be > 1 to use this function. For NuOfVar=1, use the compoisson package.")
    }
    fn <- function(pr=pr,NuOfVar=NuOfVar,lambda=lambda,nu=nu,p=p){
      i <- 0
      sum1 <- 0
      if(pr == 0)
      {
        return(0)
      }
      if(pr > 1){return("pr must be between 0 and 1.")}  
      while(sum1 < pr)
      {
        sum1 <- p.zi.scompoisson(q=i,NuOfVar=NuOfVar,lambda=lambda,nu=nu,p=p)
        i <- i + 20
      }
      i <- i - 20
      sum2 <- 0
      j <- pmax(i-20,0)
      while(sum2 < pr)
      {
        sum2 <- p.zi.scompoisson(q=j,NuOfVar=NuOfVar,lambda=lambda,nu=nu,p=p)
        j <- j + 1
      }
      return((j-1))
    }
    fn.v <- Vectorize(fn,vectorize.args="pr")
    out <- fn.v(pr=pr,NuOfVar=NuOfVar,lambda=lambda,nu=nu,p=p)
    out
  }