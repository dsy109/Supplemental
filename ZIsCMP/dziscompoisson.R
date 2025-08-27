d.zi.scompoisson <- function (x, NuOfVar, Lambda, Nu, p, log = FALSE) 
{
  fn <- Vectorize(dSCMP,vectorize.args="count")
  fx <- p * (x == 0) + (1 - p) * fn(count=x,NuOfVar=NuOfVar,Lambda=Lambda,Nu=Nu)
  if (log) 
    return(log(fx))
  else return(fx)
}
