r.zi.scompoisson <- function (n, NuOfVar, lambda, nu, p) 
{
    x <- integer(n)
    z <- rbinom(n, size = 1, prob = p)
    x[z == 1] <- 0
    x[z == 0] <- rSCMP(NuOfVar=NuOfVar, lambda=lambda, nu=nu, size = sum(z == 0))
    return(x)
}