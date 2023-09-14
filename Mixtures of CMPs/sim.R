###################################################################################
## Simulation to compare three mixture models
## sample generated from 2 components of mean-parametrized CMP
###################################################################################

source("functions.R")
set.seed(1)


## setting
n <- 200   # sample size
B <- 3     # replicates

## parameters
mu1.t <- 1      # mean parameter in component 1
mu2.t <- 15     # mean parameter in component 2
nu1.t <- 1.4    # dispersion parameter in component 1
nu2.t <- 1.5    # dispersion parameter in component 2
pi.t <- 0.3     # mixing proportion in component 1

# MLEs from mixtures
MLE.cmp <- c()
MLE.pois <- c()
MLE.nb <- c()

for(i in 1:B) {
  
  print(paste("B =",i))
  
  z <- rbinom(n,1,pi.t) #posterior  probabilities
  
  # simulated sample
  x <-   z*sample(c(0:200),n,replace=TRUE, prob=pmf(0:200,mu=mu1.t,nu=nu1.t)) +
    (1-z)*sample(c(0:200),n,replace=TRUE, prob=pmf(0:200,mu=mu2.t,nu=nu2.t))
  
  print("Poissons")
  # fit Poisson mixtures
  out.pois <- pois.mix(x,k=2)
  MLE.pois <- rbind(MLE.pois,out.pois)
  
  print("NBs")
  # fit Negative Binomial mixtures
  out.nb <- nb.mixEM(x,k=2,theta.t=NULL,eps=1e-3,maxit=1000)
  MLE.nb <- rbind(MLE.nb,tail(out.nb,1))
  
  print("CMPs")
  # fit CMP mixtures
  out.cmp <- cmpmixEM(x,k=2,
                      mu=NULL,nu=NULL,pi=NULL,
                      eps=1e-3,maxit=1000)
  MLE.cmp <- rbind(MLE.cmp,tail(out.cmp,1))
  
}




