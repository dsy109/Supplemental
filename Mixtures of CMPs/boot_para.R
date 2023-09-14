#################################################################################
#################################################################################
## Bootstrap to calculate SE for the dog mortality data with m=3 CMPs
#################################################################################
## PARAMETRIC
#################################################################################

#################################################################################
## Data
## dog mortality (Tom Lewis 2018) 
## x is the data, however, a request of the raw counts needs to be made to
## Dr. Lewis

#hist(x,breaks=seq(-1,27,1))


#################################################################################
# Non-Parametric bootstrapping
# Sampling with replacement
#################################################################################

start_time <- Sys.time()

## estimates as the starting means
estimates <- readRDS(file="bestfit.rds") ; 
estimates <- as.numeric(estimates)

## setting
n <- length(x) ;n    # sample size
B <- 1000            # bootstrap replicates

x.B <- c()
MLE.B <- c()
for(i in 1:B) {
  
  print(paste("B =",i))
  
  # parametric sampling
  z <- rmultinom(n, size=1, prob=estimates[1:3]);z  
  x.n <- z * rbind(sample(c(0:200),n,replace=TRUE, prob=pmf(0:200,mu=estimates[4],nu=estimates[7])),
                   sample(c(0:200),n,replace=TRUE, prob=pmf(0:200,mu=estimates[5],nu=estimates[8])),
                   sample(c(0:200),n,replace=TRUE, prob=pmf(0:200,mu=estimates[6],nu=estimates[9]))) 
  
  x.n <- apply(x.n, 2, sum)
  
  out <- cmpmixEM(x.n,k=3,
                  mu=c(estimates[4:6]),nu=NULL,pi=NULL,
                  eps=1e-6,maxit=1000)
  
  # store the sample and estimates
  x.B <- rbind(x.B, x.n)
  MLE.B <- rbind(MLE.B, tail(out,1))
}

end_time <- Sys.time()


## save to txt files

write(paste("Time is", end_time-start_time, attr(end_time-start_time,"units")), 
      file="Pboot.time.txt")

write.table(x.B, file="Pboot.samples.txt", sep="  ", row.names=1:B)

write.table(MLE.B, file="Pboot.MLEs.txt", sep="  ", row.names=1:B)

#################################################################################

