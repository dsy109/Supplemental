#################################################################################
#################################################################################
## Bootstrap to calculate SE for the dog mortality data with m=3 CMPs
#################################################################################
## NON-PARAMETRIC
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
estimates <- readRDS(file="bestfit.rds") 

## setting
n <- length(x) ;n    # sample size
B <- 1000            # bootstrap replicates

x.B <- c()
MLE.B <- c()
for(i in 1:B) {
  
  print(paste("B =",i))
  
  # non-parametric sampling
  x.n <- sample(x, n, replace=TRUE)
  
  out <- cmpmixEM(x.n,k=3,
                  mu=c(as.numeric(estimates[4:6])),nu=NULL,pi=NULL,
                  eps=1e-6,maxit=1000)
  
  # store the sample and estimates
  x.B <- rbind(x.B, x.n)
  MLE.B <- rbind(MLE.B, tail(out,1))
}

end_time <- Sys.time()


## save to txt files

write(paste("Time is", end_time-start_time, attr(end_time-start_time,"units")), 
      file="Nboot.time.txt")

write.table(x.B, file="Nboot.samples.txt", sep="  ", row.names=1:B)

write.table(MLE.B, file="Nboot.MLEs.txt", sep="  ", row.names=1:B)

#################################################################################


