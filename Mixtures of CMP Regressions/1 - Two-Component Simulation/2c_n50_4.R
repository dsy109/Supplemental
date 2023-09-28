###################################################################################
## Monte Carlo Simulation to test the model of CMPs mixture Regression
###################################################################################

setwd("/Users/derekyoung/Documents/CMP Research/JSTP Revision/Simulations/")
source("Dongying_functions.R")
source("Derek_Functions.R")

part <- 4
set.seed(part)   # seed


#################################################################################
## 2-component mixture of mean-parametrized CMP
## Simulate the response variable y and covariate x according to mean-parameterized CMP mixtures
## sample size n=50,100,200
#################################################################################

start_time <- Sys.time()

## setting
n <- 50      # sample size
B <- 250    # replicates

# covariate x
x <- seq(0,10,length=n)  
X <- cbind(1,x) 

k=2 # number of components
q=ncol(X) # number of columns for covariates

# beta
beta01=1.5   ;beta11=.2  # com1 - upper line 
beta02=.5   ;beta12=.4  # com2 - bottom line

beta <- matrix(c(beta01,beta11,beta02,beta12),nrow=q,ncol=k)

# dispersion
nu1 <- 0.5 
nu2 <- 5.5   

# mixing proportion
Pi <- 0.35   

# data to save  
y.B <- c()       # sample point
y.E <- c()       # points generated for non-convergent samples
MLE.cmp <- c()   # estimates

i = 0
while(i < B) {
  i <- i+1
  
  print(paste("B =",i))
  
  # simulate the component response y based on x
  
  z <- rbinom(n,1,Pi)                 #component label 
  y <- rep(0,n)
  ind.1 <- which(z==1)
  ind.2 <- which(z==0)
  mu.true <- exp(X%*%beta)
  y[ind.1] <- rcomp(length(ind.1),mu.true[ind.1,1],nu=nu1)
  y[ind.2] <- rcomp(length(ind.2),mu.true[ind.2,2],nu=nu2)
  
  # regress CMP mixtures
  tmp.test <- TRUE
  tmp.ind <- 0
  while(tmp.test){
    tmp.ind <- tmp.ind + 1
    out.cmp <- try(withTimeout(mix.mpcmp2(y,x,k=2,
                                          start=list(beta=beta,lambda=c(Pi,1-Pi),nu=c(nu1,nu2)),
                                          eps=1e-3,maxit=1000),timeout=20,onTimeout = "silent"),silent=TRUE)
    tmp.test <- ifelse(class(out.cmp)=="try-error",T,F)
    if(tmp.ind == 30){
      tmp.test <- FALSE
      out.cmp <- NULL
    }
  }
  if(is.null(out.cmp)){
    i <- i-1
    y.E <- rbind(y.E, y) #response variable y (bad samples)
  } else{
    y.B <- rbind(y.B, y) #response variable y (good samples)
    MLE.cmp <- rbind(MLE.cmp,
                     c(length(out.cmp$all.ll),tail(out.cmp$all.ll,1),out.cmp$lambda,c(out.cmp$betas),out.cmp$nus,use.names=F))
    
  if(i%in%seq(25,250,by=25)) save.image(file=paste("2c_n",n,"_",part,".RData",sep=""))
  }
}
colnames(MLE.cmp) <- c("iter","ll","Pi1","Pi2","beta01","beta11","beta02","beta12","nu1","nu2")

end_time <- Sys.time()


## save .RData
save.image(file=paste("2c_n",n,"_",part,".RData",sep=""))


## save to txt files

write(paste("Time is", end_time-start_time, attr(end_time-start_time,"units")), 
      file=paste("2c_n",n,"_time",part,".txt",sep=""))

write.table(y.B, file=paste("2c_n",n,"_samples",part,".txt",sep=""), sep="  ", row.names=1:B)
write.table(y.E, file=paste("2c_n",n,"_samples",part,".txt",sep=""), sep="  ")

write.table(MLE.cmp, file=paste("2c_n",n,"_MLEs",part,".txt",sep=""), sep="  ", row.names=1:B)





