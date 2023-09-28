###################################################################################
## Monte Carlo Simulation to test the model of CMPs mixture Regression
###################################################################################

source("/project/dyo227_uksr/JSTP/Dongying_functions.R")
source("/project/dyo227_uksr/JSTP/Derek_Functions.R")
setwd("/scratch/dyo227/")

part <- 4
set.seed(part)   # seed


#################################################################################
## 3-component mixture of mean-parametrized CMP
## Simulate the response variable y and covariate x according to mean-parameterized CMP mixtures
## sample size n=50,100,200
#################################################################################

start_time <- Sys.time()

## setting
n <- 200      # sample size
B <- 250    # replicates

# covariate x
x <- seq(0,10,length=n)  
X <- cbind(1,x) 

k=3 # number of components
q=ncol(X) # number of columns for covariates

# beta
beta01=4   ;beta11=-.2  # com1 - upper line 
beta02=.5   ;beta12=.2  # com2 - bottom line
beta03=1   ;beta13=.4  # com3 - bottom line

beta <- matrix(c(beta01,beta11,beta02,beta12,beta03,beta13),nrow=q,ncol=k)

# dispersion
nu1 <- 1 
nu2 <- .25   
nu3 <- 5

# mixing proportion
Pi1 <- 0.4
Pi2 <- 0.3
Pi3 <- 0.3

# data to save  
y.B <- c()       # sample point
y.E <- c()       # points generated for non-convergent samples
MLE.cmp <- c()   # estimates

i = 0
while(i < B) {
  i <- i+1
  
  print(paste("B =",i))
  
  # simulate the component response y based on x
  
  z <- apply(rmultinom(n,1,c(Pi1,Pi2,Pi3)),2,which.max)                #component label 
  y <- rep(0,n)
  ind.1 <- which(z==1)
  ind.2 <- which(z==2)
  ind.3 <- which(z==3)
  mu.true <- exp(X%*%beta)
  y[ind.1] <- rcomp(length(ind.1),mu.true[ind.1,1],nu=nu1)
  y[ind.2] <- rcomp(length(ind.2),mu.true[ind.2,2],nu=nu2)
  y[ind.3] <- rcomp(length(ind.3),mu.true[ind.3,3],nu=nu3)

  # regress CMP mixtures
  tmp.test <- TRUE
  tmp.ind <- 0
  while(tmp.test){
    tmp.ind <- tmp.ind + 1
    out.cmp <- try(withTimeout(mix.mpcmp2(y,x,k=3,
                                          start=list(betas=beta,lambda=c(Pi1,Pi2,Pi3),nus=c(nu1,nu2,nu3)),
                                          eps=1e-5,maxit=1000),timeout=200,onTimeout = "silent"),silent=TRUE)
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
    
  if(i%in%seq(25,250,by=25)) save.image(file=paste("3c_n",n,"_",part,".RData",sep=""))
  }
}
colnames(MLE.cmp) <- c("iter","ll","Pi1","Pi2","Pi3","beta01","beta11","beta02","beta12","beta03","beta13","nu1","nu2","nu3")

end_time <- Sys.time()


## save .RData
save.image(file=paste("3c_n",n,"_",part,".RData",sep=""))


## save to txt files

write(paste("Time is", end_time-start_time, attr(end_time-start_time,"units")), 
      file=paste("3c_n",n,"_time",part,".txt",sep=""))

write.table(y.B, file=paste("3c_n",n,"_samples",part,".txt",sep=""), sep="  ", row.names=1:B)
write.table(y.E, file=paste("3c_n",n,"_badsamples",part,".txt",sep=""), sep="  ")

write.table(MLE.cmp, file=paste("3c_n",n,"_MLEs",part,".txt",sep=""), sep="  ", row.names=1:B)





