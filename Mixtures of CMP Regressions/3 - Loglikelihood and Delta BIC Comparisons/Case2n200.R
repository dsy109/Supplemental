###################################################################################
## Monte Carlo Simulation to test the model of CMPs mixture Regression
###################################################################################

source("/project/dyo227_uksr/JSTP/Dongying_functions.R")
source("/project/dyo227_uksr/JSTP/Derek_Functions.R")
setwd("/scratch/dyo227/")

set.seed(1)   # seed
case <- 2


#################################################################################
## 2-component mixture of mean-parametrized CMP
## Simulate the response variable y and covariate x according to mean-parameterized CMP mixtures
## sample size n=50,100,200
#################################################################################

start_time <- Sys.time()

## setting
n <- 50      # sample size
B <- 1000    # replicates

# covariate x
x <- seq(0,10,length=n)  
X <- cbind(1,x) 

k=2 # number of components
q=ncol(X) # number of columns for covariates

# beta
beta01=0.5   ;beta11=-0.3  # com1 - upper line 
beta02=0.6   ;beta12=0.4  # com2 - bottom line

beta <- matrix(c(beta01,beta11,beta02,beta12),nrow=q,ncol=k)

# dispersion
nu1 <- 0.8 
nu2 <- 1.2   

# mixing proportion
Pi <- 0.5   

# data to save  
MLE.cmp <- c()   # estimates
MLE.nb <- c()
MLE.poi<- c()

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
  } else{
    MLE.cmp <- rbind(MLE.cmp,
                     c(length(out.cmp$all.ll),tail(out.cmp$all.ll,1),out.cmp$lambda,c(out.cmp$betas),out.cmp$nus,use.names=F))
    out.poi <- poisregmixEM(x=x,y=y, beta=beta,lambda=c(Pi,1-Pi), epsilon = 1e-3, k=2)
    MLE.poi <- rbind(MLE.poi,c(out.poi$lambda,c(out.poi$beta),out.poi$loglik))
    out.nb <- tail(nb.mixEMReg(x=x,y=y,theta.t=list(c(Pi,1-Pi),beta,c(1,1)), eps=1e-3, k=2),1)
    MLE.nb <- rbind(MLE.nb,tail(out.nb,1))
    
  if(i%in%seq(50,1000,by=50)) save.image(file=paste("Case_",case,"_n_",n,".RData",sep=""))
  }
}

end_time <- Sys.time()


cmp.ll <- MLE.cmp[,2]
poi.ll <- MLE.poi[,7]
nb.ll <- MLE.nb[,2]
  
ll.comp <- cbind(cmp.ll,poi.ll,nb.ll)
BIC.comp <- cbind(-2*cmp.ll+7*log(n),-2*poi.ll+5*log(n),-2*nb.ll+7*log(n))


## save .RData
save.image(file=paste("Case_",case,"_n_",n,".RData",sep=""))




