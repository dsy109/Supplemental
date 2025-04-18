###############################################################
###############################################################
# Code for RSM coverage study in Section 4
###############################################################
###############################################################

dir1 = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir1);getwd()

coverage_simulation=function(method=c("PS","CH"),P,alpha,X,x0,cm=NULL,adjust=F,betas,sigma){
  
  coeffs = betas
  model_sd = sigma

  # Get the functions to find k-factors
  source("Sim_k_factor.R")
    
  if (method=="PS"){
    TF = adjust
    k_factor = PS_k_factor(P.target=P,alpha.target=alpha,X=X,x0=x0,adjust = TF)
  } else if (method=="CH"){ # Chvostekova's method
    k_factor = C2_k_factor(P.target=P,alpha.target=alpha,X=X,x0=x0,h_sq=NULL)$kfactors
  }
  
  # EYs = EY(x0[,2])
  EYs = as.vector(x0%*%coeffs)
  
  number_of_iteration = 10000
  TI_length = matrix(NA,ncol=dim(x0)[1],nrow=number_of_iteration)
  LBs = matrix(NA,ncol=dim(x0)[1],nrow=number_of_iteration)
  UBs = matrix(NA,ncol=dim(x0)[1],nrow=number_of_iteration)
  count = 0
  
  for (i in 1:number_of_iteration){
    
    # Generate a sample from the true model:
    # Y = -19041.9 + 17930x + e where \e~N(0,130.5^2)
    
    set.seed(i)
    # obs = -19041.9 + 17930*X[,2] + rnorm(dim(X)[1],mean=0,sd=130.5)
    obs = X%*%coeffs + rnorm(dim(X)[1],mean=0,sd=model_sd)
    
    # Based on the generated sample, fit a linear model
    mod = lm(obs~X[,-1])
    
    # Based on the generated sample, find the standard error estimate
    s = summary(mod)$sigma
  
    Y.hat = x0 %*% as.matrix(mod$coefficients)
  
    TI_lower = Y.hat - k_factor*s 
    TI_upper = Y.hat + k_factor*s
    
    TI_length[i,] = TI_upper - TI_lower
    LBs[i,] = TI_lower
    UBs[i,] = TI_upper
    
    L_percentile = pnorm(TI_lower[], EYs[], model_sd)
    U_percentile = pnorm(TI_upper[], EYs[], model_sd)
    
    TIs = as.data.frame(cbind(TI_lower,TI_upper,EYs,L_percentile,U_percentile))
    colnames(TIs) = c("TI_lower","TI_upper","EYs","L_percentile","U_percentile")
    
    TIs$coverage = TIs$U_percentile - TIs$L_percentile 
    TIs$results = TIs$coverage>P
    joint_result = prod(TIs$results)
    
    count = count+joint_result
    #print(i)
  }
  
  coverage_prob = count/number_of_iteration
  
  return(list(k_factors = k_factor,coverage = coverage_prob))
              
} # end of the functio coverage_simulation()

# When m=2 Circumscribed CCD
# Table 3-4 (3-5)
# For Face-centered CCD change the X and x0 matrix.

X =　x0 = matrix(c(1,-1,-1,
             1,-1,+1,
             1,+1,-1,
             1,+1,+1,
             1,-sqrt(2),0,
             1,+sqrt(2),0,
             1,0,-sqrt(2),
             1,0,+sqrt(2),
             1,0,0,
             1,0,0,
             1,0,0,
             1,0,0,
             1,0,0),ncol=3,nrow=13,byrow=T)

# Example

coverage_simulation(method="PS",P=0.75,alpha=0.05,X=X,x0=x0,adjust=T,
                    betas=c(150,6,25),sigma=10)$coverage
# [1] 0.9633
 
# Use squared term

X =　x0 = matrix(c(1,-1,-1,
                  1,-1,+1,
                  1,+1,-1,
                  1,+1,+1,
                  1,-sqrt(2),0,
                  1,+sqrt(2),0,
                  1,0,-sqrt(2),
                  1,0,+sqrt(2),
                  1,0,0,
                  1,0,0,
                  1,0,0,
                  1,0,0,
                  1,0,0),ncol=3,nrow=13,byrow=T)
colnames(X)=c("intercept","X1","X2")

X_sq = x0_sq = cbind(X,
                     X1_sq=X[,"X1"]^2,
                     X2_sq=X[,"X2"]^2)

# Example

coverage_simulation(method="PS",P=0.75,alpha=0.05,X=X_sq,x0=x0_sq,adjust=T,
                    betas=c(150,6,25,-10,-20),sigma=10)$coverage
# 0.9802

coverage_simulation(method="PS",P=0.90,alpha=0.05,X=X_sq,x0=x0_sq,adjust=T,
                    betas=c(150,6,25,-10,-20),sigma=10)$coverage

# Use interaction term
X =　x0 = matrix(c(1,-1,-1,
                  1,-1,+1,
                  1,+1,-1,
                  1,+1,+1,
                  1,-sqrt(2),0,
                  1,+sqrt(2),0,
                  1,0,-sqrt(2),
                  1,0,+sqrt(2),
                  1,0,0,
                  1,0,0,
                  1,0,0,
                  1,0,0,
                  1,0,0),ncol=3,nrow=13,byrow=T)
colnames(X) = c("intercept","X1","X2")
X_intr = x0_intr = cbind(X,X12=X[,"X1"]*X[,"X2"])

# Interaction model: alpha=0.05

coverage_simulation(method="PS",P=0.75,alpha=0.05,X=X_intr,x0=x0_intr,adjust=T,
                    betas=c(150,6,25,-15),sigma=10)$coverage
# 0.9690

# Full model

X =　x0 = matrix(c(1,-1,-1,
                  1,-1,+1,
                  1,+1,-1,
                  1,+1,+1,
                  1,-sqrt(2),0,
                  1,+sqrt(2),0,
                  1,0,-sqrt(2),
                  1,0,+sqrt(2),
                  1,0,0,
                  1,0,0,
                  1,0,0,
                  1,0,0,
                  1,0,0),ncol=3,nrow=13,byrow=T)
colnames(X) = c("intercept","X1","X2")
X_full = x0_full = cbind(X,X1_sq=X[,"X1"]^2,X2_sq=X[,"X2"]^2,X12=X[,"X1"]*X[,"X2"])

# Full model: alpha=0.05

coverage_simulation(method="PS",P=0.75,alpha=0.05,X=X_full,x0=x0_full,adjust=T,
                    betas=c(150,6,25,-10,-20,-15),sigma=10)$coverage
# 0.9844

####################################################
####################################################
####################################################
# CCD with k=3
# Table 3-6, 3-7, 3-8

X = matrix(c(1,-1,-1,-1,
             1,-1,-1,+1,
             1,-1,+1,-1,
             1,+1,-1,-1,
             1,-1,+1,+1,
             1,+1,-1,+1,
             1,+1,+1,-1,
             1,+1,+1,+1,
             1,-2^(3/4),0,0,
             1,+2^(3/4),0,0,
             1,0,-2^(3/4),0,
             1,0,+2^(3/4),0,
             1,0,0,-2^(3/4),
             1,0,0,+2^(3/4),
             1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0),ncol=4,nrow=19,byrow=T)
colnames(X)=c("intercept","X1","X2","X3")
x0 = X  

print(xtable(X,digits=3), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)


# Main effect model
coverage_simulation(method="PS",P=0.90,alpha=0.1,
                    X=X,x0=x0,adjust=T,
                    betas=c(5,3,-0.5,-3),sigma=5)$coverage
# [1] 0.9324

# Main + sqaured term
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2)
coverage_simulation(method="PS",P=0.90,alpha=0.1,
                    X=X_sq,x0=x0_sq,adjust=T,
                    betas=c(5,3,-0.5,-3,1,1,1),sigma=5)$coverage
# [1] 0.9711

# Main + two-way interaction
X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),
                           X13=(X[,"X1"])*(X[,"X3"]),
                           X23=(X[,"X2"])*(X[,"X3"]))
coverage_simulation(method="PS",P=0.90,alpha=0.1,
                    X=X_intr,x0=x0_intr,adjust=T,
                    betas=c(5,3,-0.5,-3,1,1,1),sigma=5)$coverage
# [1] 0.9579

# Full model: Main + squared + two-way interaction
X_full = x0_full = cbind(X,
                         X1_sq=(X[,"X1"])^2,
                         X2_sq=(X[,"X2"])^2,
                         X3_sq=(X[,"X3"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),
                         X13=(X[,"X1"])*(X[,"X3"]),
                         X23=(X[,"X2"])*(X[,"X3"]))
coverage_simulation(method="PS",P=0.90,alpha=0.1,
                    X=X_full,x0=x0_full,adjust=T,
                    betas=c(5,3,-0.5,-3,1,1,1,1,1,1),sigma=5)$coverage
# [1] 0.9827


####################################################################
# Function to output the table for all P, gamma, method combination.
# Tables from 3-4 to 3.11

make_tb = function(row_name,betas,sigma,X,x0){
  
  names = as.data.frame(cbind(c(row_name,row_name,row_name),c(0.75,0.9,0.95)) ) 
  
  alpha.vector = c(0.05,0.10,0.20)
  
  for (i in 1:3){
    print(paste("i=",i))
    
    table = matrix(NA,ncol=2,nrow=3)
    alpha = alpha.vector[i]
    P.vector = c(0.75,0.9,0.95)
    
    for (j in 1:3){
      print(paste("j=",j))
      P = P.vector[j]
      table[j,1] = coverage_simulation(method="PS",P=P,alpha=alpha,
                                       X=X,x0=x0,adjust=T,
                                       betas=betas,sigma=sigma)$coverage
    } # end j
    
    for (k in 1:3){
      print(paste("k=",k))
      P = P.vector[k]
      table[k,2] = coverage_simulation(method="CH",P=P,alpha=alpha,
                          X=X,x0=x0,
                          betas=betas,sigma=sigma)$coverage
    } # end k
    
    if (i == 1){
      result1 = table
    } else if (i == 2){
      result2 = table
    } else if (i==3){
      result3 = table
    }
    
  } # end of for i
  
  result = cbind(names,as.data.frame(cbind(result1,result2,result3)))
  return(result)
}

# CCD(2) Main 
X =　x0 = matrix(c(1,-1,-1,1,-1,+1,1,+1,-1,
                  1,+1,+1,1,-sqrt(2),0,1,+sqrt(2),0,
                  1,0,-sqrt(2),1,0,+sqrt(2),
                  1,0,0,1,0,0,1,0,0,1,0,0,1,0,0),ncol=3,nrow=13,byrow=T)
main = make_tb(row_name="Main",betas=c(150,6,25),sigma=10,X=X,x0=x0)

# CCD(2) Full model
X =　x0 = matrix(c(1,-1,-1, 1,-1,+1, 1,+1,-1, 1,+1,+1,1,-sqrt(2),0,
                  1,+sqrt(2),0,1,0,-sqrt(2),1,0,+sqrt(2),
                  1,0,0,1,0,0,1,0,0,1,0,0,1,0,0),ncol=3,nrow=13,byrow=T)
colnames(X) = c("intercept","X1","X2")
X_full = x0_full = cbind(X,X1_sq=X[,"X1"]^2,X2_sq=X[,"X2"]^2,X12=X[,"X1"]*X[,"X2"])

full = make_tb(row_name="Full",betas=c(150,6,25,-10,-20,-15),sigma=10,X=X_full,x0=x0_full)
print(xtable(full,digits = 4),include.rownames=F)

# CCD(2) Face-centered Main 

X =　x0 = matrix(c(1,-1,-1,1,-1,+1,1,+1,-1,
                  1,+1,+1,1,-1,0,1,+1,0,
                  1,0,-1,1,0,+1,
                  1,0,0,1,0,0),ncol=3,nrow=10,byrow=T)

print(xtable(X,digits=1), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

main = make_tb(row_name="Main",betas=c(150,6,25),sigma=10,X=X,x0=x0)
print(xtable(main,digits = 4),include.rownames=F)

# CCD(2) Face-centered Polynomial

colnames(X) = c("intercept","X1","X2")
X_sq = x0_sq = cbind(X,X1_sq=X[,"X1"]^2,X2_sq=X[,"X2"]^2)
poly = make_tb(row_name="Polynomial",betas=c(150,6,25,-10,-20),sigma=10,X=X_sq,x0=x0_sq)
print(xtable(poly,digits = 4),include.rownames=F)

# CCD(2) Face-centered Interaction 

colnames(X) = c("intercept","X1","X2")
X_intr = x0_intr = cbind(X,X12=X[,"X1"]*X[,"X2"])
interaction = make_tb(row_name="Interaction",betas=c(150,6,25,-15),sigma=10,X=X_intr,x0=x0_intr)
print(xtable(interaction,digits = 4),include.rownames=F)

# CCD(2) Face-centered Full

colnames(X) = c("intercept","X1","X2")
X_full = x0_full = cbind(X,X1_sq=X[,"X1"]^2,X2_sq=X[,"X2"]^2,X12=X[,"X1"]*X[,"X2"])
full = make_tb(row_name="Full",betas=c(150,6,25,-10,-20,-15),sigma=10,X=X_full,x0=x0_full)
print(xtable(full,digits = 4),include.rownames=F)



### CCD(3) ###

X = matrix(c(1,-1,-1,-1,
             1,-1,-1,+1,
             1,-1,+1,-1,
             1,+1,-1,-1,
             1,-1,+1,+1,
             1,+1,-1,+1,
             1,+1,+1,-1,
             1,+1,+1,+1,
             1,-2^(3/4),0,0,
             1,+2^(3/4),0,0,
             1,0,-2^(3/4),0,
             1,0,+2^(3/4),0,
             1,0,0,-2^(3/4),
             1,0,0,+2^(3/4),
             1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0),ncol=4,nrow=19,byrow=T)
colnames(X)=c("intercept","X1","X2","X3")
x0 = X  
xtable(X)

# Main CCD3
mainccd3 = make_tb(row_name="Main",
                   betas=c(5,3,-0.5,-3),
                   sigma=5,X=X,x0=x0)
print(xtable(mainccd3,digits = 4),include.rownames=F)

# Pylinomial CCD(3)
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2)
polyccd3 = make_tb(row_name="Polynomial",
                     betas=c(5,3,-0.5,-3,1,1,1),
                     sigma=5,X=X_sq,x0=x0_sq)
print(xtable(polyccd3,digits = 4),include.rownames=F)

# Interaction CCD(3)
X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X23=(X[,"X2"])*(X[,"X3"]))
intrccd3 = make_tb(row_name="Interaction",
                 betas=c(5,3,-0.5,-3,1,1,1),
                 sigma=5,X=X_intr,x0=x0_intr)

print(xtable(intrccd3,digits = 4),include.rownames=F)

# Full model CCD(3)

X_full = x0_full = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X23=(X[,"X2"])*(X[,"X3"]))
fullccd3 = make_tb(row_name="Full",
                   betas=c(5,3,-0.5,-3,1,1,1,1,1,1),
                   sigma=5,X=X_full,x0=x0_full)
print(xtable(fullccd3,digits = 4),include.rownames=F)

########################################################################
# Main Face-centered CCD(3)

X = matrix(c(1,-1,-1,-1,1,-1,-1,+1,1,-1,+1,-1,1,+1,-1,-1,
             1,-1,+1,+1,1,+1,-1,+1,1,+1,+1,-1,1,+1,+1,+1,
             1,-1,0,0,1,+1,0,0,1,0,-1,0,1,0,+1,0,1,0,0,-1,1,0,0,+1,
             1,0,0,0,1,0,0,0),ncol=4,nrow=16,byrow=T)
colnames(X)=c("intercept","X1","X2","X3")
x0 = X  
print(xtable(X,digits=1), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

mainfcccd3 = make_tb(row_name="Main",
                     betas=c(5,3,-0.5,-3),
                     sigma=5,X=X,x0=x0)
print(xtable(mainfcccd3,digits = 4),include.rownames=F)

# Pylinomial Face-Centered CCD(3)
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2)
polyfcccd3 = make_tb(row_name="Polynomial",
                   betas=c(5,3,-0.5,-3,1,1,1),
                   sigma=5,X=X_sq,x0=x0_sq)
print(xtable(polyfcccd3,digits = 4),include.rownames=F)

# Interaction Face-Centered CCD(3)
X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X23=(X[,"X2"])*(X[,"X3"]))
intrfcccd3 = make_tb(row_name="Interaction",
                   betas=c(5,3,-0.5,-3,1,1,1),
                   sigma=5,X=X_intr,x0=x0_intr)
print(xtable(intrfcccd3,digits = 4),include.rownames=F)

# Full model Face-Centered CCD(3)

X_full = x0_full = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X23=(X[,"X2"])*(X[,"X3"]))
fullfcccd3 = make_tb(row_name="Full",
                   betas=c(5,3,-0.5,-3,1,1,1,1,1,1),
                   sigma=5,X=X_full,x0=x0_full)
print(xtable(fullfcccd3,digits = 4),include.rownames=F)

####################################################
# BBD with k=3

X = matrix(c(1,-1,-1,0,
             1,+1,-1,0,
             1,-1,+1,0,
             1,+1,+1,0,
             1,-1,0,-1,
             1,+1,0,-1,
             1,-1,0,+1,
             1,+1,0,+1,
             1,0,-1,-1,
             1,0,+1,-1,
             1,0,-1,+1,
             1,0,+1,+1,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             1,0,0,0),ncol=4,nrow=17,byrow=T)
colnames(X)=c("intercept","X1","X2","X3")
x0 = X  
print(xtable(X,digits=3), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)


# Main BBD(3)
mainbbd3 = make_tb(row_name="Main",
                   betas=c(5,3,-0.5,-3),
                   sigma=5,X=X,x0=x0)
print(xtable(mainbbd3,digits = 4),include.rownames=F)

# Pylinomial BBD(3)
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2)
polybbd3 = make_tb(row_name="Polynomial",
                   betas=c(5,3,-0.5,-3,1,1,1),
                   sigma=5,X=X_sq,x0=x0_sq)
print(xtable(polybbd3,digits = 4),include.rownames=F)

# Interaction BBD(3)
X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X23=(X[,"X2"])*(X[,"X3"]))
intrbbd3 = make_tb(row_name="Interaction",
                   betas=c(5,3,-0.5,-3,1,1,1),
                   sigma=5,X=X_intr,x0=x0_intr)
print(xtable(intrbbd3,digits = 4),include.rownames=F)

# Full model BBDD(3)

X_full = x0_full = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X23=(X[,"X2"])*(X[,"X3"]))
fullbbdd3 = make_tb(row_name="Full",
                   betas=c(5,3,-0.5,-3,1,1,1,1,1,1),
                   sigma=5,X=X_full,x0=x0_full)
print(xtable(fullbbdd3,digits = 4),include.rownames=F)

### Circumscribed CCD(4) ###

X = matrix(c(1,-1,-1,-1,-1,
             1,-1,-1,-1,+1,
             1,-1,-1,+1,-1,
             1,-1,+1,-1,-1,
             1,+1,-1,-1,-1,
             1,-1,-1,+1,+1,
             1,-1,+1,-1,+1,
             1,+1,-1,-1,+1,
             1,+1,-1,+1,-1,
             1,-1,+1,+1,-1,
             1,+1,+1,-1,-1,
             1,+1,+1,+1,-1,
             1,+1,+1,-1,+1,
             1,+1,-1,+1,+1,
             1,-1,+1,+1,+1,
             1,+1,+1,+1,+1,
             1,-2, 0, 0, 0,
             1,+2, 0, 0, 0,
             1, 0,-2, 0, 0,
             1, 0,+2, 0, 0,
             1, 0, 0,-2, 0,
             1, 0, 0,+2, 0,
             1, 0, 0, 0,-2,
             1, 0, 0, 0,+2,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0),ncol=5,nrow=29,byrow=T)
             
            
colnames(X)=c("intercept","X1","X2","X3","X4")
x0 = X  
print(xtable(X,digits=1), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

# Main CCD(4)

mainccd4 = make_tb(row_name="Main",betas=c(1,2,3,4,5),
                   sigma=5,X=X,x0=x0)
print(xtable(mainccd4,digits = 4),include.rownames=F)

# Polynomial CCD(4)
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,X4_sq=(X[,"X4"])^2)
polyccd4 = make_tb(row_name="Polynomial",betas=c(1,2,3,4,5,-2,3,-4,5),
                   sigma=5,X=X_sq,x0=x0_sq)
print(xtable(polyccd4,digits = 4),include.rownames=F)

# Interaction CCD(4)

X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X14=(X[,"X1"])*(X[,"X4"]),
                         X23=(X[,"X2"])*(X[,"X3"]),X24=(X[,"X2"])*(X[,"X4"]),X34=(X[,"X3"])*(X[,"X4"]))
intrccd4 = make_tb(row_name="Interaction",betas=c(1,2,3,4,5,1,2,-3,4,-5,6),
                   sigma=5,X=X_intr,x0=x0_intr)
print(xtable(intrccd4,digits = 4),include.rownames=F)

# Full CCD(4)

X_full = x0_full = cbind(X,
                         X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,X4_sq=(X[,"X4"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X14=(X[,"X1"])*(X[,"X4"]),
                         X23=(X[,"X2"])*(X[,"X3"]),X24=(X[,"X2"])*(X[,"X4"]),X34=(X[,"X3"])*(X[,"X4"]))
fullccd4 = make_tb(row_name="Full",betas=c(1,2,3,4,5,-2,3,-4,5,1,2,-3,4,-5,6),
                    sigma=5,X=X_full,x0=x0_full)
print(xtable(fullccd4,digits = 4),include.rownames=F)




### Face Centered CCD(4) ###

X = matrix(c(1,-1,-1,-1,-1,
             1,-1,-1,-1,+1,
             1,-1,-1,+1,-1,
             1,-1,+1,-1,-1,
             1,+1,-1,-1,-1,
             1,-1,-1,+1,+1,
             1,-1,+1,-1,+1,
             1,+1,-1,-1,+1,
             1,+1,-1,+1,-1,
             1,-1,+1,+1,-1,
             1,+1,+1,-1,-1,
             1,+1,+1,+1,-1,
             1,+1,+1,-1,+1,
             1,+1,-1,+1,+1,
             1,-1,+1,+1,+1,
             1,+1,+1,+1,+1,
             1,-1, 0, 0, 0,
             1,+1, 0, 0, 0,
             1, 0,-1, 0, 0,
             1, 0,+1, 0, 0,
             1, 0, 0,-1, 0,
             1, 0, 0,+1, 0,
             1, 0, 0, 0,-1,
             1, 0, 0, 0,+1,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0),ncol=5,nrow=26,byrow=T)
colnames(X)=c("intercept","X1","X2","X3","X4")
x0 = X  
print(xtable(X,digits=1), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

# Main Face-Centered CCD(4)
mainfcccd4 = make_tb(row_name="Main",betas=c(1,2,3,4,5),sigma=5,X=X,x0=x0)
print(xtable(mainfcccd4,digits = 4),include.rownames=F)

# Polynomial Face-Centered CCD(4)
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,X4_sq=(X[,"X4"])^2)
polyfcccd4 = make_tb(row_name="Polynomial",betas=c(1,2,3,4,5,-2,3,-4,5),sigma=5,X=X_sq,x0=x0_sq)
print(xtable(polyfcccd4,digits = 4),include.rownames=F)

# Interaction Face-Centered CCD(4)
X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X14=(X[,"X1"])*(X[,"X4"]),X23=(X[,"X2"])*(X[,"X3"]),X24=(X[,"X2"])*(X[,"X4"]),X34=(X[,"X3"])*(X[,"X4"]))
intrfcccd4 = make_tb(row_name="Interaction",betas=c(1,2,3,4,5,1,2,-3,4,-5,6),sigma=5,X=X_intr,x0=x0_intr)
print(xtable(intrfcccd4,digits = 4),include.rownames=F)

# Full Face-Centered CCD(4)
X_full = x0_full = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,X4_sq=(X[,"X4"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X14=(X[,"X1"])*(X[,"X4"]),
                         X23=(X[,"X2"])*(X[,"X3"]),X24=(X[,"X2"])*(X[,"X4"]),X34=(X[,"X3"])*(X[,"X4"]))
fullfcccd4 = make_tb(row_name="Full",betas=c(1,2,3,4,5,-2,3,-4,5,1,2,-3,4,-5,6),sigma=5,X=X_full,x0=x0_full)
print(xtable(fullfcccd4,digits = 4),include.rownames=F)


# BBD(4)

X = matrix(c(1,-1,-1,0,0,
             1,-1,+1,0,0,
             1,+1,-1,0,0,
             1,+1,+1,0,0,
             1,-1,0,-1,0,
             1,-1,0,+1,0,
             1,+1,0,-1,0,
             1,+1,0,+1,0,
             1,-1,0,0,-1,
             1,-1,0,0,+1,
             1,+1,0,0,-1,
             1,+1,0,0,+1,
             1,0,-1,-1,0,
             1,0,-1,+1,0,
             1,0,+1,-1,0,
             1,0,+1,+1,0,
             1,0,-1,0,-1,
             1,0,-1,0,+1,
             1,0,+1,0,-1,
             1,0,+1,0,+1,
             1,0,0,-1,-1,
             1,0,0,-1,+1,
             1,0,0,+1,-1,
             1,0,0,+1,+1,
             1,0,0, 0, 0,
             1,0,0, 0, 0,
             1,0,0, 0, 0,
             1,0,0, 0, 0,
             1,0,0, 0, 0),ncol=5,nrow=29,byrow=T)

colnames(X)=c("intercept","X1","X2","X3","X4")
x0 = X  
print(xtable(X,digits=1), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

# Main BBD(4)
mainbbd4 = make_tb(row_name="Main",betas=c(1,2,3,4,5),sigma=5,X=X,x0=x0)
print(xtable(mainbbd4,digits = 4),include.rownames=F)

# Polynomial BBD(4)
X_sq = x0_sq = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,X4_sq=(X[,"X4"])^2)
polybbd4 = make_tb(row_name="Polynomial",betas=c(1,2,3,4,5,-2,3,-4,5),sigma=5,X=X_sq,x0=x0_sq)
print(xtable(polybbd4,digits = 4),include.rownames=F)

# Interaction BBD(4)
X_intr = x0_intr = cbind(X,X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X14=(X[,"X1"])*(X[,"X4"]),X23=(X[,"X2"])*(X[,"X3"]),X24=(X[,"X2"])*(X[,"X4"]),X34=(X[,"X3"])*(X[,"X4"]))
intrbbd4 = make_tb(row_name="Interaction",betas=c(1,2,3,4,5,1,2,-3,4,-5,6),sigma=5,X=X_intr,x0=x0_intr)
print(xtable(intrbbd4,digits = 4),include.rownames=F)

# Full BBD(4)
X_full = x0_full = cbind(X,X1_sq=(X[,"X1"])^2,X2_sq=(X[,"X2"])^2,X3_sq=(X[,"X3"])^2,X4_sq=(X[,"X4"])^2,
                         X12=(X[,"X1"])*(X[,"X2"]),X13=(X[,"X1"])*(X[,"X3"]),X14=(X[,"X1"])*(X[,"X4"]),
                         X23=(X[,"X2"])*(X[,"X3"]),X24=(X[,"X2"])*(X[,"X4"]),X34=(X[,"X3"])*(X[,"X4"]))
fullbbd4 = make_tb(row_name="Full",betas=c(1,2,3,4,5,-2,3,-4,5,1,2,-3,4,-5,6),sigma=5,X=X_full,x0=x0_full)
print(xtable(fullbbd4,digits = 4),include.rownames=F)

# Figures: Section 4, Figure 2

library(ggplot2)
library(data.table) # to use melt() for data transpose wide to long
library(car)

## For publication

# Section 4 Figure (a)

mydata_ccd2 = read.csv("ccd2_long.csv")
mydata_ccd2$P_num     = mydata_ccd2$P
mydata_ccd2$P     = as.factor(mydata_ccd2$P)
mydata_ccd2$gamma_num = mydata_ccd2$gamma
mydata_ccd2$gamma = as.factor(mydata_ccd2$gamma)

long_ccd2 = melt(mydata_ccd2, measure.vars = c("APS", "CH"))
setnames(long_ccd2, "variable", "method")
setnames(long_ccd2, "value", "coverage")

long80_ccd2 = subset(long_ccd2,gamma_num==0.8)
p_fig_4a = ggplot(long80_ccd2, aes(x=P,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("design:",unique(long80_ccd2$design) ,"nominal coverage: 0.8"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.8) +
  scale_y_continuous(breaks=c(0.8,0.85, 0.9)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4a

ggsave("Fig_4a_CCD2_gamma80.png",
       p_fig_4a, width = 5, height = 5, dpi=300, dev='png')

# Section 4 Figure (b)

long95_ccd2 = subset(long_ccd2,gamma_num==0.95)
p_fig_4b = ggplot(long95_ccd2, aes(x=P,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("design:",unique(long95_ccd2$design) ,"nominal coverage: 0.95"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.95) +
  scale_y_continuous(breaks=c(0.95,0.96,0.97,0.98)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4b

ggsave("Fig_4b_CCD2_gamma95.png",
       p_fig_4b, width = 5, height = 5, dpi=300, dev='png')

# Section 4 Figure (c)

?read.csv
library(utils)
remove.packages(utils)

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()
list.files()

mydata_all = read.csv(file = "all_long.csv",header = T)

library(dplyr)

mydata_ccd_90_90 = mydata_all %>% filter(design == "CCD",P==0.9,gamma==0.9)

mydata_ccd_90_90$P_num     = mydata_ccd_90_90$P
mydata_ccd_90_90$P         = as.factor(mydata_ccd_90_90$P)
mydata_ccd_90_90$gamma_num = mydata_ccd_90_90$gamma
mydata_ccd_90_90$gamma     = as.factor(mydata_ccd_90_90$gamma)

p_fig_4c = ggplot(mydata_ccd_90_90, 
                  aes(x=m,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("design:",unique(mydata_ccd_90_90$design) ,"P: 0.9 nominal coverage: 0.9"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.9) +
  scale_y_continuous(breaks=c(0.9,0.925,0.95,0.975)) +
  scale_x_continuous(breaks=c(2,3,4)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4c

ggsave("Fig_4c_CCD_90_90.png",
       p_fig_4c, width = 5, height = 5, dpi=300, dev='png')

# Section 4 Figure (d)

mydata_all = read.csv("all_long.csv")

library(dplyr)
mydata_m3_90_90 = mydata_all %>% filter(m == 3,P==0.9,gamma==0.9)

mydata_m3_90_90$P_num     = mydata_m3_90_90$P
mydata_m3_90_90$P         = as.factor(mydata_m3_90_90$P)
mydata_m3_90_90$gamma_num = mydata_m3_90_90$gamma
mydata_m3_90_90$gamma     = as.factor(mydata_m3_90_90$gamma)

p_fig_4d = ggplot(mydata_m3_90_90, 
                  aes(x=design,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("m:",unique(mydata_m3_90_90$m) ,"P: 0.9 nominal coverage: 0.9"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.9) +
  scale_y_continuous(breaks=c(0.9,0.92,0.94,0.96,0.98)) +
  #scale_x_continuous(breaks=c(2,3,4)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4d

ggsave("Fig_4d_m3_90_90.png",
       p_fig_4d, width = 5, height = 5, dpi=300, dev='png')

# Beta regresssion

?as.factor

install.packages("betareg")
library(betareg)

mydata = mydata_all
mydata$model = as.factor(mydata$model)
mydata$model =  relevel(mydata$model,ref="Main")           

summary(lm())

mydata <- within(mydata, model <- relevel(model, ref = "Main"))
mydata$model <- relevel(mydata$model, ref = 3)

table(mydata$model)

summary(lm(data=mydata, coverage ~ gamma + P + method + model + m + design))

betamod = betareg(data=mydata, coverage ~ gamma + P + method + model + m + design)
summary(betamod)

beta_reg_table = round(summary(betamod)$coefficients$mean,4)

exp_est = round(exp(summary(betamod)$coefficients$mean[,"Estimate"]),2)

beta_reg_table = cbind(beta_reg_table,exp_est)

write.csv(beta_reg_table,"beta_regression_table.csv")

AIC(mod)
AIC(betamod)

# Model with interaction term

mod2 = lm(data=mydata, 
          coverage ~ gamma + P + method + model + m + design + 
            gamma*model + gamma*method + P*model + P*method + method*model + P*method*model)

summary(mod2)

betamod2 = betareg(data=mydata, 
                   coverage ~ gamma + P + method + model + m + design + 
                     gamma*model + gamma*method + P*model + P*method + method*model +
                     P*method*model)

summary(betamod2)

beta_reg_table2 = round(summary(betamod2)$coefficients$mean,4)
exp_est2 = round(exp(summary(betamod2)$coefficients$mean[,"Estimate"]),2)
beta_reg_table2 = cbind(beta_reg_table2,exp_est2)

write.csv(beta_reg_table2,"beta_regression_table2.csv")

# Drop non-sgnificant interaction

betamod3 = betareg(data=mydata, 
                   coverage ~ gamma + P + method + model + m + design + 
                     gamma*model + gamma*method + P*model + P*method + method*model)

summary(betamod3)

betamod4 = betareg(data=mydata, 
                   coverage ~ gamma + P + method + model + m + design + 
                     gamma*method + P*model + P*method + method*model)

summary(betamod4)

betamod5 = betareg(data=mydata, 
                   coverage ~ gamma + P + method + model + m + design + 
                     gamma*method + P*model + P*method)

summary(betamod5)

beta_reg_table5 = round(summary(betamod5)$coefficients$mean,4)
exp_est5 = round(exp(summary(betamod5)$coefficients$mean[,"Estimate"]),2)
beta_reg_table5 = cbind(beta_reg_table5,exp_est5)

write.csv(beta_reg_table5,"beta_regression_table_nonsig_dropped.csv")
