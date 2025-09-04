###############################################################
###############################################################
### Thyroid Data Analysis
###############################################################
###############################################################

###############################################################
### Preamble
###############################################################

Thyroid_Data <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/refs/heads/main/Reference%20Regions/Thyroid_Data_noise.txt",header=T)
library(tolerance)
library(EnvStats)
library(mvtnorm)
library(ddalpha)
library(ggplot2)
library(plotly)

BS.fn <- function(x,B=1000,P,alpha,method=c("sim","sim.cent","rect","rect.cent")){
  method <- match.arg(method)
  x.bar <- apply(x,2,mean)
  p <- length(x.bar)
  n <- nrow(x)
  S <- cov(x)
  s.ii <- diag(S)
  R <- cov2cor(S)
  x.bar.bs <- rmvnorm(B,sigma=R/n)
  S.bs <- rWishart(B,n-1,R)/(n-1)
  if(method=="sim"){
    #Alg 1 of Lucagbo and Mathew (2023)
    f1 <- function(c.rho,x.bar.bs,S.bs,P) min(pnorm(x.bar.bs+c.rho*sqrt(S.bs)))-P 
    c.rho <- sapply(1:B,function(i) uniroot(f1,interval=c(0,100),x.bar.bs=x.bar.bs[i,],S.bs=diag(S.bs[,,i]),P=P)$root)
    c1.rho <- quantile(c.rho,1-alpha)
  } else if(method=="sim.cent"){
    #Alg 2 of Lucagbo and Mathew (2023)
    z.q <- qnorm((1+P)/2)
    k.b <- sapply(1:B,function(i) max((abs(x.bar.bs[i,])+z.q)/sqrt(diag(S.bs[,,i]))))
    c1.rho <- quantile(k.b,1-alpha)
  } else if(method=="rect"){
    #Alg 4 of Lucagbo and Mathew (2023)
    f2 <- function(c.rho,x.bar.bs,S.bs,R,P) pmvnorm(lower=x.bar.bs-c.rho*sqrt(diag(S.bs)),upper=x.bar.bs+c.rho*sqrt(diag(S.bs)),corr=R)-P 
    c.rho <- sapply(1:B,function(i) uniroot(f2,interval=c(0,100),x.bar.bs=x.bar.bs[i,],S.bs=S.bs[,,i],R=R,P=P)$root)
    c1.rho <- quantile(c.rho,1-alpha)
  } else if(method=="rect.cent"){
    #New Alg of Mathew and Young
    f3 <- function(c.rho,x.bar,S,P) pmvnorm(lower=x.bar-c.rho*sqrt(diag(S)),upper=x.bar+c.rho*sqrt(diag(S)),mean=x.bar,sigma=S)-P 
    c.rho <- uniroot(f3,interval=c(0,100),x.bar=x.bar,S=S,P=P)$root
    W.b <- sapply(1:B,function(i) max((abs(x.bar.bs[i,])+c.rho)/sqrt(diag(S.bs[,,i]))) )
    c1.rho <- quantile(W.b,1-alpha)
  } 
  c1.rho
}

Liu <- function(pts,x){
  out <- apply(pts,1,depth.simplicial,data=x,exact=F,k=0.1+0.8999*((ncol(x)==2)&nrow(x)<=150)-0.0999*(ncol(x)>2|nrow(x)>150))
  out
}

#Tukey's Equivalent Blocks Code from Liu, Bretz, and Cortina-Borja (2024)
minSampSize.1d.2s <- function(rp, rconf){
  n0 = ceiling( log(1-rconf)/log(rp) ) #this is the min sample size for 1-d 1-s
  
  nk=n0
  while ( pbeta(rp, nk-1, 2, lower.tail=FALSE) < rconf ){
    nk <- nk + 1
  }
  
  return(nk)  #the minimum sample size required
}

numEquiBloc<-function(rp,rconf,nn){
  if(nn<minSampSize.1d.2s(rp,rconf)){
    print(c("The sample size is too small"))
    return("The sample size is too small")
  }else{
    nk=nn
    while ( pbeta(rp, nk, nn+1-nk, lower.tail=FALSE) > rconf ){
      nk <- nk - 1
    }
    k = nk + 1   #The  k  required
    return(k)
  } 
}

DistFreeTolRect.3d <- function(xx,rp,rgamma){
  nn <- dim(xx)[1]
  nk <- numEquiBloc(rp,rgamma,nn)
  nkD <- nn+1-nk #number of equiv blocks to delete
  nkDCycl <- floor(nkD/6) #number of cycles for striping 6 blocks R,T,F,L,B,H
  nkDRem <- nkD - 6*nkDCycl #number of blocks (<6) to strip further
  
  for (i in 1:nkDCycl) { #Strip away nkDCycl cycles of  6 blocks R,T,F,L,B,H
    
    maxRind <- which(xx[,1]==max(xx[,1])) #Indices of the obs having the max x_1 value
    maxRval <- xx[maxRind[1],1] #Record the max x_1 value
    xx <- xx[-maxRind[1],]  #Delete this observation from the sample; [1] is used
    # in case there may be ties
    maxTind <- which(xx[,2]==max(xx[,2])) #Indices of the obs having the max x_2 value
    maxTval <- xx[maxTind[1],2] #Record the max x_2 value
    xx <- xx[-maxTind[1],]  #Delete this observation from the sample
    
    maxFind <- which(xx[,3]==max(xx[,3])) #Indices of the obs having the max x_3 value
    maxFval <- xx[maxFind[1],3] #Record the max x_3 value
    xx <- xx[-maxFind[1],]  #Delete this observation from the sample
    
    minLind <- which(xx[,1]==min(xx[,1])) #Indices of the obs having the min x_1 value
    minLval <- xx[minLind[1],1] #Record the min x_1 value
    xx <- xx[-minLind[1],]  #Delete this observation from the sample
    
    minBind <- which(xx[,2]==min(xx[,2])) #Indices of the obs having the min x_2 value
    minBval <- xx[minBind[1],2] #Record the min x_2 value
    xx <- xx[-minBind[1],]  #Delete this observation from the sample
    
    minHind <- which(xx[,3]==min(xx[,3])) #Indices of the obs having the min x_3 value
    minHval <- xx[minHind[1],3] #Record the min x_3 value
    xx <- xx[-minHind[1],]  #Delete this observation from the sample
  }
  
  if(nkDRem == 0){ #Strp the remaining (<6) block(s)
    "Do nothing"
  }else if(nkDRem == 1){
    maxRind <- which(xx[,1]==max(xx[,1])) #Indices of the obs having the max x_1 value
    maxRval <- xx[maxRind[1],1] #Record the max x_1 value
    xx <- xx[-maxRind[1],]  #Delete this observation from the sample
  }else if(nkDRem == 2){
    maxRind <- which(xx[,1]==max(xx[,1])) #Indices of the obs having the max x_1 value
    maxRval <- xx[maxRind[1],1] #Record the max x_1 value
    xx <- xx[-maxRind[1],]  #Delete this observation from the sample
    
    maxTind <- which(xx[,2]==max(xx[,2])) #Indices of the obs having the max x_2 value
    maxTval <- xx[maxTind[1],2] #Record the max x_2 value
    xx <- xx[-maxTind[1],]  #Delete this observation from the sample
  }else if(nkDRem == 3){
    maxRind <- which(xx[,1]==max(xx[,1])) #Indices of the obs having the max x_1 value
    maxRval <- xx[maxRind[1],1] #Record the max x_1 value
    xx <- xx[-maxRind[1],]  #Delete this observation from the sample
    
    maxTind <- which(xx[,2]==max(xx[,2])) #Indices of the obs having the max x_2 value
    maxTval <- xx[maxTind[1],2] #Record the max x_2 value
    xx <- xx[-maxTind[1],]  #Delete this observation from the sample
    
    maxFind <- which(xx[,3]==max(xx[,3])) #Indices of the obs having the max x_3 value
    maxFval <- xx[maxFind[1],3] #Record the max x_3 value
    xx <- xx[-maxFind[1],]  #Delete this observation from the sample
  }else if(nkDRem == 4){
    maxRind <- which(xx[,1]==max(xx[,1])) #Indices of the obs having the max x_1 value
    maxRval <- xx[maxRind[1],1] #Record the max x_1 value
    xx <- xx[-maxRind[1],]  #Delete this observation from the sample
    
    maxTind <- which(xx[,2]==max(xx[,2])) #Indices of the obs having the max x_2 value
    maxTval <- xx[maxTind[1],2] #Record the max x_2 value
    xx <- xx[-maxTind[1],]  #Delete this observation from the sample
    
    maxFind <- which(xx[,3]==max(xx[,3])) #Indices of the obs having the max x_3 value
    maxFval <- xx[maxFind[1],3] #Record the max x_3 value
    xx <- xx[-maxFind[1],]  #Delete this observation from the sample
    
    minLind <- which(xx[,1]==min(xx[,1])) #Indices of the obs having the min x_1 value
    minLval <- xx[minLind[1],1] #Record the min x_1 value
    xx <- xx[-minLind[1],]  #Delete this observation from the sample
  }else if(nkDRem == 5){
    maxRind <- which(xx[,1]==max(xx[,1])) #Indices of the obs having the max x_1 value
    maxRval <- xx[maxRind[1],1] #Record the max x_1 value
    xx <- xx[-maxRind[1],]  #Delete this observation from the sample
    
    maxTind <- which(xx[,2]==max(xx[,2])) #Indices of the obs having the max x_2 value
    maxTval <- xx[maxTind[1],2] #Record the max x_2 value
    xx <- xx[-maxTind[1],]  #Delete this observation from the sample
    
    maxFind <- which(xx[,3]==max(xx[,3])) #Indices of the obs having the max x_3 value
    maxFval <- xx[maxFind[1],3] #Record the max x_3 value
    xx <- xx[-maxFind[1],]  #Delete this observation from the sample
    
    minLind <- which(xx[,1]==min(xx[,1])) #Indices of the obs having the min x_1 value
    minLval <- xx[minLind[1],1] #Record the min x_1 value
    xx <- xx[-minLind[1],]  #Delete this observation from the sample
    
    minBind <- which(xx[,2]==min(xx[,2])) #Indices of the obs having the min x_2 value
    minBval <- xx[minBind[1],2] #Record the min x_1 value
    xx <- xx[-minBind[1],]  #Delete this observation from the sample
  }    
  return(list( c(maxRval,maxTval,maxFval,minLval,minBval,minHval), c(nn,nk, nkD) ))
}

region.fn <- function(out.tr,x,alpha){
  Sigma <- cov(x) #Sigma <- cov(x)*(n+1)
  
  Mean <- apply(x, 2, mean)
  es <- eigen(Sigma)
  e1 <- es$vec %*% diag(sqrt(es$val))
  theta <- seq(0, pi, len = 100)
  phi <- seq(0, 2 * pi, len = 100)
  theta.phi <- expand.grid(theta, phi)
  r1 <- sqrt(out.tr)
  v1 <- cbind(r1 * sin(theta.phi[, 1]) * cos(theta.phi[, 
                                                       2]), r1 * sin(theta.phi[, 1]) * sin(theta.phi[, 2]), 
              r1 * cos(theta.phi[, 1]))
  pts <- t(Mean - (e1 %*% t(v1)))
  ellipse.x <- pts[, 1]
  ellipse.y <- pts[, 2]
  ellipse.z <- pts[, 3]
  out <- list(ellipse.x=ellipse.x,ellipse.y=ellipse.y,ellipse.z=ellipse.z)
  out
}


###############################################################
### Reference Interval Calculations
###############################################################

thy_comp <- Thyroid_Data[,1:3]
attach(thy_comp)
#thy_comp[,1] <- log(thy_comp[,1]) 
n <- nrow(thy_comp)
all.ints <- NULL

# Individual CIs
ind.CIs <- sapply(1:3,function(i) mean(thy_comp[,i])+c(-1,1)*qt(0.975,n-1)*sd(thy_comp[,i])/sqrt(n))
ind.CIs
all.ints <- c(ind.CIs)

# Individual PIs
ind.PIs <- sapply(1:3,function(i) as.numeric(predIntNorm(thy_comp[,i])$interval$limits))
ind.PIs
all.ints <- rbind(all.ints,c(ind.PIs))

# Individual TIs
ind.TIs <- sapply(1:3,function(i) as.numeric(normtol.int(thy_comp[,i],side=2,method="EXACT",P=0.95)[4:5]))
ind.TIs
all.ints <- rbind(all.ints,c(ind.TIs))

# Bonferroni CIs
bonf.CIs <- sapply(1:3,function(i) mean(thy_comp[,i])+c(-1,1)*qt(1-(.05/6),n-1)*sd(thy_comp[,i])/sqrt(n))
bonf.CIs
all.ints <- rbind(all.ints,c(bonf.CIs))

# Bonferroni PIs
bonf.PIs <- sapply(1:3,function(i) as.numeric(predIntNorm(thy_comp[,i],conf.level=1-(.05/6))$interval$limits))
bonf.PIs
all.ints <- rbind(all.ints,c(bonf.PIs))

# Bonferroni TIs
bonf.TIs <- sapply(1:3,function(i) as.numeric(normtol.int(thy_comp[,i],side=2,method="EXACT",P=0.95,alpha=(.05/6))[4:5]))
bonf.TIs
all.ints <- rbind(all.ints,c(bonf.TIs))

# Individual Nonparametric beta-expectation TIs
ind.npbetaTI <- sapply(1:3,function(i) as.numeric(npbetol.int(thy_comp[,i],side=2,Beta=0.95)[2:3]))
ind.npbetaTI
all.ints <- rbind(all.ints,c(ind.npbetaTI))

# Individual Nonparametric beta-content TIs
ind.npTI <- sapply(1:3,function(i) as.numeric(nptol.int(thy_comp[,i],side=2,P=0.95)[3:4]))
ind.npTI
all.ints <- rbind(all.ints,c(ind.npTI))

# Bonferroni Nonparametric beta-expectation TIs
bonf.npbetaTI <- sapply(1:3,function(i) as.numeric(npbetol.int(thy_comp[,i],side=2,Beta=1-(0.05/3))[2:3]))
bonf.npbetaTI
all.ints <- rbind(all.ints,c(bonf.npbetaTI))

# Bonferroni Nonparametric beta-content TIs
bonf.npTI <- sapply(1:3,function(i) as.numeric(nptol.int(thy_comp[,i],side=2,P=0.95,alpha=(.05/3))[3:4]))
bonf.npTI
all.ints <- rbind(all.ints,c(bonf.npTI))

# Simultaneous Normal TIs (Alg. 1 of Lucagbo and Mathew, 2023)
set.seed(1)
x.bar <- apply(thy_comp,2,mean)
S <- cov(thy_comp)
c1.rho <- BS.fn(thy_comp,B=10000,P=0.95,alpha=0.05,method="sim")
simnorm.TI <- rbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
simnorm.TI
all.ints <- rbind(all.ints,c(simnorm.TI))

# Simultaneous Central TIs for Normal Data (Alg. 2 of Lucagbo and Mathew, 2023)
set.seed(1)
x.bar <- apply(thy_comp,2,mean)
S <- cov(thy_comp)
c1.rho <- BS.fn(thy_comp,B=10000,P=0.95,alpha=0.05,method="sim.cent")
simcentnorm.TI <- rbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
simcentnorm.TI
all.ints <- rbind(all.ints,c(simcentnorm.TI))

# Rectangular Normal TIs (Alg. 4 of Lucagbo and Mathew, 2023)
set.seed(1)
x.bar <- apply(thy_comp,2,mean)
S <- cov(thy_comp)
c1.rho <- BS.fn(thy_comp,B=10000,P=0.95,alpha=0.05,method="rect")
rectnorm.TI <- rbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
rectnorm.TI
all.ints <- rbind(all.ints,c(rectnorm.TI))

# Rectangular Central Tolerance Region
set.seed(1)
x.bar <- apply(thy_comp,2,mean)
S <- cov(thy_comp)
c1.rho <- BS.fn(thy_comp,B=10000,P=0.95,alpha=0.05,method="rect.cent")
rectcentnorm.TI <- rbind(x.bar-c1.rho*sqrt(diag(S)),x.bar+c1.rho*sqrt(diag(S)))
rectcentnorm.TI
all.ints <- rbind(all.ints,c(rectcentnorm.TI))

# NP MV Beta-Content Tolerance Region (Depth-Based Approach Using Simplicial Depth)
npmvTI <- t(npmvtol.region(x = as.matrix(thy_comp), alpha = 0.05, P = 0.95, depth.fn = Liu))
npmvTI
all.ints <- rbind(all.ints,c(npmvTI))

# NP MV Beta-Expectation Tolerance Region (Depth-Based Approach Using Simplicial Depth)
npmvbetaTI <- t(npmvtol.region(x = as.matrix(thy_comp), Beta = 0.95, depth.fn = Liu))
npmvbetaTI
all.ints <- rbind(all.ints,c(npmvbetaTI))

# Tolerance Region (Tukey's Equivalent Blocks)
tukeyequiv <- DistFreeTolRect.3d(as.matrix(thy_comp),rp=0.95,rgamma=0.95)
tukeyequiv
all.ints <- rbind(all.ints,tukeyequiv[[1]][c(4,1,5,2,6,3)])
all.ints <- data.frame(all.ints)
colnames(all.ints) <- c("TSH_lower","TSH_upper","FT4_lower","FT4_upper","FT3_lower","FT3_upper")
all.ints <- cbind(all.ints,"Method"=toupper(letters[1:17]))
all.ints <- cbind(all.ints,"TSH_mid"=mean(thy_comp[,1]),"FT4_mid"=mean(thy_comp[,2]),"FT3_mid"=mean(thy_comp[,3]))

# MVN 95% Confidence Region, 95% Prediction Region, (0.95,0.95) Tolerance Region
alpha <- 0.05
p <- ncol(thy_comp)
n <- nrow(thy_comp)
out.cr <- p/(n-p)*qf(1-alpha,p,n-p)
conf.region <- region.fn(out.cr, x=thy_comp, alpha=alpha)
out.pr <- p*(n+1)/(n-p)*qf(1-alpha,p,n-p)
pred.region <- region.fn(out.pr, x=thy_comp, alpha=alpha)
set.seed(1)
out.tr <- mvtol.region(x = as.matrix(thy_comp), alpha = alpha, P = 0.95, B = 5000,
                       method = "KM")[1,1]
tol.region <- region.fn(out.tr, x=thy_comp, alpha=alpha)

region.df <- data.frame(TSH=c(conf.region$ellipse.x,pred.region$ellipse.x,tol.region$ellipse.x),
                        FT4=c(conf.region$ellipse.y,pred.region$ellipse.y,tol.region$ellipse.y),
                        FT3=c(conf.region$ellipse.z,pred.region$ellipse.z,tol.region$ellipse.z),
                        Region=factor(rep(c("95% Confidence Region","95% Prediction Region","(0.95, 0.95) Tolerance Region"),each=10000),
                                      levels=c("95% Confidence Region","95% Prediction Region","(0.95, 0.95) Tolerance Region")))

# Plots
meth.labels <- c("Individual 95% CIs",
                 "Individual 95% PIs",
                 "Individual (0.95, 0.95) TIs",
                 "Bonferroni 95% CIs",
                 "Bonferroni 95% PIs",
                 "Bonferroni (0.95, 0.95) TIs",
                 "Individual NP 0.95-Expectation TIs",
                 "Individual NP (0.95, 0.95) TIs",
                 "Bonferroni NP 0.95-Expectation TIs",
                 "Bonferroni NP (0.95, 0.95) TIs",
                 "Simultaneous Normal (0.95, 0.95) TIs",
                 "Simultaneous Central Normal (0.95, 0.95) TIs",
                 "Rectangular Normal (0.95, 0.95) TR",
                 "Rectangular Central Normal (0.95, 0.95) TR",
                 "NP Depth-Based (0.95, 0.95) TR",
                 "NP Depth-Based 0.95-Expectation TR",
                 "(0.95, 0.95) TR - Tukey's Equivalent Blocks")

ggplot(all.ints, aes(x=Method, y=TSH_mid, ymin=TSH_lower, ymax=TSH_upper))+
  geom_pointrange(size=.01,lwd=1.5,col=1:17)+
  geom_hline(yintercept = all.ints$TSH_mid[1], linetype=2)+
  coord_flip()+ggtitle("Different Statistical Intervals for TSH") +
  xlab('') + ylab("TSH (mIU/L)") + theme_bw() + annotate('rect', ymin=.4, ymax=4, xmin=0, xmax=18, alpha=.2, fill='gray50') +
  theme(text = element_text(size = 25), axis.text.y=element_blank()) + 
  geom_text(label=meth.labels,nudge_x=0.5,nudge_y=apply(all.ints[,1:2],1,mean)-all.ints$TSH_mid[1],col=1:17,size=6) + ylim(-5,12)


ggplot(all.ints, aes(x=Method, y=FT4_mid, ymin=FT4_lower, ymax=FT4_upper))+
  geom_pointrange(size=.01,lwd=1.5,col=1:17)+
  geom_hline(yintercept = all.ints$FT4_mid[1], linetype=2)+
  coord_flip()+ggtitle("Different Statistical Intervals for FT4") +
  xlab('') + ylab("FT4 (pmol/L)") + theme_bw() + annotate('rect', ymin=10, ymax=23, xmin=0, xmax=18, alpha=.2, fill='gray50') +
  theme(text = element_text(size = 25), axis.text.y=element_blank()) + 
  geom_text(label=meth.labels,nudge_x=0.5,nudge_y=apply(all.ints[,3:4],1,mean)-all.ints$FT4_mid[1],col=1:17,size=6) 


ggplot(all.ints, aes(x=Method, y=FT3_mid, ymin=FT3_lower, ymax=FT3_upper))+
  geom_pointrange(size=.01,lwd=1.5,col=1:17)+
  geom_hline(yintercept = all.ints$FT3_mid[1], linetype=2)+
  coord_flip()+ggtitle("Different Statistical Intervals for FT3") +
  xlab('') + ylab("FT3 (pmol/L)") + theme_bw() + annotate('rect', ymin=3.1, ymax=6.8, xmin=0, xmax=18, alpha=.2, fill='gray50') +
  theme(text = element_text(size = 25), axis.text.y=element_blank()) +
  geom_text(label=meth.labels,nudge_x=0.5,nudge_y=apply(all.ints[,5:6],1,mean)-all.ints$FT3_mid[1],col=1:17,size=6) 


plot_ly() %>% add_trace(x=~TSH,y=~FT4,z=~FT3,data=region.df,type="mesh3d",color=~Region,alphahull = 0, opacity = 0.25,
                        showlegend=TRUE) %>% add_markers(x = thy_comp[, 1], y = thy_comp[, 2], z = thy_comp[, 3], 
                                                         type = "scatter", marker=list(size=1),showlegend=F, colors='black') %>% layout(title = "Multivariate Normal Regions",scene=list(xaxis=list(title="TSH")))

out.all <- all.ints[,1:7]
rownames(out.all) <- NULL
out.all$Method <- meth.labels
out.all








