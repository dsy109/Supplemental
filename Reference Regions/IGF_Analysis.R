###############################################################
###############################################################
### IGF Analysis
###############################################################
###############################################################

###############################################################
### Preamble
###############################################################

library(tolerance)
library(ddalpha)
library(ggplot2)
library(mixtools)
library(ellipse)
library(plyr)
IGF <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/refs/heads/main/Reference%20Regions/IGF_Data_noise.txt",header=T)
attach(IGF)
Sex <- relevel(as.factor(Sex),ref="m")
BMI <- Weight/(Height/100)^2
matt.dat <- data.frame(IGF[,c(7:9,5)],Sex=Sex,BMI=BMI)

Liu <- function(pts,x){
  out <- apply(pts,1,depth.simplicial,data=x,exact=F,k=0.1+0.8999*((ncol(x)==2)&nrow(x)<=150)-0.0999*(ncol(x)>2|nrow(x)>150))
  out
}

out <- lm(cbind(log(IGF_I),IGFBP_2^.25,log(IGFBP_3))~I(Age-45)+I((Age-45)^2)+Sex+I(BMI-25),data=matt.dat)
n <- 427
m <- 4
R <- out$residuals
S <- (1/(n-m-1))*t(R)%*%R
R2 <- cov2cor(S)

new.X <- expand.grid(Age=seq(20,75,by=5),BMI=seq(15,40,by=5),Sex=relevel(as.factor(c("f","m")),ref="m"))
Y.hat <- predict(out,newdata=new.X)
colnames(Y.hat) <- c("log(IGF_I)","IGFBP_2^.25","log(IGFBP_3)")

matt.95.95 <- npmvtol.region(R,alpha=0.05,P=0.95,depth.fn=Liu,type="central")
matt.95.e <- npmvtol.region(R,Beta=0.95,depth.fn=Liu,type="central")
matt.f <- matt.dat[which(matt.dat$Sex=="f"),]
matt.f$Age <- as.factor(sapply(1:nrow(matt.f),function(i) seq(20,75,by=5)[which.min(abs(matt.f$Age[i]-seq(20,75,by=5)))]))

###############################################################
### Quantile-Based Reference Regions
###############################################################

# IGFBP_2 vs. IGF_I (Females)
ggplot(data = matt.f, aes(IGF_I, IGFBP_2)) + geom_rect(aes(xmin=quantile(IGF_I,0.025), xmax=quantile(IGF_I,0.975), ymin=quantile(IGFBP_2,0.025), ymax=quantile(IGFBP_2,0.975)),fill = "gray80",alpha=.1, color="black", linetype=2) +
  geom_point(size = 2, aes(color = Age)) + xlim(0,500) + ylim(0,2000) + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-2 (",mu,"g/L)"))) + theme(text = element_text(size = 20)) +
  ggtitle("Rectangular Reference Regions (All Females)",subtitle="Quantile-Based Regions")

# IGFBP_3 vs. IGF_I (Females)
ggplot(data = matt.f, aes(IGF_I, IGFBP_3)) + geom_rect(aes(xmin=quantile(IGF_I,0.025), xmax=quantile(IGF_I,0.975), ymin=quantile(IGFBP_3,0.025), ymax=quantile(IGFBP_3,0.975)),fill = "gray80",alpha=.1, color="black", linetype=2) +
  geom_point(size = 2, aes(color = Age)) + xlim(0,500) + ylim(0,7000) + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)"))) + theme(text = element_text(size = 20)) +
  ggtitle("Rectangular Reference Regions (All Females)",subtitle="Quantile-Based Regions")

# IGFBP_3 vs. IGFBP_2 (Females)
ggplot(data = matt.f, aes(IGFBP_2, IGFBP_3)) + geom_rect(aes(xmin=quantile(IGFBP_2,0.025), xmax=quantile(IGFBP_2,0.975), ymin=quantile(IGFBP_3,0.025), ymax=quantile(IGFBP_3,0.975)),fill = "gray80",alpha=.1, color="black", linetype=2) +
  geom_point(size = 2, aes(color = Age)) + xlim(0,2000) + ylim(0,7000) + 
  xlab(expression(paste("S-IGFBP-2 (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)"))) + theme(text = element_text(size = 20)) +
  ggtitle("Rectangular Reference Regions (All Females)",subtitle="Quantile-Based Regions")

###############################################################
### (0.95,0.95) Tolerance Regions
###############################################################

# IGFBP_2 vs. IGF_I (Females, by Age, fixed BMI of 20)
matt.df1 <- data.frame(xlow=exp(Y.hat[13:24,1]+matt.95.95[1,1]),
                       xup=exp(Y.hat[13:24,1]+matt.95.95[1,2]),
                       ylow=(Y.hat[13:24,2]+matt.95.95[2,1])^4,
                       yup=(Y.hat[13:24,2]+matt.95.95[2,2])^4,Age=as.factor(seq(20,75,by=5)))

ggplot(matt.df1, aes(xmin = xlow, xmax = xup, ymin = ylow, ymax = yup)) + 
  geom_rect(aes(fill = Age, linetype=Age), colour = "black", alpha=.2) + 
  geom_rect(aes(linetype=Age), colour = "black", fill=NA) + xlim(0,500) + ylim(0,2000) +
  theme(text = element_text(size = 20)) + ggtitle("Rectangular Reference Regions (Females, BMI=20)", subtitle="(0.95,0.95) Tolerance Regions") + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-2 (",mu,"g/L)")))

# IGFBP_3 vs. IGF_I (Females, by Age, fixed BMI of 20)
matt.df3 <- data.frame(xlow=exp(Y.hat[13:24,1]+matt.95.95[1,1]),
                       xup=exp(Y.hat[13:24,1]+matt.95.95[1,2]),
                       ylow=exp(Y.hat[13:24,3]+matt.95.95[3,1]),
                       yup=exp(Y.hat[13:24,3]+matt.95.95[3,2]),Age=as.factor(seq(20,75,by=5)))

ggplot(matt.df3, aes(xmin = xlow, xmax = xup, ymin = ylow, ymax = yup)) +
  geom_rect(aes(fill = Age, linetype=Age), colour = "black", alpha=.2) + 
  geom_rect(aes(linetype=Age), colour = "black", fill=NA) + xlim(0,500) + ylim(0,7000) +
  theme(text = element_text(size = 20)) + ggtitle("Rectangular Reference Regions (Females, BMI=20)", subtitle="(0.95,0.95) Tolerance Regions") + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)")))

# IGFBP_3 vs. IGFBP_2 (Females, by Age, fixed BMI of 20)
matt.df5 <- data.frame(xlow=(Y.hat[13:24,2]+matt.95.95[2,1])^4,
                       xup=(Y.hat[13:24,2]+matt.95.95[2,2])^4,
                       ylow=exp(Y.hat[13:24,3]+matt.95.95[3,1]),
                       yup=exp(Y.hat[13:24,3]+matt.95.95[3,2]),Age=as.factor(seq(20,75,by=5)))

ggplot(matt.df5, aes(xmin = xlow, xmax = xup, ymin = ylow, ymax = yup)) +
  geom_rect(aes(fill = Age, linetype=Age), colour = "black", alpha=.2) + 
  geom_rect(aes(linetype=Age), colour = "black", fill=NA) + xlim(0,2000) + ylim(0,7000) +
  theme(text = element_text(size = 20)) + ggtitle("Rectangular Reference Regions (Females, BMI=20)", subtitle="(0.95,0.95) Tolerance Regions") + 
  xlab(expression(paste("S-IGFBP-2 (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)")))

###############################################################
### 0.95-Expectation Tolerance Regions
###############################################################

# IGFBP_2 vs. IGF_I (Females, by Age, fixed BMI of 20)
matt.df1e <- data.frame(xlow=exp(Y.hat[13:24,1]+matt.95.e[1,1]),
                       xup=exp(Y.hat[13:24,1]+matt.95.e[1,2]),
                       ylow=(Y.hat[13:24,2]+matt.95.e[2,1])^4,
                       yup=(Y.hat[13:24,2]+matt.95.e[2,2])^4,Age=as.factor(seq(20,75,by=5)))

ggplot(matt.df1e, aes(xmin = xlow, xmax = xup, ymin = ylow, ymax = yup)) + 
  geom_rect(aes(fill = Age, linetype=Age), colour = "black", alpha=.2) + 
  geom_rect(aes(linetype=Age), colour = "black", fill=NA) + xlim(0,500) + ylim(0,2000) +
  theme(text = element_text(size = 20)) + ggtitle("Rectangular Reference Regions (Females, BMI=20)", subtitle="0.95-Expectation Tolerance Regions") + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-2 (",mu,"g/L)")))

# IGFBP_3 vs. IGF_I (Females, by Age, fixed BMI of 20)
matt.df3e <- data.frame(xlow=exp(Y.hat[13:24,1]+matt.95.e[1,1]),
                       xup=exp(Y.hat[13:24,1]+matt.95.e[1,2]),
                       ylow=exp(Y.hat[13:24,3]+matt.95.e[3,1]),
                       yup=exp(Y.hat[13:24,3]+matt.95.e[3,2]),Age=as.factor(seq(20,75,by=5)))

ggplot(matt.df3e, aes(xmin = xlow, xmax = xup, ymin = ylow, ymax = yup)) +
  geom_rect(aes(fill = Age, linetype=Age), colour = "black", alpha=.2) + 
  geom_rect(aes(linetype=Age), colour = "black", fill=NA) + xlim(0,500) + ylim(0,7000) +
  theme(text = element_text(size = 20)) + ggtitle("Rectangular Reference Regions (Females, BMI=20)", subtitle="0.95-Expectation Tolerance Regions") + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)")))

# IGFBP_3 vs. IGFBP_2 (Females, by Age, fixed BMI of 20)
matt.df5e <- data.frame(xlow=(Y.hat[13:24,2]+matt.95.e[2,1])^4,
                       xup=(Y.hat[13:24,2]+matt.95.e[2,2])^4,
                       ylow=exp(Y.hat[13:24,3]+matt.95.e[3,1]),
                       yup=exp(Y.hat[13:24,3]+matt.95.e[3,2]),Age=as.factor(seq(20,75,by=5)))

ggplot(matt.df5e, aes(xmin = xlow, xmax = xup, ymin = ylow, ymax = yup)) +
  geom_rect(aes(fill = Age, linetype=Age), colour = "black", alpha=.2) + 
  geom_rect(aes(linetype=Age), colour = "black", fill=NA) + xlim(0,2000) + ylim(0,7000) +
  theme(text = element_text(size = 20)) + ggtitle("Rectangular Reference Regions (Females, BMI=20)", subtitle="0.95-Expectation Tolerance Regions") + 
  xlab(expression(paste("S-IGFBP-2 (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)")))

###############################################################
### 95% Confidence Regions
###############################################################

Age=seq(20,75,by=5)
BMI=seq(15,40,by=5)

s1 <- summary(out)[[1]]$sigma
s2 <- summary(out)[[2]]$sigma
s3 <- summary(out)[[3]]$sigma

# IGFBP_2 vs. IGF_I (Females, by Age, fixed BMI of 20)
tmp1=lapply(1:12, function(i) data.frame(ellipse(R2[c(1,2),c(1,2)],scale=c(s1,s2),centre=c(Y.hat[12+i,1:2]),level=.95),Age=as.factor(Age[i])))
tmp=tmp1[[1]]
for(i in 2:12) tmp=rbind(tmp,tmp1[[i]])

tmp[,1] <- exp(tmp[,1])
tmp[,2] <- tmp[,2]^4

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls <- ddply(tmp, "Age", find_hull)

ggplot(data=tmp,aes(x=x,y=y,color=Age,fill=Age)) + 
  geom_polygon(data = hulls, aes(fill = Age, linetype=Age), colour = "black", alpha = 0.2) + 
  xlim(0,500) + ylim(0,1500) + theme(text = element_text(size = 20)) + ggtitle("Ovoidal Reference Regions (Females, BMI=20)", subtitle="Normal-Based 95% Confidence Region") + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-2 (",mu,"g/L)")))

# IGFBP_3 vs. IGF_I (Females, by Age, fixed BMI of 20)
tmp1=lapply(1:12, function(i) data.frame(ellipse(R2[c(1,3),c(1,3)],scale=c(s1,s3),centre=c(Y.hat[12+i,c(1,3)]),level=.95),Age=as.factor(Age[i])))
tmp=tmp1[[1]]
for(i in 2:12) tmp=rbind(tmp,tmp1[[i]])

tmp[,1] <- exp(tmp[,1])
tmp[,2] <- exp(tmp[,2])

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls <- ddply(tmp, "Age", find_hull)

ggplot(data=tmp,aes(x=x,y=y,color=Age,fill=Age)) + 
  geom_polygon(data = hulls, aes(fill = Age, linetype=Age), colour = "black", alpha = 0.2) + 
  xlim(0,500) + ylim(0,7000) + theme(text = element_text(size = 20)) + ggtitle("Ovoidal Reference Regions (Females, BMI=20)", subtitle="Normal-Based 95% Confidence Region") + 
  xlab(expression(paste("S-IGF-I (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)")))

# IGFBP_3 vs. IGFBP_2 (Females, by Age, fixed BMI of 20)
tmp1=lapply(1:12, function(i) data.frame(ellipse(R2[c(2,3),c(2,3)],scale=c(s2,s3),centre=c(Y.hat[12+i,c(2,3)]),level=.95),Age=as.factor(Age[i])))
tmp=tmp1[[1]]
for(i in 2:12) tmp=rbind(tmp,tmp1[[i]])

tmp[,1] <- (tmp[,1])^4
tmp[,2] <- exp(tmp[,2])

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls <- ddply(tmp, "Age", find_hull)

ggplot(data=tmp,aes(x=x,y=y,color=Age,fill=Age)) + 
  geom_polygon(data = hulls, aes(fill = Age, linetype=Age), colour = "black", alpha = 0.2) + 
  xlim(0,1500) + ylim(0,7000) + theme(text = element_text(size = 20)) + ggtitle("Ovoidal Reference Regions (Females, BMI=20)", subtitle="Normal-Based 95% Confidence Region") + 
  xlab(expression(paste("S-IGFBP-2 (",mu,"g/L)"))) + ylab(expression(paste("S-IGFBP-3 (",mu,"g/L)")))





