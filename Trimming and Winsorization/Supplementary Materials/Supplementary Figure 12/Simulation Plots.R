library(ggplot2)
setwd("/Users/kcheng/Desktop/Supplementary Figure 12/Simulation Plots")
######################################
load("Content.Normal.Reg.90 .RData")
content.vec.reg.90 <- content.vec.reg
load("Content.Normal.Tol.90 .RData")
content.vec.tol.90 <- content.vec.tol
load("Content.Normal.Reg.95 .RData")
content.vec.reg.95 <- content.vec.reg
load("Content.Normal.Tol.95 .RData")
content.vec.tol.95 <- content.vec.tol
#######################################
upper.P <- 0.03
lower.P <- 0.07

alpha <- 0.1
B <- 5000
sample.size <- seq(from=50 , to=2000 , by=20)
mean <- 0
sd <- 1
#######################################
cov.reg.90 <- cov.reg.95 <- cov.tol.90 <- cov.tol.95 <- rep(NA,length(sample.size))

for (i in 1:length(sample.size)){
  cov.reg.90[i] <- sum(content.vec.reg.90[,i] >= (1-(upper.P+lower.P)) , na.rm = TRUE)/(B-sum(is.na(content.vec.reg.90[,i])))
  cov.reg.95[i] <- sum(content.vec.reg.95[,i] >= (1-(upper.P+lower.P)) , na.rm = TRUE)/(B-sum(is.na(content.vec.reg.95[,i])))
  cov.tol.90[i] <- sum(content.vec.tol.90[,i] >= (1-(upper.P+lower.P)) , na.rm = TRUE)/(B-sum(is.na(content.vec.tol.90[,i])))
  cov.tol.95[i] <- sum(content.vec.tol.95[,i] >= (1-(upper.P+lower.P)) , na.rm = TRUE)/(B-sum(is.na(content.vec.tol.95[,i])))
}

cov.reg.90
cov.reg.95
cov.tol.90
cov.tol.95
########################################
###### Plot ######
decimal.fun <- function(x) sprintf("%.2f", x)
x.var <- sample.size
y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2
##################
df.90 <- data.frame(x.var , cov.reg.90 , cov.tol.90)
names(df.90) <- c("Size","Regular","Tolerance")

ggplot(data=df.90)+
  geom_point(aes(x=x.var , y=df.90$Regular ,
                 col=rep("Tukey",length(x.var)),
                 shape=rep("Tukey",length(x.var))),cex=7)+
  geom_line(aes(x=x.var , y=df.90$Regular ,
                col=rep("Tukey",length(x.var))),lwd=3)+
  geom_point(aes(x=x.var , y=df.90$Tolerance ,
                 col=rep("Tol",length(x.var)),
                 shape=rep("Tol",length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df.90$Tolerance ,
                col=rep("Tol",length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(sample.size),max(sample.size)), 
                     labels = seq(from=0,to=2000,by=100) , 
                     breaks = seq(from=0,to=2000,by=100)) +
  scale_y_continuous(limits=c(0,1) , labels = decimal.fun(y.lim) ,
                     breaks = y.lim) +
  geom_hline(yintercept = (1-alpha) , col="brown" , lwd=3, lty = 2) +
  theme(legend.position=c(0.5,0.075),legend.direction="horizontal",
        legend.text=element_text(size=25),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Tol'="darkblue" , 'Tukey'="red"), 
                      labels=c('Tol' = "Tolerance       " ,'Tukey'= "Tukey")) +
  scale_shape_manual(name="Method",
                     values=c('Tol'=17 , 'Tukey'=16), 
                     labels=c('Tol' = "Tolerance       " ,'Tukey'= "Tukey")) +
  theme(axis.text.x = element_text(face="bold", size=25 , angle=90) ,
        axis.title.x = element_text(face="bold", size=25) ,
        axis.text.y = element_text(face="bold", size=25) ,
        axis.title.y = element_text(face="bold", size=25))

################################################
alpha=0.05

decimal.fun <- function(x) sprintf("%.2f", x)
x.var <- sample.size
y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2
##################
df.95 <- data.frame(x.var , cov.reg.95 , cov.tol.95)
names(df.95) <- c("Size","Regular","Tolerance")

ggplot(data=df.95)+
  geom_point(aes(x=x.var , y=df.95$Regular ,
                 col=rep("Tukey",length(x.var)),
                 shape=rep("Tukey",length(x.var))),cex=7)+
  geom_line(aes(x=x.var , y=df.95$Regular ,
                col=rep("Tukey",length(x.var))),lwd=3)+
  geom_point(aes(x=x.var , y=df.95$Tolerance ,
                 col=rep("Tol",length(x.var)),
                 shape=rep("Tol",length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df.95$Tolerance ,
                col=rep("Tol",length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(sample.size),max(sample.size)), 
                     labels = seq(from=0,to=2000,by=100) , 
                     breaks = seq(from=0,to=2000,by=100)) +
  scale_y_continuous(limits=c(0,1) , labels = decimal.fun(y.lim) ,
                     breaks = y.lim) +
  geom_hline(yintercept = (1-alpha) , col="brown" , lwd=3, lty = 2) +
  theme(legend.position=c(0.5,0.075),legend.direction="horizontal",
        legend.text=element_text(size=25),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Tol'="darkblue" , 'Tukey'="red"), 
                      labels=c('Tol' = "Tolerance       " ,'Tukey'= "Tukey")) +
  scale_shape_manual(name="Method",
                     values=c('Tol'=17 , 'Tukey'=16), 
                     labels=c('Tol' = "Tolerance       " ,'Tukey'= "Tukey")) +
  theme(axis.text.x = element_text(face="bold", size=25 , angle=90) ,
        axis.title.x = element_text(face="bold", size=25) ,
        axis.text.y = element_text(face="bold", size=25) ,
        axis.title.y = element_text(face="bold", size=25))

