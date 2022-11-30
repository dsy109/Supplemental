library(ggplot2)
setwd("/Users/kedaicheng/Desktop/Winsorization/Codes/Bivariate/Independent/Exp9095/Results")
#########################################
M <- 1000
P <- 0.90
alpha <- 0.05
##########################
###### Data Loading ######
##########################
load("Parametric.Independent1 .RData")
output.1 <- output
n.vec.1 <- seq(from=300,to=400,by=20)

load("Parametric.Independent2 .RData")
output.2 <- output
n.vec.2 <- seq(from=420,to=500,by=20)

load("Parametric.Independent3 .RData")
output.3 <- output
n.vec.3 <- seq(from=520,to=600,by=20)

load("Parametric.Independent4 .RData")
output.4 <- output
n.vec.4 <- seq(from=620,to=700,by=20)

load("Parametric.Independent5 .RData")
output.5 <- output
n.vec.5 <- seq(from=720,to=800,by=20)

load("Parametric.Independent6 .RData")
output.6 <- output
n.vec.6 <- seq(from=820,to=900,by=20)

load("Parametric.Independent7 .RData")
output.7 <- output
n.vec.7 <- seq(from=920,to=1000,by=20)

load("Parametric.Independent8 .RData")
output.8 <- output
n.vec.8 <- seq(from=1020,to=1100,by=20)

load("Parametric.Independent9 .RData")
output.9 <- output
n.vec.9 <- seq(from=1120,to=1200,by=20)
############################################
###### Data Analysis for each Dataset ######
############################################
p.adj.mat.1 <- content.mat.1 <- content.adj.mat.1 <- matrix(NA,ncol=length(n.vec.1),nrow=M)

for (i in 1:M){
  p.adj.mat.1[i,] <- output.1[[i]][[1]]
  content.mat.1[i,] <- output.1[[i]][[2]]
  content.adj.mat.1[i,] <- output.1[[i]][[3]]
}

cov.vec.1 <- cov.adj.vec.1 <- rep(NA,length(n.vec.1))

for (k in 1:length(n.vec.1)){
  cov.vec.1[k] <- sum(content.mat.1[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.1[,k])))
  cov.adj.vec.1[k] <- sum(content.adj.mat.1[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.1[,k])))
}

#############################################
p.adj.mat.2 <- content.mat.2 <- content.adj.mat.2 <- matrix(NA,ncol=length(n.vec.2),nrow=M)

for (i in 1:M){
  p.adj.mat.2[i,] <- output.2[[i]][[1]]
  content.mat.2[i,] <- output.2[[i]][[2]]
  content.adj.mat.2[i,] <- output.2[[i]][[3]]
}

cov.vec.2 <- cov.adj.vec.2 <- rep(NA,length(n.vec.2))

for (k in 1:length(n.vec.2)){
  cov.vec.2[k] <- sum(content.mat.2[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.2[,k])))
  cov.adj.vec.2[k] <- sum(content.adj.mat.2[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.2[,k])))
}

#############################################
p.adj.mat.3 <- content.mat.3 <- content.adj.mat.3 <- matrix(NA,ncol=length(n.vec.3),nrow=M)

for (i in 1:M){
  p.adj.mat.3[i,] <- output.3[[i]][[1]]
  content.mat.3[i,] <- output.3[[i]][[2]]
  content.adj.mat.3[i,] <- output.3[[i]][[3]]
}

cov.vec.3 <- cov.adj.vec.3 <- rep(NA,length(n.vec.3))

for (k in 1:length(n.vec.3)){
  cov.vec.3[k] <- sum(content.mat.3[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.3[,k])))
  cov.adj.vec.3[k] <- sum(content.adj.mat.3[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.3[,k])))
}

#############################################
p.adj.mat.4 <- content.mat.4 <- content.adj.mat.4 <- matrix(NA,ncol=length(n.vec.4),nrow=M)

for (i in 1:M){
  p.adj.mat.4[i,] <- output.4[[i]][[1]]
  content.mat.4[i,] <- output.4[[i]][[2]]
  content.adj.mat.4[i,] <- output.4[[i]][[3]]
}

cov.vec.4 <- cov.adj.vec.4 <- rep(NA,length(n.vec.4))

for (k in 1:length(n.vec.4)){
  cov.vec.4[k] <- sum(content.mat.4[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.4[,k])))
  cov.adj.vec.4[k] <- sum(content.adj.mat.4[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.4[,k])))
}

#############################################
p.adj.mat.5 <- content.mat.5 <- content.adj.mat.5 <- matrix(NA,ncol=length(n.vec.5),nrow=M)

for (i in 1:M){
  p.adj.mat.5[i,] <- output.5[[i]][[1]]
  content.mat.5[i,] <- output.5[[i]][[2]]
  content.adj.mat.5[i,] <- output.5[[i]][[3]]
}

cov.vec.5 <- cov.adj.vec.5 <- rep(NA,length(n.vec.5))

for (k in 1:length(n.vec.5)){
  cov.vec.5[k] <- sum(content.mat.5[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.5[,k])))
  cov.adj.vec.5[k] <- sum(content.adj.mat.5[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.5[,k])))
}

#############################################
p.adj.mat.6 <- content.mat.6 <- content.adj.mat.6 <- matrix(NA,ncol=length(n.vec.6),nrow=M)

for (i in 1:M){
  p.adj.mat.6[i,] <- output.6[[i]][[1]]
  content.mat.6[i,] <- output.6[[i]][[2]]
  content.adj.mat.6[i,] <- output.6[[i]][[3]]
}

cov.vec.6 <- cov.adj.vec.6 <- rep(NA,length(n.vec.6))

for (k in 1:length(n.vec.6)){
  cov.vec.6[k] <- sum(content.mat.6[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.6[,k])))
  cov.adj.vec.6[k] <- sum(content.adj.mat.6[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.6[,k])))
}

#############################################
p.adj.mat.7 <- content.mat.7 <- content.adj.mat.7 <- matrix(NA,ncol=length(n.vec.7),nrow=M)

for (i in 1:M){
  p.adj.mat.7[i,] <- output.7[[i]][[1]]
  content.mat.7[i,] <- output.7[[i]][[2]]
  content.adj.mat.7[i,] <- output.7[[i]][[3]]
}

cov.vec.7 <- cov.adj.vec.7 <- rep(NA,length(n.vec.7))

for (k in 1:length(n.vec.7)){
  cov.vec.7[k] <- sum(content.mat.7[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.7[,k])))
  cov.adj.vec.7[k] <- sum(content.adj.mat.7[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.7[,k])))
}

#############################################
p.adj.mat.8 <- content.mat.8 <- content.adj.mat.8 <- matrix(NA,ncol=length(n.vec.8),nrow=M)

for (i in 1:M){
  p.adj.mat.8[i,] <- output.8[[i]][[1]]
  content.mat.8[i,] <- output.8[[i]][[2]]
  content.adj.mat.8[i,] <- output.8[[i]][[3]]
}

cov.vec.8 <- cov.adj.vec.8 <- rep(NA,length(n.vec.8))

for (k in 1:length(n.vec.8)){
  cov.vec.8[k] <- sum(content.mat.8[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.8[,k])))
  cov.adj.vec.8[k] <- sum(content.adj.mat.8[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.8[,k])))
}
#############################################
p.adj.mat.9 <- content.mat.9 <- content.adj.mat.9 <- matrix(NA,ncol=length(n.vec.9),nrow=M)

for (i in 1:M){
  p.adj.mat.9[i,] <- output.9[[i]][[1]]
  content.mat.9[i,] <- output.9[[i]][[2]]
  content.adj.mat.9[i,] <- output.9[[i]][[3]]
}

cov.vec.9 <- cov.adj.vec.9 <- rep(NA,length(n.vec.9))

for (k in 1:length(n.vec.9)){
  cov.vec.9[k] <- sum(content.mat.9[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.mat.9[,k])))
  cov.adj.vec.9[k] <- sum(content.adj.mat.9[,k] >= P , na.rm = TRUE)/(M-sum(is.na(content.adj.mat.9[,k])))
}
######################################################
######################################################
cov.org <- c(cov.vec.1 , cov.vec.2 , cov.vec.3 , cov.vec.4 , cov.vec.5,
             cov.vec.6 , cov.vec.7 , cov.vec.8 , cov.vec.9)
cov.adj <- c(cov.adj.vec.1 , cov.adj.vec.2 , cov.adj.vec.3 , cov.adj.vec.4 , cov.adj.vec.5,
             cov.adj.vec.6 , cov.adj.vec.7 , cov.adj.vec.8 , cov.adj.vec.9)
cov.org
cov.adj

N.vec <- c(n.vec.1 , n.vec.2 , n.vec.3 , n.vec.4 , n.vec.5,
           n.vec.6 , n.vec.7 , n.vec.8 , n.vec.9)

##############################################################
##############################################################
decimal.fun <- function(x) sprintf("%.2f", x)
x.var <- N.vec
y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2
###############################################
df.9095 <- data.frame(x.var , cov.org , cov.adj)
names(df.9095) <- c("Size","Unadjusted","Adjusted")

ggplot(data=df.9095)+
  geom_point(aes(x=x.var , y=df.9095$Unadjusted ,
                 col=rep("Unadj",length(x.var)),
                 shape=rep("Unadj",length(x.var))),cex=7)+
  geom_line(aes(x=x.var , y=df.9095$Unadjusted ,
                col=rep("Unadj",length(x.var))),lwd=3)+
  geom_point(aes(x=x.var , y=df.9095$Adjusted ,
                 col=rep("Adj",length(x.var)),
                 shape=rep("Adj",length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df.9095$Adjusted ,
                col=rep("Adj",length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(N.vec),max(N.vec)), 
                     labels = seq(from=min(N.vec),to=max(N.vec),by=50) , 
                     breaks = seq(from=min(N.vec),to=max(N.vec),by=50)) +
  scale_y_continuous(limits=c(0,1) , labels = decimal.fun(y.lim) ,
                     breaks = y.lim)+
  geom_hline(yintercept = (1-alpha) , col="brown" , lwd=3, lty = 2) +
  theme(legend.position=c(0.5,0.5),legend.direction="horizontal",
        legend.text=element_text(size=35),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Unadj'="red" , 'Adj'="darkblue"), 
                      labels=c("Adjusted       " , "Unadjusted")) +
  scale_shape_manual(name="Method",
                     values=c('Unadj'=16 , 'Adj'=17), 
                     labels=c("Adjusted       " , "Unadjusted")) +
  theme(axis.text.x = element_text(face="bold", size=30 , angle=90) ,
        axis.title.x = element_text(face="bold", size=30) ,
        axis.text.y = element_text(face="bold", size=30) ,
        axis.title.y = element_text(face="bold", size=30))