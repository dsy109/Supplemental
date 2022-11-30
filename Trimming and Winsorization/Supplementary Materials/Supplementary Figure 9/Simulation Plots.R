library(ggplot2)
setwd("/Users/kcheng/Desktop/Supplementary Figure 9/Simulation Plots")
#########################
###### Cauchy 9090 ######
#########################
load("CoverageCauchy9090 .RData")
Cauchy9090 <- output
P <- 0.90
alpha <- 0.10
sample.size <- seq(from=50 , to=1000 , by=10)
###############################################
decimal.fun <- function(x) sprintf("%.2f", x)

x.var <- sample.size
df<-data.frame(x.var, Cauchy9090)
names(df) <- c("N","Trim" , "Winsorize" , "Trim Tol" , "Winsorize Tol")

y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2

################
ggplot(data=df)+
  geom_point(aes(x=x.var , y=df$Trim ,
                 col=rep("Tukey" , length(x.var)),
                 shape=rep("Tukey" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$Trim ,
                col=rep("Tukey" , length(x.var))),lwd=3) +
  geom_point(aes(x=x.var , y=df$`Trim Tol` ,
                 col=rep("Tol" , length(x.var)),
                 shape=rep("Tol" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$`Trim Tol` ,
                col=rep("Tol" , length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(sample.size),max(sample.size)), 
                     labels = seq(from=50,to=1000,by=50) , 
                     breaks = seq(from=50,to=1000,by=50)) +
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

#########################
###### Cauchy 9095 ######
#########################
load("CoverageCauchy9095 .RData")
Cauchy9095 <- output
P <- 0.90
alpha <- 0.05
sample.size <- seq(from=50 , to=1000 , by=10)
###############################################
decimal.fun <- function(x) sprintf("%.2f", x)

x.var <- sample.size
df<-data.frame(x.var, Cauchy9095)
names(df) <- c("N","Trim" , "Winsorize" , "Trim Tol" , "Winsorize Tol")

y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2

################
ggplot(data=df)+
  geom_point(aes(x=x.var , y=df$Trim ,
                 col=rep("Tukey" , length(x.var)),
                 shape=rep("Tukey" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$Trim ,
                col=rep("Tukey" , length(x.var))),lwd=3) +
  geom_point(aes(x=x.var , y=df$`Trim Tol` ,
                 col=rep("Tol" , length(x.var)),
                 shape=rep("Tol" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$`Trim Tol` ,
                col=rep("Tol" , length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(sample.size),max(sample.size)), 
                     labels = seq(from=50,to=1000,by=50) , 
                     breaks = seq(from=50,to=1000,by=50)) +
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


#########################
###### Cauchy 9590 ######
#########################
load("CoverageCauchy9590 .RData")
Cauchy9590 <- output
P <- 0.95
alpha <- 0.10
sample.size <- seq(from=50 , to=1000 , by=10)
###############################################
decimal.fun <- function(x) sprintf("%.2f", x)

x.var <- sample.size
df<-data.frame(x.var, Cauchy9590)
names(df) <- c("N","Trim" , "Winsorize" , "Trim Tol" , "Winsorize Tol")

y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2

################
ggplot(data=df)+
  geom_point(aes(x=x.var , y=df$Trim ,
                 col=rep("Tukey" , length(x.var)),
                 shape=rep("Tukey" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$Trim ,
                col=rep("Tukey" , length(x.var))),lwd=3) +
  geom_point(aes(x=x.var , y=df$`Trim Tol` ,
                 col=rep("Tol" , length(x.var)),
                 shape=rep("Tol" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$`Trim Tol` ,
                col=rep("Tol" , length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(sample.size),max(sample.size)), 
                     labels = seq(from=50,to=1000,by=50) , 
                     breaks = seq(from=50,to=1000,by=50)) +
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

#########################
###### Cauchy 9595 ######
#########################
load("CoverageCauchy9595 .RData")
Cauchy9595 <- output
P <- 0.95
alpha <- 0.05
sample.size <- seq(from=50 , to=1000 , by=10)
###############################################
decimal.fun <- function(x) sprintf("%.2f", x)

x.var <- sample.size
df<-data.frame(x.var, Cauchy9595)
names(df) <- c("N","Trim" , "Winsorize" , "Trim Tol" , "Winsorize Tol")

y.lim <- round(seq(from=0,to=1,by=0.10),2)
y.size <- rep(12,length(y.lim))
y.size[which(y.lim == (1-alpha))] <- 15
y.col <- rep(1,length(y.size))
y.col[which(y.lim == (1-alpha))] <- 2

################
ggplot(data=df)+
  geom_point(aes(x=x.var , y=df$Trim ,
                 col=rep("Tukey" , length(x.var)),
                 shape=rep("Tukey" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$Trim ,
                col=rep("Tukey" , length(x.var))),lwd=3) +
  geom_point(aes(x=x.var , y=df$`Trim Tol` ,
                 col=rep("Tol" , length(x.var)),
                 shape=rep("Tol" , length(x.var))),cex=7) +
  geom_line(aes(x=x.var , y=df$`Trim Tol` ,
                col=rep("Tol" , length(x.var))),lwd=3) +
  xlab("Sample Size") + ylab("Coverage Probability") +
  scale_x_continuous(limits=c(min(sample.size),max(sample.size)), 
                     labels = seq(from=50,to=1000,by=50) , 
                     breaks = seq(from=50,to=1000,by=50)) +
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
