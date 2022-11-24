library(tolerance)
library(ggplot2)
##################################
set.seed(2020)

extreme.vec <- seq(from=-3,to=3,by=0.02)

P <- 0.90
alpha <- 0.10

N <- 1000
mean <- 0
sd <- 1

######################
###### For Mean ######
######################
exp.IF.mean.tukey.trim <- exp.IF.var.tukey.trim <- exp.IF.mean.tol.trim <- exp.IF.var.tol.trim <-
exp.IF.mean.tukey.win <- exp.IF.var.tukey.win <- exp.IF.mean.tol.win <- exp.IF.var.tol.win <- rep(NA,length(extreme.vec))

for (i in 1:length(extreme.vec)){
  extreme.value <- extreme.vec[i]
  
  set.seed(2020)
  data <- rnorm(n=N,mean=mean,sd=sd)
  data.order <- data[order(data)]
  
  x.r.tukey <- data.order[N*(1-P)/2]
  x.s.tukey <- data.order[N*(1+P)/2]
  size.tukey <- min(which(data.order >= x.s.tukey))-max(which(data.order <= x.r.tukey))+1
  
  x.r.tol <- nptol.int(data.order , alpha=alpha,P=P,side=2,method="WILKS")$"2-sided.lower"
  x.s.tol <- nptol.int(data.order , alpha=alpha,P=P,side=2,method="WILKS")$"2-sided.upper"
  size.tol <- min(which(data.order >= x.s.tol))-max(which(data.order <= x.r.tol))+1
  ### Trimmed Mean ###
  mean.trim.tukey <- sum(data.order[which(data.order >= x.r.tukey & data.order <= x.s.tukey)])/size.tukey
  mean.trim.tol <- sum(data.order[which(data.order >= x.r.tol & data.order <= x.s.tol)])/size.tol
  ### Winsorized Mean ###
  mean.win.tukey <- (which(data.order == x.r.tukey)*x.r.tukey+(N-which(data.order == x.s.tukey)+1)*x.s.tukey+
    sum(data.order[(which(data.order == x.r.tukey)+1):(which(data.order == x.s.tukey)-1)]))/N
  mean.win.tol <- (which(data.order == x.r.tol)*x.r.tol+(N-which(data.order == x.s.tol)+1)*x.s.tol+
    sum(data.order[(which(data.order == x.r.tol)+1):(which(data.order == x.s.tol)-1)]))/N
  ############
  if (extreme.value < x.r.tukey){
    exp.IF.mean.tukey.trim[i] <- x.r.tukey-mean.trim.tukey
    exp.IF.mean.tukey.win[i] <- (x.r.tukey-mean.win.tukey)-((N-size.tukey)/N)/(2*(N*(1-P)/2)/N)
    
    exp.IF.var.tukey.trim[i] <- ((x.r.tukey-mean)^2-sd^2)/(size.tukey/N)
    exp.IF.var.tukey.win[i] <- ((x.r.tukey-mean)^2-sd^2)/(size.tukey/N)+(((N-size.tukey)/N)/(2*(N*(1-P)/2)/N))^2
  } else if (extreme.value > x.s.tukey){
    exp.IF.mean.tukey.trim[i] <- x.s.tukey-mean.trim.tukey
    exp.IF.mean.tukey.win[i] <- (x.s.tukey-mean.win.tukey)+((N-size.tukey)/N)/(2*(N*(1-P)/2)/N)
    
    exp.IF.var.tukey.trim[i] <- ((x.s.tukey-mean)^2-sd^2)/(size.tukey/N)
    exp.IF.var.tukey.win[i] <- ((x.s.tukey-mean)^2-sd^2)/(size.tukey/N)+(((N-size.tukey)/N)/(2*(N*(1-P)/2)/N))^2
    
  } else {
    exp.IF.mean.tukey.trim[i] <- extreme.value-mean.trim.tukey
    exp.IF.mean.tukey.win[i] <- (extreme.value-mean.win.tukey)
    
    exp.IF.var.tukey.trim[i] <- ((extreme.value-mean)^2-sd^2)/(size.tukey/N)
    exp.IF.var.tukey.win[i] <- ((extreme.value-mean)^2-sd^2)/(size.tukey/N)
  }
  ############
  if (extreme.value < x.r.tol){
    exp.IF.mean.tol.trim[i] <- x.r.tol-mean.trim.tol
    exp.IF.mean.tol.win[i] <- (x.r.tol-mean.win.tol)-((N-size.tol)/N)/((length(which(data.order >= x.s.tol))+length(which(data.order <= x.r.tol)))/N)
    
    exp.IF.var.tol.trim[i] <- ((x.r.tol-mean)^2-sd^2)/(size.tol/N)
    exp.IF.var.tol.win[i] <- ((x.r.tol-mean)^2-sd^2)/(size.tol/N)+(((N-size.tukey)/N)/(2*(N*(1-P)/2)/N))^2
    
  } else if (extreme.value > x.s.tol){
    exp.IF.mean.tol.trim[i] <- x.s.tol-mean.trim.tol
    exp.IF.mean.tol.win[i] <- (x.s.tol-mean.win.tol)+((N-size.tol)/N)/((length(which(data.order >= x.s.tol))+length(which(data.order <= x.r.tol)))/N)
    
    exp.IF.var.tol.trim[i] <- ((x.s.tol-mean)^2-sd^2)/(size.tol/N)
    exp.IF.var.tol.win[i] <- ((x.s.tol-mean)^2-sd^2)/(size.tol/N)+(((N-size.tukey)/N)/(2*(N*(1-P)/2)/N))^2
  } else {
    exp.IF.mean.tol.trim[i] <- extreme.value-mean.trim.tol
    exp.IF.mean.tol.win[i] <- (extreme.value-mean.win.tol)
    
    exp.IF.var.tol.trim[i] <- ((extreme.value-mean)^2-sd^2)/(size.tol/N)
    exp.IF.var.tol.win[i] <- ((extreme.value-mean)^2-sd^2)/(size.tol/N)
  }
}

###################
###### Plots ######
###################
################
### For mean ###
################
decimal.fun <- function(x) sprintf("%.2f", x)

df.trim.mean <- data.frame(extreme.vec , exp.IF.mean.tukey.trim , exp.IF.mean.tol.trim)
names(df.trim.mean) <- c("Extreme","Tukey","Tolerance")

ggplot(data=df.trim.mean)+
  geom_point(aes(x=extreme.vec , y=df.trim.mean$Tukey ,
                 col=factor(rep("Tukey",length(extreme.vec))),
                 shape=factor(rep("Tukey",length(extreme.vec)))),
             cex=7,alpha=1)+
  geom_point(aes(x=extreme.vec , y=df.trim.mean$Tolerance ,
                 col=factor(rep("Tol",length(extreme.vec))),
                 shape=factor(rep("Tol",length(extreme.vec)))),
             cex=7,alpha=1) +
  xlab("X") + ylab("Empirical Influence Function") +
  scale_x_continuous(limits=c(min(extreme.vec),max(extreme.vec)), 
                     labels = decimal.fun(seq(from=-3,to=3,by=0.2)) , 
                     breaks = seq(from=-3,to=3,by=0.2)) +
  scale_y_continuous(limits=c(-2,2) , 
                     labels = decimal.fun(seq(from=-2,to=2,by=0.4)) ,
                     breaks = seq(from=-2,to=2,by=0.4))+
  theme(legend.position=c(0.75,0.1),legend.direction="horizontal",
        legend.text=element_text(size=35),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Tukey'="red" , 'Tol'="darkblue"), 
                      labels=c("Tolerance       " , "Tukey")) +
  scale_shape_manual(name="Method",
                      values=c('Tukey'=16 , 'Tol'=8), 
                      labels=c("Tolerance       " , "Tukey")) +
  theme(axis.text.x = element_text(face="bold", size=30 , angle=90) ,
        axis.title.x = element_text(face="bold", size=30) ,
        axis.text.y = element_text(face="bold", size=30) ,
        axis.title.y = element_text(face="bold", size=30))

###################################################
df.win.mean <- data.frame(extreme.vec , exp.IF.mean.tukey.win , exp.IF.mean.tol.win)
names(df.win.mean) <- c("Extreme","Tukey","Tolerance")

ggplot(data=df.win.mean)+
  geom_point(aes(x=extreme.vec , y=df.win.mean$Tukey ,
                 col=factor(rep("Tukey",length(extreme.vec))),
                 shape=factor(rep("Tukey",length(extreme.vec)))),
             cex=7,alpha=1)+
  geom_point(aes(x=extreme.vec , y=df.win.mean$Tolerance ,
                 col=factor(rep("Tol",length(extreme.vec))),
                 shape=factor(rep("Tol",length(extreme.vec)))),
             cex=7,alpha=1) +
  xlab("X") + ylab("Empirical Influence Function") +
  scale_x_continuous(limits=c(min(extreme.vec),max(extreme.vec)), 
                     labels = decimal.fun(seq(from=-3,to=3,by=0.2)) , 
                     breaks = seq(from=-3,to=3,by=0.2)) +
  scale_y_continuous(limits=c(-3,3) , 
                     labels = decimal.fun(seq(from=-3,to=3,by=0.4)) ,
                     breaks = seq(from=-3,to=3,by=0.4))+
  theme(legend.position=c(0.75,0.1),legend.direction="horizontal",
        legend.text=element_text(size=35),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Tukey'="red" , 'Tol'="darkblue"), 
                      labels=c("Tolerance       " , "Tukey")) +
  scale_shape_manual(name="Method",
                     values=c('Tukey'=16 , 'Tol'=8), 
                     labels=c("Tolerance       " , "Tukey")) +
  theme(axis.text.x = element_text(face="bold", size=30 , angle=90) ,
        axis.title.x = element_text(face="bold", size=30) ,
        axis.text.y = element_text(face="bold", size=30) ,
        axis.title.y = element_text(face="bold", size=30))

####################
### For Variance ###
####################
df.trim.var <- data.frame(extreme.vec , exp.IF.var.tukey.trim , exp.IF.var.tol.trim)
names(df.trim.var) <- c("Extreme","Tukey","Tolerance")

ggplot(data=df.trim.var)+
  geom_point(aes(x=extreme.vec , y=df.trim.var$Tukey ,
                 col=factor(rep("Tukey",length(extreme.vec))),
                 shape=factor(rep("Tukey",length(extreme.vec)))),
             cex=7,alpha=1)+
  geom_point(aes(x=extreme.vec , y=df.trim.var$Tolerance ,
                 col=factor(rep("Tol",length(extreme.vec))),
                 shape=factor(rep("Tol",length(extreme.vec)))),
             cex=7,alpha=1) +
  xlab("X") + ylab("Empirical Influence Function") +
  scale_x_continuous(limits=c(min(extreme.vec),max(extreme.vec)), 
                     labels = decimal.fun(seq(from=-3,to=3,by=0.2)) , 
                     breaks = seq(from=-3,to=3,by=0.2)) +
  scale_y_continuous(limits=c(-2,5) , 
                     labels = decimal.fun(seq(from=-2,to=5,by=0.5)) ,
                     breaks = seq(from=-2,to=5,by=0.5))+
  theme(legend.position=c(0.5,0.8),legend.direction="horizontal",
        legend.text=element_text(size=35),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Tukey'="red" , 'Tol'="darkblue"), 
                      labels=c("Tolerance       " , "Tukey")) +
  scale_shape_manual(name="Method",
                     values=c('Tukey'=16 , 'Tol'=8), 
                     labels=c("Tolerance       " , "Tukey")) +
  theme(axis.text.x = element_text(face="bold", size=30 , angle=90) ,
        axis.title.x = element_text(face="bold", size=30) ,
        axis.text.y = element_text(face="bold", size=30) ,
        axis.title.y = element_text(face="bold", size=30))

########################
df.win.var <- data.frame(extreme.vec , exp.IF.var.tukey.win , exp.IF.var.tol.win)
names(df.win.var) <- c("Extreme","Tukey","Tolerance")

ggplot(data=df.win.var)+
  geom_point(aes(x=extreme.vec , y=df.win.var$Tukey ,
                 col=factor(rep("Tukey",length(extreme.vec))),
                 shape=factor(rep("Tukey",length(extreme.vec)))),
             cex=7,alpha=1)+
  geom_point(aes(x=extreme.vec , y=df.win.var$Tolerance ,
                 col=factor(rep("Tol",length(extreme.vec))),
                 shape=factor(rep("Tol",length(extreme.vec)))),
             cex=7,alpha=1) +
  xlab("X") + ylab("Empirical Influence Function") +
  scale_x_continuous(limits=c(min(extreme.vec),max(extreme.vec)), 
                     labels = decimal.fun(seq(from=-3,to=3,by=0.2)) , 
                     breaks = seq(from=-3,to=3,by=0.2)) +
  scale_y_continuous(limits=c(-2,5) , 
                     labels = decimal.fun(seq(from=-2,to=5,by=0.5)) ,
                     breaks = seq(from=-2,to=5,by=0.5))+
  theme(legend.position=c(0.5,0.8),legend.direction="horizontal",
        legend.text=element_text(size=35),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c('Tukey'="red" , 'Tol'="darkblue"), 
                      labels=c("Tolerance       " , "Tukey")) +
  scale_shape_manual(name="Method",
                     values=c('Tukey'=16 , 'Tol'=8), 
                     labels=c("Tolerance       " , "Tukey")) +
  theme(axis.text.x = element_text(face="bold", size=30 , angle=90) ,
        axis.title.x = element_text(face="bold", size=30) ,
        axis.text.y = element_text(face="bold", size=30) ,
        axis.title.y = element_text(face="bold", size=30))
