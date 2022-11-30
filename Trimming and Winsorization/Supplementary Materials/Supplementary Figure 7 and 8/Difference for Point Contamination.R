library(ggplot2)

#########################################
###### Regular Distribution N(0,1) ######
###### Point Mass at x='outlier'   ######
#########################################
set.seed(2020)
epsilon.vec <- seq(from = 0.000 , to = 1 , by = 0.001)

N <- 1000
mean <- 0
sd <- 1
outlier <- 10 # Set outlier (contamination) at 0,2,4,10.

P <- 0.90
alpha <- 0.10
#########################################################
IF.vec.raw.win.tol <- sd.IF.vec.raw.win.tol <- mean.vec.win.tol <- mean.cont.vec.win.tol <- sd.vec.win.tol <- sd.cont.vec.win.tol <-
IF.vec.raw.win.tukey <- sd.IF.vec.raw.win.tukey <- mean.vec.win.tukey <- mean.cont.vec.win.tukey <- sd.vec.win.tukey <- sd.cont.vec.win.tukey <-
  rep(NA,length(epsilon.vec))

for (i in 1:length(epsilon.vec)){
  epsilon <- epsilon.vec[i] 
  
  data <- rnorm(n=N , mean = mean ,sd=sd)
  data.order <- data[order(data)]
  
  data.contaminate <- c(data[sample(1:N,(1-epsilon)*N,replace=FALSE)], rep(outlier , N*epsilon))
  data.contaminate.order <- data.contaminate[order(data.contaminate)]
  ################################################################
  ### Tolerance Based Winsorization Mean without Contamination ###
  ################################################################
  tol.lower <- nptol.int(data.order , alpha=alpha,P=P,side=2,method="WILKS")$"2-sided.lower"
  tol.upper <- nptol.int(data.order , alpha=alpha,P=P,side=2,method="WILKS")$"2-sided.upper"
  tol.lower.pos <- max(which(data.order <= tol.lower))
  tol.upper.pos <- min(which(data.order >= tol.upper))
  
  data.win.tol <- c(rep(tol.lower , (tol.lower.pos-1)) , 
                data.order[tol.lower.pos : tol.upper.pos] ,
                rep(tol.upper , (N-tol.upper.pos)))
  
  mean.win.tol <- mean(data.win.tol)
  sd.win.tol <- sd(data.win.tol)
  #############################################################
  ### Tolerance Based Winsorization Mean with Contamination ###
  #############################################################
  tol.lower.cont <- nptol.int(data.contaminate.order , alpha=alpha,P=P,side=2,method="WILKS")$"2-sided.lower"
  tol.upper.cont <- nptol.int(data.contaminate.order , alpha=alpha,P=P,side=2,method="WILKS")$"2-sided.upper"
  
  tol.lower.pos.cont <- max(which(data.contaminate.order <= tol.lower.cont))
  tol.upper.pos.cont <- min(which(data.contaminate.order >= tol.upper.cont))
  
  data.cont.win.tol <- c(rep(tol.lower.cont , (tol.lower.pos.cont-1)) , 
                     data.contaminate.order[tol.lower.pos.cont : tol.upper.pos.cont] ,
                     rep(tol.upper.cont , (N-tol.upper.pos.cont)))
  
  mean.cont.win.tol <- mean(data.cont.win.tol)
  sd.cont.win.tol <- sd(data.cont.win.tol)
  ########################################################
  ### Tukey's Winsorization Mean without Contamination ###
  ########################################################
  tukey.lower.pos <- ((1-P)/2)*N ### Only for even N, for simulation purpose.
  tukey.upper.pos <- N-tukey.lower.pos
  
  tukey.lower <- data.order[tukey.lower.pos+1]
  tukey.upper <- data.order[tukey.upper.pos]
  
  data.win.tukey <- c(rep(tukey.lower , tukey.lower.pos) , 
                data.order[(tukey.lower.pos+1) : tukey.upper.pos] ,
                rep(tukey.upper , tukey.lower.pos))
  
  mean.win.tukey <- mean(data.win.tukey)
  sd.win.tukey <- sd(data.win.tukey)
  #####################################################
  ### Tukey's Winsorization Mean with Contamination ###
  #####################################################
  tukey.lower.pos.cont <- ((1-P)/2)*N ### Only for even N, for simulation purpose.
  tukey.upper.pos.cont <- N-tukey.lower.pos
  
  tukey.lower.cont <- data.contaminate.order[tukey.lower.pos.cont+1]
  tukey.upper.cont <- data.contaminate.order[tukey.upper.pos.cont]
  
  data.cont.win.tukey <- c(rep(tukey.lower.cont , tukey.lower.pos.cont) , 
                  data.contaminate.order[(tukey.lower.pos.cont+1) : tukey.upper.pos.cont] ,
                  rep(tukey.upper.cont , tukey.lower.pos.cont))
  
  mean.cont.win.tukey <- mean(data.cont.win.tukey)
  sd.cont.win.tukey <- sd(data.cont.win.tukey)
  ###### Influence Function ######
  mean.vec.win.tol[i] <- mean.win.tol
  mean.cont.vec.win.tol[i] <- mean.cont.win.tol
  sd.vec.win.tol[i] <- sd.win.tol
  sd.cont.vec.win.tol[i] <- sd.cont.win.tol
  IF.vec.raw.win.tol[i] <- (mean.cont.win.tol-mean.win.tol)
  sd.IF.vec.raw.win.tol[i] <- (sd.cont.win.tol-sd.win.tol)
  
  mean.vec.win.tukey[i] <- mean.win.tukey
  mean.cont.vec.win.tukey[i] <- mean.cont.win.tukey
  sd.vec.win.tukey[i] <- sd.win.tukey
  sd.cont.vec.win.tukey[i] <- sd.cont.win.tukey
  IF.vec.raw.win.tukey[i] <- (mean.cont.win.tukey-mean.win.tukey)
  sd.IF.vec.raw.win.tukey[i] <- (sd.cont.win.tukey-sd.win.tukey)
}

################################
###### Influence Function ######
################################
### Recall that Influence function is a derivative of the difference between functionals
### of contaminated and non-contaminated date.
### Therefore, the empirical influence function is the slope of
### IF.vec.raw.win.tukey and IF.vec.raw.win.tol, respectively.

IF.vec.win.tol <- IF.vec.win.tukey <- sd.IF.vec.win.tol <- sd.IF.vec.win.tukey <- rep(NA,(length(IF.vec.raw.win.tol)-1))

for (k in 1:length(IF.vec.win.tukey)){
  IF.vec.win.tol[k] <- (IF.vec.raw.win.tol[k]-IF.vec.raw.win.tol[k+1])/epsilon.vec[1]
  IF.vec.win.tukey[k] <- (IF.vec.raw.win.tukey[k]-IF.vec.raw.win.tukey[k+1])/epsilon.vec[1]
}

# IF.vec.win.tol[c(43,916)] <- IF.vec.win.tukey[c(50,949)] <- NA

IF.raw.data <- data.frame(cbind(IF.vec.raw.win.tol , IF.vec.raw.win.tukey))
sd.IF.raw.data <- data.frame(cbind(sd.IF.vec.raw.win.tol , sd.IF.vec.raw.win.tukey))
###################################################
### Plot of Differences Between Two Functionals ###
###################################################
decimal.fun <- function(x) sprintf("%.2f", x)

ggplot()+
  geom_point(data=IF.raw.data , aes(x=rev(epsilon.vec) , y=IF.raw.data$IF.vec.raw.win.tukey , 
                                    color=factor(rep("Tukey" , length(epsilon.vec))),
                                    shape=factor(rep("Tukey" , length(epsilon.vec)))) ,
            alpha=1 , cex=5) +
  geom_point(data=IF.raw.data , aes(x=rev(epsilon.vec) , y=IF.raw.data$IF.vec.raw.win.tol ,
                                    color=factor(rep("Tol" , length(epsilon.vec))),
                                    shape=factor(rep("Tol" , length(epsilon.vec)))) ,
            alpha=0.25 , cex=5) +
  xlab("Contaminate Proportion") + ylab("Difference Between Functionals") +
  scale_x_continuous(limits=c(0,1), 
                     labels = decimal.fun(rev(seq(from=0.000 , to=1 , by=0.05))),
                     breaks = seq(from=0.000 , to=1 , by=0.05)) +
  theme(axis.text.x = element_text(face="bold", size=25 , angle=90) ,
        axis.title.x = element_text(face="bold", size=25) ,
        axis.text.y = element_text(face="bold", size=25) ,
        axis.title.y = element_text(face="bold", size=25)) +
  theme(legend.position=c(0.2,0.15),legend.direction="horizontal",
        legend.text=element_text(size=25),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                     values=c("Tukey"="red" , "Tol"="darkblue"), 
                     labels=c("Tolerance       " , "Tukey ")) +
  scale_shape_manual(name="Method",
                      values=c("Tukey"=20 , "Tol"=8), 
                      labels=c("Tolerance       " , "Tukey "))
#####################################
ggplot()+
  geom_point(data=sd.IF.raw.data , aes(x=rev(epsilon.vec) , y=sd.IF.raw.data$sd.IF.vec.raw.win.tukey , 
                                    color=factor(rep("Tukey" , length(epsilon.vec))),
                                    shape=factor(rep("Tukey" , length(epsilon.vec)))) ,
             alpha=1 , cex=5) +
  geom_point(data=sd.IF.raw.data , aes(x=rev(epsilon.vec) , y=sd.IF.raw.data$sd.IF.vec.raw.win.tol ,
                                    color=factor(rep("Tol" , length(epsilon.vec))),
                                    shape=factor(rep("Tol" , length(epsilon.vec)))) ,
             alpha=0.25 , cex=5) +
  xlab("Contaminate Proportion") + ylab("Difference Between Functionals") +
  scale_x_continuous(limits=c(0,1), 
                     labels = decimal.fun(rev(seq(from=0.000 , to=1 , by=0.05))),
                     breaks = seq(from=0.000 , to=1 , by=0.05)) +
  theme(axis.text.x = element_text(face="bold", size=25 , angle=90) ,
        axis.title.x = element_text(face="bold", size=25) ,
        axis.text.y = element_text(face="bold", size=25) ,
        axis.title.y = element_text(face="bold", size=25)) +
  theme(legend.position=c(0.5,0.15),legend.direction="horizontal",
        legend.text=element_text(size=25),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_colour_manual(name="Method",
                      values=c("Tukey"="red" , "Tol"="darkblue"), 
                      labels=c("Tolerance       " , "Tukey ")) +
  scale_shape_manual(name="Method",
                     values=c("Tukey"=20 , "Tol"=8), 
                     labels=c("Tolerance       " , "Tukey "))

