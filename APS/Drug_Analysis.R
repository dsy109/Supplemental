###############################################################
###############################################################
# Example 2: Drug Formulation Analysis
###############################################################
###############################################################

install.packages("tolerance")
install.packages("ggplot2")

source("Sim_k_factor.R")
library(tolerance)
library(xtable)
library(sqldf)
library(ggplot2)
library(metR) # To use geom_text_contour()

mydata = data.frame(
  run = 1:17,
  X1 = c(0.05,0.04,0.05,0.06,0.05,0.06,0.04,0.04,0.04,0.05,0.05,0.06,0.05,0.06,0.05,0.05,0.05),
  X2 = c(0.02,0.02,0.03,0.01,0.02,0.02,0.02,0.03,0.01,0.03,0.02,0.02,0.01,0.03,0.02,0.02,0.01),
  X3 = c(2.0,1.5,1.5,2.0,2.0,2.5,2.5,2.0,2.0,2.5,2.0,1.5,1.5,2.0,2.0,2.0,2.5),
  Y1 = c(108.00,112.60,138.23,126.27,102.33,104.60,70.80,89.08,89.72,88.55,106.67,156.37,139.87,134.07,110.13,109.70,91.10),
  PDI = c(0.221,0.220,0.218,0.207,0.220,0.211,0.208,0.225,0.205,0.223,0.220,0.210,0.237,0.256,0.231,0.238,0.213),
  Y2 = c(1.95,1.92,2.92,0.96,1.94,1.94,1.95,2.94,0.98,2.95,1.95,1.95,0.96,2.96,1.96,1.95,0.96),
  Y3 = c(97.38,96.24,97.47,95.52,97.18,96.99,97.70,97.91,97.60,98.28,97.61,97.55,95.59,98.58,98.01,97.62,96.24)
  )

coded_data = mydata
coded_data$X1 = 2*(mydata$X1 - mean(c(min(mydata$X1),max(mydata$X1))))/(max(mydata$X1)-min(mydata$X1))
coded_data$X2 = 2*(mydata$X2 - mean(c(min(mydata$X2),max(mydata$X2))))/(max(mydata$X2)-min(mydata$X2))
coded_data$X3 = 2*(mydata$X3 - mean(c(min(mydata$X3),max(mydata$X3))))/(max(mydata$X3)-min(mydata$X3))

# mod1 = lm(Y1~X1+X2+X3+I(X2^2)+I(X3^2),data=coded_data); summary(mod1)
# mod2 = lm(Y2~X1+X2+X3+X1*X2+X1*X3+I(X2^2)+I(X3^2),data=coded_data); summary(mod2)

mod3 = lm(Y3~X1+X2+X3+X1*X2+X1*X3+I(X3^2),data=coded_data)
names(anova(mod3))
anova(mod3)
s=summary(mod3)$sigma

# Original scale
summary(lm(Y3~X1+X2+X3+X1*X2+X1*X3+I(X3^2),data=mydata))
# Variable selecton
summary(lm(Y3~X1+X2+X3+X1*X2+X1*X3+X2*X3+I(X1^2)+I(X2^2)+I(X3^2),data=coded_data))
summary(lm(Y3~X1+X2+X3+X1*X2+X1*X3+X2*X3+I(X2^2)+I(X3^2),data=coded_data))
summary(lm(Y3~X1+X2+X3+X1*X2+X1*X3+I(X2^2)+I(X3^2),data=coded_data))
summary(lm(Y3~X1+X2+X3+X1*X2+X1*X3+I(X3^2),data=coded_data))

# Coded design matrix of the slected model

X1X2 = coded_data[,"X1"]*coded_data[,"X2"]
X1X3 = coded_data[,"X1"]*coded_data[,"X3"]
X3sq = coded_data[,"X3"]^2
X = cbind(intercept=1,coded_data[,c("X1","X2","X3")],X3sq,X1X2,X1X3)

print(xtable(X,digits=1), floating=FALSE, tabular.environment="bmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

X = x0 = as.matrix(X)
s=summary(mod3)$sigma

aps_k = PS_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x0,adjust = T)
ch_k = C2_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x0,h_sq=NULL)$kfactors

pw_k = rep(NA,dim(X)[1])
for (i in 1:dim(X)[1]){
  a = regtol.int(reg=mod3, new.x = as.data.frame(t(X[i,])), 
                 side = 2, alpha = 0.05, P = 0.75)
  lb = sqldf("SELECT * FROM a WHERE y IS NULL")[,"2-sided.lower"]
  ub = sqldf("SELECT * FROM a WHERE y IS NULL")[,"2-sided.upper"]
  pw_k[i] = round((ub-lb)/(2*s),2)
}

design_pts = X[,c("X1","X2","X3")]
Y_hat = round(X%*%as.matrix(mod3$coefficients),2)

pw  = cbind(pw_k = pw_k, lb= round(Y_hat - pw_k*s,2), ub= round(Y_hat + pw_k*s,2))
aps = cbind(aps_k = aps_k, lb= round(Y_hat - aps_k*s,2), ub= round(Y_hat + aps_k*s,2))
ch  = cbind(ch_k = ch_k, lb= round(Y_hat - ch_k*s,2), ub= round(Y_hat + ch_k*s,2))

table1 = cbind(design_pts,Y_hat=Y_hat,pw,aps,ch)
table = as.data.frame(table1[which(!duplicated(table1)),])
print(xtable(table,digits=c(1,0,0,0,2,2,2,2,2,2,2,2,2,2)),
      include.rownames=F)

# Contour Plot

A.vec = seq(-1,1,by=0.2)
B.vec = seq(-1,1,by=0.2)
C.vec = seq(-1,1,by=0.2)
d3 = expand.grid(X1=A.vec,X2=B.vec,X3=C.vec)


# d3$pred = predict(mod3,  data.frame(X1=d3[,1],X2=d3[,2],X3=d3[,3]))
# d3$lb   = predict(mod_lb,data.frame(X1=d3[,1],X2=d3[,2],X3=d3[,3]))　
# d3$ub   = predict(mod_ub,data.frame(X1=d3[,1],X2=d3[,2],X3=d3[,3]))


draw_contour = function(C.value){
  
  # C.value = -1
  x1 = subset(d3,X3==C.value)[,c("X1","X2","X3")]
  x1 = as.matrix(cbind(intercept=1,x1,X3sq=I(x1$X3^2),X1X2=x1$X1*x1$X2,X1X3=x1$X1*x1$X3))
  k_factors = aps_k = PS_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x1,adjust = T)
  
  
  strings = paste("select X1,X2,X3 from d3 where X3= ",C.value,sep="")
  d2 = cbind(x1,s,k_factors,pred=predict(mod3, as.data.frame(sqldf(strings))))
  d2 = as.data.frame(d2)
  d2$lb = d2$pred - d2$s * d2$k_factors     
  d2$ub = d2$pred + d2$s * d2$k_factors 
  
  d2$ub_not_bounded =  d2$ub
  d2$ub = ifelse(d2$ub > 100, 100,d2$ub)
  d2$ub = round(d2$ub,3)
  d2$lb = round(d2$lb,3)
  
  # sqldf("select * from d2 where X1=-1 and X2= -1 and X3=-1") 
  # sqldf("select * from d2 where X1=-0.8 and X2= -1 and X3=-1") 
  
  v <- ggplot(d2, aes(X1, X2, z = pred))
  v2 = v + geom_contour(colour = "black") +
    geom_raster(aes(fill = pred)) +
    scale_fill_gradientn(limits = c(93,100),colours=c("green","red")) +
    # geom_contour(colour = "black") +
    # geom_text_contour(aes(z = pred)) +
    geom_contour(breaks =seq(90,100,0.5),colour = "black") +
    geom_text_contour(breaks =seq(90,100,0.5)) +
    ggtitle(paste("Response surface X3=",C.value)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(v2)
  
  name = paste("../images/rsm_example_bbd_gupta/gupta_c_",C.value,"_rs.pdf",sep="")
  pdf(file=name, bg="transparent")
  print(v2)
  dev.off()
  
  v <- ggplot(d2, aes(X1, X2, z = lb))
  v2 = v + geom_contour(colour = "black") +
    geom_raster(aes(fill = lb)) +
    scale_fill_gradientn(limits = c(93,100),colours=c("green","red")) +
    # geom_contour(colour = "black") +
    # geom_text_contour(aes(z = lb)) +
    geom_contour(breaks =seq(90,100,0.5),colour = "black") +
    geom_text_contour(breaks =seq(90,100,0.5)) +
    ggtitle(paste("APS TI lower bound X3=",C.value)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(v2)
  
  name = paste("../images/rsm_example_bbd_gupta/gupta_c_",C.value,"_lb.pdf",sep="")
  pdf(file=name, bg="transparent")
  print(v2)
  dev.off()
  
  v <- ggplot(d2, aes(X1, X2, z = ub))
  v2 = v + geom_contour(colour = "black") +
    geom_raster(aes(fill = ub)) +
    scale_fill_gradientn(limits = c(93,100),colours=c("green","red")) +
    # geom_contour(colour = "black") +
    # geom_text_contour(aes(z = ub)) +
    geom_contour(breaks =seq(90,100,0.5),colour = "black") +
    geom_text_contour(breaks =seq(90,100,0.5)) +
    ggtitle(paste("APS TI upper bound X3=",C.value)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(v2)
  
  name = paste("../images/rsm_example_bbd_gupta/gupta_c_",C.value,"_ub.pdf",sep="")
  pdf(file=name, bg="transparent")
  print(v2)
  dev.off()
  
}　# end of draw_contour()


draw_contour(C.value=-1)
draw_contour(C.value= 0)
draw_contour(C.value= 1)
  


draw_contour_original = function(C.value){
  
  # C.value = -1
  x1 = subset(d3,X3==C.value)[,c("X1","X2","X3")]
  x1 = as.matrix(cbind(intercept=1,x1,X3sq=I(x1$X3^2),X1X2=x1$X1*x1$X2,X1X3=x1$X1*x1$X3))
  k_factors = aps_k = PS_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x1,adjust = T)
  
  
  strings = paste("select X1,X2,X3 from d3 where X3= ",C.value,sep="")
  d2 = cbind(x1,s,k_factors,pred=predict(mod3, as.data.frame(sqldf(strings))))
  d2 = as.data.frame(d2)
  d2$lb = d2$pred - d2$s * d2$k_factors     
  d2$ub = d2$pred + d2$s * d2$k_factors 
  
  d2$ub_not_bounded =  d2$ub
  d2$ub = ifelse(d2$ub > 100, 100,d2$ub)
  d2$ub = round(d2$ub,3)
  d2$lb = round(d2$lb,3)
  
  d2$X1_original = 0.01*d2$X1 + 0.05
  d2$X2_original = 0.01*d2$X2 + 0.02
  d2$X3_original = 0.5*d2$X3 + 2
  C.value_original = 0.5*C.value + 2
  C.value_original_string = gsub("\\.", "_",C.value_original)
  
  # sqldf("select * from d2 where X1=-1 and X2= -1 and X3=-1") 
  # sqldf("select * from d2 where X1=-0.8 and X2= -1 and X3=-1") 
  
  v <- ggplot(d2, aes(X1_original, X2_original, z = pred))
  v2 = v + geom_contour(colour = "black") +
    geom_raster(aes(fill = pred)) +
    scale_fill_gradientn(limits = c(93,100),colours=c("green","red")) +
    # geom_contour(colour = "black") +
    # geom_text_contour(aes(z = pred)) +
    geom_contour(breaks =seq(90,100,0.5),colour = "black") +
    geom_text_contour(breaks =seq(90,100,0.5)) +
    xlab("X1") +
    ylab("X2") +
    ggtitle(paste("Response surface X3=",C.value_original)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(v2)
  
  name = paste("/gupta_c_",C.value_original_string,"_rs.pdf",sep="")
  pdf(file=name, bg="transparent")
  print(v2)
  dev.off()
  
  v <- ggplot(d2, aes(X1_original, X2_original, z = lb))
  v2 = v + geom_contour(colour = "black") +
    geom_raster(aes(fill = lb)) +
    scale_fill_gradientn(limits = c(93,100),colours=c("green","red")) +
    # geom_contour(colour = "black") +
    # geom_text_contour(aes(z = lb)) +
    geom_contour(breaks =seq(90,100,0.5),colour = "black") +
    geom_text_contour(breaks =seq(90,100,0.5)) +
    xlab("X1") +
    ylab("X2") +
    ggtitle(paste("APS TI lower bound X3=",C.value_original)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(v2)
  
  name = paste("/gupta_c_",C.value_original_string,"_lb.pdf",sep="")
  pdf(file=name, bg="transparent")
  print(v2)
  dev.off()
  
  v <- ggplot(d2, aes(X1_original, X2_original, z = ub))
  v2 = v + geom_contour(colour = "black") +
    geom_raster(aes(fill = ub)) +
    scale_fill_gradientn(limits = c(93,100),colours=c("green","red")) +
    # geom_contour(colour = "black") +
    # geom_text_contour(aes(z = ub)) +
    geom_contour(breaks =seq(90,100,0.5),colour = "black") +
    geom_text_contour(breaks =seq(90,100,0.5)) +
    xlab("X1") +
    ylab("X2") +
    ggtitle(paste("APS TI upper bound X3=",C.value_original)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(v2)
  
  name = paste("/gupta_c_",C.value_original_string,"_ub.pdf",sep="")
  pdf(file=name, bg="transparent")
  print(v2)
  dev.off()
  
}　# end of draw_contour()

draw_contour_original(C.value=-1)
draw_contour_original(C.value= 0)
draw_contour_original(C.value= 1)
