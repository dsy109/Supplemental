###############################################################
###############################################################
# Example 1: EFH Analysis
###############################################################
###############################################################

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()
system("ls")

library(tolerance)
library(xtable)
library(sqldf)
library(ggplot2)
library(metR) # To use geom_text_contour()
source("Sim_k_factor.R")

data1 = data.frame(
  Z1 = c(8.0,8.5,7.5,8.0,8.0,8.0,8.0,7.5,7.5,8.5,8.5,8.0,8.0),
  Z2 = c(25,20,25,20,20,20,20,15,20,25,15,20,15),
  Y  = c(4526.56,6512.32,5264.65,6865.48,6795.53,6656.69,6845.12,
         5156.98,5951.18,4156.69,6965.74,6685.48,5985.56) 
)

data1$X1 = 2*(data1$Z1 - median(data1$Z1))/(max(data1$Z1)-min(data1$Z1))
data1$X2 = 2*(data1$Z2 - median(data1$Z2))/(max(data1$Z2)-min(data1$Z2))

mod = lm(Y ~ X1 + X2 + X1:X2 + I(X1^2)+I(X2^2), data=data1)
summary(mod)


# Coded design matrix of the slected model

data1

X1X2 = data1[,"X1"]*data1[,"X2"]
X1_sq = data1[,"X1"]^2
X2_sq = data1[,"X2"]^2

X = x0 = as.matrix(cbind(intercept=1,data1[,c("X1","X2")],X1_sq,X2_sq,X1X2))
s=summary(mod)$sigma

aps_k = PS_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x0,adjust = T)
ch_k = C2_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x0,h_sq=NULL)$kfactors

pw_k = rep(NA,dim(X)[1])
for (i in 1:dim(X)[1]){
  a = regtol.int(reg=mod, new.x = as.data.frame(t(X[i,])), 
                 side = 2, alpha = 0.05, P = 0.75)
  lb = sqldf("SELECT * FROM a WHERE y IS NULL")[,"2-sided.lower"]
  ub = sqldf("SELECT * FROM a WHERE y IS NULL")[,"2-sided.upper"]
  pw_k[i] = round((ub-lb)/(2*s),2)
}

pw_k
aps_k
ch_k

design_pts = X[,c("X1","X2")]

Y_hat = predict(mod,  data.frame(X1=X[,"X1"],X2=X[,"X2"]))

pw  = cbind(pw_k = pw_k, lb= round(Y_hat - pw_k*s,2), ub= round(Y_hat + pw_k*s,2))
aps = cbind(aps_k = aps_k, lb= round(Y_hat - aps_k*s,2), ub= round(Y_hat + aps_k*s,2))
ch  = cbind(ch_k = ch_k, lb= round(Y_hat - ch_k*s,2), ub= round(Y_hat + ch_k*s,2))

table1 = cbind(design_pts,Y_hat=Y_hat,pw,aps,ch)
table = as.data.frame(table1[which(!duplicated(table1)),])
print(xtable(table,digits=c(0,2,2,1,2,1,1,2,1,1,2,1,1)),
      include.rownames=F)

# Contour Plot

A.vec = seq(-1,1,by=0.1)
B.vec = seq(-1,1,by=0.1)

d3 = expand.grid(X1=A.vec,X2=B.vec)


d3$pred = predict(mod,  data.frame(X1=d3[,1],X2=d3[,2]))

x1 = d3[,c("X1","X2")]
x1 = as.matrix(cbind(intercept=1,x1,
                     X1_sq=I(x1$X1^2), X2_sq=I(x1$X2^2),
                     X1X2 = x1$X1*x1$X2))

k_factors = aps_k = PS_k_factor(P.target=0.75,alpha.target=0.05,X=X,x0=x1,adjust = T)
d2 = cbind(d3,k_factors)
d2$lb = d2$pred - d2$k_factors*s
d2$ub = d2$pred + d2$k_factors*s

d2$lb_not_bounded =  d2$lb
d2$lb = ifelse(d2$lb < 0,0,d2$lb)

# Contour 

v <- ggplot(d2, aes(X1, X2, z = pred))
v2 = v + geom_contour(bins = 1000) +
  geom_raster(aes(fill = pred)) +
  scale_fill_gradientn(limits=c(2000,10000),colours=c("green","red")) +
  # scale_fill_gradientn(limits = c(90,105),colours=c("green","red")) +
  geom_contour(colour = "black") +
  geom_text_contour(aes(z = pred)) +
  ggtitle("Response surface") +
  theme(plot.title = element_text(hjust = 0.5))
print(v2)

png(filename = "ccd_example1_contour.png")
v2
dev.off()


# upper bound    

v <- ggplot(d2, aes(X1, X2, z = ub))
v_upper = v + geom_contour(bins = 1000) +
  geom_raster(aes(fill = ub)) +
  scale_fill_gradientn(limits=c(2000,10000),colours=c("green","red")) +
  # scale_fill_gradientn(limits = c(90,105),colours=c("green","red")) +
  geom_contour(colour = "black") +
  geom_text_contour(aes(z = ub)) +
  ggtitle("APS TI upper bound") +
  theme(plot.title = element_text(hjust = 0.5))
print(v_upper)

png(filename = "ccd_example1_ub.png")
v_upper
dev.off()

# lower bound    

v <- ggplot(d2, aes(X1, X2, z = lb))
v_lower = v + geom_contour(bins = 1000) +
  geom_raster(aes(fill = lb)) +
  scale_fill_gradientn(limits=c(2000,10000),colours=c("green","red")) +
  # scale_fill_gradientn(limits = c(90,105),colours=c("green","red")) +
  geom_contour(colour = "black") +
  geom_text_contour(aes(z = lb)) +
  ggtitle("APS TI lower bound") +
  theme(plot.title = element_text(hjust = 0.5))

print(v_lower)

png(filename = "ccd_example1_lb.png")
v_lower
dev.off()

