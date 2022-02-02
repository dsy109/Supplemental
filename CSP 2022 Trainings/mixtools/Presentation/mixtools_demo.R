################################################################
################################################################
################################################################
### R Code for "Methods and Applications of Finite Mixture 
### Models, with Computing Demonstrations Using the R Package 
### mixtools" by D. S. Young
################################################################
################################################################
################################################################

library(mixtools)
library(BSDA)
library(GGally)
library(plotly)
library(mixreg)

################################################################
################################################################
### Gaussian Mixture Models
################################################################
################################################################

################################################################
### Example: Quasars Data
################################################################

quasars <- data.frame(y=read.table("http://astrostatistics.psu.edu/datasets/QSO_absorb.txt",nrows=104,skip=1)[,2])

set.seed(100)

out1 <- lapply(1:30,function(i) normalmixEM(quasars$y,epsilon=1e-06,k=2,maxit=5000))

out1 <- out1[[which.max(sapply(1:30,function(i) out1[[i]]$loglik))]]
out1[2:4]

out1.se <- boot.se(out1, B=500)
out1.se[c(2,4,6)]

m <- ggplot(quasars,aes(x=y)) + geom_histogram(aes(y=..density..),bins=9,color="black",fill="gray") +
  theme(text = element_text(size = 20)) + xlab("Velocity (km/s)") + ylab("Density") + ggtitle("Si IV 1394 Data") 
m
mix.fun <- function(x) out1$lambda[1]*dnorm(x,out1$mu[1],out1$sigma[1])+
  out1$lambda[2]*dnorm(x,out1$mu[2],out1$sigma[2])
m + stat_function(fun=mix.fun,color="red")

################################################################
### Example: Hidalgo Stamp Data
################################################################

data(Stamp)
attach(Stamp)
y <- thickness

set.seed(100)

out2 <- lapply(1:30,function(i) normalmixEM(thickness,epsilon=1e-06,k=4,maxit=5000))

out2 <- out2[[which.max(sapply(1:30,function(i) out2[[i]]$loglik))]]
out2[2:4]

out2.se <- boot.se(out2, B=500)
out2.se[c(2,4,6)]

m <- ggplot(Stamp,aes(x=thickness)) + geom_histogram(aes(y=..density..),color="black",fill="gray",bins=70) +
  theme(text = element_text(size = 20)) + xlab("Thickness (mm)") + ylab("Density") + ggtitle("Hidalgo Stamp Data") 
m
mix.fun2 <- function(x) out2$lambda[1]*dnorm(x,out2$mu[1],out2$sigma[1])+
  out2$lambda[2]*dnorm(x,out2$mu[2],out2$sigma[2])+
  out2$lambda[3]*dnorm(x,out2$mu[3],out2$sigma[3])+
  out2$lambda[4]*dnorm(x,out2$mu[4],out2$sigma[4])
m + stat_function(fun=mix.fun2,color="red")

################################################################
### Example: Diffuse Large B-Cell Lymphoma Data
################################################################

DLBCL <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/mixtools/Datasets/DLBCL.txt", header = TRUE)

set.seed(240)
DLBCLSample <- DLBCL[sample(1:nrow(DLBCL), 500), ]
out3 <- lapply(1:10,function(i) mvnormalmixEM(DLBCLSample, k = 2, verb = TRUE))

out3 <- out3[[which.max(sapply(1:10,function(i) out3[[i]]$loglik))]]
out3[2:4]

out3.se <- boot.se(out3, B=20)
out3.se[c(2,4,6)]

str(out3)
post <- as.factor(apply(out3$posterior, 1, which.max))
DLBCLSample2 <- data.frame(DLBCLSample, Component = post)

### A few visualizations of the results.

ggpairs(DLBCLSample2, columns = 1:3, mapping = aes(alpha = 0.8)) + ggtitle("DLBCL Sample") + 
  theme(text = element_text(size = 15))
ggpairs(DLBCLSample2, columns = 1:3, mapping = aes(col = Component, pch = Component, 
                                                   alpha = 0.8)) + ggtitle("Estimated Two-Component Mixture of Trivariate Normals") + 
  theme(text = element_text(size = 15))

fig <- plot_ly(DLBCLSample2, x = ~CD3, y = ~CD5, z = ~CD19, color = ~Component, colors = c("#F8766D", 
                                                                                           "#00BFC4"), size = 10)
fig <- fig %>% add_markers()
fig

################################################################
################################################################
### Mixtures of Regressions
################################################################
################################################################

################################################################
### Example: Carbon Dioxide Data
################################################################

data(CO2data)

ggplot(data = CO2data, aes(GNP, CO2, label=country)) + geom_point() + geom_label() +
  labs(title = "Carbon Dioxide Emissions", y = expression(CO[2])) + 
  theme(text = element_text(size = 13))

set.seed(100)
out4 <- lapply(1:30,function(i) regmixEM(y=CO2data$CO2, x=CO2data$GNP, k = 2))

out4 <- out4[[which.max(sapply(1:30,function(i) out4[[i]]$loglik))]]
out4[3:5]

out4.se <- boot.se(out4, B=500)
out4.se[c(2,4,6)]

CO2_df <- data.frame(CO2data,Component=as.character(apply(out4$posterior,1,which.max)))

ggplot(data = CO2_df, aes(GNP, CO2, col = Component, pch = Component)) + geom_point() + 
  labs(title = "Carbon Dioxide Emissions", y = expression(CO[2])) + 
  geom_abline(slope = out4$beta[2,1], intercept = out4$beta[1,1], col = 2) + 
  geom_abline(slope = out4$beta[2,2], intercept = out4$beta[1,2], col = 3) + 
  scale_color_manual("", values = c(2, 3), guide = F) + theme(text = element_text(size = 13), legend.position = "none")

################################################################
### Example: Aphids Data
################################################################

data(aphids)
attach(aphids)

aphids.df <- data.frame(x = aphRel, y = plntsInf)
ggplot(data = aphids.df, aes(x, y)) + geom_point() + labs(title = "Aphids Data", 
                                                          x = "Number of Aphids in Batch", y = "Number of Infected Plants") + theme(text = element_text(size = 20))

### Fit a 2-component mixture of linear regressions to the aphids data.

set.seed(100)
out5 <- regmixEM(plntsInf, aphRel, k = 2)
mix.lab <- as.factor(apply(out5$posterior, 1, which.max) + 1)
out5[3:6]

### Scatterplot of the data with the estimate from the 2-component mixture of
### regressions overlaid.

aphids.df1 <- cbind(aphids.df, group = mix.lab)
ggplot(data = aphids.df1, aes(x, y, col = group, pch = group)) + geom_point() + 
  labs(title = "Aphids Data (Parametric Mixture Fit)", 
       x = "Number of Aphids in Batch", y = "Number of Infected Plants") + 
  geom_abline(slope = out5$beta[,1][2], intercept = out5$beta[, 1][1], col = 2) + 
  geom_abline(slope = out5$beta[,2][2], intercept = out5$beta[, 2][1], col = 3) + 
  scale_color_manual("", values = c(2,3), guide = F) + 
  theme(text = element_text(size = 20), legend.position = "none")

### Estimate a semiparametric mixtures of regressions assuming only a symmetric
### error density.  The bandwidth will adapt by default.

set.seed(100)
spfit1 <- spregmix(plntsInf ~ aphRel, bw = 0.1, data = aphids, verb = TRUE)

mix.lab2 <- as.factor(apply(spfit1$posterior, 1, which.max) + 1)
aphids.df2 <- cbind(aphids.df, group = mix.lab2)
ggplot(data = aphids.df2, aes(x, y, col = group, pch = group)) + geom_point() + 
  labs(title = "Aphids Data (Semiparametric Mixture Fit; Symmetric Errors Assumption)", 
       x = "Number of Aphids in Batch", y = "Number of Infected Plants") + 
  geom_abline(slope = spfit1$beta[,  1][2], intercept = spfit1$beta[,1][1], col = 2) + 
  geom_abline(slope = spfit1$beta[,  2][2], intercept = spfit1$beta[, 2][1], col = 3) + 
  scale_color_manual("", values = c(2, 3), guide = F) + theme(text = element_text(size = 20), legend.position = "none")

### Estimate a semiparametric mixtures of regressions without assuming a symmetric
### error density.

set.seed(100)
spfit2 <- spregmix(plntsInf ~ aphRel, bw = 0.1, data = aphids, verb = TRUE, symm = FALSE)

mix.lab3 <- as.factor(apply(spfit2$posterior, 1, which.max) + 1)
aphids.df3 <- cbind(aphids.df, group = mix.lab3)
ggplot(data = aphids.df3, aes(x, y, col = group, pch = group)) + geom_point() + 
  labs(title = "Aphids Data (Semiparametric Mixture Fit; No Symmetric Errors Assumption)", 
       x = "Number of Aphids in Batch", y = "Number of Infected Plants") + 
  geom_abline(slope = spfit2$beta[, 1][2], intercept = spfit2$beta[, 1][1], col = 2) + 
  geom_abline(slope = spfit2$beta[,  2][2], intercept = spfit2$beta[, 2][1], col = 3) + 
  scale_color_manual("", values = c(2,  3), guide = F) + theme(text = element_text(size = 20), legend.position = "none")

list(out5$beta, spfit1$beta, spfit2$beta)

### Look at the nonparametric KD-based estimate of the error densities fit above,
### as well as the parametric normal density from the original fully parametric
### mixture of linear regressions fit.

z <- seq(min(c(spfit1$density.x, spfit2$density.x)), max(c(spfit1$density.x, spfit2$density.x)), 
         length = 200)
dens <- data.frame(x = c(spfit1$density.x, spfit2$density.x, z), dens = c(spfit1$density.y, 
                                                                          spfit2$density.y, dnorm(z, sd = sqrt((spfit1$np.stdev)^2 + spfit1$bandwidth^2))), 
                   Density = as.factor(c(rep("Symmetric Assumption", length(spfit1$density.x)), 
                                         rep("No Symmetric Assumption", length(spfit2$density.x)), rep("Normal Density", 
                                                                                                       length(z)))))
ggplot(data = dens, aes(x = x, y = dens, color = Density)) + geom_line() + ylab("Density") + 
  theme(text = element_text(size = 10)) + ggtitle("Kernel Density Estimates of Errors")

################################################################
################################################################
### Other Parametric Mixtures
################################################################
################################################################

################################################################
### Example: Whole Genome Duplication Data
################################################################

JLOV <- unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/mixtools/Datasets/JLOV.txt", header = FALSE))

set.seed(100)
out6 <- gammamixEM(x=JLOV,k=6,mom.start=TRUE,fix.alpha=F,epsilon=1e-3,verb=T,maxit=25000)
out6[2:3]

df <- data.frame(Dist=JLOV,Comp1=out6$lambda[1]*dgamma(JLOV,shape=out6$gamma.pars[1,1],scale=out6$gamma.pars[2,1]),
                 Comp2=out6$lambda[2]*dgamma(JLOV,shape=out6$gamma.pars[1,2],scale=out6$gamma.pars[2,2]),
                 Comp3=out6$lambda[3]*dgamma(JLOV,shape=out6$gamma.pars[1,3],scale=out6$gamma.pars[2,3]),
                 Comp4=out6$lambda[4]*dgamma(JLOV,shape=out6$gamma.pars[1,4],scale=out6$gamma.pars[2,4]),
                 Comp5=out6$lambda[5]*dgamma(JLOV,shape=out6$gamma.pars[1,5],scale=out6$gamma.pars[2,5]),
                 Comp6=out6$lambda[6]*dgamma(JLOV,shape=out6$gamma.pars[1,6],scale=out6$gamma.pars[2,6]))

p <- ggplot(df) +
  geom_histogram(aes(x = Dist, y = ..density..),
                 fill = "black", color = "black", bins=250) 
p + geom_line(data = df, aes(x = Dist, y = Comp1), color = "#E69F00",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp2), color = "#56B4E9",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp3), color = "#009E73",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp4), color = "#F0E442",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp5), color = "#0072B2",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp6), color = "#D55E00",size=1) + 
  theme(text = element_text(size = 20)) + ylab("Density") + xlab("Synonymous Distances") +
  ggtitle("Pereskia Aculeata (Cactaceae)") + scale_y_continuous(limits = c(0, 4))  

################################################################
### Example: Earthquake Data
################################################################

earthquakes <- unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/mixtools/Datasets/earthquakes.txt", header = FALSE))

earth.counts <- as.numeric(table(earthquakes))

set.seed(1)
out7 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=6))
out7 <- out7[which.max(sapply(1:20,function(i)  out7[[i]]$loglik ))][[1]]

x <- 0:210
df.e <- data.frame(Dist=earth.counts)
df <- data.frame(Dist=x,Comp1=out7$lambda[1]*dpois(x,lambda=exp(out7$beta[1])),
                 Comp2=out7$lambda[2]*dpois(x,lambda=exp(out7$beta[2])),
                 Comp3=out7$lambda[3]*dpois(x,lambda=exp(out7$beta[3])),
                 Comp4=out7$lambda[4]*dpois(x,lambda=exp(out7$beta[4])),
                 Comp5=out7$lambda[5]*dpois(x,lambda=exp(out7$beta[5])),
                 Comp6=out7$lambda[6]*dpois(x,lambda=exp(out7$beta[6])))

p <- ggplot(df.e) +
  geom_histogram(aes(x = Dist, y = ..density..),
                 fill = "black", color = "black") 
p + geom_line(data = df, aes(x = Dist, y = Comp1), color = "#E69F00",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp2), color = "#56B4E9",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp3), color = "#009E73",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp4), color = "#F0E442",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp5), color = "#0072B2",size=1) + 
  geom_line(data = df, aes(x = Dist, y = Comp6), color = "#D55E00",size=1) + 
  theme(text = element_text(size = 20)) + ylab("Density") + xlab("Number of Earthquakes") +
  ggtitle("Annual Number of Earthquakes (1900-2021)") 

################################################################
### Example: Infant Habituation Data
################################################################

infant <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/mixtools/Datasets/infant.txt", header = TRUE)

y.1 <- lapply(1:nrow(infant),function(i) infant[i,-1])
y.1 <- lapply(1:length(y.1),function(i) matrix(log(as.numeric(y.1[[i]])),ncol=1))
x.1 <- lapply(1:length(y.1), function(i) cbind(1,1:11,exp(-(1:11-1)^2)) )

set.seed(100)
out8 <- regmixEM.mixed(y.1,x.1,addintercept.random=FALSE,verb=T,epsilon=1e-3,ar.1=F)

posteriors <- apply(out8$posterior.z,1,which.max)

Infant <- data.frame(RT = log(unlist(y.1)), Time=1:11, Subject=as.factor(rep(infant$Subject,each=11)), Component=as.factor(rep(posteriors,each=11)))
Infant1 <- subset(Infant, Component==1)
Infant2 <- subset(Infant, Component==2)

ggplot(Infant1,aes(x=Time,y=RT)) + geom_point(aes(color=Subject)) + geom_line(aes(color=Subject)) +
  theme(text = element_text(size = 20),legend.position="none") +
  ylab("Reaction Time") + ggtitle("Infant Habituation Data (First Component)") +
  stat_smooth(method = 'lm', formula = y ~ x+I(exp(-x^2)), color=1, se= T) + ylim(c(1.7,2.6))

ggplot(Infant2,aes(x=Time,y=RT)) + geom_point(aes(color=Subject)) + geom_line(aes(color=Subject)) +
  theme(text = element_text(size = 20),legend.position="none") +
  ylab("Reaction Time") + ggtitle("Infant Habituation Data (Second Component)") +
  stat_smooth(method = 'lm', formula = y ~ x+I(exp(-x^2)), color=1, se= T) + ylim(c(1.7,2.6))

################################################################
################################################################
### Determining the Number of Components
################################################################
################################################################

################################################################
### Example: Earthquake Data (ctd.)
################################################################

set.seed(1)
bic2 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=2))
bic2 <- bic2[which.max(sapply(1:20,function(i)  bic2[[i]]$loglik ))][[1]]
bic3 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=3))
bic3 <- bic3[which.max(sapply(1:20,function(i)  bic3[[i]]$loglik ))][[1]]
bic4 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=4))
bic4 <- bic4[which.max(sapply(1:20,function(i)  bic4[[i]]$loglik ))][[1]]
bic5 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=5))
bic5 <- bic5[which.max(sapply(1:20,function(i)  bic5[[i]]$loglik ))][[1]]
bic6 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=6))
bic6 <- bic6[which.max(sapply(1:20,function(i)  bic6[[i]]$loglik ))][[1]]
bic7 <- lapply(1:20, function(i) poisregmixEM(earth.counts,x=rep(1,length(earth.counts)),addintercept=FALSE,k=7))
bic7 <- bic7[which.max(sapply(1:20,function(i)  bic7[[i]]$loglik ))][[1]]

-2*bic2$loglik+log(122)*3
-2*bic3$loglik+log(122)*5
-2*bic4$loglik+log(122)*7
-2*bic5$loglik+log(122)*9
-2*bic6$loglik+log(122)*11
-2*bic7$loglik+log(122)*13

################################################################
### Example: Aphids Data (ctd.)
################################################################

set.seed(100)
mod.sel <- regmixmodel.sel(x = aphids$aphRel, y = aphids$plntsInf, k = 4, type = "fixed")
mod.sel

################################################################
### Example: Carbon Dioxide Data (ctd.)
################################################################

set.seed(100)
boot.out <- boot.comp(y=CO2data$CO2, x=CO2data$GNP, max.comp=4, B=100, sig=0.05,
                      mix.type="regmix", hist=FALSE)

################################################################
################################################################
### Visualizing Estimated Mixture Models
################################################################
################################################################

################################################################
### Example: Hidalgo Stamp Data (ctd.)
################################################################

source("https://raw.githubusercontent.com/dsy109/mixtools/Updated-Functions/Updated-Functions/plotly.mixEM.R")

plot(out2, density = TRUE) #Current
plotly.mixEM(out2, density = TRUE) #Updated

################################################################
### Example: Quasars Data (ctd.)
################################################################

source("https://raw.githubusercontent.com/dsy109/mixtools/Updated-Functions/Updated-Functions/plotly.mixturegram.R")

set.seed(100)
selection=function(i,data,rep=30){
  out=replicate(30,normalmixEM(data,epsilon=1e-06,k=i,maxit=5000),simplify=F)
  counts=lapply(1:30,function(j) table(apply(out[[j]]$posterior,1,which.max)))
  counts.length=sapply(counts,length)
  counts.min=sapply(counts,min)
  counts.test=(counts.length!=i)|(counts.min<5)
  if(sum(counts.test)>0&sum(counts.test)<rep) out=out[!counts.test]
  l=unlist(lapply(out,function(x) x$loglik))
  tmp=out[[which.max(l)]]
}

all.out=lapply(2:6,selection,data=y)
pmbs=lapply(1:length(all.out),function(i) all.out[[i]]$post)
mix.out=plotly.mixturegram(y,pmbs,method="pca",all.n=FALSE,id.con=NULL,title="Mixturegram (Quasars Data)")


