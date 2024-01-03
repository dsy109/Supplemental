#### Non-Small Cell Lung Cancer Data

nscl <- read.delim("https://raw.githubusercontent.com/dsy109/Supplemental/main/ZIG/ST001269_AN002109.txt")
Stage <- read.delim("https://raw.githubusercontent.com/dsy109/Supplemental/main/ZIG/stage.txt",header=FALSE)
df1 <- nscl[1:1320,]
df2 <- nscl[1322:2641,2:9]
colnames(df2) <- as.character(nscl[1321,2:9])
nnn <- data.frame(df1,df2)
for(i in 2:96) nnn[,i] <- as.numeric(nnn[,i])
nnn[,2:96] <- matrix(nnn[,2:96])

semicont <- function(x){
  ind0 <- (x==0)
  pi0 <- mean(ind0)
  x.p <- x[!ind0]
  if(length(x.p)<=1){
    best <- "Degenerate"
    logn.theta <- gamma.theta <- c(NA, NA)
    logn.ll <- gamma.ll <- NA
  } else{
    logn.theta <- try(suppressWarnings(fitdistr(x.p,dlnorm,start=list(meanlog=mean(log(x.p)),sdlog=sd(log(x.p))))$estimate),silent=TRUE)
    gamma.theta <- try(suppressWarnings(fitdistr(x.p,dgamma,list(shape=mean(x.p)^2/var(x.p),scale=var(x.p)/mean(x.p)))$estimate),silent=TRUE)
    #  if(class(gamma.theta)=="try-error") gamma.theta <- try(suppressWarnings(fitdistr(x.p,dgamma,list(shape=mean(x.p)^2/var(x.p),scale=var(x.p)/mean(x.p)/10))$estimate),silent=TRUE)
    if(class(logn.theta)=="try-error") logn.theta <- c(meanlog=mean(log(x.p)),sdlog=sd(log(x.p)))
    if(class(gamma.theta)=="try-error") gamma.theta <- c(shape=mean(x.p)^2/var(x.p),scale=var(x.p)/mean(x.p))
    logn.ll <- sum(ind0)*log(pi0) + sum((1-pi0)*dlnorm(x.p,meanlog=logn.theta[1],sdlog=logn.theta[2],log=TRUE))
    gamma.ll <- sum(ind0)*log(pi0) + sum((1-pi0)*dgamma(x.p,shape=gamma.theta[1],scale=gamma.theta[2],log=TRUE))
    best <- ifelse(max(logn.ll,gamma.ll)==logn.ll,"ZI.logn","ZI.gamma")
  }
  out <- data.frame(meanlog=logn.theta[1],sdlog=logn.theta[2],ZI.logn=logn.ll,
                    shape=gamma.theta[1],scale=gamma.theta[2],ZI.gamma=gamma.ll,
                    Best=best)
  out
}

semicont(nnn[,2])

nnn.Ceramides <- nnn[which(nnn[,99]=="Ceramides"),2:96]
nnn.Cholesterol_Esters <- nnn[which(nnn[,99]=="Cholesterol Esters"),2:96]
nnn.Contamination <- nnn[which(nnn[,99]=="Contamination"),2:96]
nnn.Diacylglycerols <- nnn[which(nnn[,99]=="Diacylglycerols"),2:96]
nnn.Eicosanoids <- nnn[which(nnn[,99]=="Eicosanoids"),2:96]
nnn.Fatty_acids <- nnn[which(nnn[,99]=="Fatty acids"),2:96]
nnn.Glycerophospholipids <- nnn[which(nnn[,99]=="Glycerophospholipids"),2:96]
nnn.Lysoglycerophospholipids <- nnn[which(nnn[,99]=="Lysoglycerophospholipids"),2:96]
nnn.Monoacylglycerols <- nnn[which(nnn[,99]=="Monoacylglycerols"),2:96]
nnn.Sphingolipids <- nnn[which(nnn[,99]=="Sphingolipids"),2:96]
nnn.Standard <- nnn[which(nnn[,99]=="Standard"),2:96]
nnn.Triacylglycerols <- nnn[which(nnn[,99]=="Triacylglycerols"),2:96]

out.Ceramides <- sapply(1:95,function(i) semicont(nnn.Ceramides[,i]))
out.Cholesterol_Esters <- sapply(1:95,function(i) semicont(nnn.Cholesterol_Esters[,i]))
out.Contamination <- sapply(1:95,function(i) semicont(nnn.Contamination[,i]))
out.Diacylglycerol <- sapply(1:95,function(i) semicont(nnn.Diacylglycerols[,i]))
out.Eicosanoids <- sapply(1:95,function(i) semicont(nnn.Eicosanoids[,i]))
out.Fatty_acids <- sapply(1:95,function(i) semicont(nnn.Fatty_acids[,i]))
out.Glycerophospholipids <- sapply(1:95,function(i) semicont(nnn.Glycerophospholipids[,i]))
out.Lysoglycerophospholipids <- sapply(1:95,function(i) semicont(nnn.Lysoglycerophospholipids[,i]))
out.Monoacylglycerols <- sapply(1:95,function(i) semicont(nnn.Monoacylglycerols[,i]))
out.Sphingolipids <- sapply(1:95,function(i) semicont(nnn.Sphingolipids[,i]))
out.Standard <- sapply(1:95,function(i) semicont(nnn.Standard[,i]))
out.Triacylglycerols <- sapply(1:95,function(i) semicont(nnn.Triacylglycerols[,i]))

table(unlist(out.Ceramides[7,]))
table(unlist(out.Cholesterol_Esters[7,]))
table(unlist(out.Contamination[7,]))
table(unlist(out.Diacylglycerol[7,]))
table(unlist(out.Eicosanoids[7,]))
table(unlist(out.Fatty_acids[7,]))
table(unlist(out.Glycerophospholipids[7,]))
table(unlist(out.Lysoglycerophospholipids[7,]))
table(unlist(out.Monoacylglycerols[7,]))
table(unlist(out.Sphingolipids[7,]))
table(unlist(out.Standard[7,]))
table(unlist(out.Triacylglycerols[7,]))

sum(2*abs(unlist(out.Ceramides[3,])-unlist(out.Ceramides[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Cholesterol_Esters[3,])-unlist(out.Cholesterol_Esters[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Contamination[3,])-unlist(out.Contamination[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Diacylglycerol[3,])-unlist(out.Diacylglycerol[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Eicosanoids[3,])-unlist(out.Eicosanoids[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Fatty_acids[3,])-unlist(out.Fatty_acids[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Glycerophospholipids[3,])-unlist(out.Glycerophospholipids[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Lysoglycerophospholipids[3,])-unlist(out.Lysoglycerophospholipids[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Monoacylglycerols[3,])-unlist(out.Monoacylglycerols[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Sphingolipids[3,])-unlist(out.Sphingolipids[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Standard[3,])-unlist(out.Standard[6,]))<=2,na.rm=T)
sum(2*abs(unlist(out.Triacylglycerols[3,])-unlist(out.Triacylglycerols[6,]))<=2,na.rm=T)

i <- 23
nnn.profile <- nnn.Monoacylglycerols[,i]
out.profile <- out.Monoacylglycerols[,i]

XX0 <- nnn.profile[which(nnn.profile==0)]
XX <- nnn.profile[which(nnn.profile>0)]
n.X <- length(nnn.profile)

ggplot(data.frame(y=nnn.profile), aes(sample = y)) +
  stat_qq(distribution = qlnorm, dparams = out.profile[1:2], size=2) +
  stat_qq_line(distribution = qlnorm, dparams = out.profile[1:2], line.p=c(.1,.9)) +
  stat_qq(distribution = qgamma, dparams = out.profile[4:5],col=2, size=2, pch=17) +
  stat_qq_line(distribution = qgamma, dparams = out.profile[4:5], line.p=c(.1,.9),col=2) +
  ggtitle(paste("Q-Q Plot for Subject",names(nnn)[i+1])) + theme(text = element_text(size = 13)) +
  geom_point(aes(x=50000,y=70000),col=2,pch=17) +   geom_point(aes(x=50000,y=85000),col=1,pch=19) +
  annotate(geom="text",label="ZIG",x=54000,y=70000,col=2) + annotate(geom="text",label="ZILN",x=54000,y=85000,col=1) +
  ylab("Sample Quantiles") + xlab("Theoretical Quantiles")

out <- vector("list",95)
for(i in 1:95){
  set.seed(i)
  out[[i]] <- suppressWarnings(semicont.TI(nnn.Monoacylglycerols[,i],P=.95,alpha=.05,N=10000))
  print(i)
}

early <- normal <- late <- NULL
for(i in 1:95){
  cancer.stage <- Stage[which(Stage[,1]==i),2]
  if(cancer.stage=="Early"){
    early <- rbind(early,unlist(out[[i]])[-c(3,5,8,10:13)])
  } else if(cancer.stage=="Normal"){
    normal <- rbind(normal,unlist(out[[i]])[-c(3,5,8,10:13)])
  } else late <- rbind(late,unlist(out[[i]])[-c(3,5,8,10:13)])
}

all.ints <- rbind(data.frame(early,Stage="Early"),data.frame(late,Stage="Late"),data.frame(normal,Stage="Normal"))
all.ints$Stage <- as.factor(all.ints$Stage)

ggplot(all.ints,aes(x=ZIG.TI.95.,y=1:95, col=Stage, pch=Stage)) + geom_point(cex=2.2) +
  geom_segment(all.ints, mapping=aes(x=ZIG.CI.2.5.,xend=ZIG.CI.97.5.,y=1:95,yend=1:95,col=Stage),lwd=1.3) +
  geom_segment(all.ints, mapping=aes(x=0,xend=ZIG.TI.95.,y=1:95,yend=1:95,col=Stage),lwd=.1) +
  theme(text = element_text(size = 15)) + ylab("Subject") + xlab("Exosomal Lipid Profile Measurement") +
  ggtitle("(0.95,0.95) Upper Tolerance Limits and 95% Confidence Intervals (ZIG)")

out[[i]]





