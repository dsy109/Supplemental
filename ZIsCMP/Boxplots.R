library(ggplot2)

n <- 100
#n <- 250

all.ll.df <- lapply(1:12,function(i) rbind(data.frame(ll=ZIP.ll[,i],dist="ZIP",metric="loglikelihood"),
                                           data.frame(ll=ZIG.ll[,i],dist="ZIG",metric="loglikelihood"),
                                           data.frame(ll=ZINB.ll[,i],dist="ZINB",metric="loglikelihood"),
                                           data.frame(ll=ZIGP.ll[,i],dist="ZIGP",metric="loglikelihood"),
                                           data.frame(ll=ZICMP.ll[,i],dist="ZICMP",metric="loglikelihood"),
                                           data.frame(ll=ZIsCMP2.ll[,i],dist="ZIsCMP2",metric="loglikelihood"),
                                           data.frame(ll=ZIsCMP3.ll[,i],dist="ZIsCMP3",metric="loglikelihood")))

all.aic.df <- lapply(1:12,function(i) rbind(data.frame(ll=ZIP.AIC[,i],dist="ZIP",metric="AIC"),
                                           data.frame(ll=ZIG.AIC[,i],dist="ZIG",metric="AIC"),
                                           data.frame(ll=ZINB.AIC[,i],dist="ZINB",metric="AIC"),
                                           data.frame(ll=ZIGP.AIC[,i],dist="ZIGP",metric="AIC"),
                                           data.frame(ll=ZICMP.AIC[,i],dist="ZICMP",metric="AIC"),
                                           data.frame(ll=ZIsCMP2.AIC[,i],dist="ZIsCMP2",metric="AIC"),
                                           data.frame(ll=ZIsCMP3.AIC[,i],dist="ZIsCMP3",metric="AIC")))

all.bic.df <- lapply(1:12,function(i) rbind(data.frame(ll=ZIP.BIC[,i],dist="ZIP",metric="BIC"),
                                           data.frame(ll=ZIG.BIC[,i],dist="ZIG",metric="BIC"),
                                           data.frame(ll=ZINB.BIC[,i],dist="ZINB",metric="BIC"),
                                           data.frame(ll=ZIGP.BIC[,i],dist="ZIGP",metric="BIC"),
                                           data.frame(ll=ZICMP.BIC[,i],dist="ZICMP",metric="BIC"),
                                           data.frame(ll=ZIsCMP2.BIC[,i],dist="ZIsCMP2",metric="BIC"),
                                           data.frame(ll=ZIsCMP3.BIC[,i],dist="ZIsCMP3",metric="BIC")))

data.sims <- c(paste("Data Generation: ZIP (n = ",n,")",sep=""),
               paste("Data Generation: ZIB (n = ",n,")",sep=""),
               paste("Data Generation: ZINB (n = ",n,")",sep=""),
               paste("Data Generation: ZIG (n = ",n,")",sep=""),
               paste("Data Generation: ZIGP (Underdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZIGP (Overdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZICMP (Underdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZICMP (Overdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZIsCMP2 (Underdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZIsCMP2 (Overdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZIsCMP3 (Underdispersed) (n = ",n,")",sep=""),
               paste("Data Generation: ZIsCMP3 (Overdispersed) (n = ",n,")",sep=""))

ggplot(aes(y = ll, x = dist), data = all.ll.df[[i]]) + geom_boxplot() + labs(x = "Estimated Distribution",
                                                                             y = "Loglikelihood",
                                                                             title = data.sims[i])


ggplot(aes(y = ll, x = dist), data = all.aic.df[[i]]) + geom_boxplot() + labs(x = "Estimated Distribution",
                                                                             y = "AIC",
                                                                             title = data.sims[i])


ggplot(aes(y = ll, x = dist), data = all.bic.df[[i]]) + geom_boxplot() + labs(x = "Estimated Distribution",
                                                                             y = "BIC",
                                                                             title = data.sims[i])




