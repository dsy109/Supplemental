###############################################################
###############################################################
### Code for "Appendix C: Sample Size Determination Discussion"
### in "A Review of Statistical Reference Regions in Laboratory 
### Medicine: Theory and Computation" by Thomas Mathew and
### Derek S. Young
###############################################################
###############################################################

library(ggplot2)

k.n <- function(n,alpha) (1/n)*(1+qnorm(1-alpha)^2*(n-1-(2*exp(2*lgamma(n/2)-2*lgamma((n-1)/2)))))

n <- 10:500
plot(n,sqrt(k.n(n,.10)),type="l",ylab="k(n)",ylim=c(0.06,0.5))
lines(n,sqrt(k.n(n,.05)),col=2)
lines(n,sqrt(k.n(n,.01)),col=3)
abline(v=120,lty="dashed",col="gray")
legend("topright",col=1:3,text.col=1:3,lty=1,legend=c("90th","95th","99th"),title="Quantile",bty="n")

DF.q <- data.frame(n=rep(n,3),k.n=c(sqrt(k.n(n,.10)),sqrt(k.n(n,.05)),sqrt(k.n(n,.01))),Quantile=as.factor(rep(c("90th","95th","99th"),each=length(n))))
ggplot(DF.q, aes(x = n, y = k.n, group=Quantile))+ geom_vline(xintercept=120,color="black",lty="dashed") +geom_line(aes(color=Quantile,linetype=Quantile),lwd=1.3) + 
  ggtitle("Sample Size Effect on Standard Deviation of \nPoint Estimate of Percentile") + 
  theme(text = element_text(size = 15)) + ylab("k(n)") + annotate("text", x=160, y=0.5, label= "n = 120", size=5)

