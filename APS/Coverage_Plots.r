###############################################################
###############################################################
# Plots of coverage results (Section 4, Figure 2)
###############################################################
###############################################################

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir); getwd(); 
system("ls")

# install.packages("ggplot2")
library(ggplot2)
library(data.table) # to use melt() for data transpose wide to long
library(car)

## For publication

# Chap 4 Figure (a)

mydata_ccd2 = read.csv("ccd2_long.csv")
mydata_ccd2$P_num     = mydata_ccd2$P
mydata_ccd2$P     = as.factor(mydata_ccd2$P)
mydata_ccd2$gamma_num = mydata_ccd2$gamma
mydata_ccd2$gamma = as.factor(mydata_ccd2$gamma)

long_ccd2 = melt(mydata_ccd2, measure.vars = c("APS", "CH"))
setnames(long_ccd2, "variable", "method")
setnames(long_ccd2, "value", "coverage")

long80_ccd2 = subset(long_ccd2,gamma_num==0.8)
p_fig_4a = ggplot(long80_ccd2, aes(x=P,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("design:",unique(long80_ccd2$design) ,"nominal coverage: 0.8"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.8) +
  scale_y_continuous(breaks=c(0.8,0.85, 0.9)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4a

ggsave("Fig_4a_CCD2_gamma80.png",
       p_fig_4a, width = 5, height = 5, dpi=300, dev='png')

# Chap 4 Figure (b)

long95_ccd2 = subset(long_ccd2,gamma_num==0.95)
p_fig_4b = ggplot(long95_ccd2, aes(x=P,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("design:",unique(long95_ccd2$design) ,"nominal coverage: 0.95"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.95) +
  scale_y_continuous(breaks=c(0.95,0.96,0.97,0.98)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4b

ggsave("Fig_4b_CCD2_gamma95.png",
       p_fig_4b, width = 5, height = 5, dpi=300, dev='png')

# Chap 4 Figure (c)

mydata_all = read.csv("all_long.csv")

library(dplyr)
mydata_ccd_90_90 = mydata_all %>% filter(design == "CCD",P==0.9,gamma==0.9)

mydata_ccd_90_90$P_num     = mydata_ccd_90_90$P
mydata_ccd_90_90$P         = as.factor(mydata_ccd_90_90$P)
mydata_ccd_90_90$gamma_num = mydata_ccd_90_90$gamma
mydata_ccd_90_90$gamma     = as.factor(mydata_ccd_90_90$gamma)

p_fig_4c = ggplot(mydata_ccd_90_90, 
                  aes(x=m,y=coverage, group=interaction(model, method)),color=method ) + 
  ggtitle(paste("design:",unique(mydata_ccd_90_90$design) ,"P: 0.9 nominal coverage: 0.9"))+
  geom_line(aes(color=model,linetype = method)) +
  geom_point(aes(color=model,shape = method)) +
  geom_hline(yintercept = 0.9) +
  scale_y_continuous(breaks=c(0.9,0.925,0.95,0.975)) +
  scale_x_continuous(breaks=c(2,3,4)) +
  ylab("coverage") +
  theme(axis.title.y = element_text(size=14))
p_fig_4c

ggsave("Fig_4c_CCD_90_90.png",
       p_fig_4c, width = 5, height = 5, dpi=300, dev='png')

# Chap 4 Figure (d)

mydata_all = read.csv("all_long.csv")

library(dplyr)
mydata_m3_90_90 = mydata_all %>% filter(m == 3,P==0.9,gamma==0.9)

mydata_m3_90_90$P_num     = mydata_m3_90_90$P
mydata_m3_90_90$P         = as.factor(mydata_m3_90_90$P)
mydata_m3_90_90$gamma_num = mydata_m3_90_90$gamma
mydata_m3_90_90$gamma     = as.factor(mydata_m3_90_90$gamma)

p_fig_4d = ggplot(m3_9090, 
        aes(x=design,y=coverage, group=interaction(method, model)),
        color=model ) + 
   ggtitle("m:3 P:0.9 nominal coverage:0.9") +
   geom_line(aes(color=model,linetype = method)) +
   geom_point(aes(color=model,shape = method)) +
   geom_hline(yintercept = 0.9)
p_fig_4d

ggsave("Fig_4d_m3_90_90.png",
       p_fig_4d, width = 5, height = 5, dpi=300, dev='png')

