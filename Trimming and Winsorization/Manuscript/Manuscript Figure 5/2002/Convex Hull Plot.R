library(ggplot2)
setwd("/Users/kcheng/Desktop/Real Data/2002")
load("Finalized.Data.Adjusted .RData")

data <- Finalized.Data.Adj
Farm.02 <- data$X2002Farms
Acre.02 <- data$X2002Acres
data.02 <- cbind(Farm.02,Acre.02)

N <- dim(data)[1]
tol.pos.data.unadj <- 287
######################################################################
load("data.with.depth.02 .RData")

depth.02 <- data.with.depth.02[,3]

tol.pos.data.adj <- 277
#######################################################################
tukey.pos.data <- 308
########################################################################

data.02.ordered <- data.02[order(depth.02 , decreasing = FALSE),]
data.02.ordered.df <- data.frame(data.02.ordered , 
                                 depth.02[order(depth.02 , decreasing = FALSE)])
names(data.02.ordered.df) <- c("Farm" , "Acre" , "Depth")

####################################################
find_hull <- function(df) {c(chull(df[,1], df[,2]),chull(df[,1], df[,2])[1])}
tol.hull.data.02.unadj <- find_hull(data.02.ordered.df[tol.pos.data.unadj:N,])
tol.hull.data.02.points.unadj <- matrix(c(data.02.ordered.df[tol.pos.data.unadj:N,]$Farm[tol.hull.data.02.unadj],
                                          data.02.ordered.df[tol.pos.data.unadj:N,]$Acre[tol.hull.data.02.unadj]),
                                          ncol = 2 , byrow = FALSE)
tol.hull.data.02.points.unadj.df <- data.frame(tol.hull.data.02.points.unadj)
#####################################################
tol.hull.data.02.adj <- find_hull(data.02.ordered.df[tol.pos.data.adj:N,])
tol.hull.data.02.points.adj <- matrix(c(data.02.ordered.df[tol.pos.data.adj:N,]$Farm[tol.hull.data.02.adj],
                                        data.02.ordered.df[tol.pos.data.adj:N,]$Acre[tol.hull.data.02.adj]),
                                        ncol = 2 , byrow = FALSE)
tol.hull.data.02.points.adj.df <- data.frame(tol.hull.data.02.points.adj)
######################################################
tukey.hull.data.02 <- find_hull(data.02.ordered.df[tukey.pos.data:N,])
tukey.hull.data.02.points <- matrix(c(data.02.ordered.df[tukey.pos.data:N,]$Farm[tukey.hull.data.02],
                                      data.02.ordered.df[tukey.pos.data:N,]$Acre[tukey.hull.data.02]),
                                      ncol = 2 , byrow = FALSE)
tukey.hull.data.02.points.df <- data.frame(tukey.hull.data.02.points)
######################################################
# Note: hull points, groups, and color in ggplot need to match order. #
convex.hull.data.02.comb <- rbind(tol.hull.data.02.points.adj , tol.hull.data.02.points.unadj  , tukey.hull.data.02.points)
group <- c(rep(1,dim(tol.hull.data.02.points.adj)[1]),
           rep(2,dim(tol.hull.data.02.points.unadj)[1]) , 
           rep(3,dim(tukey.hull.data.02.points)[1]))
convex.hull.data.02.comb.df <- data.frame(cbind(convex.hull.data.02.comb , group))

### Plot Original Convex Hull ###
ggplot(convex.hull.data.02.comb.df) +
  geom_point(data=data , x=data$X2002Farms , y=data$X2002Acres , cex=2 , color="black") +
  geom_path(aes(x=convex.hull.data.02.comb.df$V1 , y=convex.hull.data.02.comb.df$V2 , 
                group=group , col=factor(group)) , lwd=3) +
  xlab("Number of Farms (x1000)") + ylab("Acreage (x10^6)") +
  ggtitle("2002 Census Data")+
  scale_x_continuous(limits=c((min(data$X2002Farms)-1),(max(data$X2002Farms)+1))) +
  scale_y_continuous(limits=c((min(data$X2002Acres)-1),(max(data$X2002Acres)+1))) +
  scale_colour_manual(name="Adjustment",
                      values=c('1'="#9B111E" , '2'="#F4C430" , '3'="#008080"), 
                      labels=c("Tol-Adjusted" , "Tol-Unadjusted" , "Tukey")) +
  theme(legend.position=c(0.85,0.85),legend.direction="vertical",
        legend.text=element_text(size=25),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_x_continuous(limits=c(0,max(data$X2002Farms)), 
                     labels = sprintf(seq(from=0,to=max(data$X2002Farms/1000),by=0.5) , fmt = '%#.1f') , 
                     breaks = seq(from=0,to=max(data$X2002Farms),by=500)) +
  scale_y_continuous(limits=c(0,max(data$X2002Acres)), 
                     labels = sprintf(seq(from=0,to=max(data$X2002Acres/1000000),by=1) , fmt = '%#.1f') , 
                     breaks = seq(from=0,to=max(data$X2002Acres),by=1000000)) +
  theme(axis.text.x = element_text(face="bold", size=25 , angle=90) ,
        axis.title.x = element_text(face="bold", size=25, vjust=-0.75) ,
        axis.text.y = element_text(face="bold", size=25) ,
        axis.title.y = element_text(face="bold", size=25, vjust=1.5) ,
        plot.title=element_text(family='', face='bold', size=30, 
                                vjust=0.5,hjust=0.5,margin=margin(t=30,b=-40)))

#########################################
ggplot(convex.hull.data.02.comb.df) +
  geom_point(data=data , x=data$X2002Farms , y=data$X2002Acres , cex=2 , color="black") +
  geom_path(aes(x=convex.hull.data.02.comb.df$V1 , y=convex.hull.data.02.comb.df$V2 , 
                group=group , col=factor(group)) , lwd=3) +
  xlab("Number of Farms (x1000)") + ylab("Acreage (x10^6)") +
  ggtitle("2002 Census Data")+
  scale_x_continuous(limits=c((min(data$X2002Farms)-1),(2000))) +
  scale_y_continuous(limits=c((min(data$X2002Acres)-1),(1500000))) +
  scale_colour_manual(name="Adjustment",
                      values=c('1'="#9B111E" , '2'="#F4C430" , '3'="#008080"), 
                      labels=c("Tol-Adjusted" , "Tol-Unadjusted" , "Tukey")) +
  theme(legend.position=c(0.85,0.85),legend.direction="vertical",
        legend.text=element_text(size=25),
        legend.title=element_text(size=0),
        legend.background = element_rect(color = "black",fill = "grey90", size = 2, linetype = "solid")) +
  scale_x_continuous(limits=c(0,2000), 
                     labels = sprintf(seq(from=0,to=2,by=0.5) , fmt = '%#.1f') , 
                     breaks = seq(from=0,to=2000,by=500)) +
  scale_y_continuous(limits=c(0,1500000), 
                     labels = sprintf(seq(from=0,to=1.5,by=0.5) , fmt = '%#.1f') , 
                     breaks = seq(from=0,to=1500000,by=500000)) +
  theme(axis.text.x = element_text(face="bold", size=25 , angle=90) ,
        axis.title.x = element_text(face="bold", size=25, vjust=-0.75) ,
        axis.text.y = element_text(face="bold", size=25) ,
        axis.title.y = element_text(face="bold", size=25, vjust=1.5) ,
        plot.title=element_text(family='', face='bold', size=30, 
                                vjust=0.5,hjust=0.5,margin=margin(t=30,b=-40)))

