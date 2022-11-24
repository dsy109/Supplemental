library(mixtools)
library(tolerance)
library(tidyverse)

set.seed(1)

n1 <- cbind(X1=rnorm(200),X2=rnorm(200))
d1 <- mixtools::depth(n1,n1)
n1 <- data.frame(n1)
t1 <- as.numeric(nptol.int(d1,alpha=0.10,P=0.90)[3])
r1 <- n1[d1>=t1,]
h1 <- r1 %>% slice(chull(X1, X2))

ggplot(n1,aes(x=X1,y=X2)) + geom_point() + xlim(-4,4) + ylim(-4,4) +
  geom_polygon(data = h1,alpha = 0.3,col="#00A5FF", fill="#00A5FF") +
  xlab(expression(X[1])) + ylab(expression(X[2])) 



n2 <- cbind(X1=rnorm(200),X2=rnorm(200))
d2 <- mixtools::depth(n2,n2)
n2 <- data.frame(n2)
t2 <- as.numeric(nptol.int(d2,alpha=0.10,P=0.90)[3])
r2 <- n2[d2>=t2,]
h2 <- r2 %>% slice(chull(X1, X2))

ggplot(n2,aes(x=X1,y=X2)) + geom_point() + xlim(-4,4) + ylim(-4,4) +
  geom_polygon(data = h2,alpha = 0.3,col="#00A5FF", fill="#00A5FF") +
  xlab(expression(X[1])) + ylab(expression(X[2])) 



n3 <- cbind(X1=rnorm(200),X2=rnorm(200))
d3 <- mixtools::depth(n3,n3)
n3 <- data.frame(n3)
t3 <- as.numeric(nptol.int(d3,alpha=0.10,P=0.90)[3])
r3 <- n3[d3>=t3,]
h3 <- r3 %>% slice(chull(X1, X2))

ggplot(n3,aes(x=X1,y=X2)) + geom_point() + xlim(-4,4) + ylim(-4,4) +
  geom_polygon(data = h3,alpha = 0.3,col="#00A5FF", fill="#00A5FF") +
  xlab(expression(X[1])) + ylab(expression(X[2])) 



n4 <- cbind(X1=rnorm(200),X2=rnorm(200))
d4 <- mixtools::depth(n4,n4)
n4 <- data.frame(n4)
t4 <- as.numeric(nptol.int(d4,alpha=0.10,P=0.90)[3])
r4 <- n4[d4>=t4,]
h4 <- r4 %>% slice(chull(X1, X2))

ggplot(n4,aes(x=X1,y=X2)) + geom_point() + xlim(-4,4) + ylim(-4,4) +
  geom_polygon(data = h4,alpha = 0.3,col="#00A5FF", fill="#00A5FF") +
  xlab(expression(X[1])) + ylab(expression(X[2])) 