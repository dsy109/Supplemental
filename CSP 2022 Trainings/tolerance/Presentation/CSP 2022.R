library(tolerance)
library(plotly)
library(foreign)
library(AER)
library(devtools)

source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/regtol.int2.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/nlregtol.int2.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/npregtol.int2.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly.norm.OC.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly.plottol.anova.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly.plottol.control.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly.plottol.hist.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly.plottol.multi.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly.plottol.reg.R")


########################################
###### Normal Tolerance Intervals ######
########################################
#########################
### Example: K.factor ###
#########################
K.factor(n = 10, P = 0.95, side = 2, method = "HE")
K.factor(n = 10, P = 0.95, side = 2, method = "WBE")
K.factor(n = 10, P = 0.95, side = 2, method = "EXACT", m = 50)
K.factor(n = 100, P = 0.95, side = 2, method = "HE")
K.factor(n = 100, P = 0.95, side = 2, method = "WBE")
K.factor(n = 100, P = 0.95, side = 2, method = "EXACT", m = 50)

########################
### Example: K.table ###
########################
# Order by Sample Size #
K.table(n = seq(50, 60, 10), alpha = c(0.01, 0.05, 0.10), 
        P = c(0.90, 0.95, 0.99), by.arg = "n")
# Order by Content Level #
K.table(n = seq(50, 60, 10), alpha = c(0.01, 0.05, 0.10), 
        P = c(0.90, 0.95), by.arg = "P")
# Order by Confidence Level #
K.table(n = seq(50, 60, 10), alpha = c(0.01, 0.05), 
        P = c(0.90, 0.95, 0.99), by.arg = "alpha")

#############################
### Example: K.factor.sim ###
#############################
n_sizes <- c(2:6, seq(30, 70, 10))
l_sizes <- 2:6
KM_table <- sapply(1:length(l_sizes), function(i)
  sapply(1:length(n_sizes), 
         function(j) round(K.factor.sim(n = n_sizes[j],
                                        l = l_sizes[i], side=1, alpha = 0.1,
                                        P = 0.9),3)))
dimnames(KM_table) <- list(n = n_sizes, l = l_sizes)
KM_table

############################
### Example: normtol.int ###
############################
set.seed(100)
x <- rnorm(100, 0, 0.2)
# One-Sided #
out1 <- normtol.int(x = x, alpha = 0.05, P = 0.95, side = 1,
                    method = "HE", log.norm = FALSE)
out1
plotly.plottol.hist(out1, x, side = "two",
                    x.lab = "Normal Data (One-Sided)",
                    tol.line.type = "dash")
# Two-Sided #
out2 <- normtol.int(x = x, alpha = 0.05, P = 0.95, side = 2,
                    method = "HE", log.norm = FALSE)
out2
plotly.plottol.hist(out2, x, side = "two",
                    x.lab = "Normal Data (One-Sided)",
                    tol.line.type = "dash")

########################
### Example: norm.ss ###
########################
set.seed(100)
norm.ss(alpha = 0.05, P = 0.95, side = 2, spec = c(-3, 3),
        method = "DIR", hyper.par = list(mu.0 = 0, sig2.0 = 1))

########################
### Example: norm.OC ###
########################
# By Content Level #
norm.OC(k = 4, alpha = NULL, P = c(0.90, 0.95, 0.99), 
        n = 10:20, side = 1)

plotly.norm.OC(k = 4, alpha = NULL, P = c(0.90, 0.95, 0.99),
               n = 10:20, side = 1,
               x.cex = 8, line.width = 4,
               y.lab.size = 16, x.lab.size = 16,
               x.tick.size = 16, y.tick.size = 16,
               title.size = 16, legend.size = 12)

# By Confidence Level #
norm.OC(k = 4, alpha = c(0.01, 0.05, 0.10), P = NULL, n = 10:20,
        side = 1)

plotly.norm.OC(k = 4, alpha = c(0.01, 0.05, 0.10), P = NULL,
               n = 10:20, side = 1,
               x.cex = 8, line.width = 4,
               y.lab.size = 16, x.lab.size = 16,
               x.tick.size = 16, y.tick.size = 16,
               title.size = 16, legend.size = 12)

# By Sample Size #
norm.OC(k = NULL, P = c(0.90, 0.95, 0.99),
        alpha=c(0.01,0.05,0.10), n = 10:20, side = 1)

plotly.norm.OC(k = NULL, P = c(0.90, 0.95, 0.99),
               alpha=c(0.01,0.05,0.10), n = 10:20, side = 1,
               x.cex = 8, line.width = 4,
               y.lab.size = 16, x.lab.size = 16,
               x.tick.size = 16, y.tick.size = 16,
               title.size = 16, legend.size = 12)

##################################
### Application: Silicon Wafer ###
##################################
# Reference: D. C. Montgomery (2009). "Introduction to Statistical Quality Control, 6th Edition". Wiley. Hoboken, NJ.
wafer <- c(unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/wafer.txt")))
# Histogram #
fit <- density(wafer)

decimal.fun <- function(x) {sprintf("%.3f", x)}
plot_ly() %>% 
  add_trace(x = wafer, type = "histogram", 
            marker = list(color = "#1f77b4",
                          line = list(color = "#FFFFFF",
                                      width = 1)),
            name = "Thickness" , showlegend = FALSE) %>%
  add_trace(x = fit$x, y = fit$y, type = "scatter", mode = "lines", 
            line = list(width = 8 , dash="dash" , color = "#1f77b4"),
            fill = NULL, yaxis = "y2", name = "Density" , showlegend = FALSE) %>% 
  layout(title = list(text = "Distribution of Thickness of Metal Layer",
                      x = 0.5,
                      y = 0.98,
                      font = list(size=36)),
         xaxis = list(title = "Thickness",
                      range = c(400,500),
                      tickfont = list(size = 36),
                      tickvals = seq(from=400,to=500,by=10), 
                      ticktext = seq(from=400,to=500,by=10),
                      titlefont = list(size = 36),
                      zeroline = FALSE),
         yaxis = list(title = "Counts",
                      tickfont = list(size = 36),
                      titlefont = list(size = 36)),
         bargap = 0,
         yaxis2 = list(overlaying = "y", side = "right",
                       title = "Density",
                       position=0.95,
                       range = c(0,0.03),
                       tickvals = seq(from=0,to=0.03,by=0.005), 
                       ticktext = decimal.fun(seq(from=0,to=0.03,by=0.005)),
                       tickfont = list(size = 36),
                       titlefont = list(size = 36))
  )

# 2-Sided Tolerance Limits #
wafer.tol <- normtol.int(x=wafer, alpha=0.05, P=0.95, side=2, method="EXACT") #(0.95,0.95) 2-sided normal TI

plotly.plottol.hist(wafer.tol, wafer, side = "two",
                    x.lab = "Thickness of Metal Layer",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)

############################################
###### Non-normal Tolerance Intervals ######
############################################
##################################################
### Example: Binomial Distribution: bintol.int ###
##################################################
## 85%/90% 2-sided binomial tolerance intervals for a future
## lot of 2500 when a sample of 230 were drawn from a lot of
## 1000. 
bintol.int(x = 230, n = 1000, m = 2500, alpha = 0.15, P = 0.90,
           side = 2, method = "LS")
bintol.int(x = 230, n = 1000, m = 2500, alpha = 0.15, P = 0.90,
           side = 2, method = "AS")

## Using Jeffreys' method to construct the 85%/90% 1-sided
## binomial tolerance limits. The first calculation assumes
## a prior on the proportion of defects which places greater
## density on values near 0. The second calculation assumes
## a prior on the proportion of defects which places greater
## density on values near 1.
bintol.int(x = 230, n = 1000, m = 2500, alpha = 0.15, P = 0.90, 
           side = 1, method = "JF", a1 = 2, a2 = 10)
bintol.int(x = 230, n = 1000, m = 2500, alpha = 0.15, P = 0.90, 
           side = 1, method = "JF", a1 = 5, a2 = 1)

##################################################
### Example: Poisson Distribution: poistol.int ###
##################################################
## 95%/90% 1-sided Poisson tolerance limits for future
## occurrences in a period of length 3.
poistol.int(x = 45, n = 9, m = 3, alpha = 0.05, P = 0.90,
            side = 1, method = "TAB")
poistol.int(x = 45, n = 9, m = 3, alpha = 0.05, P = 0.90,
            side = 1, method = "LS")

## 95%/90% 2-sided Poisson tolerance intervals for future
## occurrences in a period of length 15.
poistol.int(x = 45, n = 9, m = 15, alpha = 0.05, P = 0.90,
            side = 2, method = "TAB")
poistol.int(x = 45, n = 9, m = 15, alpha = 0.05, P = 0.90,
            side = 2, method = "LS")

#################################################
### Example: Weibull Distribution: exttol.int ###
#################################################
## 85%/90% 1-sided Weibull tolerance intervals for a sample
## of size 150.
set.seed(100)
WeibullData <- rweibull(150, 3, 75)
out <- exttol.int(x = WeibullData, alpha = 0.15, P = 0.90, side = 1,
                  dist = "Weibull")
out

plotly.plottol.hist(out, WeibullData , side = "two",
                    x.lab = "Random Data from Weibull Distribution",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)

## 85%/90% 2-sided Gumbel distribution for the maximum tolerance intervals 
## for a sample of size 200.
# install.packages("evd")
library(evd)
set.seed(100)
GumbelData <- rgumbel(200, loc=3, scale=5)
out <- exttol.int(x = GumbelData, alpha = 0.15, P = 0.90, side = 2,
                  dist = "Gumbel" , ext = "max")
out

plotly.plottol.hist(out, GumbelData , side = "two",
                    x.lab = "Random Data from Gumbel Distribution",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)

###############################################
### Example: Gamma Distribution: gamtol.int ###
###############################################
set.seed(100)
x <- rgamma(50, 0.30, scale = 2)
## 99%/99% 1-sided gamma tolerance intervals for a sample
## of size 50, using Howe method.
out1 <- gamtol.int(x = x, alpha = 0.01, P = 0.99, side = 1,
                   method = "HE")
out1
plotly.plottol.hist(out1, x , side = "two",
                    x.lab = "Random Data from Gamma Distribution",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)
## 95%/95% 2-sided gamma tolerance intervals for a sample
## of size 50, using EXACT method.
out2 <- gamtol.int(x = x, alpha = 0.05, P = 0.95, side = 2,
                   method = "EXACT")
out2
plotly.plottol.hist(out2, x , side = "two",
                    x.lab = "Random Data from Gamma Distribution",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)

##########################################
### Application: Binomial Distribution ###
##########################################
# Defective Chips
# Reference: NIST (2022). "Proportions Control Chart". https://www.itl.nist.gov/div898/handbook/pmc/section3/pmc332.htm
defects <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/defects.txt",header=TRUE)$Defects
binom.app<- bintol.int(x=defects, n=50*30, m=50, alpha=0.05, P=0.99, side=2, method="CP") #(0.95,0.99) 2-sided binomial TI
binom.app
plotly.plottol.control(tol.out = binom.app , x=defects,
                       x.lab = "Wafer",
                       y.lab = "Defectives",
                       x.cex = 26,
                       fit.lwd = 8,
                       tol.line.type = "dash",
                       tol.lwd = 8,
                       x.lab.size = 36, x.tick.size = 36,
                       y.lab.size = 36, y.tick.size = 36,
                       title.size = 36, title.position.y = 0.98)

#########################################
### Application: Poisson Distribution ###
#########################################
# Hospital Visits
# Reference: A. C. Cameron and P. K. Trivedi (1986). "Econometric Models Based on Count Data: Comparisons and Applications of Some Estimators and Tests." Journal of Applied Econometrics. 1(1):29-53
library(AER)
data("DoctorVisits")
visits <- DoctorVisits$visits
pois.app <- poistol.int(x=visits, n=length(visits), m=10, alpha=0.01, P=0.99, side=2, method="TAB") #(0.99,0.99) 2-sided Poisson TI for a future sample of 10 individuals
pois.app
plotly.plottol.hist(pois.app, visits , side = "two",
                    x.lab = "Doctor Visits",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)


#########################################
### Application: Weibull Distribution ###
#########################################
# Kidney Function Data
# Reference: E. Harris and J. C. Boyd (1995). "Statistical Bases of Reference Values in Laboratory Medicine." Marcel-Dekker. New York, NY.
# The variable is serum creatinine.
kidney <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/kidney.txt",header=TRUE)
weibull.app <- exttol.int(x=kidney$SCR, alpha=0.10, P=0.90, side=1, dist="Weibull") #(0.90,0.90) 1-sided upper Weibull tolerance limit
weibull.app
plotly.plottol.hist(weibull.app, kidney$SCR , side = "two",
                    x.lab = "Serum Creatinine",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)

#######################################
### Application: Gamma Distribution ###
#######################################
# Bladder Cancer
# Reference: E. T. Lee and J. Wang (2003). "Statistical Methods for Survival Data Analysis." Wiley. Hoboken, NJ.
bladder <- unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/bladder.txt",header=FALSE))
gamma.app <- gamtol.int(x=bladder, alpha=0.05, P=0.95, side=1) #(0.95,0.95) 1-sided upper gamma tolerance limit
gamma.app
plotly.plottol.hist(gamma.app, bladder , side = "upper",
                    x.lab = "Bladder Cancer",
                    tol.line.type = "dash",
                    tol.lwd = 8,
                    x.lab.size = 36, x.tick.size = 36,
                    y.lab.size = 36, y.tick.size = 36,
                    title.size = 36, title.position.y = 0.98)

##############################################
###### Nonparametric Tolerance Interval ######
##############################################
############################################
### Example: Nonparametric TI: notol.int ###
############################################
set.seed(100)
x <- rlogis(100, 5, 1)
# One-Sided TI #
out1 <- nptol.int(x = x, alpha = 0.10, P = 0.90, side = 1,
                  method = "WILKS", upper = NULL, lower = NULL)
out1
plotly.plottol.hist(out1 , x=x , side = "two" ,
                    x.lab = "X", tol.lwd = 18 , tol.line.type = "dash",
                    x.lab.size = 36, x.tick.size=36 ,
                    y.lab.size = 36 , y.tick.size = 36 ,
                    title.size = 36 , title.position.y = 0.995)
# Two-Sided TI #
out2 <- nptol.int(x = x, alpha = 0.10, P = 0.90, side = 2,
                  method = "WILKS", upper = NULL, lower = NULL)
out2
plotly.plottol.hist(out2 , x=x , side = "two" ,
                    x.lab = "X", tol.lwd = 18 , tol.line.type = "dash",
                    x.lab.size = 36, x.tick.size=36 ,
                    y.lab.size = 36 , y.tick.size = 36 ,
                    title.size = 36 , title.position.y = 0.995)

##############################################
### Example: Beta-Expected TI: npbetol.int ###
##############################################
set.seed(100)
x <- rlogis(100, 5, 1)
# One-Sided #
out1 <- npbetol.int(x = x, Beta = 0.90, side = 1,
                    upper = NULL, lower = NULL)
out1
plotly.plottol.hist(out1 , x=x , side = "two" ,
                    x.lab = "X", tol.lwd = 18 , tol.line.type = "dash",
                    x.lab.size = 36, x.tick.size=36 ,
                    y.lab.size = 36 , y.tick.size = 36 ,
                    title.size = 36 , title.position.y = 0.995)
# Two-Sided #
out2 <- npbetol.int(x = x, Beta = 0.90, side = 2,
                    upper = NULL, lower = 4)
out2
plotly.plottol.hist(out2 , x=x , side = "two" ,
                    x.lab = "X", tol.lwd = 18 , tol.line.type = "dash",
                    x.lab.size = 36, x.tick.size=36 ,
                    y.lab.size = 36 , y.tick.size = 36 ,
                    title.size = 36 , title.position.y = 0.995)

####################################################
### Example: Sample Size Determination: np.order ###
####################################################
## Only requesting the sample size.
np.order(m = 5, alpha = 0.05, P = 0.95, indices = FALSE)

## Requesting the order statistics indices as well.
np.order(m = 5, alpha = 0.05, P = 0.95, indices = TRUE)

#####################################################
### Application: Nonparametric Tolerance Interval ###
#####################################################
# Breast Cancer Remission
# Reference: M. W. A. Ramos, G. M. Cordeiro, P. R. D. Marinho, C. R. B. Dias, and G. G. Hamedani (2013). "The Zografos-Balakrishnan Log-Logistic Distribution: Properties and Applications." Journal of Statistical Theory and Applications. 13(1):65-82.
remission <- unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/remission.txt",header=FALSE))
# (0.90,0.90) 2-sided nonparametric tolerance interval (Using Wilks's Method) #
remission.WILKS <- nptol.int(x=remission, alpha=0.10, P=0.90, side=2, method="WILKS") 
remission.WILKS
plotly.plottol.hist(remission.WILKS , x=remission , side = "two" ,
                    x.lab = "Survival Time", tol.lwd = 18 , tol.line.type = "dash",
                    x.lab.size = 36, x.tick.size=36 ,
                    y.lab.size = 36 , y.tick.size = 36 ,
                    title.size = 36 , title.position.y = 0.995)
# (0.90,0.90) 2-sided nonparametric tolerance interval (Using the interpolation/extrapolation scheme of Young and Mathew (2014))
remission.YM <- nptol.int(x=remission, alpha=0.10, P=0.90, side=2, method="YM") 
remission.YM
plotly.plottol.hist(remission.YM , x=remission , side = "two" ,
                    x.lab = "Survival Time", tol.lwd = 18 , tol.line.type = "dash",
                    x.lab.size = 36, x.tick.size=36 ,
                    y.lab.size = 36 , y.tick.size = 36 ,
                    title.size = 36 , title.position.y = 0.995)

#########################################
###### Regression Tolerance Limits ######
#########################################
#################################################
### Example: Linear Regression TI: regtol.int ###
#################################################
set.seed(100)
x <- runif(100, 0, 10)
y <- 20 + 5*x + rnorm(100, 0, 3)
# One-Sided #
out1 <- regtol.int2(reg = lm(y ~ x),
                    side = 1, alpha = 0.05, P = 0.95 , new=TRUE)
out1

plotly.plottol.reg(out1 , x=x , y=y , side = "two",
                   x.cex = 22 , tol.lwd = 18, fit.lwd = 18,
                   fit.line.type = "dash", tol.line.type = "dot",
                   x.tick.size = 36 , x.lab.size = 36,
                   y.tick.size = 36 , y.lab.size = 36,
                   title.size = 36 , title.position.y = 0.995)
# Two-Sided #
out2 <- regtol.int2(reg = lm(y ~ x), new.x = data.frame(x = c(3, 6, 12)),
                    side = 2, alpha = 0.05, P = 0.95 , new=TRUE)
out2

plotly.plottol.reg(out2 , x=x , y=y , side = "two",
                   x.cex = 22 , tol.lwd = 18, fit.lwd = 18,
                   fit.line.type = "dash", tol.line.type = "dot",
                   x.tick.size = 36 , x.lab.size = 36,
                   y.tick.size = 36 , y.lab.size = 36,
                   title.size = 36 , title.position.y = 0.995)

######################################################
### Example: Nonlinear Regression TI: nlregtol.int ###
######################################################
set.seed(100)
x <- runif(50, 5, 45)
f1 <- function(x, b1, b2) b1 + (0.49 - b1)*exp(-b2*(x - 8)) +
  rnorm(50, sd = 0.01)
y <- f1(x, 0.39, 0.11)
formula <- as.formula(y ~ b1 + (0.49 - b1)*exp(-b2*(x - 8)))

# One-Sided #
out1 <- nlregtol.int2(formula = formula,
                      xy.data = data.frame(cbind(y, x)), side = 1,
                      alpha = 0.05, P = 0.95 , new = TRUE)
out1

plotly.plottol.reg(out1 , x=x , y=y , side = "two",
                   x.cex = 22 , tol.lwd = 18, fit.lwd = 18,
                   fit.line.type = "dash", tol.line.type = "dot",
                   x.tick.size = 36 , x.lab.size = 36,
                   y.tick.size = 36 , y.lab.size = 36,
                   title.size = 36 , title.position.y = 0.995)
# Two-Sided #
out2 <- nlregtol.int2(formula = formula,
                      xy.data = data.frame(cbind(y, x)),
                      x.new=cbind(c(10, 55)), side = 2,
                      alpha = 0.05, P = 0.95 , new = TRUE)
out2
plotly.plottol.reg(out2 , x=x , y=y , side = "two",
                   x.cex = 22 , tol.lwd = 18, fit.lwd = 18,
                   fit.line.type = "dash", tol.line.type = "dot",
                   x.tick.size = 36 , x.lab.size = 36,
                   y.tick.size = 36 , y.lab.size = 36,
                   title.size = 36 , title.position.y = 0.995)

#######################################################
### Example: Nonparametric Regression TI: nptol.int ###
#######################################################
set.seed(100)
x <- runif(50, 5, 45)
f1 <- function(x, b1, b2) b1+(0.49-b1)*exp(-b2*(x-8))+rnorm(50,sd = 0.01)
y <- f1(x, 0.39, 0.11)
y.hat <- loess(y~x)$fit
# One-Side #
out1 <- npregtol.int2(x = x, y = y, y.hat = y.hat, side = 1,
                      alpha = 0.05, P = 0.95, method = "WILKS",new=TRUE)
out1

plotly.plottol.reg(out1 , x=x , y=y , side = "two",
                   x.cex = 22 , tol.lwd = 18, fit.lwd = 18,
                   fit.line.type = "dash", tol.line.type = "dot",
                   x.tick.size = 36 , x.lab.size = 36,
                   y.tick.size = 36 , y.lab.size = 36,
                   title.size = 36 , title.position.y = 0.995)
# Two-Sided #
out2 <- npregtol.int2(x = x, y = y, y.hat = y.hat, side = 2,
                      alpha = 0.05, P = 0.95, method = "WILKS",new=TRUE,
                      lower = 0.38 , upper=0.49)
out2

plotly.plottol.reg(out2 , x=x , y=y , side = "two",
                   x.cex = 22 , tol.lwd = 18, fit.lwd = 18,
                   fit.line.type = "dash", tol.line.type = "dot",
                   x.tick.size = 36 , x.lab.size = 36,
                   y.tick.size = 36 , y.lab.size = 36,
                   title.size = 36 , title.position.y = 0.995)

###################################################
### Application: Linear Regression Data Example ###
###################################################
# Hospital Infections Data
# Reference: Neter, J., Kutner, M. H., Nachtsheim, C. J., Wasserman, W., et al. (1996). Applied
# linear statistical models.
library(devtools)
hospitals <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/hospitals.txt", 
                        header = TRUE)
### Calculate 2-sided pointwise (0.90, 0.90) tolerance intervals for the multiple
### regression fit of the hospital infections data.
out <- lm(InfctRsk ~ Xray + Stay, data = hospitals)
out.TI <- regtol.int2(out, new.x = data.frame(Xray = c(50, 75), Stay = c(8, 12)), 
                      alpha = 0.1, P = 0.9, side = 2, new = TRUE)
plotly.plottol.reg(out.TI, x = hospitals[c("Xray", "Stay")], y = hospitals["InfctRsk"], 
                   side = "two", rect = TRUE, x.cex = 12,
                   x.lab = "X Ray", x.lab.size = 36, x.tick.size = 18,
                   y.lab = "Stay", y.lab.size = 36, y.tick.size = 18,
                   z.lab = "Infection Risk", z.lab.size = 36, z.tick.size = 18,
                   title.size = 36 , title.position.y = 0.995)

###########################################
###### Multivariate Tolerance Region ######
###########################################
###############################
### Example: Bivariate Data ###
###############################
## 90%/90% bivariate normal tolerance region.
set.seed(100)
x1 <- rnorm(100, 0, 0.2)
x2 <- rnorm(100, 0, 0.5)
x <- cbind(x1, x2)
out1 <- mvtol.region(x = x, alpha = 0.10, P = 0.90, B = 1000,
                     method = "KM")
out1
plotly.plottol.multi(out1, x , x.lab = "X1" , y.lab = "X2",
                     x.cex = 22 , tol.lwd = 18, tol.line.type = "dot",
                     x.tick.size = 36 , x.lab.size = 36,
                     y.tick.size = 36 , y.lab.size = 36,
                     title.size = 36 , title.position.y = 0.995)

################################
### Example: Trivariate Data ###
################################
## 90%/90% trivariate normal tolerance region.
set.seed(100)
x1 <- rnorm(100, 0, 0.2)
x2 <- rnorm(100, 0, 0.5)
x3 <- rnorm(100, 5, 1)
x <- cbind(x1, x2, x3)
mvtol.region(x = x, alpha = c(0.10, 0.05, 0.01),
             P = c(0.90, 0.95, 0.99), B = 1000, method = "KM")
out2 <- mvtol.region(x = x, alpha = 0.10, P = 0.90, B = 1000,
                     method = "KM")
out2

plotly.plottol.multi(out2, x , x.cex = 22 , tol.lwd = 18,
                     x.lab = "X1" , y.lab = "X2" , z.lab = "X3",
                     title.size = 36 , title.position.y = 0.995)

#########################################################
### Example: Multivariate Regression Tolerance Region ###
#########################################################
## 95%/95% multivariate regression tolerance factors using
## a fertilizer data set presented in Anderson (2003, p. 374).
grain <- c(40, 17, 9, 15, 6, 12, 5, 9)
straw <- c(53, 19, 10, 29, 13, 27, 19, 30)
fert <- c(24, 11, 5, 12, 7, 14, 11, 18)
DF <- data.frame(grain,straw,fert)
new.x <- data.frame(fert = c(10, 15, 20))
mvreg <- lm(cbind(grain, straw) ~ fert + I(fert^2), data = DF)
set.seed(100)
out <- mvregtol.region(mvreg, new.x = new.x, alpha = 0.05,
                       P = 0.95, B = 5000)
out

########################################################
### Application: Multivariate Tolerance Regions Data ###
########################################################
# Adolescent Kidney Function Reference Regions
# Reference: D.S.Young and T. Mathew (2020), “Nonparametric Hyperrectangular Tolerance and
# Prediction Regions for Setting Multivariate Reference Regions in Laboratory Medicine.”
# Statistical Methods in Medical Research; 29(12):3569–3585.
library(mixtools)
ref.X <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/kidney_adolescents.txt", header = TRUE)
ref.males <- subset(ref.X, Gender == 1)
### A function to format the elliptical depth calculation, which is used in the
### function for constructing multivariate hyperrectangular tolerance regions.
Elliptical <- function(pts, x)
{
  out <- mixtools::depth(pts, x)
  out
}
### This will take a minute to run.
tr.male <- npmvtol.region(x = as.matrix(ref.males[, -1]), alpha = 0.05, P = 0.95, 
                          depth.fn = Elliptical, type = "semispace", adjust = "ceiling", semi.order = list(lower = NULL, 
                                                                                                           center = 2:3, upper = 1))
tr.male
### Sorting the results in order to produce pairwise scatterplots with reference
### regions overlaid.
male.df <- data.frame(C1 = (ref.males[, 2] <= tr.male[1, 2]), C2 = ((ref.males[, 
                                                                               3] <= tr.male[2, 2]) & (ref.males[, 3] >= tr.male[2, 1])), C3 = ((ref.males[, 
                                                                                                                                                           4] <= tr.male[3, 2]) & (ref.males[, 4] >= tr.male[3, 1])))
male.df <- data.frame(ref.males[, -1], C123 = as.factor(male.df$C1 & male.df$C2 & 
                                                          male.df$C3))

ggplot(data = male.df, aes(UACR, Uric.Acid, col = C123, shape = C123)) + geom_rect(aes(xmin = 0, 
                                                                                       xmax = tr.male[1, 2], ymin = tr.male[2, 1], ymax = tr.male[2, 2]), fill = "gray80", 
                                                                                   alpha = 0.1, color = "black", linetype = 2) + geom_point(size = 1.5) + theme(legend.position = "none") + 
  xlab("UACR (mg/g)") + ylab("UA (mg/dL)") + theme(text = element_text(size = 20)) + 
  ggtitle("UA vs. UACR Reference Regions (Males)") + xlim(0, 37.5) + ylim(1, 15)

ggplot(data = male.df, aes(UACR, Serum.Creatinine, col = C123, shape = C123)) + geom_rect(aes(xmin = 0, 
                                                                                              xmax = tr.male[1, 2], ymin = tr.male[3, 1], ymax = tr.male[3, 2]), fill = "gray80", 
                                                                                          alpha = 0.1, color = "black", linetype = 2) + geom_point(size = 1.5) + theme(legend.position = "none") + 
  xlab("UACR (mg/g)") + ylab("SC (mg/dL)") + theme(text = element_text(size = 20)) + 
  ggtitle("SC vs. UACR Reference Regions (Males)") + xlim(0, 37.5) + ylim(0, 5)

ggplot(data = male.df, aes(Uric.Acid, Serum.Creatinine, col = C123, shape = C123)) + 
  geom_rect(aes(xmin = tr.male[2, 1], xmax = tr.male[2, 2], ymin = tr.male[3, 1], 
                ymax = tr.male[3, 2]), fill = "gray80", alpha = 0.1, color = "black", linetype = 2) + 
  geom_point(size = 1.5) + theme(legend.position = "none") + xlab("UA (mg/dL)") + 
  ylab("SC (mg/dL)") + theme(text = element_text(size = 20)) + ggtitle("SC vs. UA Reference Regions (Males)") + 
  xlim(1.75, 15) + ylim(0, 5)
