install.packages(c("mixtools", "tolerance", "dlstats", "ggplot2", "GGally", "plotly", 
    "devtools", "mixreg"))

###########################################################
### Package Downloads
###########################################################

### This just reproduces the package download figures that I presented near the
### beginning of the talk.

library(dlstats)
library(ggplot2)

x <- cran_stats(c("mixtools", "tolerance"))

ggplot(x[1:67, ], aes(end, downloads, group = package, color = package)) + theme(legend.position = "none") + 
    ggtitle("mixtools Monthly Downloads") + geom_line(lwd = 1.7) + geom_point(aes(shape = package), 
    size = 3.5) + xlab("Month") + ylab("Downloads") + theme(text = element_text(size = 25))

for (i in 0:4) print(round(mean(x[(i * 12 + 1):((i + 1) * 12), 3]), 0))

ggplot(x[68:134, ], aes(end, downloads, group = package)) + theme(legend.position = "none") + 
    ggtitle("tolerance Monthly Downloads") + geom_line(color = "#00BFC4", lwd = 1.7) + 
    geom_point(aes(shape = package), color = "#00BFC4", size = 3.5) + xlab("Month") + 
    ylab("Downloads") + theme(text = element_text(size = 25))

for (i in 0:4) print(round(mean(x[(i * 12 + 68):((i + 1) * 12 + 68), 3]), 0))

###########################################################
### Example #1: Quasars Data
###########################################################

library(mixtools)
library(ggplot2)

y <- read.table("http://astrostatistics.psu.edu/datasets/QSO_absorb.txt", nrows = 104, 
    skip = 1)[, 2]

### Construct histogram of the quasars data.

quasars <- data.frame(y)
m <- ggplot(quasars, aes(x = y)) + geom_histogram(aes(y = ..density..), bins = 9, 
    color = "black", fill = "gray") + theme(text = element_text(size = 20)) + xlab("Velocity (km/s)") + 
    ylab("Density") + ggtitle("Si IV 1394 Data")
m

### Create a function to (1) perform multiple random starting values for the
### mixture-of-normals EM algorithm, and (2) sort the observations according to
### their posterior membership probabilities to feed into the mixturegram.

set.seed(100)
selection <- function(i, data, rep = 30)
{
    out <- replicate(30, normalmixEM(data, epsilon = 1e-06, k = i, maxit = 5000), 
        simplify = F)
    counts <- lapply(1:30, function(j) table(apply(out[[j]]$posterior, 1, which.max)))
    counts.length <- sapply(counts, length)
    counts.min <- sapply(counts, min)
    counts.test <- (counts.length != i) | (counts.min < 5)
    if (sum(counts.test) > 0 & sum(counts.test) < rep) 
        out <- out[!counts.test]
    l <- unlist(lapply(out, function(x) x$loglik))
    tmp <- out[[which.max(l)]]
}

### Order the results for each set of estimators for the different mixture models,
### and then feed those into the mixturegram.

all.out <- lapply(2:6, selection, data = y)
pmbs <- lapply(1:length(all.out), function(i) all.out[[i]]$post)
mix.out <- mixturegram(y, pmbs, method = "pca", all.n = FALSE, id.con = NULL, main = "Mixturegram (Quasars Data)")
mix.out.all <- mixturegram(y, pmbs, method = "pca", all.n = TRUE, id.con = NULL, 
    main = "Mixturegram (Quasars Data)")

### Calculate model selection criteria.

all.ll <- lapply(1:length(all.out), function(i) all.out[[i]]$loglik)
all.AIC <- c(sum(dnorm(y, mean(y), sd(y), log = T)) - 2, sapply(1:length(all.out), 
    function(i) all.ll[[i]] - (length(unlist(all.out[[i]][2:4])) - 1)))
all.BIC <- c(sum(dnorm(y, mean(y), sd(y), log = T)) - log(300), sapply(1:length(all.out), 
    function(i) all.ll[[i]] - (length(unlist(all.out[[i]][2:4])) - 1) * log(300)/2))
all.ICL <- c(all.BIC[1], sapply(1:length(all.out), function(i) all.BIC[i + 1] - sum(all.out[[i]]$lambda * 
    log(all.out[[i]]$lambda))))
all.CAIC <- c(all.BIC[1] - 1, sapply(1:length(all.out), function(i) all.BIC[i + 1] - 
    (length(unlist(all.out[[i]][2:4])) - 1)/2))
round(rbind(all.AIC, all.BIC, all.ICL, all.CAIC), 3)

### Create histogram with 2-component mixture-of-normals density overlaid.

mix.fun <- function(x) all.out[[1]]$lambda[1] * dnorm(x, all.out[[1]]$mu[1], all.out[[3]]$sigma[1]) + 
    all.out[[1]]$lambda[2] * dnorm(x, all.out[[1]]$mu[2], all.out[[1]]$sigma[2])
m + stat_function(fun = mix.fun, color = "red")

###########################################################
### Example #2: Diffuse Large B-Cell Lymphoma Data
###########################################################

library(GGally)
library(plotly)

DLBCL <- read.table("https://stat.as.uky.edu/sites/default/files/DLBCL.txt", header = TRUE)

### Estimate a 2-component mixture of multivariate normals for the DLBCL data.

set.seed(240)
DLBCLSample <- DLBCL[sample(1:nrow(DLBCL), 500), ]
out <- mvnormalmixEM(DLBCLSample, k = 2, verb = TRUE)
str(out)
post <- as.factor(apply(out$posterior, 1, which.max))
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

###########################################################
### Example #3: Aphids Data
###########################################################

library(mixreg)
data(aphids)
attach(aphids)

aphids.df <- data.frame(x = n.aphids, y = n.inf)
ggplot(data = aphids.df, aes(x, y)) + geom_point() + labs(title = "Aphids Data", 
    x = "Number of Aphids in Batch", y = "Number of Infected Plants") + theme(text = element_text(size = 20))

### Fit a 2-component mixture of linear regressions to the aphids data.

set.seed(100)
out.em <- regmixEM(n.inf, n.aphids, k = 2)
mix.lab <- as.factor(apply(out.em$posterior, 1, which.max) + 1)
out.em[3:6]

### Scatterplot of the data with the estimate from the 2-component mixture of
### regressions overlaid.

aphids.df1 <- cbind(aphids.df, group = mix.lab)
ggplot(data = aphids.df1, aes(x, y, col = group, pch = group)) + geom_point() + labs(title = "Aphids Data (Parametric Mixture Fit)", 
    x = "Number of Aphids in Batch", y = "Number of Infected Plants") + geom_abline(slope = out.em$beta[, 
    1][2], intercept = out.em$beta[, 1][1], col = 2) + geom_abline(slope = out.em$beta[, 
    2][2], intercept = out.em$beta[, 2][1], col = 3) + scale_color_manual("", values = c(2, 
    3), guide = F) + theme(text = element_text(size = 20), legend.position = "none")

### Estimate a semiparametric mixtures of regressions assuming only a symmetric
### error density.  The bandwidth will adapt by default.

set.seed(100)
spfit1 <- spregmix(n.inf ~ n.aphids, bw = 0.1, data = aphids, verb = TRUE)

mix.lab2 <- as.factor(apply(spfit1$posterior, 1, which.max) + 1)
aphids.df2 <- cbind(aphids.df, group = mix.lab2)
ggplot(data = aphids.df2, aes(x, y, col = group, pch = group)) + geom_point() + labs(title = "Aphids Data (Semiparametric Mixture Fit; Symmetric Errors Assumption)", 
    x = "Number of Aphids in Batch", y = "Number of Infected Plants") + geom_abline(slope = spfit1$beta[, 
    1][2], intercept = spfit1$beta[, 1][1], col = 2) + geom_abline(slope = spfit1$beta[, 
    2][2], intercept = spfit1$beta[, 2][1], col = 3) + scale_color_manual("", values = c(2, 
    3), guide = F) + theme(text = element_text(size = 20), legend.position = "none")

### Estimate a semiparametric mixtures of regressions without assuming a symmetric
### error density.

set.seed(100)
spfit2 <- spregmix(n.inf ~ n.aphids, bw = 0.1, data = aphids, verb = TRUE, symm = FALSE)

mix.lab3 <- as.factor(apply(spfit2$posterior, 1, which.max) + 1)
aphids.df3 <- cbind(aphids.df, group = mix.lab3)
ggplot(data = aphids.df3, aes(x, y, col = group, pch = group)) + geom_point() + labs(title = "Aphids Data (Semiparametric Mixture Fit; No Symmetric Errors Assumption)", 
    x = "Number of Aphids in Batch", y = "Number of Infected Plants") + geom_abline(slope = spfit2$beta[, 
    1][2], intercept = spfit2$beta[, 1][1], col = 2) + geom_abline(slope = spfit2$beta[, 
    2][2], intercept = spfit2$beta[, 2][1], col = 3) + scale_color_manual("", values = c(2, 
    3), guide = F) + theme(text = element_text(size = 20), legend.position = "none")

list(out.em$beta, spfit1$beta, spfit2$beta)

### Look at the nonparametric KD-based estimate of the error densities fit above,
### as well as the parameteric normal density from the original fully parametric
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

###########################################################
### Example #4: Monitoring Wells Data
###########################################################

library(devtools)

source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/ggplottol.hist.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/ggplottol.control.R")

### Fit a one-sided (0.95, 0.95) gamma tolerance interval to the wells data.

wells <- c(5.1, 2.4, 0.4, 0.5, 2.5, 0.1, 6.8, 1.2, 0.5, 0.6, 5.3, 2.3, 1.8, 1.2, 
    1.3, 1.1, 0.9, 3.2, 1, 0.9, 0.4, 0.6, 8, 0.4, 2.7, 0.2, 2, 0.2, 0.5, 0.8, 2, 
    2.9, 0.1, 4)
gam.out <- gamtol.int(x = wells, alpha = 0.05, P = 0.95, side = 1, method = "EXACT")
gam.out

### Histogram and control chart of the results.

ggplottol.hist(gam.out, x = wells, side = "upper", x.lab = paste0("Vinyl Chloride Concentrations (", 
    "&#956;", "g/L)"))
ggplottol.control(gam.out, x = wells, side = "upper", x.lab = "Index")

###########################################################
### Example #5: Hospital Infections Data
###########################################################

source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/regtol.int2.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/ggplottol.reg.R")

hospitals <- read.table("https://online.stat.psu.edu/stat462/sites/onlinecourses.science.psu.edu.stat462/files/data/infectionrisk/index.txt", 
    header = TRUE)

### Calculate pointwise (0.90, 0.90) tolerance intervals for the multiple
### regression fit of the hospital infections data.

out <- lm(InfctRsk ~ Xray + Stay, data = hospitals)
out.TI <- regtol.int2(out, new.x = data.frame(Xray = c(50, 75), Stay = c(8, 12)), 
    alpha = 0.1, P = 0.9, side = 2, new = TRUE)
out.TI

ggplottol.reg(out.TI, x = hospitals[c("Xray", "Stay")], y = hospitals["InfctRsk"], 
    side = "two", rect = TRUE)

###########################################################
### Example #6: Adolescent Kidney Function Reference Regions
###########################################################

ref.X <- read.table("https://stat.as.uky.edu/sites/default/files/refX.txt", header = TRUE)
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


