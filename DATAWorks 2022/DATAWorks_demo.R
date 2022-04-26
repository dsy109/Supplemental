################################################################
################################################################
################################################################
### R Code for "Computing Statistical Tolerance Regions Using 
### the R Package tolerance" by D. S. Young and K. Cheng
################################################################
################################################################
################################################################

library(tolerance)
library(plotly)
library(foreign)
library(AER)
library(devtools)
library(mixtools)

source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly_histtol.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly_controltol.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/regtol.int2.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly_regtol.R")
source_url("https://raw.githubusercontent.com/dsy109/tolerance/Updated-Functions/Updated%20Functions/plotly_npmvtol.R")

################################################################
################################################################
### Normal Tolerance Intervals
################################################################
################################################################

################################################################
### Example: Thickness of Silicon Wafers
################################################################

wafer <- c(unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/wafer.txt")))

wafer.tol <- normtol.int(x = wafer, alpha = 0.05, P = 0.95, side = 2, method = "EXACT") 
wafer.tol

plotly_histtol(wafer.tol, wafer, side = "two",
               x.lab = "Thickness of Metal Layer",
               tol.line.type = "dash")

plotly_controltol(wafer.tol, wafer, side = "two", 
                  x.lab = "Thickness of Metal Layer")

################################################################
################################################################
### Non-Normal Tolerance Intervals
################################################################
################################################################

################################################################
### Example: Defective Chips
################################################################

defects <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/defects.txt", header = TRUE)$Defects

defects.tol <- bintol.int(x = defects, n = 50*30, m = 50, alpha = 0.05, P = 0.99, side = 2, method = "CP") 
defects.tol

plotly_controltol(tol.out = defects.tol, x = defects,
                  x.lab = "Wafer",
                  y.lab = "Defectives",
                  tol.line.type = "dash")

################################################################
### Example: Serum Creatinine Levels
################################################################

kidney <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/kidney.txt", header = TRUE)

kidney.tol <- exttol.int(x = kidney$SCR, alpha = 0.05, P = 0.95, side = 2, dist = "Weibull") 
kidney.tol

plotly_histtol(kidney.tol, kidney$SCR, side = "two",
               x.lab = "Serum Creatinine",
               tol.line.type = "dash")

################################################################
### Example: Breast Cancer Survival Times
################################################################

distfree.est(alpha = 0.10, P = 0.90, side = 2)

bcancer <- unlist(read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/remission.txt", header = FALSE))

bcancer.WILKS <- nptol.int(x = bcancer, alpha = 0.10, P = 0.90, side = 2, method = "WILKS") 
bcancer.WILKS
plotly_histtol(bcancer.WILKS , x = bcancer , side ="two" ,
               x.lab = "Survival Time (Months)", tol.line.type = "dash")

################################################################
################################################################
### Regression Tolerance Intervals
################################################################
################################################################

################################################################
### Example: Hospital Infections
################################################################

hospitals <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/hospitals.txt", 
                        header = TRUE)

out <- lm(InfctRsk ~ Xray + Stay, data = hospitals)
out.TI <- regtol.int2(out, new.x = data.frame(Xray = 50, Stay = 12), 
                      alpha = 0.10, P = 0.90, side = 2, new = TRUE)
tail(out.TI$tol,5)

plotly_regtol(out.TI, x = hospitals[c("Xray", "Stay")], y = hospitals["InfctRsk"], 
              side = "two", rect = TRUE, x.cex = 8,
              x.lab = "X Ray", x.lab.size = 36, x.tick.size = 18,
              y.lab = "Stay", y.lab.size = 36, y.tick.size = 18,
              z.lab = "Infection Risk", z.lab.size = 36, z.tick.size = 18,
              title.size = 36 , title.position.y = 0.995)

################################################################
################################################################
### Multivariate Tolerance Regions
################################################################
################################################################

################################################################
### Example: Adolescent Kidney Function Reference Regions
################################################################

library(mixtools)
ref.X <- read.table("https://raw.githubusercontent.com/dsy109/Supplemental/main/CSP%202022%20Trainings/tolerance/Datasets/kidney_adolescents.txt", header = TRUE)
ref.males <- subset(ref.X, Gender == 1)[,-1]

Elliptical <- function(pts, x)
{
  out <- mixtools::depth(pts, x)
  out
}

tr.male <- npmvtol.region(x = as.matrix(ref.males), alpha = 0.05, P = 0.95, 
                          depth.fn = Elliptical, type = "semispace", adjust = "ceiling", semi.order = list(lower = NULL, 
                                                                                                           center = 2:3, upper = 1))
tr.male

plotly_npmvtol(tr.male, ref.males,
               var.names = colnames(ref.males), 
               title = "(0.95, 0.95) Nonparametric Rectangular Tolerance Regions",
               x.col = "#4298B5",
               x.cex = 6 ,
               x.shape = "dot",
               outlier.col = "#A6192E",
               outlier.cex = 8,
               outlier.shape = "triangle-up",
               tol.col = "#D1DDE6",
               tol.opacity = 0.4,
               x.lab.size = 12,
               x.tick.size = 12,
               y.lab.size = 12,
               y.tick.size = 12,
               title.position.x = 0.5,
               title.position.y = 0.98,
               title.size = 12,
               show.bound = TRUE , bound.type = "dot" , 
               bound.lwd = 4 , bound.col = "red")
