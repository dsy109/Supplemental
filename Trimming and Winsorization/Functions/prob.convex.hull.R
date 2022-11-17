##########################################################################################
### This is a function to calculate the cumulative probability of a given convex hull 
### There are 2 function inputs. One is a data frame of a enclosed convex hull.
### The other is a function, defined as a joint distribution of X and Y (or any other function)
### The input is a data.frame of a convex hull. The data.frame should contain two columns.
### The left column is x values, and the right column is y values.
### Note here, "Points" should be a enclosed convex hull, that is to say the first and last row
### of "Point" should be exactly the same, since it is enclosed.
##########################################################################################
library(pracma)
###################################################################################
###################################################################################
prob.convex.hull <- function(points , f.xy){
  ### Data Impute ###
  points <- data.frame(points)
  sample.vertex <- points[-dim(points)[1],]
  x.order <- sample.vertex[order(sample.vertex[,1]),]
  
  ### Find Extremes ###
  location.xmin.org <- which(sample.vertex[,1] == min(sample.vertex[,1]))  
  location.xmax.org <- which(sample.vertex[,1] == max(sample.vertex[,1]))
  
  if (length(location.xmin.org)>1){
    if(sample.vertex[location.xmin.org[1],2] > sample.vertex[location.xmin.org[2],2]){
      location.xmin <- location.xmin.org[1]
    } else {
      location.xmin <- location.xmin.org[2]
    }
  } else {location.xmin <- location.xmin.org}
  
  if (length(location.xmax.org)>1){
    if(sample.vertex[location.xmax.org[1],2] > sample.vertex[location.xmax.org[2],2]){
      location.xmax <- location.xmax.org[2]
    } else {
      location.xmax <- location.xmax.org[1]
    }
  } else {location.xmax <- location.xmax.org}
  
  ### Connect Two Extremes ###
  sample.vertex.noExtreme <- sample.vertex[-c(location.xmin,
                                              location.xmax),]
  
  line.connect.extremes <- rep(NA,2)
  line.connect.extremes[1] <- as.numeric(((sample.vertex[location.xmax,])[2]-(sample.vertex[location.xmin,])[2])/
                                           ((sample.vertex[location.xmax,])[1]-(sample.vertex[location.xmin,])[1]))
  line.connect.extremes[2] <- as.numeric((sample.vertex[location.xmax,])[2]-(line.connect.extremes[1])*(sample.vertex[location.xmax,])[1])
  
  ### Determine Upper / Lower Vertex ###
  upper.vertex <- sample.vertex[location.xmin,]
  lower.vertex <- sample.vertex[location.xmax,]
  
  for (a in 1:(dim(sample.vertex.noExtreme)[1])){
    if (sample.vertex.noExtreme[a,2] < sample.vertex.noExtreme[a,1]*line.connect.extremes[1]+line.connect.extremes[2]){
      lower.vertex <- rbind(lower.vertex,sample.vertex.noExtreme[a,])
    } else {
      upper.vertex <- rbind(upper.vertex,sample.vertex.noExtreme[a,])
    }
  }
  # Record Original Order, in terms of X #
  origin.order.upper <- order(upper.vertex[,1])
  origin.order.lower <- order(lower.vertex[,1])
  
  upper.vertex <- upper.vertex[origin.order.upper,]
  lower.vertex <- lower.vertex[origin.order.lower,]
  
  ### Determine the functions for Boundaires ###
  upper.line.mat <- matrix(NA,nrow=dim(upper.vertex)[1],ncol=2)
  lower.line.mat <- matrix(NA,nrow=dim(lower.vertex)[1],ncol=2)
  
  if (dim(upper.vertex)[1]==1){
    upper.line.mat[1,1] <- line.connect.extremes[1]
    upper.line.mat[1,2] <- line.connect.extremes[2]
  } else {
    for (b in 1:(dim(upper.vertex)[1]-1)){
      upper.line.mat[b,1] <- (upper.vertex[b+1,2]-upper.vertex[b,2])/(upper.vertex[b+1,1]-upper.vertex[b,1])
      upper.line.mat[b,2] <- upper.vertex[b,2]-upper.line.mat[b,1]*upper.vertex[b,1]
    } 
    upper.line.mat[dim(upper.vertex)[1],1] <- (sample.vertex[location.xmax,2]-upper.vertex[dim(upper.vertex)[1],2])/(sample.vertex[location.xmax,1]-upper.vertex[dim(upper.vertex)[1],1])
    upper.line.mat[dim(upper.vertex)[1],2] <- upper.vertex[dim(upper.vertex)[1],2]-upper.line.mat[dim(upper.vertex)[1],1]*upper.vertex[dim(upper.vertex)[1],1]
  }
  ####
  if (dim(lower.vertex)[1]==1){
    lower.line.mat[1,1] <- line.connect.extremes[1]
    lower.line.mat[1,2] <- line.connect.extremes[2]
  } else {
    for (c in 2:(dim(lower.vertex)[1])){
      lower.line.mat[c,1] <- (lower.vertex[c,2]-lower.vertex[c-1,2])/(lower.vertex[c,1]-lower.vertex[c-1,1])
      lower.line.mat[c,2] <- lower.vertex[c,2]-lower.line.mat[c,1]*lower.vertex[c,1]
    } 
    lower.line.mat[1,1] <- (lower.vertex[1,2]-sample.vertex[location.xmin,2])/(lower.vertex[1,1]-sample.vertex[location.xmin,1])
    lower.line.mat[1,2] <- lower.vertex[1,2]-lower.line.mat[1,1]*lower.vertex[1,1]
  }
  
  ### Determine how many times each boundary section will be used ###
  upper.rep <- rep(NA,dim(upper.vertex)[1])
  lower.rep <- rep(NA,dim(lower.vertex)[1])
  
  ###
  if (length(upper.rep)==1){
    upper.rep[1] <- dim(lower.vertex)[1]
  } else {
    for (d in 2:length(upper.rep)){
      upper.rep[d-1] <- length(which(upper.vertex[d-1,1]<lower.vertex[,1] & lower.vertex[,1] < upper.vertex[d,1]))+1
    }
    upper.rep[length(upper.rep)] <- (length(unique(sample.vertex[,1]))-1)-sum(upper.rep[-length(upper.rep)])
  }
  
  ###
  if (length(lower.rep)==1){
    lower.rep[1] <- dim(upper.vertex)[1]
  } else {
    for (e in length(lower.rep):2){
      lower.rep[e] <- length(which(lower.vertex[e-1,1]<upper.vertex[,1] & upper.vertex[,1]<lower.vertex[e,1]))+1
    }
    lower.rep[1]<- (length(unique(sample.vertex[,1]))-1)-sum(lower.rep[-1])
  }
  
  ###
  upper.section.to.use <- rep(1,upper.rep[1])
  if (length(upper.rep)==1){
    upper.section.to.use <- rep(1,upper.rep[1])
  } else {
    for (f in 2:length(upper.rep)){
      upper.section.to.use<-c(upper.section.to.use , rep(f,upper.rep[f]))
    }
  }
  
  
  lower.section.to.use <- rep(1,lower.rep[1])
  if (length(lower.rep)==1){
    lower.section.to.use <- rep(1,lower.rep[1])
  } else {
    for (g in 2:length(lower.rep)){
      lower.section.to.use<-c(lower.section.to.use , rep(g,lower.rep[g]))
    }
  }
  
  section.to.use <- matrix(c(upper.section.to.use,lower.section.to.use),byrow = TRUE,nrow=2)
  ############
  prob.hull <- 0
  
  for (h in 1:(length(unique(sample.vertex[,1]))-1)){
    x.min <- unique(x.order[,1])[h]
    x.max <- unique(x.order[,1])[h+1]
    y.min <- function(x){x*as.numeric(lower.line.mat[section.to.use[2,h],1])+as.numeric(lower.line.mat[section.to.use[2,h],2])}
    y.max <- function(x){x*as.numeric(upper.line.mat[section.to.use[1,h],1])+as.numeric(upper.line.mat[section.to.use[1,h],2])}
    ### If there is any error due to boundary problem in the "pracma" package, then skip the simulation ###
    prob.hull <- try((prob.hull + integral2(fun=f.xy , xmin=x.min , xmax=x.max , ymin=y.min , ymax=y.max)$Q) , silent=TRUE)
    if(class(prob.hull)=="try-error"){
      prob.hull <- NA
    }
    ###
  }
  print(prob.hull)
}
