### Vertex: a data frame of vertex of convex hulls. Starting and Ending Points are the same. ###
### Observations: a data frame of observations, 1st column is X values, 2nd column is Y values. ###
### Result is a matrix, 1st column is X values of projection, 2nd column is Y values of projection. ###
project.to.vertex <- function(vertex , observations){
  n.observations <- dim(observations)[1]
  projections <- matrix(NA,nrow=n.observations,ncol=2)
  vertex <- vertex[-1,]
  ######
  for (i in 1:n.observations){
    combined <- rbind(observations[i,],vertex)
    distances <- dist(combined,method="euclidean")[1:(dim(vertex)[1])]
    closest <- which(distances == min(distances))
    projections[i,] <- as.matrix(vertex[closest,])
  }
  ######
  return(data.frame(projections))
}
