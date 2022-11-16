### There are two Inputs for the function: Vertex and Observations 
### Vertex: A matrix of 2 columns. 1st column is X values and 2nd column is Y values.
### Vertex: The 1st row and the last row should be the same, to have an enclosed convex hull.
### Vertex: Observations should be in order. The path of vertex should complete a convex hull.
### Observations: A matrix of 2 columns. 1st column is X values and 2nd column is Y values.
library(LearnGeom)
##################################################################
project.to.convex.hull <- function (vertex , observations){
  ######################################################
  ### We assume observations are all outside the convex hull ###
  ######################################################
  ### Determine if project to vertex or the boundary ###
  ######################################################
  ### Step 1: Find the nearest vertex ###
  n.observations <- dim(observations)[1]
  projections <- matrix(NA,ncol=2,nrow=n.observations)
  vertex <- vertex[-1,]
  
  for (i in 1:n.observations){
    combined <- rbind(observations[i,],vertex)
    distances <- dist(combined,method="euclidean")[1:(dim(vertex)[1])]
    closest <- which(distances == min(distances))
    ###########################
    if (length(closest) == 1){
      ### Pick up the nearest vertex, along with two neighbour vertex ###
      picked.location <- rep(NA,3)
      if (closest == 1){
        picked.location <- c(dim(vertex)[1],1,2)
      } else if (closest == dim(vertex)[1]) {
        picked.location <- c(closest-1 , closest , 1)
      } else {
        picked.location <- c(closest-1 , closest, closest+1)
      }
      ### Calculate Angles ###
      angles <- rep(NA,2)
      v1 <- c((observations[i,][1]-vertex[closest,][1]),(observations[i,][2]-vertex[closest,][2]))
      v2 <- c((vertex[picked.location[1],][1]-vertex[closest,][1]),
              (vertex[picked.location[1],][2]-vertex[closest,][2]))
      v3 <- c((vertex[picked.location[3],][1]-vertex[closest,][1]),
              (vertex[picked.location[3],][2]-vertex[closest,][2]))
      angles[1] <- (v1%*%v2)/(norm(as.matrix(v1))*norm(as.matrix(v2)))
      angles[2] <- (v1%*%v3)/(norm(as.matrix(v1))*norm(as.matrix(v3)))
      
      bound.or.vertex <- NA
      if (angles[1] > 0){
        bound.or.vertex <- (-1)
      } else if (angles[2] > 0){
        bound.or.vertex <- 1
      } else {
        bound.or.vertex <- 0
      }
      ### Project to vertex or boundary ###
      if (bound.or.vertex == 0){
        projections[i,] <- vertex[closest,]
      } else if (bound.or.vertex == (-1)){
        if (vertex[picked.location[1],][1] == vertex[closest,][1]){
          projections[i,] <- c(vertex[closest,][1] , observations[i,][2])
        } else if (vertex[picked.location[1],][2] == vertex[closest,][2]) {
          projections[i,] <- c(observations[i,][1] , vertex[closest,][2])
        } else {
          P1 <- vertex[picked.location[1],]
          P2 <- vertex[closest,]
          Line <- CreateLinePoints(P1,P2)
          projections[i,] <- ProjectPoint(observations[i,] , Line)
        }
      } else {
        if (vertex[picked.location[3],][1] == vertex[closest,][1]){
          projections[i,] <- c(vertex[closest,][1] , observations[i,][2])
        } else if (vertex[picked.location[3],][2] == vertex[closest,][2]) {
          projections[i,] <- c(observations[i,][1] , vertex[closest,][2])
        } else {
          P1 <- vertex[picked.location[3],]
          P2 <- vertex[closest,]
          Line <- CreateLinePoints(P1,P2)
          projections[i,] <- ProjectPoint(observations[i,] , Line)
        }
      }
    }
    ###########################
    else {
      if (vertex[closest[1],][1] == vertex[closest[2],][1]){
        projections[i,] <- c(vertex[closest[1],][1] , observations[i,][2])
      } else if (vertex[closest[1],][2] == vertex[closest[2],][2]) {
        projections[i,] <- c(observations[i,][1] , vertex[closest[1],][2])
      } else {
        P1 <- vertex[closest[1],]
        P2 <- vertex[closest[2],]
        Line <- CreateLinePoints(P1,P2)
        projections[i,] <- ProjectPoint(observations[i,] , Line)
      }
    }
  }
  return(projections)
}