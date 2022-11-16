\name{project.to.vertex}
\title{Projection to Vertices of a Convex Hull}
\alias{project.to.vertex}
\usage{
project.to.vertex(vertex , observations)
}

\description{
\code{project.to.vertex} finds projections of any external point to the nearest vertices of the convex hull.
}

\arguments{
  \item{vertex}{A matrix of coordinates of vertices of the convex hull. The matrix should have 2 colums, where the first column is the x-coordinate and the second column is the y-coordiante of vertices. That is, each row represents for one point. The last row should be identical to the first row to enclose the convex hull.}
  
  \item{observations}{A matrix of coordinates of observations. The matrix should have 2 colums, where the first column is the x-coordinate and the second column is the y-coordiante of each observation. That is, each row represents for one observation.}
}

\value{
  \code{project.to.vertex} returns a data frame with items:
  \item{X}{x-coordinate of each projection point.}
  \item{Y}{y-coordinate of each projection point.}
  Each row represents for one projection point. Number of rows of the outcome of \code{project.to.vertex} equals to number of rows of \code{observations}. Each row represents for the corresponding projection point of the \code{observations}.
}

\seealso{
  \code{\link{vertex.to.convex.hull}}
}

\examples{### Create a convex hull ###
convex.hull <- matrix(c(4,0,
                        5,1,
                        5,2,
                        4,3,
                        3,2,
                        3,1,
                        4,0),
                      byrow=TRUE,ncol=2,nrow=7)

### Create observations ###
observations <- matrix(c(3.2,2.25,
                         6,0,
                         -7,2.4,
                         0,0),byrow=TRUE,ncol=2,nrow=4)

### Find projections from observations to vertice of the convex hull ###
project <- project.to.vertex(vertex = convex.hull ,
                             observations = observations)

### Visualization ###
library(plotly)
figure <- plot_ly() \%>\%
            add_trace(x=convex.hull[,1] , 
                      y=convex.hull[,2] , 
                      type = 'scatter' , mode = 'markers+lines' ,
                      marker = list(color = "#003DA5" , size = 15 , symbol = 'circle') ,
                      line = list(color = "#003DA5" , width = 5) ,
                      name = 'Vertices' , showlegend = TRUE) \%>\%
            add_trace(x=observations[,1] , 
                      y=observations[,2] , 
                      type = 'scatter' , mode = 'markers' ,
                      marker = list(color = "#008080" , size = 15 , symbol = 'pentagon') , 
                      name = 'Observations' , showlegend = TRUE) \%>\%
            add_trace(x=project[,1],
                      y=project[,2],
                      type = 'scatter' , mode = 'markers' ,
                      marker = list(color = "#9B111E" , size = 15 , symbol = 'diamond') , 
                      name = 'Projections' , showlegend = TRUE)
for (i in 1:(dim(observations)[1])){
  figure <- figure \%>\%
    add_segments(x=observations[i,1],xend=project[i,1],
                 y=observations[i,2],yend=project[i,2],
                 line = list(color = "#28a99e" , dash = "dash"),
                 showlegend = FALSE)
}
figure
}
\keyword{file}
