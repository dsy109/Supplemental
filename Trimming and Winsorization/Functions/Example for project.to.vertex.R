###############
### Example ###
###############
### Create a convex hull ###
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
figure <- plot_ly() %>%
            add_trace(x=convex.hull[,1] , 
                      y=convex.hull[,2] , 
                      type = 'scatter' , mode = 'markers+lines' ,
                      marker = list(color = "#003DA5" , size = 15 , symbol = 'circle') ,
                      line = list(color = "#003DA5" , width = 5) ,
                      name = 'Vertices' , showlegend = TRUE) %>%
            add_trace(x=observations[,1] , 
                      y=observations[,2] , 
                      type = 'scatter' , mode = 'markers' ,
                      marker = list(color = "#008080" , size = 15 , symbol = 'pentagon') , 
                      name = 'Observations' , showlegend = TRUE) %>%
            add_trace(x=project[,1],
                      y=project[,2],
                      type = 'scatter' , mode = 'markers' ,
                      marker = list(color = "#9B111E" , size = 15 , symbol = 'diamond') , 
                      name = 'Projections' , showlegend = TRUE)
for (i in 1:(dim(observations)[1])){
  figure <- figure %>%
    add_segments(x=observations[i,1],xend=project[i,1],
                 y=observations[i,2],yend=project[i,2],
                 line = list(color = "#28a99e" , dash = "dash"),
                 showlegend = FALSE)
}
figure

