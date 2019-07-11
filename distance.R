## distance
# 
#========================================================	
# ---
# ## title: Distance between two points
# author: Marie Gilbertson
# date: "12/19/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for computing the distance between two points when given x and y coordinates 


distance <- function(x1, y1, x2, y2) {
  dist <- sqrt(((x2-x1)^2)+((y2-y1)^2))
  return(dist)
}