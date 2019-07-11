## simul.contacts
# 
#========================================================	
# ---
# ## title: Detect simultaneous contacts
# author: Marie Gilbertson
# date: "12/19/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for identifying simultaneous contacts (simultaneous points within the distance threshold)
# Note: requires "distance.R" function

simul.contacts <- function(dataset, id1, id2, dist.thresh){
  # calculate distance between simultaneous points
  temp1 <- subset(dataset, dataset$id==id1)
  temp2 <- subset(dataset, dataset$id==id2)
  colnames(temp2) <- c("x.2", "y.2", "date.2", "id.2", "step.id.2")
  temp <- cbind(temp1, temp2)
  
  temp$dist <- distance(temp$x, temp$y, temp$x.2, temp$y.2)
  
  # extract points within distance threshold for contact definition
  cont <- subset(temp, temp$dist <= dist.thresh)
  return(cont)
}
