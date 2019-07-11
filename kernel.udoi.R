## kernel.udoi
# 
#========================================================	
# ---
# ## title: Kernel UDOI
# author: Marie Gilbertson
# date: "12/19/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for calculating kernel UDOI spatial overlap contacts
# Spatial overlap code adapted from Brandell et al (in prep)

# z = locational dataset to evaluate
# kern.perc = kernel percentage to evaluate (e.g. 95)

kernel.udoi <- function(z, kern.perc){
  z.loc <- data.frame(x=z$x, y=z$y)
  z.sp <- SpatialPointsDataFrame(z.loc, z)
  z.udoi <- kerneloverlap(z.sp[,3], meth="UDOI", percent=kern.perc, grid=200) 
  return(z.udoi)
}
