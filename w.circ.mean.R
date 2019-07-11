## w.circ.mean
# 
#========================================================	
# ---
# ## title: Weighted circular mean calculation
# author: Marie Gilbertson from Long et al 2014, Journal of Animal Ecology
# date: "08/13/2018"
#---


#Weighted circular mean calculation
w.circ.mean <- function (x,w) 
{
  sinr <- sum(w*sin(x)) 
  cosr <- sum(w*cos(x)) 
  circmean <- atan2(sinr, cosr) 
  circmean 
}
