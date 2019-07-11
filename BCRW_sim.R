## BCRW_sim
# 
#========================================================	
# ---
# ## title: BCRW Simulation Function
# author: Marie Gilbertson, modified from Long et al 2014, Journal of Animal Ecology
# date: "12/19/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates BCRW movement using a step length scaling parameter (h), bias correlation parameter (rho), bias strength parameter (b), and bias distance decay parameter (c)



BCRW_sim <- function(
  n=100,          #number of movement steps (trajectory will be subset by "minute")
  h=0.25,         #step length parameter
  rho=0,          #bias correlation parameter (0-1, where 0 -> unbiased, uncorrelated random walk, and 1 -> biased, deterministic movement)
  b=1,            #bias strength parameter
  c=0,            #bias distance decay parameter
  y0=c(0,0)       #animal starting location, also center of "home range"
){
  
  #---- Main Function ------
  y <- y0
  y.t <- y
  theta.y <- runif(1,-pi,pi)       #first direction is random
  
  
  for (i in 1:n){
    
    delta <- sqrt(sum((y0-y)^2))                              #distance to attraction point 
    psi <- atan2(y0[2]-y[2],y0[1]-y[1])                       #angle toward attraction point
    beta <- tanh(b*delta^c)                                   #bias effect     
    mu <- w.circ.mean(c(theta.y,psi),c(1-beta,beta))          #biased direction
    theta.y <- rwrpnorm(1,mu,rho)                             #"draw" actual turning angle based on "expected" angle, constrained by bias correlation parameter
    #step length from chi-squared distribution
    d.y <- h*rchi(1)
    y. <- y + c(d.y*cos(theta.y),d.y*sin(theta.y))            #calculate this "step"
    
    
    #Build the trajectory
    y.t <- rbind(y.t,y.)        
    #Save for next step!
    y <- y.
  }
  
  y.out <- data.frame(y.t,row.names=NULL)
  colnames(y.out) <- c("x","y")
  
  #add date/time to trajectory; considers date/time to be on a per-minute basis
  date <- seq(1, 60*(n+1), 60)
  y.out$date <- as_datetime(date)
  return(y.out)
}


#'----- end of BCRW_sim -----
