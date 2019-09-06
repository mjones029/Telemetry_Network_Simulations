## Cluster_Start_Simulations
# 
#========================================================	
# ---
### title: Clustered starting point simulations
# author: Marie Gilbertson
# date: "09/21/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates BCRW movement for population of 100 individuals for 90 days with clustered starting points.
# Clusters are from a negative binomial distribution, with mu set according to spacing between starting points also used for lattice configurations. 
# Number of clusters per simulation is either 1, 5, or 10, with even distribution of individuals between clusters. 
# Starting location PER CLUSTER is drawn from uniform distribution of 30,000 x 30,000 "study area" (as in random start simulations). 
#
# This script does not include contact detection or sampling; this code can be found in "Random_Start_Simulations.R"
# and is equivalent for all population spatial configurations.



##### Clear Environment
remove(list=ls())


#### Load R libraries	
library(CircStats)
library(adehabitatLT) 
library(fields)
library(reshape)
library(sp)
library(lubridate)
library(gdata) 
library(adehabitatHR) 
library(doParallel)
library(foreach)
library(plyr)
library(doRNG)


#---------------------------------------------------------
### Set simulation parameters based on array-id
#---------------------------------------------------------

# set array-id for this simulation
sim <- as.numeric(commandArgs(TRUE[1]))
sim.seed <- 32000 + sim  # ensure unique seeds for cluster start simulations



#---------------------------------------------------------------------------------------------
### Read in functions to simulate BCRW movement trajectories for a population of individuals
#---------------------------------------------------------------------------------------------

#Weighted circular mean calculation
source("w.circ.mean.R")

# BCRW simulation function
source("BCRW_sim.R")

#-------------------------------------------------------------------------------------------
### Set parameters and simulate BCRW movement trajectories for a population of individuals
#-------------------------------------------------------------------------------------------

# The following files are available in GitHub for reproducibility.
# set number of clusters
n.clust <- get(load("Number of Clusters.Rdata")) # read in dataset for number of clusters
n.clust <- n.clust[sim] # set number of clusters for this simulation

# set distances between starting locations (home range centers)
starts <- get(load("Cluster Sims Starting Locations_H15.Rdata")) # for setting mu in neg.bin distribution

# set shape of negative binomial distribution
sizes <- get(load("Cluster Start NegBin Sizes.Rdata"))

# set number of individuals per cluster (minus the starting individuals for each cluster)
# assumes total number of simulated trajectories will be 100 (for consistency, nsims should equal 100 in "Set Parameters" below)
np.clust <- NULL 
for(i in 1:n.clust){
  cl <- floor((100-n.clust)/n.clust)
  cls <- rep(cl, n.clust-1)
  cl.last <- (100-n.clust)-sum(cls)
  np.clust <- c(cls, cl.last)
}


##################################
######### Set Parameters #########
##################################
nsims <- 100 # number of individuals to simulate. 
final.duration <- 90 # duration (in days) of simualation
duration <- 7 + final.duration # account for initial "burn-in" period
x.dim <- 30000 # extent of starting study area in meters along x-axis
y.dim <- 30000 # extent of starting study area in meters along y-axis
h <- 15 # scaling parameter for step length
rho <- 0.8 # bias correlation parameter (0-1, where 0 -> unbiased, uncorrelated random walk, and 1 -> biased, deterministic movement)
b <- 0.01  # bias strength parameter (how fast bias increases or decreases)
c <- 0.3  # bias distance decay parameter (0 = no bias due to distance, >0 = bias increases w/distance, <0 = bias decreases w/distance)
##################################
##################################

# create list of starting points
all.starts <- NULL

for(i in 1:length(np.clust)){
  # distances between starting locations for each individual in cluster
  dists <- rnbinom(np.clust[i], mu = starts[sim], size = sizes[sim])
  
  # starting point for first individual in cluster
  x.start <- runif(1, min=0, max=x.dim)
  y.start <- runif(1, min=0, max=y.dim)
  xy.1 <- c(x.start, y.start)
  
  # draw direction from starting individual for starting points of remaining individuals in cluster
  # functionally, a wrapped uniform distribution
  thetas <- rwrpnorm(np.clust[i],0,0) 
  
  # calculate coordinates of starting points for remaining individuals in cluster
  x.list <- xy.1[1] + (dists*cos(thetas)) 
  y.list <- xy.1[2] + (dists*sin(thetas))
  round.starts <- cbind(data.frame(x.list), data.frame(y.list))
  round.starts <- rbind(xy.1, round.starts)
  colnames(round.starts) <- c("x", "y")
  round.starts$cluster.num <- paste(i)
  
  # save starting coordinates
  all.starts <- rbind(all.starts, round.starts)
}


# run the movement simulations, including setting up paralellization
cl <- makeCluster(3) 
registerDoParallel(cl)


clusterEvalQ(cl, {
  library(CircStats)
  library(adehabitatLT)
  library(lubridate)
  library(foreach)
})

# sim.seed ensures reproducibility and unique simulations, even with parallelization
all.sims <- foreach(n=1:nsims, .combine=rbind, .options.RNG=sim.seed) %dorng% { 
  x.sim <- all.starts$x[n]
  y.sim <- all.starts$y[n]
  simulation <- BCRW_sim(n=duration*60*24, h = h, rho = rho, y0=c(x.sim, y.sim), b = b, c = c)  
  simulation$id <- paste(n)
  simulation
}

stopCluster(cl)


# remove the first 7 days as "burn in" period
all.sims$step.id <- seq(1, 139681, 1)
remove.id <- seq(1, 10080, 1)
all.sims <- all.sims[!(all.sims$step.id %in% remove.id),]



#---------------------------------------------------------------------------------------
### Proceed with contact detection and sampling, found in "Random_Start_Simulations.R"
#---------------------------------------------------------------------------------------
