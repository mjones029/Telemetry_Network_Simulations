## Random_Start_Simulations
# 
#========================================================	
# ---
# ## title: Random Start/Distributed Simulations
# author: Marie Gilbertson
# date: "08/13/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates BCRW movement for population of 100 individuals for 90 days
# 2. Performs simulations in parallel
# 3. Detects simultaneous contacts for all 100 simulations to construct "complete" network edgelist
# 4. Performs population and individual-level sampling (proportion of population sampled, and frequency of fixes)
# 5. Calculates both space-time overlap and spatial overlap contacts
#    (spatial overlap via KDE and UDOI for q24h and q72h sampling levels) to construct sample network edgelists

# Final output: Edgelists for complete and all sample networks, contacts (locations, time, ids) from the complete network, and IDs of sampled individuals at each sampling level


# Requires several other functions:
# a. w.circ.mean.R
# b. BCRW_sim.R
# c. distance.R
# d. simul.contacts.R
# e. sample.traj_df.R
# f. detect.contacts_list.R
# g. kernel.udoi.R
# h. sample_and_contact.R (links several other functions to facilitate looping)


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
### Set simulation number and seed based on array id
#---------------------------------------------------------

# Simulations were performed as job arrays through the Minnesota Supercomputing Institute
# Seeds were set based on job array/simulation number, but for demonstration purposes, can be set more generically

# set array-id and seed for this simulation if using a job array
# sim <- as.numeric(commandArgs(TRUE[1]))
# sim.seed <- 6000 + sim

# set seed generically
sim <- 1
sim.seed <- sim


#---------------------------------------------------------------------------------------------
# ## Read in functions to simulate BCRW movement trajectories for a population of individuals
#---------------------------------------------------------------------------------------------
# These functions are modified from Long et al, 2014 Journal of Animal Ecology

#Weighted circular mean calculation
source("w.circ.mean.R")

# BCRW simulation function
source("BCRW_sim.R")

#-------------------------------------------------------------------------------------------
# ## Set parameters and simulate BCRW movement trajectories for a population of individuals
#-------------------------------------------------------------------------------------------
#
##################################
######### Set Parameters #########
##################################
nsims <- 100 # number of individuals to simulate
x.dim <- 30000 # extent of starting study area in meters along x-axis
y.dim <- 30000 # extent of starting study area in meters along y-axis
final.duration <- 90 # duration (in days) of final simualation
duration <- 7 + final.duration # account for initial "burn-in" period as all simulations start from what will be the home range center
h <- 15 # scaling parameter for step length
rho <- 0.8 # bias correlation parameter (0-1, where 0 -> unbiased, uncorrelated random walk, and 1 -> biased, deterministic movement)
b <- 0.01  # bias strength parameter (how fast bias increases or decreases)
c <- 0.3  # bias distance decay parameter (0 = no bias due to distance, >0 = bias increases w/distance, <0 = bias decreases w/distance)
##################################
##################################

# Set up paralellization
cl <- makeCluster(3) # Set number of cores to use for parallelization 
registerDoParallel(cl)

# "Export" necessary packages to other cores
clusterEvalQ(cl, {
  library(CircStats)
  library(adehabitatLT)
  library(lubridate)
  library(foreach)
})

# Run simulations based on previous parameters. 
# Note that the seed setting is necessary for reproducibility across parallel simulations

all.sims <- foreach(n=1:nsims, .combine=rbind, .options.RNG=sim.seed) %dorng% {
  x.sim <- runif(1, min=0, max=x.dim)
  y.sim <- runif(1, min=0, max=y.dim)
  simulation <- BCRW_sim(n=duration*60*24, h = h, rho = rho, y0=c(x.sim, y.sim), b = b, c = c)  
  simulation$id <- paste(n)
  simulation
}

stopCluster(cl)

# remove the first 7 days as "burn in" period
all.sims$step.id <- seq(1, 139681, 1)
remove.id <- seq(1, 10080, 1)
all.sims <- all.sims[!(all.sims$step.id %in% remove.id),]



#------------------------------------------------------------------
### Find contacts for "complete" network
## "Complete" network contacts are simultaneous points w/in a distance threshold
#------------------------------------------------------------------

##################################
######### Set Parameters #########
##################################

dist.thresh <- 100 # sets distance threshhold for a contact in meters 

##################################
##################################

# function for calculating the euclidean distance between two points
source("distance.R")

# function for identifying simultaneous contacts (simultaneous points within the distance threshold)
source("simul.contacts.R")


# prep for looping
id <- unique(all.sims$id)
contacts <- NULL # Uses dynamic allocation, which is generally not recommended for streamlining purposes. However, this allows for storing results of varying length.


# loop through whole population for pairwise contact detection
for(i in 1:(length(id)-1)){
  print(i)
  ind1 <- id[i]
  
  for(j in (i+1):length(id)){  
    ind2 <- id[j]
    
    round.contacts <- simul.contacts(all.sims, ind1, ind2, dist.thresh) 
    contacts <- rbind(contacts, round.contacts) # stores the contact ids, times, and locations for each pair of individuals
  }
}


# extract edge list from contact data and save
edge.list <- data.frame(contacts$id, contacts$id.2)
edge.list <- droplevels(unique(edge.list))
colnames(edge.list) <- c("individual.1","individual.2")
if(nrow(edge.list)==0){edge.list[1,] <- NA}
# save data about simulation/contacts for later subsetting
edge.list$value <- NA
edge.list$ind.sample <- "q1m"
edge.list$pop.sample <- "Complete"
edge.list$sim.num <- sim
edge.list$contact.type <- "space-time"

# output complete network edge list
# save according to preferred naming structure
# complete.name <- paste(<insert naming structure>, ".Rdata", sep = "")
# save(edge.list, file=complete.name)




# output contact locations/data
# save according to preferred naming structure
# contact.name <- paste(<insert naming structure>, ".Rdata", sep = "")
# save(contacts, file=contact.name)




#------------------------------------------------------------------
### Combine population and individual-level sampling of trajectories
### with space-time and spatial overlap contact identification
#------------------------------------------------------------------

# load function for subsetting dataframe of trajectories for population and individual-level sampling
source("sample.traj_df.R")

# load function for detecting space-time contacts within a dataframe of simulated trajectories
# DIFFERENT FUNCTION THAN FOR DETECTING COMPLETE NETWORK
source("detect.contacts_list.R")

# load function for calculating kernel UDOI contacts
# Spatial overlap code adapted from Ellen Brandell et al (in prep)
source("kernel.udoi.R")

# read in function that performs sampling (via sample.traj_df function) and contact detection (via detect.contacts_list and kernel.udoi functions) from sub-sampled dataset
source("sample_and_contact.R")

##################################
######### Set Parameters #########
##################################
max.time <- final.duration*24*60+1 # (in "minutes") Sets maximum time point to the last time stamp in the dataset about to be sampled.
pop.size <- length(unique(all.sims$id)) # Calculates the total population size from which to sample
dataset <- all.sims
edge.list <- edge.list
sim <- sim

pop.sample.levels <- seq(1, 0.1, -0.1) # proportions of population to sample for population-level sampling

# Uses the same distance threshold as complete network contact detection, but could use different thresholds
dist.thresh <- dist.thresh # sets distance threshhold for a contact in meters 
kern.perc <- 95 # sets percent kernel for calculating spatial overlap (KDE UDOI)
##################################
##################################



# loop through all levels of population sampling and save edgelist from each iteration
# saves at the level of population sampling; each pop. sample includes all levels of individual-level sampling and both space-time and spatial overlap definitions of contact
for(f in 1:length(pop.sample.levels)){
  print(pop.sample.levels[f])
  pop.sample <- pop.sample.levels[f]
  iter.edges <- sample_and_contact(pop.sample=pop.sample, 
                                   dataset=dataset, 
                                   max.time=max.time, 
                                   pop.size=pop.size, 
                                   edge.list=edge.list, 
                                   sim=sim, 
                                   dist.thresh=dist.thresh, 
                                   kern.perc=kern.perc)
}


# check for any warnings at end of simulation
warnings()

# Next step: Calculate complete and sample network metrics using output (edgelists) from this script.