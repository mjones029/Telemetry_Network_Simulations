## sample.traj_df
# 
#========================================================	
# ---
# ## title: Sample trajectory dataframe
# author: Marie Gilbertson
# date: "12/18/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for subsetting dataframe of trajectories for population and individual-level sampling
# (i.e. proportion of population sampled and frequency of recorded observations/locations)

# dataset = data frame of trajectories to assess
# pop.level = proportion of population to sample
# max.time = maximum simulation time point for sampling any individual trajectory
# pop.size = total population size


sample.traj_df <- function(dataset, pop.level, max.time, pop.size){
  if(pop.level <=0 | pop.level > 1){
    stop("Proportion of population sampled must be greater than 0 or less than 1")
  }
  if(class(dataset)[1] != "data.frame"){
    stop("Dataset must be an object of class data.frame")
  }
  
  # step.id should already be in dataset
  
  # population-level sampling
  pop.sample <- sort(sample(pop.size, size=paste(pop.level*pop.size), replace=F)) # randomly choose individuals to "collar"; use same set of individuals within a round of "collaring"
  sub_by.pop <- dataset[dataset$id %in% pop.sample,] # subset whole dataset to include only "collared" animals
  
  # individual-level sampling
  # all intervals will start at step.id #1
  
  # every 15min
  int.15 <- seq(1, max.time, 15) # set time step length to subset the trajectories
  sub_by.int.15 <- sub_by.pop[sub_by.pop$step.id %in% int.15,] # subset dataset of collared individuals by individual-level sampling time interval
  
  # every 60min
  int.60 <- seq(1, max.time, 60) # set time step length to subset the trajectories
  sub_by.int.60 <- sub_by.pop[sub_by.pop$step.id %in% int.60,] # subset dataset of collared individuals by individual-level sampling time interval
  
  # every 3 hours
  int.3 <- seq(1, max.time, 60*3) # set time step length to subset the trajectories
  sub_by.int.3 <- sub_by.pop[sub_by.pop$step.id %in% int.3,] # subset dataset of collared individuals by individual-level sampling time interval
  
  # every 12 hours
  int.12 <- seq(1, max.time, 60*12) # set time step length to subset the trajectories
  sub_by.int.12 <- sub_by.pop[sub_by.pop$step.id %in% int.12,] # subset dataset of collared individuals by individual-level sampling time interval
  
  # every 24 hours
  int.24 <- seq(1, max.time, 60*24) # set time step length to subset the trajectories
  sub_by.int.24 <- sub_by.pop[sub_by.pop$step.id %in% int.24,] # subset dataset of collared individuals by individual-level sampling time interval
  
  # every 72 hours
  int.72 <- seq(1, max.time, 60*72) # set time step length to subset the trajectories
  sub_by.int.72 <- sub_by.pop[sub_by.pop$step.id %in% int.72,] # subset dataset of collared individuals by individual-level sampling time interval
  
  output <- list(sub_by.int.15, sub_by.int.60, sub_by.int.3, sub_by.int.12, sub_by.int.24, sub_by.int.72, pop.sample)
  names(output) <- c("Every 15 min", "Every 60 min", "Every 3 hours", "Every 12 hours", "Every 24 hours", "Every 72 hours", "Population Sample IDs")
  return(output)
}
