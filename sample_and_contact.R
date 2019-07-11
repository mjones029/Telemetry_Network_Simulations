## sample_and_contact
# 
#========================================================	
# ---
# ## title: Sampling and Contact Detection
# author: Marie Gilbertson
# date: "12/22/2018"
#---
# ##  Preamble	
# NOTE: File naming scheme should be updated and uncommented before using (lines 102-109)

# What this code does:
# 1. Performs population and individual level sampling
# 2. Detects contacts for sample networks with both space-time and spatial overlap definitions of contact
# Note: this code combines multiple functions and processes to facilitate looping through different sampling levels in main code body.
# Requires the following functions:
# a. sample.traj_df.R
# b. detect.contacts_list.R
# c. kernel.udoi.R

# pop.sample = proportion of population to sample
# dataset = complete trajectory dataset to evaluate
# max.time = maximum trajectory time point
# pop.size = total population size
# edge.list = complete network edge list
# sim = simulation number
# dist.thresh = distance threshold for simultaneous contacts
# kern.perc = kernel percentage for spatial overlap (e.g. 95%)

sample_and_contact <- function(pop.sample=pop.sample, dataset=dataset, max.time=max.time, pop.size=pop.size, edge.list=edge.list, sim=sim, dist.thresh=dist.thresh, kern.perc=kern.perc){

  # sample from original dataset
  pop.sample <- pop.sample
  sampled.data <- sample.traj_df(dataset=dataset, pop.level=pop.sample, max.time=max.time, pop.size=pop.size)

  # calculate space-time contacts for q1m sample
  sample.ids <- unique(sampled.data[[1]]$id) # which individuals are in the population sample?
  q1m <- edge.list[edge.list$puma.1 %in% sample.ids & edge.list$puma.2 %in% sample.ids,]
  if(nrow(q1m)==0){q1m[1,] <- NA}
  q1m$ind.sample <- "q1m"
  q1m$pop.sample <- pop.sample*100
  q1m$sim.num <- sim
  q1m$contact.type <- "space-time"

  # calculate space-time contacts on all other data samples
  iter.st.contacts <- detect.contacts_list(full.dataset=sampled.data, dist.thresh=dist.thresh)
  iter.st.contacts <- na.omit(iter.st.contacts)
  if(nrow(iter.st.contacts)==0){iter.st.contacts[1,] <- NA}
  iter.st.contacts$pop.sample <- pop.sample*100
  iter.st.contacts$sim.num <- sim
  iter.st.contacts$contact.type <- "space-time"

  iter.st.contacts <- rbind(q1m, iter.st.contacts)


  # calculate spatial overlap (UDOI) contacts on all data samples
  z <- sampled.data[[5]] # start with only the q24h sampled data
  z <- z[,c(1:2,4)] # keep only the locations and individual ID for use in spatial points dataframe
  q24 <- kernel.udoi(z, kern.perc=kern.perc)
  individuals <- rownames(q24)
  y <- combn(individuals,2) # generate all pairs of individuals

  # put results into data frame
  q24.so.results <- NULL

  for(j in 1:dim(y)[2]){
    q24.so.results$puma.1[j] <- paste(y[1,j])
    q24.so.results$puma.2[j] <- paste(y[2,j])
  }
  q24.so.results <- data.frame(q24.so.results)
  q24.so.results$value <- lowerTriangle(q24, diag=F)
  q24.so.results$ind.sample <- "q24h"

  # repeat above UDOI process for q72h sample
  z <- sampled.data[[6]] # extracts only the q72h sample
  z <- z[,c(1:2,4)] # again, keep only locations and individual ID's for spatial points dataframe
  q72 <- kernel.udoi(z, kern.perc=kern.perc)
  q72.so.results <- NULL

  for(k in 1:dim(y)[2]){
    q72.so.results$puma.1[k] <- paste(y[1,k])
    q72.so.results$puma.2[k] <- paste(y[2,k])
  }
  q72.so.results <- data.frame(q72.so.results)
  q72.so.results$value <- lowerTriangle(q72, diag=F)
  q72.so.results$ind.sample <- "q72h"

  # combine spatial overlap edgelists
  iter.so.contacts <- rbind(q24.so.results, q72.so.results)

  # only keep individuals with UDOI>0 (UDOI=0 means no contact)
  iter.so.contacts <- droplevels(subset(iter.so.contacts, iter.so.contacts$value>0))
  if(nrow(iter.so.contacts)==0){iter.so.contacts[1,] <- NA}
  iter.so.contacts$pop.sample <- pop.sample*100
  iter.so.contacts$sim.num <- sim
  iter.so.contacts$contact.type <- "KDE UDOI"


  # combine space-time and spatial-overlap edgelists for this round of population sampling and save
  iter.edges <- rbind(iter.st.contacts, iter.so.contacts)

  # output simulation edgelists as one file using preferred naming scheme
  # edge.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # save(iter.edges, file=edge.name)

  # save ids of sampled/monitored individuals
  samp.ids <- sampled.data[[7]]
  # samp.ids.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # save(samp.ids, file=samp.ids.name)

}