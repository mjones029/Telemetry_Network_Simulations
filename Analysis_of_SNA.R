## Analysis of SNA
# 
#========================================================	
# ---
# ## title: Analysis of social network analysis metrics
# author: Marie Gilbertson
# date: "01/03/2019"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Calculates ranked correlation for node-level SNA metrics between complete and sample networks
# Note: correlations are using only individuals from complete network that are in corresponding sample network;
# therefore, always comparing networks of the same size.
# 2. Exports correlation files for each node-level SNA metric for each simulation
# 3. Combines complete and sample metrics for network-level metrics into one dataset and exports (to be 
# plotted later)
#
# Social network analysis metrics: Degree, Strength, Betweenness, Transitivity,
# Density, Proportion Isolates, Modularity


##### Clear Environment
remove(list=ls())


#### Load R libraries	
library(igraph)
library(dplyr) # for left_join()
library(beepr)
library(ggplot2)

#### Set simulations to analyze
# the following are for setting up reading in files and looping through different simulation variations
# the following may therefore vary pending file naming system

nsims <- 500

start.type <- "Random Start" # set starting location type. Options are: "Random Start", "Lattice Start", or "Cluster Start"

h.type <- "H15" # set step length distribution. Options are: "H15", "H34", "H60", "SAC1", "SAC3", "SAC4"

# can compare complete and sample network metrics with different contact definitions
comp.cont.type <- "100m" # set contact threshold type for the COMPLETE network. Options are: "100m", "10m", or "1m"

samp.cont.type <- "100m" # set contact threshold type for the SAMPLE networks to compare to the complete network. Options are: "100m", "10m", or "1m"



### function for only keeping KDE results from q24h and q72h individual sampling levels
fix.KDE <- function(metric.data){
  s.t <- subset(metric.data, metric.data$contact.type=="space-time")
  kde <- subset(metric.data, metric.data$contact.type=="KDE UDOI")
  kde <- subset(kde, kde$ind.sample=="q24h" | kde$ind.sample=="q72h")
  new.data <- rbind(s.t, kde)
  return(new.data)
}


############ Calculate ranked correlation for each node-level metric ##############

# set node-level metrics and sampling levels to analyze
nl.metrics <- c("Deg", "Str", "Btw", "Clust") 
ind.sample <- c("q1m", "q15m", "q60m", "q3h", "q12h", "q24h", "q72h")
pop.sample <- seq(100, 10, -10)
contact.type <- c("space-time", "KDE UDOI")


# Loop through by simulation number
for(i in 1:nsims){

  #### Set simulation number
  print(paste("Simulation", i, sep = " "))
  sim <- i

  # loop through different node-level SNA metrics
  for(j in 1:length(nl.metrics)){
    metric <- nl.metrics[j]

    # read in complete network data for given metric
    # complete.name <- paste(<insert naming structure>, ".Rdata", sep = "")
    complete.metric <- get(load(file = complete.name))

    # read in sampled network data for given metric
    # sample.name <- paste(<insert naming structure>, ".Rdata", sep = "")
    sample.metric <- get(load(file = sample.name))

    # set up empty object to store results
    # uses dynamic allocation which is less efficient but functional for these purposes
    full.cor <- NULL

    # loop through different sampling levels for given metric and calculate the ranked correlation
    for(q in 1:length(contact.type)){
      ct <- contact.type[q]
      for(r in 1:length(pop.sample)){
        ps <- pop.sample[r]
        for(s in 1:length(ind.sample)){
         is <- ind.sample[s]

          samp.met_temp <- subset(sample.metric, sample.metric$ind.sample==is & sample.metric$pop.sample==ps & sample.metric$contact.type==ct)
          # pull out the metric values from the full network that match the ids of sampled individuals
          match_cmp.met <- complete.metric[complete.metric$id %in% samp.met_temp$id,]

          # calculate ranked correlation coefficient between the sampled metric calculation and the metric for those individuals from the complete network
          smp.cor <- data.frame(cor(samp.met_temp[,2], match_cmp.met[,2], method = "spearman"))
          colnames(smp.cor) <- "cor"
          # add sampling info for tracking purposes
          smp.cor$ind.sample <- is
          smp.cor$pop.sample <- ps
          smp.cor$contact.type <- ct
          smp.cor$metric <- metric
          # save for next round
          full.cor <- rbind(full.cor, smp.cor)

        }
      }
    }

    # only keep KDE results for q24h and q72h
    full.cor <- fix.KDE(full.cor)

    # add sim.num for tracking purposes
    full.cor$sim.num <- sim

    # save correlation data
    # cor.name <- paste(<insert naming structure>, ".Rdata", sep = "")
    save(full.cor, file = cor.name)

  }
}
beep(4)


########### Combine complete and sample metrics for all network-level metrics ##############

# set network-level metrics and sampling levels to analyze
nw.metrics <- c("Dens", "Iso", "Mod")
ind.sample <- c("q1m", "q15m", "q60m", "q3h", "q12h", "q24h", "q72h")
pop.sample <- seq(100, 10, -10)
contact.type <- c("space-time", "KDE UDOI")


# Loop through by simulation number
for(i in 1:nsims){
  
  #### Set simulation number
  print(paste("Simulation", i, sep = " "))
  sim <- i
  
  # loop through different network-level SNA metrics
  for(j in 1:length(nw.metrics)){
    metric <- nw.metrics[j]
    
    # read in complete network data for given metric
    # complete.name <- paste(<insert naming strucure>, sep = "")
    complete.metric <- get(load(file = complete.name))
    
    # read in sampled network data for given metric
    # sample.name <- paste(<insert naming strucure>, ".Rdata", sep = "")
    sample.metric <- get(load(file = sample.name))
    
    if(metric!="Mod"){
      # add metric value for the complete network to the sample dataset
      sample.metric$complete <- complete.metric
      # reorder columns for ease of assessment
      sample.metric <- sample.metric[,c(5, 1:4)]
      
      # add simulation number for tracking purposes
      sample.metric$sim.num <- sim
      
    }else{
      # modularity has results for several metrics, so needs to be assessed differently
      colnames(complete.metric) <- paste("c.", colnames(complete.metric), sep="")
      # add metric value for the complete network to the sample dataset
      sample.metric <- cbind(sample.metric, complete.metric)
      
      # add simulation number for tracking purposes
      sample.metric$sim.num <- sim
    }
    
    
    # save combined data
    # conc.name <- paste(<insert naming structure>, ".Rdata", sep = "")
    save(sample.metric, file = conc.name)
    
  }
}
beep(4)