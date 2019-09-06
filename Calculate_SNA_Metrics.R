## Calculate_SNA_Metrics
# 
#========================================================	
# ---
# ## title: SNA Metrics Calculation
# author: Marie Gilbertson
# date: "08/20/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Calculates and saves all SNA metrics (see below for list) for complete and all sample networks for each simulation

# Social network analysis metrics: Degree, Strength, Betweenness, Clustering Coefficient/Transitivity,
# Density, Proportion Isolates, Modularity


##### Clear Environment #####
remove(list=ls())


#### Load R libraries	####
library(igraph)
library(dplyr) # for left_join()
library(beepr) # optional


#### Set simulations to analyze ####
# the following are for setting up reading in files and looping through different simulation variations
# the following may therefore vary pending file naming system

nsims <- 500 # set number of simulations

start.type <- "Random Start" # set starting location type (e.g. Random, Lattice, Cluster, depending upon file naming structure)

h.type <- "H15" # set step length distribution. Options are: "H15", "H34", "H60", "SAC1", "SAC3", "SAC4"

cont.type <- "100m" # set contact threshold type. Options are: "100m", "10m", or "1m"


#### load extra functions ####

# function for keeping only KDE/UDOI data from sampling frequency levels of q24h and q72h
fix.KDE <- function(metric.data){
  s.t <- subset(metric.data, metric.data$contact.type=="space-time")
  kde <- subset(metric.data, metric.data$contact.type=="KDE UDOI")
  kde <- subset(kde, kde$ind.sample=="q24h" | kde$ind.sample=="q72h")
  new.data <- rbind(s.t, kde)
  return(new.data)
}



#------------- Start of primary loop - each loop analyzes one simulation -----------#

for(z in 1:nsims){
  
  #### Set simulation number
  print(paste("Simulation", z, sep = " "))
  sim <- z
  
  
  ###### Read in Complete Edgelist Data ######
  # el.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  edge.list <- get(load(el.name))
  
  
  ###### Calculate "value" (weight) for q1m edgelist ##########
  
  # load contacts file from simulation output
  # contacts.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  contacts <- get(load(contacts.name))

  # loop through all dyads taking the number of observed contacts as the "value" for weighting edges
  for(i in 1:nrow(edge.list)){
    temp.dyad <- contacts[contacts$id==edge.list$individual.1[i] & contacts$id.2==edge.list$individual.2[i],]
    edge.list$value[i] <- nrow(temp.dyad) 
  }
  
  # convert to character class to avoid factor problems with matching/merging
  edge.list$individual.1 <- as.character(edge.list$individual.1) 
  edge.list$individual.2 <- as.character(edge.list$individual.2)
  
  
  ###### Read in Sample Edgelist and Observed Individual Data ######
  
  # loop to read in and combine all sampling edgelists, as well as the selection of individuals that were sampled in population-level sampling for each sample
  
  # "samples" below reflects the naming structure the author used as each sample edgelist is saved by proportion of the population sampled
  # "samples" will be used to load each sampling edgelist in the following loop
  samples <- c("100Samp", "90Samp", "80Samp", "70Samp", "60Samp", "50Samp", "40Samp", "30Samp", "20Samp", "10Samp")
  
  sample.el <- NULL # uses dynamic allocation which is not ideal for streamlining, but functional here
  sample.ids <- list(NULL)
  
  # loops through population sampling level to load in all sampling edgelists and id's of sampled individuals
  for(j in 1:length(samples)){   
    print(j)
    samp <- samples[j]
   
    # load sample edgelist file; sample edgelists will have been saved by the proportion of the population sampled
    # file.name <- paste(<insert naming structure>, ".Rdata", sep = "")
    el <- get(load(file.name))

    el$individual.1 <- as.character(el$individual.1) # convert to character class to avoid factor problems with matching/merging
    el$individual.2 <- as.character(el$individual.2)
    
    # assign "value" for q1m subset (draw from edge.list object)
    # q1m level had NA for value in initial simulations, so processed here differently than other space-time sampling levels
    q1m <- subset(el, el$ind.sample=="q1m")
    q1m <- left_join(q1m, edge.list, by=c("individual.1", "individual.2"))
    q1m <- q1m[,c(1,2,8,4:7)] # only keep columns: individual.1, individual.2, value.y (value from edge.list), ind.sample, pop.sample, sim.num, contact.type
    colnames(q1m) <- c("individual.1","individual.2","value","ind.sample","pop.sample","sim.num","contact.type")
    
    # bind q1m data with the rest of the sampled edgelist data
    rest <- subset(el, el$ind.sample!="q1m")
    el <- rbind(q1m, rest)
    
    # save edgelist output in dynamic allocation
    sample.el <- rbind(sample.el, el)
    
    # read in ids of those individuals included in the population sample (for identifying isolates in the network later)
    # author's file naming system again used the "samples" object to load sample IDs by proportion of the population sampled 
    # sample.file.name <- paste(<insert naming structure>, ".Rdata", sep = "")
    smp <- get(load(sample.file.name))

    # save (via dynamic allocation) the IDs of the individuals sampled at the current sampling level
    sample.ids[[j]] <- smp
  }
  
  # join complete and all sample network edgelists into one object
  full.el <- rbind(edge.list, sample.el)
  full.el <- full.el[!(is.na(full.el$individual.1)),] # make sure no blank rows were included
  
  
  # Use "value" column to assign weights to edges for subsequent social network analysis.
  # Can use direct or inverse weights
  full.el$weight <- full.el$value # need the values to be called "weights" for subsequent SNA
  
  
  ###### Calculate metrics for COMPLETE network ##########
  
  # Take only the edgelist for the complete network
  temp.el <- subset(full.el, full.el$ind.sample=="q1m" & full.el$pop.sample=="Complete" & full.el$contact.type=="space-time")
  net<-graph.data.frame(temp.el,directed=F) # convert edgelist to igraph object
  

  ### Degree ###
  c.deg <- degree(net,v=V(net), normalized = T) ## calculate degree; normalized because complete and sample networks are different sizes
  # put degree results in dataframe including individuals that had degree 0
  c.deg.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
  colnames(c.deg.df) <- "id"
  c.deg.df$deg <- c.deg
  c.deg.df <- c.deg.df[order(c.deg.df$id),]

  full.deg <- data.frame(
    id=c(1:100),
    deg=NA
  )

  have.deg <- c.deg.df$id

  # add isolates to the dataset and assign them a degree of 0
  for(m in 1:100){
    if(paste(m) %in% have.deg){
      full.deg[m,2] <- c.deg.df[c.deg.df$id==paste(m), 2]
    }
    else{
      full.deg[m,2] <- 0
    }
  }
  
  # final product: full.deg, which is a dataframe of degree values for the complete network

  
  ### Strength ###
  c.str <- strength(net, vids = V(net)) # weights automatically incorporated because of column "weights" in original data frame

  # the following lines turn the strength calculations into a workable dataframe
  c.str.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
  colnames(c.str.df) <- "id"
  c.str.df$str <- c.str
  c.str.df <- c.str.df[order(c.str.df$id),]

  full.str <- data.frame(
    id=c(1:100),
    str=NA
  )

  have.str <- c.str.df$id

  # add isolates to the dataset and assign them a strength of 0
  for(n in 1:100){
    if(paste(n) %in% have.str){
      full.str[n,2] <- c.str.df[c.str.df$id==paste(n), 2]
    }
    else{
      full.str[n,2] <- 0
    }
  }

  # final product: full.str, which is a dataframe of strength (weighted degree) values for the complete network


  ### Betweenness ###
  c.btw <- betweenness(net, v=V(net)) # If "weights" attribute present, will calculate WEIGHTED betweenness
  
  # the following lines turn the betweenness calculations into a workable dataframe
  c.btw.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
  colnames(c.btw.df) <- "id"
  c.btw.df$btw <- c.btw
  c.btw.df <- c.btw.df[order(c.btw.df$id),]
  
  full.btw <- data.frame(
    id=c(1:100),
    btw=NA
  )
  
  have.btw <- c.btw.df$id
  
  # add isolates to the dataset and assign them a betweenness of 0
  for(n in 1:100){
    if(paste(n) %in% have.btw){
      full.btw[n,2] <- c.btw.df[c.btw.df$id==paste(n), 2]
    } 
    else{
      full.btw[n,2] <- 0
    }
  }
  
  # final product: full.btw, which is a dataframe of betweenness values for the complete network
  
  
  ### Transitivity ###
  c.clust <- transitivity(net, type = c("local"), vids = V(net), isolates = c("zero")) # no normalization option; if weight included, calculates weighted metric by default
  # choice of "zero" for isolates means isolates and deg=1 individuals have transitivity of 0 (plus individuals with calculated transitivity = 0)

  # the following lines turn the transitivity calculations into a workable dataframe
  c.clust.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
  colnames(c.clust.df) <- "id"
  c.clust.df$clust <- c.clust
  c.clust.df <- c.clust.df[order(c.clust.df$id),]

  full.clust <- data.frame(
    id=c(1:100),
    clust=NA
  )

  have.clust <- c.clust.df$id

  # add isolates to the dataset and assign them a transitivity of 0
  for(n in 1:100){
    if(paste(n) %in% have.clust){
      full.clust[n,2] <- c.clust.df[c.clust.df$id==paste(n), 2]
    }
    else{
      full.clust[n,2] <- 0
    }
  }

  # final product: full.clust, which is a dataframe of transitivity values for the complete network


  ### Density ###
  full.dens <- nrow(edge.list)/choose(100,2) # if population size is greater than 100, this line would need to be altered accodingly

  # final product: full.dens, which is the density for the complete network


  ### Proportion of population which are isolates ###
  full.iso <- (nrow(subset(full.deg, full.deg$deg==0)))/100 # if population size is greater than 100, this line would need to be altered accordingly

  # final product: full.iso, which is the proportion of individuals that were isolates in the network (degree=0).

  ### Modularity ###
  # calculate several metrics of weighted and unweighted modularity
  if(gsize(net)>0){ # need to account for unconnected networks
    wtc.w <- cluster_walktrap(net, weights = E(net)$value, steps = 4, modularity = T)
    wtc.uw <- cluster_walktrap(net, weights = NULL, steps = 4, modularity = T)
    
    eb.w <- cluster_edge_betweenness(net, weights = E(net)$value, directed = F)
    eb.uw <- cluster_edge_betweenness(net, weights = NULL, directed = F)
    
    fg.w <- cluster_fast_greedy(net, weights = E(net)$value)
    fg.uw <- cluster_fast_greedy(net, weights = NULL)
    
    mod.wt.w <- modularity(net, membership(wtc.w))
    mod.wt.uw <- modularity(net, membership(wtc.uw))
    mod.eb.w <- modularity(net, membership(eb.w))
    mod.eb.uw <- modularity(net, membership(eb.uw))
    mod.fg.w <- modularity(net, membership(fg.w))
    mod.fg.uw <- modularity(net, membership(fg.uw))
    
    full.mod <- c(mod.wt.w, mod.wt.uw, mod.eb.w, mod.eb.uw, mod.fg.w, mod.fg.uw)
  } else {
    full.mod <- NA
  }
  full.mod <- data.frame(matrix(full.mod, ncol = 6))
  colnames(full.mod) <- c("mod.wt.w", "mod.wt.uw", "mod.eb.w", "mod.eb.uw", "mod.fg.w", "mod.fg.uw")
  
  # final product: full.mod, which is the weighted and unweighted modularity scores for the complete network
  
  
  # save complete network metric data
  # complete.deg.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  # complete.str.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  # complete.btw.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  # complete.clust.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  # complete.dens.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  # complete.iso.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  # complete.mod.name <- paste(<insert file naming system>, ".Rdata", sep = "")
  

  save(full.deg, file = complete.deg.name)
  save(full.str, file = complete.str.name)
  save(full.btw, file = complete.btw.name)
  save(full.clust, file = complete.clust.name)
  save(full.dens, file = complete.dens.name)
  save(full.iso, file = complete.iso.name)
  save(full.mod, file = complete.mod.name)
  
  ###### Calculate metrics for sampled networks ##########
  
  # create objects for the variations in sampling for which metrics will be calculated
  # these will be used in looping
  ind.sample <- c("q1m", "q15m", "q60m", "q3h", "q12h", "q24h", "q72h")
  pop.sample <- seq(100, 10, -10)
  contact.type <- c("space-time", "KDE UDOI")
  
  # create empty data frames to store results
  # uses dynamic allocation which is less efficient but functional for these purposes
  full_samp.deg <- NULL
  full_samp.str <- NULL
  full_samp.btw <- NULL
  full_samp.eig <- NULL
  full_samp.clust <- NULL
  full_samp.dens <- NULL
  full_samp.iso <- NULL
  
  
  # loop through each of the sampling variations
  for(q in 1:length(contact.type)){
    ct <- contact.type[q]
    for(r in 1:length(pop.sample)){
      ps <- pop.sample[r]
      for(s in 1:length(ind.sample)){
        is <- ind.sample[s]
        
        # take one sampling set edgelist at a time
        temp.el <- subset(full.el, full.el$ind.sample==is & full.el$pop.sample==ps & full.el$contact.type==ct)
        net<-graph.data.frame(temp.el,directed=F) # convert edgelist to igraph object
        
        
        ### DEGREE
        # create dataframe of degree 
        s.deg <- degree(net,v=V(net))
        s.deg.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
        colnames(s.deg.df) <- "id"
        s.deg.df$deg <- s.deg
        s.deg.df <- s.deg.df[order(s.deg.df$id),]
        
        # next, add in the individuals with degree 0 for this sampling level
        # read in list of observed individuals for this population sampling level
        observed.file <- paste(start.type, "/", h.type, "/", cont.type, "/", ps,"Samp/Sample IDs_", ps, "Sample_", "Simulation", sim, ".Rdata", sep="")
        observed <- get(load(observed.file))
        
        samp.deg <- data.frame(
          id=observed,
          deg=NA
        )

        # individuals that were not isolates (for the following loop)
        have.deg <- s.deg.df$id

        # loop to assign isolates a degree of 0
        for(l in 1:length(observed)){
          if(observed[l] %in% have.deg){
            samp.deg[l,2] <- s.deg.df[s.deg.df$id==observed[l], 2]
          }
          else{
            samp.deg[l,2] <- 0
          }
        }

        # assign sampling levels for tracking purposes
        samp.deg$ind.sample <- is
        samp.deg$pop.sample <- ps
        samp.deg$contact.type <- ct

        # save for next round
        full_samp.deg <- rbind(full_samp.deg, samp.deg)


        
        ### STRENGTH
        # repeat above steps for strength; weights must be in dataframe under "weights" or assigned separately here 
        s.str <- strength(net, vids = V(net))
        s.str.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
        colnames(s.str.df) <- "id"
        s.str.df$str <- s.str
        s.str.df <- s.str.df[order(s.str.df$id),]

        # read in list of observed individuals for this population sampling level
        samp.str <- data.frame(
          id=observed,
          str=NA
        )

        # list of individuals that were observed (for use in following loop)
        have.str <- s.str.df$id

        # loop to assign strength of 0 to isolates
        for(l in 1:length(observed)){
          if(observed[l] %in% have.str){
            samp.str[l,2] <- s.str.df[s.str.df$id==observed[l], 2]
          }
          else{
            samp.str[l,2] <- 0
          }
        }

        # assign sampling level information for tracking purposes
        samp.str$ind.sample <- is
        samp.str$pop.sample <- ps
        samp.str$contact.type <- ct

        # save for next round
        full_samp.str <- rbind(full_samp.str, samp.str)
        
        
        ### BETWEENNESS
        # repeat above steps for betweenness; can do weighted or unweighted
        s.btw <- betweenness(net, v = V(net), normalized = T, weights = NULL) # change here for weighted vs unweighted betweenness
        s.btw.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
        colnames(s.btw.df) <- "id"
        s.btw.df$btw <- s.btw
        s.btw.df <- s.btw.df[order(s.btw.df$id),]

        # read in list of observed individuals for this population sampling level
        samp.btw <- data.frame(
          id=observed,
          btw=NA
        )

        # list of individuals that were observed (for use in following loop)
        have.btw <- s.btw.df$id

        # loop to assign betweenness of 0 to isolates
        for(l in 1:length(observed)){
          if(observed[l] %in% have.btw){
            samp.btw[l,2] <- s.btw.df[s.btw.df$id==observed[l], 2]
          }
          else{
            samp.btw[l,2] <- 0
          }
        }

        # assign sampling level information for tracking purposes
        samp.btw$ind.sample <- is
        samp.btw$pop.sample <- ps
        samp.btw$contact.type <- ct

        # save for next round
        full_samp.btw <- rbind(full_samp.btw, samp.btw)

        
        ### TRANSITIVITY
        # repeat above steps for transitivity; can be weighted or unweighted
        s.clust <- transitivity(net, type = c("local"), vids = V(net), isolates = c("zero"), weights = NULL) # change here for weighted vs unweighted
        s.clust.df <- data.frame(as.numeric(unique(c(paste(temp.el$individual.1), paste(temp.el$individual.2)))))
        colnames(s.clust.df) <- "id"
        s.clust.df$clust <- s.clust
        s.clust.df <- s.clust.df[order(s.clust.df$id),]

        # read in list of observed individuals for this population sampling level
        samp.clust <- data.frame(
          id=observed,
          clust=NA
        )

        # list of individuals that were observed (for use in following loop)
        have.clust <- s.clust.df$id

        # loop to assign transitivity of 0 to isolates
        for(l in 1:length(observed)){
          if(observed[l] %in% have.clust){
            samp.clust[l,2] <- s.clust.df[s.clust.df$id==observed[l], 2]
          }
          else{
            samp.clust[l,2] <- 0
          }
        }

        # assign sampling level information for tracking purposes
        samp.clust$ind.sample <- is
        samp.clust$pop.sample <- ps
        samp.clust$contact.type <- ct

        # save for next round
        full_samp.clust <- rbind(full_samp.clust, samp.clust)
        
        
        ### DENSITY
        # calculate density by the number of edges out of all possible edges based on population sampling level
        samp.dens <- data.frame(nrow(temp.el)/choose(ps, 2))
        colnames(samp.dens) <- "dens"

        # assign sampling level information for tracking purposes
        samp.dens$ind.sample <- is
        samp.dens$pop.sample <- ps
        samp.dens$contact.type <- ct

        # save for next round
        full_samp.dens <- rbind(full_samp.dens, samp.dens)


        ### PROPORTION ISOLATES
        # save the proportion of the population that is isolates (individuals with degree of zero)
        samp.iso <- data.frame((nrow(subset(samp.deg, samp.deg$deg==0)))/ps)
        colnames(samp.iso) <- "iso"

        # assign sampling level information for tracking purposes
        samp.iso$ind.sample <- is
        samp.iso$pop.sample <- ps
        samp.iso$contact.type <- ct

        # save for next round
        full_samp.iso <- rbind(full_samp.iso, samp.iso)
        
        
        ### MODULARITY
        # calculate different weighted and unweighted metrics of modularity 
        if(gsize(net)>0){ #need to account for unconnected networks
          wtc.w <- cluster_walktrap(net, weights = E(net)$value, steps = 4, modularity = T)
          wtc.uw <- cluster_walktrap(net, weights = NULL, steps = 4, modularity = T)
          
          eb.w <- cluster_edge_betweenness(net, weights = E(net)$value, directed = F)
          eb.uw <- cluster_edge_betweenness(net, weights = NULL, directed = F)
          
          fg.w <- cluster_fast_greedy(net, weights = E(net)$value)
          fg.uw <- cluster_fast_greedy(net, weights = NULL)
          
          mod.wt.w <- modularity(net, membership(wtc.w))
          mod.wt.uw <- modularity(net, membership(wtc.uw))
          mod.eb.w <- modularity(net, membership(eb.w))
          mod.eb.uw <- modularity(net, membership(eb.uw))
          mod.fg.w <- modularity(net, membership(fg.w))
          mod.fg.uw <- modularity(net, membership(fg.uw))
          
          samp.mod <- c(mod.wt.w, mod.wt.uw, mod.eb.w, mod.eb.uw, mod.fg.w, mod.fg.uw)
        }else{
          samp.mod <- NA
        }
        samp.mod <- data.frame(matrix(samp.mod, ncol = 6))
        colnames(samp.mod) <- c("mod.wt.w", "mod.wt.uw", "mod.eb.w", "mod.eb.uw", "mod.fg.w", "mod.fg.uw")
        
        
        # assign sampling level information for tracking purposes
        samp.mod$ind.sample <- is
        samp.mod$pop.sample <- ps
        samp.mod$contact.type <- ct
        
        # save for next round
        full_samp.mod <- rbind(full_samp.mod, samp.mod)
      }
    }
  }
  
  
  # get rid of KDE data where ind.sample != q24h or q72h
  
  # DEGREE
  full_samp.deg <- fix.KDE(full_samp.deg)

  # STRENGTH
  full_samp.str <- fix.KDE(full_samp.str)

  # BETWEENNESS
  full_samp.btw <- fix.KDE(full_samp.btw)

  # TRANSITIVITY
  full_samp.clust <- fix.KDE(full_samp.clust)

  # DENSITY
  full_samp.dens <- fix.KDE(full_samp.dens)

  # PROPORTION ISOLATES
  full_samp.iso <- fix.KDE(full_samp.iso)
  
  # MODULARITY
  full_samp.mod <- fix.KDE(full_samp.mod)
  
  
  # save sample network metric data
  # sample.deg.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # sample.str.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # sample.btw.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # sample.clust.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # sample.dens.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # sample.iso.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  # sample.mod.name <- paste(<insert naming structure>, ".Rdata", sep = "")
  
  
  save(full_samp.deg, file = sample.deg.name)
  save(full_samp.str, file = sample.str.name)
  save(full_samp.btw, file = sample.btw.name)
  save(full_samp.clust, file = sample.clust.name)
  save(full_samp.dens, file = sample.dens.name)
  save(full_samp.iso, file = sample.iso.name)
  save(full_samp.mod, file = sample.mod.name)
  
}
#---------------- End of primary loop -------------#

beep(8)
