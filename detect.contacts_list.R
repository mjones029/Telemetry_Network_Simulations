## detect.contacts_list
# 
#========================================================	
# ---
# ## title: Detect contacts in a list of sampled trajectories
# author: Marie Gilbertson
# date: "12/19/2018"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for detecting space-time contacts within a dataframe of simulated (sampled) trajectories
# Note: This is a different function than that used for detecting contacts for the complete network 
# Requires simul.contacts.R and distance.R functions

# full.dataset = dataset of sampled trajectories
# dist.thresh = distance threshold for contact definition

detect.contacts_list <- function(full.dataset=full.dataset, dist.thresh=dist.thresh){
  
  # list of sampling intervals for this project
  ind.sample <- c("q15m", "q60m", "q3h", "q12h", "q24h", "q72h")
  
  # create empty dataframe for storing final results
  all.contacts <- NULL # note: uses dynamic allocation which is less than ideal for streamlining, but manages differing numbers of contacts detected 
  
  for(l in 1:length(ind.sample)){
    dataset <- full.dataset[[l]]
    print(paste(ind.sample[l],"sample"))
    # set up empty dataframes for storing results
    contacts <- NULL
    id <- unique(dataset$id)
    
    # loop through whole population for pairwise contact detection
    for(i in 1:(length(id)-1)){
      print(i)
      ind1 <- id[i]
      
      for(j in (i+1):length(id)){  
        #print(j)
        ind2 <- id[j]
        
        round.contacts <- simul.contacts(dataset, ind1, ind2) 
        contacts <- rbind(contacts, round.contacts)
      }
    }
    # count number of contacts per dyad
    edge.list <- count(contacts, c('id','id.2'))
    colnames(edge.list) <- c("puma.1", "puma.2", "value")
    if(nrow(edge.list)==0){edge.list[1,] <- NA}
    # label the edgelist by the individual-level sampling effort
    edge.list$ind.sample <- ind.sample[l]
    # save for next round of loop
    all.contacts <- rbind(all.contacts, edge.list)
  }
  return(all.contacts)
}
