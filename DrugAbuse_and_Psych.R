#Author: Arvis Sulovari, PhD


##1) What psychiatric disorders are developed as consequences of substance abuse? 
##2) What genetic mechanisms lead to drug abuse, followed by comorbid psychiatric disorders

#Discovery datasets
ssadda <- read.delim("C:/Users/arvis/[...].txt")
ssadda_AA <- ssadda[which(ssadda$A8a_RACE==4 | ssadda$A8a_RACE==5),]
ssadda_EA <- ssadda[which(ssadda$A8a_RACE==6 | ssadda$A8a_RACE==7),]

#Replication datasets
ssadda_rep_AA <- read.delim("C:/Users/arvis/[...].txt")
ssadda_rep_EA <- read.delim("C:/Users/arvis/[...].txt")
ssadda_rep <- rbind(ssadda_rep_EA,ssadda_rep_AA)
ssadda_rep_AA <- ssadda_rep[which(ssadda_rep$A8a_RACE==4 | ssadda_rep$A8a_RACE==5),]
ssadda_rep_EA <- ssadda_rep[which(ssadda_rep$A8a_RACE==6 | ssadda_rep$A8a_RACE==7),] 


#########################################################################################
#                                                                                       #
#                          INDIVIDUAL-LEVEL REPLICATIONS                                #
#                                                                                       #
#########################################################################################

#Obtain index for Age-of-Onset columns in each of the four datasets
ao_col_indeces_AA <- grep("AgeO$",colnames(ssadda_AA))
ao_col_indeces_rep_AA <- grep("AgeO$",colnames(ssadda_rep_AA))
ao_col_indeces_EA <- grep("AgeO$",colnames(ssadda_EA))
ao_col_indeces_rep_EA <- grep("AgeO$",colnames(ssadda_rep_EA))

####Clean-up non-age numbers fro mthe age of onset columns:
for(i in 1:length(ao_col_indeces_EA)){ssadda[which(as.numeric(ssadda[,ao_col_indeces_EA[i]])>70 | as.numeric(ssadda[,ao_col_indeces_EA[i]])==0),ao_col_indeces_EA[i]] <- c("UC")}
for(i in 1:length(ao_col_indeces_rep_EA)){ssadda_rep[which(as.numeric(ssadda_rep[,ao_col_indeces_rep_EA[i]])>70 | as.numeric(ssadda_rep[,ao_col_indeces_rep_EA[i]])==0),ao_col_indeces_rep_EA[i]] <- c("UC")}
for(i in 1:length(ao_col_indeces_AA)){ssadda_AA[which(as.numeric(ssadda_AA[,ao_col_indeces_AA[i]])>70 | as.numeric(ssadda_AA[,ao_col_indeces_AA[i]])==0),ao_col_indeces_AA[i]] <- c("UC")}
for(i in 1:length(ao_col_indeces_EA)){ssadda_EA[which(as.numeric(ssadda_EA[,ao_col_indeces_EA[i]])>70 | as.numeric(ssadda_EA[,ao_col_indeces_EA[i]])==0),ao_col_indeces_EA[i]] <- c("UC")}
for(i in 1:length(ao_col_indeces_rep_AA)){ssadda_rep_AA[which(as.numeric(ssadda_rep_AA[,ao_col_indeces_rep_AA[i]])>70 | as.numeric(ssadda_rep_AA[,ao_col_indeces_rep_AA[i]])==0),ao_col_indeces_rep_AA[i]] <- c("UC")}
for(i in 1:length(ao_col_indeces_rep_EA)){ssadda_rep_EA[which(as.numeric(ssadda_rep_EA[,ao_col_indeces_rep_EA[i]])>70 | as.numeric(ssadda_rep_EA[,ao_col_indeces_rep_EA[i]])==0),ao_col_indeces_rep_EA[i]] <- c("UC")}


#Function for extraction of age-of-onset variables
AO_order_fun <- function(ssadda_in=ssadda,out_mx=outmx_ao_order){
  
  #General pattern-finder for order in which Age of Onset occurred
  ao_col_indeces <- grep("AgeO$",colnames(ssadda_in))
  ssadda_ageOfOnset <- ssadda_in[,ao_col_indeces]
  
  ##Output matrix
  out_mx <- array(NA,dim=c(nrow(ssadda_ageOfOnset),length(ao_col_indeces)))
  colnames(out_mx) <- colnames(ssadda_ageOfOnset)
  
  #Populate output matrix
  for(i in 1:nrow(ssadda_ageOfOnset)) {
    #Vectorize instead of inner for-loop
    out_mx[i,] <- match(colnames(ssadda_ageOfOnset),colnames(ssadda_ageOfOnset)[order(ssadda_ageOfOnset[i,])])
    out_mx[i,which(is.na(ssadda_ageOfOnset[i,]))] <- 999
    
    if(i%%100==0){
      per_comp <- round((i/nrow(ssadda_ageOfOnset))*100,3)
      print(paste0(per_comp," %"))
    }
    else{}
  }
  
  #Order of events
  mx_overall_order <- array(NA,dim=c(135,3))
  mx_overall_order[,1] <- colnames(ssadda_ageOfOnset)
  
  for(i in 1:135){
    mx_overall_order[i,2] <- mean(out_mx[(which(out_mx[,i]<999)),i])
    #print(m)
  }
  
  #Count number of samples that support every 'direct & unidirectional' link
  ordered_out_mx <- mx_overall_order[order(as.numeric(mx_overall_order[,2])),]
  
  for(i in 1:134){
    
    factor_one <- ordered_out_mx[i,1]
    index_one <- which(colnames(ssadda_ageOfOnset)==factor_one)
    
    factor_two <- ordered_out_mx[(i+1),1]
    index_two <- which(colnames(ssadda_ageOfOnset)==factor_two)
    
    #ordered_out_mx[i,3] <- length(which(out_mx[,index_one]==i & out_mx[,index_two]==(i+1)))
    ordered_out_mx[i,3] <- length(which(out_mx[,index_one] < out_mx[,index_two]))
  }
  
  #colnames(ordered_out_mx) <- c("","","")
  return(ordered_out_mx)
}


#Run AO_order_fun() to replicate ordered events
ALL_ordered_events <- AO_order_fun(ssadda_in = ssadda)
AA_ordered_events <- AO_order_fun(ssadda_in = ssadda_AA)
EA_ordered_events <- AO_order_fun(ssadda_in = ssadda_EA)

##Reformat 
ALL_ordered_events_v2 <- cbind(ALL_ordered_events[-135,1],ALL_ordered_events[-1,1],ALL_ordered_events[-135,2],
                            ALL_ordered_events[-135,3],paste0(ALL_ordered_events[-135,1],ALL_ordered_events[-1,1]))

AA_ordered_events_v2 <- cbind(AA_ordered_events[-135,1],AA_ordered_events[-1,1],AA_ordered_events[-135,2],AA_ordered_events[-135,3],
                           paste0(AA_ordered_events[-135,1],AA_ordered_events[-1,1]))

EA_ordered_events_v2 <- cbind(EA_ordered_events[-135,1],EA_ordered_events[-1,1],EA_ordered_events[-135,2],EA_ordered_events[-135,3],
                           paste0(EA_ordered_events[-135,1],EA_ordered_events[-1,1]))

##
ALL_rep_ordered_events <- AO_order_fun(ssadda_in = ssadda_rep)
AA_rep_ordered_events <- AO_order_fun(ssadda_in = ssadda_rep_AA)
EA_rep_ordered_events <- AO_order_fun(ssadda_in = ssadda_rep_EA)

##Reformat
ALL_rep_ordered_events_v2 <- cbind(ALL_rep_ordered_events[-135,1],ALL_rep_ordered_events[-1,1],ALL_rep_ordered_events[-135,2],ALL_rep_ordered_events[-135,3],
                                paste0(ALL_rep_ordered_events[-135,1],ALL_rep_ordered_events[-1,1]))

AA_rep_ordered_events_v2 <- cbind(AA_rep_ordered_events[-135,1],AA_rep_ordered_events[-1,1],AA_rep_ordered_events[-135,2],AA_rep_ordered_events[-135,3],
                               paste0(AA_rep_ordered_events[-135,1],AA_rep_ordered_events[-1,1]))

EA_rep_ordered_events_v2 <- cbind(EA_rep_ordered_events[-135,1],EA_rep_ordered_events[-1,1],EA_rep_ordered_events[-135,2],EA_rep_ordered_events[-135,3],
                               paste0(EA_rep_ordered_events[-135,1],EA_rep_ordered_events[-1,1]))




###
#FIND THE OVERLP
##Allow for "skipping" when calculating the overlap

##NON-NETWORK SOLUTION
AA_ordered_events_v2[,1]
AA_rep_ordered_events_v2[,1]

EA_ordered_events_v2[,1]
EA_rep_ordered_events_v2[,1]

#Check if any of the paths (A->B) are replicated in all 4 datasets, allowing for "skipping"
replicated_pairs <- array(NA,dim=c(300,4))

#ONLY for the first (out of 4) loopings
counter <- 0

for(k in 1:3) {
  
  tryCatch({
    
    for(i in 1:133) {
      for(j in (i+1)) {
        
        if(k==1){
              current_pair <- c(AA_ordered_events_v2[i,1],AA_ordered_events_v2[j,1])
              current_ordered_events <- AA_ordered_events_v2
              thresh <- nrow(ssadda_AA)*0.2
        }
        
        else if(k==2) {
              current_pair <- c(AA_rep_ordered_events_v2[i,1],AA_rep_ordered_events_v2[j,1])
              current_ordered_events <- AA_rep_ordered_events_v2
              thresh <- nrow(ssadda_rep_AA)*0.2
        }
        
        else if(k==3) {
              current_pair <- c(EA_ordered_events_v2[i,1],EA_ordered_events_v2[j,1])
              current_ordered_events <- EA_ordered_events_v2
              thresh <- nrow(ssadda_EA)*0.2
        }
        
        else if(k==4) {
              current_pair <- c(EA_rep_ordered_events_v2[i,1],EA_rep_ordered_events_v2[j,1])
              current_ordered_events <- EA_rep_ordered_events_v2
              thresh <- nrow(ssadda_rep_EA)*0.2
        }
        
        
        current_ordered_events <- current_ordered_events[which(as.numeric(as.character(current_ordered_events[,4]))>=thresh),]
        
        if(length(which(AA_ordered_events_v2[,1]==current_pair[1])) + 
            length(which(AA_rep_ordered_events_v2[,1]==current_pair[1])) + 
            length(which(EA_ordered_events_v2[,1]==current_pair[1])) + 
            length(which(EA_rep_ordered_events_v2[,1]==current_pair[1])) +
            
            length(which(AA_ordered_events_v2[,1]==current_pair[2])) + 
            length(which(AA_rep_ordered_events_v2[,1]==current_pair[2])) + 
            length(which(EA_ordered_events_v2[,1]==current_pair[2])) + 
            length(which(EA_rep_ordered_events_v2[,1]==current_pair[2])) == 8) {
              
            aa_status <- (which(AA_ordered_events_v2[,1]==current_pair[1]) < 
                              which(AA_ordered_events_v2[,1]==current_pair[2]))
            
            aa_rep_status <- (which(AA_rep_ordered_events_v2[,1]==current_pair[1]) < 
                                which(AA_rep_ordered_events_v2[,1]==current_pair[2]))
            
            ea_status <- (which(EA_ordered_events[,1]==current_pair[1]) < 
                                which(EA_ordered_events[,1]==current_pair[2]))
            
            ea_rep_status <- (which(EA_ordered_events_v2[,1]==current_pair[1]) < 
                            which(EA_ordered_events_v2[,1]==current_pair[2]))
          
            if(aa_status==T & aa_rep_status==T & ea_status == T & ea_rep_status == T){
              print(c(i,j))
              counter <- counter+1
              
              replicated_pairs[counter,1] <- current_ordered_events[i,1]
              replicated_pairs[counter,2] <- current_ordered_events[j,1]
              
              replicated_pairs[counter,3] <- which((colnames(ssadda_ageOfOnset))==current_ordered_events[i,1])
              replicated_pairs[counter,4] <- which((colnames(ssadda_ageOfOnset))==current_ordered_events[j,1])
              
            }
        }
        
        }
      
    }
    
  }, error=function(e){})
  
}   

replicated_pairs_vClean <- (cbind(replicated_pairs[,1],replicated_pairs[,2],as.numeric(replicated_pairs[,3]),as.numeric(replicated_pairs[,4])))

write.csv(replicated_pairs,"Final&clean/replicated_pairs_FINAL_dec9.csv")



#########################################################################################
#                                                                                       #
#                                 GRAPH-BASED SOLUTIONS                                 #
#                                                                                       #
#########################################################################################

#Need to use Graph Theory.
##Build adjacency matrix with info on each network
require(igraph)
require(graph)
require(networkD3)


AA_adj_mx <- matrix(0,nr=135,nc=135)
#colnames(ssadda_ageOfOnset)
combined_mx_num <- cbind(match(AA_ordered_events_v2[,1],colnames(ssadda_ageOfOnset)),match(AA_ordered_events_v2[,2],colnames(ssadda_ageOfOnset)))

combined_mx_num_ALL <- cbind(c(match(AA_ordered_events_v2[,1],colnames(ssadda_ageOfOnset)),match(AA_rep_ordered_events_v2[,1],colnames(ssadda_ageOfOnset)),
                                   match(EA_ordered_events_v2[,1],colnames(ssadda_ageOfOnset)),match(EA_rep_ordered_events_v2[,1],colnames(ssadda_ageOfOnset))),
                                   c(match(AA_ordered_events_v2[,2],colnames(ssadda_ageOfOnset)),match(AA_rep_ordered_events_v2[,2],colnames(ssadda_ageOfOnset)),
                                   match(EA_ordered_events_v2[,2],colnames(ssadda_ageOfOnset)),match(EA_rep_ordered_events_v2[,2],colnames(ssadda_ageOfOnset)))
                             )


for(i in 1:nrow(combined_mx_num_ALL)) {
    #(combined_mx_num[,1]==i & combined_mx_num[,2]==j)
    AA_adj_mx[combined_mx_num_ALL[i,1],combined_mx_num_ALL[i,2]] <- 1
  }


#Plot graph of the adacency matrix
#adj.mat <- matrix(sample(c(0,1), 9, replace=TRUE), nr=3)
g <- graph.adjacency(AA_adj_mx)
plot(g)





#BEST SOLUTION
network_data <- data.frame(combined_mx_num)
network_data_AA <- data.frame(cbind(AA_ordered_events_v2[,1],AA_ordered_events_v2[,2]))
network_data_EA <- data.frame(cbind(EA_ordered_events_v2[,1],EA_ordered_events_v2[,2]))

network_data_AA_combined <- data.frame(rbind(cbind(AA_ordered_events_v2[,1],AA_ordered_events_v2[,2]),cbind(AA_rep_ordered_events_v2[,1],AA_rep_ordered_events_v2[,2])))

simpleNetwork(network_data_AA,fontSize = 10)
simpleNetwork(network_data_EA,fontSize = 10)
simpleNetwork(network_data_AA_combined,fontSize = 10)

#More complex networks
##AA
AA_network_data_forced <- data.frame(rbind(cbind(AA_ordered_events_v2[,1],AA_ordered_events_v2[,2],AA_ordered_events_v2[,4]),
                                             cbind(AA_rep_ordered_events_v2[,1],AA_rep_ordered_events_v2[,2],AA_rep_ordered_events_v2[,4])))

AA_links_data <- data.frame(cbind(match(AA_network_data_forced$X1,colnames(ssadda_ageOfOnset)),
      match(AA_network_data_forced$X2,colnames(ssadda_ageOfOnset)),
      AA_network_data_forced$X3))
colnames(AA_links_data) <- c("source","target","value")

#RUN ONLY ONCE:
#AA_links_data$source <- AA_links_data$source-1
#AA_links_data$target <- AA_links_data$target-1

AA_nodes_data <- data.frame(cbind(colnames(ssadda_ageOfOnset),c(rep(5,81),rep(2,54)),rep(20,135)))
colnames(AA_nodes_data) <- c("name","group","size")

##AA
forceNetwork(Links = AA_links_data,Nodes = AA_nodes_data,Source = "source",Target = "target",Value = "value",NodeID = "name",Group = "group",Nodesize = "size",
             fontSize = 20,charge = -400,zoom=T)


##EA
EA_network_data_forced <- data.frame(rbind(cbind(EA_ordered_events_v2[,1],EA_ordered_events_v2[,2],EA_ordered_events_v2[,4]),
                                           cbind(EA_rep_ordered_events_v2[,1],EA_rep_ordered_events_v2[,2],EA_rep_ordered_events_v2[,4])))

EA_links_data <- data.frame(cbind(match(EA_network_data_forced$X1,colnames(ssadda_ageOfOnset)),
                                  match(EA_network_data_forced$X2,colnames(ssadda_ageOfOnset)),
                                  EA_network_data_forced$X3))
colnames(EA_links_data) <- c("source","target","value")

#RUN ONLY ONCE:
#EA_links_data$source <- EA_links_data$source-1
#EA_links_data$target <- EA_links_data$target-1

EA_nodes_data <- data.frame(cbind(colnames(ssadda_ageOfOnset),c(rep(5,81),rep(2,54)),rep(20,135)))
colnames(EA_nodes_data) <- c("name","group","size")

##EA
forceNetwork(Links = EA_links_data,Nodes = EA_nodes_data,Source = "source",Target = "target",Value = "value",NodeID = "name",Group = "group",Nodesize = "size",
             fontSize = 20,charge = -400,zoom=T)



##AA+EA
BOTH_network_data_forced <- data.frame(rbind(cbind(EA_ordered_events_v2[,1],EA_ordered_events_v2[,2],EA_ordered_events_v2[,4]),
                                             cbind(EA_rep_ordered_events_v2[,1],EA_rep_ordered_events_v2[,2],EA_rep_ordered_events_v2[,4]),
                                             cbind(AA_ordered_events_v2[,1],AA_ordered_events_v2[,2],AA_ordered_events_v2[,4]),
                                             cbind(AA_rep_ordered_events_v2[,1],AA_rep_ordered_events_v2[,2],AA_rep_ordered_events_v2[,4])))

BOTH_links_data <- data.frame(cbind(match(BOTH_network_data_forced$X1,colnames(ssadda_ageOfOnset)),
                                    match(BOTH_network_data_forced$X2,colnames(ssadda_ageOfOnset)),
                                    BOTH_network_data_forced$X3))
colnames(BOTH_links_data) <- c("target","source","value")

#RUN ONLY ONCE:
#BOTH_links_data$source <- BOTH_links_data$source-1
#BOTH_links_data$target <- BOTH_links_data$target-1

BOTH_nodes_data <- data.frame(cbind(colnames(ssadda_ageOfOnset),c(rep(5,81),rep(2,54)),rep(4,135)))
colnames(BOTH_nodes_data) <- c("name","group","size")

##BOTH
forceNetwork(Links = BOTH_links_data,Nodes = BOTH_nodes_data,Source = "source",Target = "target",Value = "value",NodeID = "name",Group = "group",Nodesize = "size",
             fontSize = 40,charge = -1000,zoom=T)



#TEST NETWORK
BOTH_links_data_tmp <- BOTH_links_data
BOTH_nodes_data_tmp <- BOTH_nodes_data
BOTH_nodes_data_tmp$size <- as.numeric(as.character(BOTH_nodes_data_tmp$size))+11

#IMPORTANT!!! Revert node IDs to +1 (not zero indexed)
##RUN ONLY ONCE!
#BOTH_links_data_tmp$target <- BOTH_links_data_tmp$target+1
#BOTH_links_data_tmp$source <- BOTH_links_data_tmp$source+1

replicated_nodes_clean <- unique(cbind(na.omit(as.numeric(replicated_pairs[,3])),na.omit(as.numeric(replicated_pairs[,4]))))
replicated_nodes_clean_v2 <- cbind(replicated_nodes_clean,rep(100,123))
replicated_nodes_clean_v2 <- as.data.frame(replicated_nodes_clean_v2)
colnames(replicated_nodes_clean_v2) <- colnames(BOTH_links_data_tmp)

#Zero Index the Links df
#RUN ONCE ONLY
replicated_nodes_clean_v2$target <- replicated_nodes_clean_v2$target-1
replicated_nodes_clean_v2$source <- replicated_nodes_clean_v2$source-1


#
#for(i in 1:123){
#  indeks <- which(BOTH_links_data_tmp[,1]==replicated_nodes_clean[i,1] & BOTH_links_data_tmp[,2]==replicated_nodes_clean[i,2])
  
#  BOTH_links_data_tmp[indeks,3] <- 100 
#}

#BOTH_links_data_tmp[which(BOTH_links_data_tmp[,3]!=100),3] <- rep(0,184)

      
forceNetwork(Links = replicated_nodes_clean_v2,Nodes = BOTH_nodes_data_tmp,Source = "source",Target = "target",Value = "value",NodeID = "name",Group = "group",Nodesize = "size",
             fontSize = 40,charge = -200,zoom=T)


sankeyNetwork(Links = replicated_nodes_clean_v2,Nodes = BOTH_nodes_data_tmp,Source = "source",Target = "target",Value = "value",NodeID = "name",fontSize = 40)



###iGraph version of the D3 graph

all_adj_mx <- array(0,dim=c(122,122))
colnames(all_adj_mx) <- unique(c(replicated_nodes_clean[,1],replicated_nodes_clean[,2]))
rownames(all_adj_mx) <- unique(c(replicated_nodes_clean[,1],replicated_nodes_clean[,2]))

#Populate adjacency matrix
for(i in 1:123) {
  
  row_n <- which(rownames(all_adj_mx)==replicated_nodes_clean[i,1])
  col_n <- which(colnames(all_adj_mx)==replicated_nodes_clean[i,2])
  
  all_adj_mx[row_n,col_n] <- 1
  
}


#Cleaning-up
ssadda_ao_replicated <- ssadda_ageOfOnset[,as.numeric(as.character(unique(colnames(all_adj_mx))))]


#Re-name columns and rows of adj matrix according to actual SSADDA header name
colnames(all_adj_mx) <- colnames(ssadda_ageOfOnset)[as.numeric(as.character(colnames(all_adj_mx)))]
rownames(all_adj_mx) <- colnames(ssadda_ageOfOnset)[as.numeric(as.character(rownames(all_adj_mx)))]


#CLEAN-UP the age of onset mx
for(i in 1:ncol(ssadda_ao_replicated)){p <- which(ssadda_ao_replicated[,i]>80 | ssadda_ao_replicated[,i]<1); ssadda_ao_replicated[p,i] <- "NA"}

avg_age_arr <- array(NA,dim=c(ncol(ssadda_ao_replicated),2))
avg_age_arr[,1] <- colnames(ssadda_ao_replicated)

for(i in 1:ncol(ssadda_ao_replicated)) {
  
  avg_val <- mean(na.omit(as.numeric(as.character(ssadda_ao_replicated[,i]))))
  avg_age_arr[i,2] <- avg_val
  
}

#Order chronologically the age of onset events
avg_age_arr <- avg_age_arr[order(as.numeric(as.character(avg_age_arr[,2]))),]

##
#Remove E connecting V that are in the wrong order (according to avg_age_arr[,2])
#Working version Below!
##IMPORTANT: BEFORE RUNNING CODE BELOW, colnames(all-adj_mx) MUST BE FULL CHARACTER_TYPE NAMES
for(i in 1:ncol(all_adj_mx)) {
  
  running_node <- colnames(all_adj_mx)[i]
  next_node <- colnames(all_adj_mx)[(which(all_adj_mx[i,]==1))]
  origin_node <- colnames(all_adj_mx)[(which(all_adj_mx[,i]==1))]
  
  #Running node is the origin
  if(length(next_node)!=0){
    #print(paste0("From: ",running_node," To: ",next_node,";"))
    for(j in 1:length(next_node)) {
      if((which(avg_age_arr[,1]==toString(next_node[j]))) < (which(avg_age_arr[,1]==toString(running_node)))){
        all_adj_mx[running_node,next_node] <- 0
      }
    }
  }
  else {}
  
  
  #Running node is the destination
  if(length(origin_node)!=0){
    #print(paste0("From: ",origin_node," To: ",running_node,";"))
    for(k in 1:length(origin_node)){
      if((which(avg_age_arr[,1]==toString(origin_node[k]))) > (which(avg_age_arr[,1]==toString(running_node)))){
        all_adj_mx[running_node,next_node] <- 0
      }
    }
  }
  else {}
}


View(all_adj_mx)

###

#Add an extra column for the order in which the vertices occur (average age of onset)
renaming_mx <- cbind(colnames(all_adj_mx),match(colnames(all_adj_mx),avg_age_arr[,1]))
gsub("AgeO","",paste0(renaming_mx[,1],"_",renaming_mx[,2]))

colnames(all_adj_mx) <- gsub("AgeO","",paste0(renaming_mx[,1],"_",renaming_mx[,2]))
rownames(all_adj_mx) <- gsub("AgeO","",paste0(renaming_mx[,1],"_",renaming_mx[,2]))

#unique(c(replicated_nodes_clean[,1],replicated_nodes_clean[,2]))
#rownames(all_adj_mx) <- unique(c(replicated_nodes_clean[,1],replicated_nodes_clean[,2]))
g <- graph.adjacency(all_adj_mx)
g <- graph.adjacency(ALL_ADJ)
plot(g,layout=layout_components)
tkplot(g)
tkigraph()

############################END OF D3 NETWORKS##############################

#########################GGPLOT2 Networks(requires all_adj_mx from above)###################################
install.packages("network","sna")
install.packages("ggnetwork")
install.packages("ggrepel")
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(intergraph)
library(ggnetwork)

#Random network
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)
network.vertex.names(net) = letters[1:10]

ggnet2(net,mode="circle")

#Node colors
net %v% "phono" = ifelse(letters[1:10] %in% c("a", "e", "i"), "vowel", "consonant")
ggnet2(net, color = "phono",mode="circle",size="degree",label=1:10,directed=T) +
  theme(panel.background = element_rect(fill = "grey90"))

#Layout vriable in ggnet2 will accept 2 columns with coordinates data. It'll have as many rows as the number of nodes(i.e.122).

mynet = network(all_adj_mx,directed=T)

mycoords <- as.data.frame(cbind(as.numeric(avg_age_arr[,2]),rep(1,122)))
ggnet2(mynet,layout.par = mycoords,arrow.size = 12,arrow.gap = 0.01,edge.size = 1,edge.color = "black",label = 1:123,label.size = 12,label.color = "black")


########################End of GGPLOT2 Networks##############################



##################################ALLUVIAL NETWORK##################################
require(alluvial)
#Use ssadda_ageOfOnset
alluvial_arr <- array(NA,dim=c(13000,135))


indeks_arr <- (replicated_nodes_clean_v2$target+1)
out_arr <- array(NA,dim=c(1,2))
for(i in 1:122) {
  
  arr <- which(ssadda_ageOfOnset[,indeks_arr[i]] <= ssadda_ageOfOnset[,indeks_arr[i+1]])
  slice_arr <- cbind(rep(i,length(arr)),arr)
  out_arr <- rbind(slice_arr,out_arr)
  
  #plot(x = slice_arr[,1],y = slice_arr[,2])
  
  #out_arr <- intersect(arr,out_arr)
  #print(length(out_arr))
  
  #which(ssadda_ageOfOnset[,indeks_arr[i+1]] < ssadda_ageOfOnset[,indeks_arr[i+2]])
  #which(ssadda_ageOfOnset[,indeks_arr[i+2]] < ssadda_ageOfOnset[,indeks_arr[i+3]])
}



##############################END OF ALLUVIAL NETWORK###############################


#########################################################################################
#                                                                                       #
#                                 PAIRWISE CORRELATIONS                                 #
#                                                                                       #
#########################################################################################

###Simple merge
merge(AA_ordered_events_v2,AA_rep_ordered_events_v2,by="V5")
merge(ALL_ordered_events_v2,ALL_rep_ordered_events_v2,by="V5")


#Write into CSV files
write.csv(ALL_ordered_events,"ALL_ordered_events.csv",row.names = F)
write.csv(AA_ordered_events,"AA_ordered_events.csv",row.names = F)
write.csv(EA_ordered_events,"EA_ordered_events.csv",row.names = F)



#Define Function that takes ssadda file and does pairwise correlations for all Age of Onset phenotypes

pairwise_corr_fun <- function(ssadda_input=ssadda,pairwise_cor_output="Pairwise_summary",direction="d2p") {
  
  #Force correct variable type
  pairwise_cor_output <- as.character(pairwise_cor_output)
  direction <- as.character(direction)
  
  #Let's pick up all "Age of Onset" columns
  ssadda_ageOfOnset <- ssadda_input[,grep("AgeO$",colnames(ssadda_input))]
  ageO_n <- ncol(ssadda_ageOfOnset)
  #colnames(ssadda_ageOfOnset)
  
  #Run pairwise correlations (around 31K of them)
  summary_array <- array(NA,dim=c(choose(ageO_n,2),7))
  
  #Last column header with drug information
  lim_col_n <- max(as.numeric(grep("H2",colnames(ssadda_ageOfOnset),ignore.case = F)))  
  
  #Nested loop for pairwise correlations
  k <- 0
  for(i in 1:(ageO_n-1)) {
    for(j in (i+1):ageO_n) {
    
    #Start incrementing index
    k <- k+1
    
    ### ==> Here insert code to determine if the Age of onset (AO) of column i is < AO(j) <==
    ### calculate correlaiton test for those rows separately for AO(i) > AO(j)
    
    if(direction=="both") {
      #Find row indeces of missing or non-age data
      arr_1 <- which(is.na(ssadda_ageOfOnset[,i]) | ssadda_ageOfOnset[,i]==0 | ssadda_ageOfOnset[,i] >80)
      arr_2 <- which(is.na(ssadda_ageOfOnset[,j]) | ssadda_ageOfOnset[,j]==0 | ssadda_ageOfOnset[,j] >80)
      arr_all <- unique(c(arr_1,arr_2))
      missing_per <- ((length(arr_all)/nrow(ssadda_ageOfOnset))*100)
      
    }
    
    else if(direction=="d2p") {
      #Find row indeces of data with AO(i) < AO(j) or (!AO(i) >= AO(j))
      not_d2p <- which(ssadda_ageOfOnset[,i]>=ssadda_ageOfOnset[,j])
      arr_1 <- which(is.na(ssadda_ageOfOnset[,i]) | ssadda_ageOfOnset[,i]==0 | ssadda_ageOfOnset[,i] >80)
      arr_2 <- which(is.na(ssadda_ageOfOnset[,j]) | ssadda_ageOfOnset[,j]==0 | ssadda_ageOfOnset[,j] >80)
      arr_all <- unique(c(arr_1,arr_2,not_d2p))
      missing_per <- ((length(arr_all)/nrow(ssadda_ageOfOnset))*100)
    }
    
    else if(direction=="p2d") {
      #Find row indeces of data with AO(i) > AO(j) or (!AO(i) <= AO(j))
      not_p2d <- which(ssadda_ageOfOnset[,i]<=ssadda_ageOfOnset[,j])
      arr_1 <- which(is.na(ssadda_ageOfOnset[,i]) | ssadda_ageOfOnset[,i]==0 | ssadda_ageOfOnset[,i] >80)
      arr_2 <- which(is.na(ssadda_ageOfOnset[,j]) | ssadda_ageOfOnset[,j]==0 | ssadda_ageOfOnset[,j] >80)
      arr_all <- unique(c(arr_1,arr_2,not_p2d))
      missing_per <- ((length(arr_all)/nrow(ssadda_ageOfOnset))*100)
    }
    
    else{
      print("Unknown value for direction variable")
    }
    
    
    if(missing_per>=99) next
    
    test_sum <- cor.test(ssadda_ageOfOnset[-c(arr_all),i],ssadda_ageOfOnset[-c(arr_all),j])
    
    summary_array[k,1] <- colnames(ssadda_ageOfOnset)[i]
    summary_array[k,2] <- colnames(ssadda_ageOfOnset)[j]
    summary_array[k,3] <- test_sum$estimate
    summary_array[k,4] <- test_sum$p.value
    summary_array[k,5] <- 100-missing_per
    summary_array[k,6] <- ifelse(i<=lim_col_n,c("Drugs"),c("Other"))
    summary_array[k,7] <- ifelse(j>lim_col_n,c("Psych"),c("Other"))
  
    }
    per_comp <- round((k/choose(ageO_n,2)*100),3)
    print(paste0(per_comp," %"))
  }
  
  #Label summary array columns
  colnames(summary_array) <- c("Age_of_onset_1","Age_of_onset_2","Correlation_coefficient","Correlation_Pvalue","Available_data_%","Drugs_Age_of_onset","Psych_Age_of_onset")
  write.csv(summary_array,paste0(pairwise_cor_output,".csv"))
}


####End of Function



#Now run the function 'pairwise_corr_fun' on any ssadda files (with both)!
pairwise_corr_fun(ssadda_input = ssadda,pairwise_cor_output = "Pairwise_summary",direction="both")
pairwise_corr_fun(ssadda_input = ssadda_EA,pairwise_cor_output = "Pairwise_summary_EA",direction="both")
pairwise_corr_fun(ssadda_input = ssadda_AA,pairwise_cor_output = "Pairwise_summary_AA",direction="both")
pairwise_corr_fun(ssadda_input = ssadda_rep,pairwise_cor_output = "Pairwise_summary_rep",direction="both")
pairwise_corr_fun(ssadda_input = ssadda_rep_EA,pairwise_cor_output = "Pairwise_summary_rep_EA",direction="both")
pairwise_corr_fun(ssadda_input = ssadda_rep_AA,pairwise_cor_output = "Pairwise_summary_rep_AA",direction="both")


#D2P
pairwise_corr_fun(ssadda_input = ssadda,pairwise_cor_output = "d2p_Pairwise_summary",direction="d2p")
pairwise_corr_fun(ssadda_input = ssadda_EA,pairwise_cor_output = "d2p_Pairwise_summary_EA",direction="d2p")
pairwise_corr_fun(ssadda_input = ssadda_AA,pairwise_cor_output = "d2p_Pairwise_summary_AA",direction="d2p")
pairwise_corr_fun(ssadda_input = ssadda_rep,pairwise_cor_output = "d2p_Pairwise_summary_rep",direction="d2p")
pairwise_corr_fun(ssadda_input = ssadda_rep_EA,pairwise_cor_output = "d2p_Pairwise_summary_rep_EA",direction="d2p")
pairwise_corr_fun(ssadda_input = ssadda_rep_AA,pairwise_cor_output = "d2p_Pairwise_summary_rep_AA",direction="d2p")

#P2D
pairwise_corr_fun(ssadda_input = ssadda,pairwise_cor_output = "p2d_Pairwise_summary",direction="p2d")
pairwise_corr_fun(ssadda_input = ssadda_EA,pairwise_cor_output = "p2d_Pairwise_summary_EA",direction="p2d")
pairwise_corr_fun(ssadda_input = ssadda_AA,pairwise_cor_output = "p2d_Pairwise_summary_AA",direction="p2d")
pairwise_corr_fun(ssadda_input = ssadda_rep,pairwise_cor_output = "p2d_Pairwise_summary_rep",direction="p2d")
pairwise_corr_fun(ssadda_input = ssadda_rep_EA,pairwise_cor_output = "p2d_Pairwise_summary_rep_EA",direction="p2d")
pairwise_corr_fun(ssadda_input = ssadda_rep_AA,pairwise_cor_output = "p2d_Pairwise_summary_rep_AA",direction="p2d")





#########################################################################################
#                                                                                       #
#                                         PLOTS                                         #
#                                                                                       #
#########################################################################################
par(mfrow=c(2,2))
boxplot(as.numeric(ssadda_AA[which(ssadda_AA$aspd==1),]$F2_CocUseAgeO),as.numeric(ssadda_AA[which(ssadda_AA$aspd==2),]$F2_CocUseAgeO),ylab="Coc age of onset", xlab="ASPD(-) and ASPD(+)",main="AA discovery")
boxplot(as.numeric(ssadda_rep_AA[which(ssadda_rep_AA$aspd==1),]$F2_CocUseAgeO),as.numeric(ssadda_rep_AA[which(ssadda_rep_AA$aspd==2),]$F2_CocUseAgeO),ylab="Coc age of onset", xlab="ASPD(-) and ASPD(+)",main="AA replication")
boxplot(as.numeric(ssadda_EA[which(ssadda_EA$aspd==1),]$F2_CocUseAgeO),as.numeric(ssadda_EA[which(ssadda_EA$aspd==2),]$F2_CocUseAgeO),ylab="Coc age of onset", xlab="ASPD(-) and ASPD(+)",main="EA discovery")
boxplot(as.numeric(ssadda_rep_EA[which(ssadda_rep_EA$aspd==1),]$F2_CocUseAgeO),as.numeric(ssadda_rep_EA[which(ssadda_rep_EA$aspd==2),]$F2_CocUseAgeO),ylab="Coc age of onset", xlab="ASPD(-) and ASPD(+)",main="EA replication")


#Gap in age of onsets
hist(as.numeric(as.character(ssadda_AA$I11C_1_TkeAdvAgeO)) - as.numeric(as.character(ssadda_AA$D4D_CigAgeO)))

par(mfrow=c(2,2))
barplot(table(as.numeric(as.character(ssadda_AA$I11C_1_TkeAdvAgeO)) - as.numeric(as.character(ssadda_AA$D4D_CigAgeO))),ylab="Samples",xlab="AA age gap",main="I11C_1_TkeAdv -> D4D_CigAgeO")
barplot(table(as.numeric(as.character(ssadda_rep_AA$I11C_1_TkeAdvAgeO)) - as.numeric(as.character(ssadda_rep_AA$D4D_CigAgeO))),ylab="Samples",xlab="AA(rep) age gap")
barplot(table(as.numeric(as.character(ssadda_EA$I11C_1_TkeAdvAgeO)) - as.numeric(as.character(ssadda_EA$D4D_CigAgeO))),ylab="Samples",xlab="EA age gap")
barplot(table(as.numeric(as.character(ssadda_rep_EA$I11C_1_TkeAdvAgeO)) - as.numeric(as.character(ssadda_rep_EA$D4D_CigAgeO))),ylab="Samples",xlab="EA(rep) age gap")

par(mfrow=c(2,2))
barplot(table(as.numeric(as.character(ssadda_AA$I42A_1_DebtAgeO)) - as.numeric(as.character(ssadda_AA$F4B_CocHghDyAgeO))),ylab="Samples",xlab="AA age gap",main="I42A_1_Debt -> F4B_CocHghDy")
barplot(table(as.numeric(as.character(ssadda_rep_AA$I42A_1_DebtAgeO)) - as.numeric(as.character(ssadda_rep_AA$F4B_CocHghDyAgeO))),ylab="Samples",xlab="AA(rep) age gap")
barplot(table(as.numeric(as.character(ssadda_EA$I42A_1_DebtAgeO)) - as.numeric(as.character(ssadda_EA$F4B_CocHghDyAgeO))),ylab="Samples",xlab="EA age gap")
barplot(table(as.numeric(as.character(ssadda_rep_EA$I42A_1_DebtAgeO)) - as.numeric(as.character(ssadda_rep_EA$F4B_CocHghDyAgeO))),ylab="Samples",xlab="EA(rep) age gap")

par(mfrow=c(2,2))
barplot(table(as.numeric(as.character(ssadda_AA$I18A_2_VandalAgeO)) - as.numeric(as.character(ssadda_AA$F20C_Exp2BxAgeO))),ylab="Samples",xlab="AA age gap",main="I18A_2_Vandal -> F20C_Exp2BxAgeO")
barplot(table(as.numeric(as.character(ssadda_rep_AA$I18A_2_VandalAgeO)) - as.numeric(as.character(ssadda_rep_AA$F20C_Exp2BxAgeO))),ylab="Samples",xlab="AA(rep) age gap")
barplot(table(as.numeric(as.character(ssadda_EA$I18A_2_VandalAgeO)) - as.numeric(as.character(ssadda_EA$F20C_Exp2BxAgeO))),ylab="Samples",xlab="EA age gap")
barplot(table(as.numeric(as.character(ssadda_rep_EA$I18A_2_VandalAgeO)) - as.numeric(as.character(ssadda_rep_EA$F20C_Exp2BxAgeO))),ylab="Samples",xlab="EA(rep) age gap")





#########################################################################################
#                                                                                       #
#                                         STATS                                         #
#                                                                                       #
#########################################################################################


#Find Odds ratios for each of 4 datasets in the Left -> Right direction
t <- table(ssadda_AA[which(as.numeric(as.character(ssadda_AA$I7A_ChlngAuthAgeO))<as.numeric(as.character(ssadda_AA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_AA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_AA$aspd)[2]-as.numeric(t[2]))),nrow=2))

t <- table(ssadda_rep_AA[which(as.numeric(as.character(ssadda_rep_AA$I7A_ChlngAuthAgeO))<as.numeric(as.character(ssadda_rep_AA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_rep_AA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_rep_AA$aspd)[2]-as.numeric(t[2]))),nrow=2))

t <- table(ssadda_EA[which(as.numeric(as.character(ssadda_EA$I7A_ChlngAuthAgeO))<as.numeric(as.character(ssadda_EA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_EA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_EA$aspd)[2]-as.numeric(t[2]))),nrow=2))

t <- table(ssadda_rep_EA[which(as.numeric(as.character(ssadda_rep_EA$I7A_ChlngAuthAgeO))<as.numeric(as.character(ssadda_rep_EA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_rep_EA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_rep_EA$aspd)[2]-as.numeric(t[2]))),nrow=2))


#Find Odds ratios for each of 4 datasets i nthe Left <- Right direction
t <- table(ssadda_AA[which(as.numeric(as.character(ssadda_AA$I7A_ChlngAuthAgeO))>as.numeric(as.character(ssadda_AA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_AA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_AA$aspd)[2]-as.numeric(t[2]))),nrow=2))

t <- table(ssadda_rep_AA[which(as.numeric(as.character(ssadda_rep_AA$I7A_ChlngAuthAgeO))>as.numeric(as.character(ssadda_rep_AA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_rep_AA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_rep_AA$aspd)[2]-as.numeric(t[2]))),nrow=2))

t <- table(ssadda_EA[which(as.numeric(as.character(ssadda_EA$I7A_ChlngAuthAgeO))>as.numeric(as.character(ssadda_EA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_EA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_EA$aspd)[2]-as.numeric(t[2]))),nrow=2))

t <- table(ssadda_rep_EA[which(as.numeric(as.character(ssadda_rep_EA$I7A_ChlngAuthAgeO))>as.numeric(as.character(ssadda_rep_EA$I1B_HookyAgeO))),]$aspd)
fisher.test(matrix(c(as.numeric(t[3]),(table(ssadda_rep_EA$aspd)[3]-as.numeric(t[3])),as.numeric(t[2]),(table(ssadda_rep_EA$aspd)[2]-as.numeric(t[2]))),nrow=2))


###Log Rge for ASPD~age-of-onset
t <- na.omit(as.data.frame((cbind(as.character(ssadda_AA$I23C_1_IlglAgeO),as.character(ssadda_AA$aspd)))))
t <- t[which(t$V1!="UC"),]
t <- t[which(t$V2!=0),]
t <- cbind(t,as.numeric(as.character(t$V2))-rep(1,nrow(t)))
mod <- glm(t[,3]~as.numeric(as.character(t[,1])),family="binomial")
summary(mod)
exp(coef(mod))



#########################################################################################
#                                                                                       #
#                                         EXTRA                                         #
#                                                                                       #
#########################################################################################




#Count samples with and without progression

length(which(ssadda_AA$I7_ChlngAuth==5 & ssadda_AA$I1_Hooky == 1))
length(which(ssadda_rep_AA$I7_ChlngAuth==5 & ssadda_rep_AA$I1_Hooky == 1))
length(which(ssadda_EA$I7_ChlngAuth==5 & ssadda_EA$I1_Hooky == 1))
length(which(ssadda_rep_EA$I7_ChlngAuth==5 & ssadda_rep_EA$I1_Hooky == 1))

length(which(ssadda_AA$I2_Expell==5 & ssadda_AA$I1_Hooky == 1))
length(which(ssadda_rep_AA$I2_Expell==5 & ssadda_rep_AA$I1_Hooky == 1))
length(which(ssadda_EA$I2_Expell==5 & ssadda_EA$I1_Hooky == 1))
length(which(ssadda_rep_EA$I2_Expell==5 & ssadda_rep_EA$I1_Hooky == 1))

length(which(ssadda_AA$I7_ChlngAuth==5 & ssadda_AA$I4_StyOut == 1))
length(which(ssadda_rep_AA$I7_ChlngAuth==5 & ssadda_rep_AA$I4_StyOut == 1))
length(which(ssadda_EA$I7_ChlngAuth==5 & ssadda_EA$I4_StyOut == 1))
length(which(ssadda_rep_EA$I7_ChlngAuth==5 & ssadda_rep_EA$I4_StyOut == 1))

length(which(ssadda_AA$H10B_MJ2Prb==5 & ssadda_AA$F3_CocDaily == 1))
length(which(ssadda_rep_AA$H10B_MJ2Prb==5 & ssadda_rep_AA$F3_CocDaily == 1))
length(which(ssadda_EA$H10B_MJ2Prb==5 & ssadda_EA$F3_CocDaily == 1))
length(which(ssadda_rep_EA$H10B_MJ2Prb==5 & ssadda_rep_EA$F3_CocDaily == 1))

length(which(ssadda_AA$H10B_MJ2Prb==5 & ssadda_AA$F3_CocDaily == 1))
length(which(ssadda_rep_AA$H10B_MJ2Prb==5 & ssadda_rep_AA$F3_CocDaily == 1))
length(which(ssadda_EA$H10B_MJ2Prb==5 & ssadda_EA$F3_CocDaily == 1))
length(which(ssadda_rep_EA$H10B_MJ2Prb==5 & ssadda_rep_EA$F3_CocDaily == 1))


length(which(ssadda_AA$G1_OpiEver==5 & ssadda_AA$F1_CocEver == 1))
length(which(ssadda_rep_AA$G1_OpiEver==5 & ssadda_rep_AA$F1_CocEver == 1))
length(which(ssadda_EA$G1_OpiEver==5 & ssadda_EA$F1_CocEver == 1))
length(which(ssadda_rep_EA$G1_OpiEver==5 & ssadda_rep_EA$F1_CocEver == 1))

length(which(ssadda_AA$I18_Vandal!=1 & ssadda_AA$F20A_CocExp3 == 1))
length(which(ssadda_rep_AA$I18_Vandal!=1 & ssadda_rep_AA$F20A_CocExp3 == 1))
length(which(ssadda_EA$I18_Vandal!=1 & ssadda_EA$F20A_CocExp3 == 1))
length(which(ssadda_rep_EA$I18_Vandal!=1 & ssadda_rep_EA$F20A_CocExp3 == 1))

length(which(ssadda_AA$F20A_CocExp3==5 & ssadda_AA$E33_TrtmtProg == 1))
length(which(ssadda_rep_AA$F20A_CocExp3==5 & ssadda_rep_AA$E33_TrtmtProg == 1))
length(which(ssadda_EA$F20A_CocExp3==5 & ssadda_EA$E33_TrtmtProg == 1))
length(which(ssadda_rep_EA$F20A_CocExp3==5 & ssadda_rep_EA$E33_TrtmtProg == 1))

####
length(which(ssadda_AA$E26H_WDSymDrnk==5 & ssadda_AA$E26C_WDSym2 == 1))
length(which(ssadda_rep_AA$E26H_WDSymDrnk==5 & ssadda_rep_AA$E26C_WDSym2 == 1))
length(which(ssadda_EA$E26H_WDSymDrnk==5 & ssadda_EA$E26C_WDSym2 == 1))
length(which(ssadda_rep_EA$E26H_WDSymDrnk==5 & ssadda_rep_EA$E26C_WDSym2 == 1))

length(which(ssadda_AA$E32_SlfHlp==5 & ssadda_AA$E33_TrtmtProg == 1))
length(which(ssadda_rep_AA$E32_SlfHlp==5 & ssadda_rep_AA$E33_TrtmtProg == 1))
length(which(ssadda_EA$E32_SlfHlp==5 & ssadda_EA$E33_TrtmtProg == 1))
length(which(ssadda_rep_EA$E32_SlfHlp==5 & ssadda_rep_EA$E33_TrtmtProg == 1))

#I23C_1_IlglAgeO

length(which(ssadda_AA$I23_1_BadChk==5 &  ssadda_AA$I23_2_StolnGood==5 & ssadda_AA$I23_3_PdSex==5 & ssadda_AA$I23_4_Pimp==5 & ssadda_AA$I28_TrfkViol != 5))
length(which(ssadda_rep_AA$I23_1_BadChk==5 &  ssadda_rep_AA$I23_2_StolnGood==5 & ssadda_rep_AA$I23_3_PdSex==5 & ssadda_rep_AA$I23_4_Pimp==5 & ssadda_rep_AA$I28_TrfkViol != 5))
length(which(ssadda_EA$I23_1_BadChk==5 &  ssadda_EA$I23_2_StolnGood==5 & ssadda_EA$I23_3_PdSex==5 & ssadda_EA$I23_4_Pimp==5 & ssadda_EA$I28_TrfkViol != 5))
length(which(ssadda_rep_EA$I23_1_BadChk==5 &  ssadda_rep_EA$I23_2_StolnGood==5 & ssadda_rep_EA$I23_3_PdSex==5 & ssadda_rep_EA$I23_4_Pimp==5 & ssadda_rep_EA$I28_TrfkViol != 5))


length(which(ssadda_AA$F22_Exp3Bx==5 & ssadda_AA$F5_CocDes == 1))
length(which(ssadda_rep_AA$F22_Exp3Bx==5 & ssadda_rep_AA$F5_CocDes == 1))
length(which(ssadda_EA$F22_Exp3Bx==5 & ssadda_EA$F5_CocDes == 1))
length(which(ssadda_rep_EA$F22_Exp3Bx==5 & ssadda_rep_EA$F5_CocDes == 1))


length(which(ssadda_AA$F5_CocDes==5 & ssadda_AA$E31_AlcPro == 1))
length(which(ssadda_rep_AA$F5_CocDes==5 & ssadda_rep_AA$E31_AlcPro == 1))
length(which(ssadda_EA$F5_CocDes==5 & ssadda_EA$E31_AlcPro == 1))
length(which(ssadda_rep_EA$F5_CocDes==5 & ssadda_rep_EA$E31_AlcPro == 1))


#Run pairwise correlations and calculate p-values
require(corrplot)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p_mat <- cor.mtest(ssadda_ageOfOnset)
p_mat[is.na(p_mat)] <- 0

cor_mat <- cor(ssadda_ageOfOnset)

corrplot(cor_mat,method = "square",outline = F,p.mat = p_mat,sig.level = 0.05/136,order = "original",insig = "blank",pch = ".",pch.cex = 1.5)
