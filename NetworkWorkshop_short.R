library(statnet) #includes the libraries network, sna, and ergm
library(tnet) #includes the library igraph
library(vegan)
library(FastKNN)
library(kableExtra)
library(ggraph)
library(GGally)

# Import ceramic ceramic data into an object named ceramic. row.names=1 sets the first column values to the row names for each item.
ceramic <- read.csv(file='ceramic_clust.csv',row.names=1)
# Import attribute data into an object named ceramic.attr
ceramic.attr <- read.csv('ceramic_clust_attr.csv',row.names=1)

########################################
##### Material Co-Presence
co.p <- function (x,thresh=0.25) {
  temp <- prop.table(as.matrix(x),1) #create matrix of proportions from ceramic    
  temp[temp>=thresh] <- 1 # define anything with greater than or equal to 0.1 as present (1)
  temp[temp<1] <- 0   # define all other cells as absent (0)
  out <- temp%*%t(temp)   # matrix algebraic calculation to find co-occurence (%*% indicates matrix multiplication)
  return(out)}
ceramicP <- co.p(ceramic) # run the function
head(ceramicP)

########################################
##### Create Network Object from P/A
Pnet <- network(ceramicP,directed=F) # create network object from co-occurrence
Pnet %v% 'vertex.names' <- row.names(ceramicP) # Now let's add names for our nodes based on the row names of our original matrix
Pnet  # look at the results
plot(Pnet, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,main='co-presence network') # plot network using default layout
plot(Pnet, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,coord=ceramic.attr,main='co-presence network') # plot network using geographic coordinates

########################################
##### Brainard-Robinson Similarity
ceramic.p <- prop.table(as.matrix(ceramic), margin = 1) # This line converts the ceramic cluster frequency table to a table of proportions by row
ceramicBR <- (2-as.matrix(vegdist(ceramic.p, method='manhattan')))/2 #  Use the vegdist function to calculate the Brainard-Robinson similarity score. 
head(ceramicBR)

########################################
##### Create Network Object from BR Sim
BRnet <- network(event2dichot(ceramicBR,method='absolute',thresh=0.65),directed=F) # Define our binary network object from BR similarity
BRnet %v% 'vertex.names' <- row.names(ceramicBR) # Now let's add names for our nodes based on the row names of our original matrix
BRnet # look at the results.
plot(BRnet, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,main='BR network') # plot network using default layout
plot(BRnet, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,coord=ceramic.attr,main='BR network') # plot network using geographic coordinates

########################################
##### Chi-Squared Distance
chi.dist <- function(x) {
  rowprof <- x/apply(x,1,sum) # calculates the profile for every row
  avgprof <- apply(x,2,sum)/sum(x) # calculates the average profile
  chid <- dist(as.matrix(rowprof)%*%diag(1/sqrt(avgprof))) # creates a distance object of $\chi^{2}$ distances
  return(as.matrix(chid))}   # return the reults
# Run the script and then create the rescaled 0-1 version
ceramicX <- chi.dist(ceramic)
ceramicX01 <- ceramicX/max(ceramicX)
head(ceramicX01)

########################################
##### Create Network Object from X^2 Dist
Xnet <- network(event2dichot(1-ceramicX01,method='quantile',thresh=0.80),directed=F) # Note we use 1 minus ceramicX01 here so to convert a distance to a similarity
Xnet %v% 'vertex.names' <- row.names(ceramicX01)
Xnet
plot(Xnet, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,main='Chi-squared network') # plot network using default layout
plot(Xnet, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,coord=ceramic.attr,main='Chi-squared network') # plot network using geographic coordinates

########################################
##### K Nearest Neighbors
distMatrix <- as.matrix(dist(ceramic.attr)) # Create a distance matrix using Euclidean distances based on the site coordinates in ceramic.attr
# set k as 5 and create a function that calculates the 5 nearest neighbors for each node
k <- 5
  nrst <- lapply(1:nrow(distMatrix), function(i) k.nearest.neighbors(i, distMatrix, k = k))
  # the chunk of code below creates a symmetric matrix of 0s and then fills cells with 1 where two sites are k nearest neighbors
  dist.knn <- matrix(nrow = dim(distMatrix), ncol=dim(distMatrix),0) 
  for(i in 1:length(nrst)) for(j in nrst[[i]]) dist.knn[i,j] = 1
  # set row and column names
  row.names(dist.knn) <- row.names(ceramic.attr)
  colnames(dist.knn) <- row.names(ceramic.attr)
head(dist.knn)

########################################
##### Create Network Object from KNN
dist.net <- network(dist.knn,directed=T) # Create network object with directed ties
dist.net %v% 'vertex.names' <- row.names(ceramic.attr)
dist.net
plot(dist.net, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,main='K nearest neighbors network, k=5') # plot network using default layout
plot(dist.net, edge.col='gray', edge.lwd=0.25, vertex.cex=0.5,coord=ceramic.attr,main='K nearest neighbors network, k=5') # plot network using geographic coordinates

########################################
##### Create Weighted Network Object from P/A Data 
# create weighted network object from co-occurrence matrix by adding the ignore.eval=F argument
Pnet2 <- network(ceramicP,directed=F,ignore.eval=F,names.eval='weight')
Pnet2 %v% 'vertex.names' <- row.names(ceramicP)
Pnet2
plot(Pnet2,edge.col='weight',edge.lwd='weight',vertex.cex=0.5,vertex.col='red',main='co-presence weighted network') # plot weighted network using default layout
plot(Pnet2,edge.col='weight',edge.lwd='weight',vertex.cex=0.5,coord=ceramic.attr,vertex.col='red',main='co-presence weighted network') # plot weighted network using geographic coordinates

########################################
##### Create Weighted Network Object from BR Sim
ceramicBR2 <- ceramicBR
ceramicBR2[ceramicBR2<0.65] <- 0 # set values for similarities less than 0.65 to 0
BRnet2 <- network(ceramicBR2,directed=F,ignore.eval=F,names.eval='weight')
BRnet2 %v% 'vertex.names' <- row.names(ceramicP)
BRnet2
edge.set <- round(get.edge.value(BRnet2,'weight')*100,0)-64 #extract edge values and round, subtract 64 so that the minimum value will be 1
edge.cols <- colorRampPalette(c('gray','darkblue'))(max(edge.set)) #create color pallette for the number of values represented
plot(BRnet2,edge.col=edge.cols[edge.set],edge.lwd='weight',vertex.cex=0.5,vertex.col='red',main='co-presence weighted network') # plot weighted network using default layout
plot(BRnet2,edge.col=edge.cols[edge.set],edge.lwd='weight',vertex.cex=0.5,coord=ceramic.attr,vertex.col='red',main='co-presence weighted network') # plot weighted network using geographic coordinates

########################################
##### Calculate basic network statistics for BRnet
## degree centrality
dg <- as.matrix(sna::degree(BRnet,gmode='graph'))
## eigenvector centrality
eg <- as.matrix(sna::evcent(BRnet)) 
eg <- sqrt((eg^2)*length(eg)) # standardized to start network as is frequently seen in the literature
## betweenness centrality
bw <- sna::betweenness(BRnet,gmode='graph') 

########################################
##### Calculate weighted network statistics for BRnet
# calculate weighted degree as the sum of weights - 1
dg.wt <- as.matrix(rowSums(ceramicBR)-1) 
# calculate weighted eigenvector centrality and rescale
eg.wt <- as.matrix(sna::evcent(ceramicBR)) 
eg.wt <- sqrt((eg.wt^2)*length(eg.wt)) # standardize to start network
# calculate weighted betweenness from the tnet package (we use the suppressWarnings package to avoid notifications)
bw.wt <- suppressWarnings(betweenness_w(ceramicBR,directed=F))[,2] 

########################################
##### Define functions for calculating and tabulating all centrality metrics
# Calculate centrality scores for binary networks
net.stats <- function(y){ 
  # calculate degree centrality
  dg <- as.matrix(sna::degree(y,gmode='graph')) 
  # calculate and scale eigenvector centrality
  eg <- as.matrix(sna::evcent(y)) 
  eg <- sqrt((eg^2)*length(eg)) 
  # calculate betweenness centrality
  bw <- sna::betweenness(y,gmode='graph') 
  # combine centrality scores into matrix
  output <- cbind(dg,eg,bw) 
  rownames(output) <- rownames(as.matrix(y)) 
  colnames(output) <- c('dg','eg','bw') 
  return (output)} # return results of this function

# Calculate centrality scores for weighted networks (similarity matrices)
net.stats.wt <- function(y){ 
  # calculate weighted degree as the sum of weights - 1
  dg.wt <- as.matrix(rowSums(y)-1) 
  # calculate weighted eigenvector centrality and rescale
  eg.wt <- as.matrix(sna::evcent(y)) 
  eg.wt <- sqrt((eg.wt^2)*length(eg.wt))
  # calculate weighted betweenness from the tnet package (we use the suppressWarnings package to avoid notifications)
  bw.wt <- suppressWarnings(betweenness_w(y,directed=F))[,2] 
  output <- cbind(dg.wt,eg.wt,bw.wt) 
  rownames(output) <- rownames(as.matrix(y))
  colnames(output) <- c('dg.wt','eg.wt','bw.wt') 
  return (output)} # return results of this function

########################################
##### Calculate net stats for BRnet and ceramicBR
BR.stats <- net.stats(BRnet)
BR.stats
BR.w.stats <- net.stats.wt(ceramicBR)
BR.w.stats

########################################
##### Plotting network objects (lots of options)
## Using SNA/network package options
?plot.network
plot.network(BRnet,displaylabels=T,label.cex=0.5,edge.col='gray',
             edge.lwd=0.25,vertex.col='blue')

plot.network(BRnet,displaylabels=T,label.cex=0.5,edge.col='gray',
             edge.lwd=0.25,vertex.col='blue',interactive=T)

plot.network(BRnet,displaylabels=T,label.cex=0.5,edge.col='gray',
             edge.lwd=0.25,vertex.col='blue',vertex.cex=dg/6)


## Using GGally and ggraph
BRnet %>% 
  ggraph(layout = 'fr') + 
  geom_edge_link(color='gray') + 
  geom_node_point(aes(size = bw.wt, color = bw.wt)) + 
  scale_color_continuous(guide = 'legend') + 
  theme_graph()

BRnet %>% 
  ggraph(layout = 'manual',  node.positions = ceramic.attr[,1:2]) + 
  geom_edge_link(color='gray') + 
  geom_node_point(aes(size = bw.wt, color = bw.wt)) + 
  scale_color_continuous(guide = 'legend') + 
  theme_graph()

edge_weight <- get.edge.value(BRnet2,'weight')
BRnet2 %>% 
  ggraph(layout = 'fr') + 
  geom_edge_link(aes(alpha=edge_weight), color='black') + 
  geom_node_point(aes(size = bw.wt, color = bw.wt)) + 
  geom_node_text(aes(label = name), repel=T)+
  scale_color_continuous(guide = 'legend') + 
  theme_graph()
