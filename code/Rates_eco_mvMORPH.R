
###########################################
#                                         #
#  Rates per bone/eco category - mvMORPH  #
#                                         #
###########################################

install.packages('devtools')
library(devtools)
install_github("JClavel/mvMORPH")
library(mvMORPH)

# phylogenetic tree
treeMAND <- read.nexus("treeMAND2.nexus")

#All of my species data is already rearranged to the same order as the tree
species_data <- read.csv('species_data_phylo_no_mysts.csv', row.names = 2, fileEncoding="UTF-8-BOM")

#Pull out the factor to paint the tree with below

#Eco cats e.g. dentition
fact_diet <- species_data$dentition

# Body size
fact_size <- read.csv("Csize.csv", row.names = 1)

# Morphology
data <- read.csv("jaw_var.csv", row.names = 1, fileEncoding="UTF-8-BOM")

# Prepare data for analyses and check that tips match the data labels
dat <- list(Y=as.matrix(data[treeMAND$tip.label,]),
            size=as.numeric(fact_size[treeMAND$tip.label,]))
#Add the function to reconstruct from a sample of trees
paintAllTree <- function(tree, ancestral, tips){
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }
  return(treebis)
}
#What category to paint your simmap tree with
cat_suborder = as.factor(fact_suborder)
names(cat_suborder) = treeMAND$tip.label

# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates.
#Ideally this should be performed on a sample of stochastic maps.
#Or by using a mapping constructed from the ML reconstruction for instance

simm_tr<- make.simmap(treeMAND, model="ARD", cat_suborder, nsim=100)
a_single_simmap_tree <- paintAllTree(treeMAND, describe.simmap(simm_tr), as.character(cat_suborder))

# fit linear model > here I include size in the model as a covariate. This should be equivalent to estimating rates on size corrected dat
fit <- mvgls(Y~size, data=dat, tree=a_single_simmap_tree, model="BMM", method = "PL", error=TRUE) # if you remove the error argument, then it is assumed that there's no ME nor intraspecific variance that may blur the results.

# Save results
results <- list(fit=fit$param, tree=a_single_simmap_tree)
save(results, file="suborder_PCscores_pan.Rdata")
results$fit
  
