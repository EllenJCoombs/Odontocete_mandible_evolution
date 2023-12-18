
###################################
#                                 #
#    BayesTraits  MCMC plots      #
#                                 #
###################################


#If zip package doesn't work - need:
devtools::install_github("hferg/btrtools", force = TRUE)

#Plot caterpillar plots
                                                        
library(geiger)
library(phytools)
library(BTRTools) #load from packages - install - package archive - folder must be zipped 
library(coda)
#Codes for mcmc chains diagnostics:
tracePlots <- function(file, burnin=0, thinning=1, plot=TRUE, display=c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")){
  require(BTRTools)
  require(coda)
  rjout <- loadRJ(file, burnin = burnin, thinning = thinning)
  chain_out <- type.convert(rjout$rj_output)
  rownames(chain_out) = chain_out[,"It"]
  chain_out = chain_out[,-1]
  # Retrieve numerical
  index <- sapply(chain_out,function(x) is.numeric(x))
  chain <- mcmc(chain_out[,index])
  # plot the trace
  if(plot){
    plot(chain[,display])
  }
  # Just compute some statistics (autocorrelation...)
  cat("Effective sample size:","\n")
  print(effectiveSize(chain[,display]))
  # return results
  invisible(chain)
}
setwd("D:/BayesTraits/Lambda_results")
#### EXAMPLE ####
test = tracePlots(file=paste("whole_var.txt.VarRates.txt"), burnin=0)
#tracePlots(file=paste("./BT_adj_shape_skull_95/Tree_173species_Adj_shape_95_pPCA_tree_order_BM/Adj_shape_95_pPCA_tree_order_BM.txt.VarRates.txt"), burnin=0)
#(This will return the ESS (effective sample size) and trace plot)
#  If you need to compare several independent runs (e.g. using general mcmc diagnostic functions in "coda" package)
test2 = tracePlots(file=paste("whole2_var.txt.VarRates.txt"), burnin=0) 
#("./BT_adj_shape_skull_95/Tree_173species_Adj_shape_95_pPCA_tree_order_BM_2/Adj_shape_95_pPCA_tree_order_BM_2.txt.VarRates.txt"),burnin=0) ##another dataset
my_list_of_chains = mcmc.list(list(test[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")],
                                   test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
#  Then simply use any functions from "coda" that works with objects of class "mcmc.list"
# e.g.
plot(my_list_of_chains)
gelman.diag(my_list_of_chains, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)
