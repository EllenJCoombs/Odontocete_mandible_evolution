
#############################################
#                                           #
#  Code to import, resample, and slide LMs  #
#                                           #
#############################################

#load packages 

rm(list=ls()) #to clear workspace 
library(tidyverse)
library(broom) #needed with tidyverse
library(latticeExtra) #qgraph
library(readr)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph) #load in packages bit 
library(ape)
library(geiger)
library(abind)
library(tibble)
library(dplyr)
library(tidyr)
library(rgl)
library("devtools")
devtools::install_github("rnfelice/SURGE")
library(SURGE)
library(RColorBrewer) # for color palettes
library(magick)


#set path to PTS final LHS - this is where all of the curves are 

#import table defining curves
#different curves table for the odontocetes 
curve_table <- read_csv('new_curves.csv')

#identify the folder where your pts files are
ptsfolder <- "D:/Smithsonian Postdoc/Year 1/Mandible scans/LMs and curves/pts"


#import the pts file names
ptslist <- dir(ptsfolder, pattern='.pts', recursive=F)

my_curves <- create_curve_info(curve_table, n_fixed = 36)
setwd(ptsfolder)
subsampled.lm <- import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE)

#capture the output if it's too big to read in the console 
capture.output(import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE), file = "odonts.txt", append = T)
#if you have any missing points, Checkpoint will set them to x,y,z=9999
#this makes for bad plots in checkLM, so switch them to NA
subsampled.lm[subsampled.lm == 9999] <- NA
subsampled.lm[subsampled.lm == Inf] <- NA
subsampled.lm[subsampled.lm == -Inf] <- NA


#SET WORKING DIRECTORY TO ASCII PLY
#check to make sure your curves look okay on each specimen
checkLM(subsampled.lm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)


newpts <- subsampled.lm

#Create missing list 
misslist<-createMissingList(dim(newpts)[3])
for (j in 1:dim(newpts)[[3]]){
  misslist[[j]]<-which(is.na(newpts[,1,j]))
} 
newpts2<-fixLMtps(newpts)
#check that your pts and ply match
plyfolder <- "D:/Smithsonian Postdoc/Year 1/Mandible scans/ply ASCII"
ptslist2<-gsub(pattern="\\.pts$","",ptslist)
plylist <-  dir(plyfolder, pattern='.ply', recursive=F)
plylist2<-gsub(pattern="\\.ply$","",plylist)
setdiff(plylist2,ptslist2) #should be zero if all matching up OK
#######NOTE
#####my plys are awkwardly in a folder called  "C:/Data/whale_ply/ply/ply"
##### I am setting the working directory ONE FOLDER UP from where the meshes are

#mc.cores = 1 for WINDOWS. mc.cores = 3 for MAC
####so for you it should be 

#################################################
#                                               #
#        RUN THE SLIDER 3D_2 CODE               #
#                                               #
#################################################

##################################
#                                #
#    STEP SIZE - 0.2 or 0.002    #
#                                #
##################################


setwd("D:/Smithsonian Postdoc/Year 1/Mandible scans/ply ASCII/ply")
{
  slided4.all <- slider3d_2(newpts2$out, SMvector= my_curves$Sliding.LMs,
                            outlines = my_curves$Curves, sur.path = "./ply", sur.name = NULL, 
                            meshlist = paste("./ply/",dimnames(newpts2$out)[[3]],".ply",sep=""), ignore = NULL,
                            sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                            recursive = TRUE, iterations = 3, initproc = TRUE,
                            pairedLM = 0, mc.cores = 1, bending=TRUE,
                            fixRepro = FALSE,stepsize=0.02,
                            missingList=misslist)
  dimnames(slided4.all[["dataslide"]])[3]<-dimnames(newpts2$out)[3]
  #save(slided4.all,file="~/Google Drive/NHM/crocs/data/slid.crocs.all.apr25.R")
  #load("~/Google Drive/NHM/crocs/data/slid.crocs.sept23.R")
}



#re-estimate missing post sliding
slided4.all$dataslide[which(is.na(newpts))]<-NA
#Fix missing landmarks
slid.lms<-fixLMtps(slided4.all$dataslide)
#the landmarks ready for mirroring are in an object called slid.lms$out


slidedlms <- slid.lms$out

checkLM(slidedlms,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)

save(slidedlms, file = 'Final_dataset_LMs_102.R')
