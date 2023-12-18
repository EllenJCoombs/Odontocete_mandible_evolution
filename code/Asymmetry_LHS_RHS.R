
###############################################
#                                             #
#  Asymmetry - LHS to RHS reflection          #
#  LandvR calculation of Euclidean distances  #
#                                             #
###############################################

#LHS LMs and curves 
LHS <- final_mands[c(1:18, 37:326),,]
#RHS LMs and curves 
RHS <- final_mands[c(19:36, 327:616),,]


#rearrange LHS so we can match LHS to RHS
LHS_arranged=abind::abind(LHS[c(1:9),,], LHS[14:18,,],
                          LHS[c(10:13),,], 
                          LHS[c(19:308),,],
                          along= 1)

#If RHS needs to be rearranged too (in older data)
RHS_arranged=abind::abind(RHS[c(1:8),,],
                          RHS[c(15:18),,], 
                          RHS[c(9:14),,], 
                          #RHS[c(13:16),,],
                          #RHS[c(19:308),,],
                          along= 1)

#Overwrite the naming system if needed - wanted names/numbers first, current names/numbers second 
RHS_arranged[c(11, 12, 13), , ] <- RHS_arranged[c(13, 11, 12), , ]

#Hereafter 'LHS' refers to the reflected LHS landmarks (reflected to the right)

#reflect the LHS to the RHS
arr.refl <-(-1)*(LHS_arranged[1:length(1:nrow(LHS_arranged)),,])

#Bind the data set and Procrustes together 

#Bind the 2 datasets 
final_dataset_LHS_ref_RHS=abind::abind(RHS_arranged,
                           arr.refl,
                           along = 3) 


species_data <- read.csv('RHS_and_reflected.csv')
#pull the phylo names from the species data 
Full_names_phylo=species_data$phylo_name

#makes the species phylo.names the names of the array 
dimnames(final_dataset_LHS_ref_RHS)[[3]]<-Full_names_phylo 


#Procrustes the data 

final_data_LHS_ref_RHS_PROC=geomorph::gpagen(final_dataset_LHS_ref_RHS) #Remove non-shape aspects 
final_data_asymm_COORDS=final_data_LHS_ref_RHS_PROC$coords #Subset out the coords 


#PCA
pca_res <- geomorph::gm.prcomp(final_data_asymm_COORDS)

pca_summary <- summary(pca_res)
#mutate the species data to add to the matrix

plot(pca_res)



#pull out the pairwise distances 

distances_LHS <- dist(two.d.array(final_data_asymm_COORDS[,,1:100])) %>%
  broom::tidy(.) %>%
  mutate(., tog.or.sep = "LHS")

distances_RHS <- dist(two.d.array(final_data_asymm_COORDS[,,101:200])) %>%
  broom::tidy(.) %>%
  mutate(., tog.or.sep = "RHS")

#combine datasets

tibble2 <- tibble(Distances.LHS = distances_LHS$distance, Distances.RHS = distances_RHS$distance)

#Pairwise distance plot 

distdistplot <- ggplot(data = tibble2, aes(x = Distances.LHS, y = Distances.RHS)) +
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "Pairwise Distances- LHS landmarks", y = "Pairwise Distances- RHS landmarks") +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical"
  ) +
  ggtitle("Pairwise Distance vs Pairwise Distance") +
  stat_smooth(
    mapping = aes(x = Distances.LHS, y = Distances.RHS),
    se = TRUE,
    data = tibble2,
    method = lm, fill = "grey",
    formula = y ~ x, inherit.aes = FALSE
  )
distdistplot

#t-test 

t.test(x = tibble2$Distances.LHS, y = tibble2$Distances.RHS)

#do as above for archs and odonts 


#Now make the density plots 

cbbPalette <- c("#000000", #black
                "#E69F00", #orange
                "#56B4E9", #lightblue
                "#009E73", #green
                "#F0E442", #yellow
                "#0072B2", #darkblue
                "#D55E00", #red
                "#CC79A7") #pink


comparison <- rep(c("LHS","RHS"),each=4950)

t1 <- tibble(distances = c(tibble2$Distances.LHS,tibble2$Distances.RHS))
t1 <- t1 %>% mutate(label = comparison)


g2<- ggplot(t1, aes(x = distances, fill = label)) +
  geom_density(alpha= .7) +
  theme_bw()+
  scale_fill_manual(name = NULL, values=c(cbbPalette[2],cbbPalette[7]), labels = c("LHS", "RHS"))+
  xlab("Euclidean Distance") +
  ylab("Density")+
  
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"))
g2


#For landvR - plotting Eucidean distances by colour

#(was x)
LHS_dists <- final_data_asymm_COORDS[,,c(1:100)] 
#was y)
RHS_dists <- final_data_asymm_COORDS[,,c(101:200)]


if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

############# LANDVR ###################
#Check the values
MirroredAC[1,,2] == Proc_manual[1,,2]
# FALSE FALSE FALSE
## This is due to rounding, in fact they have the same 9 digits - round them 
round(MirroredAC[1,,2], digits = 9) == round(Proc_manual[1,,2], digits = 9)
# TRUE TRUE TRUE

#I’ve updated landvR to version 0.3 where the coordinates.difference function now have a tolerance optional argument.
#You can use the following to get the 0 difference results:

differences_between_lms <- coordinates.difference(coordinates = LHS_dists[,,1],
                                                  reference = RHS_dists[,,1],
                                                  type = "spherical",
                                                  rounding = 9)

#Remove errornous missing landmarks (these should be zero because they are static)
differences_between_lms[[1]][1:18, 1:3] <- c(0.000000, 0.000000, 0.000000)


#Ellen's own colour function 
colfunc <- colorRampPalette(c("red", "yellow", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

get.col.spectrum <- landvR::procrustes.var.plot(LHS_dists[,,85], RHS_dists[,,85], col.val = differences_between_lms[[1]][,1], col = colfunc)

test=differences_between_lms[[1]][,1] #this is a test for specimen 1 to look at the differences between lms 
test

##### LOOKING AT AN AVERAGE SPECIMEN ######
N=308 #number of landmarks 
specs=100 #number of specimens 
all_combined=array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)

i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = x[,,i],
                                            reference = y[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}


#55, 56, 57, 59, 60 are all missing data and should be zero 
all_combined[1:18, 1:3, 1:94] <- c(0.000000, 0.000000, 0.000000)
#write.csv(all_combined, file = 'all_combined.csv')

radii=all_combined[,1,] #looking at the second column (usually x,y,z) but here it is the radii, aziumuth, and polar 
radii_mean=apply(radii, c(1), mean) #c(1) looking at the first column which is the radii 
#test=all_combined[[1]][,,1] #this is a test for specimen 1 to look at the differences between lms 
#test

radii=all_combined[,1,] #second column of whole dataset with just the radii [,1,]


############################################################################
#Reararrange LHS and mandibles - these code fed into the above 
#rearrange LHS 
LHS_arranged=abind::abind(LHS[c(1:9),,], LHS[14:18,,],
                             LHS[c(10:13),,], 
                             LHS[c(19:308),,],
                             along= 1)


arranged_mands=abind::abind(final_mands[c(1:26),,], final_mands[33:34,,],
                             final_mands[c(27),,],
                             final_mands[c(35:36),,],
                             final_mands[c(28),,],
                             final_mands[c(29:32),,],
                             final_mands[c(37:896),,],
                             along= 1)


LHS_arranged <- arranged_mands[c(1:18, 37:326),,]
RHS_arranged <- arranged_mands[c(19:36, 327:616),,]

spheres3d(x[,,89], col = 'green', radius = 0.0004)
rglwidget()
spheres3d(y[,,89], col = 'red', radius = 0.0004)
rglwidget()
