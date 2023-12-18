#Code from Andy Knapp - ecomorphospace 

library(phytools)
library(ggplot2)
library(musculusColors)
library(gridExtra)
library(geomorph)

#classifer.whales is the species data
#Add species data 

classifier.whales <- read.csv('species_data_phylo_with_mysts.csv')

#load shape data 
load('final_mands_102_tree.R')

#load tree 
treeMAND <- read.nexus('treeMAND.nexus')

# Firstly, perform Procrustes on shape data
Y.gpa=gpagen(final_mands_102)
Y.gpa <-Y.gpa$coords
PCA.3D.W <- geomorph::gm.prcomp( Y.gpa )

PCW <-data.frame(cbind(PCA.3D.W$x[,c(1:14)],classifier.whales,dimnames(Y.gpa)[[3]]))

###############################################################################
#### Get tree data and convert to geom_segment to create phylomorphospace #####
###############################################################################

## create matrix for use in phytools::fastAnc()
mat1 <- cbind(eval(substitute(PCW$Comp1), PCW),eval(substitute(PCW$Comp2), PCW))
mat2 <- cbind(eval(substitute(PCW$Comp1), PCW),eval(substitute(PCW$Comp3), PCW))
rownames(mat1) <- eval(substitute(PCW$dimnames.Y.gpa...3..), PCW)
rownames(mat2) <- eval(substitute(PCW$dimnames.Y.gpa...3..), PCW)
stopifnot(length(setdiff(treeMAND$tip.label, rownames(mat1))) == 0)
stopifnot(length(setdiff(treeMAND$tip.label, rownames(mat2))) == 0)

xAnc <- fastAnc(treeMAND, mat1[,1])
yAnc <- fastAnc(treeMAND, mat1[,2])
zAnc <- fastAnc(treeMAND, mat2[,2])

all_node_coords.1 <-
  data.frame(
    #put PC values for all nodes and tips in a dataframe
    #tips go first in order of tip labels, then numerical order for nodes
    x=c(mat1[treeMAND$tip.label,1], xAnc),
    y=c(mat1[treeMAND$tip.label,2], yAnc),
    nodeid=1:(treeMAND$Nnode + length(treeMAND$tip.label))
  )
all_node_coords.2 <-
  data.frame(
    #put PC values for all nodes and tips in a dataframe
    #tips go first in order of tip labels, then numerical order for nodes
    x=c(mat2[treeMAND$tip.label,1], xAnc),
    y=c(mat2[treeMAND$tip.label,2], zAnc),
    nodeid=1:(treeMAND$Nnode + length(treeMAND$tip.label))
  )

#get edge list from tree object
edges.1 <- data.frame(treeMAND$edge)
names(edges.1) <- c("node1", "node2")
#translate tip/node numbers into PC coordinates in all_node_coords dataframe
edgecoords.1 <- merge(
  merge(edges.1, all_node_coords.1, by.x="node1", by.y="nodeid"),
  all_node_coords.1, by.x="node2", by.y="nodeid")

#get edge list from tree object
edges.2 <- data.frame(treeMAND$edge)
names(edges.2) <- c("node1", "node2")
#translate tip/node numbers into PC coordinates in all_node_coords dataframe
edgecoords.2 <- merge(
  merge(edges.2, all_node_coords.2, by.x="node1", by.y="nodeid"),
  all_node_coords.2, by.x="node2", by.y="nodeid")

PhM1 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2), colour = "black", shape = 21,  size=3.5)+
  theme_bw() 
  
  plot(PhM1)

#Coloured by eco variable 
 
PhM1 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = diet, fill = diet), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 5))+
  theme_bw()
plot(PhM1) 


PhM2 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = dentition, fill = dentition), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 4))+
  theme_bw()
plot(PhM2) 


PhM3 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = FM, fill = FM), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 3))+
  theme_bw()
plot(PhM3) 


PhM4 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = echolocation, fill = echolocation), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 4))+
  theme_bw()
plot(PhM4) 


PhM5 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = habitat, fill = habitat), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 4))+
  theme_bw()
plot(PhM5) 


PhM6 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = freq_type, fill = freq_type), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 7))+
  theme_bw()
plot(PhM6) 


PhM7 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = age, fill = age), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 5))+
  theme_bw()
plot(PhM7) 


PhM8 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = age, fill = dive_type), colour = "black",shape = 21,  size=3)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 5))+
  theme_bw()
plot(PhM8) 


#Morphospace by suborder

PhM9 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = suborder, fill = suborder), colour = "black", shape = 21,  size= 4.5)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 3))+
  theme_bw()
plot(PhM9) 

grid.arrange(PhM7, PhM2, PhM1, PhM8, PhM4, PhM3, PhM6, PhM5, ncol=2)

library(gtable)
library(grid)
g2 <- ggplotGrob(PhM1)
g3 <- ggplotGrob(PhM2)
g4 <- ggplotGrob(PhM3)
g5 <- ggplotGrob(PhM4)
g6 <- ggplotGrob(PhM5)
g7 <- ggplotGrob(PhM6)
g8 <- ggplotGrob(PhM7)
g <- rbind(g2, g3,g4, g5, g6, g7, g8)
g$widths <- unit.pmax(g2$widths, g3$widths, g4$widths,g5$widths, g6$widths, g7$widths, g8$widths)
grid.newpage()
grid.draw(g)


#plot colours, split out families, add phylogeny

# Create a named vector of colors for specific families
my_colors <- c("Kogiidae" = 'green', "Delphinidae" = 'pink', "Phocoenidae" = 'black', 
               "Ziphiidae" = 'purple', "Iniidae" = 'red', "Eurhinodelphinidae" = 'yellow', "Basilosauridae" = 'orange')

# Set the default color for all other families
default_color <- "green"

#plot colours, split out families, 
PhM1 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), linewidth=0.7, alpha=0.5) +
  geom_point(data = PCW,aes(x=Comp1, y = Comp2, color=family, fill=family, shape=suborder), size=6)+
  scale_fill_manual(values = c(my_colors, default_color)) +
  scale_shape_manual(values = c(24, 22, 21)) +
  labs(fill = "suborder", shape= "suborder") +
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide="none")+
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  theme_bw() 

plot(PhM1)
