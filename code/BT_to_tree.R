
################################################
#                                              #
#  Plotting outputs from BTs onto phylogeny    #
#  Code: Ryan Felice                           #
#                                              #
################################################

install.packages("devtools")
library(devtools)
install_github("hferg/BTprocessR")
library(BTprocessR)

library(viridis)
library(ape)
library(phytools)

mytreebybranch<-function (tree, x, mode = c("edges", "tips", "nodes"), palette = "rainbow", nbin=15, leg.cex=.04,
                          legend = TRUE, xlims = NULL, ...)
{
  mode <- mode[1]
  if (!inherits(tree, "phylo"))
    stop("tree should be an object of class \"phylo\".")
  if (mode == "tips") {
    x <- c(x[tree$tip.label], fastAnc(tree, x))
    names(x)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  else if (mode == "nodes") {
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  if (hasArg(tol))
    tol <- list(...)$tol
  else tol <- 1e-06
  if (hasArg(prompt))
    prompt <- list(...)$prompt
  else prompt <- FALSE
  if (hasArg(type))
    type <- list(...)$type
  else type <- "phylogram"
  if (hasArg(show.tip.label))
    show.tip.label <- list(...)$show.tip.label
  else show.tip.label <- TRUE
  if (hasArg(show.node.label))
    show.node.label <- list(...)$show.node.label
  else show.node.label <- FALSE
  if (hasArg(edge.width))
    edge.width <- list(...)$edge.width
  else edge.width <- 4
  if (hasArg(edge.lty))
    edge.lty <- list(...)$edge.lty
  else edge.lty <- 1
  if (hasArg(font))
    font <- list(...)$font
  else font <- 3
  if (hasArg(cex))
    cex <- list(...)$cex
  else cex <- par("cex")
  if (hasArg(adj))
    adj <- list(...)$adj
  else adj <- NULL
  if (hasArg(srt))
    srt <- list(...)$srt
  else srt <- 0
  if (hasArg(no.margin))
    no.margin <- list(...)$no.margin
  else no.margin <- TRUE
  if (hasArg(root.edge))
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (hasArg(label.offset))
    label.offset <- list(...)$label.offset
  else label.offset <- 0.01 * max(nodeHeights(tree))
  if (hasArg(underscore))
    underscore <- list(...)$underscore
  else underscore <- FALSE
  if (hasArg(x.lim))
    x.lim <- list(...)$x.lim
  else x.lim <- NULL
  if (hasArg(y.lim))
    y.lim <- list(...)$y.lim
  else y.lim <- if (legend && !prompt && type %in% c("phylogram",
                                                     "cladogram"))
    c(1 - 0.06 * length(tree$tip.label), length(tree$tip.label))
  else NULL
  if (hasArg(direction))
    direction <- list(...)$direction
  else direction <- "rightwards"
  if (hasArg(lab4ut))
    lab4ut <- list(...)$lab4ut
  else lab4ut <- NULL
  if (hasArg(tip.color))
    tip.color <- list(...)$tip.color
  else tip.color <- "black"
  if (hasArg(plot))
    plot <- list(...)$plot
  else plot <- TRUE
  if (hasArg(rotate.tree))
    rotate.tree <- list(...)$rotate.tree
  else rotate.tree <- 0
  if (hasArg(open.angle))
    open.angle <- list(...)$open.angle
  else open.angle <- 0
  if (is.function(palette))
    cols <- palette(n = nbin)
  else {
    if (palette == "heat.colors")
      cols <- heat.colors(n = nbin)
    if (palette == "gray")
      cols <- gray(nbin:1/nbin)
    if (palette == "rainbow")
      cols <- rainbow(nbin, start = 0.7, end = 0)
  }
  if (is.null(xlims))
    xlims <- range(x) + c(-tol, tol)
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  par(lend = 2)
  xx <- plot.phylo(tree, type = type, show.tip.label = show.tip.label,
                   show.node.label = show.node.label, edge.color = colors,
                   edge.width = edge.width, edge.lty = edge.lty, font = font,
                   cex = cex, adj = adj, srt = srt, no.margin = no.margin,
                   root.edge = root.edge, label.offset = label.offset, underscore = underscore,
                   x.lim = x.lim, y.lim = y.lim, direction = direction,
                   lab4ut = lab4ut, tip.color = tip.color, plot = plot,
                   rotate.tree = rotate.tree, open.angle = open.angle, lend = 2,
                   new = FALSE)
  if (legend == TRUE && is.logical(legend))
    legend <- round(0.5 * max(nodeHeights(tree)), 2)
  if (legend) {
    if (hasArg(title))
      title <- list(...)$title
    else title <- "trait value"
    if (hasArg(digits))
      digits <- list(...)$digits
    else digits <- 1
    if (prompt)
      add.color.bar(legend, cols, title, xlims, digits, fsize=leg.cex,
                    prompt = TRUE)
    else add.color.bar(legend, cols, title, xlims, digits, fsize=leg.cex,
                       prompt = FALSE, 
                       x = par()$usr[1] + 0.05 * (par()$usr[2] - par()$usr[1]), 
                       y = par()$usr[3] + 0.05 * (par()$usr[4] - par()$usr[3])) }
  invisible(xx)
}



#Plot your tree 

BT_rates_tree <-
  BTprocessR::rjpp(
    # rates results (.VarRates.txt file)
    rjlog = 
      "lambda_var_whole_var.txt.VarRates.txt",
    # posterior trees (.Output.trees file)
    rjtrees = 
      "lambda_var_whole_var.txt.Output.trees",
    # tree that you used for the BayesTrauts input
    tree = treeMAND2
  )


color4<-colorRampPalette(viridis(n = 100 , option = "plasma")[-c(1:10)])


mytreebybranch(BT_rates_tree$tree_summary$tree_summaries$original_tree,
               title="log(Mean Rate)", 
               x=log(BT_rates_tree$scalars$meanRate[-1]),
               mode="edges",
               type="fan",
               leg.cex=.6,
               cex=c(.3),
               edge.width=1.5,
               palette=color4,
               nbin=6,
               no.margin=FALSE)


#Andy's code for standardising the scale 

#create function for ordering all values from 0 to 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# work out range of all branch rates and save them as a single vector, running 0-1 (change BT_rates_tree to match your data)
ratesRange <-range01(c(BT_rates_tree1$scalars$meanRate,BT_rates_tree2$scalars$meanRate))

ratesRange2 <-(ratesRange*100)+1

# plot two trees, in same order as worked out above, using first half of branch rates vector to colour branches
mytreebybranch(BT_rates_tree1$tree_summary$tree_summaries$original_tree,
               title="log(Mean Rate)",
               x=log(ratesRange2[1:(length(ratesRange2)/2)][-1]),
               mode="edges",
               type="fan",
               leg.cex=0.5,
               cex=c(0.75),
               edge.width=3,
               palette=color4,
               nbin=6,
               no.margin=FALSE)

# plot second tree, changing value for x to incorporate second half of branch rates vector
mytreebybranch(BT_rates_tree2$tree_summary$tree_summaries$original_tree,
               title="log(Mean Rate)",
               x=log(ratesRange2[(length(ratesRange2)/2):length(ratesRange2)][-1]),
               mode="edges",
               type="fan",
               leg.cex=0.5,
               cex=c(0.75),
               edge.width=3,
               palette=color4,
               nbin=6,
               no.margin=FALSE)
