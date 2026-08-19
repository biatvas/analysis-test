#============#
#### HEAD ####
#============#

setwd("C:/Users/moniq/OneDrive/doctorate/data/disparity/disparityAnalyzes")
setwd("~/Library/CloudStorage/OneDrive-Personal/doctorate/data/disparity/disparityAnalyzes") #mac

if (!require(librarian)) install.packages("librarian"); library("librarian")

librarian::shelf(tidyverse, FactoMineR, phytools, data.table)

tree <- read.tree("dataset/mimosa_clean-VASCONCELOS2020.txt")
treeData <- read.csv("dataset/treeData.csv", header = T)
phyloData <- read.csv("dataset/Mimosa_data_updated_140520-VASCONCELOS2020 - Mimosa_data_updated_140520-VASCONCELOS2020.csv")
#listCharac_old <- readRDS(file="list/listCharac.RData")
listCharac <- readRDS(file="list/listCharac_20_05_2025.RData")
prunedTree <- keep.tip(tree, row.names(listCharac$flower))

#===========#
#### PCA ####
#===========#

## PCA data
dataPCA <- list(flower = PCA(listCharac$flower, scale.unit=F),
                calyx = PCA(listCharac$calyx, scale.unit = F),
                corolla = PCA(listCharac$corolla, scale.unit = F),
                androecium = PCA(listCharac$androecium, scale.unit = F),
                gynoecium = PCA(listCharac$gynoecium, scale.unit = F)
                )
#saveRDS(dataPCA, "list/results/pca/dataPCA.RDS")

row.names(dataPCA$flower$var$coord) <- c("coW","pW","pL","pTL",
                                  "caW","spW","spTL","spL",
                                  "fL","aW","aL",
                                  "oL","oW","sL")

colGroups <- rep(0, 14)
colGroups[row.names(dataPCA$flower$var$coord) %in% c("caW","spW","spTL","spL")] <- "#F8B621" #"calyx  
colGroups[row.names(dataPCA$flower$var$coord) %in% c("coW","pW","pL","pTL")] <- "#1BA3C6" #corolla 
colGroups[row.names(dataPCA$flower$var$coord) %in% c("fL","aW","aL")] <- "#FC719E" #"androecium
colGroups[row.names(dataPCA$flower$var$coord) %in% c("oL","oW","sL")] <- "#33A65C" #"gynoecium

variablesContribution <- plot(dataPCA$flower, choix = "var", title = "", col.var = colGroups) +
  labs(x = "PC1 (54.24%)", y = "PC2 (19.73%)") +
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid = element_blank()) +   
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)

variablesContribution_PC3 <- plot(dataPCA$flower, choix = "var", title = "", 
                                  col.var = colGroups, axes = 2:3) +
  labs(x = "PC2 (19.73%)", y = "PC3 (10.13%)") +  
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid = element_blank()) + 
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)

variables_contribution_PCA <- as.data.frame(dataPCA$flower$var$contrib)

variables_contribution_PCA[nrow(variables_contribution_PCA) + 1, ] =
  as.numeric(format(dataPCA$flower$eig[1:5, 2], scientific = FALSE, digits = 3))

variables_contribution_PCA <- format(variables_contribution_PCA, digits = 2)

write.csv2(variables_contribution_PCA, "list/results/pca_and_phylomorphospace/variables_contribution.csv")

## Total phylomorphospace
#phylomorphospace(prunedTree, dataEigs[,1:2], node.size=c(0,1.2))
tip.cols<-rep(c("black"),Ntip(prunedTree))
names(tip.cols) <- prunedTree$tip.label
cols <- c(tip.cols[prunedTree$tip.label],rep("white",prunedTree$Nnode))
names(cols)<-1:(length(prunedTree$tip)+prunedTree$Nnode)

#pdf("figures/phylomorphospace/phylo_plot_total_2.pdf",
#    width = 8.24,
#    height = 6.18,
#    paper = "a4",
#    family="ArialMT",
#    bg = "transparent") #6.18 inches by 8.24 inches 

legend <- c(order(abs(dataPCA$flower$ind$coord[,1]), decreasing = TRUE)[1:10], 
            order(abs(dataPCA$flower$ind$coord[,2]), decreasing = TRUE)[1:10])

plot(dataPCA$flower)
phylomorphospace(prunedTree, dataPCA$flower$ind$coord[,1:2], xlab = "PC 1 (54.24%)", 
                 ylab = "PC 2 (19.73%)", label = "off",
                 control=list(col.node=cols),
                 lwd = 2,
                 node.size = c(0.5,1.25), 
                 node.by.map=T) -> phylomorphospace

# grid(lty = 1, col = "darkgray", lwd = 1)
text(dataPCA$flower$ind$coord[legend, 1], dataPCA$flower$ind$coord[legend, 2], 
     labels = rownames(dataPCA$flower$ind$coord)[legend], 
     pos = 3, offset = 0.5, cex = 0.7, col = "black")
#dev.off()

legend <- c(order(abs(dataPCA$flower$ind$coord[,2]), decreasing = TRUE)[1:10], 
            order(abs(dataPCA$flower$ind$coord[,3]), decreasing = TRUE)[1:10])

phylomorphospace(prunedTree, dataPCA$flower$ind$coord[,2:3], xlab = "PC 2 PC 2 (19.73%)", 
                 ylab = "PC 3 ", label = "off",
                 control=list(col.node=cols),
                 lwd = 2,
                 node.size = c(0.5,1.25), 
                 node.by.map=T) -> phylomorphospace

# grid(lty = 1, col = "darkgray", lwd = 1)
text(dataPCA$flower$ind$coord[legend, 2], dataPCA$flower$ind$coord[legend, 3], 
     labels = rownames(dataPCA$flower$ind$coord)[legend], 
     pos = 3, offset = 0.5, cex = 0.7, col = "black")
#dev.off()

## Phylomorphospace calyx
tip.cols_calyx<-rep(c("#F8B621"),Ntip(prunedTree))
names(tip.cols_calyx) <- prunedTree$tip.label
cols_calyx <- c(tip.cols_calyx[prunedTree$tip.label],rep("white",prunedTree$Nnode))
names(cols_calyx)<-1:(length(prunedTree$tip)+prunedTree$Nnode)

legend_calyx <- c(order(abs(dataPCA$calyx$ind$coord[,1]), decreasing = TRUE)[1:10], 
                  order(abs(dataPCA$calyx$ind$coord[,2]), decreasing = TRUE)[1:10])

plot(dataPCA$calyx)

phylomorphospace(prunedTree, dataPCA$calyx$ind$coord[,1:2], xlab = "PC 1 (87.56%)", 
                 ylab = "PC 2 (9.95%)", label = "off",
                 control=list(col.node=cols_calyx),
                 lwd = 2,
                 node.size = c(0.5,1.25), 
                 node.by.map=TRUE) -> phylomorphospace_calyx

text(dataPCA$calyx$ind$coord[legend_calyx, 1], dataPCA$calyx$ind$coord[legend_calyx, 2], 
     labels = rownames(dataPCA$calyx$ind$coord)[legend_calyx], 
     pos = 3, offset = 0.5, cex = 0.7, col = "black")

## Phylomorphospace corolla
tip.cols_corolla<-rep(c("#1BA3C6"),Ntip(prunedTree))
names(tip.cols_corolla) <- prunedTree$tip.label
cols_corolla <- c(tip.cols_corolla[prunedTree$tip.label],rep("white",prunedTree$Nnode))
names(cols_corolla)<-1:(length(prunedTree$tip)+prunedTree$Nnode)

legend_corolla <- c(order(abs(dataPCA$corolla$ind$coord[,1]), decreasing = TRUE)[1:10], 
                    order(abs(dataPCA$corolla$ind$coord[,2]), decreasing = TRUE)[1:10])

plot(dataPCA$corolla)
phylomorphospace(prunedTree, dataPCA$corolla$ind$coord[,1:2], xlab = "PC 1 (69.13%)", 
                 ylab = "PC 2 (27.88%)", label = "off",
                 control=list(col.node=cols_corolla),
                 lwd = 2,
                 node.size = c(0.5,1.25), node.by.map=TRUE) -> phylomorphospace_corolla
  #theme(panel.grid.major = element_blank(),
  #      axis.title = element_text(size = 1, color = "black"),
  #      panel.grid = element_blank()) +
  text(dataPCA$corolla$ind$coord[legend_corolla, 1], dataPCA$corolla$ind$coord[legend_corolla, 2], 
       labels = rownames(dataPCA$corolla$ind$coord)[legend_corolla], 
       pos = 3, offset = 0.5, cex = 0.7, col = "black")

## Phylomorphospace androecium
tip.cols_androecium<-rep(c("#FC719E"),Ntip(prunedTree))
names(tip.cols_androecium) <- prunedTree$tip.label
cols_androecium <- c(tip.cols_androecium[prunedTree$tip.label],rep("white",prunedTree$Nnode))
names(cols_androecium)<-1:(length(prunedTree$tip)+prunedTree$Nnode)

legend_androecium <- c(order(abs(dataPCA$androecium$ind$coord[,1]), decreasing = TRUE)[1:10], 
                       order(abs(dataPCA$androecium$ind$coord[,2]), decreasing = TRUE)[1:10])
plot(dataPCA$androecium)
phylomorphospace(prunedTree, dataPCA$androecium$ind$coord[,1:2], xlab = "PC 1 (77.08%)", 
                 ylab = "PC 2 (18.07%)", label = "off",
                 control=list(col.node=cols_androecium),
                 lwd = 2,
                 node.size = c(0.5,1.5), node.by.map=TRUE) -> phylomorphospace_androecium
 # theme(panel.grid.major = element_blank(),
 #       axis.title = element_text(size = 1, color = "black"),
 #       panel.grid = element_blank()) +
  text(dataPCA$androecium$ind$coord[legend_androecium, 1], dataPCA$androecium$ind$coord[legend_androecium, 2], 
      labels = rownames(dataPCA$androecium$ind$coord)[legend_androecium], 
      pos = 3, offset = 0.5, cex = 0.7, col = "black")

## Phylomorphospace gynoecium
tip.cols_gynoecium<-rep(c("#33A65C"),Ntip(prunedTree))
names(tip.cols_gynoecium) <- prunedTree$tip.label
cols_gynoecium <- c(tip.cols_gynoecium[prunedTree$tip.label],rep("white",prunedTree$Nnode))
names(cols_gynoecium)<-1:(length(prunedTree$tip)+prunedTree$Nnode)
plot(dataPCA$gynoecium)
legend_gynoecium <- c(order(abs(dataPCA$gynoecium$ind$coord[,1]), decreasing = TRUE)[1:10], 
                      order(abs(dataPCA$gynoecium$ind$coord[,2]), decreasing = TRUE)[1:10])

phylomorphospace(prunedTree, dataPCA$gynoecium$ind$coord[,1:2], xlab = "PC 1 (62.98%)", 
                 ylab = "PC 2 (31.06%)", label = "off",
                 control=list(col.node=cols_gynoecium),
                 lwd = 2,
                 node.size = c(0.5,1.5), node.by.map=TRUE) -> phylomorphospace_gynoecium
  #theme(panel.grid.major = element_blank(),
  #      axis.title = element_text(size = 1, color = "black"),
  #      panel.grid = element_blank()) +
  text(dataPCA$gynoecium$ind$coord[legend_gynoecium, 1], dataPCA$gynoecium$ind$coord[legend_gynoecium, 2], 
       labels = rownames(dataPCA$gynoecium$ind$coord)[legend_gynoecium], 
       pos = 3, offset = 0.5, cex = 0.7, col = "black")


#### Phylomorphospace with ggplot ####

phylo_data <- list(
  total = data.frame(
    xstart = phylomorphospace$xx[phylomorphospace$edge[, 1]],
    ystart = phylomorphospace$yy[phylomorphospace$edge[, 1]],
    xstop = phylomorphospace$xx[phylomorphospace$edge[, 2]],
    ystop = phylomorphospace$yy[phylomorphospace$edge[, 2]],
    nodestart = phylomorphospace$edge[, 1],
    nodestop = phylomorphospace$edge[, 2]),
  calyx = data.frame(
    xstart = phylomorphospace_calyx$xx[phylomorphospace_calyx$edge[, 1]],
    ystart = phylomorphospace_calyx$yy[phylomorphospace_calyx$edge[, 1]],
    xstop = phylomorphospace_calyx$xx[phylomorphospace_calyx$edge[, 2]],
    ystop = phylomorphospace_calyx$yy[phylomorphospace_calyx$edge[, 2]],
    nodestart = phylomorphospace_calyx$edge[, 1],
    nodestop = phylomorphospace_calyx$edge[, 2]),
  corolla = data.frame(
    xstart = phylomorphospace_corolla$xx[phylomorphospace_corolla$edge[, 1]],
    ystart = phylomorphospace_corolla$yy[phylomorphospace_corolla$edge[, 1]],
    xstop = phylomorphospace_corolla$xx[phylomorphospace_corolla$edge[, 2]],
    ystop = phylomorphospace_corolla$yy[phylomorphospace_corolla$edge[, 2]],
    nodestart = phylomorphospace_corolla$edge[, 1],
    nodestop = phylomorphospace_corolla$edge[, 2]),
  androecium = data.frame(
    xstart = phylomorphospace_androecium$xx[phylomorphospace_androecium$edge[, 1]],
    ystart = phylomorphospace_androecium$yy[phylomorphospace_androecium$edge[, 1]],
    xstop = phylomorphospace_androecium$xx[phylomorphospace_androecium$edge[, 2]],
    ystop = phylomorphospace_androecium$yy[phylomorphospace_androecium$edge[, 2]],
    nodestart = phylomorphospace_androecium$edge[, 1],
    nodestop = phylomorphospace_androecium$edge[, 2]),
  gynoecium = data.frame(
    xstart = phylomorphospace_gynoecium$xx[phylomorphospace_gynoecium$edge[, 1]],
    ystart = phylomorphospace_gynoecium$yy[phylomorphospace_gynoecium$edge[, 1]],
    xstop = phylomorphospace_gynoecium$xx[phylomorphospace_gynoecium$edge[, 2]],
    ystop = phylomorphospace_gynoecium$yy[phylomorphospace_gynoecium$edge[, 2]],
    nodestart = phylomorphospace_gynoecium$edge[, 1],
    nodestop = phylomorphospace_gynoecium$edge[, 2]))

node_coords <- list(
  total = data.frame(x = phylomorphospace$xx[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                     y = phylomorphospace$yy[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                     node = (Ntip(prunedTree) + 1):(Ntip(prunedTree) +  prunedTree$Nnode)),
  calyx = data.frame(x = phylomorphospace_calyx$xx[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                     y = phylomorphospace_calyx$yy[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                     node = (Ntip(prunedTree) + 1):(Ntip(prunedTree) +  prunedTree$Nnode)),
  corolla = data.frame(x = phylomorphospace_corolla$xx[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                       y = phylomorphospace_corolla$yy[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                       node = (Ntip(prunedTree) + 1):(Ntip(prunedTree) +  prunedTree$Nnode)),
  androecium = data.frame(x = phylomorphospace_androecium$xx[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                          y = phylomorphospace_androecium$yy[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                          node = (Ntip(prunedTree) + 1):(Ntip(prunedTree) +  prunedTree$Nnode)),
  gynoecium = data.frame(x = phylomorphospace_gynoecium$xx[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                         y = phylomorphospace_gynoecium$yy[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                         node = (Ntip(prunedTree) + 1):(Ntip(prunedTree) +  prunedTree$Nnode))
)

library(ggtree)
phylomorphospace_total <-
  ggplot() +
  geom_segment(
    data = phylo_data$total,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop),
    linewidth = 0.25,
    color = "black"
  ) +
  geom_point(data = node_coords$total,
             aes(x = x, y = y),
             shape = 21,
             fill = "white",
             color = "black",
             size = 0.07) +
  geom_point(data = dataPCA$flower$ind$coord,
             aes(
               x = Dim.1,
               y = Dim.2,
               label = row.names(dataPCA$flower$ind$coord)),
             size = 3) +
  coord_fixed(
    ratio = 1,
    expand = TRUE,    
    xlim = c(-3.5,3),
    ylim = c(-2,2),
    clip = "on") +
  labs(x = "PC1 (54.24%)", y = "PC2 (19.74%)") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

phylomorphospace_total

#interactive phylomorphospace. 
#ggplotly(phylomorphospace_total, tooltip = "label") #animated phylomorphospace

phylomorphospace_calyx_2 <-
  ggplot() +
  geom_segment(
    data = phylo_data$calyx,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop),
    linewidth = 0.25,
    color = "#F8B621",
    alpha = 0.4
  ) +
  geom_point(data = node_coords$calyx,
             aes(x = x, y = y),
             shape = 21,
             fill = "white",
             color = "#F8B621",
             size = 0.03) +
  geom_point(data = dataPCA$calyx$ind$coord,
             aes(
               x = Dim.1,
               y = Dim.2,
               label = row.names(dataPCA$calyx$ind$coord)),
             shape = 21,
             fill = "#F8B621",
             color = "#F8B621",
             size = 2) +
  coord_fixed(
    ratio = 1,
    #expand = TRUE,
    xlim = c(-3,2.5),
    ylim = c(-1.5,1.5),
    clip = "on") +
  labs(x = "PC1 (87.56%)", y = "PC2 (9.95%)") +
  theme_bw() +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

phylomorphospace_calyx_2

phylomorphospace_corolla_2 <- 
  ggplot() +
  geom_segment(
    data = phylo_data$corolla,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop),
    linewidth = 0.25,
    color = "#1BA3C6",
    alpha = 0.4) +
  geom_point(data = node_coords$corolla,
             aes(x = x, y = y),
             shape = 21,
             fill = "white",
             color = "#1BA3C6",
             size = 0.03) +
  geom_point(data = dataPCA$corolla$ind$coord,
             aes(
               x = Dim.1,
               y = Dim.2,
               color = "#1BA3C6",
               label = row.names(dataPCA$corolla$ind$coord)),
             fill = "#1BA3C6",
             color = "#1BA3C6",
             size = 2) +
  coord_fixed(
    ratio = 1,
    expand = TRUE,
    xlim = c(-3,2.5),
    ylim = c(-1.5,1.5),
    clip = "on") +
  labs(x = "PC1 (69.13%)", y = "PC2 (27.88%)") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

phylomorphospace_corolla_2

## Androecium
phylomorphospace_androecium_2 <-
  ggplot() +
  geom_segment(
    data = phylo_data$androecium,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop),
    linewidth = 0.25,
    color = "#FC719E",
    alpha = 0.4
  ) +
  geom_point(data = node_coords$androecium,
             aes(x = x, y = y),
             shape = 21,
             fill = "white",
             color = "#FC719E",
             size = 0.03) +
  geom_point(data = dataPCA$androecium$ind$coord,
             aes(
               x = Dim.1,
               y = Dim.2,
               label = row.names(dataPCA$androecium$ind$coord)),
             fill = "#FC719E",
             color = "#FC719E",
             size = 2) +
  coord_fixed(
    ratio = 1,
    expand = TRUE,
    xlim = c(-3,2.5),
    ylim = c(-1.5,1.5),
    clip = "on") +
  labs(x = "PC1 (77.08%)", y = "PC2 (18.07%)") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

phylomorphospace_androecium_2

## Gynoecium
phylomorphospace_gynoecium_2 <-
  ggplot() +
  geom_segment(
    data = phylo_data$gynoecium,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop),
    linewidth = 0.25,
    color = "#33A65C",
    alpha = 0.4
  ) +
  geom_point(data = node_coords$gynoecium,
             aes(x = x, y = y),
             shape = 21,
             fill = "white",
             color = "#33A65C",
             size = 0.03) +
  geom_point(data = dataPCA$gynoecium$ind$coord,
             aes(
               x = Dim.1,
               y = Dim.2,
               label = row.names(dataPCA$gynoecium$ind$coord)),
             fill = "#33A65C",
             color = "#33A65C",
             size = 2) +
  coord_fixed(
    ratio = 1,
    expand = TRUE,
    xlim = c(-3,2.5),
    ylim = c(-1.5,1.5),
    clip = "on") +
  labs(x = "PC1 (62.97%)", y = "PC2 (31.07%)") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

phylomorphospace_gynoecium_2

pdf("figures/phylomorphospace/phylomorphospace_total.pdf",
    width = 7, height = 7)
plot(phylomorphospace_total)
dev.off()

pdf("figures/phylomorphospace/phylomorphospace_calyx.pdf",
    width = 7, height = 7)
plot(phylomorphospace_calyx_2)
dev.off()

pdf("figures/phylomorphospace/phylomorphospace_corolla.pdf",
    width = 7, height = 7)
plot(phylomorphospace_corolla_2)
dev.off()

pdf("figures/phylomorphospace/phylomorphospace_androecium.pdf",
    width = 7, height = 7)
plot(phylomorphospace_androecium_2)
dev.off()

pdf("figures/phylomorphospace/phylomorphospace_gynoecium.pdf",
    width = 7, height = 7)
plot(phylomorphospace_gynoecium_2)
dev.off()

library(ggpubr)
phylomorphospace_whorls <- ggarrange(phylomorphospace_calyx_2, phylomorphospace_corolla_2, 
                                     phylomorphospace_androecium_2, phylomorphospace_gynoecium_2,
                                     ncol = 2, nrow = 2)


pdf("figures/phylomorphospace/phylomorphospace_whorls.pdf",
    width = 7, height = 7)
plot(phylomorphospace_whorls)
dev.off()

variables_contribution_all <- ggarrange(variablesContribution, variablesContribution_PC3, ncol = 2, nrow = 1)

pdf("figures/phylomorphospace/variables_contribution_all.pdf",
    width = 7, height = 7)
plot(variables_contribution_all)
dev.off()

#============================#
##### Ploting PC2 and PC3 ####

rm(list = setdiff(ls(), c("dataPCA", "prunedTree")))

## Total phylomorphospace

tip.cols<-rep(c("black"),Ntip(prunedTree))
names(tip.cols) <- prunedTree$tip.label
cols <- c(tip.cols[prunedTree$tip.label],rep("white",prunedTree$Nnode))
names(cols)<-1:(length(prunedTree$tip)+prunedTree$Nnode)

#pdf("figures/phylomorphospace/phylo_plot_total_2.pdf",
#    width = 8.24,
#    height = 6.18,
#    paper = "a4",
#    family="ArialMT",
#    bg = "transparent") #6.18 inches by 8.24 inches 

legend <- c(order(abs(dataPCA$flower$ind$coord[,2]), decreasing = TRUE)[1:10], 
            order(abs(dataPCA$flower$ind$coord[,3]), decreasing = TRUE)[1:10])

plot(dataPCA$flower, axes = 2:3)
phylomorphospace(prunedTree, dataPCA$flower$ind$coord[,2:3], xlab = "PC 2 (19.73%)", 
                 ylab = "PC 3 (10.13%)", label = "off",
                 control=list(col.node=cols),
                 lwd = 2,
                 node.size = c(0.5,1.25), 
                 node.by.map=T) -> phylomorphospace

# grid(lty = 1, col = "darkgray", lwd = 1)
text(dataPCA$flower$ind$coord[legend, 2], dataPCA$flower$ind$coord[legend, 3], 
     labels = rownames(dataPCA$flower$ind$coord)[legend], 
     pos = 3, offset = 0.5, cex = 0.7, col = "black")
#dev.off()

#### Phylomorphospace with ggplot ####

phylo_data <- list(
  total = data.frame(
    xstart = phylomorphospace$xx[phylomorphospace$edge[, 1]],
    ystart = phylomorphospace$yy[phylomorphospace$edge[, 1]],
    xstop = phylomorphospace$xx[phylomorphospace$edge[, 2]],
    ystop = phylomorphospace$yy[phylomorphospace$edge[, 2]],
    nodestart = phylomorphospace$edge[, 1],
    nodestop = phylomorphospace$edge[, 2]))

node_coords <- list(
  total = data.frame(x = phylomorphospace$xx[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                     y = phylomorphospace$yy[(Ntip(prunedTree) + 1):(Ntip(prunedTree) + prunedTree$Nnode)],
                     node = (Ntip(prunedTree) + 1):(Ntip(prunedTree) +  prunedTree$Nnode)))

library(ggtree)
phylomorphospace_total_2 <-
  ggplot() +
  geom_segment(
    data = phylo_data$total,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop),
    linewidth = 0.25,
    color = "black"
  ) +
  geom_point(data = node_coords$total,
             aes(x = x, y = y),
             shape = 21,
             fill = "white",
             color = "black",
             size = 0.07) +
  geom_point(data = dataPCA$flower$ind$coord,
             aes(
               x = Dim.2,
               y = Dim.3,
               label = row.names(dataPCA$flower$ind$coord)),
             size = 3) +
  coord_fixed(
    ratio = 1,
    expand = TRUE,    
    xlim = c(-2,2),
    ylim = c(-1.5,1.5),
    clip = "on") +
  labs(x = "PC2 (19.73%)", y = "PC3 (10.13%)") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

phylomorphospace_total_2


pdf("figures/phylomorphospace/phylomorphospace_total_PC2_and_PC3.pdf",
    width = 7, height = 7)
plot(phylomorphospace_total_2)
dev.off()