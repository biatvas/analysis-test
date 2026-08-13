setwd("C:/Users/Eulália/Desktop/bia/Labis/analysis_test")
if (!require(librarian)) install.packages("librarian"); library("librarian")

library(librarian)
librarian::shelf(tidyverse, FactoMineR, phytools, data.table, cluster, ape)
tree <- read.nexus("pruned_tree-mimosa.nex")
dat_quali <- read.csv2(file="2.dados/mimosa_data_quali.csv")
str(dat_quali)
rownames(dat_quali) <- dat_quali$X
data <- dat_quali[, -1]
data[] <- lapply(data, as.factor)
# Calcular a distância 
gower_dist <- daisy(data, metric = "gower")
# PCoA com correção Cailliez para matrizes não euclidianas
pcoa_result <- pcoa(gower_dist, correction = "cailliez")
summary(pcoa_result)
plot(pcoa_result$vectors[,1], pcoa_result$vectors[,2], xlab="PCoA1", ylab="PCoA2")
#ajustes finais do resultado
text(pcoa_result$vectors[,1], pcoa_result$vectors[,2], 
     labels = rownames(pcoa_result$vectors), 
     pos = 4, cex = 0.7)
#ajustando nomes de colunas
library(ggplot2)

pcoa_df <- as.data.frame(pcoa_result$vectors[,1:2])
pcoa_df$species <- rownames(pcoa_df)

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, label = species)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  labs(x = "PCoA1", y = "PCoA2") +
  theme_minimal()

summary(pcoa_result)
View(pcoa_df)

