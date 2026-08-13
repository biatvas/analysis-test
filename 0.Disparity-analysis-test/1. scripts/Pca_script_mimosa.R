setwd("C:/Users/Eulália/Desktop/bia/Labis/Dados_treino")
if (!require(librarian)) install.packages("librarian"); library("librarian")

librarian::shelf(tidyverse, FactoMineR, phytools, data.table)
dat_flower <- read.csv2(file="mimosa_data_example.csv")
str(dat_flower)
rownames(dat_flower) <- dat_flower[,1]
dat_flower <- dat_flower[,-1]
# Nomes dos dados
data_names <- rownames(dat_flower)
# Corrigir dados que são texto
dat_flower[] <- lapply(dat_flower, function(x) as.numeric(gsub(",", ".", x)))
# Padronizar
dados_padronizados <- log(dat_flower) 
# a função scale gera um dado ordenado para informações que tem escalas diferentes 
# PCA
pca_resultado <- prcomp(dados_padronizados)
summary(pca_resultado)
biplot(pca_resultado)
scores <- as.data.frame(pca_resultado$x[, 1:2])
# Adicionar os nomes das espécies como coluna "especie"
scores$especie <- rownames(scores)
# Verifique se a coluna foi criada corretamente
head(scores)
library(ggplot2)
ggplot(scores, aes(x = PC1, y = PC2, label = especie)) +
  geom_point(color = "red") +
  geom_text(vjust = -0.8, size = 3) +
  theme_minimal() +
  labs(title = "PCA - Espécies de Mimosa",
       x = "PC 1",
       y = "PC 2")
scores <- pca_resultado$x[, 1:2]

install.packages("factorextra")
library("factorextra")
#fviz_contrib(pca_resultado, choice = "ind", axes = 1:2)
## Phylomorphospace ========================================================
#prune tree phylogeny
all(rownames(scores) %in% prunedTree$tip.label)
library(phytools)
prunedTree <- keep.tip(tree, rownames(dat_flower))
tree <- read.nexus("pruned_tree-mimosa.nex")
phylomorphospace(prunedTree, scores, label = "horizontal", node.size = c(0, 1))

