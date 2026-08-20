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

##notas de outras analises 
#Now, all variables except filamentColor:stamensNumber are numeric. Below, I will transform the quali variables
# into factors (0,1,2 etc)

#Color of filaments
unique(morpho_data_2$filamentColor) # the states are pink, white and yellow. 

morpho_data_2$filamentColor <- as.numeric(factor(morpho_data_2$filamentColor, 
                                                 levels = c('white', 'yellow', 'pink'))) - 1
# Now, for the variable filamentColor, 0 is white, 1 is yellow and 2 is pink

#Number of stamens
unique(morpho_data_2$stamensNumber) # the states are 8, 10 and 4
as.numeric(factor(morpho_data_2$stamensNumber, levels = c('10', '8', '4'))) - 1

#Maybe is better to define them as haplo or diplo 
morpho_data_2 %>%
  mutate(stamensNumber = case_when(
    stamensNumber == numberPetals ~ "haplostemonous",
    stamensNumber == 2 * numberPetals ~ "diplostemonous",
    TRUE ~ as.character(stamensNumber)
  )) -> morpho_data_2

#Now, 0 if they are diplo and 1 if they are haplo
morpho_data_2$stamensNumber <- as.numeric(factor(morpho_data_2$stamensNumber, 
                                                 levels = c('diplostemonous', 'haplostemonous'))) - 1

unique(morpho_data_2$numberPetals) # the states are 4 and 5
morpho_data_2$numberPetals <- as.numeric(factor(morpho_data_2$numberPetals, 
                                                levels = c('4', '5'))) - 1
# Now, for the variable numberPetals, 0 is 4-merous 1 is 5-merous.

head(morpho_data_2) #now, all variables are presented as states. 
