if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(tidyverse, FactoMineR, phytools, data.table,vegan, cluster, readxl, ggplot2)
mimosa_data_quali <- read_excel("C:/Users/Eulália/Desktop/bia/Labis/Dados_treino/mimosa_data_quali.xlsx")
mimosa_quali_sel <- mimosa_data_quali %>% 
  select(infOrganization, infType, infShape, flowerMerosity, calyxShape, corollaShape) %>%
  mutate(across(everything(), as.factor))  

gower_dist <- daisy(mimosa_quali_sel, metric = "gower")

nmds_results <- monoMDS(gower_dist, trymax = 100)

plot(nmds_results)

nmds_points <- as.data.frame(scores(nmds_results))
nmds_points$species <- mimosa_data_quali[[1]]  # nomes da primeira coluna
## Adicionar os nomes das espécies como coluna "especie"
#scores$especie <- rownames(scores)
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color = "red", size = 3) +
  geom_text(aes(label = species), vjust = -0.5, size = 3) +
  theme_minimal()

summary(nmds_results)

#para saber qual variável tem maior peso na contribuição dos 
##componentes principais
nmds_results$stress
install.packages("phytools", "dispRity","FactoMineR", "ape")
library(phytools)
library(FactoMineR)
library(dispRity)
library(ape)
##Phylomorfospace 
tree <- read.nexus("C:/Users/Eulália/Desktop/bia/Labis/Dados_treino/pruned_tree-mimosa.nex")
# Colocar a primeira coluna como rownames
# Extrair coordenadas NMDS + adicionar nomes das espécies
scores_nmds2 <- scores(nmds_results)[, 1:2] %>%
  as.data.frame() %>%
  mutate(species = mimosa_data_quali[[1]])

mimosa_data_quali2 <- mimosa_data_quali
rownames(mimosa_data_quali2) <- mimosa_data_quali2[[1]]

# Verificar nomes que não batem
setdiff(rownames(mimosa_data_quali2), tree$tip.label)  # nos dados, mas não na árvore
setdiff(tree$tip.label, rownames(mimosa_data_quali2))  # na árvore, mas não nos dados
# Podar árvore para manter apenas espécies que existem nos dados
prunedTree <- keep.tip(tree, rownames(mimosa_data_quali2))
scores_for_plot <- scores_nmds2[, 1:2]
rownames(scores_for_plot) <- scores_nmds2$species  # garantir que os nomes batem com a árvore
phylomorphospace(prunedTree, scores_for_plot, label = "horizontal", node.size = c(0, 1))
