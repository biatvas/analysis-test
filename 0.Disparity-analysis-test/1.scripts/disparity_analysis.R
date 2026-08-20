### download packages =====
library(readr)
library(stats)
library(dispRity)
library(phytools)
library(ggplot2)
library(vegan)
library(FactoMineR) #PCA
library(ade4) #tem a funcao que trabalha com dados mistos!
library(dplyr)
library(cluster)
library(ape)

#read morphological data 
mimosa_data <- read.csv("Documents/GitHub/analysis-test/0.Disparity-analysis-test/2.dados/mimosa_data_example.csv", sep = ";")

rownames(mimosa_data) = mimosa_data$X
traits <- mimosa_data[, 2:15]

#distancia entre taxon
matriz_distancia <- dist(traits, method = "euclidean")
print(matriz_distancia)

matriz_distancia_objeto <- as.matrix(matriz_distancia)
print(matriz_distancia_objeto)


#plot PCA for quantitative data 
pca <- cmdscale(matriz_distancia, k = 14, eig = TRUE)
pca_values <- pca$values
pca_points <- pca$points

plot(pca_points, type = "p", main = "PCA - Mimosa") +
  geom_point(color = "red") +
  labs(x = "PC2", y = "PC1") +
  xlim(-0.5, 0.5) +
  ylim(-0.5,0.5)
  
# padronização dos dados
dados_padronizados <- scale(traits) #ou uso log?
str(dados_padronizados)

#matriz de correlação dos dados
matriz_cor <- cor(dados_padronizados)
print(matriz_cor)

# PCA com prcomp =======
# pca dos dados
pca_resultado <- prcomp(dados_padronizados)
summary(pca_resultado)
pca_resultado$rotation  # variáveis por componente
#para ver o resultado
pca_resultado$rotation[order(pca_resultado$rotation[,1],decreasing = T),]
pca_resultado$sdev

biplot(pca_resultado, scale = 0)
# Extrai os scores (coordenadas dos indivíduos no novo espaço)
scores <- as.data.frame(pca_resultado$x)
# Adiciona a coluna com os nomes das espécies
scores$especie <- mimosa_data$X
summary(pca_resultado)

ggplot(scores, aes(x = PC1, y = PC2, label = especie )) +
  geom_point(color = "red") +
  geom_text(vjust = -0.5, size = 2) +
  theme_minimal() +
  labs(title = "PCA - Espécies de Mimosa",
       x = "PC 1",
       y = "PC 2")

#para saber qual variável tem maior peso na contribuição dos 
##componentes principais

help("PCA")
#PCA(X, scale.unit = TRUE, ncp = 5, ind.sup = NULL, quanti.sup = NULL, 
# quali.sup = NULL, row.w = NULL, col.w = NULL, graph = TRUE, axes = c(1,2))

## PCoA =====
dat_quali <- read.csv("Documents/GitHub/analysis-test/0.Disparity-analysis-test/2.dados/mimosa_data_quali.csv",
                                           sep = ";")
rownames(dat_quali) = dat_quali$species
traits <- dat_quali[, 2:21]

#definir categorias dos dados quali? e talvez avaliar como fazer input de dados

str(dat_quali)
data <- dat_quali[, -1] #exclui a coluna com nomes de especies
data[] <- lapply(data, as.factor)

# Calcular a distância 
gower_dist <- daisy(data, metric = "gower")

# PCoA com correção Cailliez para matrizes não euclidianas
pcoa_result <- pcoa(gower_dist, correction = "cailliez")
summary(pcoa_result)
pcoa_result$values

#plot simples
plot(pcoa_result$vectors[,1], pcoa_result$vectors[,2], xlab="PCoA1", ylab="PCoA2")

text(pcoa_result$vectors[,1], pcoa_result$vectors[,2], 
     labels = rownames(pcoa_result$vectors), 
     pos = 4, cex = 0.7) 

#ajustando nomes de colunas
pcoa_df <- as.data.frame(pcoa_result$vectors[,1:2])
pcoa_df$species <- rownames(pcoa_df)

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, label = species)) +
  geom_point(color = "yellow3") +
  geom_text(vjust = -0.5, size = 3) +
  labs(x = "PCoA1", y = "PCoA2") +
  theme_minimal() +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)

summary(pcoa_result)
View(pcoa_df)
 
# phylomorphospace
tree <- read.nexus("Documents/GitHub/analysis-test/0.Disparity-analysis-test/4.trees/pruned_tree-mimosa.nex")

# Verificar nomes
setdiff(rownames(data), tree$tip.label)  # nos dados, mas não na árvore
setdiff(tree$tip.label, rownames(data))
prunedTree <- keep.tip(tree, rownames(data))

phylomorphospace(prunedTree, pcoa_df[,1:2], label = "horizontal", 
                 xlab="PCoA2",
                 ylab="PCoA1",
                 node.size = c(0.5, 1.25),
                 node.by.map=T)
title(main="Phylomorphospace of flower morphology in Mimosa", font.main=3)


#### ainda nao conferi essa etapa do script ======
#======== NMDS ============================
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


#======== Phylomorphospace =====
tree <- read.nexus("pruned_tree-mimosa.nex")


###========================================================###
###=== Disparity analysis and phylomorphospace for Inga =====
## testing with subset for disparity metrics

filamentlenght = c()
calyxlength = c()
corollalength = c()
stamencount = c()
inflorescencetype = c()

Inga = cbind(filamentlength,calyxlength,corollalength,stamencount,inflorescencetype)
rownames(Inga) = c()


