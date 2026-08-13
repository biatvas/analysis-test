library(readxl)
mimosa_data_example <- read.csv("C:/Users/Eulália/Desktop/bia/Labis/mimosa_data_example.csv", sep = ";")
View(mimosa_data_example)
str(mimosa_data_example)
# padronização dos dados
dados_numericos <- mimosa_data_example[, -1]  # remove a coluna 1 (nomes)
dados_padronizados <- scale(dados_numericos)
str(dados_padronizados)
#matriz de correlação dos dados
matriz_cor <- cor(dados_padronizados)
print(matriz_cor)
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
scores$especie <- mimosa_data_example[, 1] 
library(ggplot2)
ggplot(scores, aes(x = PC1, y = PC2, label = especie)) +
  geom_point(color = "red") +
  geom_text(vjust = -0.5, size = 2) +
  theme_minimal() +
  labs(title = "PCA - Espécies de Mimosa",
       x = "PC 1",
       y = "PC 2")
#para saber qual variável tem maior peso na contribuição dos 
##componentes principais
pca_resultado$rotation
install.packages("phytools")
library(phytools)
install.packages("dispRity")
install.packages("FactoMineR")
library(FactoMineR)
library(dispRity)
str(mimosa_data_example)
help("PCA")
#PCA(X, scale.unit = TRUE, ncp = 5, ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL, row.w = NULL, col.w = NULL, graph = TRUE, axes = c(1,2))
help("dispRity")
