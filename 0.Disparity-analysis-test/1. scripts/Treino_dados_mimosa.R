library(readxl)
mimosa_data_example <- read_excel("C:/Users/Eulália/Desktop/bia/dados para treino/mimosa_data_example.csv")
View(mimosa_data_example)
install.packages(stats)
library(stats)
dados = mimosa_data_example
matriz_distancia <- dist(dados, method = "euclidean")

print(matriz_distancia)
matriz_distancia_objeto <- as.matrix(matriz_distancia)
print(matriz_distancia_objeto)
library(ggplot2)
library(vegan)
pcoa <- cmdscale(matriz_distancia, k = 14, eig = TRUE)
pcoa_values <- pcoa$values
pcoa_points <- pcoa$points
plot(pcoa_points, type = "p", main = "PCoA Plot")
