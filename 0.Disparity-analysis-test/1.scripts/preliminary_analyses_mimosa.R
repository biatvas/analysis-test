#==============================================================================#
# 1. Head ----------------------------------------------------------------------

#defining work directory setwd("")

if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, purrr, factoextra, readr, stats, tidyr, stringr,
                 dispRity, ggplot2, vegan, FactoMineR, ade4, cluster, ape, ticyverse) #installing and/or loading packages

morpho_data <- read.csv("Documents/GitHub/bvasconcelos-IC-disparidade-floral/1.datasets/mimoseae_subset_clean.csv") #reading morphological dataset
  
#==============================================================================#
# 2. Checking morpho dataset ---------------------------------------------------
#==============================================================================#
# 2.1 Select columns and correct data -----------------------------------------
str(morpho_data) #checking the type of each variable

traits <- cbind("taxon" = morpho_data$taxon, morpho_data[,6:83])
#244 obs and 79 variables

# Dados contínuos em traits nem sempre estão corretamente organizados nos seus respectivos 
# min, low, high e max. São as colunas contínuas:
cols <- names(traits)
range_traits <- unique(sub("(.+)_min(_|$).*$", "\\1", cols[grepl("_min(_|$)", cols)]))

min_col <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                               "_min(_|$)"),  cols, value = TRUE)))
low_col <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                               "_low(_|$)"),  cols, value = TRUE)))
high_col <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                "_high(_|$)"), cols, value = TRUE)))
max_col  <- unlist(lapply(range_traits,  function(x) grep(paste0("^", x, 
                                                                 "_max(_|$)"), cols, value = TRUE)))

all.equal(length(min_col), length(low_col), length(high_col), length(max_col))

continuous_col <- c(min_col, low_col, high_col, max_col)
continuous_col[!continuous_col %in% colnames(traits)] #checando

# Podem estar ter valores imputados como:
unique(unlist(unname(traits[continuous_col]))) #logo, antes de transformações precisam ser corrigidos

# Gerando um outro banco de dados para trabalharmos sem alterar o original
traits_2 <- traits

## 3.1 Removing spaces ####
#=========================#

# Primeiro, vamos remover espaços das colunas contínuas. 
no_space <- as.data.frame(sapply(traits_2[continuous_col], function(x) gsub("\\s+", "", x)))

#verificando se nada foi alterado no nome das colunas ou se NA foram gerados
all(names(no_space) == names(traits_2[continuous_col])) 
sum(is.na(traits_2[continuous_col])) #121684
sum(is.na(no_space)) #121684

#já que tudo parece estar certo, trocando as colunas contínuas para as sem espaços em traits_2
traits_2[continuous_col] <- no_space
sum(is.na(traits)) #166507
sum(is.na(traits_2)) #166507

#============================#
## 3.2 Removing ≤ symbols ####
#============================#

# Abaixo, faremos com que todas as ocorrencias de números iniciadas com ≤ , por exemplo "≤8", 
# sejam copiadas nas colunas de nome high e o simbolo sejá removido

#Novamente, criando um traits_3 para não modificar o que já estava certo antes
traits_3 <- traits_2

for (col in continuous_col) {
  
  root <- sub("_(min|low|high|max)$", "", col)
  
  col_high <- paste0(root, "_high") #isso aqui é porque os valores serão adicionados em high
  
  if (!(col_high %in% names(traits_3))) next
  
  vals <- traits_3[[col]]
  high_vals <- traits_3[[col_high]]
  
  idx <- grepl("^≤", vals) & (is.na(high_vals) | high_vals == "") #entao se tem esse simbolo em alguma coluna e os valores high são NA, então transfere pra high
  
  traits_3[[col_high]][idx] <- sub("^≤", "", vals[idx])
  traits_3[[col]][idx] <- NA_character_ #a coluna que tinha o símbolo e não era high vira NA
}

#verificando se deu certo:
##antes era:
traits_2$inflorescence_length_min[grepl("^≤", traits_2$inflorescence_length_min)]
which(grepl("^≤", traits_2$inflorescence_length_min) == T)
## agora:
traits_3$inflorescence_length_min[which(grepl("^≤", traits_2$inflorescence_length_min) == T)] #o que pedimos para virar NA
traits_3$inflorescence_length_high[which(grepl("^≤", traits_2$inflorescence_length_min) == T)] #os valores foram adionados em high

sum(is.na(traits_2)) #166507
sum(is.na(traits_3)) #166507

# Já que está tudo certo com traits_3, podemos substituir o dataset
traits_2 <- traits_3
remove(traits_3)

#============================#
## 3.3 Removing ± symbols ####
#============================#

#Se tem o símbolo ±, então o valor é adicionado em low e em high. Onde ocorria antes é removido. 

#criando um dataset para trabalharmos
traits_3 <- traits_2

for (col in continuous_col) {
  
  root <- sub("_(min|low|high|max)$", "", col)
  
  col_low  <- paste0(root, "_low") #porque incluiremos em low
  col_high <- paste0(root, "_high") #porque incluiremos em high
  
  if (!(col_low %in% names(traits_3)) || !(col_high %in% names(traits_3))) next #para colunas que não são low e high, aplica o definido abaixo
  
  vals <- traits_3[[col]]
  low_vals  <- traits_3[[col_low]]
  high_vals <- traits_3[[col_high]]
  
  #seleciona quem tem ± no início
  idx <- grepl("^±\\s*", vals) &
    (is.na(low_vals) | low_vals == "") &
    (is.na(high_vals) | high_vals == "")
  
  #tira o simbolo
  clean_vals <- sub("^±\\s*", "", vals[idx])
  
  #copia em low e high
  traits_3[[col_low]][idx]  <- clean_vals
  traits_3[[col_high]][idx] <- clean_vals
  
  #apaga onde estavam as ocorrencias de +/-
  traits_3[[col]][idx] <- NA_character_
}

#verificando se deu tudo certo usando inflorescence_length_min como exemplo 
traits_2$inflorescence_length_min[grepl("^±", traits_2$inflorescence_length_min)]
which(grepl("^±", traits_2$inflorescence_length_min) == T)

traits_3$inflorescence_length_min[which(grepl("^±", traits_2$inflorescence_length_min) == T)]
traits_3$inflorescence_length_high[which(grepl("^±", traits_2$inflorescence_length_min) == T)]

sum(is.na(traits_3)) #166506
sum(is.na(traits_2)) #166506

#td certo, então subtituindo os dataset
traits_2 <- traits_3
remove(traits_3)

#===============================#
## 3.4 Removing "()" symbols ####
#===============================#

# agora resolvendo problemas dos que estão entre parenteses, atribuindo esses valores para
# low, min, max ou high, de acordo com a posição dos valores

#dataset para trabalharmos
traits_3 <- traits_2

for (col in continuous_col) {
  
  root <- sub("_(min|low|high|max)$", "", col)
  
  col_min  <- paste0(root, "_min")
  col_low  <- paste0(root, "_low")
  col_high <- paste0(root, "_high")
  col_max  <- paste0(root, "_max")
  col_check <- paste0(root, "_check") #pra checar se deu certo
  
  vals <- traits_3[[col]]
  
  # selecionando quais colunas tem casos de parenteses pra deixar elas na coluna check, caso seja preciso checar
  idx <- grepl("\\(", vals)
  traits_3[[col_check]] <- NA_character_
  traits_3[[col_check]][idx] <- vals[idx]
  
  x <- vals[idx]
  
  #extraindo os numeros dos parenteses no inicio do numero medio, ou seja, valroes que serão incluidos no min
  
  has_min <- grepl("^\\(", x) 
  min_vals <- rep(NA_character_, length(x))
  
  min_vals[has_min] <- sub("^\\(([^\\)]+)\\).*", "\\1", x[has_min])
  min_vals <- sub("-$", "", min_vals)
  
  traits_3[[col_min]][idx] <- min_vals
  
  #extraindo os numeros dos parenteses no final do numero medio, ou seja, valroes que serão incluidos no max
  
  has_max <- grepl("\\)$", x)
  max_vals <- rep(NA_character_, length(x))
  
  max_vals[has_max] <- sub(".*\\(([^\\)]+)\\)$", "\\1", x[has_max])
  max_vals <- sub("^-", "", max_vals)
  
  traits_3[[col_max]][idx] <- max_vals
  
  #removendo parenteses dos que trocamos  
  core <- gsub("\\([^\\)]+\\)", "", x)
  core <- gsub("--+", "-", core)
  core <- gsub("^-|-$", "", core)
  
  #definindo quem é low e high
  low_vals  <- rep(NA_character_, length(core))
  high_vals <- rep(NA_character_, length(core))
  
  has_range <- grepl("-", core) #agora dos dos que tem parenteses, que tbm estao entre traços, ou seja, low e high 
  
  low_vals[has_range]  <- sub("^(\\d+\\.?\\d*).*", "\\1", core[has_range])
  high_vals[has_range] <- sub(".*-(\\d+\\.?\\d*)$", "\\1", core[has_range])
  #outras inconsistencias com hifen
  no_range <- !has_range & core != ""
  low_vals[no_range] <- core[no_range]
  
  traits_3[[col_low]][idx]  <- low_vals
  traits_3[[col_high]][idx] <- high_vals
  
}

#verificando se deu certo:
sum(grepl("\\(", traits_2$inflorescence_length_min)) #tinha 7 ocorrencias com ()
sapply(traits_2[continuous_col], function(x) sum(grepl("\\(", x))) #todas na coluna
#inflorescence_length_min

traits_2$inflorescence_length_min[grepl("\\)", traits_2$inflorescence_length_min)]
which(grepl("\\(", traits_2$inflorescence_length_min) == T)

#valores que antes estavam entre parenteses e antes do número médio, agora estão em min
traits_3$inflorescence_length_min[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]
#valores que antes estavam entre parenteses e depois do número médio, agora estão em max
traits_3$inflorescence_length_max[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]

# valores que não estavam entre parenteses e que se referiam ao low, agora estçai em low
traits_3$inflorescence_length_low[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]

# valores que não estavam entre parenteses e que se referiam ao high, agora estão em high
traits_3$inflorescence_length_high[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]

#de traits_3 selecionando apenas as colunas que tbm estão em traits_2, ou seja, removendo os check
traits_3 <- select(traits_3, -c((names(traits_3[!colnames(traits_3) %in% colnames(traits_2)]))))

#verificando se tá td certo
all(names(traits_3) %in% names(traits_2))
#é pra ter diminuido o numero de NA, já que dividimos uma coluna em outras
sum(is.na(traits_3)) #166489
sum(is.na(traits_2)) #166506

#tudo parece estar certo, substituindo os dataset
traits_2 <- traits_3
remove(traits_3)

#verificando o que falta arrumar
unique(unlist(unname(traits_2[continuous_col])))

#==============================#
## 3.5 Removing "-" symbols ####
#==============================#

#Tratando os numeros entre "-", atribuindo os valores para low ou para high

#dataset para trabalharmos
traits_3 <- traits_2

for (col in continuous_col) {
  
  root <- sub("_(min|low|high|max)$", "", col)
  
  col_low  <- paste0(root, "_low")
  col_high <- paste0(root, "_high")
  
  vals <- traits_3[[col]]
  
  idx <- grepl("^\\d+\\.?\\d*-\\d+\\.?\\d*$", vals)
  
  low_vals  <- sub("^(\\d+\\.?\\d*)-.*", "\\1", vals[idx])
  high_vals <- sub(".*-(\\d+\\.?\\d*)$", "\\1", vals[idx])
  
  traits_3[[col_low]][idx]  <- low_vals
  traits_3[[col_high]][idx] <- high_vals
  
  traits_3[[col]][idx] <- NA_character_
  
}

#verificando
traits_2$inflorescence_length_min[grepl("-", traits_2$inflorescence_length_min)]
which(grepl("-", traits_2$inflorescence_length_min) == T)

traits_3$inflorescence_length_low[which(grepl("-", traits_2$inflorescence_length_min) == T)]
traits_3$inflorescence_length_high[which(grepl("-", traits_2$inflorescence_length_min) == T)]

#novamente, é pra diminuir o número de NA, porque estão dividindo o que tinha em apenas uma coluna em duas
sum(is.na(traits_3)) #166440
sum(is.na(traits_2)) #166489

#porque aparentemente está tudo certo, substituindo os datasets
traits_2 <- traits_3
remove(traits_3)
all(names(traits) %in% names(traits_2))

#========================#
# 2. Merging  columns ####
#========================#

# função para incluir valores de min e max em low e high, respectivamente

update_trait_values <- function(traits, min_col, low_col, high_col, max_col) {
  
  min_to_low <- which(is.na(traits[[low_col]]) & !is.na(traits[[min_col]])) #se é NA em low e não é em min, min é transferido pra low 
  max_to_high <- which(is.na(traits[[high_col]]) & !is.na(traits[[max_col]])) #se é NA em max e não é em high, high é transferido pra low
  
  if (length(min_to_low)) {
    traits[[low_col]][min_to_low] <- traits[[min_col]][min_to_low] #transfere min para low
    traits[[min_col]][min_to_low] <- NA #exclui em min
  }
  if (length(max_to_high)) {
    traits[[high_col]][max_to_high] <- traits[[max_col]][max_to_high] #transfere max para min
    traits[[max_col]][max_to_high] <- NA #exclui em max
  }
  
  return(traits)
}

cols <- names(traits)
range_traits <- unique(sub("(.+)_min(_|$).*$", "\\1", cols[grepl("_min(_|$)", cols)]))

#gerando noto dataset para trabalharmos
traits_2 <- traits

for (root in range_traits) {
  min_col  <- grep(paste0("^", root, "_min(_|$)"),  cols, value = TRUE) 
  low_col  <- grep(paste0("^", root, "_low(_|$)"),  cols, value = TRUE)
  high_col <- grep(paste0("^", root, "_high(_|$)"), cols, value = TRUE)
  max_col  <- grep(paste0("^", root, "_max(_|$)"),  cols, value = TRUE)
  
  n <- min(length(min_col), length(low_col), length(high_col), length(max_col))
  
  if (n > 0) {
    for (i in seq_len(n)) { #looping para fazer para todas as ocorrencias
      traits_2 <- update_trait_values(
        traits_2,
        min_col[i],
        low_col[i],
        high_col[i],
        max_col[i] 
      )
    }
  }
}

sum(is.na(traits_2)) #167229
sum(is.na(traits)) #167229

#===========================================#
## 2.1 Mean values for continuous traits ####
#===========================================#

# Verificando se existe algo que tem informação em min mas não tem em low ou se tem em max mas não tem em high
min_cols <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                "_min(_|$)"),  cols, value = TRUE)))

low_cols <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                "_low(_|$)"),  cols, value = TRUE)))

high_cols <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                 "_high(_|$)"), cols, value = TRUE)))

max_cols  <- unlist(lapply(range_traits,  function(x) grep(paste0("^", x, 
                                                                  "_max(_|$)"), cols, value = TRUE)))

any(sapply(seq_along(min_cols), function(i) {
  idx <- !is.na(traits_2[[min_cols[i]]]) & traits_2[[min_cols[i]]] != "" &
    (is.na(traits_2[[low_cols[i]]]) | traits_2[[low_cols[i]]] == "")
  any(idx)
})) #verifica se há algum caso que tem NA ou é vazio na coluna low, mas tem dados na coluna min

any(sapply(seq_along(max_cols), function(i) {
  idx <- !is.na(traits_2[[max_cols[i]]]) & traits_2[[max_cols[i]]] != "" &
    (is.na(traits_2[[high_cols[i]]]) | traits_2[[high_cols[i]]] == "")
  any(idx)
})) #verifica se há algum caso que tem NA ou é vazio na coluna high, mas tem dados na coluna max

# já que não há nada que tenha em min e max que não tenha valores em low e high (ou seja, FALSE foi retornado), 
# vou remover as colunas min e max e tirar a média entre low e high (se só houver apenas um valor, 
# ele será usado)

traits_2 #3172 obs, 79 variables
sum(is.na(traits_2)) #167229

for (i in seq_along(low_cols)) {
  
  low  <- as.numeric(traits_2[[low_cols[i]]])
  high <- as.numeric(traits_2[[high_cols[i]]])
  
  mean_col <- sub("_low$", "_mean", low_cols[i])
  
  traits_2[[mean_col]] <- rowMeans(
    cbind(low, high),
    na.rm = TRUE
  )
}

#Verificando
sum(is.na(traits_2[colnames(traits)])) #167229, o mesmo que antes, então não foram gerados NAs ao estimar a média

all.equal(traits_2$height_mean, rowMeans(cbind(as.numeric(traits_2$height_low),
                                               as.numeric(traits_2$height_high)),na.rm = TRUE))

#===================================#
## 2.2 Keeping only mean columns ####
#===================================#

# removendo colunas com low, min, max, high (manter só mean)
continuous_col <- c(min_cols, low_cols, high_cols, max_cols)
traits_3 <- traits_2

traits_3 <- traits_3[!colnames(traits_3) %in% continuous_col]
#43 variaveis e 3172 obs
sum(is.na(traits_3)) #69858

#============================#
# 3. Unit standardization ####
#============================#
# traits_2_backup <- traits_2

traits_2 <- traits_3
remove(traits_3)
sum(is.na(traits_2)) #69858

#==========================#
## 3.1 Correcting typos ####
#==========================#
traits_3 <- traits_2

#corrigindo um erro na escrita
colnames(traits_3) <- sub("calyx_lobe_length_unit.", "calyx_lobe_length_unit", colnames(traits_3))

cols <- names(traits_3)
range_traits <- unique(sub("(.+)_mean(_|$).*$", "\\1", cols[grepl("_mean(_|$)", cols)]))

unit_col <- paste(range_traits, "unit", sep = "_") #colunas com unit
all(unit_col %in% colnames(traits_3)) #todas colunas de unit_col está em traits

unique(unname(unlist(lapply(traits_3[unit_col], function (x) unique(x))))) #tom? 

which(apply(traits_3[unit_col], 1, function(row) {
  any(row == "tom", na.rm = TRUE)
}) == T)

traits_3[165,]$taxon

#corrigindo
traits_3[165,"pedicel_width_unit"] <- "mm"
#sum(is.na(traits_3)) #69858

unique(unname(unlist(lapply(traits_3[unit_col], function (x) unique(x)))))

traits_2 <- traits_3
remove(traits_3)

#=======================================#
## 3.2 Applying unit standardization ####
#=======================================#
traits_3 <- traits_2

mean_cols <- names(traits_3)[grepl("_mean$", names(traits_3))]

for (var in mean_cols) {
  
  unit_var <- sub("_mean$", "_unit", var)
  
  if (!unit_var %in% names(traits_3)) next
  
  traits_3[[unit_var]] <- as.character(traits_3[[unit_var]])
  traits_3[[unit_var]] <- trimws(tolower(traits_3[[unit_var]]))
  
  idx_m  <- traits_3[[unit_var]] %in% "m"
  idx_dm <- traits_3[[unit_var]] %in% "dm"
  idx_mm <- traits_3[[unit_var]] %in% "mm"
  
  traits_3[[var]][idx_m]  <- traits_3[[var]][idx_m] * 100
  traits_3[[var]][idx_dm] <- traits_3[[var]][idx_dm] * 10
  traits_3[[var]][idx_mm] <- traits_3[[var]][idx_mm] / 10
  
  idx_valid <- !is.na(traits_3[[unit_var]]) & traits_3[[unit_var]] != ""
  traits_3[[unit_var]][idx_valid] <- "cm"
}

#tem que ter a mesma soma de NA
sum(is.na(traits_2)) #69858
sum(is.na(traits_3)) #69858

which(traits_2$inflorescence_length_unit == "mm")[3]
traits_2[10,"inflorescence_length_unit"]
traits_2[10,"inflorescence_length_mean"] #1.75
#precisa ser 1.75/10
traits_3[10,"inflorescence_length_mean"]

which(traits_2$inflorescence_length_unit == "cm")[7]
traits_2[95,"inflorescence_length_unit"]
traits_2[95,"inflorescence_length_mean"]
#precisa ser 3
traits_3[95,"inflorescence_length_mean"]

#checando se algum NA foi introduzido
which(is.na(traits_3[mean_cols]) & !is.na(traits_2[mean_cols]) == T)

#pelos testes, parece estar tudo ok. podemos remover as colunas descrevendo as unidades
traits_2 <- traits_3
remove(traits_3)

#traits contain 31 variables
traits_3 <- traits_3[, !names(traits_3) %in% unit_col]
sum(is.na(traits_3)) #54868

#=========================================#
# 4. Saving continuous cleaned dataset ####
#=========================================#

write.csv(traits_3, "3.outputs/continuous_data_cleaned-20260410.csv", row.names = F,
          fileEncoding = "UTF-8")

#============================#
# 5. Dataset completeness ####
#============================#

cleaned_traits <- read.csv("3.outputs/continuous_data_cleaned-20260410.csv")

#checking traits with less than 15% completeness
traits_percent <- colMeans(!is.na(cleaned_traits)) * 100
names(traits_percent[traits_percent < 15])
traits_percent_original <- colMeans(!is.na(traits)) * 100
names(traits_percent_original[traits_percent_original < 15])
trait_completeness <- colMeans(!is.na(cleaned_traits[,-1])) * 100

#checking less than 30% 
names(traits_percent[traits_percent < 30])

# acho que podemos revisar os dados não contínuos tbm e depois a gnt faz a filtragem

# traits_filtered <- cleaned_traits[, c(TRUE, trait_completeness >= 30)]

#now i will keep 30% filtered columns and remove taxa with less than 60% completeness
# species_completeness <- rowMeans(!is.na(traits_filtered[,-1])) * 100
# traits_filtered$species_completeness <- species_completeness
#write.csv(traits_filtered, "4.outputs/traits_filtered.csv", row.names = F)

#remove species less than 60% completeness
# traits_final <- traits_filtered %>%
#   dplyr::filter(species_completeness >= 60) %>%
#   dplyr::select(-species_completeness)

#now we have 1292 species

#write.csv(traits_final,"final_species_filtered_data.csv", row.names = F)

# traits_final <- traits_final %>% 
#   separate(taxon, into = c("genus", "epithet"),
#            sep = "_")
# unique(traits_final$genus)
# traits_final %>% distinct() %>% count(genus, name = "traits_final")

#==============================================================================#
# 2.4 Check quali data ---------------------------------------------------------



#==============================================================================#
# 3. Morphospaces --------------------------------------------------------------
#==============================================================================#

#==============================================================================#
# 3.1 PCA ----------------------------------------------------------------------
#==============================================================================#
# Morphospaces are PCA, PCoA, NMDS... I will conduct a PCA first, so only with quantitative variables
row.names(morpho_data_2) <- morpho_data_2$X
morpho_data_2 <- select(morpho_data_2, -X)

pca <- prcomp(morpho_data_2[,1:14], scale. = F, center = T)

#doing a PCA by myself
X <- morpho_data_2[,1:14] #the quanti data
X_centered <- scale(X, center = T, scale = F) #centering the data but not scaling them bc all variables have the same nature (measurements in mm)
C <- cov(X_centered) #calculating the covariance matrix
eig <- eigen(C) #computing eigenvalues and eigenvectors of the covariance matrix
eig$vectors

PCA_scores <- X_centered %*% eig$vectors #calculating the PCA scores for each species
colnames(PCA_scores) <- paste0("PC", 1:14)

plot(PCA_scores[,1:2]) #the same as plot(pca$x)

# the variables are not log transformed. This can affect our analyses.
fviz_contrib(pca, choice = "ind", axes = 1:2)
fviz_contrib(pca, choice = "var", axes = 1:2)

# filamentLength and styleLength are much more variable compared with other traits, 
# bc they are bigger. We need to normalize the variables. We can do that using a log transformation.
morpho_data_2 %>% 
  mutate_at(vars(colWidth:styleLength), log) -> morpho_data_2

# PCA for log-transformed data
pca_2 <- prcomp(morpho_data_2[,1:14], scale. = F, center = T)

# Comparing PCA without log-transformed data and log-transformed data
plot(pca_2$x)
plot(pca$x)

#ind variation (i.e., species scores)
fviz_contrib(pca_2, choice = "ind", axes = 1:2)
fviz_contrib(pca, choice = "ind", axes = 1:2)

#variable variation (i.e., traits scores)
fviz_contrib(pca_2, choice = "var", axes = 1:2)
fviz_contrib(pca, choice = "var", axes = 1:2)

#Log-transformed is much better.
pca <- pca_2

plot(pca$x)

#We can remove everything we dont need now and clean our global environment.
remove(C,eig,morpho_data,PCA,X,X_centered,PCA_scores,pca_2)

# Next, we will do a PCoA and a NMDS

#==============================================================================#
# 4. Phylomorphospaces ---------------------------------------------------------
#==============================================================================#

#==============================================================================#
# 5. Phylo signal --------------------------------------------------------------
#==============================================================================#

#==============================================================================#
# 6. Disparity analyses --------------------------------------------------------
#==============================================================================#

#==============================================================================#
# 6. Checking differences between groups ---------------------------------------
#==============================================================================#

