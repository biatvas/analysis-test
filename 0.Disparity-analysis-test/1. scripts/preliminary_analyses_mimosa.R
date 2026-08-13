#==============================================================================#
# 1. Head ----------------------------------------------------------------------
#==============================================================================#

setwd("Documents/GitHub/analysis-test/0.Disparity-analysis-test/") #defining work directory

if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, purrr, factoextra) #installing and/or loading packages

morpho_data <- read.csv("1.datasets/mimosa_data_quali.csv", sep = ";") #reading morphological dataset
morpho_data_2 <- read.csv("1.datasets/mimosa_data_example.csv", sep= ",")
  
#==============================================================================#
# 2. Checking morpho dataset ---------------------------------------------------
#==============================================================================#

#==============================================================================#
# 2.1 Checking morphological dataset -------------------------------------------
#==============================================================================#

str(morpho_data) #checking the type of each variable

morpho_data %>% 
  select(infOrganization:corollaShape) %>%
  map(unique) #checking the states for each quali variable

#There are NA for all quali variables. But how much?
sum(is.na(morpho_data[,2:21])) / prod(dim(morpho_data[,2:21])) * 100 #calculating the percentage of NAs. 
#I am not selecting the first column bc it is the species name. Result: ~6.74% is NA. Not too bad.

# Now, I will adjust the other dataframe
str(morpho_data_2) #All variable are a chr. I will transform some of them

morpho_data_2 %>%
  mutate_at(vars(colWidth:styleLength),
            ~ gsub(",", ".", .x)) %>%
  mutate_at(vars(colWidth:styleLength),
                as.numeric) -> morpho_data_2

str(morpho_data_2)

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

# Next, we will do a PCOA and a NMDS

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

