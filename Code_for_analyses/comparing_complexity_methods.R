# Title: Comparing complexity methods

# Name: Heather White

# Date created: 17/07/19

# Last modified: 21/04/20

# License: MIT license

# PCA analysis (individuals, variables), method trajectories, clustering, corrplots, PC loadings

##################################################################################################

rm(list = ls())

library(mvMORPH)
library(tidyverse)
library(ade4)
library(factoextra)
library(corrplot)
library(corrgram)
library(MASS)



# Read in the csv containing the results of all the methods
method_results <- read.csv('Analysis/all_complexity_results.csv')
# Extract the necessary data into a variable
data <- dplyr::select(method_results, 2:6)



# The number of specimens
  # Julien used a phylogney here to illustrate a hierarchical structure
clust <- pbtree(n=79)
# Plot the phylogeny
plotTree(clust)
# random dependencies between n number of metrics (in this 2) - the number of different methods I am comparing
cov_var <- crossprod(matrix(rnorm(3*3), 3))



# Get PCA data using ade4
res.pca <- dudi.pca(data,
                    scannf = FALSE,   # hide the scree plot
                    nf = 3            # Number of axes to keep
)




####################################################################################################

# Plotting PCAs to compare the different methods - using factoextra package

# This method used for a figure
# Plot a PCA with just the variables - ie the different methods compared
# This plot also highlights the different method arrows in different colours to show how much each one contributes to the variability
fviz_pca_var(res.pca,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             axes=c(2,3) # change this to 1 and 2 to get PC1 and PC2 or 2 and 3 to get PC2 and PC3
)


# This code below plots the method trajectories with the specimens with the eigenvalues for each PC.
scatter(res.pca, xax = 1, yax = 2)


# To see the eigenvalue contribution of each PC:
biplot(res.pca)


# Accessing eigenvalues for the method trajectory comparison PCA
summary(res.pca)




# Plot a PCA with both the specimens (n=79) and the variables (methods used)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969"  
)


###################################################################################################


# Comparing methods by plotting PCAs - using ade4 package 
# Similar to factoextra package


# PCA data calculated higher up, using ade4 - dudi.pca function
res.pca


# PCA circle - method contributions
s.corcircle(res.pca$co)


# A scatter plot with all the specimens and the method used represented as arrows
scatter(res.pca,
        posieig = "none", # Cacher le scree plot
        clab.row = 0      # Caché l'annotation de slignes
)



#####################################################################################################

# Clustering plotting to compare the methods

# Now we can also use a clustering approach (eg hierarchical cluster) or kmean to define distincts groupings (e.g. the 2 methods tested)
clusters <- kmeans(data, 3)


# Using ade4 package
# Clustering plot without labelling specimens and not on a PCA plot
groups <- as.factor(clusters$cluster)
s.class(res.pca$li,
        fac = groups,  # colorer par groupes
        col = 1:3#c("#00AFBB",  "#FC4E07")
)

# This method used for a figure
# Using factoextra package
# Clustering plot with labelled specimens on a PCA with the method trajectories
fviz_pca_biplot(res.pca,
                col.ind = groups, # colorer par groupes
                palette = 1:3,#c("#00AFBB",  "#FC4E07"),
                addEllipses = T, # Ellipses de concentrations
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = T,
)


# Clustering plot with labelled specimens on a PCA without the method trajectories
fviz_pca_ind(res.pca,
              col.ind = groups, # colorer par groupes
              palette = 1:3,#c("#00AFBB",  "#FC4E07"),
              addEllipses = TRUE, # Ellipses de concentrations
              ellipse.type = "confidence",
              legend.title = "Groups",
              repel = TRUE
)



###########################################################################################


# Corrplot


# Corrplot allows a comparison between the methods, shows if they are similarly related or not.
# This is better than doing regression analysis as regression analysis relies on which you choose as predictor or dependent


# Input data 'data' = matrix with specimens listed and every method result for each specimen
# cor function calculates the correlation between x and y.
cor <- cor(data)
# Calculating the p-value and CI for the methods using Pearson correlation - can also use Spearman's
res1 <- cor.mtest(data, method = 'pearson', conf.level = .95)
res2 <- cor.mtest(data, method = 'pearson', conf.level = .95)
# To produce the corrplot - a graphical display of a correlation matrix
corrplot(cor, method = "circle")
corrplot(cor, method = "color")
corrplot(cor, method = "number")
# Producing a corrplot with crosses on comparisons not significant at the level defined
corrplot(cor, p.mat = res1$p, sig.level = .001) # Everything significant p>0.01
# Producing a corrplot with the starred p-values
corrplot(cor, p.mat = res1$p, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")
corrplot.mixed(cor, lower = "number", upper = "circle")


# Can also calculate the covariances of x and y
cov <- cov(data)
# Can't plot a corrplot for this as the values are not within the range of -1, 1


# Can also calculate the variance of x
var <- var(data)
# Can't plot a corrplot for this as the values are not within the range of -1, 1


write.csv(cor, "/Users/heatherwhite/Documents/PhD/PhD/R_Projects/suture_methods_comparison/Analysis/Method_correlations.csv")
write.csv(cov, "/Users/heatherwhite/Documents/PhD/PhD/R_Projects/suture_methods_comparison/Analysis/Method_covariances.csv")
# Saving method variances in R - not a matrix so can't save a csv
save(var, file='Analysis/Method_variances.R')



# Example of a corrplot for the dimensions of a PCA.
var <- get_pca_var(res.pca)
corrplot(corrgram(var$cos2), is.corr = FALSE) # Using the corrgram function converts the 'var' variable to a matrix needed for corrplot




###########################################################################################

# PC loadings

# PC loadings shows the relationship between the variable (ie suture method) and PC axis - instead of doing bivariate analysis between PC axis and method


# Extracting PC loadings for first 3 PCs
pc_loadings <- res.pca$c1 # c1 = loadings of variables
write.csv(pc_loadings, "/Users/heatherwhite/Documents/PhD/PhD/R_Projects/suture_methods_comparison/Analysis/PC_loadings.csv")
pc_loadings2 <- read.csv('Analysis/PC_loadings.csv') # To get the method headings in the dataframe

# The variables' coordinates - normed to the square roots of the eigenvalues
# c1 is calculated from these co values
res.pca$co


# Bar plot of method loadings on PC1
ggplot(pc_loadings2, aes(x = X, y = CS1))+ 
  geom_bar(stat = 'identity', fill='cadetblue3', colour ='black', size =.4)+
  ggtitle('PC1 - 65.4% of total variance')+
  xlab("Complexity Methods")+
  ylab("Loadings")+
  geom_hline(yintercept = 0)+
  theme_test()

# Bar plot of method loadings on PC1 - reordered ascending
ggplot(pc_loadings2, aes(x = reorder(X, CS1, sum), y = CS1))+ 
  geom_bar(stat = 'identity', fill='cadetblue3', colour ='black', size =.4)+
  ggtitle('PC1 - 65.4% of total variance')+
  xlab("Complexity Methods")+
  ylab("Loadings")+
  geom_hline(yintercept = 0)+
  theme_test()


# Bar plot of method loadings on PC2 
ggplot(pc_loadings2, aes(x = X, y = CS2))+ 
  geom_bar(stat = 'identity', fill='cadetblue3', colour ='black', size =.4)+
  ggtitle('PC2 - 17.8% of total variance')+
  xlab("Complexity Methods")+
  ylab("Loadings")+
  geom_hline(yintercept = 0)+
  theme_test()

# Bar plot of method loadings on PC2 - reordered ascending
ggplot(pc_loadings2, aes(x = reorder(X, CS2, sum), y = CS2))+ 
  geom_bar(stat = 'identity', fill='cadetblue3', colour ='black', size =.4)+
  ggtitle('PC2 - 17.8% of total variance')+
  xlab("Complexity Methods")+
  ylab("Loadings")+
  geom_hline(yintercept = 0)+
  theme_test()


# Bar plot of method loadings on PC3
ggplot(pc_loadings2, aes(x = X, y = CS3))+ 
  geom_bar(stat = 'identity', fill='cadetblue3', colour ='black', size =.4)+
  ggtitle('PC3 - 13.7% of total variance')+
  xlab("Complexity Methods")+
  ylab("Loadings")+
  geom_hline(yintercept = 0)+
  theme_test()

# Bar plot of method loadings on PC3 - reordered ascending
ggplot(pc_loadings2, aes(x = reorder(X, CS3, sum), y = CS3))+ 
  geom_bar(stat = 'identity', fill='cadetblue3', colour ='black', size =.4)+
  ggtitle('PC3 - 13.7% of total variance')+
  xlab("Complexity Methods")+
  ylab("Loadings")+
  geom_hline(yintercept = 0)+
  theme_test()
  



# Method contributions to each PC axis

# Method contributions to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 5) # var is selecting the variables from the PCA

# Method contributions to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 5)

# Method contributions to PC3
fviz_contrib(res.pca, choice = "var", axes = 3, top = 5)



