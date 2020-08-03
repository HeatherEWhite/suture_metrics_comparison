# Title: Principal component analysis

# Name: Heather White

# Date created: 22/01/20

# Last modified: 24/05/20

# License: MIT license



# PCA for 2D landmarks of test dataset

#######################################################################################################


rm(list = ls())

library(tidyverse)
library(geomorph)
library(Morpho)
library(rgl)
library(ape)
library(paleomorph)
library(RRPP)
library(arrayhelpers)
library(corrplot)
library(svd)



#######################################################################################################

# Reading in the landmarks and adding to an array for later analysis


# Read in the landmarks
folder <- 'Data/Resampled_landmarks'


# Number of specimens
nspecimens<-79
# Read in the file names to a list variable:
csvlist <- list.files(path = folder, pattern = "*.csv", full.names = T)
# Set the size of the array:
csvarray<-array(dim=c(500,2,nspecimens))
#dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
# Add the landmarks for each specimen to the array:
for(i in 1:length(csvlist))
{
  csvarray[,,i]<-as.matrix(read.csv(file=csvlist[i],header=T,sep=",",row.names=1))
}
# To call on one of the specimens from the array use CSVarray[,,1] for example in the console to access the first specimen



# Associating a specimen name with each array (setting the dimnames of each specimen)
dimnames(csvarray)[3]<-list(
  substr(dir("./",pattern=".csv"),1,(nchar(dir("./",pattern=".csv"))-4)))###gives specimens names based on file names in ply folder
lm_array<-csvarray
lm_array
dim(lm_array)


# Reading in the species names .csv - this needs to be in the same order as the specimens
species_names <- read.csv('Data/Species_names_ordered.csv', header = F)
dimnames(lm_array)[3]=species_names


####################################################################################################

# Analysis of the landmarks - PCA


# Procrustes analysis - remove non-shape aspects
Y.gpa=gpagen(lm_array)
test=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Raw size of specimens 

plot(test)



# Plot tangent space (PCA) for PC1 and PC2
# This PCA is of shape data for the semilandmarks with specimen labels
PCA=plotTangentSpace(test, axis1=3, axis2=4, warpgrids = F, label = dimnames(test)[[3]])
# This PCA is of shape data for the semilandmarks without specimen labels
PCA=plotTangentSpace(test, axis1=1, axis2=2, warpgrids = F)


# Plotting the PCA with a heatmap for complexity scores


# Read in the complexity data
complexity <- read.csv('Analysis/all_complexity_results.csv')
# Extract the data for each method
FD_box <- as.numeric(complexity[,2]) 
FD_mad <- as.numeric(complexity[,3]) 
PSD <- as.numeric(complexity[,4]) 
SI <- as.numeric(complexity[,5])
SCI <- as.numeric(complexity[,6]) 


# PCA heatmap 
# replace on first line the variable with each of the 5 method variables above
# Create 10 bins of the complexity data
bins=cut(FD_mad, breaks=10, ordered_result = TRUE) 
# Create a colour scale for the heatmap
colfac=colorRampPalette(c("blue", "red")) 
# Divide the colour scale into 10
col.gp=colfac(10)
# Associate the data with colour bins
names(col.gp) <- levels(bins)
col.gp <- col.gp[match(bins, names(col.gp))] # col.gp must NOT be a factor
# Plot shape data with groups coloured with these colour bins
# Red is the most complex suture and blue is the least complex
PCA=plotTangentSpace(test, groups= col.gp, axis1=2, axis2=3, warpgrids = T) 
# Plot PCA with specimen labels to check that the complexity scores have been associated correctly
PCA=plotTangentSpace(test, groups= col.gp, axis1=1, axis2=2, warpgrids = F, label = dimnames(test)[[3]]) 




##############################################################################################

# Accessing PC scores from the PCA


# Accessing the PC scores for all the specimens
PC_scores <- PCA$pc.scores

# Writing PC score results for each specimen to a csv
write.csv(PC_scores, "Y:/xxx/xxx/xxx.csv")


######################################################################################################

# Other useful information from the PCA


# Accessing the minimum PC shape for PC1
PCA$pc.shapes$PC1min

# Summary of contributions from each PC axis
PCA$pc.summary


######################################################################################################

# PCA using svd package

# Obtain Procrustes coordinates
proc <- gpagen(lm_array)
# Convert coordinates to two-dimensional matrix
coords2d <- two.d.array(proc$coords)
# Calculate consensus and flatten to single vectors
consensus <- apply(proc$coords, c(1,2), mean)
consensusvec <- apply(coords2d, 2, mean)
# Calculate Procrustes residuals (Procrustes coordinates - consensus)
resids <- t(t(coords2d)-consensusvec)
# Calculate covariance matrix
P <- cov(resids)
# Calculate eigenvector and eigenvalues with SVD
pca.stuff <- svd(P)
eigenvalues <- pca.stuff$d
eigenvectors <- pca.stuff$u
# Calculate PCA scores
scores <- resids%*%eigenvectors


######################################################################################################

# Relationship between complexity score and PC score from shape data

# Using Corrplot


# Read in all the complexity results and PC score results
data <- read.csv('Analysis/all_complexity_results_with_PC_scores.csv', header = T)
data


# Get the necessary data
# Five method results and PC1, PC2, PC3, PC4 (PCs with >5% variation)
data2 <- dplyr::select(data, 2:10)

# Input data 'data' = matrix with specimens listed and every method result for each specimen
# cor function calculates the correlation between x and y.
cor <- cor(data2)
cor2 <- cor[1:5, 6:9]
# Significance testing
res1 <- cor.mtest(data2, method = 'pearson', conf.level = .95)
# To produce the corrplot - a graphical display of a correlation matrix
corrplot(cor2, method = "circle")
corrplot(cor2, method = "color")
# Producing a corrplot with crosses on comparisons not significant at the level defined
corrplot(cor, p.mat = res1$p, sig.level = .05)


write.csv(cor2, "Y:/xxx/xxx/xxx.csv")







