# Title: Principal component analysis

# Name: Heather White

# Date created: 22/01/20

# Last modified: 24/04/20

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


# Reading in the landmarks
folder <- 'Data/Resampled_landmarks'


# Number of specimens I have
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


# Reading in the species names .csv - this is in the same order as the specimens
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
# Shows if there is grouping based on high complexity and morphology


# Read in complexity data
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
write.csv(PC_scores, "Y:/Heather/R_Projects/suture_methods_comparison/Analysis/PC_scores.csv")


######################################################################################################

# Other useful information from the PCA


# Accessing the minimum PC shape for PC1 - the extreme end, so I can plot this and the max end for on PCA figure
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


# Corrplot allows a comparison between the methods, shows if they are similarly related or not.
# This is better than doing regression analysis as regression analysis relies on which you choose as predictor or dependent

# Get the necessary data
# Five method results and PC1, PC2, PC3, PC4 (all >5% variation)
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


write.csv(cor2, "/Users/heatherwhite/Documents/PhD/PhD/R_Projects/suture_methods_comparison/Analysis/complexity_vs_PC_score_correlations.csv")




#############################################################################################

# Bivariate analysis - relationship between method and PC scores from shape data

# This was the analysis done for the original manuscript


# Plotting PC1 scores alongside SI result
ggplot(data, aes(x = PC1, y = SI))+ 
  geom_point() +
  xlab("PC1 score") +
  ylab("Sinuosity Index")+
  theme_bw()

# Fit the general linear model
model1 <- lm(PC1 ~ SI, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model1) 
mtext("Model assumptions (PC1 and SI)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model1)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model1)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC1, y = SI))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC1 score") +
  ylab("Sinuosity Index")+
  theme_bw(base_size = 14)



#################################################


# Plotting PC2 scores alongside SI result
ggplot(data, aes(x = PC2, y = SI))+ 
  geom_point() +
  xlab("PC2 score") +
  ylab("Sinuosity Index")+
  theme_bw()

# Fit the general linear model
model2 <- lm(PC2 ~ SI, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model2) 
mtext("Model assumptions (PC2 and SI)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model2)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model2)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC2, y = SI))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC2 score") +
  ylab("Sinuosity Index")+
  theme_bw(base_size = 14)


#################################################


# Plotting PC1 scores alongside SCI result
ggplot(data, aes(x = PC1, y = SCI))+ 
  geom_point() +
  xlab("PC1 score") +
  ylab("Suture Complexity Index")+
  theme_bw()

# Fit the general linear model
model3 <- lm(PC1 ~ SCI, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model3) 
mtext("Model assumptions (PC1 and SCI)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model3)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model3)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC1, y = SCI))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC1 score") +
  ylab("Suture Complexity Index")+
  theme_bw(base_size = 14)




#################################################


# Plotting PC2 scores alongside SCI result
ggplot(data, aes(x = PC2, y = SCI))+ 
  geom_point() +
  xlab("PC2 score") +
  ylab("Suture Complexity Index")+
  theme_bw()

# Fit the general linear model
model4 <- lm(PC2 ~ SCI, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model4) 
mtext("Model assumptions (PC2 and SCI)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model4)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model4)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC2, y = SCI))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC2 score") +
  ylab("Suture Complexity Index")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC1 scores alongside FD box result
ggplot(data, aes(x = PC1, y = FD_box))+ 
  geom_point() +
  xlab("PC1 score") +
  ylab("Fractal Dimension (box counting)")+
  theme_bw()

# Fit the general linear model
model5 <- lm(PC1 ~ FD_box, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model5) 
mtext("Model assumptions (PC1 and FD box)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model5)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model5)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC1, y = FD_box))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC1 score") +
  ylab("Fractal Dimension (box counting)")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC2 scores alongside FD box result
ggplot(data, aes(x = PC2, y = FD_box))+ 
  geom_point() +
  xlab("PC2 score") +
  ylab("Fractal Dimension (box counting)")+
  theme_bw()

# Fit the general linear model
model6 <- lm(PC2 ~ FD_box, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model6) 
mtext("Model assumptions (PC2 and FD_box)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model6)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model6)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC2, y = FD_box))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC2 score") +
  ylab("Fractal Dimension (box counting)")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC1 scores alongside FD madogram result
ggplot(data, aes(x = PC1, y = FD_madogram))+ 
  geom_point() +
  xlab("PC1 score") +
  ylab("Fractal Dimension (madogram)")+
  theme_bw()

# Fit the general linear model
model7 <- lm(PC1 ~ FD_madogram, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model7) 
mtext("Model assumptions (PC1 and FD madogram)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model7)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model7)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC1, y = FD_madogram))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC1 score") +
  ylab("Fractal Dimension (madogram)")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC2 scores alongside FD madogram result
ggplot(data, aes(x = PC2, y = FD_madogram))+ 
  geom_point() +
  xlab("PC2 score") +
  ylab("Fractal Dimension (madogram)")+
  theme_bw()

# Fit the general linear model
model8 <- lm(PC2 ~ FD_madogram, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model8) 
mtext("Model assumptions (PC2 and FD madogram)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model8)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model8)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC2, y = FD_madogram))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC2 score") +
  ylab("Fractal Dimension (madogram)")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC1 scores alongside PSD result
ggplot(data, aes(x = PC1, y = PSD))+ 
  geom_point() +
  xlab("PC1 score") +
  ylab("Power Spectrum Density")+
  theme_bw()

# Fit the general linear model
model9 <- lm(PC1 ~ PSD, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model9) 
mtext("Model assumptions (PC1 and PSD)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model9)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model9)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC1, y = PSD))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC1 score") +
  ylab("Power Spectrum Density")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC2 scores alongside PSD result
ggplot(data, aes(x = PC2, y = PSD))+ 
  geom_point() +
  xlab("PC2 score") +
  ylab("Power Spectrum Density")+
  theme_bw()

# Fit the general linear model
model10 <- lm(PC2 ~ PSD, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model10) 
mtext("Model assumptions (PC2 and PSD)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model10)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model10)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC2, y = PSD))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC2 score") +
  ylab("Power Spectrum Density")+
  theme_bw(base_size = 14)




#################################################


# Plotting PC3 scores alongside SI result
ggplot(data, aes(x = PC3, y = SI))+ 
  geom_point() +
  xlab("PC3 score") +
  ylab("Sinuosity Index")+
  theme_bw()

# Fit the general linear model
model11 <- lm(PC3 ~ SI, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model11) 
mtext("Model assumptions (PC3 and SI)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model11)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model11)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC3, y = SI))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC3 score") +
  ylab("Sinuosity Index")+
  theme_bw(base_size = 14)



#################################################

# Plotting PC3 scores alongside SCI result
ggplot(data, aes(x = PC3, y = SCI))+ 
  geom_point() +
  xlab("PC2 score") +
  ylab("Suture Complexity Index")+
  theme_bw()

# Fit the general linear model
model12 <- lm(PC3 ~ SCI, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model12) 
mtext("Model assumptions (PC3 and SCI)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model12)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model12)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC3, y = SCI))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC3 score") +
  ylab("Suture Complexity Index")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC3 scores alongside FD box result
ggplot(data, aes(x = PC3, y = FD_box))+ 
  geom_point() +
  xlab("PC3 score") +
  ylab("Fractal Dimension (box counting)")+
  theme_bw()

# Fit the general linear model
model13 <- lm(PC3 ~ FD_box, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model13) 
mtext("Model assumptions (PC3 and FD_box)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model13)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model13)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC3, y = FD_box))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC3 score") +
  ylab("Fractal Dimension (box counting)")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC3 scores alongside FD madogram result
ggplot(data, aes(x = PC3, y = FD_madogram))+ 
  geom_point() +
  xlab("PC3 score") +
  ylab("Fractal Dimension (madogram)")+
  theme_bw()

# Fit the general linear model
model14 <- lm(PC3 ~ FD_madogram, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model14) 
mtext("Model assumptions (PC3 and FD madogram)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model14)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model14)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC3, y = FD_madogram))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC3 score") +
  ylab("Fractal Dimension (madogram)")+
  theme_bw(base_size = 14)




#################################################

# Plotting PC3 scores alongside PSD result
ggplot(data, aes(x = PC3, y = PSD))+ 
  geom_point() +
  xlab("PC3 score") +
  ylab("Power Spectrum Density")+
  theme_bw()

# Fit the general linear model
model15 <- lm(PC3 ~ PSD, data = data)

# Evaulating the assumptions of the model
par(mfrow = c(2,2))
plot(model15) 
mtext("Model assumptions (PC3 and PSD)", outer=T, cex=1, line=-2)
par(mfrow = c(1,1))

# To get a summary of the model
summary(model15)
# To get the anova summary (Pr(>F) = p-value) 
# I don't think this anova bit is necessary
anova(model15)

# Add the model interpretation to the scatter plot (ie the line)
ggplot(data, aes(x = PC3, y = PSD))+ 
  geom_point(colour = "cornflowerblue") +
  geom_smooth(method = 'lm', colour = "black", size = 1) +
  xlab("PC3 score") +
  ylab("Power Spectrum Density")+
  theme_bw(base_size = 14)



