# Title: Calculating fractal dimension

# Name: Heather White

# Date created: 29/05/19

# Last modified: 13/06/19

# License: MIT license



# Script to calculate the fractal dimension of each specimen

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(fractaldim)
library(data.table)
library(pcaPP)
library(wavelets)
library(geomorph)



#######################################################################################################


# Reading in the .csv files into an array:

folder <- 'Data/Resampled_landmarks'
# Specify the number of taxa in the folder:
ntaxa<-79 
# Read in the file names to a list variable:
csvlist <- list.files(path = folder, pattern = "*.csv", full.names = T)
# Set the size of the array:
csvarray<-array(dim=c(500,2,ntaxa))
#dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
# Add the landmarks for each specimen to the array:
for(i in 1:length(csvlist))
{
  csvarray[,,i]<-as.matrix(read.csv(file=csvlist[i],header=T,sep=",",row.names=1))
}
# To call on one of the specimens from the array use CSVarray[,,1] for example in the console to access the first specimen



# Associating a specimen name with each array
dimnames(csvarray)[3]<-list(
  substr(dir("./",pattern=".csv"),1,(nchar(dir("./",pattern=".csv"))-4)))###gives specimens names based on file names in ply folder
lm_array<-csvarray



# Checking the array:

# Check the dimensions of lm_array
dim(lm_array)

# Check the class of lm_array - for fractal dimension this has to be matrix or dataframe
class(lm_array)
# Check random specimens in the lm_array
lm_array[,,1]
lm_array[,,6]
# Check the class of each specimen within the lm_array
class(lm_array[,,1])



# Reading in the species names .csv - this is in the same order as the specimens
species_names <- read.csv('Data/Species_names_ordered.csv', header = F)




#######################################################################################################

# Fractal dimension method test


# 2D Methods: 'transect.var', 'transect.incr1', 'isotropic', 'squareincr' and 'filter1' can be used for 2D matrices or dataframes
# 2D matrices have to be square matrices

# 1D Methods: 'boxcount', 'hallwood', 'variogram', 'madogram', 'rodogram', 'variation', 'incr1', 'genton', 'periodogram', 'wavelet', 'dctII'
# ????????????????????????????
# ??? These methods seem to use one column of data so I am not sure how it is calculating fractal dimension for the 2 coordinate columns???


# Test to see how FD uses the matrices (1 or both sets of coordinates):

# First specimen matrix
a<-lm_array[,,1]
# x-coordinates of first specimen
b<-lm_array[,,1][,1]
# y-coordinates of first specimen
c<-lm_array[,,1][,2]

# Setting up the figure structure
par(mfrow=c(1,3))
# FD of first specimen, D = 1.46
a1<-fd.estim.madogram(a, nlags = "all", plot.loglog = TRUE,
                   plot.allpoints = TRUE)
# Gives summary of a1 - mean FD but no standard deviation
summary(a1)
# FD of x-coordinates of first specimen, D = 1.04
b1<-fd.estim.madogram(b, nlags = "all", plot.loglog = TRUE,
                      plot.allpoints = TRUE)
# Gives summary of b1 - mean FD but no standard deviation
summary(b1)
# FD of y-coordinates of first specimen, D = 1.54
c1<-fd.estim.madogram(c, nlags = "all", plot.loglog = TRUE,
                      plot.allpoints = TRUE)
# Gives summary of c1 - mean FD but no standard deviation
summary(c1)


# All three give different values of D, I assume this means that using the method for the specimen as a whole uses both coordinate sets



##############################################################################################################

# Selecting which FD method to use

# Setting up the figure structure
par(mfrow=c(2,4))
# Testing some of the different methods
fd.estim.variogram (a, nlags = 20, plot.loglog = TRUE)
fd.estim.variation (a, nlags = 20, plot.loglog = TRUE)
fd.estim.variogram (a, nlags = 3, plot.loglog = TRUE,
                    plot.allpoints = TRUE)
fd.estim.variation (a, plot.loglog = TRUE, plot.allpoints = TRUE)
fd.estim.hallwood (a, nlags = 10, plot.loglog = TRUE)
fd.estim.boxcount (a, nlags = "all", plot.loglog = TRUE,
                   plot.allpoints = TRUE)
fd.estim.periodogram (a, plot.loglog = TRUE)
fd.estim.dctII (a, plot.loglog = TRUE)


# All of the different methods - not tested
fd.estim.boxcount (a, plot.loglog = FALSE, nlags = "auto",
                   shift.up=TRUE, plot.allpoints = FALSE, legend.type = 's', debuglevel = 0)
fd.estim.hallwood (a, plot.loglog = FALSE, nlags = "auto",
                   plot.allpoints = FALSE, legend.type = 's', debuglevel = 0)
fd.estim.variogram (a)
fd.estim.madogram (a)
fd.estim.rodogram (a)
fd.estim.variation (a, p.index = 1)
fd.estim.incr1(a, p.index=2)
fd.estim.genton (a)
fd.estim.periodogram (a, plot.loglog = FALSE, nlags = "auto")
fd.estim.wavelet (a, plot.loglog=FALSE, plot.allpoints = FALSE,
                  filter = "haar", J1 = max(1,floor(log2(length(data))/3-1)),
                  J0 = floor(log2(length(data))), legend.type = 's', debuglevel = 0)
fd.estim.dctII (a, plot.loglog = FALSE, nlags = "auto")


##########################################################################################################

# Scaling using Procrustes'


# Performing Procrustes'
Procrustes<-gpagen(lm_array)
# Saving coordinate data to a new variable once normalised
lm_array_corrected<-Procrustes$coords


# Scaling using Procrustes' seems to have worked, but this is why the final FD values are in a very small range of each other
# Because, it is scaling them to a similar thing.
# Julian thinks this method should be okay for scaling and length should not be a problem.
# But he will look at if there is a way possible to do scaling to the length - although, there will be biases for the suture length.


# Test to see how the distances compare when they have been corrected and not corrected for Procrustes'
# Below is before Procrustes' plotted against the corrected FD with Procrustes, large distance differences
dist_lm = sapply(1:79, function(i) sqrt(sum((lm_array[1,,i]-lm_array[500,,i])^2)))
plot(FD_array, dist_lm) # Might want to to change the FD_array variable here to before it has been corrected to get the plot
# Below is after Procrustes' correction, showing very small distance differences
# Hence why Julian thinks this should be enough.
dist_lm2 = sapply(1:79, function(i) sqrt(sum((lm_array_corrected[1,,i]-lm_array_corrected[500,,i])^2)))
plot(FD_array, dist_lm2)


##########################################################################################################

# Calculating fractal dimension for all specimens - madogram method

par(mfrow=c(2,4))
FD_array <- array(dim=c(79,1))
for(i in 1:length(csvlist))
{
  run<-fd.estim.madogram(lm_array_corrected[,,i], nlags = "all", plot.loglog = TRUE,
                                   plot.allpoints = TRUE, main = i)
  FD_array[i]<-run$fd
}



# To get a summary of the FD_array variable
FD_array
# To save the FD_array as an R object
save(FD_array, file="Analysis/FD_madogram_array.R")

# In FD_array get rownames to save as species names
# Have to use the function dimnames not rownames because it is an array
# The [1] accesses the first position of the array, when you do dim(FD_array) this gives you 79, 1. So using [1] renames the species names
dimnames(FD_array)[1]=species_names


# Writing FD summary results for each specimen to a csv
write.csv(FD_array, "Y:/Heather/R_Projects/2D_landmarks_sutures/Analysis/FD_results_madogram.csv")



##########################################################################################################

# Calculating fractal dimension for all specimens - boxcounting method

par(mfrow=c(2,4))
FD_array_box <- array(dim=c(79,1))
for(i in 1:length(csvlist))
{
  run_box<-fd.estim.boxcount(lm_array_corrected[,,i], nlags = "all", plot.loglog = TRUE,
                         plot.allpoints = TRUE, main = i)
  FD_array_box[i]<-run_box$fd
}



# To get a summary of the FD_array variable
FD_array_box
# To save the FD_array as an R object
save(FD_array_box, file="Analysis/FD_box_array.R")


# In FD_array get rownames to save as species names
# Have to use the function dimnames not rownames because it is an array
# The [1] accesses the first position of the array, when you do dim(FD_array) this gives you 79, 1. So using [1] renames the species names
dimnames(FD_array_box)[1]=species_names


# Writing FD summary results for each specimen to a csv
write.csv(FD_array_box, "Y:/Heather/R_Projects/2D_landmarks_sutures/Analysis/FD_results_box.csv")







