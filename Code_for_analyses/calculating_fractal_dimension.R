# Title: Calculating fractal dimension

# Name: Heather White

# Date created: 29/05/19

# Last modified: 13/02/20

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


# Read in the LM files into an array:

folder <- 'Data/Resampled_landmarks'
# Specify the number of taxa:
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



# Associate a specimen name with each array
dimnames(csvarray)[3]<-list(
  substr(dir("./",pattern=".csv"),1,(nchar(dir("./",pattern=".csv"))-4)))###gives specimens names based on file names in ply folder
lm_array<-csvarray



# Checking the array:

# Check the dimensions of lm_array
dim(lm_array)

# Check the class of lm_array - for fractal dimension this has to be matrix or dataframe
class(lm_array)




# Read in the species names .csv - makre sure this is in the same order as the specimens
species_names <- read.csv('Data/Species_names_ordered.csv', header = F)




#######################################################################################################

# Fractal dimension method test


# Scaling using Procrustes'
# Performing Procrustes'
Procrustes<-gpagen(lm_array)
# Saving coordinate data to a new variable once normalised
lm_array_corrected<-Procrustes$coords



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
# Save the FD_array as an R object
save(FD_array, file="Y:/xxx/xxx/xxx.R")

# Associate species names with each FD madogram complexity score
dimnames(FD_array)[1]=species_names


# Write FD madogram summary results for each specimen to a .csv
write.csv(FD_array, "Y:/xxx/xxx/xxx.csv")



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
# Save the FD_array as an R object
save(FD_array_box, file="Analysis/FD_box_array.R")


# Associate species names with each FD box counting complexity score
dimnames(FD_array_box)[1]=species_names


# Write FD box counting complexity scores for each specimen to a .csv
write.csv(FD_array_box, "Y:/xxx/xxx/xxx.csv")







