# Title: Calculating STFT and PSD

# Name: Heather White

# Date created: 05/06/19

# Last modified: 20/02/19

# License: MIT license



# Script to calculate the short-time Fourier transform and the subsequent power spectrum density

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(e1071)
library(psd)
library(geomorph)
library(Momocs)



#######################################################################################################

# Read in the .csv LM files into an array:


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



# Associating a specimen name with each array
dimnames(csvarray)[3]<-list(
  substr(dir("./",pattern=".csv"),1,(nchar(dir("./",pattern=".csv"))-4)))###gives specimens names based on file names in ply folder
lm_array<-csvarray



# Checking the array:

# Check the dimensions of lm_array
dim(lm_array)

# Check the class of lm_array - for fractal dimension this has to be matrix or dataframe
class(lm_array)



# Read in the species names .csv - make sure this is in the same order as the specimens
species_names <- read.csv('Data/Species_names_ordered.csv', header = F)



#######################################################################################################

# Scaling the data (to centroid size and line length)


# Performing Procrustes'
Procrustes<-gpagen(lm_array)
# Saving coordinate data to a new variable once normalised
lm_array_corrected<-Procrustes$coords


#######################################################################################################

# STFT Analysis


STFT_array <- array(dim=c(39,64,79)) 
# 64 is the number of columns, determined by the number of coefficients and 39 is the number of rows for each coefficient
#FD_summary <- array(79,9)
for(i in 1:length(csvlist))
{
  run<-stft(lm_array_corrected[,,i], win=min(80,floor(length(lm_array_corrected[,,i])/10)), 
            inc=min(24, floor(length(lm_array_corrected[,,i])/30)), coef=64, wtype="hanning.window")
  STFT_array[,,i]<-run$values
}
# STFT results:
# rows = windows used in STFT; columns = Fourier coefficients (64) for each window


# Checking STFT_array has worked for a random specimen, number 78 in this case
STFT_array[,,78]


# Associating the species name with the correct matrix
dimnames(STFT_array)[3]=species_names


# Save STFT_array in R
save(STFT_array, file='Y:/xxx/xxx/xxx.R')


# Write STFT_array to a csv
write.csv(STFT_array, 'Y:/xxx/xxx/xxx.csv')


# Convert and save the array to a 2D matrix
STFT_2D_array <- two.d.array(STFT_array, sep='.')
save(STFT_2D_array, file='Y:/xxx/xxx/xxx.R')



#############################################################################################################

# Power Spectrum Density (PSD)


# Calculate the average of the squared STFT coefficients over each frequency/local transforms
# This calculates the Power value for each harmonic 
STFT_average <- array(dim=c(39,1,79)) 
for(i in 1:length(csvlist))
{
  STFT_average[,,i] <- rowMeans(STFT_array[,,i]^2, na.rm = FALSE, dims = 1)
}

# Viewing the results
STFT_average



# Sum the Power values at each harmonic to get a power value for each suture
PSD_array <- array(dim = c(79,1))
for(i in 1:length(csvlist))
{
  PSD_array[i,] <- sum(STFT_average[,,i])
}

# View the PSD values for each suture
PSD_array

# Associate the species name with the correct matrix
dimnames(PSD_array)[1]=species_names


# Save PSD results to a .csv
write.csv(PSD_array, "Y:/xxx/xxx/xxx.csv")




