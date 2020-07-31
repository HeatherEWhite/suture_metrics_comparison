# Title: Sinuosity index

# Name: Heather White

# Date created: 24/10/19

# Last modified: 31/10/19

# License: MIT license



# Script to calculate the sinuosity index of each specimen

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


# Sinuosity index

# Formula: length of suture following the curve / suture length start to end


######################################################################################################

# Part 1: This below part of the code calculates the top line of the sinuosity index equation


# Setting up the array for containing the landmark to landmark to landmark ... suture length measurements
ABlength_array<-array(dim=c(499,ntaxa))
# Calculating suture lengths (landmark to landmark) for all 79 taxa and adding these to an array 
b = 1
for(i in 1:length(csvlist))
{
  for(b in 1:499)
  {
    start <- csvarray[b,,]
    end <- csvarray[c(b+1),,]
    ABlength_array[c(b),i] <-  sqrt(((start[1,i] - end[1,i])^2) + (start[2,i] - end[2,i])^2)
    b= b + 1
  }
}



# Setting up an empty array for the sum of landmark to landmark length measurements
lengths_sum<-array(dim=c(1,ntaxa))
# Summing up the length values for between all the landmarks to get top hald of the SI equation
for(i in 1:length(csvlist))
{
  lengths_sum[,i] <- sum(ABlength_array[,i])
    
}


# Checking the summing up for specimen number 1
spec1 <- ABlength_array[,1]
spec1_sum <- sum(spec1)



#######################################################################################################

# Part 2: This below part of the code calculates the bottom line of the sinuosity index equation


# Creating arrays of 2D landmarks for the start and end landmarks of each of the 79 sutures
# Array in position 1 specifies x and y coordinates (therefore 2D landmarks)
# Array in position 2 specifies the specimen number (1 to 79)
start2 <- csvarray[1,,]
end2 <- csvarray[500,,]


# Calculating suture start to end length of suture 1 - test
# I am starting from landmark 3 not 1, as the sutures 
ABlength2 <- sqrt(((start2[1,1] - end2[1,1])^2) + (start2[2,1] - end2[2,1])^2)


# Setting up the array for containing A to B suture length measurement
ABlength_array2<-array(dim=c(1,ntaxa))
# Calculating suture lengths (start to end) for all 79 taxa and adding these to an array 
for(i in 1:length(csvlist))
{
  ABlength_array2[,i] <- sqrt(((start2[1,i] - end2[1,i])^2) + (start2[2,i] - end2[2,i])^2)
}




#####################################################################################################


# Calculating sinuosity index based by: Part1/Part2
SI <- lengths_sum/ABlength_array2

# Transposing the SI variable so that columns become rows
SI2 <- t(SI)


# Associating the species name with the correct matrix
# The number 1 is used to say associated the species names to the 3rd part of the array - in this case this is each matrix
dimnames(SI2)[1]=species_names


# Writing SI summary results for each specimen to a csv
write.csv(SI2, "Y:/Heather/R_Projects/2D_landmarks_sutures/Analysis/SI_results.csv")




