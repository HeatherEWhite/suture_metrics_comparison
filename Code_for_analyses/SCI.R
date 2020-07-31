# Title: Suture complexity index

# Name: Heather White

# Date created: 01/11/19

# Last modified: 01/11/19

# License: MIT license

# Calculating suture complexity index using sinuosity index calculations

##################################################################################################


rm(list = ls())

library(tidyverse)


##################################################################################################

# Loading in data


# Read in Sinuosity index csv
SI_raw <- read.csv('Analysis/SI_results.csv', header = T)
# Selecting the SI column
SI <- select(SI_raw, V1)


# Read in the CF data
CF_raw <- read.csv('Analysis/CF_results.csv', header = T)
# Selecting the CF column
CF <- select(CF_raw, CF)


# Reading in the species names .csv - this is in the same order as the specimens
species_names <- read.csv('Data/Species_names_ordered.csv', header = F)


################################################################################################

# Calculating SCI

SCI <- SI * CF

# In SCI varaible get rownames to save as species names
# Have to use the function dimnames not rownames because it is an array
# The [1] accesses the first position of the array, when you do dim(FD_array) this gives you 79, 1. So using [1] renames the species names
dimnames(SCI)[1]=species_names[1]


# Writing FD summary results for each specimen to a csv
write.csv(SCI, "Y:/Heather/R_Projects/2D_landmarks_sutures/Analysis/SCI_results.csv")



