# Title: Suture complexity index

# Name: Heather White

# Date created: 01/11/19

# Last modified: 01/02/20

# License: MIT license

# Calculating suture complexity index using sinuosity index calculations

##################################################################################################


rm(list = ls())

library(tidyverse)


##################################################################################################

# Load in data


# Read in .csv for Sinuosity index results
SI_raw <- read.csv('Analysis/SI_results.csv', header = T)
# Selecting the SI column
SI <- select(SI_raw, V1)


# Read in complexity factor (CF) data
CF_raw <- read.csv('Analysis/CF_results.csv', header = T)
# Select the CF column
CF <- select(CF_raw, CF)


# Read in the species names .csv - make sure this is in the same order as the specimens
species_names <- read.csv('Data/Species_names_ordered.csv', header = F)


################################################################################################

# Calculating SCI

SCI <- SI * CF

# Save the species names with each SCI result
dimnames(SCI)[1]=species_names[1]


# Write SCI complexity score results to a csv
write.csv(SCI, "Y:/xxx/xxx/xxx.csv")



