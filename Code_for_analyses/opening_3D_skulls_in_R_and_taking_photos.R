# Title: Opening 3D skulls in R

# Name: Heather White

# Date created: 15/05/19

# Last modified: 15/05/19

# License: MIT license

# Script to load 3D skulls into R and export .png photos


#####################################################################################################

rm(list=ls())

library(Morpho)
library(rgl)

# Load file
Giraffa_camelopardalis=ply2mesh(file="file:///Y:/Heather/Results/Suture Methods/Suture Samples Test Dataset ASCII/Giraffa_camelopardalis.ply")

# Bring up the 3D visualiser
shade3d(Giraffa_camelopardalis,col=bone1)
# Once I have the 3D visualiser up, zoom in using mouse scroller to region of interest

# Run this line of code to get the printscreen of the region of interest as .png format
rgl.snapshot("Y:/Heather/Results/Suture Methods/Suture Samples Test Dataset Photos/Giraffa_camelopardalis.png")

