# Title: Opening 3D skulls in R

# Name: Heather White

# Date created: 15/05/19

# Last modified: 15/11/19

# License: MIT license

# Script to load 3D skulls into R and export .png photos


#####################################################################################################

rm(list=ls())

library(Morpho)
library(rgl)

# Load .ply file
specimen_name=ply2mesh(file="file:///Y:/xxx/xxx/xxx/specimen_name.ply")

# Bring up the 3D visualiser
shade3d(specimen_name,col=bone1)
# Once open in the 3D visualiser, zoom in to the region of interest

# Then run this line of code to get the printscreen of the region of interest as .png format
rgl.snapshot("Y:/xxx/xxx/xxx/specimen_name.png")

