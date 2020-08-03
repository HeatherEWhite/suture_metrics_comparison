# Title: StereoMorph: collecting 2D landmarks

# Name: Heather White

# Date created: 21/05/19

# Last modified: 23/11/19

# License: MIT license

# StereoMorph for collecting 2D curves from jpeg suture photos

##################################################################################################

rm(list=ls())

library(StereoMorph)
library(geomorph)


# Initialising the app for positioning landmarks
digitizeImages(image.file = 'Raw_Data/Photos', shapes.file = 'Raw_Data/Original_Landmarks', landmarks.ref = 'Raw_Data/Landmarks.txt', curves.ref = 'Raw_Data/Curves.txt')

# Save the original landmark data
shapedat<-readShapes("./Raw_Data/Original_Landmarks")

# Save the resampled shape data - each suture resampled to have 500 semi-landmarks along its path
shapedat_resampled<-readland.shapes(shapedat, c(500))

# Save the resampled landmarks to a variable for each specimen
# This has to be run each time with every specimen name
Tragelaphus_scriptus_resampled <- shapedat_resampled[["landmarks"]][["specimen_name"]]

# Write the resampled (500) landmark data to a .csv file
write.csv(specimen_name_resampled, file = 'Y:/xxx/xxx/xxx.csv')



