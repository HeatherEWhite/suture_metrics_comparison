# Title: StereoMorph: collecting 2D landmarks

# Name: Heather White

# Date created: 21/05/19

# Last modified: 23/05/19

# License: MIT license

# StereoMorph for collecting 2D curves from jpeg suture photos

##################################################################################################

rm(list=ls())

library(StereoMorph)
library(geomorph)


# Initialising the app for positioning landmarks
digitizeImages(image.file = 'Raw_Data/Photos', shapes.file = 'Raw_Data/Original_Landmarks', landmarks.ref = 'Raw_Data/Landmarks.txt', curves.ref = 'Raw_Data/Curves.txt')

# Saving the original landmark data in the Analysis folder
shapedat<-readShapes("./Raw_Data/Original_Landmarks")

# Saving the resampled shape data, so there are the same number of semilandmarks for each curve of each suture
shapedat_resampled<-readland.shapes(shapedat, c(500))

# Saving the resampled landmarks to a variable for each specimen
# This has to be run each time with every specimen name
Tragelaphus_scriptus_resampled <- shapedat_resampled[["landmarks"]][["Tragelaphus_scriptus"]]

# Writing the resampled (500) landmark data to a csv file and saving it in the 'Data' folder rather than the raw data folder
write.csv(Tragelaphus_scriptus_resampled, file = 'Data/Resampled_landmarks/Tragelaphus_scriptus_resampled.csv')
