#-------------------------------------------
# Useful packages to install
#-------------------------------------------

#For time series
install.packages("zoo")

#For using ARIMA model
install.packages("forecast")

#Curvas ROC
install.packages("ROCR")

#imputation anf filling of missing values
install.packages("mice")

#For using sample.split()
install.packages("caTools")

#for plots
install.packages("ggplot2")
install.packages("ggmap")

#for CART
install.packages("rpart")
install.packages("rpart.plot")
install.packages("randomForest")
install.packages("class")
install.packages("caret")
install.packages("e1071")

#text analytics
install.packages("tm")
install.packages("SnowballC")

#K-MEANS testin
install.packages("flexclust")

#MCA Analysis
install.packages("FactoMineR")

#JSON
install.packages("jsonlite")

#GET (URL)
install.packages("RCurl")

#Recommenderlab
install.packages("recommenderlab")

#Web Services
install.packages("Rook")
install.packages("googleVis")

#SHINY APP
install.packages('rsconnect')
install.packages("RSQLite")
install.packages(c('shiny', 'ggvis', 'dplyr', 'RSQLite'))


#Data Manipulation
install.packages("plyr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("reshape")
install.packages("data.table")
install.packages("splitstackshape")
install.packages("reshape2")


#SQL Package
install.packages("sqldf")

#Biclustering
install.packages("biclust")

#recommender package
install.packages("recommenderlab")

#EXCEL
install.packages("xlsx")

#GRadient Boosting
install.packages("xgboost")

#instalar devtools
install.packages("devtools")

#githhub packages
library(devtools)
install_github("ggbiplot", "vqv")

#Data mining with R. SMOTE
install.packages("DMwR")

#KNN
install.packages("FNN")


#XGBoost
install.packages("readr")
install.packages("drat", repos="https://cran.rstudio.com")
drat:::addRepo("dmlc")
install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")


#R Google Analytics API 
#install devtools package for downloading packages from github
install.packages("devtools")
library(devtools)
#install curl for easier use
install.packages("curl")
library(curl)
#installing rga package from github
install_github("skardhamar/rga")

#Excel files
install.packages("XLConnect")

#XML Files
install.packages("XML")

##MySQL
install.packages("RmySQL")

#Bioconductor
source("https://www.bioconductor.org/biocLite.R")
biocLite("rhdf5")
biocLite("impute")

#Access a lot of foreign data types
install.packages("foreign")

#clustering categorical variables
install.packages("klaR")
