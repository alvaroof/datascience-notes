########################################################################################################
########################################################################################################
### Title: Summary of commands in R (Coursera)
### Date: 20/06/2017
########################################################################################################
########################################################################################################

.libPaths("C:/Datos/RStudio/library")

#HELP
?rnorm
help(rnorm)
help.search("rnorm")
args(rnorm)
rnorm
search() #list elemnts into the environment

#OS
getwd()
setwd()
dir() #list all files in wdundeb            
ls() #list objects in R workingspace
list.files("Directory name") #return a list with all names fot files in that folder
file.exist("Directory name") #returns TRUE or FALSE
dir.create("Directory Name") #create directory in that path

#DATATYPES
vector("numeric", length = 10) #create vectors
a <- list("a", 2:5, 1+01)
attributes(a)
matrix(1:6, nrow = 2, ncol = 3)
#"another way of creating a matrix is editing el attrbute dimension of a vector
a <- 1:10
dim(a) <- c(2,5)
#Other standar method is using cbind() and rbind()
fact <- factor(c("yes","yes","no", "maybe"), levels = c("maybe","yes","no")) #you can set the order of factors
unclass(fact)

is.na()
is.nan()
typeof()
class()
#LOAD DATA
read.table()
read.csv()
readLines()
source()
dget()
load()
unserialize()

a <- array(rnorm(2*2*10), c(2,2,10)) #generates 10 2x2 matrices

#to automatically initialize the classes of each column (optimizes reading time)
initial <- read.table("data.txt", nrows = 10)
classes <- sapply(initial, class)
final <- read.table("data.txt", colClasses = classes)

#To roughly calculate the size of a dataset in memory (RAM) -- rowsxcolumnsx(storage/datatype) [8 bytes/numeric]
#As a rule of thumb, it will be required twice as much memory
#e.g. 1.5m rows 120columns numeric is gonna require 1.34GBx2 = 2.7Gb of RAM

#EXPORT DATA
write.table()
writeLines() #write each element to a line in a text file
dump() #same as dump(), with several objects. preserve metadata, output is a text file editable, less affected by corrruption, 
       #better to track with version control software
dput() #same as dump(). export an R object
save()
serialize()

#CONNECTIONS 
file()
gzfile()
bzfile()
url()

#same as read.csv("foo.txt")
con <- file("foo.txt","r")
data <- read.csv(con)
close(con)

#reading HTML
con <- url("https://www.jhsph.edu/", "r")
html <- readLines(con)
head(html)

#PASTE AND CONCATENATE
paste()
cat() #paste strings a put them out to a file

#DATES AND TIMES
x <- as.Date("1970-01-01") #day Zero in R
p <- Sys.time()
strptime() #Parse dates in other format

#LOOP FUNCTIONS
x <- 1:5
lapply(x, runif, min = 0, max = 10)
sapply(x, runif, min = 0, max = 10) #does the same as laaply but return 
#something simplified (not in this specific case...)
tapply() #over subset of an object

#All this apply a function over the elements of the same object
mapply() #this does that over different objects

split() 

#split in conjunction with sapply is really useful
library(datasets)
s <- split(airquality, airquality$Month)
sapply(s, function(x){colMeans(x)})

#gl() generate factos to use with split
f <- gl(3,10)

#interaction() serves for creating sublevels
f1 <- gl(2,5)
f2 <- gl(5,2)
interaction(f1,f2)
str(split(x, list(f1,f2), drop = TRUE)) #split doesn'0t need interaction to be called
#drop = TRUE takes out the emptly combind levels

##DEBUGGING
lm(x~y)
traceback() #tell you how different functions where called when you execute one
    #it gives information about the most recent error
debug(lm) #typical debugging, let you go through the code of a function executing it one step at a time
lm(x~y)
# Browse[2]> n 
#for goint to he next line

browser() #same as debug() but from any point inside the function 
trace() #let you debug something without editing it

options(error = recover)
read.csv("mouchfile") #you can change the default error behaviour (throwing a message) to
        #freeze the execution wherever the error ocurred


##SIMULATION
#useful functions for normal distributions (analog to other dist e.g. poisson)
rnorm() #random generation
dnorm() #density 
pnorm() #cumulative
qnorm() #quantile

set.seed(144) #in simulation is important to set seeds to we can reproduce results!

sample(1:10, replace = TRUE)

##PROFILING
system.time(str(mtcars)) #gives the elapsed execution time for the input expression
Rprof() #starts Rprofiler
summaryRprof() #summarizes output Rprofiler by.self is the best option

##DPLYR::::package to manipulate dataframes
select(iris, Sepal.Length) #return subset of columns
iris2 <- filter(iris, Sepal.Length > 5) #return subset of rows
iris2 <- arrange(iris, Sepal.Length, desc(Petal.Length)) #reorder rows
iris2 <- rename(iris, SLength = Sepal.Length, PLength = Petal.Length) #rename some columns
iris2 <- mutate(iris2, sdLen = (SLength - PLength)) #create new columns
iris2 <- mutate(iris2, Sepallong = factor(1 * (SLength > 5),labels = c("short", "long"))) #also useful for defining
#labeled subgroups
species <- group_by(iris2, Species)
summarize(species, SLength = mean(SLength)) #Reduces multiple values down to a single value

#Also, pipeline operator (%>%) is pretty handy to perform operations in one line
iris %>% mutate(L = Sepal.Length - Petal.Length, W = Sepal.Width - Petal.Width) %>% 
    group_by(Species) %>% summarize(LDev = sum(L - mean(L)), WDev = sum(W - mean(W)))


#INDEXES
match() #returns happeing index of an element in within an array



###GETTING AND CLEANING DATA
##DOWNLOAD FILES and REad FILES
download.file() #downloads a file from the internet. Useful for downloading csv 
#files, excel files, and so on.
#method = "curl" for https
dateDownloaded <- date() ##good to know when the file was downloaded
read.xlsx()
read.xlsx2() #faster than normal version, but unstable to read subset of rows
write.xlsx() 

#XML Files
library(XML)
fileUrl <- "URL for XML file"
doc <- xmlTreeParse(fileUrl, useInternal=TRUE)
rootNode <- xmlRoot(doc)
xmlName(rootNode) @shows
names(rootNode) #elements nested right below the wrapper
rootNode[[1]] #you can access element like in a list
rootNode[[1]][[1]]
xmlSApply(rootNode,xmlValue) #Sapply for XML (extracts aLLelements)
xmlSApply(rootNode,"//node1",XmlValue) #extracts node value for all elements
xmlSApply(rootNode,"//node2",XmlValue) 
xmlSApply(doc, "//li[@class='team-name']", xmlValue)

#JSON files
library(jsonlite)
jSon <- "http://by0nja.de.bayer.cnb/roi/services/generalServices/getDashboardDetails.php?dashboardName=PDCM%20Dashboard%20New"
jsonData <- fromJSON(jSon)
names(jsonData$data)

#convert files to JSON
myjson <- toJSON(iris, pretty = TRUE)
cat(myjson)

#data.table datatables (really fast, should give it a try when working with large files)
library(data.table) #works as data.frame, but faster. creates objects of 
#class table

##EXPLORATORY DATA ANALYSIS
##GRAPHS
summary(iris)

boxplot(iris$Sepal.Length ~iris$Species, col = "green")
abline(h = 6.75)

hist(iris$Sepal.Length, col = "red", breaks = 10)
abline(v = mean(iris$Sepal.Length), lwd = 5, lty = 2)
rug(iris$Sepal.Length)

par("lty") #~says the default value for that characteristic
colors()
par(las = 1, bg = "red", mar = c(4,4,2,1), oma = c(0,0,0,0), mfrow = c(2,1))
barplot(table(iris$Species), col = "blue")


par(mfrow = c(3,1), mar = c(4,4,2,1))
hist(subset(iris, iris$Species == "setosa")$Sepal.Length)
hist(subset(iris, iris$Species == "versicolor")$Sepal.Length)
hist(subset(iris, iris$Species == "virginica")$Sepal.Length)
par(mfrow = c(1,1), mar = c(4,4,2,1))

plot(iris$Sepal.Length, iris$Sepal.Width, col = iris$Species, pch = 12)
legend("topright", pch = 1, col = c("black", "red", "green"),legend = unique(iris$Species))

lines(iris$Petal.Length, iris$Petal.Width)
points(iris$Petal.Length, iris$Petal.Width)
title()
text()
mtext()
axis()

example(points) #to see examples of use

##PLOTTING
#base plotting system (ABOVE)
#plot regression line
iris2 <- subset(iris, iris$Species == "virginica")
plot(iris2$Sepal.Width, iris2$Sepal.Length, pch = 20)
model <- lm(iris2$Sepal.Length ~ iris2$Sepal.Width)
abline(model)

#multiple plots
par(mfrow = c(1,3), mar = c(4,4,2,1), oma = c(0,0,2,0))
with(airquality, {
    plot(Wind, Ozone, main = "Ozone ~ Wind")
    plot(Solar.R, Ozone, main = "Ozone ~ Solar.R")
    plot(Temp, Ozone, main = "Ozone ~ Temp")
    mtext("Ozone and Weather in NY", outer = TRUE)
})

#lattice system
library(lattice)
xyplot(iris$Sepal.Length ~ iris$Sepal.Width | iris$Species, layout = c(1,3))
bwplot(iris$Sepal.Length ~ iris$Species)
#ggplot2
library(ggplot2)
qplot(iris$Sepal.Length, iris$Sepal.Width)

##GRAPHIC DEVICES
#types: bitmap (for images), vector (for graphs), svg (supports interactiviry, webbased)
#bitmap (png (line drawing, solid colors), jpeg (natural scenes), tiff, bmp (icons))
#to send to a pdf
pdf(file = "myplot.pdf") 
par(mfrow = c(1,3), mar = c(4,4,2,1), oma = c(0,0,2,0))
with(airquality, {
plot(Wind, Ozone, main = "Ozone ~ Wind")
plot(Solar.R, Ozone, main = "Ozone ~ Solar.R")
plot(Temp, Ozone, main = "Ozone ~ Temp")
mtext("Ozone and Weather in NY", outer = TRUE)
})
dev.off()


dev.cur() #currently open graphic devices 


##MYSQL
library(RMySQL)
library(sqldf)
#see databases in server
ucscDb <- dbConnect(MySQL(),user="genome",host="genome-mysql.cse.ucsc.edu")
result <- dbGetQuery(ucscDb, "show databases;"); 
dbDisconnect(ucscDb)

#connect to specific database
hg19 <- dbConnect(MySQL(),user="genome", db="hg19", host="genome-mysql.cse.ucsc.edu")
allTables <- dbLisTables(hg19) 
length(allTables)

#see the fields of a particular table
dbListFields(hg19, "table")

#convert table to a dataframe
affyData <- dbReadTable(hg19, "table")
head(affyData)

#select a specfici subset
query <- dbSendquery(hg19, "select * from affyU133Plus2 where blablabla")
affyMis <- fetch(query); quantile(affyMis$misMatches)
affyMisSmall <- fetch(query, n=10)
dbClearResult(query)

##ACCESING WEBS
library(httr)
library(XML)
url <- "https://www.bioconductor.org/install/"
html2 <- GET(url)
content2 <- content(html2, as="text")
parsedHtml <- htmlParse(content2, asText=TRUE)
xpathSApply(parsedHtml, "///title", xmlValue)

#for authentications
html2 <- GET(url, authenticate("user","passwd"))


##READING IMAGES
install.packages("png")
install.packages("jpg")
install.packages("readbitmap")

##READING GIS
install.packages("rdgal")
install.packages("rgeos")
install.packages("raster")

##READING MUSIC DATA
install.packages("tuneR")
install.packages("seewave")


##visualization with lattice
#panel function serves to personalize what is shown
xyplot(airquality$Ozone ~ airquality$Wind | airquality$Month, panel = function(x,y,...){
    panel.xyplot(x,y,...)
    panel.lmline(x,y, col = 2)
})

##GITHUB APIs
library(httr)
# 1. Find OAuth settings for github:
#    http://developer.github.com/v3/oauth/
oauth_endpoints("github")
# 2. To make your own application, register at at
#    https://github.com/settings/applications. Use any URL for the homepage URL
#    (http://github.com is fine) and  http://localhost:1410 as the callback url
#
#    Replace your key and secret below.

myapp <- oauth_app("github",
                   key = "cbdbff315748875caa0f",
                   secret = "d622f9e28f3262b3d0e633be2a10c1c2093db573")

# 3. Get OAuth credentials
github_token <- oauth2.0_token(oauth_endpoints("github"), myapp)

# 4. Use API
gtoken <- config(token = github_token)
req <- GET("https://api.github.com/users/jtleek/repos", gtoken)
stop_for_status(req)
content(req)

# OR:
req <- with_config(gtoken, GET("https://api.github.com/rate_limit"))
stop_for_status(req)
content(req)

##COMPARISONS
identical(a1r, a1s)

##read files with width fixed
x <- read.fwf(
  file=url("http://www.cpc.ncep.noaa.gov/data/indices/wksst8110.for"),
  skip=4,
  widths=c(12, 7, 4, 9, 4, 9, 4, 9, 4))

##VISUALIZACION CON GGPLOT2
#qplot with two variables is similar to plot
qplot(iris$Sepal.Length, iris$Sepal.Width, 
      color = iris$Species, 
      geom = c("point", "smooth"))

#qplot with one variable is similar to hist
qplot(iris$Sepal.Length, fill = iris$Species)
qplot(iris$Sepal.Length, geom = "density", color = iris$Species)

#facets is similar to the lattice spliting by factor
qplot(iris$Sepal.Length, iris$Sepal.Width, 
      data = iris,
      facets = . ~ iris$Species)

#to plot linear models
qplot(iris$Sepal.Length, iris$Sepal.Width, 
      color = iris$Species) + 
geom_smooth(method = "lm")

