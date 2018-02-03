########################################################################################################
########################################################################################################
### Title: Summary of commands in R (Coursera)
### Date: 20/06/2017
########################################################################################################
########################################################################################################

.libPaths("C:/Datos/RStudio/library")
setwd("C:/datos")
source("C:/Datos/Coursera/get_info.R")

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
object.size(iris)

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

sample(1:10, replace = TRUE) #create random vector from a sample

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


##########################################################
###GETTING AND CLEANING DATA
#############################################################

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

table(iris$Sepal.Length, useNA = "ifany") #count number of occurences
xt <- xtabs(iris$Sepal.Length ~ iris$Sepal.Width + iris$Petal.Width ) #make crosstabs
ftable(xt) #better visualization


##check for missing values
sum(is.na(iris$Petal.Width))
any(is.na(iris$Petal.Width))
all(iris$Sepal.Length>0)
colSums(is.na(iris))
all(colSums(is.na(iris)))
table(iris$Species %in% c("versicolor", "virginica"))


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
install.packages("jpeg")
install.packages("readbitmap")

library("jpeg")
download.file("https://d396qusza40orc.cloudfront.net/getdata%2Fjeff.jpg", destfile = "image.jpeg")
img <- readJPEG("image.jpeg", native = TRUE) ##native = TRUE reads the jpeg in 
#nativeRaster format


##READING GIS
install.packages("rdgal")
install.packages("rgeos")
install.packages("raster")

##READING MUSIC DATA
install.packages("tuneR")
install.packages("seewave")



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


#Know the size of a file
object.size(iris)

##CREAITNG NEW VAR
#creating sequences
s1 <- seq(1,10,by=2)
s2 <- seq(1,100,length = 4)

##creating binary variables
iris$large <- ifelse(iris$Sepal.Length>mean(iris$Sepal.Length), TRUE,FALSE)
table(iris$large)

#creating categorical variables
cutpoints <- quantile(iris$Sepal.Length, seq(0,1,length=5), na.rm = TRUE)
iris$groups <- cut(iris$Sepal.Length, cutpoints)
table(iris$groups)


#RESHAPING DATASETS

#melting dataframes: simplify datasets setting which var are ids(keys) and which are measures
library(reshape2)
mtcars$carname <- rownames(mtcars)
carMelt <- melt(mtcars, id=c("carname", "gear", "cyl"), measure.vars = c("mpg", "hp")) #melt reshape the dataframe
head(carMelt, n=3)

#casting dataframes: reshape dataframes to summarize one var in functions of others. By default summarize by length
cylData <- dcast(carMelt, cyl ~ variable, mean)

#averaging values
tapply(InsectSprays$count, InsectSprays$spray, sum)

#split into a list
l <- split(iris$Sepal.Length, iris$Species)
lapply(l, sum)
unlist(l)

#intersect of two group
intersect(names(iris), c("Species"))

#merging data: like doing joins in SQL
library(dplyr)
colClasses = c("character", "character", "character", "character",
               "character", "character", "character", "character", "character", "character")

gdp <- read.csv(url("https://d396qusza40orc.cloudfront.net/getdata%2Fdata%2FGDP.csv"), 
                skip = 4, colClasses = colClasses)
edu <- read.csv(url("https://d396qusza40orc.cloudfront.net/getdata%2Fdata%2FEDSTATS_Country.csv"))
gdp <- rename(gdp, Id = X)
edu <- rename(edu, Id = CountryCode)

df <- merge(gdp, edu, by = "Id") #slower than join, but more featured
df2 <- merge(gdp, edu, by = "Id") #better for multiple dfs

df <- rename(df, gdp = X.4)
df <- mutate(df, gdp = as.numeric(gsub(",","",gdp)))
df <- select(df, gdp, Income.Group)
gdp_groups <- group_by(df, Income.Group)
summarize(gdp_groups, mean = mean(gdp, na.rm = TRUE))

df2$rank <- as.numeric(df2$X.1)
df2 <- select(df2, rank, Income.Group)
rank_groups <- group_by(df2, Income.Group)
summarize(rank_groups, mean = mean(rank, na.rm = TRUE))

cutpoints <- quantile(df2$rank, seq(0,1,length=5), na.rm = TRUE)
df2$quantile <- cut(df2$rank, cutpoints)

##EDITING TEXT VARIABLES
txt <- c("AlVAro", "FOO", "F Fo fef ,,,//", "Variable.Name1")
tolower(txt) #make all letters lowercase
toupper(txt) #make all letters uppercase
strsplit(txt, "\\.") #split strings by pattern
sub("F","_therewasanFhere_",txt) #replace first ocurrence
gsub("F","_therewasanFhere_",txt) #replace all matches
grep("FOO", txt) #return position of elements where match occurs
grep("FOO", txt, value=TRUE) #returns the actual elements
grepl("FOO", txt) #return boolean vector
nchar(txt) #number of characters
substr(txt, 1,3) #three first letters
paste(txt[1],txt[2],sep = "_pasting_") #combines two strings
paste0(txt[1],txt[2]) #pastes without space
trimws(txt) #removes leading and trailing spaces

##REGULAR EXPRESSIONS
# ^:beginning of line
# $:end of line
# [Bb][Uu][Ss][Hh]: set of characters
# [0-9][a-z][A-Z][0-9a-z]: range of characters
# [^a-z]: NOT a-z
# .:any single character
# blood|fire: or character (match blood or fire)
# ( and ): as in programming
# [a-z]?: this expression is 'optional'
# (.*): any number of repetitions of any character, including none (*? is ungreedy)
# .+: same as before but at least one repetition
#{m,n}: restrict matches to m to n occurences
# {m,} at least m matches

##WORKING WITH DATES
library(lubridate) #convert any numbers to date. there are functions for times as well
dat <- "14::03::12"
ydm(dat)
dmy(dat)
#ultimately yiou want your dates to be of class POSIXlt or pOSIXct




#####################################################################
##EXPLORATORY DATA ANALYSIS
#####################################################################

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

pairs(iris)
plot(iris[,1:4]) #multi plot


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

##visualization with lattice
#panel function serves to personalize what is shown
xyplot(airquality$Ozone ~ airquality$Wind | airquality$Month, panel = function(x,y,...){
    panel.xyplot(x,y,...)
    panel.lmline(x,y, col = 2)
})

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

##how to use generic ggplot

g <- ggplot(iris, aes(Sepal.Length, Sepal.Width))

g + geom_point(aes(color = Species), size = 5) +
    labs(title = "ggplot2 example") +
    labs(x = expression("log " * PM[2.5]), y = "Width") +
    geom_smooth(method = "lm") + 
    facet_grid(. ~ Species) +
    theme_bw() +
    coord_cartesian(ylim = c(2.5,4.0))


cutpoints <- quantile(iris$Petal.Length, seq(0,1,length=5), na.rm = TRUE)
iris$Petal.factor <- cut(iris$Petal.Length, cutpoints)
g <- ggplot(iris, aes(Sepal.Length, Sepal.Width))
g + geom_point(aes(color = Species)) +
    facet_grid(Petal.factor ~ Species) +
    theme_bw() 

#boxplot
library(UsingR)
data(InsectSprays)
g <- ggplot(data=InsectSprays, aes(y = InsectSprays$count, x = InsectSprays$spray, color=InsectSprays$spray))
g + geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

#for density histograms
xvalues <- seq(min(sdistribution), max(sdistribution), length=50)
yvalues <- dnorm(xvalues, mean=1/lambda, sd=(1/lambda/sqrt(n)))
d1 <- data.frame(x = sdistribution)
d2 <- data.frame(x = xvalues, y= yvalues)
g <- ggplot()
#here is the important line ..density..
g <- g + geom_histogram(aes(x=x, y=..density..), color="black", fill="antiquewhite2", binwidth=0.5, data = d1)
g <- g + geom_line(aes(x=x, y=y), colour="red", data=d2)
g <- g + labs(title = "Comparison to a normal distribution with same mean and sd") + theme_bw() 


##HIERARCHYCAL CLUSTERING
#let's generate and annotate random data points 
set.seed(1234)
par(mar = c(0,0,0,0))
x <- rnorm(12, mean = rep(1:3, each = 4), sd = 0.2)
y <- rnorm(12, mean = rep(c(1,2,1), each = 4), sd = 0.2)
plot(x,y, col="blue", pch=19,cex=2)
text(x+0.05,y+0.05, labels = as.character(1:12))

dataFrame <- data.frame(x=x, y=y)
distxy <- dist(dataFrame) #by default is the euclidean distance metric

#hierarchical clustering
#good thing about hclust is that it doesn't tell you how many cluster there are
#this is something that you define 'cutting' the tree
#when calculating distance between cluster there are two common methods
#average linkage: distance between center of gravity
#complete linkage: distance is the distance between farthest points.
#advice: try to scale variables before applying
hClustering <- hclust(distxy)
plot(hClustering)

#heatmap function for visualizing matrix data
#what heatmap basically does is to run a hierarchical clustering
#both in the rows and in the columns of the matrix
set.seed(144)
dataMatrix <- as.matrix(dataFrame)[sample(1:12),]
heatmap(dataMatrix)


##K-MEANS CLUSTERING
#contrary to hierarchical clustering, k menas needs to know the number
#of clusters in advance.
#let's generate and annotate random data points 
set.seed(1234)
par(mar = c(0,0,0,0))
x <- rnorm(12, mean = rep(1:3, each = 4), sd = 0.2)
y <- rnorm(12, mean = rep(c(1,2,1), each = 4), sd = 0.2)
plot(x,y, col="blue", pch=19,cex=2)
text(x+0.05,y+0.05, labels = as.character(1:12))
dataFrame <- data.frame(x=x, y=y)

kmeansObj <- kmeans(dataFrame, centers = 3)
names(kmeansObj)
#clusters has the clustering of the points
plot(x,y,col = kmeansObj$cluster, pch = 19, cex = 2)
points(kmeansObj$centers, col = 1:3, pch = 3, cex = 3, lwd = 3)

#we can also see clustering from k-means in a heatmap-style
#using image() function
set.seed(144)
dataMatrix <- as.matrix(dataFrame)[sample(1:12),]
kmeansObj2 <- kmeans(dataMatrix, centers = 3)
par(mfrow = c(1,2), mar = c(2,4,0.1,0.1))
image(t(dataMatrix)[,nrow(dataMatrix):1]) #ordered as in the matrix
image(t(dataMatrix)[,order(kmeansObj2$cluster)]) #ordered witihn clusters


##DIMENSION REDUCTION
#exploratory pattern analysis
set.seed(144)
par(mar = rep(4,4))
dataMatrix <- matrix(rnorm(400), nrow = 40)
image(1:10,1:40, t(dataMatrix))[,nrow(dataMatrix):1] 
heatmap(dataMatrix) #shows heatmap

#reorder rows according to clusters
hh <- hclust(dist(dataMatrix))
dataMatrixOrdered <- dataMatrix[hh$order, ]
par(mfrow = c(1,3))
image(1:10,1:40, t(dataMatrixOrdered))[,nrow(dataMatrixOrdered):1] 
plot(rowMeans(dataMatrixOrdered), 40:1,xlab = "Row Mean", ylab = "Row", pch = 19)
plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column mean", pch = 19)
#we don't see any pattern as the numbers were uniformly randomly generated

#PRINCIPAL COMPONENT ANALAYSIS // SINGULAR VARIABLE DECOMPOSITION
#find a new set of fewer uncorrelated variables that explains
#as much variance as possible of the original data
#SVD is finding a decomposition of matrix X as: X = U*D*V(transpose)
#where U is called left singular vectors
#V is called right singular vectors
#D is singular values
#PCA is the same as SVD, but before decomposing the mean is substracted from
#each column and then each one is divided by the standar deviation
set.seed(144)
dataMatrix <- matrix(rnorm(400), nrow = 40)
svd1 <- svd(scale(dataMatrix))
par(mfrow = c(1,3))
image(t(dataMatrix)[,nrow(dataMatrix):1])
plot(svd1$u[,1], ylab = "1st left singular vector", xlab = "Rows", pch = 19)
plot(svd1$v[,1], ylab = "1st right singular vector", xlab = "Columns", pch = 19)
#D matrix in decomposition comprehends the amount of variance explained
#by each component.
par(mfrow = c(1,1))
plot(svd1$d^2/sum(svd1$d^2), xlab = "Columns or Components",
     ylab = "Percentage of total variance explained", pch = 19)

#Calculate PCA1 is similar to 1st right singular vector ans so on
set.seed(144)
dataMatrix <- matrix(rnorm(400), nrow = 40)
pca1 <- prcomp(dataMatrix, scale = TRUE)

#Creating an artifical pattern (actually two) to prove the previous methods
set.seed(144)
for(i in 1:40){
    #flip a coin
    flip1 <- rbinom(1, size = 1, prob = 0.5)
    flip2 <- rbinom(1, size = 1, prob = 0.5)
    #if coin is head, add a common pattern to that row
    if(flip1){
        dataMatrix[i,] <- dataMatrix[i,] + rep(c(0,5), each = 5)
    }
    if(flip2){
        dataMatrix[i,] <- dataMatrix[i,] + rep(c(0,5), each = 5)
    }
}
hh <- hclust(dist(dataMatrix))
dataMatrixOrdered <- dataMatrix[hh$order,]


#handling missing values in svd or pca
library(impute)
set.seed(144)
dataMatrixNA <- matrix(rnorm(400), nrow = 40)
dataMatrixNA[sample(1:100, size = 40, replace = FALSE)] <- NA
#use k-nearest neighbors to do the imputation of NA values in rows
dataMatrix <- impute.knn(dataMatrixNA, k = 5)$data 


##USING COLORS IN R
#good for interpolate between colors
palette <- colorRamp(c("red", "green")) #◙extremes of the defined palette
palette(seq(0,1,len=10))
palette2 <- colorRampPalette(c("red", "green"))
palette2(10)

#good for choosing colors
require("RColorBrewer")
newcolors <- brewer.pal(3, "BuGn") #look for palette names in help page
newcolors
newpalette <- colorRampPalette(newcolors)
image(volcano, col = newpalette(5))

#for plotting many many points
#option1
x <- rnorm(10000)
y <- rnorm(10000)
smoothScatter(x,y)
#option2
plot(x,y,col=rgb(0,0,0,0.2)) #last parameter is alpha parameter


#rgb function allows you to choose exactly the color you want


#####################################################################
##REPRODUCIBLE RESEARCH
#####################################################################

sessionInfo() #this is good for keeping reproducibility across PCs
set.seed() #mandatory for random simulations


#####################################################################
##STATISTICAL INFERENCE
#####################################################################

#CONFIDENCE INTERVALS
library(datasets); data("ChickWeight"); library(reshape2)
wideCW <- dcast(ChickWeight, Diet + Chick ~ Time, value.var = "weight")

wideCW$gain <- wideCW['21'] - wideCW['0']

wideCW14 <- subset(wideCW, wideCW$Diet %in% c(1,4))

#calculate confidence interval
t.test(gain ~ Diet, paired = FALSE, var.equal = TRUE, data = wideCW14)$conf
#paired represents if the subject are the same individuals in the two groups
#var.equal represents if the variance is the same in both groups

#HYPOTHESIS TESTING
data(father.son)
t.test(father.son$sheight-father.son$fheight)
#to perform a t hypothesis test, in this case, that the height of both father 
#and son is the same

#to caluclate p-values,since we want the upper tails we use the lower.tail=FALSE
pbinom(6,8,0.5,lower.tail = FALSE)

#different and independent groups via T test (equal variances)
nx=10
ny=10
sp <- sqrt((9*0.60^2+9*0.68^2)/(10+10-2))
5-3+ c(-1,1)*qt(0.975,18)*sp*(1/10+1/10)^0.5

#different and independent groups via T test (unequal variances)
6-4+c(-1,1)*qnorm(0.975)*(0.5^2/100+2^2/100)^0.5


##POWER AND TYPE II ERROR
library(manipulate)
power.t.test(delta = 2,sd = 1, sig.level = 0.05, 
             power = 0.95, type = "one.sample", alt = "one.sided") #if you omit one argument it solves for it
#╢in power, the only parameter important is the effecgt size = (mu2-mu1)/sd


##MULTIPLE COMPARISONS
#we need to adjust out significance level or our pvalues if we do multiple tests
#like in the comic from xkcd about jelly beans and p<0.05...
#in this script we adjust the p values. for notes about adjusting the alpha values 
#please refer to the course slides


#we generate 1000 random datasets
set.seed(42)
pValues <- rep(NA,1000)
for(i in 1:1000){
    y <- rnorm(20)
    x <- rnorm(20)
    pValues[i] <- summary(lm(y~x))$coeff[2,4]
}
alpha <- 0.05

#we can see that we have 57 false positives (about 5%)
sum(pValues < alpha)
#and now we test different kinds of control types

#Controls FWER - Bonferroni correction - 
#very conservative - no false positives but some false negatives
sum(p.adjust(pValues, method = "bonferroni")<alpha)

#Controls FDR - Benjamini Hochberg correction
#less conservative than bonferroni - less false positives
sum(p.adjust(pValues, method = "BH")<alpha)


##RESAMPLING
#bootstrap
library(UsingR)
data(father.son)
x <- father.son$sheight
n <- length(x)
B <- 10000
resamples <- matrix(sample(x,n*B, replace = TRUE), B, n)
resampledMedians <- apply(resamples, 1, median)
#if we plot the last object we will see the distribution for the empirical median
g = ggplot(data.frame(medians = resampledMedians), aes(x= medians))
g = g + geom_histogram(color="black", fill="lightblue", binwidth = 0.05)
g 
sd(resampledMedians)
quantile(resampledMedians, c(0.025,0.975)) #to calculate the confidence interval
#there are some corrections to this in the package bootstrap


##permutation tests - group comparisons
##there are several variations for permutation tests
library(UsingR)
data(InsectSprays)
subdata <- InsectSprays[InsectSprays$spray %in% c("B","C"),]
y <- subdata$count
group <- as.character(subdata$spray)
testStat <- function(w,g){mean(w[g=="B"])-mean(w[g=="C"])}
observedStat <- testStat(y,group)
permutations <- sapply(1:10000, function(i) testStat(y, sample(group)))
observedStat
mean(permutations > observedStat) #as this is zero we reject the null hypothesis

#####################################################################
##REGRESSION 
#####################################################################
#least squares for the function
myfun <- function(mu){
    x <- c(0.18, -1.54, 0.42, 0.95)
    w <- c(2, 1, 3, 1)
    sum(w*(x-mu)^2)
}
optim(0.5, myfun)$par

#regression through the origin
x <- c(0.8, 0.47, 0.51, 0.73, 0.36, 0.58, 0.57, 0.85, 0.44, 0.42)
y <- c(1.39, 0.72, 1.55, 0.48, 1.19, -1.59, 1.23, -0.65, 1.49, 0.05)
lm(I(y) ~ 0+ I(x)) ##the 0 or minus 1 is used in R formulas to declare that the intercept is cero
#This I() statement inside the lm function serves to do arithmetical
#operations insede the function call


#residuals
library(UsingR)
data(diamond)
y <- diamond$price; x<-diamond$carat; n<- length(y)
fit <- lm(y~x)
e <- resid(fit) ##same as y-predict(fit)
yhat <- predict(fit)
sum(e) ##sum of residuals is zero (if model has intercept)
sum(e*x) ##has also to be zero

##plot points, fit and residuals
plot(x,y,bg="lightblue", col="black", pch=21, frame=FALSE)
abline(fit, lwd=2)
for(i in 1:n){
    lines(c(x[i],x[i]),c(y[i],yhat[i]), col="red", lwd=2)
}

#let's plot residuals versus x in horizontal axis, easier to evaluate
plot(x,e,bg="lightblue", col="black", pch=21, frame=FALSE)
abline(h=0)
for(i in 1:n){
    lines(c(x[i],x[i]),c(e[i],0), col="red", lwd=2)
}

##and now with ggplot
g=ggplot(data.frame(x=x,y=resid(lm(y~x))), aes(x=x,y=y))
g=g+geom_hline(yintercept = 0, size=2)
g=g+geom_point(size=7, colour="black", alpha=0.4)
g=g+geom_point(size=5, colour="red", alpha=0.4)
g=g+xlab("X")+ylab("Residual")
g #→it is good for residuals to appear to be randomly distributed . If we detect a pattern
#, some tendency, it means our model doesn't account for some variations of the dependent var


#♥here we compare the variation in residuals against the mean, and agains a lm model
e = c(resid(lm(y~1)), resid(lm(y~x)))
fit = factor(c(rep("Int", length(y)), rep("int,slope", length(y))))
g <- ggplot(data.frame(e=e,fit=fit), aes(y=e,x=fit,fill=fit))
g=g+geom_dotplot(binaxis = "y", size=2, stackdir = "center", binwidth = 30)
g


#variability in residuals
#if we include intercept the sum of residuals has to be zero
#standar deviation of residuals is 
sd <- sqrt((1/(n-2))*sum(resid(lm(y~x))^2))
#or
fit <- lm(y~x)
summary(fit)$sigma


#Inference in regression: how to extract the coefficients
library(UsingR)
data(diamond)
y <- diamond$price; x<-diamond$carat; n<- length(y)
fit <- lm(y~x)
coefs <- summary(fit)$coefficients
#build the confidence interval for the intercept and the slope
coefs[1,1]+c(-1,1)*qt(.975,df=fit$df)*coefs[1,2]
coefs[2,1]+c(-1,1)*qt(.975,df=fit$df)*coefs[2,2]

##PREDICTION
data(mtcars)
n <- nrow(mtcars)
fit <- lm(mpg~wt, data=mtcars)
sigmahat <- sqrt((1/(n-2))*sum(resid(fit)^2))
pred <- predict(fit, newdata=data.frame(wt=meanwt), interval="prediction")
conf <- predict(fit, newdata=data.frame(wt=meanwt), interval="confidence")
# conf is for the confidence interval of the regression line at that point
# pred is for thje conf interval of the predicted value itself
# in this example is calculed at the mean


###REGRESSION MODELS
##harmonic: Fourier Basis just for sins
notes4 <- c(261.63,293.66,329.63,349.23,392.00,440.00,493.88,523.25)
t <- seq(0,2, by=.001); n<- length(t)
c4 <- sin(2*pi*notes4[1]*t)
e4 <- sin(2*pi*notes4[3]*t)
g4 <- sin(2*pi*notes4[6]*t)
chord <- c4+e4+g4+rnorm(n,0,0.3)
x <- sapply(notes4,function(freq) sin(2*pi*freq*t))
fit <- lm(chord ~ x -1)
plot(notes4,fit$coefficients^2, type = "l", col = "red") #we can see the notes used
abline(v=notes4)
#fourier discrete transform: both for sins and cosins
#this shows the spectrum for the chord
a <- fft(chord); plot(Re(a)^2, type = "l")

#odds ratio of using autolander comparing headwind to tailwind
library(MASS)
shuttle$auto <- as.integer(shuttle$use == "auto")
fit <- glm(auto ~ wind - 1 , "binomial", shuttle)
cf <- exp(coef(fit))
oddsratio <- cf[1]/cf[2]
oddsratio

#poisson regression
fit <- glm(count ~spray-1, poisson, InsectSprays)
cf <- exp(fit$coefficients)
cf[1]/cf[2]

#poisson regression with an offset
n <- nrow(InsectSprays)
deltat <- 0.1
t <- seq(0.1,(n-1)*deltat+0.1,by = deltat)
fit <- glm(count ~spray+offset(log(t)), poisson, InsectSprays)
fit2 <- glm(count ~spray+offset(log(10)+log(t)), poisson, InsectSprays)

#knot methodology // spline
knots <- c(0)
splineTerms <- sapply(knots, function(knot) (x>knot)*(x-knot))
xMat <- cbind(1, x, splineTerms)
fit <- lm(y ~ xMat - 1)
yHat <- predict(fit)

#stimated slope
maxPos <- which.max(x)
minPos <- which(x == 0)
len <- maxPos - minPos
(yHat[maxPos] - yHat[minPos]) / len

