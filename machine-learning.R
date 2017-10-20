#comments

Sys.setlocale("LC_ALL", "C") #Elimina problemas con lenguaje Inglés

tapply(data$ChildMortality, data$Region, min, na.rm = TRUE)

setwd("C:/Users/GGORS")
mvt <- read.csv("mvtWeek1.csv")

#DEFINING FUNCTIONS

myfunction <- function(arg1){
  table(arg)[1]/(table(arg)[2]+table(arg)[1])*100
  return(object)
}

#function for reading all the .csv in a folder

import.all.csv <- function(folder){
  setwd(folder)
  temp = list.files(pattern="*.csv")
  data <- lapply(temp, read.csv)
  for (i in seq(1:length(data))) {
    assign(temp[i], data[[i]])
  }
  return()
}




#De factor a Date con formato correcto
DateConvert = as.Date(strptime(mvt$Date, "%m/%d/%y %H:%M")) 
IBM$Date = as.Date(IBM$Date, "%m/%d/%y")

mvt$Month = months(DateConvert)
mvt$Weekday = weekdays(DateConvert)
mvt$Date = DateConvert
hist(mvt$Date, breaks=100)

Top5 <- subset(mvt, mvt$LocationDescription == "PARKING LOT/GARAGE(NON.RESID.)" | mvt$LocationDescription == "ALLEY" | mvt$LocationDescription == "STREET"| mvt$LocationDescription == "GAS STATION" | mvt$LocationDescription == "DRIVEWAY - RESIDENTIAL")
Top5$LocationDescription = factor(Top5$LocationDescription)
table(Top5$LocationDescription)


setwd("C:/Users/GGORS")
IBM <- read.csv("IBMStock.csv")
GE <- read.csv("GEStock.csv")
ProcterGamble <- read.csv("ProcterGambleStock.csv")
CocaCola <- read.csv("CocaColaStock.csv")
Boeing <- read.csv("BoeingStock.csv")
IBM$Date = as.Date(IBM$Date, "%m/%d/%y")
GE$Date = as.Date(GE$Date, "%m/%d/%y")
CocaCola$Date = as.Date(CocaCola$Date, "%m/%d/%y")
ProcterGamble$Date = as.Date(ProcterGamble$Date, "%m/%d/%y")
Boeing$Date = as.Date(Boeing$Date, "%m/%d/%y")

plot(CocaCola$Date, CocaCola$StockPrice, col="red")
lines(ProcterGamble$Date, ProcterGamble$StockPrice, col = "blue", lty=2)
abline(v=as.Date(c("1997-09-01")), lwd=2)
abline(h=100, lwd=2)
colors()
plot(CocaCola$Date[301:432], CocaCola$StockPrice[301:432], type="l", col="red", ylim=c(0,210))
lines(IBM$Date[301:432], IBM$StockPrice[301:432], col = "blue", lty=2)
lines(GE$Date[301:432], GE$StockPrice[301:432], col = "green", lty=2)
lines(Boeing$Date[301:432], Boeing$StockPrice[301:432], col = "black", lty=2)
lines(ProcterGamble$Date[301:432], ProcterGamble$StockPrice[301:432], col = "purple", lty=2)



tapply(IBM$StockPrice, months(IBM$Date), mean)


CPS = merge(CPS, MetroAreaMap, by.x="MetroAreaCode", by.y="Code", all.x=TRUE)
table(tapply(CPS$Race == "Asian", CPS$MetroArea, mean) > 0.20)




model1 <- lm(Price ~ AGST + HarvestRain, data=wine) #Método para regresión lineal
summary(model1) 
model1$residuals #Vector con los residuos
sum(model1$residuals^2) #SSE

cor(wine$WinterRain, wine$Price) #Correlacion entre dos variables
cor(wine) #Matriz de correlaciones (todas las variables han de ser num)
abs(cor(wine)) > 7 #Para ver las variables independientes con alta correlación


#Realizar ciones
model4 <- lm(Price ~ AGST + HarvestRain + WinterRain + Age, data = wine)
predictTest <- predict(model4, newdata=wineTest) #PRedicciones
predictTest

memory.limit() #checking max RAM usable according to conf

SSE <- sum((wineTest$Price-predictTest)^2)
SST <- sum((wineTest$Price-mean(wine$Price))^2)





PointsReg <- lm(PTS ~ X2PA + X3PA + FTA + AST + ORB + TOV + STL + BLK, data=NBA)
RMSE <- (sum(PointsReg$residuals^2/nrow(NBA))^(1/2))


#We can automatically simplify a starting model with the coman step()
model <- lm(Temp ~ MEI + CO2 + CH4 + N2O + CFC.11 + CFC.12 + TSI + Aerosols, data=train)
modelstep <- step(model)


#Determining which columns has Na's values
NAS <- function(x){ sum(is.na(x))}
data <- (apply(train, 2, NAS) > 0)
names(data[data == TRUE])

#remove NAs
train <- na.omit(train)
test <- na.omit(test)

#setting a diferent reference level in unordered factros

train$raceeth <- relevel(train$raceeth, "White")
test$raceeth <- relevel(test$raceeth, "White")


#comandos para el estudio del Google Flu Trends
FluTrain <- read.csv("FluTrain.csv")
hist(FluTrain$ILI, breaks = 100)

#It is very interesting to use log() when a variable is too SKEWED
plot(FluTrain$Queries, log(FluTrain$ILI))

#To see if it is skewed use an histogram
hist(FluTrain$ILI)

#when coming back from a prediction using log() we must use exp() function
PredTest1 <- exp(predict(FluTrend1, newdata=FluTest))

#relative error (Observed - Estimated)/Observed

#when working with factors is better to subset
data <- subset(FluTest, FluTest$Week == "2012-03-11 - 2012-03-17")


#working with time series
ILILag2 = lag(zoo(FluTrain$ILI), -2, na.pad=TRUE) #apply lag too zoo object
FluTrain$ILILag2 = coredata(ILILag2) #extract real values from zoo object

#For Displaying multiple plots at once
par(mfrow=c(1,2))


#do a split in a dataset into train and test sets
library(caTools)
set.seed(88)
split <- sample.split(framingham)
split <- sample.split(framingham$TenYearCHD, SplitRatio = 0.65)
train <- subset(framingham, split == TRUE)
framinghamLog <- glm(TenYearCHD ~., data=train, family = binomial) #family = binomial for logistic regression
predictTest <- predict(framinghamLog, type = "response", newdata=test)

#if the dependent variable is not discrete (e.g continous)
spl = sample(1:nrow(data), size=0.7 * nrow(data))
train = data[spl,]
test = data[-spl,]



#AUC value is the prob of differentiating between a randomly chosen positive or negative point 
library(ROCR)
ROCRpred <- prediction(predictTest, test$TenYearCHD)
as.numeric(performance(ROCRpred, "auc")@y.values) #AUC (Area under ROC Curve) value


#for plotting the ROC curve
ROCRpred <- prediction(predictTrain, train$PoorCare)
ROCRperf <- performance(ROCRpred, "tpr", "fpr")
plot(ROCRperf) #Different versions of the plot
plot(ROCRperf, colorize = TRUE)
plot(ROCRperf, colorize = TRUE, print.cutoffs.at = seq(0,1,0.1))
plot(ROCRperf, colorize = TRUE, print.cutoffs.at = seq(0,1,0.1), text.adj = c(-0.2,1.7))

#AUC measures how good is a model. 0.5 tells us that is no better than choosing randomly
#The AUC of a model has the following nice interpretation: given a random patient from 
#the dataset who actually received poor care, and a random patient from the dataset who actually 
#received good care, the AUC is the perecentage of time that our model will classify which is which
#correctly. 


#filling Nas
library(mice)
polling = read.csv("PollingData.csv")
simple <- polling[c("Rasmussen", "SurveyUSA", "DiffCount", "PropR")]
set.seed(144)
imputed <- complete(mice(simple))
polling$Rasmussen <- imputed$Rasmussen
polling$SurveyUSA <- imputed$SurveyUSA


#sign function for setting a baseline model in logistic regression
table(sign(train$Rasmussen))

#removing NAs
#split into train and test
#defining baselinemodel
#looking for multicollinearity

mod1 <- glm(Republican ~ PropR, data=train, family="binomial")
pred1 <- predict(mod1, type="response") #since we are predicting into train newdata is not required.
#type = "response" for probabilities
LogModel2 = glm(voting ~ sex + control + sex:control, data=gerber, family="binomial") #we can add
#the itneraction terms between the independent variables as shown above

#for not using some non-numerical variables into the models we can define a vector
nonvars = c("year", "songtitle", "artistname", "songID", "artistID")
train = train[ , !(names(train) %in% nonvars) ]
test = test[ , !(names(test) %in% nonvars) ]
mod2 = glm(Top10 ~ . - loudness, data=train, family=binomial) #if you want to substract numerical variables

#quita una variable de un set
vars.for.imputation = setdiff(names(loans), "not.fully.paid")

#Imputation 
library(mice)
set.seed(144)
vars.for.imputation = setdiff(names(loans), "not.fully.paid")
imputed = complete(mice(loans[vars.for.imputation]))
loans[vars.for.imputation] = imputed


#CART method
library(rpart)
library(rpart.plot)
StevensTree <- rpart(Reverse ~ Circuit + Issue + Petitioner 
                     + Respondent + LowerCourt + Unconst, data=train, method = "class", minbucket = 25)


#for knowing the first appearance of an element into a vector use match
summer = c(4,3,2,5,7)
match(4, summer)


#RANDOM FOREST
library(randomForest)
stevens <- read.csv("stevens.csv")
split <- sample.split(stevens$Reverse, SplitRatio = 0.7)
Train <- subset(stevens, split == TRUE)
Test <- subset(stevens, split == FALSE)
Train$Reverse <- as.factor(Train$Reverse) #the varible that is intented to be forecasted has to be a factor
Test$Reverse <- as.factor(Test$Reverse)
stevensForest <- randomForest(as.factor(Reverse) ~Circuit + Issue + 
Petitioner + Respondent + LowerCourt + Unconst, data=Train, nodesize = 25, ntree = 200, importance = TRUE)
PredictForest <- predict(stevensForest, newdata = Test)
a <- table(Test$Reverse, PredictForest)
acc <- (a[2,2]+a[1,1])/(a[2,2]+a[1,1]+a[1,2]+a[2,1])
acc

#RandomForest loses interpretability with respect to CART models, but we
#can gain insight on which variables were used thrpough the following

importance(model)
vu = varUsed(MODEL, count=TRUE)
vusorted = sort(vu, decreasing = FALSE, index.return = TRUE)
dotchart(vusorted$x, names(MODEL$forest$xlevels[vusorted$ix]))

#wecan also measure the "impurity" in our model
varImpPlot(MODEL)


#CROSS VALIDATION
library(caret)
library(e1071)
set.seed(2)
numFolds <- trainControl(method = "cv", number = 10)
cpGrid <- expand.grid(.cp=seq(0.01,0.5,0.01))
tr <- train(Reverse ~ Circuit + Issue + Petitioner + Respondent + LowerCourt + Unconst, data=Train, method="rpart", trControl = numFolds, tuneGrid = cpGrid)
bestcp <- tr$bestTune
best.tree <- tr$finalModel
plot(tr$results[,1], tr$results[,2]) #plots the RMSE - cp graph
StevensTreeCV <- rpart(Reverse ~ Circuit + Issue + Petitioner + Respondent + LowerCourt + Unconst, data=Train, method = "class", cp = bestcp)
PredictTestCV <- predict(StevensTreeCV, newdata = Test, type="class")


#CART with Penalty Matrix
claimsTree <- rpart(bucket2009 ~ age + arthritis + alzheimers + cancer + copd + depression + diabetes + heart.failure + ihd + kidney + osteoporosis + stroke + bucket2008 + reimbursement2008, data=train, method="class", cp=0.00005, parms = list(loss=PenaltyMatrix))



#solid point 
plot(boston$LON, boston$LAT)
points(boston$LON[boston$CHAS == 1],boston$LAT[boston$CHAS == 1], col = "red")
points(boston$LON[boston$CHAS == 1],boston$LAT[boston$CHAS == 1], col = "blue")
points(boston$LON[boston$CHAS == 1],boston$LAT[boston$CHAS == 1], col = "blue", pch = 19)



#Error in plot.new() : figure margins too large
par("mar") #should show [1] 5.1 4.1 4.1 2.1
par(mar=c(1,1,1,1))


#Randomly selects a sample from a train set

set.seed(1)
trainSmall = train[sample(nrow(train), 2000), ]


#getting all the dataframes from the workspace and put them into a list

l.df <- lapply(ls(), function(x) if (class(get(x)) == "data.frame") get(x))


#TEXT ANALYTICS
#1) Convert the title variable to corpusTitle and the abstract variable to corpusAbstract.
#2) Convert corpusTitle and corpusAbstract to lowercase. After performing this step, remember to run the lines:
#3) Remove the punctuation in corpusTitle and corpusAbstract.
#4) Remove the English language stop words from corpusTitle and corpusAbstract.
#5) Stem the words in corpusTitle and corpusAbstract (each stemming might take a few minutes).
#6) Build a document term matrix called dtmTitle from corpusTitle and dtmAbstract from corpusAbstract.
#7) Limit dtmTitle and dtmAbstract to terms with sparseness of at most 95% (aka terms that appear in at least 5% of documents).
#8) Convert dtmTitle and dtmAbstract to data frames (keep the names dtmTitle and dtmAbstract).

library(tm)
library(SnowballC)

tweets <- read.csv("tweets.csv", stringsAsFactors = FALSE) #use always this parameter
corpus <- Corpus(VectorSource((tweets$Tweet))) #define corpus with the text
corpus <- tm_map(corpus, tolower) #to lowercase
corpus = tm_map(corpus, PlainTextDocument)
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, c(stopwords("english")))
corpus <- tm_map(corpus, stemDocument)
corpus[[1]]$content
frequencies <- DocumentTermMatrix(corpus)
inspect(frequencies[1000:1005,505:515]) #inspecting the frequency matrix
findFreqTerms(frequencies, lowfreq = 20) #ranking of the wrods used more than 20 times
sparse <- removeSparseTerms(frequencies, 0.995) #keep terms that appear in more than 0.5% of documents
tweetsSparse <- as.data.frame(as.matrix(sparse)) #convert to DF
colnames(tweetsSparse) = make.names(colnames(tweetsSparse)) #make sure that names are apporpiate names. not numbers
tweets$Negative <- as.factor(tweets$Avg <= -1) #labeling tweets as negative
tweetsSparse$Negative <- tweets$Negative #adding new variable negative to DF
library(caTools) #splitting in Train and Test sets
set.seed(123)
split <- sample.split(tweetsSparse$Negative, SplitRatio = 0.7)
trainSparse <- subset(tweetsSparse, split == TRUE)
testSparse <- subset(tweetsSparse, split == FALSE)


strwrap(emails$email[1]) #shows a long string dividided per lines

#The grepl function returns TRUE if a string is found in another string, e.g.
grepl("cat","dogs and cats",fixed=TRUE) # TRUE
grepl("cat","dogs and rats",fixed=TRUE) # FALSE

#for declaring logical variables running tests on previous data
wikiWords2$HTTP <- ifelse(grepl("http", wiki$Added, fixed = TRUE), 1, 0)

#counting number of characters
which.min(nchar(trials$title))

#summing along rows or columns
colSums()
rowSums()


# Remove unnecessary variables
movies$ID = NULL
movies$ReleaseDate = NULL
movies$VideoReleaseDate = NULL
movies$IMDB = NULL

# Remove duplicates
movies = unique(movies)


#Obtener centroides. funciones sobre listas y grupos
spl = split(movies[2:20], clusterGroups)
lapply(spl, colMeans)

#HIERARCHICAL CLUSTERING with WARD and euclidean. WARD finds compact and spherical clusters
flower <- read.csv("flower.csv", header = FALSE)
flowerMatrix <- as.matrix(flower)
flowerVector <- as.vector(flowerMatrix)
distance <- dist(flowerVector, method = "euclidean")
clusterIntensity <- hclust(distance, method = "ward.D")
plot(clusterIntensity) #the taller the line, the more difference between clusters
rect.hclust(clusterIntensity, k=3) #remark 3 clusters
flowerClusters <- cutree(clusterIntensity, k=3) #define number of clusters. assgins each element to a cluster
#we are going to plot the image
dim(flowerClusters) <- c(50,50) #redim flower image as matrix
image(flowerClusters, axes = FALSE) # plot the image
image(flowerMatrix, axes = FALSE, col = grey(seq(0,1, length=256))) #see the original grayscale image

#K-MEANS CLUSTERING
k <- 10
set.seed(200)
KMC <- kmeans(householdsNorm, centers = k, iter.max = 1000)
str(KMC)
healthyClusters <- KMC$cluster
dim(healthyClusters) <- c(nrow(healthyMatrix), ncol(healthyMatrix))
image(healthyClusters, axes = FALSE, col = rainbow(k))
tumor <- read.csv("tumor.csv", header = FALSE) #reads test set
tumorMatrix <- as.matrix(tumor)
tumorVector <- as.vector(tumorMatrix)

#For obtaining a train and test set 
library(flexclust) #for test model on test set
KMC.kcca <- as.kcca(KMC, healthyVector) #kmeansclusteranalysis object
tumorClusters <- predict(KMC.kcca, newdata = tumorVector)
dim(tumorClusters) <- c(nrow(tumorMatrix),ncol(tumorMatrix))


table(articleClusters, KMC$cluster) #comparing HierClusters with KMCClusters

#normalizing variables to clustering
preproc = preProcess(limitedTrain)
normTrain = predict(preproc, limitedTrain)
normTest = predict(preproc, limitedTest) #remember that we use the same preproc for both!

#Visualization
library(ggplot2)
scatterplot <- ggplot(WHO, aes(x=GNI,y=FertilityRate))
scatterplot + geom_point()
scatterplot + geom_point(col = "blue", size=3, shape=17)
scatterplot + geom_point(col = "blue", size=3, shape=17) + ggtitle("Graphic")
FertilityGNIplot <- scatterplot + geom_point(col = "blue", size=3, shape=17) + ggtitle("Graphic")
pdf("pdf.pdf")
print(FertilityGNIplot)
dev.off()
gplot(WeekDayCounts, aes(x = Var1, y = Freq)) + geom_line(aes(group=1))
geom_line(aes(group=1), linetype = 12, alpha=0.3)


#split variable by color and size
FertilityGNIplot <- ggplot(WHO, aes(x=GNI, y=FertilityRate, color=Region, size=LifeExpectancy)) + geom_point() + ggtitle("Graphic")
FertilityGNIplot <- ggplot(WHO, aes(x=GNI, y=FertilityRate, color=Region)) + geom_point() + ggtitle("Graphic") 

#plotting a linear regression
ggplot(WHO, aes(x=log(FertilityRate), y=Under15)) + geom_point() + ggtitle("Graphic") + stat_smooth(method = "lm")

#more easy to see colors
scale_color_brewer(palette="Dark2")

#convert dates to suitable formats
mvt$Date <- strptime(mvt$Date, format="%m/%d/%y %H:%M")
mvt$Weekday <-weekdays(mvt$Date)

#Reordering the values of a variable into a data frame
WeekDayCounts$Var1 <- factor(WeekDayCounts$Var1, ordered = TRUE, levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))

#HEATMAP
ggplot(DayHourCounts, aes(x=Hour,y=Var1)) + geom_tile(aes(fill=Freq))
                                                      



                                                      
----------------------------------
  

#READ from file que multiple byte separator
dat <- readLines('test.csv')
dat <- gsub("sep", " ", dat)
dat <- textConnection(dat)
dat <- read.table(dat)
                                                    
                                                      
                                                    
#Join tables
library(plyr)
df <- join(movies, ratings, type = "inner")

#Creating Matrix with variables and values (mapped) 
ratingmat <- dcast(ratings, UserID~MovieID, value.var = "Rating", na.rm=TRUE)
ratingmat[is.na(ratingmat)] <- 0
ratingmat <- as.matrix(ratingmat)

#aggregate data frames, summarizing some variable
wordmat <- aggregate(tags$tag, by = list(tags$movieId), paste, collapse = " ") #aggregate all the tags per movie
names(wordmat) <- c("movieID", "tag")

#only correlating numeric variable
cor(train[sapply(train, is.numeric)])


#PCA Principal Component analysis
# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
split <- sapply(train,is.numeric)
train_numeric <- train[,split]
ir.pca <- prcomp(train_numeric,
                 center = TRUE,
                 scale. = TRUE) 
# print method
print(ir.pca)
# plot method
plot(ir.pca, type = "l")
# summary method
summary(ir.pca)
# Predict PCs
predict(ir.pca, 
        newdata=tail(log.ir, 2))

#PCA plots
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

#PCA con caret
require(caret)
trans = preProcess(train_numeric, 
                   method=c("BoxCox", "center", 
                            "scale", "pca"))
PC = predict(trans, train_numeric)

#PCA con princomp
arc.pca1 <- princomp(train_numeric, scores=TRUE, cor=TRUE)
summary(arc.pca1)
plot(arc.pca1)
biplot(arc.pca1)
arc.pca1$scores #print the projection

#PLOT pairwise scatterplots. exploratory data analysis
pairs(data)

#Exploratory data analysis. categorical. numerical
library(ggplot2)
qplot(quantity, data=train, geom="bar", fill=status_group) + 
  theme(legend.position = "top") + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0))

ggplot(train, aes(x = construction_year)) + 
  geom_histogram(bins = 20) + 
  facet_grid( ~ status_group)
#CARET preprocessing

require(caret)
trans = preProcess(data, 
                   c("BoxCox", "center", "scale"))
predictorsTrans = data.frame(
  trans = predict(trans, data))

#SMOTE with DMwR. SMOTE is for ovresamplying categories with less 15%
library(DMwR)
train$status_group <- as.factor(train$status_group)
train <- SMOTE(status_group ~ ., train, perc.over = , perc.under=100)

#Confusion Matrix
library(caret)
confusionMatrix(pred_forest_train, train$status_group)


#compare variables two by two, to see relations. margin 1 rows margin 2 columns
prop.table(table(train$payment, train$status_group), margin = 1)

#join data into one variaable
train$latlong <- paste(round(train$latitude,2), round(train$longitude, 2), sep = ":")

#relates a bin number to each entry
.bincode(train$Purchase, c(0,5000,10000,15000,20000, 30000))

#PCA 

#PCA
train.numeric <- subset(train, select = sapply(train,is.numeric))
test.numeric <- subset(test, select = sapply(test,is.numeric))

pca.train <- prcomp(train.numeric, scale. = T)
pca.test <- prcomp(test.numeric, scale. = T)

#compute standard deviation of each principal component
pca.data <- pca.test
std_dev <- pca.data$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained",type = "b")
plot(cumsum(prop_varex), xlab = "Principal Component",ylab = "Cumulative Proportion of Variance Explained",type = "b")
