#####################################################################
##CARET
##date: 23/01/2018
##author: Alvaro O.
##description: general information and command examples for the caret
##library, extracted from the related coursera course
#####################################################################

library(caret)
library(e1071)

#predict functions
#lda :==: predict(obj)
#glm :==: predict(obj, type="response")
#gbm :==: predict(obj, type="response", n.trees)
#mda :==: predict(obj, type="posterior")
#rpart :==: predict(obj, type="prob")
#Weka :==: predict(obj, type="probability")
#LogitBoost :==: predict(obj, type="raw", nIter)

library(kernlab); data(spam)
inTrain <- createDataPartition(y=spam$type, p=0.75, list= FALSE)
training <- spam[inTrain,]
testing <- spam[-inTrain,]
dim(training)
set.seed(144)
modelFit <- train(type ~., data=training, method="glm")
modelFit
modelFit$finalModel #choose the best model of the multiple ones created

predTest <- predict(modelFit, newdata=testing)

confusionMatrix(predTest, testing$type)



##DAY3

##COVARIATES
#create dummy covariates:: from factor to indicator variables

library(ISLR);library(caret); data(Wage);
inTrain <- createDataPartition(y=Wage$wage, p=0.75, list= FALSE)
training <- Wage[inTrain,]
testing <- Wage[-inTrain,]

dummies <- dummyVars(wage ~ jobclass, data=training)
head(predict(dummies, newdata=training))

##removing zero covariates or covariates with low variability
nsv <- nearZeroVar(training, saveMetrics = TRUE)
nsv

#splines basis: fit curvy lines: polinomials
library(splines)
bsBasis <- bs(training$age,df=3) #polinomial grade 3
lm1 <- lm(wage ~ bsBasis, data=training)
plot(training$age, training$wage,pch=19,cex=0.5)
points(training$age, 
       predict(lm1, newdata=training), 
       col="red")

#to create the basis on test set you use
predict(bsBasis, age=testing$age)

##PCA

library(kernlab); data(spam)
inTrain <- createDataPartition(y=spam$type, p=0.75, list= FALSE)
training <- spam[inTrain,]
testing <- spam[-inTrain,]

#look for correlation
M <- abs(cor(training[,-58]))
diag(M) <- 0
which(M>0.8, arr.ind = T)

#PCA with caret
preProc <- preProcess(log10(spam[,-58]+1), 
                      method = "pca", pcaComp = 2)
spamPC <- predict(preProc, log10(spam[,-58]+1))
plot(spamPC[,1], spamPC[,2], col=spam$type)

trainPC <- predict(preProc, log10(training[,-58]+1))
modelFit <- train(training$type ~., 
                  method="glm", data=trainPC)

#we need to apply the preProp object that calculated with the train set
testPC <- predict(preProc, log10(testing[,-58]+1))


#diagnosis of models
inTrain <- createDataPartition(y=spam$type, p=0.75, list= FALSE)
training <- spam[inTrain,]
testing <- spam[-inTrain,]
dim(training)
set.seed(144)
modelFit <- train(type ~., data=training, method="glm")
finMod <- modelFit$finalModel

plot(finMod,1,pch=19,cex=0.5,col="red")
qplot(finMod$fitted.values, finMod$residuals, data=training)
plot(finMod$residuals,pch=19) #shouldn't be a trend here!!
#if there is a trend that means there is some relation with
#time or other variable not measured

##show several graphs at once. compare pairs
featurePlot(x = train[,vars], y = train$classe, plot="pairs")

#####################################################################
##CARET
#####################################################################

#For caret documentation please go to the script caret.R
library(caret)


#how to do model stacking
library(caret)
library(gbm)
set.seed(3433)
library(AppliedPredictiveModeling)
data(AlzheimerDisease)
adData = data.frame(diagnosis,predictors)
inTrain = createDataPartition(adData$diagnosis, p = 3/4)[[1]]
training = adData[ inTrain,]
testing = adData[-inTrain,]
set.seed(62433)

#three different models
mod_rf <- train(diagnosis ~ ., data = training, method = "rf")
mod_gbm <- train(diagnosis ~ ., data = training, method = "gbm")
mod_lda <- train(diagnosis ~ ., data = training, method = "lda")
pred_rf <- predict(mod_rf, testing)
pred_gbm <- predict(mod_gbm, testing)
pred_lda <- predict(mod_lda, testing)

#stacked using RF
predDF <- data.frame(pred_rf, pred_gbm, pred_lda, diagnosis = testing$diagnosis)
combModFit <- train(diagnosis ~ ., method = "rf", data = predDF)
combPred <- predict(combModFit, predDF)
confusionMatrix(combPred, testing$diagnosis)$overall[1]


##regularized models: LASSO
library(elasticnet)
set.seed(3523)
library(AppliedPredictiveModeling)
data(concrete)
inTrain = createDataPartition(concrete$CompressiveStrength, p = 3/4)[[1]]
training = concrete[ inTrain,]
testing = concrete[-inTrain,]
set.seed(233)
mod_lasso <- train(CompressiveStrength ~ ., data = training, method = "lasso")
library(elasticnet)
plot.enet(mod_lasso$finalModel, xvar = "penalty", use.color = TRUE)


#TRAINING TIME SERIES
library(forecast)
mod_ts <- bats(tstrain)
fcast <- forecast(mod_ts, level = 95, h = dim(testing)[1])
sum(fcast$lower < testing$visitsTumblr & testing$visitsTumblr < fcast$upper) / 
    dim(testing)[1]


##support vector machine
set.seed(3523)
library(AppliedPredictiveModeling)
data(concrete)
inTrain = createDataPartition(concrete$CompressiveStrength, p = 3/4)[[1]]
training = concrete[inTrain, ]
testing = concrete[-inTrain, ]
set.seed(325)
library(e1071)
mod_svm <- svm(CompressiveStrength ~ ., data = training)
pred_svm <- predict(mod_svm, testing)
accuracy(pred_svm, testing$CompressiveStrength)

##near zero variance
library(caret)
zeroVar <- nearZeroVar(train,saveMetrics=TRUE)
train <- train[,!zeroVar$nzv]
test <- test[,!zeroVar$nzv]

##compair pairs, plot
featurePlot(x = train[,vars], y = train$classe, plot = "pairs")

