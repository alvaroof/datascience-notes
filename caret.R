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