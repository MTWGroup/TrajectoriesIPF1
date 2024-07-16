#Project: IPF new project
#
# Purpose: Profile project: imputations RF MCMC Cens2021May with visit 4
# Version:5F
# Date: 24/06/2021
# Author:  HPF
#
# Input: 
# Main data bases:MDataH_CensMay21_230621_415oX32v
#
# Output: Monte carlo imputations  
#        imputations All visits with visits 1 and 3 and 4
#         Detection of error in data collection
#         Data Preparation and Reprocessing
#         Split the dataset into training and validation
#         Descriptive statistics
#         Training and train control
#         Make imputations on the test data
#         Compute the prediction error RMSE and deviation
#         Graphs RMSE and errorDev
#         Tables
#        
# variables to impute: Vs1==ppFVC @ baseline
#                      Vs3==ppFVC @ 90 days   
#                      Vs4==ppFVC @ 180 days
#                      Vs5==ppFVC @ 365 days
#                      Vs6==ppFVC @ 730 days
#                      Vs7==ppFVC @ 1095 days
# 
# Dependencies:@/F:PROFILE FINAL RF SC 130421
#
# Notes: use the new dataset, could be a problem with the dates.
#
# = 1 working space ============================================================

setwd("F:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
list.files("F:/")
setwd("F:/S3.1 PROFILE FINAL RF_MCMC 040621")
list.files("F:/S3.1 PROFILE FINAL RF_MCMC 040621")

# = 2 packages need ============================================================

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("sva"))
#install.packages("eeptools") 

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)

#Formatting/documentation packages
library(dplyr)
library(stringr)
library(tidyr)
library(plyr)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(Metrics)

#  K-Nearest Neighbors Essentials
library(caret)
library(skimr)
library(RANN)
library(DMwR)
library(FNN)
library(ranger)
library(tidyverse)
library(e1071)
library(rpart)
library(rpart.plot)
library(doSNOW)
library(randomForest)
library(tree)
library(MASS)
library(doParallel)

# = 3 open database ============================================================

MDW<- read.csv("MDataH_CensMay21_230621_415oX32v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
rownames(MDW)<-MDW$patientid
MDW$X=NULL

SIM<-read.csv("MDataV_Simulated_030621.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(SIM)
dim(SIM)
SIM$X=NULL

# = 4 data preparation and selection Vs4 =======================================

MDWVs4<-data.frame(MDW[ ,c("Vs1","Vs3","Vs4")])
sapply(MDWVs4, function(x) sum(is.na(x)))

MDWTemp<-na.omit(MDWVs4)
sapply(MDWTemp, function(x) sum(is.na(x)))

MDWVs3<-MDWVs4%>%filter(is.na(Vs3)) 
MDWVs3<-rbind(MDWVs3,MDWTemp)
MDWVs3<-mutate(MDWVs3,NAVs3Imp = !is.na(MDWVs3$Vs3))
summary(MDWVs3)

MDWVs4<-MDWVs4%>%filter(is.na(Vs4)) 
MDWVs4<-rbind(MDWVs4,MDWTemp)
MDWVs4<-mutate(MDWVs4,NAVs4Imp = !is.na(MDWVs4$Vs4))

SIM_Vs3<-data.frame(SIM[ ,c("Vs1","Vs3","Vs4","Simulated")])
SIM_Vs3<-SIM_Vs3%>%filter(Simulated=="SimX100")
SIM_Vs3$Simulated=NULL

# = 6 split the dataset into training and validation visit 3 ===================

fvcTrain<-SIM_Vs3

fvcTest<-MDWVs3
fvcTest$NAVs3Imp=NULL

# = 7 Computing forest model classifier FVC ====================================

cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs3~., 
                data=fvcTrain, 
                method="rf",
                tuneGrid = tunegrid,
                trControl=control)

#Shutdown cluster
stopCluster(cl)

print(rf_fit)
plot(rf_fit)
fvcImp <- varImp(rf_fit, scale = FALSE)
plot(fvcImp, top = 3)

FVCTrain_pred<- data.frame(rf_fit %>% predict(fvcTrain))
fvcTrain<-data.frame(fvcTrain,FVCTrain_pred$rf_fit.....predict.fvcTrain.)
names(fvcTrain)[4]<-paste("Vs3pred")

FVCTest_pred <- data.frame(rf_fit %>% predict(fvcTest))
fvcTest<-data.frame(fvcTest,FVCTest_pred$rf_fit.....predict.fvcTest.)
names(fvcTest)[4]<-paste("Vs3pred")

# = 8 Errors in the predictions and graphs =====================================

TrainVs3<-fvcTrain %>% dplyr:: summarise(n = n(),max=max(Vs3),min=min(Vs3),sd= sd(Vs3))
TestVs3<-fvcTest %>% dplyr:: summarise(n = n(),max=max(Vs3, na.rm=TRUE),min=min(Vs3,na.rm=TRUE),sd= sd(Vs3, na.rm=TRUE))

FVCTrain.RMSEVs3<-round(rmse(fvcTrain$Vs3pred,fvcTrain$Vs3)/(TrainVs3$sd),2)

fvcTestTemp<-na.omit(fvcTest)
FVCTest.RMSEVs3<-round(rmse(fvcTestTemp$Vs3pred ,fvcTestTemp$Vs3)/(TestVs3$sd),2)

MDWVs3<-rbind(fvcTest)
MDWVs3<-mutate(MDWVs3,NAVs3Imp = !is.na(MDWVs3$Vs3))
MDWVs3<- mutate(MDWVs3, Vs3= ifelse(is.na(Vs3), Vs3pred, Vs3))
rownames(MDWVs3)->MDWVs3$patientid

MDWVs3$ErrorDev<-round(((MDWVs3$Vs3pred - MDWVs3$Vs3)/(MDWVs3$Vs3)*100),2)

# = 9 Computing forest model classifier FVC ====================================

fvcTrain<-SIM_Vs3
fvcTest<-MDWVs4
fvcTest$NAVs4Imp=NULL

cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs4~., 
                data=fvcTrain, 
                method="rf",
                tuneGrid = tunegrid,
                trControl=control)

#Shutdown cluster
stopCluster(cl)

print(rf_fit)
plot(rf_fit)
fvcImp <- varImp(rf_fit, scale = FALSE)
plot(fvcImp, top = 3)

FVCTrain_pred<- data.frame(rf_fit %>% predict(fvcTrain))
fvcTrain<-data.frame(fvcTrain,FVCTrain_pred$rf_fit.....predict.fvcTrain.)
names(fvcTrain)[4]<-paste("Vs4pred")

FVCTest_pred <- data.frame(rf_fit %>% predict(fvcTest))
fvcTest<-data.frame(fvcTest,FVCTest_pred$rf_fit.....predict.fvcTest.)
names(fvcTest)[4]<-paste("Vs4pred")


# = 10 Errors in the predictions and graphs ====================================

TrainVs4<-fvcTrain %>% dplyr:: summarise(n = n(),max=max(Vs4),min=min(Vs4),sd= sd(Vs4,na.rm=TRUE))
TestVs4<-fvcTest %>% dplyr:: summarise(n = n(),max=max(Vs4, na.rm=TRUE),min=min(Vs4,na.rm=TRUE),sd= sd(Vs4, na.rm=TRUE))

FVCTrain.RMSEVs4<-round(rmse(fvcTrain$Vs4pred,fvcTrain$Vs4)/(TrainVs4$sd),2)

fvcTestTemp<-na.omit(fvcTest)
FVCTest.RMSEVs4<-round(rmse(fvcTestTemp$Vs4pred ,fvcTestTemp$Vs4)/(TestVs4$sd),2)

MDWVs4<-rbind(fvcTest)
rownames(MDWVs4)->MDWVs4$patientid

MDWVs4<-mutate(MDWVs4,NAVs4Imp = !is.na(MDWVs4$Vs4))
MDWVs4<- mutate(MDWVs4, Vs4= ifelse(is.na(Vs4), Vs4pred, Vs4))
MDWVs4$ErrorDev<-round(((MDWVs4$Vs4pred - MDWVs4$Vs4)/(MDWVs4$Vs4)*100),2)

# = 11 Join Vs3 and Vs4 ========================================================

MDWVs3T<-data.frame(MDWVs3[ ,c("patientid","Vs1","Vs3","Vs4")])
MDWVs3TT<-data.frame(MDWVs3[ ,c("patientid","Vs3pred","NAVs3Imp","ErrorDev")])

MDWVs4T<-data.frame(MDWVs4[ ,c("patientid","Vs1","Vs3","Vs4")])
MDWVs4TT<-data.frame(MDWVs4[ ,c("patientid","Vs4pred","NAVs4Imp","ErrorDev")])

MDWVs4T<-rbind(MDWVs3T,MDWVs4T)
MDWVs4T<-(MDWVs4T%>%filter(!duplicated(patientid)))

MDWVs4T<-plyr::join(MDWVs4T,MDWVs3TT, by="patientid")
MDWVs4T<-plyr::join(MDWVs4T,MDWVs4TT, by="patientid")
names(MDWVs4T)[7]<-paste("ErrorDev3")
names(MDWVs4T)[10]<-paste("ErrorDev4")

MDWVs4<-(MDWVs4T%>%mutate(ErrorDev3=replace_na(ErrorDev3,0))%>%mutate(ErrorDev4=replace_na(ErrorDev4,0))
          %>%mutate(NAVs3Imp=replace_na(NAVs3Imp,"TRUE"))%>%mutate(NAVs4Imp=replace_na(NAVs4Imp,"TRUE"))
          %>%mutate(Vs3pred= ifelse(is.na(Vs3pred), Vs3, Vs3pred))
          %>%mutate(Vs4pred= ifelse(is.na(Vs4pred), Vs4, Vs4pred)))

MDWVs4<-data.frame(MDWVs4[ ,c("patientid","Vs1","Vs3","Vs4","Vs3pred","Vs4pred","NAVs3Imp","NAVs4Imp","ErrorDev3","ErrorDev4")])


# = 12 start imputation visit V5 ===============================================
# = 12.1 data preparation and selection Vs5 ====================================

MDWVs5<-data.frame(MDW[ ,c("patientid","Vs5")])
MDWVs5<-plyr::join(MDWVs4,MDWVs5, by="patientid")
MDWVs5<-data.frame(MDWVs5[ ,c("patientid","Vs1","Vs3","Vs4","Vs5")])
sapply(MDWVs5, function(x) sum(is.na(x)))
rownames(MDWVs5)<-MDWVs5$patientid

MDWVs5$patientid=NULL

SIM_Vs5<-data.frame(SIM[ ,c("Vs1","Vs3","Vs4","Vs5","Simulated")])
SIM_Vs5<-SIM_Vs5%>%filter(Simulated=="SimX100")
SIM_Vs5$Simulated=NULL

# = 13 split the dataset into training and validation ==========================

fvcTrain<-SIM_Vs5
fvcTest <-rbind(MDWVs5)

# = 14 Computing forest model (using ranger) classifier FVC ====================
cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs5~., 
                data=fvcTrain, 
                method="rf",
                tuneGrid = tunegrid,
                trControl=control)

#Shutdown cluster
stopCluster(cl)

print(rf_fit)
plot(rf_fit)
fvcImp <- varImp(rf_fit, scale = FALSE)
plot(fvcImp, top = 3)

FVCTrain_pred<- data.frame(rf_fit %>% predict(fvcTrain))
fvcTrain<-data.frame(fvcTrain,FVCTrain_pred$rf_fit.....predict.fvcTrain.)
names(fvcTrain)[5]<-paste("Vs5pred")

FVCTest_pred <- data.frame(rf_fit %>% predict(fvcTest))
fvcTest<-data.frame(fvcTest,FVCTest_pred$rf_fit.....predict.fvcTest.)
names(fvcTest)[5]<-paste("Vs5pred")

# = 15 Errors in the predictions and graphs Vs5 ================================

TrainVs5<-fvcTrain %>% dplyr:: summarise(n = n(),max=max(Vs5),min=min(Vs5),sd= sd(Vs5))
TestVs5<-fvcTest %>% dplyr:: summarise(n = n(),max=max(Vs5, na.rm=TRUE),min=min(Vs5, na.rm=TRUE),sd= sd(Vs5, na.rm=TRUE))

FVCTrain.RMSEVs5<-round(rmse(fvcTrain$Vs5pred,fvcTrain$Vs5)/(TrainVs5$sd),2)

fvcTestTemp<-na.omit(fvcTest)
FVCTest.RMSEVs5<-round(rmse(fvcTestTemp$Vs5pred ,fvcTestTemp$Vs5)/(TestVs5$sd),2)

MDWVs5<-rbind(fvcTest)
rownames(MDWVs5)->MDWVs5$patientid

MDWVs5<-mutate(MDWVs5,NAVs5Imp = !is.na(MDWVs5$Vs5))
MDWVs5<- mutate(MDWVs5, Vs5= ifelse(is.na(Vs5), Vs5pred, Vs5))
MDWVs5$ErrorDev5<-round(((MDWVs5$Vs5pred - MDWVs5$Vs5)/(MDWVs5$Vs5)*100),2)

MDWVs5<-plyr::join(MDWVs4,MDWVs5, by="patientid")
MDWVs5[11]=NULL
MDWVs5[11]=NULL
MDWVs5[11]=NULL
MDWVs5<-data.frame(MDWVs5[ ,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs3pred","Vs4pred","Vs5pred",
                                "NAVs3Imp","NAVs4Imp","NAVs5Imp","ErrorDev3","ErrorDev4","ErrorDev5")])

MDWVs5T<-data.frame(MDWVs5[ ,c("NAVs3Imp","NAVs4Imp","NAVs5Imp")])

MDWVs5T <-MDWVs5T %>% mutate_if(is.character, as.logical)
summary(MDWVs5T)

MDWVs5<-data.frame(MDWVs5[ ,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs3pred","Vs4pred","Vs5pred","ErrorDev3","ErrorDev4","ErrorDev5")])
MDWVs5<-cbind(MDWVs5,MDWVs5T)
summary(MDWVs5)

MDWVSite<-data.frame(MDW[ ,c("patientid","site")])
MDWVs5<-plyr::join(MDWVSite,MDWVs5, by="patientid")

# = 16 start imputation visit Vs6 ==============================================
# = 16.1 data preparation and selection Vs6 ====================================

MDWVs6<-data.frame(MDW[ ,c("patientid","Vs6")])
MDWVs6<-plyr::join(MDWVs5,MDWVs6, by="patientid")
MDWVs6<-data.frame(MDWVs6[ ,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs6")])
sapply(MDWVs6, function(x) sum(is.na(x)))

rownames(MDWVs6)<-MDWVs6$patientid
MDWVs6$patientid=NULL

SIM_Vs6<-data.frame(SIM[ ,c("Vs1","Vs3","Vs4","Vs5","Vs6","Simulated")])
SIM_Vs6<-SIM_Vs6%>%filter(Simulated=="SimX100")
SIM_Vs6$Simulated=NULL

# = 17 split the dataset into training and validation ==========================

fvcTrain<-SIM_Vs6
fvcTest <-rbind(MDWVs6)

# = 18 Computing forest model (using ranger) classifier FVC ====================
cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs6~., 
                data=fvcTrain, 
                method="rf",
                tuneGrid = tunegrid,
                trControl=control)

#Shutdown cluster
stopCluster(cl)

print(rf_fit)
plot(rf_fit)
fvcImp <- varImp(rf_fit, scale = FALSE)
plot(fvcImp, top = 3)

FVCTrain_pred<- data.frame(rf_fit %>% predict(fvcTrain))
fvcTrain<-data.frame(fvcTrain,FVCTrain_pred$rf_fit.....predict.fvcTrain.)
names(fvcTrain)[6]<-paste("Vs6pred")

FVCTest_pred <- data.frame(rf_fit %>% predict(fvcTest))
fvcTest<-data.frame(fvcTest,FVCTest_pred$rf_fit.....predict.fvcTest.)
names(fvcTest)[6]<-paste("Vs6pred")



# = 19 Errors in the predictions and graphs Vs6 ================================

TrainVs6<-fvcTrain %>% dplyr:: summarise(n = n(),max=max(Vs6),min=min(Vs6),sd= sd(Vs6))
TestVs6<-fvcTest %>% dplyr:: summarise(n = n(),max=max(Vs6, na.rm=TRUE),min=min(Vs6, na.rm=TRUE),sd=sd(Vs6, na.rm=TRUE))

FVCTrain.RMSEVs6<-round(rmse(fvcTrain$Vs6pred,fvcTrain$Vs6)/(TrainVs6$sd),2)

fvcTestTemp<-na.omit(fvcTest)
FVCTest.RMSEVs6<-round(rmse(fvcTestTemp$Vs6pred ,fvcTestTemp$Vs6)/(TestVs6$sd),2)

MDWVs6<-data.frame(fvcTest)
rownames(MDWVs6)->MDWVs6$patientid

MDWVs6<-mutate(MDWVs6,NAVs6Imp = !is.na(MDWVs6$Vs6))
MDWVs6<-mutate(MDWVs6, Vs6= ifelse(is.na(Vs6), Vs6pred, Vs6))
MDWVs6$ErrorDev6<-round(((MDWVs6$Vs6pred - MDWVs6$Vs6)/(MDWVs6$Vs6)*100),2)

MDWVs6<-plyr::join(MDWVs5,MDWVs6, by="patientid")
MDWVs6[16]=NULL
MDWVs6[16]=NULL
MDWVs6[16]=NULL
MDWVs6[16]=NULL

MDWVs6<-data.frame(MDWVs6[ ,c("patientid","site","Vs1","Vs3","Vs4","Vs5","Vs6",
                               "Vs3pred","Vs4pred","Vs5pred","Vs6pred",
                               "NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp",
                               "ErrorDev3","ErrorDev4","ErrorDev5","ErrorDev6")])
summary(MDWVs6)


# = 20 start imputation visit Vs7 ==============================================
# = 20.1 data preparation and selection Vs7 ====================================

MDWVs7<-data.frame(MDW[ ,c("patientid","Vs7")])
MDWVs7<-plyr::join(MDWVs6,MDWVs7, by="patientid")
MDWVs7<-data.frame(MDWVs7[ ,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])
sapply(MDWVs7, function(x) sum(is.na(x)))

rownames(MDWVs7)<-MDWVs7$patientid
MDWVs7$patientid=NULL

SIM_Vs7<-data.frame(SIM[ ,c("Vs1","Vs3","Vs4","Vs5","Vs6","Vs7","Simulated")])
SIM_Vs7<-SIM_Vs7%>%filter(Simulated=="SimX100")
SIM_Vs7$Simulated=NULL

# = 21 split the dataset into training and validation ==========================

fvcTrain<-SIM_Vs7
fvcTest<-data.frame(MDWVs7)

# = 22 Computing forest model (using ranger) classifier FVC ====================

cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs7~., 
                data=fvcTrain, 
                method="rf",
                tuneGrid = tunegrid,
                trControl=control)

#Shutdown cluster
stopCluster(cl)

print(rf_fit)
plot(rf_fit)
fvcImp <- varImp(rf_fit, scale = FALSE)
plot(fvcImp, top = 3)

FVCTrain_pred<- data.frame(rf_fit %>% predict(fvcTrain))
fvcTrain<-data.frame(fvcTrain,FVCTrain_pred$rf_fit.....predict.fvcTrain.)
names(fvcTrain)[7]<-paste("Vs7pred")

FVCTest_pred <- data.frame(rf_fit %>% predict(fvcTest))
fvcTest <-data.frame(fvcTest,FVCTest_pred$rf_fit.....predict.fvcTest.)
names(fvcTest)[7]<-paste("Vs7pred")

# = 23 Errors in the predictions and graphs Vs7 ================================

TrainVs7<-fvcTrain %>% dplyr:: summarise(n = n(),max=max(Vs7),min=min(Vs7),sd= sd(Vs7, na.rm=TRUE))

TestVs7<-fvcTest %>%dplyr::summarise(n = n(),max=max(Vs7, na.rm=TRUE),min=min(Vs7, na.rm=TRUE),sd= sd(Vs7, na.rm=TRUE))

FVCTrain.RMSEVs7<-round(rmse(fvcTrain$Vs7pred,fvcTrain$Vs7)/(TrainVs7$sd),2)

fvcTestTemp<-na.omit(fvcTest)
FVCTest.RMSEVs7<-round(rmse(fvcTestTemp$Vs7pred ,fvcTestTemp$Vs7)/(TestVs7$sd),2)

MDWVs7<-rbind(fvcTest)
rownames(MDWVs7)->MDWVs7$patientid

MDWVs7<-mutate(MDWVs7,NAVs7Imp = !is.na(MDWVs7$Vs7))
MDWVs7<- mutate(MDWVs7, Vs7= ifelse(is.na(Vs7), Vs7pred, Vs7))
MDWVs7$ErrorDev7<-round(((MDWVs7$Vs7pred - MDWVs7$Vs7)/(MDWVs7$Vs7)*100),2)

MDWVs7<-plyr::join(MDWVs6,MDWVs7, by="patientid")
MDWVs7[20]=NULL
MDWVs7[20]=NULL
MDWVs7[20]=NULL
MDWVs7[20]=NULL
MDWVs7[20]=NULL

MDWVs7<-data.frame(MDWVs7[ ,c("patientid","site","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7",
                                "Vs3pred","Vs4pred","Vs5pred","Vs6pred","Vs7pred",
                                "NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp",
                                "ErrorDev3","ErrorDev4","ErrorDev5","ErrorDev6","ErrorDev7")])
summary(MDWVs7)

# = 24 Data aggregations =======================================================

MDW<-data.frame(MDW[ ,c("patientid","Imp","birthyear","YearDOD","d_v1date","VisitDate_F",
                         "ageV1","ageDOD","ageNow","DeadFromV1","Prog","gender","smoking","comorbidcount","v1weight",
                         "Status","DiedDuringStudy","Diedin5Years")])

MDW<-plyr::join(MDW,MDWVs7,by="patientid")

ErrorsFVC<-data.frame(c(FVCTest.RMSEVs3,FVCTest.RMSEVs4,FVCTest.RMSEVs5,FVCTest.RMSEVs6,FVCTest.RMSEVs7))

names(ErrorsFVC)[1]<-paste("NRMSE")
ErrorsFVC$Visits<-(c("Vs3","Vs4","Vs5","Vs6","Vs7"))
ErrorsFVC$Dataset<-(c("Test","Test","Test","Test","Test"))

f1 <- ggplot(ErrorsFVC, aes(x=Visits, y=NRMSE,fill=Dataset))+geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("NRMSE Random-Forest (Repeated Cross-Validation - Imputation)")
f1<-f1+ theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
f1

WDEV3<-data.frame(MDW[ ,c("patientid","site","NAVs3Imp","ErrorDev3")])
WDEV3<-WDEV3%>%mutate(visit_N="Vs3")
WDEV4<-data.frame(MDW[ ,c("patientid","site","NAVs4Imp","ErrorDev4")])
WDEV4<-WDEV4%>%mutate(visit_N="Vs4")
WDEV5<-data.frame(MDW[ ,c("patientid","site","NAVs5Imp","ErrorDev5")])
WDEV5<-WDEV5%>%mutate(visit_N="Vs5")
WDEV6<-data.frame(MDW[ ,c("patientid","site","NAVs6Imp","ErrorDev6")])
WDEV6<-WDEV6%>%mutate(visit_N="Vs6")
WDEV7<-data.frame(MDW[ ,c("patientid","site","NAVs7Imp","ErrorDev7")])
WDEV7<-WDEV7%>%mutate(visit_N="Vs7")
names(WDEV3)[3]<-paste ("NAVImp")
names(WDEV4)[3]<-paste ("NAVImp")
names(WDEV5)[3]<-paste("NAVImp")
names(WDEV6)[3]<-paste("NAVImp")
names(WDEV7)[3]<-paste("NAVImp")
names(WDEV3)[4]<-paste ("ErrorDev")
names(WDEV4)[4]<-paste ("ErrorDev")
names(WDEV5)[4]<-paste("ErrorDev")
names(WDEV6)[4]<-paste("ErrorDev")
names(WDEV7)[4]<-paste("ErrorDev")

ErroDevCom<-rbind(WDEV3,WDEV4,WDEV5,WDEV6,WDEV7)
ErroDevCom<-ErroDevCom%>%filter(NAVImp=="TRUE")

f2<-ggplot(ErroDevCom, aes(x=visit_N, y=ErrorDev, color=as.factor(site))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="Error Dev % ", x="Visit and Site", title="Random Forest Imputation - Error Dev in %",caption = "RiskStats") +
  scale_color_manual(values=c("Darkorange","Darkgreen"))+ 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
f2<-f2+ theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))
f2<- f2+ ylim(-100, 120)
f2         

dev.off()
dpi=100
tiff("Random Forest imputations Errors 070621.tiff",  res=dpi, height=10*dpi, width=25*dpi)
grid.arrange(f1,f2,ncol=2)
dev.off()

# = 25 imputation effect =======================================================


MDataH<- read.csv("MData_M4C_2502oX5v_130521.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataH)
dim(MDataH)
summary(MDataH)
MDataH$X=NULL
MDataH<-MDataH%>%mutate(Imputation="Imputed_Natural")
MDataH<-data.frame(MDataH[ ,c("patientid","ppFVC","visit_Idx","Imputation")])

MDWV<- read.csv("MDataV_CensMay21_230621_2490oX33v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))

MDWV<- data.frame(MDWV[,c("patientid","ppFVC","visit_Idx")])

MDWVx<-MDWV%>%mutate(Imputation="No_Imputed")
MDWVx<-(MDWVx%>%mutate(visit_Idx=case_when(visit_Idx==1 ~  1, visit_Idx==2 ~ 2, visit_Idx==3 ~ 3,
                                           visit_Idx==4 ~ 5, visit_Idx==5 ~ 10, visit_Idx==6 ~ 15)))

MDW1<- data.frame(MDW[,c("patientid","Vs1")])
MDW1<-MDW1%>%mutate(visit_Idx=1)
names(MDW1)[2]<-paste("ppFVC")

MDW3<- data.frame(MDW[,c("patientid","Vs3")])
MDW3<-MDW3%>%mutate(visit_Idx=2)
names(MDW3)[2]<-paste("ppFVC")

MDW4<- data.frame(MDW[,c("patientid","Vs4")])
MDW4<-MDW4%>%mutate(visit_Idx=3)
names(MDW4)[2]<-paste("ppFVC")

MDW5<- data.frame(MDW[,c("patientid","Vs5")])
MDW5<-MDW5%>%mutate(visit_Idx=5)
names(MDW5)[2]<-paste("ppFVC")

MDW6<- data.frame(MDW[,c("patientid","Vs6")])
MDW6<-MDW6%>%mutate(visit_Idx=10)
names(MDW6)[2]<-paste("ppFVC")

MDW7<- data.frame(MDW[,c("patientid","Vs7")])
MDW7<-MDW7%>%mutate(visit_Idx=15)
names(MDW7)[2]<-paste("ppFVC")

MDWImp<-rbind(MDW1,MDW3,MDW4,MDW5, MDW6, MDW7)
MDWImp<-MDWImp%>%arrange(patientid)

MDWImp<-MDWImp%>%mutate(Imputation="Impute Syntetic")

MData<-rbind(MDWVx,MDWImp,MDataH)


FImp<-ggplot(MData, aes(x=visit_Idx, y=ppFVC, color=as.factor(Imputation))) + 
  geom_jitter(position=position_jitter(0.1), cex=3, alpha = 0.1)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="ppFVC post-imputation",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name="Imputed",values=c("#00798c","#edae49","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FImp<-FImp + scale_x_discrete(name ="Time (Days)", limits=c("00","70-177","154-350","","322-581","","","","","564-823","","","","","1017-1251"))
FImp<-FImp+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), size= 0.3)   
#FImp<-FImp+ theme(legend.position = "none")  
FImp


MWSym1<- data.frame(SIM[,c("patientid","Vs1","Simulated")])
MWSym1<-MWSym1%>%mutate(visit_Idx=1)
names(MWSym1)[2]<-paste("ppFVC")

MWSym3<- data.frame(SIM[,c("patientid","Vs3","Simulated")])
MWSym3<-MWSym3%>%mutate(visit_Idx=2)
names(MWSym3)[2]<-paste("ppFVC")

MWSym4<- data.frame(SIM[,c("patientid","Vs4","Simulated")])
MWSym4<-MWSym4%>%mutate(visit_Idx=3)
names(MWSym4)[2]<-paste("ppFVC")

MWSym5<- data.frame(SIM[,c("patientid","Vs5","Simulated")])
MWSym5<-MWSym5%>%mutate(visit_Idx=5)
names(MWSym5)[2]<-paste("ppFVC")

MWSym6<- data.frame(SIM[,c("patientid","Vs6","Simulated")])
MWSym6<-MWSym6%>%mutate(visit_Idx=10)
names(MWSym6)[2]<-paste("ppFVC")

MWSym7<- data.frame(SIM[,c("patientid","Vs7","Simulated")])
MWSym7<-MWSym7%>%mutate(visit_Idx=15)
names(MWSym7)[2]<-paste("ppFVC")

MWSym<-rbind(MWSym1,MWSym3,MWSym4,MWSym5, MWSym6, MWSym7)
MWSym100<-MWSym%>%filter(Simulated=="SimX100")
MWSym100x<-data.frame(MWSym100[ ,c("patientid","ppFVC","visit_Idx")])
MWSym100x<-MWSym100x%>%mutate(Imputation="Syntetic")

G2A<-ggplot(MWSym100, aes(x =visit_Idx, y = ppFVC, group=patientid,color=(Simulated))) +
  ggtitle("ppFVC")+ geom_line(size = 0.5)+ geom_point(size = 0.9) + theme_bw()+
  labs(y="ppFVC", title="Synthetic database",caption = "Cohort Profile") +
  theme(plot.title = element_text(size=25, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333", size=12, angle=0))+
  scale_y_continuous(breaks=seq(0,180,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)
G2A<-G2A + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))  
G2A<-G2A +scale_color_manual(name="Site",values=c("#2e4057","#d1495b"))
G2A<-G2A + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
G2A


MDWV<- read.csv("MDataV_CensMay21_230621_2490oX33v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
MDWV<-MDWV%>%group_by(patientid)%>%filter(Imp==0)
MDWV<- data.frame(MDWV[,c("patientid","ppFVC","visit_Idx")])

MDWVx<-MDWV%>%mutate(Imputation="No_Imputed")
MDWVx<-(MDWVx%>%mutate(visit_Idx=case_when(visit_Idx==1 ~  1, visit_Idx==2 ~ 2, visit_Idx==3 ~ 3,
                                           visit_Idx==4 ~ 5, visit_Idx==5 ~ 10, visit_Idx==6 ~ 15)))

MDataX<-rbind(MDWImp,MWSym100x,MDWVx)


FImpX<-ggplot(MDataX, aes(x=visit_Idx, y=ppFVC, color=as.factor(Imputation))) + 
  geom_jitter(position=position_jitter(0.1), cex=3, alpha = 0.5)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="ppFVC post-imputation",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name="Imputed",values=c("#00798c","#edae49","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FImpX<-FImpX + scale_x_discrete(name ="Time (Days)", limits=c("00","70-177","154-350","","322-581","","","","","564-823","","","","","1017-1251"))
FImpX<-FImpX + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), size= 0.3)   
#FImp<-FImp+ theme(legend.position = "none")  
FImpX


dev.off()
dpi=100
tiff("Random Forest Syntectic Imputations 220621.tiff",  res=dpi, height=6*dpi, width=22*dpi)
grid.arrange(FImp,FImpX,G2A,ncol=3)
dev.off()

write.csv(ErrorsFVC,file ="MSynImp_Errors_220621_15ovX3v.csv",na="")
write.csv(MDW,file ="MDataImp_Syntetic_220621_415oX40v.csv",na="")