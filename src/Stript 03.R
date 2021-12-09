#Project: IPF new project
#
# Purpose: Profile project: vs3 imputations RF
# Version:2
# Date: 26/04/2020
# Author: HPF
#
# Input: MData_PROFILE_Sims_260421_85oX13v
# Main data base: MData_PROFILE_Sims_260421_85oX13v
#       
#
# Output: imputations visits 1 and 7
#         RF imputations
#         RMSE 
#         % Differences
#         DataBase
#
# Dependencies:@ F:/S2 PROFILE FINAL TEST SC 210421
#
# Notes: use the new dataset, could be a problem with the dates.
#
# =    1 working space =========================================================
setwd("F:/")
list.files("F:/")
#if (! dir.exists("Project_IPF_SOM")) dir.create("Project_IPF_SOM")
setwd("F:/S2 PROFILE FINAL TEST SC 210421")
list.files("F:/S2 PROFILE FINAL TEST SC 210421")
getwd()

# = 2 packages need ============================================================

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("sva"))
#install.packages("") 

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)

#Formatting/documentation packages
library(dplyr)
library(tidyr)
library(plyr)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(Metrics)
library(stringr)

#  K-Nearest Neighbors Essentials
library(caret)
library(skimr)
library(RANN)
library(DMwR2)
library(FNN)
library(missForest)
library(tree)
library(randomForest)
library(doSNOW)

# = 3 open database ============================================================

MWD<- read.csv("MData_PROFILE_Sims_260421_85oX13v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD)
dim(MWD)
MWD$X=NULL

# = 4 Random Forest (RF) imputation Vs4 ========================================

MWD_Vs4<-data.frame(MWD[ ,c("patientid","Vs1","Vs3","Vs4Imp","NAVs4Imp")])
MWD_Vs4_C<-MWD_Vs4%>%filter(NAVs4Imp==0)
MWD_Vs4_I<-MWD_Vs4

fvcTrain<-data.frame(MWD_Vs4_C[ ,c("Vs1","Vs3","Vs4Imp")])
fvcTest <-data.frame(MWD_Vs4_I[ ,c("Vs1","Vs3","Vs4Imp")])

# = 5 Computing forest model classifier FVC ====================================

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs4Imp~., 
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

# = 6 RF error model stimation Vs4 =============================================

MWD_RFTestA<-data.frame(fvcTest)
MWD_RFTestA<- MWD_RFTestA%>%mutate(Vs4Imp= ifelse(is.na(Vs4Imp),Vs4pred, Vs4Imp))
MWD$patientid->MWD_RFTestA$patientid

MWD_RFTestA<-data.frame(MWD_RFTestA[ ,c("patientid","Vs4Imp","Vs4pred")])
MWD_RF<-MWD
MWD_RF$Vs4Imp=NULL
MWD_RF<-plyr::join(MWD_RF,MWD_RFTestA, by="patientid")
MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site")])

MWD_RF$ErrorPer4<-round(((MWD_RF$Vs4Imp - MWD_RF$Vs4)/(MWD_RF$Vs4)*100),2)
MWD_RFvs4<-MWD_RF%>%filter(ErrorPer4>0 | ErrorPer4<0)

TestVs4<-MWD_RF%>% dplyr:: summarise(n = n(),max=max(Vs4, na.rm=TRUE),min=min(Vs4,na.rm=TRUE),sd= sd(Vs4, na.rm=TRUE))
RMSEvs4RF<-(round(rmse(MWD_RFvs4$Vs4Imp,MWD_RFvs4$Vs4)/(TestVs4$sd),2))

# = 7 Random Forest (RF) imputation Vs3 ========================================

MWD_Vs3<-data.frame(MWD_RF[ ,c("patientid","Vs1","Vs3Imp","Vs4Imp","NAVs3Imp")])
MWD_Vs3_C<-MWD_Vs3%>%filter(NAVs3Imp==0)
MWD_Vs3_I<-MWD_Vs3

fvcTrain<-data.frame(MWD_Vs3_C[ ,c("Vs1","Vs3Imp","Vs4Imp")])
fvcTest <-data.frame(MWD_Vs3_I[ ,c("Vs1","Vs3Imp","Vs4Imp")])

# = 8 Computing forest model classifier FVC ====================================

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs3Imp~., 
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

# = 9 RF error model estimation Vs3 =============================================

MWD_RFTestA<-data.frame(fvcTest)
MWD_RFTestA<- MWD_RFTestA%>%mutate(Vs3Imp= ifelse(is.na(Vs3Imp),Vs3pred, Vs3Imp))
MWD$patientid->MWD_RFTestA$patientid

MWD_RFTestA<-data.frame(MWD_RFTestA[ ,c("patientid","Vs3Imp","Vs3pred")])
MWD_RF$Vs3Imp=NULL
MWD_RF<-plyr::join(MWD_RF,MWD_RFTestA, by="patientid")
MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site","ErrorPer4")])

MWD_RF$ErrorPer3<-round(((MWD_RF$Vs3Imp - MWD_RF$Vs3)/(MWD_RF$Vs3)*100),2)
MWD_RFvs3<-MWD_RF%>%filter(ErrorPer3>0 | ErrorPer3<0)

TestVs3<-MWD_RF%>% dplyr:: summarise(n = n(),max=max(Vs3, na.rm=TRUE),min=min(Vs3,na.rm=TRUE),sd= sd(Vs3, na.rm=TRUE))
RMSEvs3RF<-(round(rmse(MWD_RFvs3$Vs3Imp,MWD_RFvs3$Vs3)/(TestVs3$sd),2))

MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site",
                              "ErrorPer3","ErrorPer4")])


# = 11 Random Forest (RF) imputation Vs5 =======================================

MWD_Vs5<-data.frame(MWD_RF[ ,c("patientid","Vs1","Vs3Imp","Vs4Imp","Vs5Imp","NAVs5Imp")])
MWD_Vs5_C<-MWD_Vs5%>%filter(NAVs5Imp==0)
MWD_Vs5_I<-MWD_Vs5

fvcTrain<-data.frame(MWD_Vs5_C[ ,c("Vs1","Vs3Imp","Vs4Imp","Vs5Imp")])
fvcTest <-data.frame(MWD_Vs5_I[ ,c("Vs1","Vs3Imp","Vs4Imp","Vs5Imp")])


# = 12 Computing forest model classifier FVC ====================================

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs5Imp~., 
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

# = 13 RF error model estimation Vs5 ===========================================

MWD_RFTestA<-data.frame(fvcTest)
MWD_RFTestA<- MWD_RFTestA%>%mutate(Vs5Imp= ifelse(is.na(Vs5Imp),Vs5pred, Vs5Imp))
MWD$patientid->MWD_RFTestA$patientid

MWD_RFTestA<-data.frame(MWD_RFTestA[ ,c("patientid","Vs5Imp","Vs5pred")])
MWD_RF$Vs5Imp=NULL
MWD_RF<-plyr::join(MWD_RF,MWD_RFTestA, by="patientid")
MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site",
                              "ErrorPer3","ErrorPer4")])

MWD_RF$ErrorPer5<-round(((MWD_RF$Vs5Imp - MWD_RF$Vs5)/(MWD_RF$Vs5)*100),2)
MWD_RFvs5<-MWD_RF%>%filter(ErrorPer5>0 | ErrorPer5<0)

TestVs5<-MWD_RF%>% dplyr:: summarise(n = n(),max=max(Vs5, na.rm=TRUE),min=min(Vs5,na.rm=TRUE),sd= sd(Vs5, na.rm=TRUE))
RMSEvs5RF<-(round(rmse(MWD_RFvs5$Vs5Imp,MWD_RFvs5$Vs5)/(TestVs5$sd),2))

MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site",
                              "ErrorPer3","ErrorPer4","ErrorPer5")])


# = 14 Random Forest (RF) imputation Vs6 =======================================

MWD_Vs6<-data.frame(MWD_RF[ ,c("patientid","Vs1","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp","NAVs6Imp")])
MWD_Vs6_C<-MWD_Vs6%>%filter(NAVs6Imp==0)
MWD_Vs6_I<-MWD_Vs6

fvcTrain<-data.frame(MWD_Vs6_C[ ,c("Vs1","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp")])
fvcTest <-data.frame(MWD_Vs6_I[ ,c("Vs1","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp")])

# = 15 Computing forest model classifier FVC ====================================

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs6Imp~., 
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

# = 16 RF error model estimation Vs6 ===========================================

MWD_RFTestA<-data.frame(fvcTest)
MWD_RFTestA<- MWD_RFTestA%>%mutate(Vs6Imp= ifelse(is.na(Vs6Imp),Vs6pred, Vs6Imp))
MWD$patientid->MWD_RFTestA$patientid

MWD_RFTestA<-data.frame(MWD_RFTestA[ ,c("patientid","Vs6Imp","Vs6pred")])
MWD_RF$Vs6Imp=NULL
MWD_RF<-plyr::join(MWD_RF,MWD_RFTestA, by="patientid")
MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site",
                              "ErrorPer3","ErrorPer4","ErrorPer5")])


MWD_RF$ErrorPer6<-round(((MWD_RF$Vs6Imp - MWD_RF$Vs6)/(MWD_RF$Vs6)*100),2)
MWD_RFvs6<-MWD_RF%>%filter(ErrorPer6>0 | ErrorPer6<0)

TestVs6<-MWD_RF%>% dplyr:: summarise(n = n(),max=max(Vs6, na.rm=TRUE),min=min(Vs6,na.rm=TRUE),sd= sd(Vs6, na.rm=TRUE))
RMSEvs6RF<-(round(rmse(MWD_RFvs6$Vs6Imp,MWD_RFvs6$Vs6)/(TestVs6$sd),2))

# = 17 Random Forest (RF) imputation Vs7 =======================================

MWD_Vs7<-data.frame(MWD_RF[ ,c("patientid","Vs1","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp","Vs7Imp","NAVs7Imp")])
MWD_Vs7_C<-MWD_Vs7%>%filter(NAVs7Imp==0)
MWD_Vs7_I<-MWD_Vs7

fvcTrain<-data.frame(MWD_Vs7_C[ ,c("Vs1","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp","Vs7Imp")])
fvcTest <-data.frame(MWD_Vs7_I[ ,c("Vs1","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp","Vs7Imp")])

# = 15 Computing forest model classifier FVC ====================================

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3, 
                        search='grid')

tunegrid <- expand.grid(.mtry = (1:15)) 


set.seed(2000)
rf_fit <- train(Vs7Imp~., 
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
fvcTest<-data.frame(fvcTest,FVCTest_pred$rf_fit.....predict.fvcTest.)
names(fvcTest)[7]<-paste("Vs7pred")

# = 16 RF error model estimation Vs6 ===========================================

MWD_RFTestA<-data.frame(fvcTest)
MWD_RFTestA<- MWD_RFTestA%>%mutate(Vs7Imp= ifelse(is.na(Vs7Imp),Vs7pred, Vs7Imp))
MWD$patientid->MWD_RFTestA$patientid

MWD_RFTestA<-data.frame(MWD_RFTestA[ ,c("patientid","Vs7Imp","Vs7pred")])
MWD_RF$Vs7Imp=NULL
MWD_RF<-plyr::join(MWD_RF,MWD_RFTestA, by="patientid")

MWD_RF<-data.frame(MWD_RF[ ,c("patientid","Vs1", "Vs3", "Vs4","Vs5","Vs6","Vs7", "Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp",
                              "Vs7Imp","NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp","site",
                              "ErrorPer3","ErrorPer4","ErrorPer5","ErrorPer6")])



MWD_RF$ErrorPer7<-round(((MWD_RF$Vs7Imp - MWD_RF$Vs7)/(MWD_RF$Vs7)*100),2)
MWD_RFvs7<-MWD_RF%>%filter(ErrorPer7>0 | ErrorPer7<0)

TestVs7<-MWD_RF%>% dplyr:: summarise(n = n(),max=max(Vs7, na.rm=TRUE),min=min(Vs7,na.rm=TRUE),sd= sd(Vs7, na.rm=TRUE))
RMSEvs7RF<-(round(rmse(MWD_RFvs7$Vs7Imp,MWD_RFvs7$Vs7)/(TestVs7$sd),2))


# = 17 saving data =============================================================


write.csv(MWD_RF, file="MDataTest_PROFILE_RF_260421_83oX24v.csv", na="")

#= 18 Edit the databases for graphs ============================================
  
#General
MWD_RF1<- data.frame(MWD_RF[,c("patientid","Vs1")])
MWD_RF1<-MWD_RF1%>%mutate(visit_Idx=1)
names(MWD_RF1)[2]<-paste("ppFVC")
MWD_RF1<-MWD_RF1%>%mutate(Imp="FALSE")

MWD_RF2<-MWD_RF%>%filter(NAVs3Imp==0)
MWD_RF2<- data.frame(MWD_RF2[,c("patientid","Vs3")])
MWD_RF2<-MWD_RF2%>%mutate(visit_Idx=2)
names(MWD_RF2)[2]<-paste("ppFVC")
MWD_RF2<-MWD_RF2%>%mutate(Imp="FALSE")

MWD_RF3<-MWD_RF%>%filter(NAVs4Imp==0)
MWD_RF3<- data.frame(MWD_RF3[,c("patientid","Vs4")])
MWD_RF3<-MWD_RF3%>%mutate(visit_Idx=3)
names(MWD_RF3)[2]<-paste("ppFVC")
MWD_RF3<-MWD_RF3%>%mutate(Imp="FALSE")

MWD_RF4<-MWD_RF%>%filter(NAVs5Imp==0)
MWD_RF4<- data.frame(MWD_RF4[,c("patientid","Vs5")])
MWD_RF4<-MWD_RF4%>%mutate(visit_Idx=5)
names(MWD_RF4)[2]<-paste("ppFVC")
MWD_RF4<-MWD_RF4%>%mutate(Imp="FALSE")

MWD_RF5<-MWD_RF%>%filter(NAVs6Imp==0)
MWD_RF5<- data.frame(MWD_RF[,c("patientid","Vs6")])
MWD_RF5<-MWD_RF5%>%mutate(visit_Idx=10)
names(MWD_RF5)[2]<-paste("ppFVC")
MWD_RF5<-MWD_RF5%>%mutate(Imp="FALSE")

MWD_RF6<-MWD_RF%>%filter(NAVs7Imp==0)
MWD_RF6<- data.frame(MWD_RF[,c("patientid","Vs7")])
MWD_RF6<-MWD_RF6%>%mutate(visit_Idx=15)
names(MWD_RF6)[2]<-paste("ppFVC")
MWD_RF6<-MWD_RF6%>%mutate(Imp="FALSE")

MWD_RFI2<-MWD_RF%>%filter(NAVs3Imp==1)
MWD_RFI2<- data.frame(MWD_RFI2[,c("patientid","Vs3Imp")])
MWD_RFI2<-MWD_RFI2%>%mutate(visit_Idx=2)
names(MWD_RFI2)[2]<-paste("ppFVC")
MWD_RFI2<-MWD_RFI2%>%mutate(Imp="TRUE")

MWD_RFI3<-MWD_RF%>%filter(NAVs4Imp==1)
MWD_RFI3<- data.frame(MWD_RFI3[,c("patientid","Vs4Imp")])
MWD_RFI3<-MWD_RFI3%>%mutate(visit_Idx=3)
names(MWD_RFI3)[2]<-paste("ppFVC")
MWD_RFI3<-MWD_RFI3%>%mutate(Imp="TRUE")

MWD_RFI4<-MWD_RF%>%filter(NAVs5Imp==1)
MWD_RFI4<- data.frame(MWD_RFI4[,c("patientid","Vs5Imp")])
MWD_RFI4<-MWD_RFI4%>%mutate(visit_Idx=5)
names(MWD_RFI4)[2]<-paste("ppFVC")
MWD_RFI4<-MWD_RFI4%>%mutate(Imp="TRUE")

MWD_RFI5<-MWD_RF%>%filter(NAVs6Imp==1)
MWD_RFI5<- data.frame(MWD_RFI5[,c("patientid","Vs6Imp")])
MWD_RFI5<-MWD_RFI5%>%mutate(visit_Idx=10)
names(MWD_RFI5)[2]<-paste("ppFVC")
MWD_RFI5<-MWD_RFI5%>%mutate(Imp="TRUE")

MWD_RFI6<-MWD_RF%>%filter(NAVs7Imp==1)
MWD_RFI6<- data.frame(MWD_RFI6[,c("patientid","Vs7Imp")])
MWD_RFI6<-MWD_RFI6%>%mutate(visit_Idx=15)
names(MWD_RFI6)[2]<-paste("ppFVC")
MWD_RFI6<-MWD_RFI6%>%mutate(Imp="TRUE")


MWD_RF_ALL<-rbind(MWD_RF1,MWD_RF2,MWD_RF3,MWD_RF4,MWD_RF5,MWD_RF6,
                    MWD_RFI2,MWD_RFI3,MWD_RFI4,MWD_RFI5,MWD_RFI6)

# = 19 Graphs ==================================================================

G1A<-ggplot(MWD_RF_ALL, aes(x=visit_Idx, y=ppFVC, color=as.factor(Imp))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Imputation RF",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333", size=12, angle=0))+
  scale_color_manual(values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
G1A<- G1A + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
G1A<- G1A + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
G1A

ERRORVs3<-data.frame(c(RMSEvs3RF))
ERRORVs3$Imputation<-(c("RF"))
ERRORVs3$ML<-(c("YES"))
ERRORVs3$Visit<-(c("90"))
names(ERRORVs3)[1]<-paste("RMSE")

ERRORVs4<-data.frame(c(RMSEvs4RF))
ERRORVs4$Imputation<-(c("RF"))
ERRORVs4$ML<-(c("YES"))
ERRORVs4$Visit<-(c("180"))
names(ERRORVs4)[1]<-paste("RMSE")

ERRORVs5<-data.frame(c(RMSEvs5RF))
ERRORVs5$Imputation<-(c("RF"))
ERRORVs5$ML<-(c("YES"))
ERRORVs5$Visit<-(c("365"))
names(ERRORVs5)[1]<-paste("RMSE")

ERRORVs6<-data.frame(c(RMSEvs6RF))
ERRORVs6$Imputation<-(c("RF"))
ERRORVs6$ML<-(c("YES"))
ERRORVs6$Visit<-(c("730"))
names(ERRORVs6)[1]<-paste("RMSE")

ERRORVs7<-data.frame(c(RMSEvs7RF))
ERRORVs7$Imputation<-(c("RF"))
ERRORVs7$ML<-(c("YES"))
ERRORVs7$Visit<-(c("1095"))
names(ERRORVs7)[1]<-paste("RMSE")

ERRORS<-rbind(ERRORVs3,ERRORVs4,ERRORVs5,ERRORVs6,ERRORVs7)
str(ERRORS)
ERRORS$Visit<-sapply(ERRORS$Visit, as.numeric)

f1 <- ggplot(ERRORS, aes(x=Imputation, y=RMSE, fill=as.factor(Visit)))+
  geom_bar(stat="identity", position="dodge")+ theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y="NRMSD",x="Imputation" ,title="NRMSD - Error RF",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(values=c("#00798c","#d1495b"))+
  geom_hline(yintercept = 0.13 , color="Black", linetype="dashed")+
  geom_hline(yintercept = 0.34 , color="Black", linetype="dashed")+
  geom_hline(yintercept = 0.43, color="Black", linetype="dashed")+
  geom_hline(yintercept = 0.48 , color="purple", linetype="dashed")+
  geom_hline(yintercept = 0.60 , color="purple", linetype="dashed")
f1<-f1+ theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))            
f1<-f1 + scale_fill_manual(name="Days", values=c("#edae49","#d1495b","#66a182","#2e4057","#8d96a3"))
f1

# = 20 ERROR in four visits ====================================================

write.csv(ERRORS, file="PROFILE_testIMP_RF_260421.csv", na="")

MWD_RF_ALL<-MWD_RF_ALL%>%mutate(Method="RF")

write.csv(MWD_RF_ALL, file="PROFILE_RFV_260421.csv", na="")


dev.off()
dpi=100
tiff("Test RF imputations.tiff", res=dpi, height=12*dpi, width=30*dpi)
grid.arrange(G1A,f1,ncol=2)
dev.off()

