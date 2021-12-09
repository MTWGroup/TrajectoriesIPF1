#Project: IPF new project
#
# Purpose: Profile project: vs3 imputations Linear vs Knn and RF
# Version:1
# Date: 22/04/2020
# Author: HPF
#
# Input: 
# Main data base: MDataH_Nottingham_120421_242oX27v 
#       
#
# Output: imputations visits 1 and 3
#         linear and poly models
#         Detection of error in data collection
#         Data Preparation and Preprocessing
#         Make imputations on the test data
#         Method to test : LOCF, Linear, 10% decline, kknn, and RF
#         Compute the prediction error RMSE and deviation
#         Graphs RMSE and errorDev
#         Tables
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
library(DMwR)
library(FNN)
library(missForest)
library(tree)
library(randomForest)
library(doSNOW)

# = 3 open database ============================================================

MWDNotts<- read.csv("MDataH_Nottingham_120421_242oX27v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWDNotts)
dim(MWDNotts)
MWDNotts$X=NULL


MWDBromp<- read.csv("MDataH_Brompton_120421_175oX31v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWDBromp)
dim(MWDBromp)
MWDBromp$X=NULL

# = 4 data prepatation and selection ===========================================

MWDNotts_fvcS<-MWDNotts%>%filter(Imp==0)
MWDNotts_fvcS <- data.frame(MWDNotts_fvcS[,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])

MWDBromp_fvcS<-MWDBromp%>%filter(Imp==0)
MWDBromp_fvcS <- data.frame(MWDBromp_fvcS[,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])

MWD_fvcS<-rbind(MWDNotts_fvcS,MWDBromp_fvcS)

# = 6 Simulartion of missing data ==============================================

set.seed(123)
ind <- sample(2, nrow(MWDNotts_fvcS), replace = TRUE, prob = c(0.35, 0.65))
fvcTrainNotts <- MWDNotts_fvcS[ind==1,]
fvcTestNotts <- MWDNotts_fvcS[ind==2,]

set.seed(123)
ind <- sample(2, nrow(MWDBromp_fvcS), replace = TRUE, prob = c(0.35, 0.65))
fvcTrainBromp <- MWDBromp_fvcS[ind==1,]
fvcTestBromp <- MWDBromp_fvcS[ind==2,]

fvcTrain<-rbind(fvcTrainNotts,fvcTrainBromp)

createNAs <- function (x, pctNA = 0.94) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}

MWDNotts_fvcSNA<-data.frame(fvcTestNotts[ ,c("Vs7")])
set.seed(2000)
MWDNotts_fvcSNA<-data.frame(createNAs(MWDNotts_fvcSNA))
names(MWDNotts_fvcSNA)[1]<-paste("Vs7na")  
fvcL1NA_Notts<-cbind(fvcTestNotts,MWDNotts_fvcSNA)


MWDBromp_fvcSNA<-data.frame(fvcTestBromp[ ,c("Vs7")])
set.seed(2000)
MWDBromp_fvcSNA<-data.frame(createNAs(MWDBromp_fvcSNA))
names(MWDBromp_fvcSNA)[1]<-paste("Vs7na")  
fvcL1NA_Bromp<-cbind(fvcTestBromp,MWDBromp_fvcSNA)


createNAs <- function (x, pctNA = 0.63) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}

fvcL2NA_Notts<-data.frame(fvcTestNotts[ ,c("Vs6")])
set.seed(1211)
fvcL2NA_Notts<-data.frame(createNAs(fvcL2NA_Notts))
names(fvcL2NA_Notts)[1]<-paste("Vs6na")  
fvcL1NA_Notts<-cbind(fvcL1NA_Notts,fvcL2NA_Notts)


fvcL2NA_Bromp<-data.frame(fvcTestBromp[ ,c("Vs6")])
set.seed(1211)
fvcL2NA_Bromp<-data.frame(createNAs(fvcL2NA_Bromp))
names(fvcL2NA_Bromp)[1]<-paste("Vs6na")  
fvcL1NA_Bromp<-cbind(fvcL1NA_Bromp,fvcL2NA_Bromp)


createNAs <- function (x, pctNA = 0.38) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}


fvcL3NA_Notts<-data.frame(fvcTestNotts[ ,c("Vs5")])
set.seed(2000)
fvcL3NA_Notts<-data.frame(createNAs(fvcL3NA_Notts))
names(fvcL3NA_Notts)[1]<-paste("Vs5na")  
fvcL1NA_Notts<-cbind(fvcL1NA_Notts,fvcL3NA_Notts)

fvcL3NA_Bromp<-data.frame(fvcTestBromp[ ,c("Vs5")])
set.seed(2000)
fvcL3NA_Bromp<-data.frame(createNAs(fvcL3NA_Bromp))
names(fvcL3NA_Bromp)[1]<-paste("Vs5na")  
fvcL1NA_Bromp<-cbind(fvcL1NA_Bromp,fvcL3NA_Bromp)


createNAs <- function (x, pctNA = 0.26) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}


fvcL4NA_Notts<-data.frame(fvcTestNotts[ ,c("Vs4")])
set.seed(2000)
fvcL4NA_Notts<-data.frame(createNAs(fvcL4NA_Notts))
names(fvcL4NA_Notts)[1]<-paste("Vs4na")  
fvcL1NA_Notts<-cbind(fvcL1NA_Notts,fvcL4NA_Notts)

fvcL3NA_Notts<-data.frame(fvcTestNotts[ ,c("Vs3")])
names(fvcL3NA_Notts)[1]<-paste("Vs3na")  
fvcL1NA_Notts<-cbind(fvcL1NA_Notts,fvcL3NA_Notts)


createNAs <- function (x, pctNA = 0.50) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}

MWD_fvcSNA_Notts <- data.frame(fvcL1NA_Notts[,c("patientid","Vs1","Vs3","Vs4na","Vs5na","Vs6na","Vs7na")])
names(MWD_fvcSNA_Notts)[4]<-paste("Vs4")
names(MWD_fvcSNA_Notts)[5]<-paste("Vs5")
names(MWD_fvcSNA_Notts)[6]<-paste("Vs6")
names(MWD_fvcSNA_Notts)[7]<-paste("Vs7")


fvcL3NA_Bromp<-data.frame(fvcTestBromp[ ,c("Vs3")])
set.seed(2000)
fvcL3NA_Bromp<-data.frame(createNAs(fvcL3NA_Bromp))
names(fvcL3NA_Bromp)[1]<-paste("Vs3na")  
fvcL1NA_Bromp<-cbind(fvcL1NA_Bromp,fvcL3NA_Bromp)

fvcL4NA_Bromp<-data.frame(fvcTestBromp[ ,c("Vs4")])
names(fvcL4NA_Bromp)[1]<-paste("Vs4na")  
fvcL1NA_Bromp<-cbind(fvcL1NA_Bromp,fvcL4NA_Bromp)


MWD_fvcSNA_Bromp <- data.frame(fvcL1NA_Bromp[,c("patientid","Vs1","Vs3na","Vs4na","Vs5na","Vs6na","Vs7na")])
names(MWD_fvcSNA_Bromp)[3]<-paste("Vs3")
names(MWD_fvcSNA_Bromp)[4]<-paste("Vs4")
names(MWD_fvcSNA_Bromp)[5]<-paste("Vs5")
names(MWD_fvcSNA_Bromp)[6]<-paste("Vs6")
names(MWD_fvcSNA_Bromp)[7]<-paste("Vs7")


MWD_fvcSNA<-rbind(MWD_fvcSNA_Notts,MWD_fvcSNA_Bromp)
MWD_fvcSNA<-rbind(fvcTrain,MWD_fvcSNA)

# = 7 imputation matrix ========================================================

MWD_fvcSNA<-mutate(MWD_fvcSNA,NAVs3Imp = !is.na(MWD_fvcSNA$Vs3))
MWD_fvcSNA<-(MWD_fvcSNA%>%mutate(NAVs3Imp=case_when(NAVs3Imp=="TRUE" ~ 0 ,NAVs3Imp=="FALSE" ~ 1 )))

MWD_fvcSNA<-mutate(MWD_fvcSNA,NAVs4Imp = !is.na(MWD_fvcSNA$Vs4))
MWD_fvcSNA<-(MWD_fvcSNA%>%mutate(NAVs4Imp=case_when(NAVs4Imp=="TRUE" ~ 0 ,NAVs4Imp=="FALSE" ~ 1 )))

MWD_fvcSNA<-mutate(MWD_fvcSNA,NAVs5Imp = !is.na(MWD_fvcSNA$Vs5))
MWD_fvcSNA<-(MWD_fvcSNA%>%mutate(NAVs5Imp=case_when(NAVs5Imp=="TRUE" ~ 0 ,NAVs5Imp=="FALSE" ~1 )))

MWD_fvcSNA<-mutate(MWD_fvcSNA,NAVs6Imp = !is.na(MWD_fvcSNA$Vs6))
MWD_fvcSNA<-(MWD_fvcSNA%>%mutate(NAVs6Imp=case_when(NAVs6Imp=="TRUE" ~ 0 ,NAVs6Imp=="FALSE" ~1 )))

MWD_fvcSNA<-mutate(MWD_fvcSNA,NAVs7Imp = !is.na(MWD_fvcSNA$Vs7))
MWD_fvcSNA<-(MWD_fvcSNA%>%mutate(NAVs7Imp=case_when(NAVs7Imp=="TRUE" ~ 0 ,NAVs7Imp=="FALSE" ~ 1 )))

MWD_fvcSNA$SUMimp<-rowSums(MWD_fvcSNA[8:12])

names(MWD_fvcSNA)[3]<-paste("Vs3Imp")
names(MWD_fvcSNA)[4]<-paste("Vs4Imp")
names(MWD_fvcSNA)[5]<-paste("Vs5Imp")
names(MWD_fvcSNA)[6]<-paste("Vs6Imp")
names(MWD_fvcSNA)[7]<-paste("Vs7Imp")

MWD_fvcSNA <- data.frame(MWD_fvcSNA[,c("patientid","Vs3Imp","Vs4Imp","Vs5Imp","Vs6Imp","Vs7Imp",
                                       "NAVs3Imp","NAVs4Imp","NAVs5Imp","NAVs6Imp","NAVs7Imp","SUMimp")])
MWD_fvcSim<-plyr::join(MWD_fvcS,MWD_fvcSNA, by="patientid")

MWDBromp<-data.frame(MWDBromp[ ,c("patientid","site")])
MWDNotts<-data.frame(MWDNotts[ ,c("patientid","site")])
MWDSite<-rbind(MWDBromp,MWDNotts)

MWD_fvcSim<-plyr::join(MWD_fvcSim,MWDSite, by="patientid")


# = 8 Save simulations =========================================================

write.csv(MWD_fvcSim, file="MData_PROFILE_Sims_260421_85oX13v.csv", na="")