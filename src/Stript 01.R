#Project: IPF new project
#
# Purpose: Profile project: verical data base with 7 visits
# Version: Gisly imputs (Brompton)
# Date: 12/04/2021
# Author:  HPF
#
# Input: Main data base: PROFILE_compcsv_iain2020.csv
#
#
# Output: consolidation of the database add ppFVC only Brompton
#         patientid
#         Site : Brompton
#         visit_N
#         visit_Idx
#         TimeV1toVfinal
#         visit_F
#         d_gen
#         d_dob
#         DOD
#         ageV1
#         ageNow
#         age
#         Status
#         DeadFromV1
#         FVCna (true missing)
#         incppFVF
#Add: mising data from visits : 3&4
#
# Dependencies:@/F:PROFILE FINAL IMP SC 120421/Stage 1 Data Set
#
# Notes: use the new dataset, cuold be a problem with the dates.
#
# =    1 working space =========================================================
setwd("F:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
list.files("F:/")
setwd("F:/S1 PROFILE FINAL DataEDIT SC 120421")
list.files("F:/S1 PROFILE FINAL DataEDIT SC 120421")

# = 2 packages need ============================================================

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("sva"))
#install.packages("eeptools") 

#Plotting and color options packages
library(gplots)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(lattice)

#Formatting/documentation packages
library(plyr)
library(doBy)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyverse)
library(eeptools)
library(zoo)
library(data.table)

# = 3 open database ============================================================

MData <- read.csv("PROFILE_linked_comp_updMay21A.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MData)
dim(MData)

# = 4 Select need it data from MData ===========================================
MData<-data.frame(MData[ ,c("patientid","site","popipf","lfvisitN","visit_date","d_v1date" ,"age","smoking",
                            "birthyear","DOD","Cens_date","ppFVC","gender","v1weight","comorbidcount")])
MData<-MData%>%filter(popipf==1)%>%filter(site=="Brompton"|site=="Nottingham")
sapply(MData, function(x) sum(is.na(x)))
names(MData)[4]<-paste("visit_N")
names(MData)[3]<-paste("d_ipf")
MData<-data.frame(MData[ ,c("patientid","site","visit_N","visit_date","d_v1date","age","birthyear",
                            "DOD","Cens_date","ppFVC","gender","smoking","v1weight","comorbidcount","d_ipf")])

# = 5 Repeating rows of data.frame =============================================

U_MData<-data.frame(MData[ ,c("patientid","visit_N")])
U_MData<-(U_MData%>%filter(!duplicated(patientid))) # total N= 276 IPF patients
U_MData<-data.frame(U_MData[ ,c("patientid")])
names(U_MData)[1]<-paste("patientid")
U_MData<-join(U_MData,MData, by="patientid") # garanty of not repite patients

U_MData<-data.frame(U_MData[ ,c("patientid","visit_N")])
U_MData<-(U_MData%>%group_by(patientid)%>%dplyr::summarise(n=n())) #elimination of patients with only one visit
U_MData<-(U_MData%>%filter(n>1)) # new total 271 from notts

U_MData<-data.frame(U_MData[-2])
U_MData<-(U_MData %>%do( data.frame(patientid = rep(.$patientid, each = 7), stringsAsFactors = FALSE))) # prototype of the working frame
U_MData <-(U_MData %>% group_by(patientid) %>% dplyr::mutate(visit_N = 1:n())) # add of the visits indexes

U_MData1<-(U_MData%>%filter(visit_N==1))
U_MData2<-(U_MData%>%filter(visit_N>2))
U_MData<-rbind(U_MData1,U_MData2) 
U_MData<-U_MData%>%arrange(patientid) # final working framework, no including visit two

U_MData<-data.frame(U_MData%>%unite(patientid, c("patientid","visit_N")))
U_MData<-(U_MData%>%filter(!duplicated(patientid)))

# = 6 Sync up MData to get 7 visits ============================================

MData<-data.frame(MData%>%unite(patientid, c("patientid","visit_N")))
MData<-join(U_MData,MData, by="patientid")
MData<-(MData%>%filter(!duplicated(patientid)))

MData<-(MData %>% dplyr::rowwise() %>% dplyr::mutate(visit_N = strsplit(patientid, split="_")[[1]][2]))  
MData<-(MData %>% dplyr::rowwise() %>% dplyr::mutate(patientid = strsplit(patientid, split="_")[[1]][1]))

MData<-data.frame(MData %>% group_by(patientid) %>% fill(site,birthyear,DOD,Cens_date,gender,smoking,v1weight,d_v1date,comorbidcount,d_ipf) 
                  %>%fill(site,birthyear,DOD,Cens_date,gender,smoking,v1weight,d_v1date,comorbidcount,d_ipf, .direction = "up") %>%ungroup())

sapply(MData, function(x) sum(is.na(x)))

# = 9 elimination of the patient without baseline and visit 3 or 4 ppFVC =======

MData1<-data.frame(MData[ ,c("patientid","visit_N", "ppFVC")]) # elimination of missing visit 1
MData1<-(MData1%>%filter(visit_N==1))
sapply(MData1, function(x) sum(is.na(x)))
MData1<-na.omit(MData1) # new total 267 patients
MData1<-data.frame(MData1[ ,c("patientid")])
names(MData1)[1]<-paste("patientid") 
MData1<-join(MData1,MData, by="patientid")
MData1<-data.frame(MData1[ ,c("patientid","visit_N", "ppFVC")]) 

MData3<-(MData1%>%filter(visit_N==3)) # elimination of missing visit 3
sapply(MData3, function(x) sum(is.na(x)))
MData3<-na.omit(MData3) # new total 215 patients

MData4<-(MData1%>%filter(visit_N==4)) # elimination of missing visit 4
sapply(MData4, function(x) sum(is.na(x)))
MData4<-na.omit(MData4) # new total 207 patients

MData1<-(MData1%>%filter(visit_N==1))
MData134<-rbind(MData1,MData3,MData4)
MData134<-(MData134%>%group_by(patientid)%>%dplyr::summarise(n=n()))
MData134<-MData134%>%filter(n>1)
MData134<-data.frame(MData134[ ,c("patientid")]) 

MData<-join(MData134,MData, by="patientid")

# = 10 format the dates ========================================================
# = 10.1 calulate the dates of the visits ======================================

MData<-(MData%>%mutate(FromV1Theo=case_when(visit_N==1 ~  0, visit_N==3 ~ 90, visit_N==4 ~ 180,
                                            visit_N==5 ~ 365, visit_N==6 ~ 730, visit_N==7 ~ 1095)))

# = 10.2 format visit 1 ========================================================

DoVisit1<- data.frame(MData[,c("patientid","d_v1date","visit_date","visit_N","FromV1Theo")])

DoVisit1<-(DoVisit1 %>% dplyr::rowwise() %>% dplyr::mutate(Day = strsplit(d_v1date, split="-")[[1]][1]))
DoVisit1<-(DoVisit1 %>% dplyr::rowwise() %>% dplyr::mutate(Month = strsplit(d_v1date, split="-")[[1]][2]))
DoVisit1<-(DoVisit1 %>% dplyr::rowwise() %>% dplyr::mutate(Year = strsplit(d_v1date, split="-")[[1]][3]))            

DoVisit1<-(DoVisit1%>%unite(d_v1date, c("Day", "Month","Year")))
DoVisit1<-(DoVisit1%>%mutate(d_v1date=str_replace_all(d_v1date,"_","/")))

DoVisit1$d_v1date <- as.Date(DoVisit1$d_v1date, format = "%d/%b/%y")

DoVisit1<-(DoVisit1 %>% dplyr::rowwise() %>% dplyr::mutate(Day = strsplit(visit_date, split="-")[[1]][1]))
DoVisit1<-(DoVisit1 %>% dplyr::rowwise() %>% dplyr::mutate(Month = strsplit(visit_date, split="-")[[1]][2]))
DoVisit1<-(DoVisit1 %>% dplyr::rowwise() %>% dplyr::mutate(Year = strsplit(visit_date, split="-")[[1]][3]))            

DoVisit1<-(DoVisit1%>%unite(visit_date, c("Day", "Month","Year")))
DoVisit1<-(DoVisit1%>%mutate(visit_date=str_replace_all(visit_date,"_","/")))

DoVisit1$visit_date <- as.Date(DoVisit1$visit_date, format = "%d/%b/%y")
DoVisit1$visit_dateTheo<-(DoVisit1$d_v1date + DoVisit1$FromV1Theo)
str(DoVisit1)
DoVisit1<-(DoVisit1%>%mutate(visit_date=ifelse(is.na(visit_date),visit_dateTheo,visit_date)))
str(DoVisit1)
DoVisit1$visit_date<-as.Date(DoVisit1$visit_date)
str(DoVisit1)
DoVisit1$FromV1<- floor(age_calc(DoVisit1$d_v1date,DoVisit1$visit_date ,units = "days"))
DoVisit1<-DoVisit1%>%group_by(patientid)%>%dplyr::mutate(visit_Idx=1:n())%>%ungroup()
DoVisit1<-data.frame(DoVisit1[,c("patientid","visit_N","visit_Idx","d_v1date","visit_date","visit_dateTheo","FromV1","FromV1Theo")])

# = 10.3 format DOD ============================================================

DoDeath<-data.frame(MData[,c("patientid","DOD","birthyear","visit_N","Cens_date")])
str(DoDeath)

DoDeath  <- mutate(DoDeath, DOD = ifelse(is.na(DOD),Cens_date, DOD))

DoDeath<-(DoDeath %>% dplyr::rowwise() %>% dplyr::mutate(Day = strsplit(DOD, split="-")[[1]][1]))
DoDeath<-(DoDeath %>% dplyr::rowwise() %>% dplyr::mutate(Month  = strsplit(DOD, split="-")[[1]][2]))
DoDeath<-(DoDeath %>% dplyr::rowwise() %>% dplyr::mutate(Year = strsplit(DOD, split="-")[[1]][3]))            
DoDeath<-(DoDeath %>% dplyr::rowwise() %>% dplyr::mutate(YearDOD = strsplit(DOD, split="-")[[1]][3]))  

DoDeath$Year[DoDeath$Year==40]<-20
DoDeath$YearDOD[DoDeath$YearDOD==40]<-20
sapply(DoDeath, function(x) sum(is.na(x)))

DoDeath<-(DoDeath%>%arrange(patientid)%>%unite(DOD, c("Day","Month", "Year")))

DoDeath<-(DoDeath%>%mutate(DOD=str_replace_all(DOD,"_","/")))
DoDeath$DOD <- as.Date(DoDeath$DOD, format = "%d/%b/%y")
str(DoDeath)
DoDeath$YearDOD<-sapply(DoDeath$YearDOD, as.numeric)
DoDeath<-DoDeath%>%mutate(YearDOD=YearDOD+2000)
DoDeath$birthyear<-sapply(DoDeath$birthyear, as.numeric)
str(DoDeath)

#= 10.4 joining dates ========================================================= 
  
DoDeath<-(DoDeath%>%unite(patientid, c("patientid","visit_N")))
DoVisit1<-(DoVisit1%>%unite(patientid, c("patientid","visit_N")))
Dates<-join(DoDeath,DoVisit1, by="patientid")

Dates<-(Dates %>% dplyr::rowwise() %>% dplyr::mutate(visit_N = strsplit(patientid, split="_")[[1]][2]))  
Dates<-(Dates %>% dplyr::rowwise() %>% dplyr::mutate(patientid = strsplit(patientid, split="_")[[1]][1]))

Dates<-Dates%>%mutate(ageDOD=(YearDOD-birthyear))

Dates<-Dates%>%mutate(ageNow=(2021-birthyear))

Dates$Prog<- floor(age_calc(Dates$d_v1date , Dates$DOD ,units = "days"))
Dates<-(Dates%>%mutate(DeadFromV1=(Prog/365))%>%dplyr::mutate_at(vars(DeadFromV1), funs(round(., 2))))
Dates<-(Dates%>%unite(patientid, c("patientid","visit_N")))
Dates$d_dob<-(Dates$DOD - (Dates$ageDOD)*365)  # only for Brompton

# = 11 Join date and all the data ==============================================

MData<-data.frame(MData[ ,c("patientid","site","visit_N","ppFVC","gender","smoking","comorbidcount","v1weight")])       
str(MData)
MData<-(MData%>%unite(patientid, c("patientid","visit_N")))
MData<-join(MData, Dates, by= "patientid")
MData<-(MData%>%filter(!duplicated(patientid)))

MData<-(MData %>% dplyr::rowwise() %>% dplyr::mutate(visit_N = strsplit(patientid, split="_")[[1]][2]))  
MData<-(MData %>% dplyr::rowwise() %>% dplyr::mutate(patientid = strsplit(patientid, split="_")[[1]][1]))

# = 12 detection of the last ppFVC value and visit of a patient ================

MData1<-data.frame(MData[ ,c("patientid","visit_N", "ppFVC","visit_date")])
MData1<-na.omit(MData1) 

MData1<-(MData1%>% group_by(patientid) %>% top_n(1, visit_N)%>%ungroup())
MData<-join(MData,MData1, by="patientid")

names(MData)[24]<-paste("visit_F")
names(MData)[25]<-paste("ppFVC_F")
names(MData)[26]<-paste("VisitDate_F")

MData$FVCna<-is.na(MData$ppFVC)

MData<-MData%>%mutate(ageV1=ageDOD-DeadFromV1)
MData<-(MData%>%mutate(Status=case_when(ageNow>ageDOD ~ 1 ,ageNow==ageDOD ~0 )))
MData<-(MData%>%mutate(DiedDuringStudy=case_when(Prog>1094 ~ 0, Prog<1095 ~1 )))
MData<-(MData%>%mutate(Diedin5Years=case_when(Prog>1824 ~ 0, Prog<1825 ~1 )))


sapply(MData, function(x) sum(is.na(x)))
MData<-MData%>%dplyr::mutate(smoking=replace_na(smoking,"UNK"))%>%dplyr::mutate(comorbidcount=replace_na(comorbidcount,0))
sapply(MData, function(x) sum(is.na(x)))

# = 13 selection of patient with 7 visits and imputation per patient ===========

MData1<-(MData%>%filter(FVCna=="FALSE"))
MData1<-(MData1%>%group_by(patientid)%>%dplyr::summarise(n=n())%>%mutate(Imp=6-n))
MData<-join(MData,MData1, by="patientid")

# = 14 horizontal data base ====================================================

MDataH <- data.frame(MData [,c("patientid","ppFVC","visit_N")])
MDataH <- spread(MDataH, patientid, ppFVC)

MDataH<-(MDataH%>%dplyr::mutate(visit_N=str_replace(visit_N,"1","Vs1"))
         %>%dplyr::mutate(visit_N=str_replace(visit_N,"3","Vs3"))
         %>%dplyr::mutate(visit_N=str_replace(visit_N,"4","Vs4"))
         %>%dplyr::mutate(visit_N=str_replace(visit_N,"5","Vs5"))
         %>%dplyr::mutate(visit_N=str_replace(visit_N,"6","Vs6"))
         %>%dplyr::mutate(visit_N=str_replace(visit_N,"7","Vs7")))

rownames(MDataH)<-MDataH$visit_N
MDataH$visit_N=NULL
MDataH<-data.frame(t(MDataH))
str(MDataH)

MDataH$patientid<-rownames(MDataH)

MDataH<-data.frame(MDataH[ ,c("patientid","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])
MData1<-(MData%>%filter(visit_N==1))

MData1$ppFVC=NULL

MDataH<-join(MData1,MDataH,by="patientid")

# = 15 saving data =============================================================

write.csv(MData, file="MDataV_CensMay21_230621_2490oX33v.csv", na="")
write.csv(MDataH,file ="MDataH_CensMay21_230621_415oX32v.csv",na="")
