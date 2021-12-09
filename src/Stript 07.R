#Project: IPF new project
#
# Purpose: Anova PROFILE  
# Version: 1
# Date: 24/04/2021
# Author:  HPF
#
# Input: The Oneway Analysis of Variance (ANOVA)
#
#
# Output: Comparsion visit by visits 
#
#
# Dependencies:@/F:PROFILE FINAL STATS SC 1504021
#
# Notes: use the new dataset, cuold be a problem with the dates.
#
# =    1 working space =========================================================
setwd("F:/")
#if (! dir.exists("Project_IPF_SOM")) dir.create("Project_IPF_SOM")
setwd("F:/T.8 PROFILE STATS MCMC 24 vs MCMC")
list.files("F:/T.8 PROFILE STATS MCMC 24 vs MCMC")
getwd()

# = 2 packages need ============================================================
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("impute"))

#Analytical packages
library("limma")

#Formatting/documentation packages
library("tidyverse")
library("plyr")
library("doBy")
library("reshape2")
library("stringr")

#Graphic packages
library("ggplot2")
library("gridExtra")
library("RColorBrewer")

#Tables packeges
library("knitr")
library("kableExtra")
library("lsr")
library("psych")
options(width = 130)
options(knitr.table.format = "html") 

library(cowplot)
library(patchwork)
library(lattice)
library(grid)
library(ggpubr)
# = 3 open database ============================================================

MData<- read.csv("MDataH_CensMay21_230621_415oX32v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MData)
dim(MData)
summary(MData)
MData$X=NULL
MData<-MData%>%mutate(ImpType="NO_IMP")

MData_H<- read.csv("MDataV_CensMay21_230621_2490oX33v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MData_H)
dim(MData_H)
summary(MData_H)
MData_H$X=NULL
MData_H<-MData_H%>%mutate(ImpType="NO_IMP")
MData_H<-data.frame(MData_H[ ,c("patientid","visit_N","ppFVC","ImpType")]) 
MData_H<-MData_H%>%mutate(visit_Idx=case_when(visit_N==1~1,visit_N==3~2,visit_N==4~3,visit_N==5~5,visit_N==6~10,visit_N==7~15))
MData_H<-data.frame(MData_H[ ,c("patientid","visit_Idx","ppFVC","ImpType")])


MDataMCMC<- read.csv("MDataSoMV_3YAll_MCMC_415oX49v_110821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataMCMC)
dim(MDataMCMC)
summary(MDataMCMC)
MDataMCMC$X=NULL
MDataMCMC$DeadIn2Y=NULL
MDataMCMC<-MDataMCMC%>%mutate(ImpType="MCMC")

MDataMCMC_H<- read.csv("MDataSoMH_3YAll_MCMC_2490oX7v_110821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataMCMC_H)
dim(MDataMCMC_H)
MDataMCMC_H$X=NULL
MDataMCMC_H$DeadIn2Y=NULL
MDataMCMC_H<-MDataMCMC_H%>%mutate(ImpType="MCMC")
MDataMCMC_H<-data.frame(MDataMCMC_H[ ,c("patientid","visit_Idx","ppFVC","ImpType")]) 


MDataMCMC24<- read.csv("MDataSoMV_3Y_EXTR_415oX48v_050821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataMCMC24)
dim(MDataMCMC24)
summary(MDataMCMC24)
MDataMCMC24$X=NULL
MDataMCMC24<-MDataMCMC24%>%mutate(ImpType="MCMC24")

MDataMCMC24_H<- read.csv("MDataSoMH_3Y_EXTR_2490oX6v_050821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataMCMC24_H)
dim(MDataMCMC24_H)
MDataMCMC24_H$X=NULL
MDataMCMC24_H<-MDataMCMC24_H%>%mutate(ImpType="MCMC24")
MDataMCMC24_H<-data.frame(MDataMCMC24_H[ ,c("patientid","visit_Idx","ppFVC","ImpType")]) 

# = 4 Join date and all the data ===============================================

MDataALL_H<-rbind(MData_H,MDataMCMC_H,MDataMCMC24_H)
MDataALL_H<-MDataALL_H%>%mutate(ImpType=case_when(ImpType=="NO_IMP"~ "No Imputed Values", ImpType=="MCMC"~"RF Naive", 
                                                  ImpType=="MCMC24" ~ "RF Theoretic"))


ImpppFVC <- MDataALL_H%>% group_by(visit_Idx,ImpType) %>% dplyr:: summarise(n = n(),
                                                                          mean = mean(ppFVC,  na.rm=TRUE),
                                                                          median = median(ppFVC,  na.rm=TRUE),
                                                                          sd = sd(ppFVC,  na.rm=TRUE)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                                    CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                                    CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)


names(ImpppFVC)[4]<-paste("ppFVC")

ALLM<-ggplot(MDataALL_H, aes(x=visit_Idx, y=ppFVC, color=as.factor(ImpType)))+
             geom_jitter(position=position_jitter(0.1), cex=1, alpha= 0.3)+ theme_bw()+
             geom_errorbar(aes(ymin = ppFVC-sd, ymax = ppFVC+sd), data =ImpppFVC, width =1, size= 1, position = position_dodge(0.1))+
            #stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
             labs(y="ppFVC ", title=" ",caption = "Cohort PROFILE", tag = "Supp F1")+
             theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
                   axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
                   axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
            geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
            theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
            axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
            scale_color_manual(name=" " ,values=c("#00798c","#edae49","#d1495b"))+
            scale_y_continuous(breaks=seq(0,120,20))+
            scale_x_continuous(breaks=seq(0,12,4))+
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
            scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))+
            stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), size= 1)+
            geom_segment(aes(x=9.8,xend=15.2,y=135,yend=135))+
            geom_segment(aes(x=9.8,xend=9.8,y=135,yend=134))+
            geom_segment(aes(x=15.2,xend=15.2,y=135,yend=134))+
            annotate("text", label = "A,B***", x = 12.5, y = 137 , size = 5, colour = "darkred", fontface = "bold")+
            annotate("text", label = "A,C***", x = 12.5, y = 141 , size = 5, colour = "darkred", fontface = "bold")+
            geom_segment(aes(x=14.8,xend=15.2,y=128,yend=128))+
            geom_segment(aes(x=14.8,xend=14.8,y=128,yend=127))+
            geom_segment(aes(x=15.2,xend=15.2,y=128,yend=127))+ 
            annotate("text", label = "A,B,C***", x = 15, y = 130 , size = 5, colour = "darkred", fontface = "bold")+
            theme(legend.text=element_text(size=18)) 
            #theme(legend.position = "none")  
ALLM



textsup1 <- paste("Shown are mean (±SE) changes from baseline (dashed line) in lung function over", 
                  "the three year period. The data shown are for patients with available data without", 
                  "imputation (green) as well as for the results of two imputation methods used to account", 
                  "for missing data (naïve- yellow and theoretic -red). Dots represent the values -", 
                  "non imputed (green) and imputed (yellow and red) for each patient assessed at that ",
                  "time point. Different superscript indicates significant differences between non ",
                  "imputed values (A) and both imputations methods  naïve (B) or theoretic (C) (***p<0.005).", sep = " ")

textsup11 <- ggparagraph(text = textsup1 ,size = 20, color = "black")

dev.off()
dpi=300
tiff("Supplementary Figure 1.tiff", res=dpi, height=15*dpi, width=15*dpi)
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(6,6)))

define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(ALLM, vp=define_region(1:4,1:6))
print(textsup11, vp=define_region(6:5,1:6))
dev.off()

# = 5 Join the vertical data ===================================================

MDataMCMC<-data.frame(MDataMCMC[ ,c("patientid","ImpType","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])              
MDataMCMC24<-data.frame(MDataMCMC24[ ,c("patientid","ImpType","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])  
MData<-data.frame(MData[ ,c("patientid","ImpType","Vs1","Vs3","Vs4","Vs5","Vs6","Vs7")])  
#MData[is.na(MData)] <- 0

MData<-rbind(MData,MDataMCMC,MDataMCMC24)
MDataA<-data.frame(MData%>%unite(patientid, c("patientid","ImpType")))

#= 6 explore data ==============================================================

Table1<-describeBy(MData$Vs1,group = MData$ImpType, mat = TRUE) %>% #create dataframe
  select(Cluster=group1, Patients=n, ppFVC=mean, SD=sd, Median=median, Min=min, Max=max, 
         Skew=skew, Kurtosis=kurtosis, SEM=se) %>% 
  kable(align=c("lrrrrrrrr"), digits=2, row.names = FALSE,
        caption="ppFVC Baseline Descriptive Statistics Vs1") %>% 
  kable_styling(bootstrap_options=c("bordered", "responsive","striped"), full_width = FALSE)
Table1


Table3<-describeBy(MData$Vs3,group = MData$ImpType, mat = TRUE) %>% #create dataframe
  select(Cluster=group1, Patients=n, ppFVC=mean, SD=sd, Median=median, Min=min, Max=max, 
         Skew=skew, Kurtosis=kurtosis, SEM=se) %>% 
  kable(align=c("lrrrrrrrr"), digits=2, row.names = FALSE,
        caption="ppFVC Baseline Descriptive Statistics Vs3") %>% 
  kable_styling(bootstrap_options=c("bordered", "responsive","striped"), full_width = FALSE)
Table3


Table4<-describeBy(MData$Vs4,group = MData$ImpType, mat = TRUE) %>% #create dataframe
  select(Cluster=group1, Patients=n, ppFVC=mean, SD=sd, Median=median, Min=min, Max=max, 
         Skew=skew, Kurtosis=kurtosis, SEM=se) %>% 
  kable(align=c("lrrrrrrrr"), digits=2, row.names = FALSE,
        caption="ppFVC Baseline Descriptive Statistics Vs4") %>% 
  kable_styling(bootstrap_options=c("bordered", "responsive","striped"), full_width = FALSE)
Table4

Table5<-describeBy(MData$Vs5,group = MData$ImpType, mat = TRUE) %>% #create dataframe
  select(Cluster=group1, Patients=n, ppFVC=mean, SD=sd, Median=median, Min=min, Max=max, 
         Skew=skew, Kurtosis=kurtosis, SEM=se) %>% 
  kable(align=c("lrrrrrrrr"), digits=2, row.names = FALSE,
        caption="ppFVC Baseline Descriptive Statistics Vs5") %>% 
  kable_styling(bootstrap_options=c("bordered", "responsive","striped"), full_width = FALSE)
Table5

Table6<-describeBy(MData$Vs6,group = MData$ImpType, mat = TRUE) %>% #create dataframe
  select(Cluster=group1, Patients=n, ppFVC=mean, SD=sd, Median=median, Min=min, Max=max, 
         Skew=skew, Kurtosis=kurtosis, SEM=se) %>% 
  kable(align=c("lrrrrrrrr"), digits=2, row.names = FALSE,
        caption="ppFVC Baseline Descriptive Statistics Vs6") %>% 
  kable_styling(bootstrap_options=c("bordered", "responsive","striped"), full_width = FALSE)
Table6

Table7<-describeBy(MData$Vs7,group = MData$ImpType, mat = TRUE) %>% #create dataframe
  select(Cluster=group1, Patients=n, ppFVC=mean, SD=sd, Median=median, Min=min, Max=max, 
         Skew=skew, Kurtosis=kurtosis, SEM=se) %>% 
  kable(align=c("lrrrrrrrr"), digits=2, row.names = FALSE,
        caption="ppFVC Baseline Descriptive Statistics Vs7") %>% 
  kable_styling(bootstrap_options=c("bordered", "responsive","striped"), full_width = FALSE)
Table7


#= 5 analysis of the time ppFVC ================================================

MDataT<-t(MDataA[-1])
colnames(MDataT)<-c(MDataA$patientid)
MDataT<-data.frame(MDataT)
MDataT<-MDataT %>% mutate_if(is.character,as.numeric)

# = 6 eBayer model set up ======================================================

treatments_Cl<-c("NO_IMP","MCMC","MCMC24")

assay_names_Cl<-c(rep("NO_IMP",415),rep( "MCMC",415), rep("MCMC24",415))

design_Cl<-model.matrix(~0 + factor(assay_names_Cl, levels= treatments_Cl))
colnames(design_Cl)<-treatments_Cl

contrasts_Cl<-makeContrasts(NO_IMP-MCMC,NO_IMP-MCMC24,
                            MCMC-MCMC24,levels=treatments_Cl) 

fit_Cl=lmFit(MDataT, design_Cl)
fit1_Cl=contrasts.fit(fit_Cl, contrasts=contrasts_Cl)
fit1_Cl=eBayes(fit1_Cl)


# = 7 Contracts ================================================================

StatsppFVC <- MDataALL_H %>% group_by(visit_Idx,ImpType) %>% dplyr:: summarise(n = n(),
                                                                       mean = mean(ppFVC, na.rm=T),
                                                                       median = median(ppFVC, na.rm=T),
                                                                       sd = sd(ppFVC, na.rm=T)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                                 CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                                 CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

StatsppFVC<-data.frame(StatsppFVC)
StatsppFVC<-StatsppFVC%>%mutate(visit_Idx=case_when(visit_Idx=="1" ~ "Vs1", visit_Idx=="2" ~ "Vs3", 
                                                    visit_Idx=="3" ~ "Vs4", visit_Idx=="5" ~ "Vs5",
                                                    visit_Idx=="10" ~ "Vs6", visit_Idx=="15" ~ "Vs7"))


Cl1<-StatsppFVC%>%filter(ImpType=="NO_IMP")
names(Cl1)[1]<-paste("visit")
names(Cl1)[4]<-paste("ppFVC.NO_IMP")
names(Cl1)[7]<-paste("sem.NO_IMP")
Cl1<- data.frame(Cl1[,c("visit","ppFVC.NO_IMP","sem.NO_IMP")])             


Cl2<-StatsppFVC%>%filter(ImpType=="MCMC")
names(Cl2)[1]<-paste("visit")
names(Cl2)[4]<-paste("ppFVC.MCMC")
names(Cl2)[7]<-paste("sem.MCMC")
Cl2<- data.frame(Cl2[,c("visit","ppFVC.MCMC","sem.MCMC")])             


Cl3<-StatsppFVC%>%filter(ImpType=="MCMC24")
names(Cl3)[1]<-paste("visit")
names(Cl3)[4]<-paste("ppFVC.MCMC24")
names(Cl3)[7]<-paste("sem.MCMC24")
Cl3<- data.frame(Cl3[,c("visit","ppFVC.MCMC24","sem.MCMC24")])

Cl_list_A1<-topTable(fit1_Cl, number=nrow(MDataT))

Cl_list_B1G<-topTable(fit1_Cl, coef="NO_IMP - MCMC", number = Inf)
Cl_list_B1G<-Cl_list_B1G%>%mutate(adj.P.Val.Bon=(p.adjust(Cl_list_B1G$P.Value, method ="bonferroni", n = 12)))%>%
  mutate(contrast="NO_IMP-MCMC")%>%mutate(P.Value.Star=case_when(adj.P.Val.Bon <0.05 ~ "*", adj.P.Val.Bon >0.05 ~ "NS"))%>%
  mutate(visit=rownames(Cl_list_B1G))
Cl_list_B1G<-join(Cl2,Cl_list_B1G, by="visit")
Cl_list_B1G<-join(Cl1,Cl_list_B1G, by="visit")
Cl_list_B1G$AveExpr=NULL
Cl_list_B1G$t=NULL
Cl_list_B1G$B=NULL
Cl_list_B1G<-Cl_list_B1G%>%mutate(visit=case_when(visit=="Vs1"~ 0.0,visit=="Vs3"~ 90, visit=="Vs4"~ 180,  visit=="Vs5"~ 365, 
                                                  visit=="Vs6"~ 730, visit=="Vs7" ~ 1095))


Cl_list_B2G<-topTable(fit1_Cl, coef="NO_IMP - MCMC24", number = Inf)
Cl_list_B2G<-Cl_list_B2G%>%mutate(adj.P.Val.Bon=(p.adjust(Cl_list_B2G$P.Value, method ="bonferroni", n = 12)))%>%
  mutate(contrast="NO_IMP-MCMC24")%>%mutate(P.Value.Star=case_when(adj.P.Val.Bon <0.05 ~ "*", adj.P.Val.Bon >0.05 ~ "NS"))%>%
  mutate(visit=rownames(Cl_list_B2G))
Cl_list_B2G<-join(Cl3,Cl_list_B2G, by="visit")
Cl_list_B2G<-join(Cl1,Cl_list_B2G, by="visit")
Cl_list_B2G$AveExpr=NULL
Cl_list_B2G$t=NULL
Cl_list_B2G$B=NULL
Cl_list_B2G<-Cl_list_B2G%>%mutate(visit=case_when(visit=="Vs1"~ 0.0,visit=="Vs3"~ 90, visit=="Vs4"~ 180,  visit=="Vs5"~ 365, 
                                                  visit=="Vs6"~ 730, visit=="Vs7" ~ 1095))


Cl_list_B3G<-topTable(fit1_Cl, coef="MCMC - MCMC24", number = Inf)
Cl_list_B3G<-Cl_list_B3G%>%mutate(adj.P.Val.Bon=(p.adjust(Cl_list_B3G$P.Value, method ="bonferroni", n = 12)))%>%
  mutate(contrast="MCMCVsMCMC24")%>%mutate(P.Value.Star=case_when(adj.P.Val.Bon <0.05 ~ "*", adj.P.Val.Bon >0.05 ~ "NS"))%>%
  mutate(visit=rownames(Cl_list_B3G))
Cl_list_B3G<-join(Cl2,Cl_list_B3G, by="visit")
Cl_list_B3G<-join(Cl3,Cl_list_B3G, by="visit")
Cl_list_B3G$AveExpr=NULL
Cl_list_B3G$t=NULL
Cl_list_B3G$B=NULL
Cl_list_B3G<-Cl_list_B3G%>%mutate(visit=case_when(visit=="Vs1"~ 0.0,visit=="Vs3"~ 90, visit=="Vs4"~ 180,  visit=="Vs5"~ 365, 
                                                  visit=="Vs6"~ 730, visit=="Vs7" ~ 1095))



# = 12 Save the dataframes =====================================================

write.csv(Cl_list_A1, file="Analysis of Variance Imputation ppFVC 415 080921.csv",na="")
write.csv(Cl_list_B1G, file="post hoc comparison Imputation ppFVC 415 NO_IMP vs MCMC 080921.csv", na="")
write.csv(Cl_list_B2G, file="post hoc comparison Imputation ppFVC 415 NO_IMP vs MCMC24 080921.csv", na="")
write.csv(Cl_list_B3G, file="post hoc comparison Imputation ppFVC 415 MCMC vs MCMC24 080921.csv", na="")



dev.off()
dpi=100
tiff("Stats PROFILE 080921.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ALLM
dev.off()



MDataMCMC_H<- read.csv("MDataSoMH_3YAll_MCMC_2490oX7v_110821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataMCMC_H)
dim(MDataMCMC_H)
MDataMCMC_H$X=NULL
MDataMCMC_H$DeadIn2Y=NULL
MDataMCMC_Cl<-data.frame(MDataMCMC_H[ ,c("patientid","Cl4_1095")])
MDataMCMC_Cl<-(MDataMCMC_Cl%>%filter(!duplicated(patientid)))
MDataMCMC_H<-MDataMCMC_H%>%mutate(ImpType="MCMC")
MDataMCMC_H<-data.frame(MDataMCMC_H[ ,c("patientid","visit_Idx","ppFVC","ImpType","DiedDuringStudy")]) 

MDataMCMC24_H<- read.csv("MDataSoMH_3Y_EXTR_2490oX6v_050821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataMCMC24_H)
dim(MDataMCMC24_H)
MDataMCMC24_H$X=NULL
MDataMCMC24_H<-MDataMCMC24_H%>%mutate(ImpType="MCMC24")
MDataMCMC24_H<-data.frame(MDataMCMC24_H[ ,c("patientid","visit_Idx","ppFVC","ImpType","DiedDuringStudy")]) 

MDataALL_H<-rbind(MDataMCMC_H,MDataMCMC24_H)
MDataALL_H<-join(MDataALL_H,MDataMCMC_Cl,by="patientid")


# = 9.1 graphs Cluster 3 M4 ====================================================

ClM4<- MDataALL_H%>% group_by(Cl4_1095, visit_Idx,ImpType) %>% dplyr:: summarise(n = n(),
                                                                                          mean = mean(ppFVC),
                                                                                          median = median(ppFVC),
                                                                                          sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                                                    CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                                                    CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

Cl1M4<-MDataALL_H%>%filter(Cl4_1095==1)

FCl1M4D<-ggplot(Cl1M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(ImpType))) + 
  geom_jitter(position=position_jitter(0.1),alpha= 0.3, cex=2)+ 
  theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="",caption = "Cohort PROFILE", tag = "A") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#A4DE02"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#  annotate("text", label ="67/140(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
#  annotate("text", label ="73/140(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl1M4D<-FCl1M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4D<-FCl1M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl1M4D<-FCl1M4D+ theme(legend.position = "none")  
FCl1M4D

FCl1M4D<-FCl1M4D + annotate("text", label ="  0.0%", x = 1, y = 145 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label =" -1.8%", x = 2, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label =" -3.9%", x = 3, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label ="-12.8%",x = 5, y = 145 , size = 4, colour  ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label ="-20.8%", x = 10, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label ="-23.8%", x = 15, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label ="NS", x = 2, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label ="NS", x = 3, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label =" 4.3e-05", x = 5, y = 150 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label =" 2.8e-12", x = 10, y = 150 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D<-FCl1M4D + annotate("text", label =" 1.2e-15", x = 15, y = 150 , size = 4, colour ="#993333", fontface = "bold")
FCl1M4D

Cl2M4<-MDataALL_H%>%filter(Cl4_1095==2)

FCl2M4D<-ggplot(Cl2M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(ImpType))) + 
  geom_jitter(position=position_jitter(0.1),alpha= 0.3, cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title=" ",caption = "Cohort PROFILE", tag = "B") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#d1495b","#F6BDC0"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#  annotate("text", label ="67/140(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
#  annotate("text", label ="73/140(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl2M4D<-FCl2M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4D<-FCl2M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl2M4D<-FCl2M4D+ theme(legend.position = "none")  
FCl2M4D

FCl2M4D<-FCl2M4D + annotate("text", label ="  0.0%", x = 1, y = 145 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label ="  5.3%", x = 2, y = 145 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label ="  7.2%", x = 3, y = 145 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label ="  4.8%", x = 5, y = 145 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" -2.4%", x = 10, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" -7.9%", x = 15, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" NS ", x = 2, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" NS ", x = 3, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" NS ", x = 5, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" NS ", x = 10, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D<-FCl2M4D + annotate("text", label =" NS ", x = 15, y = 150 , size = 4, colour ="#8d96a3", fontface = "bold")
FCl2M4D


Cl3M4<-MDataALL_H%>%filter(Cl4_1095==3)

FCl3M4D<-ggplot(Cl3M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(ImpType))) + 
  geom_jitter(position=position_jitter(0.1),alpha= 0.3, cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title=" ",caption = "Cohort PROFILE", tag = "C") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#edae49","#ffe5b4"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#  annotate("text", label ="67/140(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
#  annotate("text", label ="73/140(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl3M4D<-FCl3M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4D<-FCl3M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl3M4D<-FCl3M4D+ theme(legend.position = "none")  
FCl3M4D

FCl3M4D<-FCl3M4D + annotate("text", label ="0.0%", x = 1, y = 145 , size = 4, colour     ="#8d96a3", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" -12.5%", x = 2, y = 145 , size = 4, colour  ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" -16.9%", x = 3, y = 145 , size = 4, colour  ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" -19.1%", x = 5, y = 145 , size = 4, colour  ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" -21.7%", x = 10, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" -23.8%", x = 15, y = 145 , size = 4, colour ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" 8.5e-04", x = 2, y = 150 , size = 4, colour ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" 2.0e-06", x = 3, y = 150 , size = 4, colour ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" 4.9e-08", x = 5, y = 150 , size = 4, colour ="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" 3.9e-10", x = 10, y = 150 , size = 4, colour="#993333", fontface = "bold")
FCl3M4D<-FCl3M4D + annotate("text", label =" 6.2e-12", x = 15, y = 150 , size = 4, colour="#993333", fontface = "bold")
FCl3M4D

Cl4M4<-MDataALL_H%>%filter(Cl4_1095==4)

FCl4M4D<-ggplot(Cl4M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(ImpType))) + 
  geom_jitter(position=position_jitter(0.1),alpha= 0.3, cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title=" ",caption = "Cohort PROFILE", tag = "D") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#000080","#87CEEB"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#  annotate("text", label ="67/140(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
#  annotate("text", label ="73/140(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl4M4D<-FCl4M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4D<-FCl4M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl4M4D<-FCl4M4D+ theme(legend.position = "none")  
FCl4M4D

FCl4M4D<-FCl4M4D + annotate("text", label ="  0.0%", x = 1, y = 145 , size = 4, colour  = "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label =" -5.3%", x = 2, y = 145 , size = 4, colour  = "#993333", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label =" -1.9%", x = 3, y = 145 , size = 4, colour  = "#993333", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="  2.2%", x = 5, y = 145 , size = 4, colour  = "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="  1.4%", x = 10, y = 145 , size = 4, colour=  "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label =" -3.9%", x = 15, y = 145 , size = 4, colour=  "#993333", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="NS", x = 2, y = 150 , size = 4, colour = "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="NS", x = 3, y = 150 , size = 4, colour = "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="NS", x = 5, y = 150 , size = 4, colour = "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="NS", x = 10, y = 150 , size = 4, colour= "#8d96a3", fontface = "bold")
FCl4M4D<-FCl4M4D + annotate("text", label ="NS", x = 15, y = 150 , size = 4, colour= "#8d96a3", fontface = "bold")
FCl4M4D






textsup1 <- paste("Figure 3: Description of Self Organising Maps-derived clusters in PROFILE",
                 "Average lung functions of the four clusters obtained from Self Organising Maps of ",
                 "three years' spirometric data. Solid dark colour lines in each graph represent the",
                 "mean trajectory ppFVC imputed by the naive RF MCMC imputation model; the light ",
                 "solid colour lines are the imputed mean ppFVC trajectory obtained by the theoretical",
                 "RF MCMC imputation model. Dots: dark colour is the values obtained by the naïve",
                 "model RF MCMC and light colour represent the imputed values obtained by the",
                 "theoretical model on each patient. Each box blot represent mean (±SE) changes on",
                 "each time point on the same colour patterns. Graphs: A) Classical linear decline",
                 "(CL1, n=140 - green tone), B) Sinusoidal trajectory (CL2, n=100- red tone), C) Early",
                 "decline trajectory (CL3, n=113 -yellow tone) and D) Newtonian trajectory (CL4, n=62 -",
                 "- blue tone) - the top line in each graph indicates the total decline in annual FVC in %.", 
                 "The line below indicates the statistical significance of each time point ppFVC to baseline and", 
                 "the indicate p-value is Bonferroni corrected in each case." ,sep= " ")

textsup11 <- ggparagraph(text = textsup1 ,size = 15, color = "black")


dev.off()
dpi=300
tiff("Fig3.tiff", res=dpi, height=15*dpi, width=20*dpi)
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(7,8)))

define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(FCl1M4D, vp=define_region(1:3,1:4))
print(FCl2M4D, vp=define_region(1:3,5:8))
print(FCl3M4D, vp=define_region(4:6,1:4))
print(FCl4M4D, vp=define_region(4:6,5:8))
print(textsup11, vp=define_region(7:7,1:8))
dev.off()




dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Imputation effect 270921 Years All.tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4D,FCl2M4D,FCl3M4D,FCl4M4D, ncol=2)
dev.off()