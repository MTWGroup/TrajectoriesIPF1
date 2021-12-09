#
# Purpose: Profile project: verical data base with 7 visits
# Version: Gisly imputs (Brompton)
# Date: 12/06/2021
# Author:  HPF
#
# Input: 
#
#
# Output: limitiing the data to 1 year alive + imputation. SOM
# Add: mising data from visits : 3&4
#
# Dependencies:@/F:PROFILE FINAL IMP SC 120421/Stage 1 Data Set
#
# Notes: use the new dataset, cuold be a problem with the dates.
#
# =    1 working space =========================================================
setwd("F:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
list.files("F:/")
setwd("F:/T.1 PROFILE  12 month alive  Clusters 120621")
list.files("F:/T.1 PROFILE  12 month alive  Clusters 120621")

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
library(tidyverse)

# Machine learing package
library(som)

# = 3 open database ============================================================

MWD<- read.csv("MDataSoMV_OptimModels_Syntetic_415oX48v_020721.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD)
dim(MWD)
sapply(MWD, function(x) sum(is.na(x)))
rownames(MWD)<-MWD$patientid
MWD$X=NULL
MWD<-(MWD[1:40])

# = 4 Select need it data from MData ===========================================

MWD1095A<-MWD%>%filter(DeadFromV1>=3)
MWD1095A1<-data.frame(MWD1095A[1],MWD1095A[20:25])
rownames(MWD1095A1)<-MWD1095A1$patientid
MWD1095A1$patientid=NULL

# = 5 Normalize the data =======================================================

RSoM1095A<- som::normalize(MWD1095A1)

RSoM1095A<-data.frame(RSoM1095A)
RSoM1095A<-RSoM1095A%>%mutate(patientid=rownames(RSoM1095A))
RSoM1095A<-join(MWD1095A,RSoM1095A, by="patientid")
names(RSoM1095A)[41]<-paste("Vs1N")
names(RSoM1095A)[42]<-paste("Vs3N")
names(RSoM1095A)[43]<-paste("Vs4N")
names(RSoM1095A)[44]<-paste("Vs5N")
names(RSoM1095A)[45]<-paste("Vs6N")
names(RSoM1095A)[46]<-paste("Vs7N")
rownames(RSoM1095A)<-RSoM1095A$patientid

RSoM1095A<-(RSoM1095A%>%mutate(Diedin4Years=case_when(DeadFromV1>4 ~ 0, DeadFromV1<4 ~ 1, DeadFromV1==4 ~0)))

# = 6 data preparation and selection ===========================================


set.seed(222)
MAPSoMCl4 <- som(RSoM1095A[41:46], xdim=2, ydim=2 ,alpha = c(0.05, 0.01), radius = 1 , topol="rect", neigh="bubble")

dev.off()
dpi=300
tiff("SoMCl41095A.tiff", res=dpi, height=20*dpi, width=35*dpi)
plot(MAPSoMCl4,yadj=0.15, main="FVC profiles from visits 1 to 7 obtained by Kohonen Neural Network Cl4 complete", 
     xlab="visits: Vs1,Vs3,Vs4,Vs5,Vs6,Vs7")
dev.off()

RSoMCl4<-data.frame(MAPSoMCl4$visual)
RSoMCl4<-RSoMCl4%>%mutate(patientid=RSoM1095A$patientid)
names(RSoMCl4)[3]<-paste("qErrorCl4")
RSoMCl4<-(RSoMCl4%>%mutate(Cl4=case_when( y==0 & x==1 ~ 2,  y==0 & x==0 ~ 1, y==1 & x==1 ~ 4,  y==1 & x==0 ~ 3)))
RSoMCl4<-data.frame(RSoMCl4[,c("patientid","Cl4","qErrorCl4")])
RSoM1095A<-join(RSoM1095A,RSoMCl4, by="patientid")

# =  7 graphical cluster analysis ==============================================

#General
RSoMVs1<- data.frame(RSoM1095A[,c("patientid","Vs1","Vs1N","Cl4","Diedin4Years")])
RSoMVs1<-RSoMVs1%>%mutate(visit_Idx=1)
names(RSoMVs1)[2]<-paste("ppFVC")
names(RSoMVs1)[3]<-paste("NppFVC")

RSoMVs3<- data.frame(RSoM1095A[,c("patientid","Vs3","Vs3N","Cl4","Diedin4Years")])
RSoMVs3<-RSoMVs3%>%mutate(visit_Idx=2)
names(RSoMVs3)[2]<-paste("ppFVC")
names(RSoMVs3)[3]<-paste("NppFVC")

RSoMVs4<- data.frame(RSoM1095A[,c("patientid","Vs4","Vs4N","Cl4","Diedin4Years")])
RSoMVs4<-RSoMVs4%>%mutate(visit_Idx=3)
names(RSoMVs4)[2]<-paste("ppFVC")
names(RSoMVs4)[3]<-paste("NppFVC")


RSoMVs5<- data.frame(RSoM1095A[,c("patientid","Vs5","Vs5N","Cl4","Diedin4Years")])
RSoMVs5<-RSoMVs5%>%mutate(visit_Idx=5)
names(RSoMVs5)[2]<-paste("ppFVC")
names(RSoMVs5)[3]<-paste("NppFVC")


RSoMVs6<- data.frame(RSoM1095A[,c("patientid","Vs6","Vs6N","Cl4","Diedin4Years")])
RSoMVs6<-RSoMVs6%>%mutate(visit_Idx=10)
names(RSoMVs6)[2]<-paste("ppFVC")
names(RSoMVs6)[3]<-paste("NppFVC")

RSoMVs7<- data.frame(RSoM1095A[,c("patientid","Vs7","Vs7N","Cl4","Diedin4Years")])
RSoMVs7<-RSoMVs7%>%mutate(visit_Idx=15)
names(RSoMVs7)[2]<-paste("ppFVC")
names(RSoMVs7)[3]<-paste("NppFVC")

MWDSoM417<-rbind(RSoMVs1,RSoMVs3,RSoMVs4,RSoMVs5, RSoMVs6, RSoMVs7)
MWDSoM417<-MWDSoM417%>%arrange(patientid,visit_Idx)

# = 9.1 graphs Cluster 3 M4 ====================================================

ClM4<- MWDSoM417 %>% group_by(Cl4, visit_Idx, Diedin4Years) %>% dplyr:: summarise(n = n(),
                                                                                  mean = mean(ppFVC),
                                                                                  median = median(ppFVC),
                                                                                  sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                                            CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                                            CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

qErrorCl4<- RSoMCl4 %>% group_by(Cl4) %>% dplyr:: summarise(n = n(),
                                                            mean = mean(qErrorCl4),
                                                            median = median(qErrorCl4),
                                                            sd = sd(qErrorCl4)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                          CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                          CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)



Cl1M4<-MWDSoM417%>%filter(Cl4==1)

FCl1M4N<-ggplot(Cl1M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 1 Model 4 structure - 3 Years ",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl1M4N<-FCl1M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4N<-FCl1M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl1M4N <- FCl1M4N + annotate("text", label ="n=89/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl1M4N <- FCl1M4N + annotate("text", label ="qError (0.882±0.03)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl1M4N

FCl1M4D<-ggplot(Cl1M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin4Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 1 Model 4 - ppFVC Died < 4 Yrs - 3 Years",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="65/89(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="24/89(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl1M4D<-FCl1M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4D<-FCl1M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl1M4D<-FCl1M4D+ theme(legend.position = "none")  
FCl1M4D


FCl1M4M<-ggplot(Cl1M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 1 Model 4 - 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl1M4M<-FCl1M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4M<-FCl1M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl1M4D<- FCl1M4D+ theme(legend.position = "none")  
FCl1M4M


# = 9.2 graphs Cluster 2 M4 ====================================================

Cl2M4<-MWDSoM417%>%filter(Cl4==2)

FCl2M4N<-ggplot(Cl2M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 2 Model 4 structure - 3 Years ",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl2M4N<-FCl2M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4N<-FCl2M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl2M4N <- FCl2M4N + annotate("text", label ="n=59/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl2M4N <- FCl2M4N + annotate("text", label ="qError (1.521±0.05)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl2M4N

FCl2M4D<-ggplot(Cl2M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin4Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 2 Model 4 - ppFVC Died < 4 Yrs - 3 Years ",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="51/59(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="08/59(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl2M4D<-FCl2M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4D<-FCl2M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl2M4D<- FCl2M4D+ theme(legend.position = "none")  
FCl2M4D

FCl2M4M<-ggplot(Cl2M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 2 Model 4- 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl2M4M<-FCl2M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4M<-FCl2M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl2M4D<- FCl2M4D+ theme(legend.position = "none")  
FCl2M4M


# = 9.3 graphs Cluster 3 M4 ====================================================

Cl3M4<-MWDSoM417%>%filter(Cl4==3)

FCl3M4N<-ggplot(Cl3M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 3 Model 4 structure - 3 Years ",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl3M4N<-FCl3M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4N<-FCl3M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl3M4N <- FCl3M4N + annotate("text", label ="n=61/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl3M4N <- FCl3M4N + annotate("text", label ="qError(1.228±0.061)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl3M4N


FCl3M4D<-ggplot(Cl3M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin4Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 3 Model 4 - ppFVC Died < 4 Yrs- 3 Years",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="54/61(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="07/61(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl3M4D<-FCl3M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4D<-FCl3M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl3M4D<- FCl3M4D+ theme(legend.position = "none")  
FCl3M4D

FCl3M4M<-ggplot(Cl3M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 3 Model 4 - 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#edae49"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl3M4M<-FCl3M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4M<-FCl3M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl3M4D<- FCl3M4D+ theme(legend.position = "none")  
FCl3M4M

# = 9.4 graphs Cluster 4 M4 ====================================================

Cl4M4<-MWDSoM417%>%filter(Cl4==4)

FCl4M4N<-ggplot(Cl4M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 4 Model 4 structure - 3 Years",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl4M4N<-FCl4M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4N<-FCl4M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl4M4N <- FCl4M4N + annotate("text", label ="n=33/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl4M4N <- FCl4M4N + annotate("text", label ="qError (2.161±0.058)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl4M4N

FCl4M4D<-ggplot(Cl4M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin4Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 4 Model 4 - ppFVC Died < 3 Yrs - 3 Years",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="26/33(4Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="07/33(4Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl4M4D<-FCl4M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4D<-FCl4M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl4M4D<- FCl4M4D+ theme(legend.position = "none")  
FCl4M4D

FCl4M4M<-ggplot(Cl4M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 4 Model 4 - 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#66a182"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl4M4M<-FCl4M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4M<-FCl4M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl4M4D<- FCl4M4D+ theme(legend.position = "none")  
FCl4M4M

dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Raw 200721 3 Years .tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4N,FCl2M4N,FCl3M4N,FCl4M4N, ncol=2)
dev.off()

dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Mortality 200721 3 Years .tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4D,FCl2M4D,FCl3M4D,FCl4M4D, ncol=2)
dev.off()

dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Mean 200721 3 Years .tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4M,FCl2M4M,FCl3M4M,FCl4M4M, ncol=2)
dev.off()

write.csv(MWDSoM417, file="MDataSoMH_3Y_MCMC_492oX7v_120721.csv", na="")
write.csv(RSoM1095A, file="MDataSoMV_3Y_MCMC_242oX49v_120721.csv", na="")


# = 9.1 graphs Cluster 3 M4 ====================================================

ClM4A<- MWDSoM417 %>% group_by(Cl4A, visit_Idx, Diedin5Years) %>% dplyr:: summarise(n = n(),
                                                                                    mean = mean(ppFVC),
                                                                                    median = median(ppFVC),
                                                                                    sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                                              CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                                              CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

qErrorCl4A<- RSoMCl4 %>% group_by(Cl4) %>% dplyr:: summarise(n = n(),
                                                             mean = mean(qErrorCl4),
                                                             median = median(qErrorCl4),
                                                             sd = sd(qErrorCl4)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                           CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                           CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)



Cl1M4<-MWDSoM417%>%filter(Cl4A==1)

FCl1M4N<-ggplot(Cl1M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 1 Model 4 structure - 3 Years ",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl1M4N<-FCl1M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4N<-FCl1M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl1M4N <- FCl1M4N + annotate("text", label ="n=68/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl1M4N <- FCl1M4N + annotate("text", label ="qError (0.882±0.03)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl1M4N

FCl1M4D<-ggplot(Cl1M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin5Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 1 Model 4 - ppFVC Died < 5 Yrs - 3 Years",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="37/68(5Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="31/68(5Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl1M4D<-FCl1M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4D<-FCl1M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl1M4D<-FCl1M4D+ theme(legend.position = "none")  
FCl1M4D


FCl1M4M<-ggplot(Cl1M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4A))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 1 Model 4 - 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl1M4M<-FCl1M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl1M4M<-FCl1M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl1M4D<- FCl1M4D+ theme(legend.position = "none")  
FCl1M4M


# = 9.2 graphs Cluster 2 M4 ====================================================

Cl2M4<-MWDSoM417%>%filter(Cl4A==2)

FCl2M4N<-ggplot(Cl2M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 2 Model 4 structure - 3 Years ",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl2M4N<-FCl2M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4N<-FCl2M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl2M4N <- FCl2M4N + annotate("text", label ="n=75/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl2M4N <- FCl2M4N + annotate("text", label ="qError (1.521±0.05)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl2M4N

FCl2M4D<-ggplot(Cl2M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin5Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 2 Model 4 - ppFVC Died < 5 Yrs - 3 Years ",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="41/75(5Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="34/75(5Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl2M4D<-FCl2M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4D<-FCl2M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl2M4D<- FCl2M4D+ theme(legend.position = "none")  
FCl2M4D

FCl2M4M<-ggplot(Cl2M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4A))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 2 Model 4- 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl2M4M<-FCl2M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl2M4M<-FCl2M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl2M4D<- FCl2M4D+ theme(legend.position = "none")  
FCl2M4M


# = 9.3 graphs Cluster 3 M4 ====================================================

Cl3M4<-MWDSoM417%>%filter(Cl4A==3)

FCl3M4N<-ggplot(Cl3M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 3 Model 4 structure - 3 Years ",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl3M4N<-FCl3M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4N<-FCl3M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl3M4N <- FCl3M4N + annotate("text", label ="n=40/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl3M4N <- FCl3M4N + annotate("text", label ="qError(1.228±0.061)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl3M4N


FCl3M4D<-ggplot(Cl3M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin5Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 3 Model 4 - ppFVC Died < 5 Yrs- 3 Years",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="27/40(5Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="13/40(5Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl3M4D<-FCl3M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4D<-FCl3M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl3M4D<- FCl3M4D+ theme(legend.position = "none")  
FCl3M4D

FCl3M4M<-ggplot(Cl3M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4A))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 3 Model 4 - 3 Years",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#edae49"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl3M4M<-FCl3M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl3M4M<-FCl3M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl3M4D<- FCl3M4D+ theme(legend.position = "none")  
FCl3M4M

# = 9.4 graphs Cluster 4 M4 ====================================================

Cl4M4<-MWDSoM417%>%filter(Cl4A==4)

FCl4M4N<-ggplot(Cl4M4,aes(x =visit_Idx, y = NppFVC, colour = factor(patientid))) +
  geom_line(size = 1, alpha=0.2)+ geom_point(size = 0.5, alpha=0.2) + theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.position = "none")+
  labs(y="NppFVC ", title="Cluster 4 Model 4 structure - 3 Years",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))
FCl4M4N<-FCl4M4N+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4N<-FCl4M4N+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), colour = "darkred")   
FCl4M4N <- FCl4M4N + annotate("text", label ="n=59/242", x = 11, y = -2 , size = 5, colour = "darkgreen", fontface = "bold")
FCl4M4N <- FCl4M4N + annotate("text", label ="qError (2.161±0.058)", x = 11, y = -2.2 , size = 5, colour = "darkred", fontface = "bold")
FCl4M4N

FCl4M4D<-ggplot(Cl4M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Diedin5Years))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 4 Model 4 - ppFVC Died < 5 Yrs - 3 Years",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", label ="44/59(5Yr>)", x = 11, y = 30 , size = 5, colour = "#00798c", fontface = "bold")+
  annotate("text", label ="15/59(5Yr<)", x = 11, y = 26 , size = 5, colour = "#d1495b", fontface = "bold")
FCl4M4D<-FCl4M4D + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4D<-FCl4M4D + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
FCl4M4D<- FCl4M4D+ theme(legend.position = "none")  
FCl4M4D

FCl4M4M<-ggplot(Cl4M4, aes(x=visit_Idx, y=ppFVC, color=as.factor(Cl4A))) + 
  geom_jitter(position=position_jitter(0.1), cex=2)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="Cluster 4 Model 4",caption = "Cohort PROFILE") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name= " " ,values=c("#66a182"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCl4M4M<-FCl4M4M+ scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FCl4M4M<-FCl4M4M+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE))   
#FCl4M4D<- FCl4M4D+ theme(legend.position = "none")  
FCl4M4M

dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Raw 120721 3 Years M4C .tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4N,FCl2M4N,FCl3M4N,FCl4M4N, ncol=2)
dev.off()

dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Mortality 120721 3 Years MC4 .tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4D,FCl2M4D,FCl3M4D,FCl4M4D, ncol=2)
dev.off()

dev.off()
dpi=100
tiff("SOM MCMC visit 1 to 7 PROFILE Model 4 Cluster Mean 120721 3 Years M4C .tiff",  res=dpi, height=14*dpi, width=21*dpi)
grid.arrange(FCl1M4M,FCl2M4M,FCl3M4M,FCl4M4M, ncol=2)
dev.off()