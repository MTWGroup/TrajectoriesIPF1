#Project: IPF new project
#
# Purpose: Profile project: verical data base with 7 visits
# Version: Gisly imputs (Brompton)
# Date: 12/06/2021
# Author:  HPF
#
# Input: 
#
#
# Output: limitiing the data to 1 year alive + imputation.
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
setwd("F:/T.5 PROFILE MCMC BOTH 160721")
list.files("F:/T.5 PROFILE MCMC BOTH 160721")

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
library(cowplot)
library(ggpubr)

#Formatting/documentation packages
library(plyr)
library(doBy)
library(reshape2)
library(stringr)
library(tidyverse)

# Analytic packages
library(factoextra)
library(NbClust)

# = 3 open database ============================================================

MWD<- read.csv("MDataSoMV_OptimModels_415oX11v_220621.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MWD)
dim(MWD)
MWD$X=NULL

# = 4 Select need it data from MData ===========================================

MWD365<-data.frame(MWD[1],MWD[20:23])
rownames(MWD365)<-MWD365$patientid
MWD365$patientid=NULL

MWD730<-data.frame(MWD[1],MWD[20:24])
rownames(MWD730)<-MWD730$patientid
MWD730$patientid=NULL

MWD1095<-data.frame(MWD[1],MWD[20:25])
rownames(MWD1095)<-MWD1095$patientid
MWD1095$patientid=NULL

MWDNOIMP<-MWD%>%filter(Imp==0)
MWDNOIMP<-data.frame(MWDNOIMP[1],MWDNOIMP[20:25])
rownames(MWDNOIMP)<-MWDNOIMP$patientid
MWDNOIMP$patientid=NULL

# = 5 Data preparation =========================================================

MWD365ASC <- scale(MWD365)
MWD730ASC <- scale(MWD730)
MWD1095ASC <- scale(MWD1095)
MWDNOIMPASC <- scale(MWDNOIMP)

# = 6 Analisis =================================================================


Elbow365<-fviz_nbclust(MWD365ASC, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  geom_hline(yintercept = 163, linetype = 2)+
  geom_hline(yintercept = 183, linetype = 2)+
  geom_hline(yintercept = 196, linetype = 2)+
  geom_hline(yintercept = 207, linetype = 2)+
  geom_hline(yintercept = 232, linetype = 2)+
  geom_hline(yintercept = 285, linetype = 2)+
  geom_hline(yintercept = 387, linetype = 2)+
  annotate("text", label ="WSS gap between clusters 3 vs 4 = 102", x = 7, y = 750 , size = 4, colour ="#993333", fontface = "bold")+
  labs(subtitle = " ", tag = "C", title = " ")

Elbow365

Elbow730<-fviz_nbclust(MWD730ASC, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  geom_hline(yintercept = 225, linetype = 2)+
  geom_hline(yintercept = 237, linetype = 2)+
  geom_hline(yintercept = 257, linetype = 2)+
  geom_hline(yintercept = 280, linetype = 2)+
  geom_hline(yintercept = 312, linetype = 2)+
  geom_hline(yintercept = 391, linetype = 2)+
  geom_hline(yintercept = 509, linetype = 2)+
  annotate("text", label ="WSS gap between clusters 3 vs 4 = 118", x = 7, y = 750 , size = 4, colour ="#993333", fontface = "bold")+
  labs(subtitle = " ", tag = "B", title = " ")

Elbow730


Elbow1095<-fviz_nbclust(MWD1095ASC, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  geom_hline(yintercept = 250, linetype = 2)+
  geom_hline(yintercept = 270, linetype = 2)+
  geom_hline(yintercept = 300, linetype = 2)+
  geom_hline(yintercept = 340, linetype = 2)+
  geom_hline(yintercept = 395, linetype = 2)+
  geom_hline(yintercept = 437, linetype = 2)+
  geom_hline(yintercept = 600, linetype = 2)+
  annotate("text", label ="WSS gap between clusters 3 vs 4 = 167", x = 7, y = 750 , size = 4, colour ="#993333", fontface = "bold")+
  labs(subtitle = " ", tag = "A", title = " ")

Elbow1095


ElbowNOIMP<-fviz_nbclust(MWDNOIMPASC, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  geom_hline(yintercept = 55, linetype = 2)+
  geom_hline(yintercept = 65, linetype = 2)+
  geom_hline(yintercept = 68, linetype = 2)+
  geom_hline(yintercept = 71, linetype = 2)+
  geom_hline(yintercept = 75, linetype = 2)+
  geom_hline(yintercept = 85, linetype = 2)+
  geom_hline(yintercept = 114, linetype = 2)+
  annotate("text", label ="WSS gap between clusters 3 vs 4 = 29", x = 7, y = 400 , size = 4, colour ="#993333", fontface = "bold")+
  labs(subtitle = " ", tag = "D", title = " ")

ElbowNOIMP

dev.off()
dpi=100
tiff("Clusters 270921 All Patients.tiff",  res=dpi, height=12*dpi, width=20*dpi)
grid.arrange(Elbow1095,Elbow730,Elbow365,ElbowNOIMP,ncol=2)
dev.off()


textsup1 <- paste("Supplementary Figure 2: Within-Sum of Square inertia  (or  variance) graphs for the cluster number estimation generated", 
                  "by the elbow algorithm.", 
                  "Within-Sum of Square inertia (or variance) graphs for the cluster data partition generated by  the  elbow  algorithm,",
                  "from  1  to 10 clusters for the whole PROFILE lung function cohort. Each graph represents a different lung function data partition:",  
                  "A) 1095 days, B) 730 days, C) 360 days or D) 1095 days data of not imputed and complete spirometric records.",
                  "The estimated curvature is displayed as a dashed solid line, and shows in each data partition case a stabilisation decline (or inertia)", 
                  "at 4 clusters indicated by a dashed vertical line.", sep= " ")

textsup11 <- ggparagraph(text = textsup1 ,size = 20, color = "black")


dev.off()
dpi=300
tiff("Supp Fig 2.tiff", res=dpi, height=20*dpi, width=20*dpi)
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(8,8)))

define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(Elbow1095, vp=define_region(1:3,1:4))
print(Elbow730, vp=define_region(1:3,5:8))
print(Elbow365, vp=define_region(4:6,1:4))
print(ElbowNOIMP, vp=define_region(4:6,5:8))
print(textsup11, vp=define_region(8:8,1:8))
dev.off()
