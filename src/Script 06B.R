setwd("F:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
list.files("F:/")
setwd("F:/T.12 lmer PROFIL internal validation 271021")
list.files("F:/T.12 lmer PROFIL internal validation 271021")

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
library(ggpubr)
library(ggplot2)


# Editing packages
library(dplyr)
library(plyr)
library(jaccard)
library(qvalue)
library(cluster)
library(factoextra)

# = 3 open database ============================================================

MDataP <- read.csv("MDataSoMV_3YAll_MCMC_415oX49v_200721.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataP)
dim(MDataP)
MDataP$X=NULL
MDataP<-MDataP%>%filter(Imp==0)

MDataS <- read.csv("MDataSoMV_3Y_MCMC_242oX49v_120721.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataS)
dim(MDataS)
MDataS$X=NULL
MDataS<-MDataS%>%filter(Imp==0)

MDataN <- read.csv("MDataSoMV_NoIMP_MCMC_82oX49v_130721.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataN)
dim(MDataN)
MDataN$X=NULL
MDataN<-MDataN%>%filter(Imp==0)


#define Jaccard Similarity function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


MPCl1<-MDataP%>%filter(Cl4_1095==1)
MSCl1<-MDataS%>%filter(Cl4==1)
MNCl1<-MDataN%>%filter(Cl4_NoIMP==1)

Cl1_PvN<-jaccard(MPCl1$patientid, MNCl1$patientid)
Cl1_PvS<-jaccard(MPCl1$patientid, MSCl1$patientid)
Cl1_SvN<-jaccard(MSCl1$patientid, MNCl1$patientid)


MPCl2<-MDataP%>%filter(Cl4_1095==2)
MSCl2<-MDataS%>%filter(Cl4A==2)
MNCl2<-MDataN%>%filter(Cl4M_NoIMP==2)

Cl2_PvN<-jaccard(MPCl2$patientid, MNCl2$patientid)
Cl2_PvS<-jaccard(MPCl2$patientid, MSCl2$patientid)
Cl2_SvN<-jaccard(MSCl2$patientid, MNCl2$patientid)

MPCl3<-MDataP%>%filter(Cl4_1095==3)
MSCl3<-MDataS%>%filter(Cl4==3)
MNCl3<-MDataN%>%filter(Cl4_NoIMP==3)

Cl3_PvN<-jaccard(MPCl3$patientid, MNCl3$patientid)
Cl3_PvS<-jaccard(MPCl3$patientid, MSCl3$patientid)
Cl3_SvN<-jaccard(MSCl3$patientid, MNCl3$patientid)

MPCl4<-MDataP%>%filter(Cl4_1095==4)
MSCl4<-MDataS%>%filter(Cl4==4)
MNCl4<-MDataN%>%filter(Cl4_NoIMP==4)

Cl4_PvN<-jaccard(MPCl4$patientid, MNCl4$patientid)
Cl4_PvS<-jaccard(MPCl4$patientid, MSCl4$patientid)
Cl4_SvN<-jaccard(MSCl4$patientid, MNCl4$patientid)



Jarc<-data.frame(c(Cl1_PvN,Cl1_PvS,Cl1_SvN,
                   Cl2_PvN,Cl2_PvS,Cl2_SvN,
                   Cl3_PvN,Cl3_PvS,Cl3_SvN,
                   Cl4_PvN,Cl4_PvS,Cl4_SvN))

names(Jarc)[1]<-paste("Jaccard_Coef")
str(Jarc)


Jarc2<-data.frame(c("Cl1_PvT","Cl1_PvS","Cl1_SvT",
                    "Cl2_PvT","Cl2_PvS","Cl2_SvT",
                    "Cl3_PvT","Cl3_PvS","Cl3_SvT",
                    "Cl4_PvT","Cl4_PvS","Cl4_SvT"))

names(Jarc2)[1]<-paste("Clusters")
str(Jarc)

Jarc<-cbind(Jarc,Jarc2)




f1 <- ggplot(Jarc, aes(x=Clusters, y=Jaccard_Coef, fill=as.factor(Clusters)))+
  geom_bar(stat="identity", position="dodge")+ theme_bw()+
  theme(legend.position = "none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y="Jaccard / Tanimoto coefficient",x="Cluster comparison" ,title="Jaccard / Tanimoto cluster similarity coefficient",caption = "Cohort PROFILE")+
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=6, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  geom_vline(xintercept = 3.5 , color="Red", linetype="dashed",size = 2)+
  geom_vline(xintercept = 6.5 , color="Blue", linetype="dashed",size = 2)+
  geom_vline(xintercept = 12.5 , color="Darkgreen", linetype="dashed",size = 2)+
  geom_vline(xintercept = 9.5 , color="Darkorange", linetype="dashed",size = 2)+
  geom_hline(yintercept = 0.5 , color="Black", linetype="dashed")
f1


dev.off()
dpi=100
tiff("Jaccard graphs 291021 .tiff",  res=dpi, height=14*dpi, width=21*dpi)
f1
dev.off()



