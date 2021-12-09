#Project: IPF new project
#
# Purpose: Profile project: MCMC simulations for random forests
# Version: test
# Date: 03/06/2021
# Author: HPF
#
# Input: 
# Main data base:   
#
# Output: Work with the clusters
#         1) run a Monte Carlo (MC) simulation of future ppFVC
#         3) Find Conf intervals at 95% and 05% from 1 and 0
#         4) Find the SD, MAX, MIN and MEANs
#         5) Charactetion of the mean function of the cluster
#         6) MC simulation in base of basic population model.
#         6.1) N[t] = (N[t-1] * propDeclineppFCV(t-1) * rlnorm(1, 0, sd_decline(t))) * (1 - ppFVClossRate)
#         7) Graphical comparation between the fuctions, histograms, and line graphs
#         9) validation of the Simulation x1, x10 x100 by Generalized Linear Models with boottrap
#         10) posibilites of decline below 80 ppFVC in the simulation
#     
# Dependencies:@/home/svzhpf/IPF Project 150720
#
# Notes: missing biomarkers
#
#


# =    1 working space =======================================================================

setwd("F:/")
#if (! dir.exists("Project_IPF_SOM")) dir.create("Project_IPF_SOM")
setwd("F:/S2.2 PROFILE TEST MCMC 030621")
list.files("F:/S2.2 PROFILE TEST MCMC 030621")
getwd()

# = 2 packages need ===========================================================================

#Formatting/documentation packages
library(plyr)
library(doBy)
library(tidyverse)
library(broom)
library(LearnBayes)
library(polynom)
library(caret)
library(gridExtra)
library(RColorBrewer)
library(Metrics)
library(zoo)
library(MASS)
library(mc2d)

# = 3 open database ============================================================

MDataN<- read.csv("MDataV_Nottingham_120421_1452oX29v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataN)
dim(MDataN)
summary(MDataN)
MDataN$X=NULL
sapply(MDataN, function(x) sum(is.na(x)))

MDataB<- read.csv("MDataV_Brompton_120421_1050oX25v.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MDataB)
dim(MDataB)
summary(MDataB)
MDataB$X=NULL
sapply(MDataB, function(x) sum(is.na(x)))


# = 4 editing of the dataframe =================================================

MDataN<-data.frame(MDataN[ ,c("patientid","visit_N","ppFVC","Imp")])
MDataB<-data.frame(MDataB[ ,c("patientid","visit_N","ppFVC","Imp")])

MData<-rbind(MDataB,MDataN)
MData<-data.frame(MData%>%group_by(patientid)%>%filter(Imp==0)%>%ungroup())
MData<-(MData%>%mutate(visit_Idx=case_when(visit_N==1 ~  1, visit_N==3 ~ 2, visit_N==4 ~ 3,visit_N==5 ~ 5, visit_N==6 ~ 10, visit_N==7 ~ 15)))
MData<-MData%>%mutate(Imputation="No_Imputed")

FImp<-ggplot(MData, aes(x=visit_Idx, y=ppFVC, color=as.factor(Imputation))) + 
  geom_jitter(position=position_jitter(0.1), cex=3, alpha = 0.1)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="ppFVC No imputation",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name="Imputed",values=c("#00798c","#d1495b"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FImp<-FImp + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
FImp<-FImp+ stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), size= 0.3)   
FImp<-FImp+ theme(legend.position = "none")  
FImp


G1<-ggplot(MData, aes(ppFVC)) + geom_density(alpha = 0.2) +facet_wrap(~visit_Idx)   
G1

# = 5 characteristics complete dataset  ========================================


AlStats <- MData %>% group_by(visit_Idx) %>% dplyr:: summarise(n = n(),
                                                               mean = mean(ppFVC),
                                                               median = median(ppFVC),
                                                               max=max(ppFVC),
                                                               min=min(ppFVC),
                                                               sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                         CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

# = 6 Log Distribution of simulation parameters  ===============================

AlStats<-AlStats%>%mutate(ppFVC0=AlStats$mean[1])%>%mutate(ppFVC1=AlStats$mean[2])%>%mutate(ppFVC2=AlStats$mean[3])
AlStats<-AlStats%>%mutate(ppFVC3=AlStats$mean[4])%>%mutate(ppFVC4=AlStats$mean[5])%>%mutate(ppFVC5=AlStats$mean[6])
AlStats<-(AlStats%>%mutate(Grow0=mean/ppFVC0)%>%mutate(Grow1=mean/ppFVC1)%>%mutate(Grow2=mean/ppFVC2)%>%mutate(Grow3=mean/ppFVC3)
          %>%mutate(Grow4=mean/ppFVC4)%>%mutate(Grow5=mean/ppFVC5))

AlStats<-(AlStats%>%mutate(U0=(((mean-ppFVC0)/ppFVC0)/10)*-1)%>%mutate(U1=(((mean-ppFVC1)/ppFVC1)/10)*-1)%>%mutate(U2=(((mean-ppFVC2)/ppFVC2)/10)*-1)
          %>%mutate(U3=(((mean-ppFVC3)/ppFVC3)/10)*-1)%>%mutate(U4=(((mean-ppFVC4)/ppFVC4)/10)*-1)%>%mutate(U5=(((mean-ppFVC5)/ppFVC5)/10)*-1))  

AlStats<-AlStats%>%mutate(maxR=max/mean)%>%mutate(minR=min/mean)%>%mutate(medianR=median/mean)


# = 8 Simulation x1  ==============================================================


simulations <-(83*1)
# randomly generation

f_fvc <- function(){
  ppFVC <- c(AlStats$mean[1])
  RLoss0 <- c(AlStats$Grow0[1])
  max0<-c(AlStats$maxR[1])
  min0<-c(AlStats$minR[1])
  median0<-c(AlStats$medianR[1])
  Dist0  <- rpert(1, min = min0, mode = median0, max = max0, shape =1)
  ppFVC0 <- (ppFVC * (RLoss0) * Dist0)*(1)                                                          # first end-of-ppFVC
  RLoss1 <- c(AlStats$Grow0[2])                                                                     # Rate of loss ppFVC by 90 days (visit 3)
  RLoss2 <- c(AlStats$Grow1[3])
  RLoss3 <- c(AlStats$Grow2[4])
  RLoss4 <- c(AlStats$Grow3[5])
  RLoss5 <- c(AlStats$Grow4[6])
  PLoss1 <- c(AlStats$U0[2])                                                                        # Fix of FVC loss by 90 days (visit 3)
  PLoss2 <- c(AlStats$U1[3])
  PLoss3 <- c(AlStats$U2[4])
  PLoss4 <- c(AlStats$U3[5])
  PLoss5 <- c(AlStats$U4[6])
  max1<-c((AlStats$max[2]/(1-PLoss1))/(AlStats$mean[1]*RLoss1))
  max2<-c((AlStats$max[3]/(1-PLoss2))/(AlStats$mean[2]*RLoss2))
  max3<-c((AlStats$max[4]/(1-PLoss3))/(AlStats$mean[3]*RLoss3))
  max4<-c((AlStats$max[5]/(1-PLoss4))/(AlStats$mean[4]*RLoss4))
  max5<-c((AlStats$max[6]/(1-PLoss5))/(AlStats$mean[5]*RLoss5))
  min1<-c((AlStats$min[2]/(1-PLoss1))/(AlStats$mean[1]*RLoss1))
  min2<-c((AlStats$min[3]/(1-PLoss2))/(AlStats$mean[2]*RLoss2))
  min3<-c((AlStats$min[4]/(1-PLoss3))/(AlStats$mean[3]*RLoss3))
  min4<-c((AlStats$min[5]/(1-PLoss4))/(AlStats$mean[4]*RLoss4))
  min5<-c((AlStats$min[6]/(1-PLoss5))/(AlStats$mean[5]*RLoss5))
  Dist1 <- rpert(1, min = min1, mode =0.99, max = max1, shape =20)                            # log normalised STD on 90 days (visit 3) 
  Dist2 <- rpert(1, min = min2, mode =1, max = max2, shape =40)  
  Dist3 <- rpert(1, min = min3, mode =1, max = max3, shape =40)
  Dist4 <- rpert(1, min = min4, mode =1, max = max4, shape =40)  
  Dist5 <- rpert(1, min = min5, mode =1, max = max5, shape =40) 
  Sim1<- (ppFVC0 * (RLoss1) * Dist1)*(1-PLoss1)                                                    # Simulation (ppFVC) 90 days (visit 3)
  Sim2<- (Sim1 * (RLoss2) * Dist2)*(1-PLoss2)
  Sim3<- (Sim2 * (RLoss3) * Dist3)*(1-PLoss3)
  Sim4<- (Sim3 * (RLoss4) * Dist4)*(1-PLoss4)
  Sim5<- (Sim4 * (RLoss5) * Dist5)*(1-PLoss5)
  return(list(
    ppFVC0= (ppFVC0),
    Sim1= Sim1,
    Sim2= Sim2,
    Sim3= Sim3,
    Sim4= Sim4,
    Sim5= Sim5))
}

# Monte Carlo simulations
seeds <- c(308,17773,16882,8694,36702,65634)

fvc_dfAx1 <- data.frame(ppFVC0=rep(NA, simulations),
                        Sim1=rep(NA, simulations),
                        Sim2=rep(NA, simulations),
                        Sim3=rep(NA, simulations),
                        Sim4=rep(NA, simulations),
                        Sim5=rep(NA, simulations))

for (seed in seeds) {
  set.seed(seed)
  for (i in seq(simulations)){
    my_simulation <- f_fvc()
    fvc_dfAx1$ppFVC0[i] <- my_simulation$ppFVC0
    fvc_dfAx1$Sim1[i] <- my_simulation$Sim1
    fvc_dfAx1$Sim2[i] <- my_simulation$Sim2
    fvc_dfAx1$Sim3[i] <- my_simulation$Sim3
    fvc_dfAx1$Sim4[i] <- my_simulation$Sim4
    fvc_dfAx1$Sim5[i] <- my_simulation$Sim5
    
  }
}


fvc_dfAx1<-fvc_dfAx1%>%dplyr::mutate(patientid= 1:n())
fvc_dfAx1G<-fvc_dfAx1%>% gather(Simulation,ppFVC,ppFVC0,Sim1,Sim2,Sim3,Sim4,Sim5)
fvc_dfAx1G<-data.frame(fvc_dfAx1G[,c("patientid","Simulation","ppFVC")])
fvc_dfAx1G<-(fvc_dfAx1G%>%mutate(visit_Idx=case_when(Simulation=="ppFVC0" ~1 ,Simulation=="Sim1" ~ 2,
                                                     Simulation=="Sim2" ~ 3, Simulation=="Sim3" ~5,
                                                     Simulation=="Sim4" ~10, Simulation=="Sim5" ~15))%>%arrange(patientid))

Sim1All<-fvc_dfAx1G%>%mutate(Simulated="SimX1")
Sim1All<-data.frame(Sim1All%>%unite(patientid, c("patientid","Simulated")))
Sim1All<-Sim1All%>%mutate(Simulated="SimX1")
#Sim1All<-data.frame(Sim1All%>%group_by(patientid)%>%filter(ppFVC<139.14413)%>%ungroup())
#Sim1All<-data.frame(Sim1All%>%group_by(patientid)%>%filter(ppFVC>34.13491)%>%ungroup())

Sim1StatsAll <-Sim1All%>% group_by(visit_Idx) %>% dplyr:: summarise(n = n(),
                                                                    mean = mean(ppFVC),
                                                                    median = median(ppFVC),
                                                                    max=max(ppFVC),
                                                                    min=min(ppFVC),
                                                                    sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                              CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                              CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)


MDataA<-data.frame(MData[ ,c("patientid","visit_Idx","ppFVC")])
MDataA<-MDataA%>%mutate(Simulated="Actual")
Sim1All<-data.frame(Sim1All[ ,c("patientid","visit_Idx","ppFVC")])
Sim1All<-Sim1All%>%mutate(Simulated="SimX1")

GSim<-rbind(MDataA,Sim1All)

GSimDG<-ggplot(GSim, aes(ppFVC, fill =as.factor(Simulated))) + geom_density(alpha = 0.2) +facet_wrap(~visit_Idx)   
GSimDG



# = 2 Simulation x100  ==============================================================


simulations <-(83*100)
# randomly generation

f_fvc <- function(){
  ppFVC <- c(AlStats$mean[1])
  RLoss0 <- c(AlStats$Grow0[1])
  max0<-c(AlStats$maxR[1])
  min0<-c(AlStats$minR[1])
  median0<-c(AlStats$medianR[1])
  Dist0  <- rpert(1, min = min0, mode = median0, max = max0, shape =1)
  ppFVC0 <- (ppFVC * (RLoss0) * Dist0)*(1)                                                          # first end-of-ppFVC
  RLoss1 <- c(AlStats$Grow0[2])                                                                     # Rate of loss ppFVC by 90 days (visit 3)
  RLoss2 <- c(AlStats$Grow1[3])
  RLoss3 <- c(AlStats$Grow2[4])
  RLoss4 <- c(AlStats$Grow3[5])
  RLoss5 <- c(AlStats$Grow4[6])
  PLoss1 <- c(AlStats$U0[2])                                                                        # Fix of FVC loss by 90 days (visit 3)
  PLoss2 <- c(AlStats$U1[3])
  PLoss3 <- c(AlStats$U2[4])
  PLoss4 <- c(AlStats$U3[5])
  PLoss5 <- c(AlStats$U4[6])
  max1<-c((AlStats$max[2]/(1-PLoss1))/(AlStats$mean[1]*RLoss1))
  max2<-c((AlStats$max[3]/(1-PLoss2))/(AlStats$mean[2]*RLoss2))
  max3<-c((AlStats$max[4]/(1-PLoss3))/(AlStats$mean[3]*RLoss3))
  max4<-c((AlStats$max[5]/(1-PLoss4))/(AlStats$mean[4]*RLoss4))
  max5<-c((AlStats$max[6]/(1-PLoss5))/(AlStats$mean[5]*RLoss5))
  min1<-c((AlStats$min[2]/(1-PLoss1))/(AlStats$mean[1]*RLoss1))
  min2<-c((AlStats$min[3]/(1-PLoss2))/(AlStats$mean[2]*RLoss2))
  min3<-c((AlStats$min[4]/(1-PLoss3))/(AlStats$mean[3]*RLoss3))
  min4<-c((AlStats$min[5]/(1-PLoss4))/(AlStats$mean[4]*RLoss4))
  min5<-c((AlStats$min[6]/(1-PLoss5))/(AlStats$mean[5]*RLoss5))
  Dist1 <- rpert(1, min = min1, mode =0.99, max = max1, shape =20)                            # log normalised STD on 90 days (visit 3) 
  Dist2 <- rpert(1, min = min2, mode =1, max = max2, shape =40)  
  Dist3 <- rpert(1, min = min3, mode =1, max = max3, shape =40)
  Dist4 <- rpert(1, min = min4, mode =1, max = max4, shape =40)  
  Dist5 <- rpert(1, min = min5, mode =1, max = max5, shape =40) 
  Sim1<- (ppFVC0 * (RLoss1) * Dist1)*(1-PLoss1)                                                    # Simulation (ppFVC) 90 days (visit 3)
  Sim2<- (Sim1 * (RLoss2) * Dist2)*(1-PLoss2)
  Sim3<- (Sim2 * (RLoss3) * Dist3)*(1-PLoss3)
  Sim4<- (Sim3 * (RLoss4) * Dist4)*(1-PLoss4)
  Sim5<- (Sim4 * (RLoss5) * Dist5)*(1-PLoss5)
  return(list(
    ppFVC0= (ppFVC0),
    Sim1= Sim1,
    Sim2= Sim2,
    Sim3= Sim3,
    Sim4= Sim4,
    Sim5= Sim5))
}

# Monte Carlo simulations
seeds <- c(308,17773,16882,8694,36702,65634)

fvc_dfAx1 <- data.frame(ppFVC0=rep(NA, simulations),
                        Sim1=rep(NA, simulations),
                        Sim2=rep(NA, simulations),
                        Sim3=rep(NA, simulations),
                        Sim4=rep(NA, simulations),
                        Sim5=rep(NA, simulations))

for (seed in seeds) {
  set.seed(seed)
  for (i in seq(simulations)){
    my_simulation <- f_fvc()
    fvc_dfAx1$ppFVC0[i] <- my_simulation$ppFVC0
    fvc_dfAx1$Sim1[i] <- my_simulation$Sim1
    fvc_dfAx1$Sim2[i] <- my_simulation$Sim2
    fvc_dfAx1$Sim3[i] <- my_simulation$Sim3
    fvc_dfAx1$Sim4[i] <- my_simulation$Sim4
    fvc_dfAx1$Sim5[i] <- my_simulation$Sim5
    
  }
}


fvc_dfAx1<-fvc_dfAx1%>%dplyr::mutate(patientid= 1:n())
fvc_dfAx1G<-fvc_dfAx1%>% gather(Simulation,ppFVC,ppFVC0,Sim1,Sim2,Sim3,Sim4,Sim5)
fvc_dfAx1G<-data.frame(fvc_dfAx1G[,c("patientid","Simulation","ppFVC")])
fvc_dfAx1G<-(fvc_dfAx1G%>%mutate(visit_Idx=case_when(Simulation=="ppFVC0" ~1 ,Simulation=="Sim1" ~ 2,
                                                     Simulation=="Sim2" ~ 3, Simulation=="Sim3" ~5,
                                                     Simulation=="Sim4" ~10, Simulation=="Sim5" ~15))%>%arrange(patientid))



Sim100All<-fvc_dfAx1G%>%mutate(Simulated="SimX100")
Sim100All<-data.frame(Sim100All%>%unite(patientid, c("patientid","Simulated")))
Sim100All<-Sim100All%>%mutate(Simulated="SimX100")


Sim100StatsAll <-Sim100All%>% group_by(visit_Idx) %>% dplyr:: summarise(n = n(),
                                                                    mean = mean(ppFVC),
                                                                    median = median(ppFVC),
                                                                    max=max(ppFVC),
                                                                    min=min(ppFVC),
                                                                    sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                              CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                              CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)



Sim100All<-data.frame(Sim100All[ ,c("patientid","visit_Idx","ppFVC")])
Sim100All<-Sim100All%>%mutate(Simulated="SimX100")

GSim<-rbind(MDataA,Sim1All,Sim100All)

GSimDG<-ggplot(GSim, aes(ppFVC, fill =as.factor(Simulated))) + geom_density(alpha = 0.2) +facet_wrap(~visit_Idx)   
GSimDG



# = 2 Simulation x1000  ==============================================================


simulations <-(83*10)
# randomly generation

f_fvc <- function(){
  ppFVC <- c(AlStats$mean[1])
  RLoss0 <- c(AlStats$Grow0[1])
  max0<-c(AlStats$maxR[1])
  min0<-c(AlStats$minR[1])
  median0<-c(AlStats$medianR[1])
  Dist0  <- rpert(1, min = min0, mode = median0, max = max0, shape =1)
  ppFVC0 <- (ppFVC * (RLoss0) * Dist0)*(1)                                                          # first end-of-ppFVC
  RLoss1 <- c(AlStats$Grow0[2])                                                                     # Rate of loss ppFVC by 90 days (visit 3)
  RLoss2 <- c(AlStats$Grow1[3])
  RLoss3 <- c(AlStats$Grow2[4])
  RLoss4 <- c(AlStats$Grow3[5])
  RLoss5 <- c(AlStats$Grow4[6])
  PLoss1 <- c(AlStats$U0[2])                                                                        # Fix of FVC loss by 90 days (visit 3)
  PLoss2 <- c(AlStats$U1[3])
  PLoss3 <- c(AlStats$U2[4])
  PLoss4 <- c(AlStats$U3[5])
  PLoss5 <- c(AlStats$U4[6])
  max1<-c((AlStats$max[2]/(1-PLoss1))/(AlStats$mean[1]*RLoss1))
  max2<-c((AlStats$max[3]/(1-PLoss2))/(AlStats$mean[2]*RLoss2))
  max3<-c((AlStats$max[4]/(1-PLoss3))/(AlStats$mean[3]*RLoss3))
  max4<-c((AlStats$max[5]/(1-PLoss4))/(AlStats$mean[4]*RLoss4))
  max5<-c((AlStats$max[6]/(1-PLoss5))/(AlStats$mean[5]*RLoss5))
  min1<-c((AlStats$min[2]/(1-PLoss1))/(AlStats$mean[1]*RLoss1))
  min2<-c((AlStats$min[3]/(1-PLoss2))/(AlStats$mean[2]*RLoss2))
  min3<-c((AlStats$min[4]/(1-PLoss3))/(AlStats$mean[3]*RLoss3))
  min4<-c((AlStats$min[5]/(1-PLoss4))/(AlStats$mean[4]*RLoss4))
  min5<-c((AlStats$min[6]/(1-PLoss5))/(AlStats$mean[5]*RLoss5))
  Dist1 <- rpert(1, min = min1, mode =0.99, max = max1, shape =20)                            # log normalised STD on 90 days (visit 3) 
  Dist2 <- rpert(1, min = min2, mode =1, max = max2, shape =40)  
  Dist3 <- rpert(1, min = min3, mode =1, max = max3, shape =40)
  Dist4 <- rpert(1, min = min4, mode =1, max = max4, shape =40)  
  Dist5 <- rpert(1, min = min5, mode =1, max = max5, shape =40) 
  Sim1<- (ppFVC0 * (RLoss1) * Dist1)*(1-PLoss1)                                                    # Simulation (ppFVC) 90 days (visit 3)
  Sim2<- (Sim1 * (RLoss2) * Dist2)*(1-PLoss2)
  Sim3<- (Sim2 * (RLoss3) * Dist3)*(1-PLoss3)
  Sim4<- (Sim3 * (RLoss4) * Dist4)*(1-PLoss4)
  Sim5<- (Sim4 * (RLoss5) * Dist5)*(1-PLoss5)
  return(list(
    ppFVC0= (ppFVC0),
    Sim1= Sim1,
    Sim2= Sim2,
    Sim3= Sim3,
    Sim4= Sim4,
    Sim5= Sim5))
}

# Monte Carlo simulations
seeds <- c(308,17773,16882,8694,36702,65634)

fvc_dfAx1 <- data.frame(ppFVC0=rep(NA, simulations),
                        Sim1=rep(NA, simulations),
                        Sim2=rep(NA, simulations),
                        Sim3=rep(NA, simulations),
                        Sim4=rep(NA, simulations),
                        Sim5=rep(NA, simulations))

for (seed in seeds) {
  set.seed(seed)
  for (i in seq(simulations)){
    my_simulation <- f_fvc()
    fvc_dfAx1$ppFVC0[i] <- my_simulation$ppFVC0
    fvc_dfAx1$Sim1[i] <- my_simulation$Sim1
    fvc_dfAx1$Sim2[i] <- my_simulation$Sim2
    fvc_dfAx1$Sim3[i] <- my_simulation$Sim3
    fvc_dfAx1$Sim4[i] <- my_simulation$Sim4
    fvc_dfAx1$Sim5[i] <- my_simulation$Sim5
    
  }
}


fvc_dfAx1<-fvc_dfAx1%>%dplyr::mutate(patientid= 1:n())
fvc_dfAx1G<-fvc_dfAx1%>% gather(Simulation,ppFVC,ppFVC0,Sim1,Sim2,Sim3,Sim4,Sim5)
fvc_dfAx1G<-data.frame(fvc_dfAx1G[,c("patientid","Simulation","ppFVC")])
fvc_dfAx1G<-(fvc_dfAx1G%>%mutate(visit_Idx=case_when(Simulation=="ppFVC0" ~1 ,Simulation=="Sim1" ~ 2,
                                                     Simulation=="Sim2" ~ 3, Simulation=="Sim3" ~5,
                                                     Simulation=="Sim4" ~10, Simulation=="Sim5" ~15))%>%arrange(patientid))



Sim10All<-fvc_dfAx1G%>%mutate(Simulated="SimX1k")
Sim10All<-data.frame(Sim10All%>%unite(patientid, c("patientid","Simulated")))
Sim10All<-Sim10All%>%mutate(Simulated="SimX1k")


Sim1KStatsAll <-Sim10All%>% group_by(visit_Idx) %>% dplyr:: summarise(n = n(),
                                                                        mean = mean(ppFVC),
                                                                        median = median(ppFVC),
                                                                        max=max(ppFVC),
                                                                        min=min(ppFVC),
                                                                        sd = sd(ppFVC)) %>%mutate(sem = sd / sqrt(n - 1),
                                                                                                  CI_lower = mean + qt((1-0.95)/2, n - 1) *  sem, 
                                                                                                  CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)



Sim10All<-data.frame(Sim10All[ ,c("patientid","visit_Idx","ppFVC")])
Sim10All<-Sim10All%>%mutate(Simulated="SimX10")

GSim<-rbind(MDataA,Sim1All,Sim100All,Sim10All)

# = 10.1  Graphical characterisation of cluster Simulations ====================

GSimDG<-ggplot(GSim, aes(ppFVC, fill =as.factor(Simulated))) + geom_density(alpha = 0.2) +facet_wrap(~visit_Idx)   
GSimDG

Sim1G<-ggplot(GSim, aes(x=visit_Idx, y=ppFVC, color=as.factor(Simulated))) + 
  geom_jitter(position=position_jitter(0.1), cex=3, alpha = 0.1)+ theme_bw()+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="crossbar", width=0.5 )+
  labs(y="ppFVC ", title="ppFVC Simulation",caption = "Cohort PROFILE") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=15, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=15, hjust = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", color = "red", size=0.5)+
  theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=12, angle=0))+
  scale_color_manual(name="Simulation",values=c("#00798c","#d1495b","#edae49","#66a182"))+
  scale_y_continuous(breaks=seq(0,120,20))+
  scale_x_continuous(breaks=seq(0,12,4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
Sim1G<-Sim1G + scale_x_discrete(name ="Time (Days)", limits=c("00","90","180","","365","","","","","730","","","","","1095"))
Sim1G<- Sim1G + stat_smooth(method = "lm", se = TRUE, fill = NA, formula = y ~ poly(x, 3, raw = TRUE), size= 0.3)   
#Sim1G<-Sim1G  + theme(legend.position = "none")  
Sim1G



# = 15 Validations =================================================================================
# = 15.1 Regression with resampling GLM ============================================================


# configure caret training parameters to 200 bootstrap samples
fitControl <- trainControl(method = "boot", number = 20)

GSimM0<-GSim%>%filter(visit_Idx==0)
GSimM1<-GSim%>%filter(visit_Idx==1)
GSimM2<-GSim%>%filter(visit_Idx==2)
GSimM3<-GSim%>%filter(visit_Idx==3)
GSimM5<-GSim%>%filter(visit_Idx==5)
GSimM10<-GSim%>%filter(visit_Idx==10)
GSimM15<-GSim%>%filter(visit_Idx==15)

fit <- train(ppFVC ~ Simulated, method="glm",data=GSimM15, trControl = fitControl)

# print output object
fit

fit$resample[1:10,]


GSimAll0<-GSim%>%filter(visit_Idx==15)
observed_F_value<- anova(lm(GSimAll0$ppFVC ~ GSimAll0$Simulated))
observed_F_value

# = 10 save files =========================================================================


dev.off()
dpi=100
tiff("MC Global simuation 030621.tiff",  res=dpi, height=08*dpi, width=15*dpi)
grid.arrange(Sim1G,GSimDG ,ncol=2)
dev.off()

GSimT<- data.frame(GSim [,c("patientid","ppFVC","visit_Idx")])
GSimH <- spread(GSimT, patientid, ppFVC)

rownames(GSimH)<-GSimH$visit_Idx
GSimH$visit_Idx=NULL
GSimH<-data.frame(t(GSimH))
str(GSimH)

GSimH$patientid<-rownames(GSimH)
GSimH<-data.frame(GSimH[ ,c("patientid","X1","X2","X3","X5","X10","X15")])
GSim1<-(GSim%>%filter(visit_Idx==1))

GSimH<-join(GSim1,GSimH,by="patientid")
GSimH<-data.frame(GSimH[ ,c("patientid","X1","X2","X3","X5","X10","X15","Simulated")])

names(GSimH)[2]<-paste("Vs1")
names(GSimH)[3]<-paste("Vs3")
names(GSimH)[4]<-paste("Vs4")
names(GSimH)[5]<-paste("Vs5")
names(GSimH)[6]<-paste("Vs6")
names(GSimH)[7]<-paste("Vs7")

write.csv(GSimH, file="MDataV_Simulated_030621.csv", na="")