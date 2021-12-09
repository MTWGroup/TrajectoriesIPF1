#Project: IPF new project
#
# Purpose: Data surv Notts SOM  
# Version: 1
# Date: 22/06/2021
# Author:  HPF
#
# Input: Main data base:
#
# Output: survival analysis of the clusters generated for all profile
#         general analysis
#         poc hod 
#         Testing proportional Hazards assumption
#         
#
#
# Dependencies:@/F:S5.1 PROFILE MCMC SURV SC 220521
#
# Notes: use the new dataset, could be a problem with the dates.
#
# =    1 working space =========================================================
setwd("F:/")
#if (! dir.exists("Project_IPF_SOM")) dir.create("Project_IPF_SOM")
setwd("F:/T.5 PROFILE MCMC BOTH 160721")
list.files("F:/T.5 PROFILE MCMC BOTH 160721")
getwd()

# = 2 packages need ============================================================
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("impute"))

#Analytical packages
library("RColorBrewer")
library("affycoretools")
library("survminer")
library("FactoMineR")
library("factoextra")
library("survival")  

#Formatting/documentation packages
library("plyr")
library("doBy")
library("reshape2")
library("stringr")
library("tidyverse")
library("gridExtra")

library(lattice)
library(cowplot)
library(ggpubr)
library(grid)
# = 3 open database ============================================================

MData<- read.csv("MDataSoMV_3YAll_MCMC_415oX49v_110821.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
str(MData)
dim(MData)

MData<-(MData%>%mutate(DeadFromIn5Y=case_when(DeadFromV1>5 ~ 5, DeadFromV1<5 ~ DeadFromV1, DeadFromV1==5 ~5.00)))
# = 4 Survival Analysis and Visualization between sites ========================

fit <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ site, data = MData)


ggsurvplot(fit, data = MData)
ggsurvplot(fit, data = MData, censor.shape="|", censor.size = 4)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(face="bold", size=13, hjust = 0.0, vjust = 1),
      axis.title.x = element_text(face="bold", size=13, hjust = 0.5, vjust = 0.4),
      axis.title.y = element_text(face="bold", size= 13, hjust = 0.5, angle=90, vjust = 3),
      axis.text.x = element_text(face="bold", color="#993333",size=10, angle=0),
      axis.text.y = element_text(face="bold", color="#993333",size=10, angle=0)
    )
}

ggsurvA <- ggsurvplot(
  fit,                  # survfit object with calculated statistics.
  data =  MData,        # data used to fit survival curves.
  title = "PROFILE Site Survival Analysis",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#00798c","#d1495b"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Brompton", "Nottingham")    # change legend labels.
)

ggsurvA

# = 4.1 fit a cox proportional hazards model and plot ==========================

res.coxS <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ site, data = MData)
res.coxS
test.phS <- cox.zph(res.coxS)
test.phS
ggcoxzph(test.phS)


# = 4.2 Save graphs ============================================================

dev.off()
dpi=100
tiff("SURV All Site PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA
dev.off()

# = 5 Survival Analysis and Visualization between sites ========================

fit <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4_1095, data = MData)
fit1 <- survfit(Surv(DeadFromV1,Status) ~ Cl4_1095, data = MData)


ggsurvplot(fit, data = MData)
ggsurvplot(fit, data = MData, censor.shape="|", censor.size = 4)

ggsurvA <- ggsurvplot(
  fit,                  # survfit object with calculated statistics.
  data =  MData,        # data used to fit survival curves.
  title = "A",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#00798c","#d1495b","#edae49","#000080"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = F,# show bars instead of names in text annotations
  risk.table.title="B ",

   # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c("Cluster 1", "Cluster 2","Cluster 3", "Cluster 4")# change legend labels.
  )

ggsurvA

ggsurvA$table <- ggsurvA$table +
                 theme(plot.caption = element_text(hjust = 0, size = 13) ) +
                labs(caption = "Figure 4: Kaplan-Meir survival graphs.
Kaplan-Meier estimates of survival of IPF patients based on A) cluster allocation. B) the number of deaths 
is shown for every cluster during the following five years.Survival for cluster 4 was significantly 
different from that all the other clusters: vs Cl1-p<0.0001, vs Cl2 -p=0.033 and vs Cl3 -p<0.0001.") 

  



ggsurvA


dev.off()
dpi=300
tiff("Fig4.tiff", res=dpi, height=10*dpi, width=9*dpi)
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(7,7)))

define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(ggsurvA, vp=define_region(1:7,1:7))
dev.off()


res.coxS <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ site, data = MData)
res.coxS
test.phS <- cox.zph(res.coxS)
test.phS
ggcoxzph(test.phS)



# = 5.2 Save graphs ============================================================

dev.off()
dpi=100
tiff("SURV All Clusters PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA
dev.off()


# = 6.1 Test Between cluster 4 vs 6 ============================================

MDataS1A <-MData%>%filter(Cl4==1)
MDataS4A <-MData%>%filter(Cl4==4)
MDataS1vs4<-rbind(MDataS1A,MDataS4A)

fit1Vs4 <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS1vs4)

ggsurvplot(fit1Vs4, data = MDataS1vs4)
ggsurvplot(fit1Vs4, data = MDataS1vs4, censor.shape="|", censor.size = 4)


ggsurvA1Vs4 <- ggsurvplot(
  fit1Vs4,                  # survfit object with calculated statistics.
  data =  MDataS1vs4,        # data used to fit survival curves.
  title = "PROFILE Survival Analysis Clusters 1 and 4",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#00798c","#000080"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Cluster 1", "Cluster 4")    # change legend labels.
)

ggsurvA1Vs4

# = 6.2 fit a cox proportional hazards model and plot ==========================

res.coxCl1VsCl4 <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS1vs4)
res.coxCl1VsCl4
test.phCl1VsCl4 <- cox.zph(res.coxCl1VsCl4)
test.phCl1VsCl4
ggcoxzph(test.phCl1VsCl4)

# = 6.3 Save graphs ============================================================

dev.off()
dpi=100
tiff("SURV 1 vs 4 Cluster PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA1Vs4
dev.off()

# = 7.1 Test Between cluster 3 vs 4 ============================================

MDataS3A <-MData%>%filter(Cl4==3)
MDataS4A <-MData%>%filter(Cl4==4)
MDataS3vs4<-rbind(MDataS3A,MDataS4A)

fit3Vs4 <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS3vs4)

ggsurvplot(fit3Vs4, data = MDataS3vs4)
ggsurvplot(fit3Vs4, data = MDataS3vs4, censor.shape="|", censor.size = 4)


ggsurvA3Vs4 <- ggsurvplot(
  fit3Vs4,                   # survfit object with calculated statistics.
  data =  MDataS3vs4,        # data used to fit survival curves.
  title = "PROFILE Survival Analysis Clusters 3 and 4",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#edae49","#000080"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Cluster 3", "Cluster 4")    # change legend labels.
)

ggsurvA3Vs4

res.coxCl3VsCl4 <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS3vs4)
res.coxCl3VsCl4
test.phCl3VsCl4 <- cox.zph(res.coxCl3VsCl4)
test.phCl3VsCl4
ggcoxzph(test.phCl3VsCl4)


# = 7.3 Save graphs ============================================================

dev.off()
dpi=100
tiff("SURV 3 vs 4 Cluster PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA3Vs4
dev.off()

# = 8.1 Test Between cluster 2 vs 6 ============================================

MDataS2A <-MData%>%filter(Cl4==2)
MDataS4A <-MData%>%filter(Cl4==4)
MDataS2vs4<-rbind(MDataS2A,MDataS4A)

fit2Vs4 <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS2vs4)

ggsurvplot(fit2Vs4, data = MDataS2vs4)
ggsurvplot(fit2Vs4, data = MDataS2vs4, censor.shape="|", censor.size = 4)


ggsurvA2Vs4 <- ggsurvplot(
  fit2Vs4,                   # survfit object with calculated statistics.
  data =  MDataS2vs4,        # data used to fit survival curves.
  title = "PROFILE Survival Analysis Clusters 2 and 4",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#d1495b","#000080"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Cluster 2", "Cluster 4")    # change legend labels.
)

ggsurvA2Vs4

# = 8.2 fit a cox proportional hazards model and plot ==========================

res.coxCl2VsCl4 <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS2vs4)
res.coxCl2VsCl4
test.phCl2VsCl4 <- cox.zph(res.coxCl2VsCl4)
test.phCl2VsCl4
ggcoxzph(test.phCl2VsCl4)

# = 8.3 Save graphs ============================================================

dev.off()
dpi=100
tiff("SURV 2 vs 4 Cluster PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA2Vs4
dev.off()

# = 9.1 Test Between cluster 5 vs 6 ============================================

MDataS1A <-MData%>%filter(Cl4==1)
MDataS2A <-MData%>%filter(Cl4==2)
MDataS1vs2<-rbind(MDataS1A,MDataS2A)

fit1Vs2 <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS1vs2)

ggsurvplot(fit1Vs2, data = MDataS1vs2)
ggsurvplot(fit1Vs2, data = MDataS1vs2, censor.shape="|", censor.size = 4)


ggsurvA1Vs2 <- ggsurvplot(
  fit1Vs2,                   # survfit object with calculated statistics.
  data =  MDataS1vs2,        # data used to fit survival curves.
  title = "PROFILE Survival Analysis Clusters 1 and 2",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#00798c","#d1495b"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.heiht = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Cluster 1", "Cluster 2")    # change legend labels.
)

ggsurvA1Vs2


# = 9.2 fit a cox proportional hazards model and plot ==========================

res.coxCl1VsCl2 <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS1vs2)
res.coxCl1VsCl2
test.phCl1VsCl2 <- cox.zph(res.coxCl1VsCl2)
test.phCl1VsCl2
ggcoxzph(test.phCl1VsCl2)


# = 9.3 Save graphs ============================================================

dev.off()
dpi=100
tiff("SURV 1 vs 2 Cluster PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA1Vs2
dev.off()

# = 10.1 Test Between cluster 2 vs 3 ===========================================

MDataS2A <-MData%>%filter(Cl4==2)
MDataS3A <-MData%>%filter(Cl4==3)
MDataS2vs3<-rbind(MDataS2A,MDataS3A)

fit2Vs3 <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS2vs3)

ggsurvplot(fit2Vs3, data = MDataS2vs3)
ggsurvplot(fit2Vs3, data = MDataS2vs3, censor.shape="|", censor.size = 4)


ggsurvA2Vs3 <- ggsurvplot(
  fit2Vs3,                   # survfit object with calculated statistics.
  data =  MDataS2vs3,        # data used to fit survival curves.
  title = "PROFILE Survival Analysis Clusters 2 and 3",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#d1495b","#edae49"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Cluster 2", "Cluster 3")    # change legend labels.
)

ggsurvA2Vs3


# = 10.2 fit a cox proportional hazards model and plot =========================

res.coxCl2VsCl3 <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS2vs3)
res.coxCl2VsCl3
test.phCl2VsCl3 <- cox.zph(res.coxCl2VsCl3)
test.phCl2VsCl3
ggcoxzph(test.phCl2VsCl3)

# = 10.3 Save graphs ===========================================================

dev.off()
dpi=100
tiff("SURV 2 vs 3 Cluster PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA2Vs3
dev.off()

# = 11.1 Test Between cluster 1 vs 3 ===========================================

MDataS1A <-MData%>%filter(Cl4==1)
MDataS3A <-MData%>%filter(Cl4==3)
MDataS1vs3<-rbind(MDataS1A,MDataS3A)

fit1Vs3 <- survfit(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS1vs3)

ggsurvplot(fit1Vs3, data = MDataS1vs3)
ggsurvplot(fit1Vs3, data = MDataS1vs3, censor.shape="|", censor.size = 4)


ggsurvA1Vs3 <- ggsurvplot(
  fit1Vs3,                   # survfit object with calculated statistics.
  data =  MDataS1vs3,        # data used to fit survival curves.
  title = "PROFILE Survival Analysis Clusters 1 and 3",
  ggtheme=custom_theme(),
  risk.table = TRUE,    # show risk table.
  pval = TRUE,          # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#00798c","#edae49"),
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Cluster 1", "Cluster 3")    # change legend labels.
)

ggsurvA1Vs3

# = 11.2 fit a cox proportional hazards model and plot =========================

res.coxCl1VsCl3 <- coxph(Surv(DeadFromIn5Y,Diedin5Years) ~ Cl4, data = MDataS1vs3)
res.coxCl1VsCl3
test.phCl1VsCl3 <- cox.zph(res.coxCl1VsCl3)
test.phCl1VsCl3
ggcoxzph(test.phCl1VsCl3)

# = 13.3 Save graphs ===========================================================

dev.off()
dpi=100
tiff("SURV 1 vs 3 Cluster PROFILE 220621.tiff",  res=dpi, height=10*dpi, width=15*dpi)
ggsurvA1Vs3
dev.off()


