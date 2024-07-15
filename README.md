# Analysis of Forced Vital Capacity (FVC) trajectories in Idiopathic Pulmonary Fibrosis (IPF) identifies four distinct clusters of disease behaviour.
This repository included the codes associated with the Analysis of Forced Vital Capacity (FVC) trajectories in Idiopathic Pulmonary Fibrosis (IPF)

## Overview of the project codes and workflow included in the repository
### Script 01
1) Database selection of variables
2) Time frame formatting 
3) Visit selection baseline, 90 days or 180 days
4) Estimation of time between spirometric visits
5) Anthropometric variables

### Script 02
Missing visits simulation

### Script 03
Imputation test sample (random forest)
Error imputation calculation: 
1) normalised root-mean-squared deviation
2) % Error (no use in the final publication)
3) Graphs (ggplot)

### Script 04
Monte Carlo Markov Chain Monte-Carlo simulation (Naive)

### Script 05 
Self Organizing Maps

### Script 06A
Elbow analysis

### Script 06B
Jaccard indexes

### Script 07 
Basic statistical analysis

### Script 08
Survival analysis of the clusters generated for all profile

These scripts were executed and designed using R 4.0.3. 
Functions and arguments used in R 4.0.3 might be deprecated or behave differently in newer versions. 
Review the R changelog and package news for deprecated functions or changed behaviours.

** This repository does not include the data (GPDR - protected) **
