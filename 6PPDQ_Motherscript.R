# TO DO
# • Run univariate stats on the WHOLE dataset

# Preamble ====
set.seed(123)

library(pacman)
p_load(
  bench,
  Cairo,
  cowplot,
  data.table,
  doParallel,
  foreach,
  HDANOVA, 
  imputeLCMD, 
  janitor,
  mixOmics,
  multcompView,
  parallel,
  readr,
  rstudioapi, 
  tidyverse, 
  vegan,
  WRS2)

setwd(dirname(getActiveDocumentContext()$path))

pal <- c(
  Control = "black",
  Solvent = "#bebebe",
  Low = "#278192",
  High = "#3f3f7b")

# Source required scripts for analyses ====
# Loads all of the functions I've written for this pipeline
source(file = "Metabolomics_Funcs.R")

# Gets the data and wrangles it into the right format
source(file = "Import_6PPDQ.R")

# Caclulates unviariate statistics for tadpoles' physical traits
source(file = "Univariate_Tadpole_Physical.R")

# Calculates relative standard deviations and stores in a nice list
# source(file = "RSD_6PPDQ.R")

# Pretreats QToF metabolite data
source(file = "QToF_Metabo_6PPDQ.R")

# Pretreats QToF lipid data
source(file = "QToF_Lipids_6PPDQ.R")

# Pretreats central carbon metabolite data
source(file = "CCM_6PPDQ.R")

# Merges all data into one dataframe for further analyses. Final pretreatments:
#     • 80% rule
#     • Median normalise
#     • Pareto scale
source(file = "Merge_6PPDQ.R")

# Calculates the means of each metabolite by tank to prevent pseudoreplication:
source(file = "Tank_Means_6PPDQ.R")

# Univariate statistical analyses
source(file = "Univariate_6PPDQ.R")

# Univariate statistics for control-solvent comparisons to investigate solvent
# effects.
source(file = "Univariate_Control_vs_Solvent_6PPDQ.R")

# Univariate statistics after Gosner Stage stratification to investigate
# developmental effects.
source(file = "Univariate_Gosner_Stratified_6PPDQ.R")

# PERMANOVA model to see if Gosner stage had any effect
source(file = "Gosner_PERMANOVA.R")

# Main multivariate statistical analyses
source(file = "Multivariate_6PPDQ.R")

# Multivariate Control-Solvent contrast analyses
source(file = "Multivariate_Control_vs_Solvent_6PPDQ.R")

# Plots for individual metabolites
source(file = "Plot_Univariate_6PPDQ.R")

# Perform a PERMANOVA to test if Gosner Stage has an effect

# Print debugging
print("6PPDQ_Motherscript.R ran all modules successfully!")
