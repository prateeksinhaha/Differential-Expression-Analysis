## INSTALLATION OF PACKAGES ##

# All objects from the current workspace (R memory) were removed.

rm(list=ls()) 

# Working Directory was set.

setwd("C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/")

# The following packages were installed into R.

install.packages("BiocManager")
install.packages("forcats")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("writexl")
install.packages("rio")
install.packages("magrittr")
install.packages("utils")

library(BiocManager)

BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("edgeR")
