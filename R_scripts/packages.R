###############################################################################
################ Install packages for the RNA seq analyses ###################


# clear the environment 
rm(list = ls())


# set teh working directory 
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/UOB documents_/RNA Seq/CSHL")


#install.packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db", force = TRUE)
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("dplyr", dependencies = TRUE)
install.packages("Seurat")
install.packages("ggrepel")


# load libraries 
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)



