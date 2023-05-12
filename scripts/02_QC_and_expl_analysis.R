
#Part 1: cellranger pipeline was done in UNIX. now comes analysis in R.

#Part 2: Quality control and exploratory analysis 
#We will now check the quality of the data further:



# 1. Mapping quality
# 2. Cell counts
# 3. Distribution of keys quality metrics

#We will then:

# 1. Filter genes with very low expression
# 2. Identify low-quality cells
# 3. Filter and/or mark low quality cells

library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)

path = '/scratch/cube/sango/sarcoidosis_project'
setwd(path)


#loading a single sample
