# Part 6: Cluster marker genes

# Load packages ----
library(scater)
library(scran)
library(pheatmap)
library(tidyverse) # always load tidyverse after other packages

path = "/scratch/cube/sango/sarcoidosis_project"
# data corrected using batch integration with mutual nearest neighbours
sce <- readRDS(paste0(path, "/results/sarcoidosis_postSct_clust.Rds"))
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)

genes = as.data.frame(rownames(rowData(sce)))

reducedDimNames(sce)[3] = 'UMAP'
plotUMAP(sce, 
         colour_by = "louvain", 
         text_by = "louvain")

#Our objective is to identify genes that distinguish between these clusters. 
#For example genes such as the “CST3” gene, which is a known monocyte marker:

plotTSNE(sce, 
         colour_by = "Cst3", 
         text_by = "louvain", 
         by_exprs_values = "sctrans_norm")


# Marker gene identification ----

# identify marker genes based on mean expression differences
# default options do not need to be specified, but shown here for illustration
markers_default <- findMarkers(
  sce, 
  groups = factor(sce$louvain), # clusters to compare
  block = sce$Phenotype,    # covariates in statistical model
  test.type = "t",   # t-test (default)
  direction = "any", # test for either higher or lower expression (default)
  lfc = 0, # null hypothesis log-fold-change = 0 (default)
  pval.type = "any" # ranking of p-values based on any comparison (default)
)

# returns a list of length equal to the number of clusters
markers_default

# check the result of a particular cluster
markers_default[[22]]

# extract results for one of the clusters
c22_markers_default <- markers_default[[22]]
c22_markers_default[1:10, c(1, 5:25)]

#We can then use the Top field to identify a set of genes that
#is guaranteed to distinguish cluster 8 from any other cluster.

# identify set of genes in the top 3 p-value ranking
c22_markers_default[c22_markers_default$Top <= 3, ]

#Cst3 and Tyrobp are here in cluster 8: suggest these are monocytes
rownames(c8_markers_default[c8_markers_default$Top <= 10, ])
#As we suspected, cluster 8 is likely to contain monocytes, based 
#on the expression of CST3 and TYROBP, for example Each DataFrame 
#also contains several other statistics that may be of interest. The
#summary.logFC field provides a convenient summary of the direction and
#effect size for each gene, and is defined here as the log-fold change
#from the comparison with the lowest p-value. The p.value field contains
#the combined p-value that is obtained by applying Simes’ method to the
#pairwise p-values for each gene and represents the evidence against the
#joint null hypothesis, i.e., that the gene is not DE between cluster 8 
#and any other cluster. Examination of these statistics permits a quick
#evaluation of the suitability of a candidate marker; if both of these
#metrics are poor (small log-fold change, large p-value), the gene
#can most likely be dismissed.


plotTSNE(sce, 
         colour_by = "Cst3", 
         text_by = "louvain", 
         by_exprs_values = "sctrans_norm")

# visualise the expression of the gene on the uncorrected values
plotExpression(sce, 
               features = "Cst3", 
               x = "louvain")
plotExpression(sce, 
               features = "Hba-a1", #not a marker bc uninteresting LFC
               x = "louvain")


# modify the previous call to findMarkers to focus on genes that are up-regulated
markers_up <- findMarkers(
  sce, 
  groups = factor(sce$louvain), # clusters to compare
  block = sce$Phenotype,    # covariates in statistical model
  test.type = "t",   # t-test (default)
  direction = "up", # test for up-regulated genes only
  lfc = 0, # null hypothesis log-fold-change = 0 (default)
  pval.type = "any" # ranking of p-values based on any comparison (default)
)

c8_markers_up <- markers_up[[8]]
rownames(c8_markers_up)[c8_markers_up$Top <= 3]

# [1] "Fcer1g"   "Arhgap15" "Cst3"     "Ctss"     "Lyn"      "Gimap6"   "Tmsb10"   "Ifitm3"   "Ets1"    
# [10] "Gpx1"     "Rpsa"     "Actg1"    "Nr4a1"    "H2-Ab1"   "Ptprcap"  "Malat1"   "Cytip"    "Mbnl1"   
# [19] "Klf4"     "Cd52"     "Tyrobp"   "Mir142hg" "Rpl10a"   "H2-Aa"    "Ablim1"   "S100a10"  "Bach2"   
# [28] "Dennd4a"  "Psap"     "Srgn"     "Fyb"      "Rpl3"     "Pcbp2"    "H2-Eb1"  

c8_markers_up["Hba-a1", ]


#These two settings yield a more focused set of candidate 
#marker genes that are up-regulated in cluster 8.


# testing for the alternative hypothesis that LFC > 1
markers_up_lfc1 <- findMarkers(
  sce, 
  groups = factor(sce$louvain), # clusters to compare
  block = sce$Phenotype,    # covariates in statistical model
  test.type = "t",   # t-test (default)
  direction = "up", # test for up-regulated genes only
  lfc = 1, # null hypothesis log-fold-change = 1
  pval.type = "any" # ranking of p-values based on any comparison (default)
)


# fetching top markers for cluster 8
c8_markers_up_lfc1 <- markers_up_lfc1[[8]]
c8_markers_up_lfc1[c8_markers_up_lfc1$Top <= 3, ]

#Finding cluster-specific markers
#By default, findMarkers() will give a high ranking to genes that are 
#differentially expressed in any pairwise comparison. This means that
#a gene only needs a very low p-value in a single pairwise comparison
#to achieve a low Top value.


#Take the Cst3 gene, which does seem indeed more highly expressed
#in cluster 8, but only compared to a couple of other clusters (7 and maybe 9).

#While this gene is partially contributing to the distinction
#between clusters, it is not the most diagnostic gene for cluster
#8 (if that is what we were interested in).

#A more stringent approach would only consider genes that are
#differentially expressed in all pairwise comparisons involving
#the cluster of interest. To achieve this, we set pval.type="all"
#in findMarkers() to use an intersection-union test (Berger and Hsu
#1996) where the combined p-value for each gene is the maximum of the
#p-values from all pairwise comparisons. A gene will only achieve a low
#combined p-value if it is strongly DE in all comparisons to other clusters.

# ranking based on the maximum p-value across all pairwise comparisons
markers_up_all <- findMarkers(
  sce, 
  groups = factor(sce$louvain), # clusters to compare
  block = sce$Phenotype,    # covariates in statistical model
  test.type = "t",   # t-test (default)
  direction = "up", # test for up-regulated genes only
  lfc = 0, # null hypothesis log-fold-change = 1
  pval.type = "all" # ranking of p-values based on all comparisons
)


#In this case, the resulting tables do not include a 
#Top column any more, as the ranking is simply based
#on the maximum p-value observed across all comparisons.
#The table is now simply ranked from low-to-high p-value.

# fetching top markers for cluster 8
c8_markers_up_all <- markers_up_all[[8]]
c8_markers_up_all[1:10, ]

plotExpression(sce, 
               features = "Cst3",
               
               x = "louvain")

plotTSNE(sce, colour_by = "Fabp5", 
         text_by = "louvain", 
         by_exprs_values = "sctrans_norm")

plotTSNE(sce, colour_by = "S100a9", 
         text_by = "louvain", 
         by_exprs_values = "sctrans_norm")

# select some top genes for cluster 8
c8_top10 <- c8_markers_up_all[c8_markers_up_all$p.value <= 1e-2, ]


# Alternative testing strategies ----

# Wilcoxon rank-sum test
markers_wilcox_up <- findMarkers(
  sce, 
  groups = sce$louvain, # clusters to compare
  block = sce$Phenotype,    # covariates in statistical model
  test.type = "wilcox",   # t-test (default)
  direction = "up",
  pval.type = "all"
)

rownames(markers_wilcox_up[[22]])[1:50]

plotTSNE(sce, colour_by = "Lgals3", 
         text_by = "louvain", 
         by_exprs_values = "sctrans_norm")



                                                                