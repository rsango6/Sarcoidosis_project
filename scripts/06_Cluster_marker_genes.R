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

# [1] "Fcer1g"   "Arhgap15" "Cst3"     "Ctss"     "Lyn"      "Gimap6"   "Tmsb10"   "Ifitm3"   "Ets1"    
# [10] "Gpx1"     "Rpsa"     "Actg1"    "Nr4a1"    "H2-Ab1"   "Ptprcap"  "Malat1"   "Cytip"    "Mbnl1"   
# [19] "Klf4"     "Cd52"     "Tyrobp"   "Mir142hg" "Rpl10a"   "H2-Aa"    "Ablim1"   "S100a10"  "Bach2"   
# [28] "Dennd4a"  "Psap"     "Srgn"     "Fyb"      "Rpl3"     "Pcbp2"    "H2-Eb1"  


#Finding cluster-specific markers
#By default, findMarkers() will give a high ranking to genes that are 
#differentially expressed in any pairwise comparison. This means that
#a gene only needs a very low p-value in a single pairwise comparison
#to achieve a low Top value.


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

# pick clusters: 6, 7, 10, 18, 19, 21, 22 (based on )

c10 = markers_up_all[[10]]
c22 = markers_up_all[[22]]
c18 = markers_up_all[[18]]
c7 = markers_up_all[[7]]
c2 = markers_up_all[[2]]
c21 = markers_up_all[[21]]
c11 = markers_up_all[[11]]
c20 = markers_up_all[[20]]
c13 = markers_up_all[[13]]
c16 = markers_up_all[[16]]
c14 = markers_up_all[[14]]



cluster_list = list()

process_clusters = function(c) {
  c = c %>%
    as.data.frame() %>%
    dplyr::arrange(desc(summary.logFC)) %>%
    dplyr::select(c(1:3)) %>%
    dplyr::mutate(Gene = rownames(c)) #%>%
    #dplyr::top_n(20)
  
}

for (i in 1:length(markers_up_all)) {
  
  cluster_list[[i]] = markers_up_all[[i]]
  
  cluster_list[[i]] = process_clusters(cluster_list[[i]])
  
  names(cluster_list)[i] = paste0("Cluster", i)
  
  #setwd("/scratch/cube/sango/sarcoidosis_project/results")
  #write.csv(cluster_list[[i]], paste0(names(cluster_list)[i], ".csv"))
  
}

c10 = process_clusters(c10) # Lef1, Tcf7
c22 = process_clusters(c22) #Fabp5, Fth1, Ftl1, Lgals3, Cstb, S100a1, Atp6v0c, Gpnmb
c18 = process_clusters(c18) #Ccl5, Gzmk
c7 = process_clusters(c7)  #Ctsc, Rbpj, Mrc1, C1qa, Psd3, Dab2, Pf4, (Pid1?)
c2 = process_clusters(c2) #Ccl4, Ccl3, Il6, Mcpt8
c21 = process_clusters(c21) #Tpsb2, Cpa3, Mcpt4, Cma1
c11 = process_clusters(c11) #S100a9, S100a8, Cd44, Mcl1, Retnlg, Il1b, Srgn, Slpi, Csf3r, Clec4d
c20 = process_clusters(c20) #Fabp4, Ly6c1, Flt1, Cavin2, Cav1
c13 = process_clusters(c13) #Gzma, Klrb1c, Ncr1
c16 = process_clusters(c16) #Lyz2, Fn1, Plcb1
c14 = process_clusters(c14) #Dcn, Gsn, Mgp, Igfbp7, C3, Serping1, Fbln1, Col3a1, Clec3b, Lum

top20genes_postSct = read.csv(paste0(path, "/results/03_transformed_top_20_genes.csv"))

quick_tsne = function(gene) {
  plotTSNE(sce, 
           by_exprs_values = "sctrans_norm",
           colour_by = gene,
           text_by = "louvain") +
    ggtitle(gene)
}

top20genes_tsne = list()

for (i in 1:length(top20genes_postSct$Symbol)) {
  top20genes_tsne[[i]] = quick_tsne(top20genes_postSct$Symbol[i])
}
names(top20genes_tsne) = top20genes_postSct$Symbol

for (i in names(top20genes_tsne)) {
  setwd("/scratch/cube/sango/sarcoidosis_project/Plots/06_ClusterMarkerGenes")
  
  ggsave(filename = paste("TSNE", i, ".pdf", sep = "_"), plot = top20genes_tsne[[i]],
         device = "pdf", width = 8, height = 8, units = "in")
  
}

quick_tsne('Fabp5')
  
  
  
  
  

                                                                