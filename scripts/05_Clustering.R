#Part 5: Clustering

library(tidyverse) # data wrangling and plotting (ggplot2)
library(scater) # scRnaSeq QC
library(scran) # scRnaSeq normalisation
library(bluster) # scRnaSeq clustering
library(dynamicTreeCut) # tree cutting in clustering
library(cluster) # for silhouette
library(igraph) # for graph-based clustering and plotting networks 
library(reticulate) # for graph-based clustering with leiden
library(leiden) # for community detection (clustering)
library(pheatmap) # for heatmap plotting
library(patchwork) # to combine plots


# Read mnn.out object in:
sce = readRDS(paste0(path, "/results/sarcoidosis_postSct_dimRed.Rds"))

#1.2.1.1 Clustering

#Here we will apply hierarchical clustering on the Euclidean distances
#between cells, using the Ward D2 criterion to minimize the total variance within each cluster.

# get PCs
pcs <- reducedDim(sce, "PCA")
# compute distance
pc.dist <- dist(pcs)
# derive tree
hc.tree <- hclust(pc.dist, method="ward.D2")
hcd <- as.dendrogram(hc.tree)

#The dendrogram below shows each cell as a leaf.

plot(hcd, type = "rectangle", ylab = "Height", leaflab = "none")


# identify clusters by cutting branches, requesting a minimum cluster size of 20 cells.
hc.clusters <- unname(cutreeDynamic(hc.tree,
                                    distM = as.matrix(pc.dist),
                                    minClusterSize = 20,
                                    verbose = 0))
#Cell counts for each cluster (rows) and each sample group (columns):
# cells per sample group
table(hc.clusters, sce$Phenotype)

# hc.clusters Control Sarcoidosis
#     1     4712        3849
#     2     3194        3232
#     3      952        1639
#     4     1684         816
#     5     1710         707
#     6      257         949
#     7      702         452
#     8      161         868
#     9      401         233
#     10     258         337


#Cell counts for each cluster (rows) and each sample (columns):
  
# per sample
table(hc.clusters, sce$SampleName)

# hc.clusters FP19 FP21 FP24 FP25 FP27 FP28
#         1  1823 1527  957 2077 1365  812
#         2  1299 1128 1334  877  770 1018
#         3   393  841  401  160  397  399
#         4   799  396  130  293  290  592
#         5   278  385  208  102  114 1330
#         6   130  487  369   45   93   82
#         7   378  189  144  174  119  150
#         8     6  249  248   65  371   90
#         9    19   29  151  282   53  100
#         10   34   26  257   93   54  131

#This data (cell number) may also be shown as heatmap.

tmpTab <- table(hc.clusters, sce$SampleName)
rownames(tmpTab) = paste("hc", rownames(tmpTab), sep = "_")

# columns annotation with cell name:
mat_col <- colData(sce) %>% data.frame() %>% dplyr::select(SampleName, Phenotype) %>% dplyr::distinct()

rownames(mat_col) <- mat_col$SampleName
mat_col$SampleName <- NULL

# Prepare colours for clusters:
colourCount = length(unique(mat_col$Phenotype))
getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

mat_colors <- list(Phenotype = getPalette(colourCount))
names(mat_colors$Phenotype) <- unique(mat_col$Phenotype)

#Heatmap, with samples ordered as in sample sheet:
  
  # without column clustering
  pheatmap(tmpTab,
           #cluster_cols = FALSE,
           #cluster_rows = FALSE,
           annotation_col = mat_col,
           annotation_colors = mat_colors)
  
  
  
  #If clusters mostly include cells from one sample or the other,
  #it suggests that the samples differ, and/or the presence of batch effect.
  
  #Let us show cluster assignments on the t-SNE, with cells shaped
  #by cell type, colored by cluster and sized total UMI counts (‘sum’).
  
  # store cluster assignment in SCE object:
  sce$cluster <- factor(hc.clusters)
  # make, store and show TSNE plot:
  g <- plotTSNE(sce, colour_by = "cluster", shape_by="Phenotype")
  g

  
  # split by sample and show:
  g <- g + facet_wrap(. ~ sce$Phenotype)
  g

  
  
  #We can also check the contribution of each sample by painting cells by ‘SampleName’:
    
  g <- plotTSNE(sce, colour_by = "SampleName", shape_by="Phenotype") +
    facet_wrap(. ~ sce$Phenotype)
  g

  
  
  #1.2.1.2 Separatedness
  
  #The congruence of clusters may be assessed by computing the
  #sillhouette for each cell. The larger the value the closer
  #the cell to cells in its cluster than to cells in other clusters.
  #Cells closer to cells in other clusters have a negative value.
  #Good cluster separation is indicated by clusters whose cells 
  #have large silhouette values.
  
  #We first compute silhouette values.

  sil <- silhouette(hc.clusters, dist = pc.dist)
  
  
  #We then plot silhouettes with one color per cluster and cells
  #with a negative silhouette with the color of their closest
  #cluster. We also add the average silhouette for each cluster and all cells.
  
  # prepare colours:
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
  sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
  # plot:
  plot(sil,
       main = paste(length(unique(hc.clusters)), "clusters"),
       border = sil.cols,
       col = sil.cols,
       do.col.sort = FALSE) 

  
  
  #1.2.2 k-means
  # define clusters with kmeans()
  # because results depend on the initial cluster centers,
  # it is usually best to try several times,
  # by setting 'nstart' to say 20,
  # kmeans() will then retain the run with the lowest within cluster variation.
  kclust <- kmeans(pcs, centers = 10, nstart = 20) 
  
  # compute silhouette
  sil <- silhouette(kclust$cluster, dist(pcs))
  
  # plot silhouette
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
  sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
  plot(sil, main = paste(length(unique(kclust$cluster)), "clusters"), 
       border=sil.cols, col=sil.cols, do.col.sort=FALSE)
  
  #tSNE
  sce$k10 <- factor(kclust$cluster)
  # make, store and show TSNE plot:
  g <- plotTSNE(sce, colour_by = "k10", shape_by="Phenotype")
  g = g + facet_wrap(~ sce$Phenotype)
  g
  
  
  
  
#To find the most appropriate number of clusters, one performs the analysis for a series
  #of values for k and for each k value computes a ‘measure of fit’ of the clusters
  #defined. Several metrics exist. We will look at:
    
  #1. within-cluster sum-of-squares
  #2. silhouette
  #3. gap statistic
  
  #1.2.2.3 Within-cluster sum-of-squares
  
 # The within-cluster sum-of-squares is the sum of the squared deviations
  #from each observation and the cluster centroid. This metric measures
  #the variability of observations within a cluster. It decreases as k
  #increases, by an amount that decreases with k. Indeed, for low k values
  #increasing k usually improves clustering, as shown by a sharp drop in
  #the within-cluster sum-of-squares. For higher k values, increasing k
  #does not reduce the within-cluster sum-of-squares much: dividing 
  #clusters further only slightly reduces variablility inside clusters.
  #On a plot of the within-cluster sum-of-squares against k, the curve
  #shows two parts: the first with a steep negative slope and the second
  #with a small slope. The ‘elbow’ of the curve indicates the most
  #appropriate number of clusters.
  
  
  library(broom)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # get PCA matrix (i.e. rotated values)
  points <- as_tibble(pcs)
  
  # define clusters for different number of clusters
  # from 1 to 20:
  kclusts <- tibble(k = 1:20) %>%
    dplyr::mutate(
      kclust = map(k, ~kmeans(points, .x)), # define clusters
      tidied = map(kclust, tidy), # 'flatten' kmeans() output into tibble
      glanced = map(kclust, glance), # convert model or other R object to convert to single-row data frame
      augmented = map(kclust, augment, points) # extract per-observation information
    )
  
  # get cluster assignments
  # unnest a list column with unnest(),
  # i.e. make each element of the list its own row.
  clusters <- kclusts %>%
    unnest(tidied)
  
  # get assignments
  assignments <- kclusts %>% 
    unnest(augmented)
  
  table(assignments$k)
  
  # get clustering outcome
  clusterings <- kclusts %>%
    unnest(glanced)
  
  #We now plot the total within-cluster sum-of-squares and decide on k.
  
  ggplot(clusterings, aes(k, tot.withinss)) +
    geom_point() +
    geom_line()
  
  clusterings %>% 
    mutate(tot.withinss.diff = tot.withinss - lag(tot.withinss)) %>%
    arrange(desc(tot.withinss.diff))
  
  
  # get the smallest negative drop
  k_neg <- clusterings %>% 
    mutate(tot.withinss.diff = tot.withinss - lag(tot.withinss)) %>%
    dplyr::filter(tot.withinss.diff < 0) %>%
    slice_max(tot.withinss.diff, n = 1) %>%
    pull(k)
  k_neg <- k_neg -1 
  #k = 15
  # get the first positive diff
  k_pos <- clusterings %>% 
    mutate(tot.withinss.diff = tot.withinss - lag(tot.withinss)) %>%
    dplyr::filter(tot.withinss.diff >= 0) %>%
    slice_min(k, n = 1) %>%
    pull(k)
  k_pos <- k_pos -1
  #k = 17
  
  
  #1.2.2.4 Silhouette
  
  #We now compute the Silhouette values for each set of clusters defined above for the series of values of k.
  
  # compute pairwise distance between cells in the PC space,
  # as it is needed to compute silhouette:
  pcDist <- dist(pcs)
  # compute silhouette for each set of clusters defined above for a series of values of k:
  Ks=sapply(2:20,
            function(i) {
              tmpClu <- as.numeric(kclusts$augmented[[i]]$.cluster)
              #table(tmpClu)
              sil <- silhouette(tmpClu, pcDist)
              summary(sil)$avg.width
            }
  )
  # plot average width against k:
  plot(2:20, Ks,
       xlab="k",
       ylab="av. silhouette",
       type="b",
       pch=19)
  
  
  #1.2.2.5 Gap statistic
  
  #Because the variation quantity decreases as the number of clusters
  #increases, a more reliable approach is to compare the observed variation
  #to that expected from a null distribution. With the gap statistic,
  #the expected variation is computed: by generating several reference
  #data sets using the uniform distribution (no cluster),
  #for each such set and the series of k values, by 
  #clustering and computing the within-cluster variation.
  
  #The gap statistic measures for a given k the difference between
  #the expected and observed variation. The most appropriate number
  #of clusters is that with the higest gap value.
  
  #We will use cluster::clusGap() to compute the gap statistic for k between 1 and 20.
  
  set.seed(123)
  gaps <- cluster::clusGap(
    x = pcs,
    FUNcluster = kmeans,
    K.max = 20,
    nstart = 5, # low for expediency here but should use higher
    B = 10 # low for expediency here but should use higher
  )
  
  ## find the "best" k
  best.k <- cluster::maxSE(gaps$Tab[, "gap"], gaps$Tab[, "SE.sim"])
  # in case the the best k was '1',
  # skip first row of output table:
  gapsTab <- gaps$Tab[2:nrow(gaps$Tab),]
  best.k <- cluster::maxSE(gapsTab[, "gap"], gapsTab[, "SE.sim"])
  best.k
  #The “optimal” k value is 13
  
  
  # copy kbest
  df <- as.data.frame(assignments)
  sce$kmeans.best <- as.numeric(df[df$k == best.k, ".cluster"])
  
  # plot sihouette
  clust.col <- colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(best.k)
  sil <- silhouette(sce$kmeans.best, dist = pc.dist)
  sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
  sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
  plot(sil, main = paste(length(unique(sce$kmeans.best)), "clusters"), 
       border=sil.cols, col=sil.cols, do.col.sort=FALSE)
  
  
  #1.2.3 Graph-based clustering
  
  #We will:
  
  # 1. build the graph,
  # 2. define clusters,
  # 3. check membership across samples,
  # 4. show membership on a t-SNE plot,
  # 5. assess clustering quality.
  
  #We now build the shared nearest neighbour (SNN) graph, using scran’s buildSNNGraph() with:
    
  # 1. the reduced and denoised data set (PCs)
  # 2. the default number of neighbours (k=10)
  # 3. the default type of edge weight (type=“rank”)
  
  snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
  
  
  # subset graph down to 1000 cells
  # https://igraph.org/r/doc/subgraph.html
  
  # Add vertices (nodes. ie cells) annotation
  V(snn.gr)$SampleName <- colData(sce)$SampleName
  V(snn.gr)$Phenotype <- as.character(colData(sce)$Phenotype)
  
  # pick 1000 nodes randomly
  edgesToGet <- sample(nrow(snn.gr[]), 7000)
  
  # subset graph for these 1000 ramdomly chosen nodes
  snn.gr.subset <- subgraph(snn.gr, edgesToGet)
  g <- snn.gr.subset
  
  # set colors for clusters
  if(length(unique(V(g)$Phenotype)) <= 2)
  {
    cols <- colorspace::rainbow_hcl(length(unique(V(g)$Phenotype)))
  } else {
    cols <- RColorBrewer::brewer.pal(n = length(unique(V(g)$Phenotype)), name = "RdYlBu")
  }
  names(cols) <- unique(V(g)$Phenotype)
  
  # plot graph
  plot.igraph(
    g, layout = layout_with_fr(g),
    vertex.size = 3, vertex.label = NA,
    vertex.color = cols[V(g)$Phenotype],
    frame.color = cols[V(g)$Phenotype],
    main = "default parameters"
  )
  
  # add legend
  legend('topright',
         legend=names(cols),
         pch=21,
         col=cols, # "#777777",
         pt.bg=cols,
         pt.cex=1, cex=.6, bty="n", ncol=1)
  
  
#1.2.3.2 Walktrap method
  
#The walktrap method relies on short random walks (a few steps) through '
#the network. These walks tend to be ‘trapped’ in highly-connected regions
#of the network. Node similarity is measured based on these walks. Nodes are
#first each assigned their own community. Pairwise distances are computed and
#the two closest communities are grouped. These steps are repeated a given 
#number of times to produce a dendrogram. Hierarchical clustering is then
#applied to the distance matrix. The best partition is that with the highest modularity.
  
# identify clusters with walktrap
# default number of steps: 4
cluster.out <- cluster_walktrap(snn.gr)


# cluster assignments are stored in the membership slot
wt.clusters <- cluster.out$membership

table(wt.clusters)

table(wt.clusters, sce$Phenotype)

tmpTab <- table(wt.clusters, sce$SampleName)
rownames(tmpTab) = paste("4-step", rownames(tmpTab), sep = "_")

# columns annotation with cell name:
mat_col <- colData(sce) %>% data.frame() %>% dplyr::select(SampleName, Phenotype) %>% distinct()
rownames(mat_col) <- mat_col$SampleName
mat_col$SampleName <- NULL

# Prepare colours for clusters:
colourCount = length(unique(mat_col$Phenotype))
getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

mat_colors <- list(Phenotype = getPalette(colourCount))
names(mat_colors$Phenotype) <- unique(mat_col$Phenotype)

pheatmap(tmpTab,
         annotation_col    = mat_col,
         annotation_colors = mat_colors)

#The table below shows the cluster distribution across samples. 
#Most clusters comprise cells from several replicates of a same
#sample type, while few are observed in only one sample.


#With 15-step walks:

cluster.out.s15 <- cluster_walktrap(snn.gr, steps = 15)
wt.clusters.s15 <- cluster.out.s15$membership
table(wt.clusters.s15)

#We now apply the Louvain approach, store its outcome in the SCE object and show cluster sizes.


ig.louvain <- igraph::cluster_louvain(snn.gr)
cl <- ig.louvain$membership
cl <- factor(cl)
# store membership
sce$louvain <- cl
# show cluster sizes:
table(sce$louvain)

sce

V(snn.gr)$SampleName <- colData(sce)$SampleName
V(snn.gr)$SampleGroup <- as.character(colData(sce)$Phenotype)
V(snn.gr)$louvain <- as.character(colData(sce)$louvain)

# once only
edgesToGet <- sample(nrow(snn.gr[]), 10000)

snn.gr.subset <- subgraph(snn.gr, edgesToGet)
g <- snn.gr.subset

tmpLayout <- layout_with_fr(g)

rgb.palette <- colorRampPalette(c("purple","yellow"), space="rgb") # Seurat-like
cols <- rgb.palette(length(unique(V(g)$louvain)))
names(cols) <- as.character(1:length(unique(V(g)$louvain)))

# Plot the graph, color by cluster assignment
igraph::plot.igraph(
  g, layout = tmpLayout,
  vertex.color = cols[V(g)$louvain],
  frame.color = cols[V(g)$louvain],
  vertex.size = 3, vertex.label = NA, main = "Louvain"
)

# add legend
legend('bottomright',
       legend=names(cols),
       pch=21,
       col=cols, # "#777777",
       pt.bg=cols,
       pt.cex=1, cex=.6, bty="n", ncol=2)



fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
# show clusters on TSNE
reducedDimNames(sce)[2] = "TSNE"
p <- plotTSNE(sce, colour_by="louvain") + fontsize
p + facet_wrap(~ sce$Phenotype)



#1.2.3.4 Leiden method


cluToUse <- "louvain"
plotTSNE(sce, colour_by = cluToUse,
         text_by = "louvain") +
  ggtitle("louvain clusters")


rownames(sce) <- scater::uniquifyFeatureNames(
  rowData(sce)$ID,
  rowData(sce)$Symbol)
p1 <- plotTSNE(sce, by_exprs_values = "sctrans_norm", colour_by = "Ms4a1") +
  ggtitle("B cells: Ms4a1")
p2 <- plotTSNE(sce, by_exprs_values = "sctrans_norm", colour_by = "Cd79a") +
  ggtitle("B cells: Cd79a")
gridExtra::grid.arrange(p1, p2, ncol=2)


rm(p1, p2)

quick_tsne = function(gene) {
  plotTSNE(sce, by_exprs_values = "sctrans_norm", colour_by = gene) +
    ggtitle(gene)
}


top20genes_postSct = read.csv(paste0(path, "/results/03_transformed_top_20_genes.csv"))

top20genes_tsne = list()
for (i in 1:length(top20genes_postSct$Symbol)) {
  top20genes_tsne[[i]] = quick_tsne(top20genes_postSct$Symbol[i])
  names(top20genes_tsne) = top20genes_postSct$Symbol
}

p1 <- plotExpression(sce,
                     exprs_values = "sctrans_norm",
                     x=cluToUse,
                     colour_by=cluToUse,
                     features= "Ms4a1") +
  ggtitle("B cells")
p2 <- plotExpression(sce,
                     exprs_values = "sctrans_norm",
                     x=cluToUse,
                     colour_by=cluToUse,
                     features= "Cd79a") +
  ggtitle("B cells")
gridExtra::grid.arrange(p1, p2, ncol=2)

markGenesFull <- list(
  "Naive CD4+ T cells"=c("Il7r", "Ccr7"),
  "Memory CD4+ T cells"=c("Il7r", "S100a4"),
  "B cells"=c("Ms4a1", "Cd79a"),
  "CD8+ T cells"=c("Cd8a"),
  "NK cells"=c("Nkg7"),
  "CD14+ Monocytes"=c("Cd14", "Lyz1"),
  "FCGR3A+ Monocytes"=c("Ms4a7"),
  "Dendritic Cells"=c("Fcer1a", "Cst3"),
  "Platelets"=c("Ppbp")
)
markGenesAvail <- lapply(markGenesFull,
                         function(x){
                           x[x %in% rownames(sce)]
                         })


plotExpression2 <- function(sce, ctx, markGenesAvail) {
  if(length(markGenesAvail[[ctx]])>0) {
    plotExpression(sce, exprs_values = "sctrans_norm",
                   x=cluToUse, colour_by=cluToUse,
                   features=markGenesAvail[[ctx]]) + ggtitle(ctx)
  }
}

# Naive CD4+ T
ctx <- "Naive CD4+ T cells"
plotExpression2(sce, ctx, markGenesAvail) 

# Memory CD4+
ctx <- "Memory CD4+ T cells"
plotExpression2(sce, ctx, markGenesAvail)

# B cells
ctx <- "B cells"
plotExpression2(sce, ctx, markGenesAvail)

# CD8+ T
ctx <- "CD8+ T cells"
plotExpression2(sce, ctx, markGenesAvail)

# NK
ctx <- "NK cells"
plotExpression2(sce, ctx, markGenesAvail)

# CD14+ Mono
ctx <- "CD14+ Monocytes"
plotExpression2(sce, ctx, markGenesAvail)

# FCGR3A+ Mono
ctx <- "FCGR3A+ Monocytes"
plotExpression2(sce, ctx, markGenesAvail)

# DC
ctx <- "Dendritic Cells"
plotExpression2(sce, ctx, markGenesAvail)

# Platelet
ctx <- "Platelets"
plotExpression2(sce, ctx, markGenesAvail)

quick_tsne = function(gene, cells) {
  plotTSNE(sce, 
           by_exprs_values = "sctrans_norm",
           colour_by = gene,
           text_by = "louvain") +
    ggtitle(paste(gene, cells, sep = " - "))
}

KnownMarkGenesListTSNE = list()
for (i in 1:length(unlist(markGenesAvail))) {
  KnownMarkGenesListTSNE[[i]] = quick_tsne(unlist(markGenesAvail)[i],
                                           names(unlist(markGenesAvail)[i]))
}

names(KnownMarkGenesListTSNE) = unlist(markGenesAvail)

quick_tsne('Lyz2', 'general Myeloid marker')

saveRDS(sce, file=paste0(path, "/results/sarcoidosis_postSct_clust.Rds"))

