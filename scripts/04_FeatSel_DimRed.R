
#Part 4: Feature Selection and Dimensionality Reduction
library(scater) 
library(scran)
library(PCAtools)
library(tidyverse) # always load tidyverse after other packages

# read data ----
sce <- readRDS(paste0(path, "/results/sarcoidosis_postSct.Rds"))
sce

#To make some of our plots later on easier to interpret, we 
#will replace the rownames of the object (containing Ensembl
#gene IDs) with the common gene names. Sometimes it happens
#that there is no common gene name, or different genes have
#the same common name. In both cases, we should instead keep
#the Ensembl ID, to guarantee unique gene names in our object.
#A safe way to do this, is to use the uniquifyFeatureNames() function:


# use common gene names instead of Ensembl gene IDs
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)

#3.1 Quantifying per-gene variation
#Assuming that, for most genes, the observed variance across
#cells is due to technical noise, we can assess technical
#variation by fitting a trend line between the mean-variance
#relationship across all genes. Genes that substantially
#deviate from this trend may then be considered as highly-variable,
#i.e. capturing biologically interesting variation.

# Mean-variance model ----

# fit the model; the output is a DataFrame object
gene_var <- modelGeneVar(sce)

# plot the variance against mean of expression with fitted trend line
gene_var %>% 
  # convert to tibble for ggplot
  as_tibble() %>% 
  # make the plot
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")


#3.2 Selecting highly variable genes
#Once we have quantified the per-gene variation, the next
#step is to select the subset of HVGs to use in downstream
#analyses. A larger subset will reduce the risk of discarding
#interesting biological signal by retaining more potentially
#relevant genes, at the cost of increasing noise from irrelevant
#genes that might obscure said signal. It is difficult to
#determine the optimal trade-off for any given application
#as noise in one context may be useful signal in 
#another. Commonly applied strategies are:

#1. take top X genes with largest (biological) variation
#2. based on significance
#3. keeping all genes above the trend
#4. selecting a priori genes of interest

#In our example, we will define ‘HVGs’ as the top 10% of genes
#with the highest biological component. This is a fairly arbitrary
#choice. A common practice is to pick an arbitrary threshold 
#(either based on number of proportion) and proceed with the
#rest of the analysis, with the intention of testing other choices
#later, rather than spending much time worrying about
#obtaining the “optimal” value.

# identify the top 10% most variable genes
hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs) # check how many genes we have: 1067

hvgs[1:10]

# "Cd74"   "Lyz2"   "Igkc"   "S100a9" "Apoe"   "Ccl5"   "H2-Aa"  "S100a8" "Ebf1"   "H2-Eb1"

#The result is a vector of gene IDs ordered by their biological variance
#(i.e. highest deviation from the trend line shown above). We can use
#this with functions that accept a list of genes as option to restrict
#their analysis to that subset of genes (e.g. when we do PCA later on).

#For example, we may want to visualise the expression of the top most-variable
#genes determined in this way. We can do this with a violin plot for each gene,
#using the plotExpression() function:

plotExpression(sce, features = hvgs[1:20], point_alpha = 0.05, jitter = "jitter")

#4.1.1 Running PCA
# PCA ----
#The runPCA() function can be used to run PCA on a SCE object, 
#and returns an updated version of that object with the PCA 
#result added to the reducedDim slot. Importantly, we can also
#restrict the PCA to use only some of the features (rows) of the
#object, which in this case we do by using the highly variable
#genes we identified earlier.

sce <- runPCA(sce, subset_row = hvgs)
sce

reducedDim(sce, "PCA")[1:10, 1:5]

#By default, runPCA() returns the first 50 PCs, but you can change this number by specifying the ncomponents option.
## extract variance explained
pca_pct_variance <- data.frame(variance = attr(reducedDim(sce, "PCA"), "percentVar"))
pca_pct_variance$PC <- 1:nrow(pca_pct_variance)
# visualise percentage variance explained by PCs (scree plot)
pca_pct_variance %>% 
  ggplot(aes(PC, variance)) +
  geom_col() +
  labs(y = "Variance explained (%)")

# visualise PC plot
plotReducedDim(sce, dimred = "PCA", colour_by = "SampleName")

# visualise multiple PCs
plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "SampleName")

# more custom visualisations with ggcells (e.g. add facets)
ggcells(sce, aes(x = PCA.1, y = PCA.2, colour = SampleName)) +
  geom_point(size = 0.5) +
  facet_wrap(~ SampleName) +
  labs(x = "PC1", y = "PC2", colour = "SampleName")


#4.1.2 PCA Diagnostics


ah <- AnnotationHub()
ens.mm.102 <- query(ah, c("Mus musculus", "EnsDb", 102))[[1]]

genes <- rowData(sce)$ID
gene_annot <- AnnotationDbi::select(ens.mm.102, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME")) %>%
  set_names(c("ID", "Chromosome"))

rowData(sce)$Chromosome = gene_annot$Chromosome[match(rowData(sce)$ID, gene_annot$ID)]

is.mito <- which(rowData(sce)$Chromosome=="MT")

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)

# extract correlations between different variables and our PC scores
explan_pcs <- getExplanatoryPCs(sce,
                                variables = c(
                                  "sum",
                                  "detected",
                                  "Phenotype",
                                  "SampleName",
                                  "subsets_Mito_percent"
                                  ))
        

plotExplanatoryPCs(explan_pcs/100)

# distribution of correlations between each gene's expression and our variables of interest
plotExplanatoryVariables(sce,
                         variables = c(
                           "sum",
                           "detected",
                           "Phenotype",
                           "SampleName",
                           "subsets_Mito_percent"
                         ))
#This analysis indicates that individual and detection rate have the highest 
#explanatory power for many genes, and we don’t see technical covariates
#having as high correlations. If that were the case, we might need to 
#repeat the normalization step while conditioning out for these covariates,
#or we would include them in downstream analysis.

#4.1.3 Chosing the number of PCs

#how many PCs have percent variation contribution bigger than 1?
table(pca_pct_variance$variance > 1)
#FALSE  TRUE 
#40    10 


#4.1.3.1 Elbow point
# identify elbow point from explained variances
chosen_elbow <- findElbowPoint(pca_pct_variance$variance)
chosen_elbow #7

# scree plot (PC vs variance) with elbow highlighted
pca_pct_variance %>% 
  ggplot(aes(PC, variance)) +
  geom_point() +
  geom_vline(xintercept = chosen_elbow)

#4.1.3.2 Denoising PCA

#The assumption of this method is that the biology drives most
#of the variance and hence should be captured by the first few PCs,
#while technical noise affects each gene independently, hence it
#should be captured by later PCs. Therefore, our aim in this approach
#is to find the minimum number of PCs that explains more variance
#than the total technical variance across genes (estimated
#from our mean-variance trend).

# run denoise PCA step
sce <- denoisePCA(sce, technical = gene_var)
# check dimensions of the "denoised" PCA
ncol(reducedDim(sce, "PCA")) #5


#4.2 t-SNE: t-Distributed Stochastic Neighbor Embedding
# we run t-SNE based on the PCA we ran previously

set.seed(123) # set a random seed to ensure reproducibility
sce <- runTSNE(sce, 
               name = "TSNE_perplex50",
               perplexity = 50, 
               dimred = "PCA")

# Make a custom visualisation using ggcells
ggcells(sce, aes(x = TSNE_perplex50.1, y = TSNE_perplex50.2, 
                 colour = SampleName)) +
  geom_point()


## run the UMAP with 50 neighbours
set.seed(123) # set seed for reproducibility
sce <- runUMAP(sce, 
               name = "UMAP_neighbors50",
               dimred = "PCA")

ggcells(sce, aes(x = UMAP_neighbors50.1, y = UMAP_neighbors50.2, 
                 colour = SampleName)) +
  geom_point()

saveRDS(sce, paste0(path, "/results/sarcoidosis_postDeconv_dimRed.Rds"))

