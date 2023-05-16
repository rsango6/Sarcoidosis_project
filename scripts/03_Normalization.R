
#PART 3: Normalisation
#Why normalise?
#Systematic differences in sequencing coverage between libraries 
#occur because of low input material, differences in cDNA capture
#and PCR amplification. Normalisation removes such differences so
#that differences between cells are not technical but biological,
#allowing meaningful comparison of expression profiles between cells.
#Normalisation and batch correction have different aims. Normalisation
#addresses technical differences only, while batch correction considers
#both technical and biological differences.

library(scran)
library(sctransform)
bp.params = MulticoreParam(workers = 7)


#We will load the R object created after QC and check its content (class, dimensions, assays, …)

path = "/scratch/cube/sango/sarcoidosis_project"
sce = readRDS(paste0(path, "/results/Sarcoidosis_filtered_genes.rds"))

sce

dd = colData(sce) %>%
  data.frame() %>%
  dplyr::rename(SampleName=Sample) %>%
  DataFrame()

colData(sce) = dd

#We can also count the number of cells for each sample:

colData(sce) %>%
  # colData() returns a DFrame
  # that we need to convert to a data.frame for parsing
  data.frame() %>%
  # group by some columns only: SampleName, SampleId, SampleGroup
  # (could do with SampleName only but we would miss SampleId, SampleGroup later)
  group_by(SampleName, Phenotype) %>%
  # count cells for each group
  summarise(nbCells=n()) %>%
  # display output table
  DT::datatable(rownames = FALSE,
                options = list(dom="tpl", pageLength = 15))



# have new list of cell barcodes for each sample
vec.bc <- colData(sce) %>%
  data.frame() %>%
  group_by(SampleName) %>%
  pull(Barcode)

#Library size normalization
#For each cell, the library size factor is proportional 
#to the library size such that the average size factor 
#across cell is one.

lib.sf <- librarySizeFactors(sce)
summary(lib.sf)

#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.08331  0.49517  0.74477  1.00000  1.20332 12.92843 
#Size factor distribution: wide range, typical of scRNA-seq data.


dd <- data.frame("log10libSf"=log10(lib.sf))
ggplot(dd, aes(x=log10libSf)) + 
  geom_histogram(bins=50)


#1.4.5 Deconvolution
#The method below increases read counts by pooling cells into groups, computing size 
#factors within each of these groups and scaling them so they are comparable across 
#clusters. This process is repeated many times, changing pools each time to collect
#several size factors for each cell, from which is derived a single value for that cell.

#Cluster cells
#The table below show the number and size of clusters found:
set.seed(100) # clusters with PCA from irlba with approximation
clust <- quickCluster(sce, BPPARAM=bp.params) # slow with all cells.
table(clust)


#Compute size factors
sce <- computePooledFactors(sce,
                            clusters = clust,
                            min.mean = 0.1, # min.mean
                            # A numeric scalar specifying the minimum (library size-adjusted)
                            # average count of genes to be used for normalization.
                            BPPARAM = bp.params)
deconv.sf <- sizeFactors(sce)
summary(deconv.sf)

#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.04279  0.48846  0.74597  1.00000  1.19288 13.59290 


#Plot deconvolution size factors against library size factors:

sce <- addPerFeatureQC(sce, BPPARAM = bp.params) # PATCH
sce <- addPerCellQC(sce, BPPARAM = bp.params)


colData(sce)$cell_sparsity <- 1 - (colData(sce)$detected / nrow(sce))
rowData(sce)$gene_sparsity <- (100 - rowData(sce)$detected) / 100

deconvDf <- data.frame(lib.sf, deconv.sf,
                       "source_name" = sce$Phenotype,
                       "sum" = sce$sum,
                       #"mito_content" = sce$subsets_Mito_percent,
                       "cell_sparsity" = sce$cell_sparsity)

# colour by sample type
sp <- ggplot(deconvDf, aes(x=lib.sf, y=deconv.sf, col=source_name)) +
  geom_point()
sp

# colour by cell sparsity
sp <- ggplot(deconvDf, aes(x=lib.sf, y=deconv.sf, col=cell_sparsity)) +
  geom_point()
sp


#Apply size factors
#For each cell, raw counts for genes are divided by the size factor for that
#cell and log-transformed so downstream analyses focus on genes with strong relative differences. We use scater::logNormCounts().


sce <- logNormCounts(sce) # adds logcounts
# check list of assays stored:
print(assayNames(sce))

saveRDS(sce, paste0(path, "/results/sarcoidosis_postDeconv.Rds"))

###################################################
# this was deconvolution normalisation strategy. ##
# now second type: sctransform                   ##
###################################################


#sctransform: 
#With scaling normalisation a correlation remains between the mean 
#and variation of expression (heteroskedasticity). This affects 
#downstream dimensionality reduction as the few main new dimensions
#are usually correlated with library size. sctransform addresses the
#issue by regressing library size out of raw counts and providing
#residuals to use as normalized and variance-stabilized expression
#values in some downstream analyses, such as dimensionality reduction.
#We will use the sctransform vignette.

#We will first obtain the raw counts matrix:

# keep raw counts in a 'counts' variable:
counts <- counts(sce)
# check the class of the object
# expect a 'dgCMatrix': Compressed, sparse, column-oriented numeric matrices
# the “standard” class for sparse numeric matrices in the Matrix package
print(class(counts))
print(dim(counts)) #18202 27113
## name columns (cells) with barcodes
colnames(counts) <- colData(sce)$Barcode



#Inspect data
#We will now calculate some properties and visually inspect the data.
#Our main interest is in the general trends not in individual outliers.
#Neither genes nor cells that stand out are important at this step;
#we focus on the global trends.

#Derive gene and cell attributes from the UMI matrix
#Gene attributes include for each gene:
  
#mean UMI count across cells
#number of cells where the gene is detected
#variance of UMI counts across cells
#the mean and variance above on the log10 scale
#Cells attributes include for each cell:
  
#total UMI count across genes (library size)
#number of genes detected (with at least 1 UMI)

# gene attributes:
# prepare a data frame named e.g. 'gene_attr' to keep gene attributes, inc:
gene_attr <- data.frame(mean = rowMeans(counts), 
                        detection_rate = rowMeans(counts > 0),
                        var = rowVars(counts))
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
# name rows of the 'gene_attr' data frame:
rownames(gene_attr) <- rownames(counts)

# cell attributes:
cell_attr <- data.frame(n_umi = colSums(counts),
                        n_gene = colSums(counts > 0))
rownames(cell_attr) <- colnames(counts)

dim(gene_attr) #18202     5
head(gene_attr)
dim(cell_attr) #27113     2
head(cell_attr)

#Mean-variance relationship
#For the genes, on the log10 scale we can see that up to a mean 
#UMI count of 0 the variance follows the line through the origin
#with slop one, i.e. variance and mean are roughly equal as expected 
#under a Poisson model. However, genes with a higher average UMI
#count show overdispersion compared to Poisson.

ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(linewidth = 0.3) +
  geom_abline(intercept = 0, slope = 1, color='red')

#Mean-variance relationship
#For the genes, on the log10 scale we can see that up to a mean UMI count of 0 
#the variance follows the line through the origin with slop one, i.e. variance and 
#mean are roughly equal as expected under a Poisson model. However, genes with a 
#higher average UMI count show overdispersion compared to Poisson.

# add the expected detection rate under Poisson model
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x,
                            detection_rate = 1 - dpois(0, lambda = 10^x))
ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_line(data=poisson_model, color='red') +
  theme_gray(base_size = 8)


#Cells attributes
#The plot below show the relationship between the to cell attributes computed:
#library size (n_umi) and number of genes detected (n_gene).

ggplot(cell_attr, aes(n_umi, n_gene)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3)

#Transformation
#Method
#“Based on the observations above, which are not unique to this particular data set,
#we propose to model the expression of each gene as a negative binomial random variable
#with a mean that depends on other variables. Here the other variables can be used to
#model the differences in sequencing depth between cells and are used as independent
#variables in a regression model. In order to avoid overfitting, we will first fit
#model parameters per gene, and then use the relationship between gene mean and
#parameter values to fit parameters, thereby combining information across genes.
#Given the fitted model parameters, we transform each observed UMI count into a
#Pearson residual which can be interpreted as the number of standard deviations
#an observed count was away from the expected mean. If the model accurately describes
#the mean-variance relationship and the dependency of mean and latent factors, then
#the result should have mean zero and a stable variance across the range of expression.” sctransform vignette.


#In short:

# 1. expression of a gene is modeled by a negative binomial random variable with a mean that depends on library size

# 2. use library size as independent variable in regression model

# 3. it model for each gene, then combine data across genes to fit parameters
# 4. convert UMI counts to residuals akin to the number of standard deviations away from the expected mean.

#Assumptions:
  
#accurate model of the mean-variance relationship
#accurate model of the dependency of mean and latent factors
#Outcome:
  
# 1. mean zero
# 2. stable variance across expression range

#Estimation and transformation
#We will now estimate model parameters and transform data

#The vst function estimates model parameters and performs the variance stabilizing transformation.

#Here we use the log10 of the total UMI counts of a cell as variable for 
#sequencing depth for each cell. After data transformation we plot the model
#parameters as a function of gene mean (geometric mean).

print(dim(counts)) #18202 27113


# We use the Future API for parallel processing;
# set parameters here
future::plan(strategy = 'multicore', workers = 7)
options(future.globals.maxSize = 10 * 1024 ^ 3)

# transform counts:
set.seed(44) # for reproducibility
vst_out <- sctransform::vst(
  counts, # A matrix of UMI counts with genes as rows and cells as columns
  latent_var = c('log_umi'), # The independent variables to regress out as a character vector
  return_gene_attr = TRUE, # Make cell attributes part of the output
  return_cell_attr = TRUE, # Calculate gene attributes and make part of output
  verbosity = 0 # An integer specifying what to show (0: nothing, 1: messages, 2: + progress bar)
)


#Parameters plots
# diagnostic plots: estimated and fitted model parameters
# by default parameters shown are:
# - intercept
# - latent variables, here log_umi
# - overdispersion factor (od_factor)
sctransform::plot_model_pars(
  vst_out, # The output of a vst run
  verbosity = 2 # Messages only, no progress bar
)

print(vst_out$model_str)

#y ~ log_umi

#The distribution of residual mean is cetered around 0:
ggplot(vst_out$gene_attr, aes(residual_mean)) +
  geom_histogram(binwidth=0.01)

#The distribution of residual variance is centered around 1:
ggplot(vst_out$gene_attr, aes(residual_variance)) +
  geom_histogram(binwidth=0.1) +
  geom_vline(xintercept=1, color='red') +
  xlim(0, 10)

#The following plot of the residual variance against the mean: after transformation
#there is no relationship between gene mean and variance.
ggplot(vst_out$gene_attr,
       aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3)

#Check genes with large residual variance. These genes would be markers of expected cell populations. 
#Note how they represent a great range of mean UMI and detection rate values.
dd <- vst_out$gene_attr %>%
  arrange(-residual_variance) %>%
  slice_head(n = 20) %>%
  dplyr::mutate(across(where(is.numeric), round, 2))

dd %>% tibble::rownames_to_column("ID") %>%
  left_join(as.data.frame(rowData(sce))[,c("ID", "Symbol")], "ID") %>%
  DT::datatable(rownames = FALSE)


print(dim(vst_out$y))
vst_out$y[1:10,1:5]
sce
print(assayNames(sce))

# genes that are expressed in fewer than 5 cells are not used and not returned
# so to add vst_out$y as an assay we need to ditch the missing genes completely.
# https://github.com/ChristophH/sctransform/issues/27

geneOverlap <- rownames(sce) %in% rownames(vst_out$y)
if(!all(geneOverlap)) #if not all TRUE
{
  table(rownames(sce) %in% rownames(vst_out$y))
  tmpInd <- which(rownames(sce) %in% rownames(vst_out$y))
  sce <- sce[tmpInd,]
  assayNames(sce)
}

sce #18202 27113 # all good
## reading 10X data with vector above adds a prefix to sce colnames
# so we will not pass vstMat colnames when copying it in a assay slot,
# but must first check that barcodes are indeed in the same order
# in sce and vstMat. == this is comment to the line 376
vstMat <- as(vst_out$y[rownames(sce),], "dgCMatrix")

all(colnames(vstMat) == sce$Barcode)
all(rownames(vstMat) == rownames(sce))


assay(sce, "sctrans_norm", withDimnames=FALSE) <- vstMat

saveRDS(sce, paste0(path, "/results/sarcoidosis_postSct.Rds"))

sce
