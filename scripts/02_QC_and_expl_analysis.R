
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
library(xlsx)

path = '/scratch/cube/sango/sarcoidosis_project'
setwd(path)

#sample meta data

samplesheet = read.xlsx('Sample_annotation_sheet_v2.xlsx',
                        sheetIndex = 3)[1:10, ]
samplesheet = samplesheet %>%
  dplyr::filter(!Sample.Description == "MMP12-inhibitor treated") %>%
  dplyr::filter(!Sample.Name == "FP22")


bp.params = MulticoreParam(workers = 7)

#loading a single sample

sample.path = 'COUNT/FP19_transcriptome/filtered_feature_bc_matrix/'
sce.sing = read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)
sce.sing

#The counts matrix
#Compared to bulk RNA-seq, Single-cell RNA-seq data is sparse, i.e. there 
#many missing values or zeroes. This is particularly true with 
#droplet-based methods such as 10X, mostly because:
  
# - any given cell does not express each gene
# - the library preparation does not capture all transcript the cell does express
# - the sequencing depth per cell is far lower and so fewer of the expressed genes are detected
# - We can access the counts matrix with counts. Given the large number of droplets in a sample, count matrices can be large.


dim(counts(sce.sing))
#32285  5416 (genes, cells)


rowData(sce.sing) #Details about the “features” (in this case genes)

colData(sce.sing) #The rows of this table correspond to the data in the columns of the count matrix; 
#the row names of this table will match the column names of the counts 
#matrix - currently these are the droplet barcodes.

#if we look at the number of genes detected in each cell, we can see that this ranges from 20 to 8686, with a median of 1851
genesPerCell = colSums(counts(sce.sing) > 0)
summary(genesPerCell)

plot(density(genesPerCell), main="", xlab="Genes per cell")

#Total UMI for a gene versus the number of times detected
tmpCounts <- counts(sce.sing)[,1:1000]

plot(rowSums(tmpCounts),
     rowMeans(tmpCounts > 0),
     log = "x",
     xlab="total number of UMIs",
     ylab="proportion of cells expressing the gene"
)

rm(tmpCounts)

#Distribution of counts for a gene across cells

#We could also look at the distribution of counts for individual genes across all cells.
#The plot below shows this distribution for the top 20 genes detected.

rel_expression <- t( t(counts(sce.sing)) / colSums(counts(sce.sing))) * 100
rownames(rel_expression) <- rowData(sce.sing)$Symbol
most_expressed <- sort(rowSums( rel_expression ),T)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)

#Load multiple samples

samples_list <- paste0(samplesheet$Sample.Name, "_transcriptome")
list_of_files = str_c(path, 
                      "/COUNT/",
                      samples_list, 
                      "/filtered_feature_bc_matrix")
names(list_of_files) <- samplesheet$Sample.Name
list_of_files

bp.params = MulticoreParam(workers = 2)

sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)

sce
colData(sce)

#Modify the droplet annotation

colData(sce) <- colData(sce) %>% 
  as.data.frame() %>%
  rownames_to_column("RowName") %>% 
  mutate(SampleNum = str_extract(RowName, "^[0-9]+")) %>%
  mutate(Barcode = str_replace(Barcode, "1$", SampleNum)) %>%
  left_join(samplesheet, by=c(Sample="Sample.Name")) %>%
  #rename(SampleId=Sample) %>% 
  #rename(Sample=SampleName) %>%    
  #     mutate(Sample = case_when(
  #         SampleId == "SRR9264351" ~ str_c(Sample, "a"),
  #         SampleId == "SRR9264352" ~ str_c(Sample, "b"),
  #         TRUE ~ Sample)) %>% 
  column_to_rownames("RowName") %>% 
  dplyr::select(Sample, Barcode, Phenotype, Sample.Description) %>% #SampleGroup=Phenotype
  DataFrame()


#Undetected genes
#Although the count matrix has 36601 genes, many of these will not have been detected in any droplet.

detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)

#detected_genes
#FALSE  TRUE 
#6932 25353 , 21% of genes not detected in any droplet, let's remove them:

sce <- sce[detected_genes,]


#Annotate genes
#In order to assess the percentage of mitochondrial UMIs, we will need to be able
#to identify mitochondrial genes. The simplest way to do this is to annotate 
#the genes with their chromosome of origin.


ah <- AnnotationHub()
ens.mm.102 <- query(ah, c("Mus musculus", "EnsDb", 102))[[1]]

genes <- rowData(sce)$ID
gene_annot <- AnnotationDbi::select(ens.mm.102, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME")) %>%
  set_names(c("ID", "Chromosome"))

rowData(sce)$Chromosome = gene_annot$Chromosome[match(rowData(sce)$ID, gene_annot$ID)]

rowData(sce)

na.omit(rowData(sce))

unique(rowData(sce)$Chromosome)

#Add per cell QC metrics
#We can now add the per cell QC metrics to the droplet annotation using the 
#function addPerCellQC. In order to get the metrics for the subset of mitochondrial
#genes, we need to pass the function a vector indicating which genes are mitochondrial.


is.mito <- which(rowData(sce)$Chromosome=="MT")

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)

colData(sce)

#The function has added six columns to the droplet annotation:

# sum: total UMI count
# detected: number of features (genes) detected
# subsets_Mito_sum: number of UMIs mapped to mitochondrial transcripts
# subsets_Mito_detected: number of mitochondrial genes detected
# subsets_Mito_percent: percentage of UMIs mapped to mitochondrial transcripts
# total: also the total UMI count


#QC metric distribution

#Before moving on to do the actual cell filtering, it is always a good idea
#to explore the distribution of the metrics across the droplets.


plotColData(sce, x="Sample", y="sum",other_fields="Phenotype") + 
  facet_wrap(~Phenotype, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Total count of UMIs")

plotColData(sce, x="Sample", y="detected", other_fields="Phenotype") + 
  facet_wrap(~Phenotype, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Detected features (genes)")

plotColData(sce, x="Sample", y="subsets_Mito_percent", other_fields="Phenotype") + 
  facet_wrap(~Phenotype, nrow=1, scales = "free_x") +
  ggtitle("Mito percent")


#scatter plot shows the extent to which library size and numbers of genes detected are correlated.

colData(sce) %>% 
  as.data.frame() %>% 
  arrange(subsets_Mito_percent) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mito_percent > 10)) + 
  facet_wrap(vars(Phenotype))


#Identification of low-quality cells with adaptive thresholds


#One could use hard threshold for the library size, number of genes detected and 
#mitochondrial content based on the distributions seen above. These would need 
#vary across runs and the decision making process is somewhat arbitrary. It may
#therefore be preferable to rely on outlier detection to identify cells that 
#markedly differ from most cells.

#We saw above that the distribution of the QC metrics is close to Normal. Hence,
#we can detect outliers using the median and the median absolute deviation (MAD)
#from the median (not the mean and the standard deviation which both are sensitive to outliers).

#For a given metric, an outlier value is one that lies over some number of MADs
#away from the median. A cell will be excluded if it is an outlier in the part 
#of the range to avoid, for example low gene counts, or high mitochondrial content.
#For a normal distribution, a threshold defined with a distance of 3 MADs from the
#median retains about 99% of values.

#The scater function isOutlier can be used to detect outlier cells based on any 
#metric in the colData table. It returns a boolean vector that identifies outliers.
#By default it will mark any cell that is 3 MADS in either direction from the median as an outlier.


#Library size
#With library size we wish to identify outliers that have very low 
#library sizes, this indicates that the droplets either contain poor
#quality cells, perhaps damaged or dying, or do not contain a cell at all.

low_lib_size <- isOutlier(sce$sum, log=TRUE, type="lower")
table(low_lib_size)

attr(low_lib_size, "thresholds")

# lower   higher => threshold values for library size. none of the cell have less than 400 UMIs, which
#400.2822      Inf

colData(sce)$low_lib_size <- low_lib_size

plotColData(sce, 
            x="Sample", 
            y="sum",
            other_fields="Phenotype", 
            colour_by = "low_lib_size") + 
  facet_wrap(~Phenotype, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Total count", title = "Total count") +
  guides(colour=guide_legend(title="Discarded"))


#Number of genes
#As with the library size, we will log tranform the number of genes detected 
#prior to filtering using the median absolute deviation.

low_n_features <- isOutlier(sce$detected, log=TRUE, type="lower")
table(low_n_features)

#FALSE  TRUE 
#29967  1359 #This has excluded out 1359 cells. The threshold value was:

attr(low_n_features, "thresholds")[1]
#309.1619


colData(sce)$low_n_features <- low_n_features
plotColData(sce, 
            x="Sample", 
            y="detected",
            other_fields="Phenotype", 
            colour_by = "low_n_features") + 
  facet_wrap(~Phenotype, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Genes detected", title = "Genes detected") +
  guides(colour=guide_legend(title="Discarded"))



#Mitochondrial content
#For the mitochondrial content the exclusion zone is in the higher part of the distribution.
#For this reason we do not need to worry about log transforming the data as want to 
#remove the long right hand tail anyway.


high_Mito_percent <- isOutlier(sce$subsets_Mito_percent, type="higher")
table(high_Mito_percent)

#FALSE  TRUE 
#27285  4041  #This has removed 4041 cells in total. The upper threshold value:

attr(high_Mito_percent, "thresholds")[2]

#7.569205

colData(sce)$high_Mito_percent <- high_Mito_percent
plotColData(sce,  
            x="Sample",
            y="subsets_Mito_percent",
            other_fields="Phenotype",
            colour_by = "high_Mito_percent") + 
  facet_wrap(~Phenotype, nrow=1, scales = "free_x") + 
  labs(y = "Percentage mitochondrial UMIs",
       title = "Mitochondrial UMIs") +
  guides(colour=guide_legend(title="Discarded"))


#Summary of discarded cells
#Having applied each of the three thresholds separately, we can now 
#combine them to see how many droplets in total we will be excluding.

data.frame(`Library Size` = sum(low_lib_size),
           `Genes detected` = sum(low_n_features),
           `Mitochondrial UMIs` = sum(high_Mito_percent),
           Total = sum(low_lib_size | low_n_features | high_Mito_percent))

#All three filter steps at once
#The three steps above may be run in one go using the quickPerCellQC function.
#This creates a DataFrame with 4 columns containing TRUE/FALSE - one for each
#filter metric and one called “discard” that combined the three logicals.

cell_qc_results <- quickPerCellQC(colData(sce),
                                  percent_subsets=c("subsets_Mito_percent"))
colSums(as.data.frame(cell_qc_results))

sce$discard <- cell_qc_results$discard

#Filtering out poor quality cells

sce.filtered <- sce[, !sce$discard]

colData(sce.filtered) <- colData(sce.filtered)[,1:3]
sce.filtered <- addPerCellQC(sce.filtered, BPPARAM = bp.params)

colData(sce.filtered) #27113 cells instead of 31326


#Mitochondrial content versus library size
#A useful diagnostic plot for assessing the impact of the filtering is to do a 
#scatter plot of the mitochondrial content against the library size. We can 
#overlay our final filter metric using the point colour.

plotColData(sce, 
            x="sum", 
            y="subsets_Mito_percent", 
            other_fields="Sample",
            colour_by="discard") +
  facet_wrap(~Sample, ncol=5, scale="free_x")


#QC and Filtering based on sparsity
#The approach above identified poor-quality using thresholds on the 
#number of genes detected and mitochondrial content. We will here 
#specifically look at the sparsity of the data, both at the gene and cell levels.


#Sparsity plots
#We will compute:
  
#the cell sparsity: for each cell, the proportion of genes that are not detected
#the gene sparsity: for each gene, the proportion of cells in which it is not detected

#To help calculate the gene sparsity we can generate QC metrics for genes with 
#addPerFeatureQC. This adds two columns to the gene annotation (rowData):
  
#mean - the mean UMI count for the gene across all cells
#detected - the percentage of cells in which the gene was detected

sce.filtered <- addPerFeatureQC(sce.filtered, BPPARAM = bp.params)
rowData(sce.filtered)

#Now we can calculate sparsity using the “detected” columns in the colData and the rowData.

colData(sce.filtered)$cell_sparsity <- 1 - (colData(sce.filtered)$detected / nrow(sce.filtered))
rowData(sce.filtered)$gene_sparsity <- (100 - rowData(sce.filtered)$detected) / 100

#The cell sparsity plot shows that most cells have between 85% and 99% 0’s, which is typical.
hist(sce.filtered$cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity", main="")

#The gene sparsity plot shows that a large number of genes are almost never detected, which is also regularly observed.
hist(rowData(sce.filtered)$gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity", main="")

#Filter by sparsity
#We could remove cells with sparsity higher than 0.99, and/or mitochondrial content higher than 10%.

#Genes detected in a few cells only are unlikely to be informative and would hinder 
#normalization. We will remove genes that are expressed in fewer than 20 cells.

sparse.cells <- sce.filtered$cell_sparsity > 0.99 #all FALSE, meaning all cells from prev. filtering step have sparsity less this
mito.cells <- sce.filtered$subsets_Mito_percent > 10 # NULL since that is already filtered in the previous step

min.cells <- 1 - (20 / ncol(sce.filtered))
sparse.genes <- rowData(sce.filtered)$gene_sparsity > min.cells

#Number of genes removed:

table(sparse.genes)

## FALSE  TRUE 
## 18202  7151


#Number of cells removed:

table(sparse.cells)
table(mito.cells)

sum(sparse.genes)
sce.filtered <- sce.filtered[!sparse.genes, ]

rowData(sce.filtered) #18202, additional 7151 genes were filtered out bc they have high sparsity
colData(sce.filtered) #27113 cells

#save sce object for next step: 03_Normalization

saveRDS(sce.filtered, paste0(path,"/results/Sarcoidosis_filtered_genes.rds"))
