
# Part 6:Differential expression and abundance between conditions

library(edgeR)
library(ggplot2)

colLabels(sce) <- sce$louvain
tab <- table(colLabels(sce), sce$Phenotype)
tab

tab <- table(colLabels(sce), sce$SampleName)

pheatmap::pheatmap(tab,
                   border_color      = NA,
                   drop_levels       = TRUE,
                   cluster_rows      = !FALSE,
                   cluster_cols      = !FALSE
)

p1 <- plotTSNE(sce, colour_by="Phenotype", text_by="label", point_size=0.3)
p2 <- plotTSNE(sce, colour_by="SampleName", point_size=0.3)
gridExtra::grid.arrange(p1, p2+facet_wrap(~colData(sce)$SampleName), ncol=2)


p1 <- plotTSNE(sce, colour_by="Phenotype", text_by="label", point_size = 0.3)
p2 <- plotTSNE(sce, colour_by="SampleName", point_size = 0.3) +
  facet_wrap(~colData(sce)$Phenotype)

p1 + p2 + plot_layout(widths = c(1, 2))



# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
columnsToUse <- c("SampleName", "Phenotype", "louvain")
colData(sce) <- colData(sce) %>% as.data.frame() %>% dplyr::select(all_of(columnsToUse)) %>% DataFrame
summed <- aggregateAcrossCells(sce, 
                               id = DataFrame(
                                 label = sce$louvain,
                                 sample = sce$SampleName
                               )
)
colData(summed) %>% head(3)

#1.3.2 Performing the DE analysis
#1.3.2.1 Introduction
#The DE analysis will be performed using quasi-likelihood (QL)
#This uses a negative binomial generalized linear model (NB GLM)
#to handle overdispersed count data in experiments with limited 
#replication. In our case, we have biological variation with few
#replicates per sample group, so edgeR (or its contemporaries)
#is a natural choice for the analysis.

#We do not use all labels for GLM fitting as the strong DE between
#labels makes it difficult to compute a sensible average abundance
#to model the mean-dispersion trend. Moreover, label-specific batch
#effects would not be easily handled with a single additive term in
#the design matrix for the batch. Instead, we arbitrarily pick one
#of the labels to use for this demonstration.


labelToGet <- "7"
current <- summed[,summed$label==labelToGet]
colData(current)

# Creating up a DGEList object for use in edgeR:
countsToUse <- counts(current)
colnames(countsToUse) <- colData(current)$SampleName
y <- DGEList(countsToUse, samples=colData(current))
y


#1.3.2.2 Pre-processing
#A typical step in bulk RNA-seq data analyses is to remove samples with very
#low library sizes due to failed library preparation or sequencing. The very
#low counts in these samples can be troublesome in downstream steps such as
#normalization or for some statistical approximations used in the DE analysis.
#In our situation, this is equivalent to removing label-sample combinations
#that have very few or lowly-sequenced cells. The exact definition of “very low” 
#will vary, but in this case, we remove combinations containing fewer than 20 cells

discarded <- current$ncells < 20
y <- y[,!discarded]
summary(discarded)
##    Mode   FALSE 
## logical       6

#Another typical step in bulk RNA-seq analyses is to remove genes that are
#lowly expressed. This reduces computational work, improves the accuracy
#of mean-variance trend modelling and decreases the severity of the multiple
#testing correction. Genes are discarded if they are not expressed above
#a log-CPM threshold in a minimum number of samples (determined from the
#size of the smallest treatment group in the experimental design).

keep <- filterByExpr(y, group=current$Phenotype)
y <- y[keep,]
summary(keep)
#Mode      FALSE    TRUE 
#logical   10236    7966 -> went from 18202 to 7966 genes to keep


#Finally, we correct for composition biases by computing normalization
#factors with the trimmed mean of M-values method (Robinson and Oshlack 2010).
#Counts for our pseudo-bulk samples are large enough to apply bulk normalization methods.

y <- calcNormFactors(y)
y$samples

#As part of the usual diagnostics for a bulk RNA-seq DE analysis, we generate a
#mean-difference (MD) plot for each normalized pseudo-bulk profile. This should
#exhibit a trumpet shape centered at zero indicating that the normalization 
#successfully removed systematic bias between profiles. Lack of zero-centering
#or dominant discrete patterns at low abundances may be symptomatic of deeper
#problems with normalization, possibly due to insufficient cells/reads/UMIs
#composing a particular pseudo-bulk profile.

par(mfrow=c(2,3))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}

#We also generate a multi-dimensional scaling (MDS) plot for the pseudo-bulk
#profiles. This is closely related to PCA and allows us to visualize the
#structure of the data. Here, the aim is to check whether samples separate
#by our known factors of interest. Strong separation foreshadows a large 
#number of DEGs in the subsequent analysis.

pca <- prcomp(t(y$counts), center = T, scale = T)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pcaresults <- summary(pca)

scree.data <- as.data.frame(pcaresults$importance) # eigenvalues / standard deviations
score.data <- as.data.frame(pcaresults$x) # coordinates of the samples (i.e., observations)
loading.data <- as.data.frame(pcaresults$rotation) # loading scores for individual variables

score.data <- score.data[, 1:2]
score.data$Group <- y$samples$Phenotype
#rownames(score.data) = df$X[df$Day == day]
perc.var = 100 * pcaresults[["importance"]][2, ]

scores_plot <- ggplot(score.data, aes(x = PC1, y = PC2, label = rownames(score.data))) +
  geom_text(aes(colour = factor(Group))) +
  #geom_label_repel(label = rownames(score.data)) +
  xlab(paste("PC1:", perc.var[1], "%")) +
  ylab(paste("PC2:", perc.var[2], "%"))
  #my.theme +

scores_plot


#1.3.2.3 Statistical modelling
#Our aim is to test whether the log-fold change between sample
#groups is significantly different from zero.

design <- model.matrix(~factor(Phenotype), y$samples)
design

#We estimate the negative binomial (NB) dispersions with estimateDisp().
#The role of the NB dispersion is to model the mean-variance trend, which
#is not easily accommodated by QL dispersions alone due to the quadratic
#nature of the NB mean-variance trend.

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)


#We also estimate the quasi-likelihood dispersions with glmQLFit() (Chen, 
#Lun, and Smyth 2016). This fits a GLM to the counts for each gene and 
#estimates the QL dispersion from the GLM deviance. We set robust=TRUE
#to avoid distortions from highly variable clusters (Phipson et al. 2016).
#The QL dispersion models the uncertainty and variability of the per-gene
#variance - which is not well handled by the NB dispersions, so the two
#dispersion types complement each other in the final analysis.

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

summary(fit$df.prior)

#QL dispersion estimates for each gene as a function of abundance. Raw estimates
#(black) are shrunk towards the trend (blue) to yield squeezed estimates (red).

plotQLDisp(fit)

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))

topTags(res)$table

#1.3.3 Putting it all together
#Now that we have laid out the theory underlying the DE analysis,
#we repeat this process for each of the labels. This is conveniently
#done using the pseudoBulkDGE function from scran, which will loop
#over all labels and apply the exact analysis described above to each
#label. To prepare for this, we filter out all sample-label 
#combinations with insufficient cells.

summed.filt <- summed[,summed$ncells >= 20]

#We construct a common design matrix that will be used in the
#analysis for each label. Recall that this matrix should have
#one row per unique sample (and named as such), reflecting the
#fact that we are modelling counts on the sample level instead of the cell level.

# Pulling out a sample-level 'targets' data.frame:
targets <- colData(sce)[!duplicated(sce$SampleName),] %>%
  data.frame() %>%
  dplyr::select(-louvain)

# Constructing the design matrix:
design <- model.matrix(~factor(Phenotype), data=targets)
rownames(design) <- targets$SampleName


#We then apply the pseudoBulkDGE function to obtain a list of DE genes
#for each label. This function puts some additional effort into automatically
#dealing with labels that are not represented in all sample groups, for which
#a DE analysis between conditions is meaningless; or are not represented in
#a sufficient number of replicate samples to enable modelling of biological variability.

summed.filt$Phenotype <- factor(summed.filt$Phenotype)

de.results <- pseudoBulkDGE(summed.filt, 
                            label = summed.filt$label,
                            design = ~Phenotype,
                            coef = "PhenotypeSarcoidosis",
                            condition = summed.filt$Phenotype
)


#We examine the numbers of DEGs at a FDR of 5% for each label using the
#decideTestsPerLabel function. Note that genes listed as NA were either
#filtered out as low-abundance genes for a given label’s analysis, or
#the comparison of interest was not possible for a particular label,
#e.g., due to lack of residual degrees of freedom or an absence
#of samples from both conditions.

is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)

#     -1  0   1    NA
# 1   0  8538   3  9661
# 10 53  7498 190 10461
# 11  0  3699   0 14503
# 12  0  9459   2  8741
# 13  3  9537  25  8637
# 14  0  4742   3 13457
# 15  8  2763   9 15422
# 16  8  8466   4  9724
# 17  0  6944  11 11247
# 18 18  8682  32  9470
# 19  0  5481   0 12721
# 2   0  8640   1  9561
# 20  0  6350   0 11852
# 21  0  3607   0 14595
# 3   0  7882   4 10316
# 4   0  2866   0 15336
# 5   0  8976   4  9222
# 6  14  9829  19  8340
# 7   0  7965   1 10236
# 8   0 11168   0  7034
# 9   1  8552   0  9649

#For each gene, we compute the percentage of cell types in which
#that gene is upregulated or downregulated. (Here, we consider a
#gene to be non-DE if it is not retained after filtering.).

# Upregulated across most cell types in Sarcoidosis
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

#Plcb4         Fabp5         Itgad 4931431B13Rik        Lgals3         Ier5l         Gpnmb        Abcb1a         Eomes 
#0.3333333     0.2380952     0.2380952     0.2380952     0.2380952     0.1904762     0.1904762     0.1428571     0.1428571 
#Gzmk 
#0.1428571 

# Downregulated across cell types in Sarcoidosis
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)


#We further identify label-specific DE genes that are significant
#in our label of interest yet not DE in any other label. As hypothesis
#tests are not typically geared towards identifying genes that are not
#DE, we use an ad hoc approach where we consider a gene to be consistent
#with the null hypothesis for a label if it fails to be detected even
#at a generous FDR threshold of 50%.

remotely.de <- decideTestsPerLabel(de.results, threshold=0.5)
not.de <- remotely.de==0 | is.na(remotely.de)

# first cluster in is.de
cx <- colnames(is.de)[2]
cx


other.labels <- setdiff(colnames(not.de), cx)
unique.degs <- is.de[,cx]!=0 & rowMeans(not.de[,other.labels])==1
unique.degs <- names(which(unique.degs))
head(unique.degs)


# Choosing the top-ranked gene for inspection:
de.inspec <- list()
de.inspec[[cx]] <- de.results[[cx]] 
de.inspec[[cx]] <- de.inspec[[cx]][order(de.inspec[[cx]]$PValue),]
de.inspec[[cx]] <- de.inspec[[cx]][rownames(de.inspec[[cx]]) %in% unique.degs,]

sizeFactors(summed.filt) <- NULL
plotExpression(logNormCounts(summed.filt), 
               features=rownames(de.inspec[[cx]])[1],
               x="Phenotype", colour_by="Phenotype", 
               other_fields="label") + 
  facet_wrap(~label) +
  ggtitle(glue::glue("{cx}: {rownames(de.inspec[[cx]])[1]}"))



#1.4 Differential abundance between conditions
#1.4.1 Overview
#In a DA analysis, we test for significant changes in per-label
#cell abundance across conditions. This will reveal which cell
#types are depleted or enriched upon treatment, which is arguably
#just as interesting as changes in expression within each cell type.
#The DA analysis has a long history in flow cytometry (Finak et al.
#2014; Lun, Richard, and Marioni 2017) where it is routinely used
#to examine the effects of different conditions on the composition
#of complex cell populations. By performing it here, we effectively
#treat scRNA-seq as a “super-FACS” technology for defining relevant
#subpopulations using the entire transcriptome.

#We prepare for the DA analysis by quantifying the number of cells assigned to each label (or cluster).


abundances <- table(sce$louvain, sce$SampleName) 
abundances <- unclass(abundances) 
head(abundances)

#     FP19 FP21 FP24 FP25 FP27 FP28
# 1   384  442  375  452  275  331
# 2    92  330  194   18   56   77
# 3   163  116  115   89  192  159
# 4    19   29  151  282   53   99
# 5   509  568  277  830  588  194
# 6   376  191  141  175  118  146
# 7   185   67   64  144   47  429
# 8   352  229   98  114   87  348
# 9   204  445  157   56   93  221
# 10  505   48   71  535   49  283
# 11  134  497  371   47   96   82
# 12  817  361  243  473  330  172
# 13  492  383  554  134  259  263
# 14  392  177    3   20   22   22
# 15   34   26  257   92   55  133
# 16  151  508  158   44  150  123
# 17   75  116   46  287  171   58
# 18   90  247  549  153  367  230
# 19   51   99   41   22   28 1104
# 20  109   55   47  162  196  160
# 21   25   84   56   38   63   70
# 22    0  239  231    1  331    0

# Attaching some column metadata.
extra.info <- colData(sce)[match(colnames(abundances), sce$SampleName),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

#We filter out low-abundance labels as previously described. This avoids cluttering the result
#table with very rare subpopulations that contain only a handful of cells. For a DA analysis
#of cluster abundances, filtering is generally not required as most clusters will not be of
#low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).

keep <- filterByExpr(y.ab, group=y.ab$samples$Phenotype)
y.ab <- y.ab[keep,]
summary(keep)
#   Mode    TRUE 
#logical      22


#Unlike DE analyses, we do not perform an additional normalization step with
#calcNormFactors(). This means that we are only normalizing based on the
#“library size”, i.e., the total number of cells in each sample. Any changes
#we detect between conditions will subsequently represent differences 
#in the proportion of cells in each cluster.

#Here, the log-fold change in our model refers to the change in cell abundance
#between sample groups, rather than the change in gene expression.

design <- model.matrix(~factor(Phenotype), y.ab$samples)

#We use the estimateDisp() function to estimate the NB dispersion for each cluster. 
#We turn off the trend as we do not have enough points for its stable estimation


y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)

plotQLDisp(fit.ab, cex=1)


res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))

topTags(res)

#topTags(res)
# Coefficient:  factor(Phenotype)Sarcoidosis 
# logFC   logCPM         F       PValue          FDR
# 22  9.3135663 14.98256 68.238970 1.388919e-06 3.055621e-05
# 10 -2.8546676 15.74002 16.892984 1.186161e-03 1.304777e-02
# 7  -1.9785876 15.06108  8.232980 1.296132e-02 9.504968e-02
# 11  1.9348839 15.43908  7.052183 1.955460e-02 1.075503e-01
# 2   1.6886052 14.75203  4.705114 4.884979e-02 2.149391e-01
# 18  1.4359550 15.96638  4.161066 6.184181e-02 2.150901e-01
# 16  1.4112406 15.31237  3.934742 6.843776e-02 2.150901e-01
# 19 -2.7764404 15.55434  3.760059 9.230080e-02 2.538272e-01
# 8  -0.8983275 15.41975  1.753077 2.079019e-01 5.082048e-01
# 21  0.6938517 13.66760  1.073228 3.187749e-01 6.885652e-01
