## DIFFERENTIAL EXPRESSION ANALYSIS ##


# The limma (Linear models for microarray) package is used to perform 
# differential expressions to compare the sample groups. 
# Example: SampleInfo$group. 

# A design matrix of 0 and 1 is created, one row for each sample and one column 
# for each sample group. The column names can be renamed as required.

# The presence of any lowly-expressed genes, that will affect the quality of DE 
# analysis was checked in the expression data.

# A problem in performing statistical analysis like limma is the inference of 
# type 1 statistical errors (false positive). One simple way to reduce 
# this is by filtering data. For example, since not all genes are expressed
# in all tissues and many genes will not be expressed in any sample, these non-
# expressed genes can be removed.

# The cut off can be at the median of the expression values, which means to 
# consider around 50% of the genes that will not be expressed.

# The expressed genes that are present in more than 2 samples were kept. It was 
# seen that around half of the genes are not qualified as an "expressed" 
# gene since the cut-off is the median value.


library(limma)

design <- model.matrix(~group + 0, GEOSeries)  #~0 means no intercept

design  # Optionally View(design) can be used to view it as data frame
View(design)

# The column names were renamed

colnames(design) <- levels(gs)
design

View(design)
write.csv(design, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/design.csv")

# The median expression level (cutoff) was calculated

# cutoff <- median(exprs(GEOSeries))

# cutoff  # The cutoff was printed

# TRUE or FALSE for whether each gene is "expressed" in each sample:

# is_expressed <- exprs(GEOSeries) > cutoff

# Identified genes expressed in MORE than 2 samples

# keep <- rowSums(is_expressed) > 2   #checks how many TRUE (1-Sum) in each gene

# Genes that are removed / retained can be checked by

# table(keep)

# Only those genes which were expressed were selected

# GEOSeries <- GEOSeries[keep,]

# GEOSeries


## OUTLIER MANAGEMENT ##

# This has to be done carefully so the filtered data won't be too biased. 
# 'weights' were calculated to define the reliability of each sample. 
# The 'arrayWeights' function (limma) assigned a score to each sample, with a 
# value of 1 implying equal weight. Samples with score less than 1 were 
# down-weighed, or else up-weighed.

# array_weight <- arrayWeights(exprs(GEOSeries) , design)

# Now by having a design matrix, it was needed to estimate the coefficients. 
# For this design, the replicate arrays were essentially averaged for each 
# sample level. In addition, standard deviations for each gene, and the average 
# intensity for the genes across all microarrays were calculated.

# There is a need to inform the package which pairwise contrasts has to be made. 
# For this experiment, treatment (two types of texane drugs) and control in each 
# serum type is to be contrasted. So there are 4 contrasts to be specified.


# To do the statistical comparisons, Limma uses Bayesian statistics to minimize 
# type 1 error. The eBayes function performs the tests. To summarize the results 
# of the statistical test, 'topTable' adjusts the p-values and return the top
# genes that meet the cutoffs that are supplied as arguments; while 
# 'decideTests'makes calls for DEGs by adjusting the p-values and applying a 
# logFC cutoff similar to topTable.

# Fitting the coefficients

fit <- lmFit(GEOSeries, design) 

# This fits linear model for each gene in a given series of arrays


# The function makeContrasts expresses contrasts between a set of parameters 
# as a numeric matrix. The parameters are usually the coefficients from a 
# linear model fit, so the matrix specifies which comparisons between the 
# coefficients are to be extracted from the fit.

contrasts <- makeContrasts(AD - CTRL, levels = design)

write.csv(contrasts, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/contrasts.csv")

fit2 <- contrasts.fit(fit, contrasts) # contrasts.fit function computes 
# contrasts from linear model fit

fit2 <- eBayes(fit2, 0.01)  # empirical Bayes statistics for differential expression
topTable <- topTable(fit2, adjust="fdr", sort.by="logFC", number=250)
topTable <- subset(topTable, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","GENE_SYMBOL","GENE_NAME","SEQUENCE"))

# Table of top genes from linear model fit

# To see the results of the second contrast (if it exists)
# topTable(fit2, coef=2) and etc

topTable2 <- topTable(fit2, adjust="fdr", sort.by="logFC", number=Inf)


# To know how many genes are differentially expressed overall, decideTest 
# function is used:

library(edgeR)

summary(decideTests(fit2, adjust.method="fdr", p.value=0.05))
decideTests(fit2, adjust.method="fdr", p.value=0.05)
table(decideTests(fit2, adjust.method="fdr", p.value=0.05))



de <- decideTests(fit2, method="global")
up <- which(de[,1] == 1)     # Up Regulated Genes
upGenes <- fit2[up, ]
upTable <- topTable(upGenes, n=Inf)
dim(upTable)

down <- which(de[,1] == -1)  # Down Regulated Genes
downGenes <- fit2[down,]
downTable <- topTable(downGenes, n=Inf)
dim(downTable)

notSig <- which(de[,1] == 0) # Not Significant Genes
notSigGenes <- fit2[notSig,]
notSigTable <- topTable(notSigGenes, n=Inf)
dim(notSigTable)


#------------------------------------------------------------------------------#
