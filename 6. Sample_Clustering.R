## SAMPLE CLUSTERING ##


library(pheatmap)  # Visualizes the correlations between the samples by
                   # hierarchical clustering

# cor() function is used to compute the correlations of microarray and gene 
# expression profiles for a set of genes. It calculates correlation on the scale
# of 0 to 1, in a pairwise fashion between all samples, then visualize on 
# heat map.

corMatrix <- cor(exprs(GEOSeries), use = "c")  # use = "c" stops an error if 
z   <- corMatrix                               # there are  missing data points
write.csv(z, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/corMatrix.csv")

pheatmap(corMatrix)

rownames(SampleInfo)  # Printed the row names to check if it matched the 
                      # correlation matrix

colnames(corMatrix)

# If it does not match, the row names were forced to match the columns
# rownames(SampleInfo) <- colnames(corMatrix)

png(filename = "Sample Clustering of Correlation Matrix (Expression Values).png"
    ,res=720, width = 10, height = 6, units = "in")
pheatmap(corMatrix, annotation_col = SampleInfo, 
         main = "Sample Clustering of Correlation Matrix (Expression Values)",
         fontsize = 8, fontsize_row = 9, fontsize_col = 9)
dev.off()

#------------------------------------------------------------------------------#

