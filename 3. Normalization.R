## NORMALIZATION CHECK AND USAGE OF SCALES ##


# To check if the expression data was RMA (Robust Multiarray Average) 
# normalized or not (presence of low-expressing genes or not):

pData(GEOSeries)$data_processing[1]

# For printing the distribution of expression levels:

library(utils)
summary(exprs(GEOSeries))
x <- summary(exprs(GEOSeries))
write.csv(x, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/summary.csv")

# For visualization and statistical analysis, the data was inspected to 
# know what scale the data were presented in. The methods used assumes that
# the data were on a log2 scale in the range of 0 to 20.

# If the expression values do not lie in the range of 0 to 20, log2 
# transformation has to be performed as follows:

# exprs(GEOSeries) <- log2(exprs(GEOSeries))
# summary(exprs(GEOSeries))      # To check the summary again

# Here expression values lie in the range of 0 to 20,
# so no log2 transformation has to be done

#------------------------------------------------------------------------------#
