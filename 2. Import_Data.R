## IMPORTING THE DATA ##

# Required packages were loaded.

library(GEOquery)            # Loaded Biobase and BiocGenerics packages.

AccId <- "GSE138260"         # Assigned a variable to GSE series dataset.

GEOSeries <- getGEO(AccId)   # Assigned a variable for getGEO function. 
                             # in GEOquery package to direct download
                             # GSE series.

View(GEOSeries)              # Spreadsheet-style data viewer on R 
                             # object was viewed to understand the overview.

# The number of platforms used were checked.

length(GEOSeries)           # If length > 1, it meant that some datasets on GEO 
                            # were derived from different microarray platforms.

# The [[ form subsets by selecting single element using integer or character 
# indices.

GEOSeries <- GEOSeries[[1]]

# If more than one platform is present, they can be analysed by 
# changing the number inside the [[...]] e.g. GEOSeries <- GEOSeries[[2]].


# pData accesses information on recorded experimental phenotypes:
# phenotypic data (co-variates) and meta-data (description of co-variates).

# fData accesses feature data and meta-data - gene annotation.
# exprs accesses matrix of expression data.
# Example: pData(GEOSeries)[1,] displays phenotypic data of 1st row through
# every column.


#------------------------------------------------------------------------------#





