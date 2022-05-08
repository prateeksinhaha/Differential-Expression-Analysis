##  HEATMAP FOR N TOP GENES  ##

# Here the top n genes from the differential expression analysis were taken
# and heatmap was produced. The default colour palette goes from low expression 
# in blue to high expression in red, which is a good alternative to the 
# traditional red/green heatmaps which are not suitable for those with forms of 
# colour-blindness.

library(dplyr)
library(pheatmap)
library(tibble)

# Subset of top genes

topNgenes <- filtered_Data$ID[1:20]


# Heatmap was plotted
png(filename = "Expression Values Heatmap of Top 20 Genes.png"
    ,res=720, width = 10, height = 6, units = "in")

pheatmap(exprs(GEOSeries)[topNgenes,], annotation_col = SampleInfo, 
         labels_row = filtered_Data$GENE_SYMBOL, fontsize = 10, 
         main = "Expression Values Heatmap of Top 20 Genes",
         fontsize_row = 9, fontsize_col = 9)
dev.off()

#------------------------------------------------------------------------------#

