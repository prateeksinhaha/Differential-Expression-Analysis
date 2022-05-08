# SELECTING TOP GENES #

# For selecting top 250 up regulated and top 250 down regulated genes in
# increasing or decreasing order of their logFC, dplyr package can be used.

library(dplyr)
library(writexl)

# Selected columns like ID, GB_ACC, GENE_SYMBOL, P.Value, logFC and GENE_TYPE

data <- select(filtered_Data, ID, GB_ACC, GENE_SYMBOL, P.Value, logFC, GENE_TYPE)

table(data$GENE_TYPE) # This give no. of up and down regulated genes.

up <- arrange(data, desc(data$logFC))
down <- arrange(data, data$logFC)

# Selected only Up Regulated Gene Type

upR <- filter(up, GENE_TYPE == "Up Regulated")

# Selected only Down Regulated Gene Type

downR <- filter(down, GENE_TYPE == "Down Regulated")

# The P.Value of these upR and downR data frame is already in increasing order
# Now only top 250 genes of each were selected and stored.

upR <- upR[1:250,] # Top 250 Genes
downR <- downR[1:250,] # Top 250 Genes

# Write these data into xlsx format

write_xlsx(upR, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/upR.xlsx")
write_xlsx(downR, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/downR.xlsx")

#------------------------------------------------------------------------------#