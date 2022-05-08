## UPREGULATED OR DOWNREGULATED HAVING GENE_SYMBOLS ##

# Filtering only those top genes having Gene Symbols with GB_ACC

# Filter only those rows with no blanks or empty blanks
# Especially for GENE_SYMBOL and GB_ACC

filtered_Data <- full_resultsX[!(!is.na(full_resultsX$GENE_SYMBOL) & full_resultsX$GENE_SYMBOL == ""), ]
dim(filtered_Data)
write_xlsx(filtered_Data,"C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/filtered_Data.xlsx")

length(filtered_Data$GENE_SYMBOL)  # So total of 3077 gene symbols were obtained

