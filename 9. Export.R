##  EXPORTING DATA FRAME TO EXCEL  ##

# The writexl is a package that contains a function write_xlsx() function 
# which is used to write a data frame to an Excel (.xlsx) file.

# Package loaded
library("writexl") 

# Saved the data frame at the specified path
write_xlsx(full_results1,"C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/full_results1.xlsx")
write_xlsx(full_resultsX,"C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/full_resultsX.xlsx")

#------------------------------------------------------------------------------#