## INSPECTING THE CLINICAL VARIABLES ##


# GEO submitted data contain sample labels and information about processing 
# protocol which can be extracted by pData function. For a given data, it has to
# be decided which columns will be useful in the analysis. Example: In the given
# dataset characteristics_ch1 is useful. The columns can be renamed and 
# selected by using "rename" and "select" function using dplyr package.

library(dplyr)  

SampleInfo <- pData(GEOSeries)  # Assigned a variable to phenotypic data of 
                                # GEOSeries

head(SampleInfo)                # To display first 6 rows (by default)

table(SampleInfo$characteristics_ch1)  # $ - extracts or subsets a specific 
                                       # part of a data object. 

# table() performs categorical tabulation of data with the variable and its 
# frequency

# Columns having factors needed for analysis were selected

SampleInfo <- select(SampleInfo, characteristics_ch1)

# Columns can be renamed as follows 

SampleInfo <- rename(SampleInfo, sample = characteristics_ch1)
head(SampleInfo)
dim(SampleInfo)                     # Displayed the dimensions of an object
SampleInfo$sample


# "stringr" is used to make a column of simplified group names for each sample. 
# A new column is created, named "group". The function 
# "str_detect" detects the presence of the words, and then fills the row 
# accordingly. It depends on dataset to make these categories in the new 
# columns. These commands can be modified for the dataset of interest.

library(stringr)  
SampleInfo$group <- ""
for(i in 1:nrow(SampleInfo)){
        if(str_detect(SampleInfo$sample[i], "disease state: AD"))
        {SampleInfo$group[i] <- "AD"}
        
        if(str_detect(SampleInfo$sample[i], "disease state: Control"))
        {SampleInfo$group[i] <- "CTRL"}
        
}

SampleInfo
View(SampleInfo)
y <- head(SampleInfo)
write.csv(y, "C:/Users/Prateek/R/GEOData/GSE138260/Differential-Expression-Analysis/SampleInfoHead.csv")
#------------------------------------------------------------------------------#
