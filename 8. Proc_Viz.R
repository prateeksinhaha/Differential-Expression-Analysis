## PROCESSING AND VISUALIZATION OF DIFFERENTIAL EXPRESSION RESULTS ##

# The gene name associated with the gene ID can be found out. The annotation 
# data can be retrieved with the 'fData' function. The ID, GB_ACC (GenBank 
# Accession ID) needs to be selected. And then can be added into fit2 table.

# The "Volcano Plot" function can be used to visualize the results of a DE 
# analysis. The x axis shows the log-fold change and the y axis is some measure 
# of statistical significance, which in this case is the log-odds, or 
# "B" statistic. The color of those genes with p value cutoff more than 0.05, 
# and fold change cut off more than 1 can be changed.

library(dplyr)

annotation <- fData(GEOSeries)
head(annotation)

annotation <- dplyr::select(annotation, ID, GB_ACC, GENE_SYMBOL, GENE_NAME, SEQUENCE)  # To select only ID and GB_ACC
head(annotation)


fit2$genes <- annotation # Assigning the annotation data to $genes 
# (newly created) in fit2

topTable(fit2)

full_results1 <- topTable(fit2, coef=1, number=Inf)

# The Volcano Plot can now be created

library(ggplot2)

png(filename = "Volcano Plot (logFC vs B).png"
    ,res=720, width = 12, height = 8, units = "in")

ggplot(full_results1, 
       aes(x = logFC, y=B), ) + geom_point() + 
  labs(title = "Volcano Plot (logFC vs B)") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12))

dev.off()

# Created a simple volcano plot 

png(filename = "Volcano Plot (logFC vs adj.P.Val).png"
    ,res=720, width = 12, height = 8, units = "in")

vol_plot <- full_results1 %>%
  ggplot(aes(x = logFC,
             y = -log10(adj.P.Val))) + 
  geom_point() + labs(title = "Volcano Plot (logFC vs adj.P.Val)")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12)) 
vol_plot
dev.off()

# The flexibility of ggplot2 allows to automatically label points on the plot 
# that might be of interest. Example: genes that meet an adjusted p-value 
# and log fold-change cut-off. 

# To label genes into the following groups:

# 1. Genes with logFC > 1 and adj.P.Val < 0.05 labelled as Highly Up Regulated.
# 2. Genes with logFC < -1 and adj.P.Val < 0.05 labelled as Highly Down Regulated.
# 3. All other genes labelled as Non-Significant i.e. non-significant.

# Create new categorical column 

full_results1 <- full_results1 %>%
  mutate(GENE_TYPE = case_when(logFC >= 1 & adj.P.Val <= 0.05 ~ "Highly Up Regulated",
                               logFC <= -1 & adj.P.Val <= 0.05 ~ "Highly Down Regulated",
                               TRUE ~ "Non-Significant"))  

# Add colour, size to volcano plot 

full_results1 %>%
  ggplot(aes(x = logFC,
             y = -log10(adj.P.Val),
             fill = GENE_TYPE,    
             size = GENE_TYPE)) + 
  labs(title = "Volcano Plot (logFC vs adj.P.Val) With Gene Type") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        legend.position = "bottom") +
  geom_point(shape = 21, # Specify shape and color as fixed local parameters    
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point color
  scale_size_manual(values = sizes) + # Modify point size
  scale_x_continuous(breaks = c(seq(-2.5, 2.5, 0.5)),       
                     limits = c(-2.5, 2.5))
dev.off()


# Furthermore, the identity of some genes can be labelled. A limit can be set
# of the top "N" genes that has to be labelled, and label each gene according 
# to it's Symbol.

library(ggrepel)
#p_cutoff <- 0.05
#fc_cutoff <- 0
#topN <- 20

#full_results1 %>% 
#     mutate(Significant = adj.P.Val < p_cutoff, logFC > fc_cutoff ) %>% 
#    mutate(Rank = 1:n(), Label = ifelse(Rank < topN, GENE_SYMBOL,"")) %>% 
#   ggplot(aes(x = logFC, y = -log10(adj.P.Val), col=Significant,label=Label)) + 
#  geom_point() + geom_text_repel(col="black") +
# labs(title = "Volcano Plot (logFC vs adj.P.Val) With Top 20 Genes")


## FILTERING AND EXPORTING THE RESULT TABLE

# The filter function from dplyr gives a convenient way to interrogate the 
# table of results.

# To look into the fold change data of a selected gene, whether it is 
# significantly differential expressed or not:

# The results for particular gene of interest were obtained
# GB_ACC for Nkx3-1 is NM_001256339 or NM_006167
# no NM_001256339 in this data


# filter(full_results1, GB_ACC == "NM_006167") #imp


# New result is formed which satisfies the cutoff

full_results1 <- full_results1 %>%
  mutate(GENE_TYPE = case_when(logFC > 0 & adj.P.Val <= 0.05 ~ "Up Regulated",
                               logFC < 0 & adj.P.Val <= 0.05 ~ "Down Regulated",
                               TRUE ~ "Non-Significant"))  



p_cutoff <- 0.05
fc_cutoff <- 0

full_resultsX <- full_results1 %>% 
  mutate(Significant = adj.P.Val < p_cutoff)

head(full_resultsX)

table(full_resultsX$Significant)

# A new colunm is formed named Significant 
# which have binary values of either TRUE or FALSE in each 559136 genes.

# TRUE Values (97457 genes) are that genes which satisfies logFC cut i.e. > 0
# which could be either Up regulated or Down regulated


# Therefore, it can be seen that 97457 genes satisfies the cutoff.
# These genes can be filtered out as well and stored.

full_resultsX <- filter(full_resultsX, Significant==TRUE)
dim(full_resultsX)

#------------------------------------------------------------------------------#





