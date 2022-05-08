## DISTRIBUTION OF EXPRESSION LEVELS USING BOX PLOT ##


# Grouped membership for all samples

gsms <- "000000000000000001111111111111111111"
sml <- strsplit(gsms, split="")[[1]]

# Assigned samples to groups and set up design matrix

gs <- factor(sml)
groups <- make.names(c("AD","CTRL"))
levels(gs) <- groups
GEOSeries$group <- gs
dev.new(width = 3 + ncol(GEOSeries)/6, height=5)

# Box Plot

ord <- order(gs)  # Ordered samples by group

png(filename = "Distribution of Expression Levels (GSE138260).png"
    ,res=720, width = 8, height = 8, units = "in")
palette(c("#66A61E", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#1B9E77", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))

title <- paste("Distribution of Expression Levels (GSE138260)")
boxplot(exprs(GEOSeries)[,ord], boxwex = 0.6, outline = FALSE, main = title, 
        ylab = "Gene Expression of Each Sample", las = 2, col=gs[ord])  

# Boxplot was plotted
# Orientation can be changed by 
# changing las=2

legend("topright", groups, fill=palette(), bty="n")

dev.off()


#------------------------------------------------------------------------------#