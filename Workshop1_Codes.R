library(DESeq2)
library(tidyverse)


####################################################################################
# Reading data
####################################################################################
Raw_Counts <- readRDS(file = "~/CSTAT/Workshops/Workshop1_Bulk_RNA_seq/Workshop1_Bulk_RNA_seq/Raw_counts.rds")

Sample_Info <- read.table(file = "~/CSTAT/Workshops/Workshop1_Bulk_RNA_seq/Workshop1_Bulk_RNA_seq/samplesheet.tsv",
                          sep = "\t", header = TRUE)


# Is the order of the samples in the Sample_Info file the same as the order of the samples in the data?
# Check visually or use code
colnames(Raw_Counts)
Sample_Info$SampleName

# use R to check each one
colnames(Raw_Counts)==Sample_Info$SampleName
# use R to check all
all(colnames(Raw_Counts)==Sample_Info$SampleName)

####################################################################################
# Filtering out low quality genes
# We keep genes with a threshold of total expression, and non-zero expression in at least three samples.
# Reducing number of genes that are being tested will also help with multiple testing correction.We will discuss multiple testing later.
####################################################################################
rowSums(Raw_Counts)

# filter the genes with total measured counts of less than 50
keep <- rowSums(Raw_Counts) > 5
table(keep)
Raw_Counts_filtered <- Raw_Counts[keep,]

# filter the genes that are non-zero in less than 3 samples (to be able to calculate SD)
keep <- rowSums(Raw_Counts_filtered > 0) >= 3
table(keep)
Raw_Counts_filtered <- Raw_Counts_filtered[keep,]



####################################################################################
# Data exploration and analysis
####################################################################################
summary(Raw_Counts_filtered)

# Problem 1) Outlier measurements
# few outliers affect distribution visualization
boxplot(Raw_Counts_filtered, main='Raw counts', las=2)

# Problem 2) Dispersion
# Raw counts mean expression Vs standard Deviation (SD)
plot(rowMeans(Raw_Counts_filtered), rowSds(Raw_Counts_filtered), 
     main='Raw Counts: SD vs Mean')

plot(rowMeans(Raw_Counts_filtered), rowSds(Raw_Counts_filtered), 
     main='Raw Counts: SD vs Mean', 
     xlim=c(0,10000),
     ylim=c(0,5000))

# Problem 3) Different library sizes
# the total number of reads per sample (library size), highly influenced by the sample preparation and PCR rate
colSums(Raw_Counts_filtered)
# 24M vs 34M, it is a 40% difference in total number of reads
# Question: Will the total number of reads affect the measured abundance of a gene? YES.

####################################################################################
# Data transformation:
# Reduce the effect of outlier genes
# Reduce dependence of variance to mean
# Reduce dependence of a gene's expression to the library size
# log2, VST (Variance Stablizing Transformation), rlog (regularized log)
# https://link.springer.com/article/10.1186/s13059-019-1874-1#:~:text=The%20primary%20goal,shallowly%20sequenced%20cells.
####################################################################################

# log2 transformation
logcounts <- log2(Raw_Counts_filtered + 1)
# adding 1 to the counts, they are called pseudo-counts in the business

# boxplot of the log transformed counts: No Outliers anymore
boxplot(logcounts, main='Log Transformed Raw Counts', las=2)

# Log counts Mean expression vs Standard Deviation (SD): Much less dispersion
plot(rowMeans(logcounts), rowSds(logcounts), 
     main='Log2 Counts: SD vs Mean')
# looks pretty good and the strong dependence of the SD to Mean expression has significantly reduced

# how does the transformed library size look?
colSums(logcounts)
# The difference is now 144k vs 152k, only 6% difference between the largest and smallest


# Let's style the boxplot of the log2 transformed data to look better
# make a colour vector
statusCols <- case_when(Sample_Info$Status=="Infected" ~ "red", 
                        Sample_Info$Status=="Uninfected" ~ "green")

# Check distributions of samples using boxplots
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts+1)",
        las=2,
        col=statusCols,
        main="Log2(Counts+1)")
# Let's add a blue horizontal line that corresponds to the overall median
abline(h=median(logcounts), col="blue")



# Variance Stablizing Transformation from DESeq2 (wrapper for VST from sctransform)
# developed for scRNA-seq, but it works well for the case of Bulk RNA-seq as well.
# VST fits a regularized Negative Binomial Model to a subset of the genes to determine the parameters of the transformation.
# Pearson residuals of the fitted model are used.
vst_counts <- vst(Raw_Counts_filtered)

# Check distributions of samples using boxplots
boxplot(vst_counts, 
        xlab="", 
        ylab="VST counts",
        las=2,
        col=statusCols,
        main="Variance Stabilizing Transformation")
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(vst_counts), col="blue")

# VST counts standard deviation (sd) vs mean expression
plot(rowMeans(vst_counts), rowSds(vst_counts), 
     main='VST counts: SD vs Mean')
# it looks even flatter than log2 transformed counts

# how does the transformed library size look?
colSums(vst_counts)
# The difference is less than 1%! Pretty good.


############
###*****##########
# Exercise: use the rlog function from the package DESeq2 to transform the data. Check the SD vs Mean and boxplots of the data after transformation.
# Give them 5-10 minutes and walk around to help them with this.
###*****##########
############

rlog_counts <- rlog(Raw_Counts_filtered)


# Check distributions of samples using boxplots
boxplot(rlog_counts, 
        xlab="", 
        ylab="rlog counts",
        las=2,
        col=statusCols,
        main="Variance Stabilizing Transformation")
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(rlog_counts), col="blue")

# rlog counts standard deviation (sd) vs mean expression
plot(rowMeans(rlog_counts), rowSds(rlog_counts), 
     main='rlog counts: SD vs Mean')

# how does the transformed library size look?
colSums(rlog_counts)
# The difference is less than 1%! Pretty good.



####################################################################################
# Dimensional Reduction using PCA:
# Explain PCA from 3 dimensions going to 2 dimensions. Then from 20k dimensions to 4.
# Explain scree plot and how to choose the necessary number of dimensions.
####################################################################################
library(ggfortify) # install.packages("ggfortify")

dim(rlog_counts)

# run PCA
pcDat <- prcomp(t(rlog_counts))
# Notes: centering the data is performed by default, but scaling is advised and better performed. However, in this data the scaling is making it worse!
# The same applies when using VST normalized and log2 normalized values.
# Scaling the data, possibly puts the low expressed genes and high expressed genes at the same level with each other. That means, low expressed genes
# will have the same effect as highly expressed genes. 
# https://www.biostars.org/p/442601/



# plot principal components (commonly known as a scree plot)
plot(1:12, pcDat$sdev, type = "b")
# it looks like 3 or 4 PCs are enough

# plot PCA
autoplot(pcDat)

# let's make a better looking plot of PCA
autoplot(pcDat,
         data = Sample_Info, 
         colour="Status", 
         shape="TimePoint",
         size=5)


# some notes on PCA:
summary(pcDat)
# to calculate the percentage of variances, divide the square of the sdev to sum of all sdev
pcDat[["sdev"]]^2/sum(pcDat$sdev^2)

############
###*****##########
# Exercise: Make autoplot of 1st and 3rd PCs. Make autoplot of 2nd and 3rd PCs.
# Give them 5-10 minutes and walk around to help them with this.
# Consult the help page for autoplot.prcomp or type ?autoplot.prcomp in the console
###*****##########
############
autoplot(pcDat,
         x=1, y=3,
         data=Sample_Info, 
         colour="Status", 
         shape="TimePoint",
         size=5)

autoplot(pcDat,
         x=2, y=3,
         data=Sample_Info, 
         colour="Status", 
         shape="TimePoint",
         size=5)



# visualize pairs of PCs simultaneously (can apply pairs() from the base or use ggpairs from GGally)
library(GGally) # install.packages("GGally")

PCA_df <- data.frame(pcDat$x[,1:4])
PCA_df$Status <- Sample_Info$Status
PCA_df$TimePoint <- Sample_Info$TimePoint
ggpairs(PCA_df, aes(color=Status, shape=TimePoint))



####################################################################################
# Hierachical clustering
####################################################################################
library(ggdendro) # install.packages("ggdendro")

hclDat <-  t(rlog_counts) %>%
  dist(method = "euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate=TRUE)
# This clustering is using the entire 20k genes. So the clustering is performed on the high dimensional data.
# Generally it is best to perform HC on the dimensional reduced data.
# However, for clean data sets that might not be necessary as we are about to see in the following.

# let's add better labels to the plot
hclDat2 <- hclDat # save the original and change a copy
hclDat2$labels <- paste0(Sample_Info$Status, ":", Sample_Info$TimePoint)
ggdendrogram(hclDat2, rotate=TRUE)


hclDat <-  pcDat$x[,1:4] %>%
  dist(method = "euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate=TRUE)
# This clustering is using the first 4 Principal Directions (distances measured in the reduced dimensional space) 

# let's add better labels to the plot
hclDat2 <- hclDat # save the original and change a copy
hclDat2$labels <- paste0(Sample_Info$Status, ":", Sample_Info$TimePoint)
ggdendrogram(hclDat2, rotate=TRUE)

# Since this is a very clean and high quality data, using the hierarchical clustering using the original data and the PCA gave the
# same result. When using general data, use the PC embedding. Originql data might have too much noise (biological or technical noise or both)




#######################################################################################
# DEG Analysis using DESeq2 package
#######################################################################################

# First we save the column data and transform the columns to factors
# We can choose the levels so that Uninfected is first level and the times are in chronological order (for later convenience)
column_Data <- data.frame(Status=factor(Sample_Info$Status, levels = c("Uninfected", "Infected")),
                          TimePoint=factor(Sample_Info$TimePoint, levels = c("d11", "d33")))

rownames(column_Data) <- Sample_Info$SampleName

# We also need a design model for the fixed effects
# As a simple model, first assume that the TimePoint variable has no effect on the gene expressions
# We design a model for Status only
model_Status <- as.formula(~ Status)

# This model as a matrix looks like the following:
model.matrix(model_Status, data = column_Data)


# Create DESeq Object from the counts matrix
DESeq_Obj <- DESeqDataSetFromMatrix(Raw_Counts_filtered,
                                    colData = column_Data,
                                    design = model_Status)


# DESeq command does three tasks:
# 1) Estimates Size Factors: Calculate the “median ratio” normalisation size factors for each sample and adjust for average transcript length on a per gene per sample basis.
### DESeq2 has calculated a normalizsation factor for each gene for each sample.
# 2) Estimates Dispersions: Estimates the approximate relationship between Mean and Variance for each gene.
# 3) Applies Negative Binomial GLM fitting and calculate Wald statistics.

DESeq_Obj <- DESeq(DESeq_Obj)

results_Model1 <- results(DESeq_Obj, alpha=0.05)
results_Model1

table(results_Model1$padj < 0.05)

# Histograms of p-values and adj p-values
hist(results_Model1$pvalue)
hist(results_Model1$padj)

# Ordering the results by the adj p-value
Model1_df <- data.frame(results_Model1@listData)
rownames(Model1_df) <- rownames(results_Model1)
Model1_df <- Model1_df %>% arrange(padj)


# log2FC = log2(Treatment/Control)
# Exercise: How many genes are up-regulated and padj<0.05? How many are down-regulated and padj<0.05?
# Exercise: Repeat the above analysis for a model based on TimePoint only, ignoring the Status.






# Both + Interaction model
model_Interaction <- as.formula(~ Status + TimePoint + Status:TimePoint)
# Create DESeq Object from the counts matrix
DESeq_Obj_Int <- DESeqDataSetFromMatrix(Raw_Counts_filtered,
                                    colData = column_Data,
                                    design = model_Interaction)


DESeq_Obj_Int <- DESeq(DESeq_Obj_Int)

# contrasts
resultsNames(DESeq_Obj_Int)


results_Model_Int <- results(DESeq_Obj_Int, alpha=0.05)
results_Model_Int

table(results_Model_Int$padj < 0.05)

# Histograms of p-values and adj p-values
hist(results_Model_Int$pvalue)
hist(results_Model_Int$padj)

# Ordering the results by the adj p-value
Model_Int_df <- data.frame(results_Model_Int@listData)
rownames(Model_Int_df) <- rownames(results_Model_Int)
Model_Int_df <- Model_Int_df %>% arrange(padj)




###############################################################################
# Gene Set Enrichment Analysis
###############################################################################
library(gprofiler2)


Sig_Genes_Int <- Model_Int_df %>%
  filter(padj<0.05) %>%
  rownames()

GSEA_Int <- gost(query = Sig_Genes_Int, organism = "mmusculus", ordered_query = F, significant = F,
                             multi_query = F, measure_underrepresentation = F, evcodes = F,
                             domain_scope = "custom", custom_bg = rownames(Raw_Counts_filtered), user_threshold = 0.05,
                             correction_method = "g_SCS", numeric_ns = "", as_short_link = F)
GSEA_Int$result <- GSEA_Int$result %>% arrange(p_value)
sig_pathways_Int <- GSEA_Int$result[c(1:7, 9,10),]
sig_pathways_Int$log10p <- -log10(sig_pathways_Int$p_value)
sig_pathways_Int$Title <- "GSEA Interaction Model"
sig_pathways_Int <- data.frame(sig_pathways_Int)
GSEA_plot(sig_pathways_Int)


######################################
GSEA_plot <- function(datt){
  p <- ggplot(data=datt, aes(x=forcats::fct_reorder(term_name, log10p), y = log10p, fill=source, label = ifelse(significant == TRUE, "*", ""))) +
    geom_bar(stat="identity", position = position_dodge(), color = NA, width = 0.8) +
    scale_fill_manual(values = c("#f6bd60", "#f5cac3", "#bfcc94", "#f28482", "#a1234f")) +
    labs(title=paste0(datt$Title), x = paste0(""), y = "-log10(p value)") +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 2) +
    geom_text(vjust = 1, nudge_y = 0.5, size = 15) +
    theme(plot.title = element_text(size = 30, color = "black"),
          text = element_text(size = 30, color = "black"),
          axis.text.x=element_text(size = 30, color = "black") ,
          legend.position = "none")
  return(p)
}
#############################################
