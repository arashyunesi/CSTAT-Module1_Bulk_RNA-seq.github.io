library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


# load the data
# option 1 using Read10X from Seurat, use this option if each sample is in a
# separate folder and not mixed with other samples
counts_data <- Read10X("~/CSTAT/Workshops/Workshop3_Spatial_RNA-seq/Data/2D/outs/filtered_feature_bc_matrix/")

# option 2 using ReadMtx from Seurat
# use this option if you need to specify the counts and barcodes inside a folder shared
# with other samples
Data_path <- "~/CSTAT/Workshops/Workshop3_Spatial_RNA-seq/Data/2D/outs"
counts_data <- ReadMtx(mtx = file.path(Data_path, "filtered_feature_bc_matrix/matrix.mtx.gz"),
                  cells = file.path(Data_path, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
                  features = file.path(Data_path, "filtered_feature_bc_matrix/features.tsv.gz"))

# what is dgCMatrix format?
# take a quick look at this compressed matrix format
# I explained this in detail in the second workshop of this series.

# create Seurat object from this data
sample_obj <- CreateSeuratObject(counts = counts_data, project = "Spatial_Visium",
                                 assay = "Spatial")
DefaultAssay(sample_obj) <- "Spatial"
sample_obj$slice <- 1
sample_obj$region <- "2D"

# load the high resolution image and attach it to the Seurat object
image_data <- Read10X_Image(file.path(Data_path, "spatial"),
                            image.name = "tissue_lowres_image.png")
image_data <- image_data[colnames(sample_obj)] # put the spots in the same order as the RNA-seq barcodes
sample_obj[["image"]] <- image_data # store the image into the Seurat object

SpatialFeaturePlot(sample_obj, features = "nCount_Spatial")
SpatialFeaturePlot(sample_obj, features = "nFeature_Spatial")
VlnPlot(sample_obj, features = "nCount_Spatial")

# we can put a variable transparency that increases with the nCount to see the tissue in the background
SpatialFeaturePlot(sample_obj, features = "nCount_Spatial", alpha = c(0.2,1))

# Let us visualize the TLS structures
# B cells
SpatialFeaturePlot(sample_obj, features = c("CD19", "CD22", "CD79A", "CD79B"), slot = "counts")

# T cells
SpatialFeaturePlot(sample_obj, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"), slot = "counts", alpha = c(0.3,1))


plot(rowMeans(counts_data), apply(counts_data, 1, sd),
     xlim = c(0,200), ylim = c(0,200))
# There is strong relationship between mean and sd for spatial data. In case of Bulk RNA-seq and scRNA-seq we transformed the data (log2, vst, and LogNormalize)
# In scRNA-seq and Bulk RNA-seq, we performed data transformation to correct for technical variance between observations (cells or samples), i.e. library size differences
# In Spatial RNA-seq data the variance between different spots (observations) depends on technical variation plus the variance in
# tissue structure. There are hotspots e.g. here due to immune system activity around cancer regions
# Look at the TLS structures in the sample
# It is common to perform SCTransform on spatial data
# refer to https://link.springer.com/article/10.1186/s13059-019-1874-1 for details of SCTransform
# In a nutshell: builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance.

# for SCTransform to run faster, install glmGamPoi package, otherwise it will run very slow
# BiocManager::install('glmGamPoi')
sample_obj <- SCTransform(sample_obj, assay = "Spatial", verbose = TRUE)





########################################################################
# Dimensionality reduction, clustering, and visualization
########################################################################
# run Principal Component Analysis on the SCT data
sample_obj <- RunPCA(sample_obj, assay = "SCT", verbose = TRUE) # default number of calcuated PCs is 50, can be set higher or lower

# how many PCs should we use for downstream analysis?
plot(1:50, sample_obj@reductions$pca@stdev) # 15 to 20 PCs are enough

# perform kNN on the PCA embedded data
sample_obj <- FindNeighbors(sample_obj, reduction = "pca", dims = 1:20)

# Run clustering analysis on the graph of nearest neighbors from the previous step
sample_obj <- FindClusters(sample_obj, verbose = TRUE, resolution = 0.5, random.seed = 1223)
# Let's discuss how to choose the resolution here.
# Changing the resolution might change the clustering pattern.
# Large resolution usually finds more clusters and small resolution finds less clusters.
# A trick to decide the resolution would be using Clustering Stability
# e.g. change the resolution around the selected value and compare the new clustering with the original clustering using ARI.
# We can work on this if there is enough time.


# run UMAP to visualize the clusters
sample_obj <- RunUMAP(sample_obj, reduction = "pca", dims = 1:20)
DimPlot(sample_obj, reduction = "umap", label = TRUE)
SpatialDimPlot(sample_obj, label = TRUE, label.size = 3) # set the  pt.size.factor = 1 to see the tissue in the background


# The plot above is too messy and difficult to see. Can we separate the clusters in the plots?
SpatialDimPlot(sample_obj, cells.highlight = CellsByIdentities(object = sample_obj, idents = 0:10), facet.highlight = TRUE, ncol = 3)

# Interactive plots would be very helpful
SpatialDimPlot(sample_obj, label = TRUE, label.size = 3, interactive = TRUE)
# Interactive plots over the genes
SpatialFeaturePlot(sample_obj, interactive = TRUE, features = "CD79A")



###############################################################################################
# Finding Spatially Variable Genes
###############################################################################################
# Option 1:
# What did we do in the case scRNA-seq? We found the genes that were statistically different in expression between the clusters.
# We can do the same type of finding variable genes between the clusters, this time using the spatial clustering performed.
Markers_Clusters <- FindAllMarkers(sample_obj, logfc.threshold = 1, min.pct = 0.05, min.diff.pct = 0.1)
# for faster implementation of finding markers, install.packages("Rfast)

Markers_Clusters <- Markers_Clusters %>%
  filter(p_val_adj < 0.01)


# What are the issues with this type of finding spatially variable genes?
# 1) It is reliant on the clustering performed, not very good clustering would give ...
# 2) Some gene patterns might not follow the cluster patterns exactly
# 3) The length scale of change in some genes might be larger or smaller than the cluster sizes


# Option 2:
# Use built-in Spatial Variable Gene selection method, Seurat uses Moran's I
sample_obj <- FindSpatiallyVariableFeatures(sample_obj, assay = "SCT", selection.method = "moransi",
                                              features = VariableFeatures(sample_obj))

View(sample_obj@assays[["SCT"]]@meta.features)
# As you can see this method is very slow, even for a subset of genes (3000 genes from the original 30k genes).
# Moran's I is a very primitive method of finding spatial patterns in data.
# Does not take into account the nature of RNA-seq data.
# Does not distinguish different length scales of spatial variation.
# In this example has produced NAs for most of the 3000 genes that we tested.


# Option 3: Use other packages for SVG detection
###########################################################################
# SPARK and SPARK-X for Spatially Variable Gene selection
###########################################################################
# if you don't have devtools first install.packages('devtools')
devtools::install_github('xzhoulab/SPARK')


library(SPARK)

# First let us read the 
Locations <- read.table(file = file.path(Data_path, "spatial/tissue_positions_list.csv"),
                        sep = ",", row.names = 1)
# What are the different columns of the locations data?
# look at the locations, there are 4992 spots
# look at the sample_obj, how many spots are there? 3963 spots only! Why the discrepancy?
# remove the extra spots from the locations data.
spots <- intersect(rownames(Locations), colnames(sample_obj))
counts_data <- counts_data[, spots] # putting the counts data at the same order with the spot names and locations
Locations <- Locations[spots,] # removing the extra locations and setting the order of the spots same as the counts data

# create a SPARK object from the counts and locations data
SPARK_obj <- CreateSPARKObject(counts_data, data.frame(Locations[,4:5]), min_total_counts = 10, percentage = 0.05)
SPARK_obj@lib_size <- apply(SPARK_obj@counts, 2, sum)

# the theoretical model of SPARK and SPARKX in a nutshell: https://xzhoulab.github.io/SPARK/
# fit a gaussian model to the counts data
SPARK_obj <- spark.vc(SPARK_obj, covariates = NULL, lib_size = SPARK_obj@lib_size,
                      num_core = parallelly::availableCores(), fit.model = "gaussian", verbose = F)
# note that here I am using the parallelization option in SPARK for model fitting. I use the "parallelly" package to find the cpu cores
# that I have reserved. Do not use DetectCores from the parallel package on a shared system!

SPARK_obj <- spark.test(SPARK_obj, check_positive = T, verbose = F) # testing function is not parallelized, so it is slow.
SPARK_analysis <- SPARK_obj@res_mtest
# SPARK uses Cauchy Combination rule to combine p-values from different tests on the same gene
# Multiple testing adjustment (adjusted p-value) is produced using Benjamini-Yekutieli (BY) procedure
# BY is a more conservative procedure than Benjamini-Hochberg (BH), e.g. stronger Type-I error control rate





###############################################################
# Subset anatomical regions
###############################################################
# Let us subset the two TLS structures in the bottom left corner
# We can use coordinates to subset the data
# To find the coordinates, we can use an interactive plot

SpatialDimPlot(sample_obj, interactive = TRUE, pt.size.factor = 0.5)
# on this plot, find the coordinates of a box containing the two TLS 
# c(3852, 5015) c(5011, 4131) are the opposite corners

TLS_locations <- Locations %>%
  filter(V5>3852 & V5<5011) %>%
  filter(V6>4131 & V6<5015)

# first, we add a new metadata to include the spot names
sample_obj$spot_barcode <- colnames(sample_obj)
# then subset the Seurat object using the barcodes of the TLS hotspots
TLS_obj <- subset(sample_obj, subset = spot_barcode %in% rownames(TLS_locations))

SpatialDimPlot(TLS_obj, crop = TRUE, label = TRUE, pt.size.factor = 4, label.size = 3)
SpatialDimPlot(TLS_obj, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 0)




###############################################################
# Integration of Spatial Visium Data with scRNA-seq data
###############################################################
library(SeuratData) # install with "devtools::install_github('satijalab/seurat-data')"


# find an available reference dataset from SeuratData package
AvailableData()
# if they don't have a reference dataset for the organ you need, you can download an annotated dataset online, e.g. from Gene Expression Omnibus

# if it is your first time using the reference dataset, you need to install it first. This only applies to datasets provided by Seurat
InstallData("lungref")
# check the installation
InstalledData()

# load the refrence dataset
LoadData("lungref", "azimuth")

# I could not make the lungref from SeuratData work! After an hour of trying, I stopped.
# I decided to go to Human Lung Cell Atlas and download their Core dataset (500k cells annotated from many donors and studies)
# BiocManager::install("anndataR") if you do not have the package
# library(anndataR)
# 
# HLCA_core <- read_h5ad("~/CSTAT/Workshops/Workshop3_Spatial_RNA-seq/Data/HLCA_core.h5ad", as = "InMemoryAnnData")
# 
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# 
# library(SeuratDisk)
# Convert("~/CSTAT/Workshops/Workshop3_Spatial_RNA-seq/Data/HLCA_core.h5ad", dest = "h5seurat", overwrite = TRUE)
# HLCA_core <- LoadH5Seurat("~/CSTAT/Workshops/Workshop3_Spatial_RNA-seq/Data/HLCA_core.h5seurat")



####################################################
# not working!! Let's change gears.
####################################################
# We follow Seurat's vignettes to process mouse brain Visium and annotation procedure
# You can see the vignette used for this integration at:
# https://satijalab.org/seurat/articles/spatial_vignette

InstallData("stxBrain")

brain <- LoadData("stxBrain", type = "anterior1")

# preprocessing steps
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = TRUE)

# visualizations
plot <- SpatialFeaturePlot(brain, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
                                                               legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
# if you want to save the plot use the next lines:
# jpeg(filename = "../output/images/spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
# print(plot)
# dev.off()

# dimensional reduction, clustering, ...
brain <- RunPCA(brain, assay = "SCT", verbose = TRUE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = TRUE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain,idents = c(2, 1, 4, 3, 5, 8)),
               facet.highlight = TRUE, ncol = 3)


# Integration with scRNA-seq data
allen_reference <- readRDS("~/CSTAT/Workshops/Workshop3_Spatial_RNA-seq/Data/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = TRUE) %>%
  RunPCA(verbose = TRUE) %>%
  RunUMAP(dims = 1:30)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)


# Subset the mouse brain sample to clusters of interest
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = TRUE) %>%
  RunPCA(verbose = TRUE)


anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4", "Astro"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)



