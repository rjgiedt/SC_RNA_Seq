# Clear Environment
closeAllConnections()
rm(list=ls())

# Load Libraries 
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)

# Load dataset (10X reads in ) 
test.data <- Read10X("/Users/jimmy/Documents/Large Data/SC RNA-Seq-Mouse/Rawdata/FC_03294/Unaligned_10xc_PF_TenX_mm1/cellranger_count/Garren1_ABH_aPD1_SI-GA-A3/outs/raw_gene_bc_matrices/mm10/")

#Examine the object size
sparse.size <- object.size(x = pbmc.data)
sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
test <- CreateSeuratObject(raw.data = test.data, min.cells = 3, min.genes = 200, project = "10X_test")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = test@data), value = TRUE)
percent.mito <- Matrix::colSums(test@raw.data[mito.genes, ])/Matrix::colSums(test@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
test <- AddMetaData(object = test, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = test, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = test, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene")

#### Normalizing Data
test <- NormalizeData(object = test, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Detection of variable genes across the single cells
test <- FindVariableGenes(object = test, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = test@var.genes)

test <- ScaleData(object = test, vars.to.regress = c("nUMI", "percent.mito"))

test <- RunPCA(object = test, pc.genes = test@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

PrintPCA(object = test, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = test, pcs.use = 1:2)
PCAPlot(object = test, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
test <- ProjectPCA(object = test, do.print = FALSE)

PCHeatmap(object = test, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = test, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

# Cluster the cells
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
test <- FindClusters(object = test, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = test)

test <- RunTSNE(object = test, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = test)