# Clear Environment
closeAllConnections()
rm(list=ls())

# Load Libraries 
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)

# Load dataset (10X reads in ) 
test.data <- Read10X("/Users/jimmy/Documents/Large Data/SC RNA-Seq-Mouse/Rawdata/FC_03294/Unaligned_10xc_PF_TenX_mm1/cellranger_count/Garren1_Control_SI-GA-A1/outs/filtered_gene_bc_matrices/mm10")

# see https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/dgTMatrix-class.html
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes

dim(test.data)
# 27998 genes and 4390 single cells

# Check out the first six genes and cells
test.data[1:6, 1:6]

# check how many genes have at least one transcript in each cell
summary(colSums(as.matrix(test.data))) # As matrix to eliminate matrix sparsity

# check how many genes have at least one transcript in each cell
at_least_one <- apply(test.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
# Median number of detected genes among the single cells is ~ 1300

hist(colSums(as.matrix(test.data)),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
# Median sum of expression among the single cells is ~4,000

# We will filter out genes and single cells before we continue with the analysis. 
# The tutorial has arbitrary values of keeping genes expressed in three or more cells 
# and keeping cells with at least 200 detected genes.

# manually check the number of genes detected in three or more cells
tmp <- apply(test.data, 1, function(x) sum(x>0))
table(tmp>=3)
# a lot of genes are not detected in 3 or more cells (14,726 vs. 13272)

# all cells have at least 200 detected genes
keep <- tmp>=3
tmp <- test.data[keep,]
at_least_one <- apply(tmp, 2, function(x) sum(x>0))
summary(at_least_one)

dim(tmp)

# Create Seurat Class
test <- CreateSeuratObject(raw.data = test.data,
                           min.cells = 3,
                           min.genes = 200,
                           project = "10X_PBMC")

# see ?seurat for more information on the class
class(test)
test
# an object of class seurat in project 10X_PBMC 
# 13714 genes across 2700 samples following filtering for min cells/ min genese.

slotNames(test)

# The tutorial states that "The number of genes and UMIs (nGene and nUMI) are 
# automatically calculated for every object by Seurat." The nUMI is calculated as 
# num.mol <- colSums(object.raw.data), i.e. each transcript is a unique molecule. 
# The number of genes is simply the tally of genes with at least 1 transcript; 
# num.genes <- colSums(object.raw.data > is.expr) where is.expr is zero.

# A common quality control metric is the percentage of transcripts from the mitochondrial 
# genome. According to the paper "Classification of low quality cells from single-cell 
# RNA-seq data" the reason this is a quality control metric is because if a single cell 
# is lysed, cytoplasmic RNA will be lost apart from the RNA that is enclosed in the 
# mitochondria, which will be retained and sequenced.

# mitochondria genes conveniently start with MT
mito.genes <- grep(pattern = "^mt-", x = rownames(x = test@data), value = TRUE)
length(mito.genes)

percent.mito <- Matrix::colSums(test@raw.data[mito.genes, ]) / Matrix::colSums(test@raw.data)

# check out the meta data
head(test@meta.data)

# add some more meta data
test <- AddMetaData(object = test,
                    metadata = percent.mito,
                    col.name = "percent.mito")

head(test@meta.data)

# plot number of genes, UMIs, and % mitochondria
VlnPlot(object = test,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)

#The GenePlot() function can be used to visualise gene-gene relationships as well as any 
# columns in the seurat object. Below we use the plotting function to spot cells that have 
# a high percentage of mitochondrial RNA and to plot the relationship between the number 
# of unique molecules and the number of genes captured.

par(mfrow = c(1, 2))
GenePlot(object = test, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')

# Next we'll use the FilterCells() function to subset the pbmc object based on the number 
# of genes detected in each cell and by the percent mitochondria. Two thresholds need are 
# specified for each filter, a low and a high; -Inf and Inf are used if you only want to 
# specify a lower or upper bound respectively. Cells that have less than 200 genes and 
# more than 2,500 genes and over 5% mitochondrial content are filtered out.

# manual check; I already know all cells have >200 genes
table(test@meta.data$percent.mito < 0.05 & test@meta.data$nGene<2500)

# perform the filtering using FilterCells()
test <- FilterCells(object = test,
                    subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf),
                    high.thresholds = c(2500, 0.05))

# No cells are filtered out; numbers consistent with above
test

# The next step is to normalise the data, so that each cell can be compared against each other. 
# At the time of writing, the only normalisation method implemented in Seurat is by log 
# normalisation. Gene expression measurements for each cell are normalised by its total 
# expression, scaled by 10,000, and log-transformed.

hist(colSums(as.matrix(test@data)),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

# currently there is only log-normalisation
# according to the documentation more methods to be included soon
test <- NormalizeData(object = test,
                      normalization.method = "LogNormalize",
                      scale.factor = 1e4)
  
hist(colSums(as.matrix(test@data)),
       breaks = 100,
       main = "Total expression after normalisation",
       xlab = "Sum of expression")

# Once the data is normalised, the next step is to find genes which vary between single cells; 
# genes that are constant among all cells have no distinguishing power. 
# The FindVariableGenes() function calculates the average expression and dispersion for each 
# gene, places these genes into bins, and then calculates a z-score for dispersion within 
# each bin. I interpret that as take each gene, get the average expression and variance of 
# the gene across the 2,638 cells, categorise genes into bins (default is 20) based on their 
# expression and variance, and finally normalise the variance in each bin. This was the 
# same approach in Macosko et al. and new methods for detecting genes with variable expression
# patterns will be implemented in Seurat soon (according to the tutorial). The parameters 
# used below are typical settings for UMI data that is normalised to a total of 10,000 
# molecules and will identify around 2,000 variable genes. The tutorial recommends that 
# users should explore the parameters themselves since each dataset is different.

# the variable genes slot is empty before the analysis
test@var.genes
logical(0)

# refer to ?FindVariableGenes
test <- FindVariableGenes(object = test,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)

# vector of variable genes
head(test@var.genes)

# number of variable genes
length(test@var.genes)

# mean and variance of genes are stored pbmc@hvg.info
head(test@hvg.info)

#Seurat constructs linear models to predict gene expression based on user-defined variables to 
#help remove unwanted sources of variation. The idea is that confounding factors, e.g. batch 
# effects and cell cycle stage, affect the observed gene expression patterns and one should 
# adjust for these factors to infer the "correct" gene expression pattern. Buettner et al. 
# demonstrate how correcting for confounding factors improved their downstream analyses.

# The example provided in the tutorial used the number of detected molecules per cell and 
# the percentage mitochondrial RNA to build a linear model. The scaled z-scored residuals, 
# i.e. how much the actual expression differs from the linear model, are stored in the 
# scale.data slot, which are used for dimensionality reduction and clustering.

# slot is empty before running ScaleData()
test@scale.data

# build linear model using nUMI and percent.mito
test <- ScaleData(object = test,
                  vars.to.regress = c("nUMI", "percent.mito"))

class(test@scale.data)
test@scale.data[1:6, 1:6]

# Many genes and single cells are not 
# interesting because they don't vary a lot. One goal of Principal Component Analysis (PCA) 
# is to find the direction/s (usually the first two principal components) in which there is the 
# most variance. The RunPCA() function performs the PCA on genes in the @var.genes slot by 
# default and this can be changed using the pc.genes parameter.

test <- RunPCA(object = test,
               pc.genes = test@var.genes,
               do.print = TRUE,
               pcs.print = 1:5,
               genes.print = 5)

PrintPCAParams(test)

# The PrintPCA() function outputs a set of genes that most strongly define a set of principal 
# components.

PrintPCA(object = test, pcs.print = 1:2, genes.print = 5, use.full = FALSE)

# Furthermore, Seurat has various functions for visualising the cells and genes that define 
# the principal components.

# visualise top genes associated with principal components
VizPCA(object = test, pcs.use = 1:2)

# The PCAPlot() function plots the principal components from a PCA; cells are coloured by their 
# identity class according to test@ident.

PCAPlot(object = test, dim.1 = 1, dim.2 = 2)

# However, the PCA was only performed on the most variable genes, which is a subset of the 
# dataset. The ProjectPCA step scores each gene in the dataset based on their correlation 
# with the calculated components. This is useful because there may be genes that were not 
# detected as variable genes in the variable gene selection step, which are still strongly 
# correlated with cellular heterogeneity.

# the results of the projected PCA can be explored by setting use.full=TRUE in the functions above
test <- ProjectPCA(object = test, do.print = FALSE)

# The PCHeatmap() function produces a heatmap based on the PCA; by default the function uses 
# the first principal component and plots 30 genes across the number of cells specified in 
# cells.use. Setting cells.use to a number plots the "extreme" cells on both ends of the 
# spectrum, which dramatically speeds plotting for large datasets.

PCHeatmap(object = test,
          pc.use = 1,
          cells.use = 500,
          do.balanced = TRUE,
          label.columns = FALSE)
# Can safely ignore warnings

# Plotting all cells

PCHeatmap(object = test,
          pc.use = 1,
          do.balanced = TRUE,
          label.columns = FALSE)

# Plotting on 12 principal components.
PCHeatmap(object = test,
          pc.use = 1:12,
          cells.use = 500,
          do.balanced = TRUE,
          label.columns = FALSE)

# The next steps are to determine how many principal components to use in downstream analyses, 
# which is an important step for Seurat. The tutorial goes through two methods: one uses a 
# statistical test based on a random null model, which is time-consuming for large datasets 
# due to the resampling and may not return a clear cutoff and the other is a commonly used 
# heuristic. The statistical test is carried out by the JackStraw() function, which randomly 
# permutes a subset of data, and calculates projected PCA scores for these "random" genes. 
# Then compares the PCA scores for the "random" genes with the observed PCA scores to 
# determine statistical significance. End result is a p-value for each gene's association 
# with each principal component. This resampling test was inspired by the jackstraw procedure.
# Significant principal components are those with a strong enrichment of low p-value genes.

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
system.time(
  test <- JackStraw(object = test,
                    num.replicate = 100,
                    do.print = FALSE)
)

# The JackStrawPlot() function provides a visualisation tool for comparing the distribution 
# of p-values for each principal component with a uniform distribution (dashed line). 
# "Significant" principal components will show a strong enrichment of genes with low p-values 
# (solid curve above the dashed line).

JackStrawPlot(object = test, PCs = 1:12)
# In this case, it appears that principal components 1 to 10 are significant.

# Another approach for deciding how many principal components to use is to examine the 
# standard deviations of the principle components, which is performed by the PCElbowPlot() 
# function. A cutoff can be drawn where there is a clear elbow in the graph.

PCElbowPlot(object = test)

# It looks like an elbow would fall around principle component 9.

# Seurat includes a graph-based clustering approach, which is quite technical. The 
# approach was heavily inspired by recent work that applied graph-based clustering approaches 
# to scRNA-seq data, namely SNN-Cliq and PhenoGraph. Below is the technical description of the 
# approach:
  
# Briefly, these methods embed cells in a graph structure (e.g. a K-nearest neighbour (KNN) 
# graph) with edges drawn between cells with similar gene expression patterns, and then attempt 
# to partition this graph into highly interconnected "quasi-cliques" or "communities." As in 
# PhenoGraph, we first construct a KNN graph based on the Euclidean distance in PCA space, and 
# refine the edge weights between any two cells based on the shared overlap in their local 
# neighbourhoods (Jaccard distance). To cluster the cells, we apply modularity optimisation 
# techniques SLM, to iteratively group cells together, with the goal of optimising the standard 
# modularity function.

# The FindClusters() function implements the procedure above. The resolution parameter adjusts 
# the granularity of the clustering with higher values leading to more clusters, i.e. higher 
# granularity. According to the authors of Seurat, setting resolution between 0.6 - 1.2 
# typically returns good results for datasets with around 3,000 cells. The clusters are saved 
# in the @ident slot of the Seurat object.

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
test <- FindClusters(object = test, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

# use PrintFindClustersParams() to print summary
# of parameters used to FindClusters()
PrintFindClustersParams(object = test)

# the authors suggest using the same PCs as input to the clustering analysis
# although computing the tSNE based on scaled gene expression
# is also supported using the genes.use argument

test <- RunTSNE(object = test,
                dims.use = all(),
                do.fast = TRUE)

TSNEPlot(object = test, do.label = TRUE)

# The FindMarkers() function identifies positive and negative markers by comparing genes in 
# cells of one cluster against genes in all other cells. The FindAllMarkers() function 
# automates this process for all clusters, but you can also test groups of clusters vs. 
# each other, or against all cells. From the tutorial:
  
# The min.pct argument requires a gene to be detected at a minimum percentage in either of 
# the two groups of cells, and the thresh.test argument requires a gene to be differentially 
# expressed (on average) by some amount between the two groups. You can set both of these to 
# 0, but with a dramatic increase in time - since this will test a large number of genes that 
# are unlikely to be highly discriminatory. As another option to speed up these computations, 
# max.cells.per.ident can be set. This will downsample each identity class to have no more 
# cells than whatever this is set to. While there is generally going to be a loss in power, 
# the speed increases can be significant and the most highly differentially expressed genes 
# will likely still rise to the top.

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = test,
                                ident.1 = 5,
                                min.pct = 0.25)

head(cluster1.markers)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = test,
                                ident.1 = 5,
                                ident.2 = c(2,4),
                                min.pct = 0.25)
head(cluster5.markers)

# find markers for every cluster compared to all remaining cells, report only the positive ones
test.markers <- FindAllMarkers(object = test,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               thresh.use = 0.25)

head(test.markers)
test.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# Seurat has four tests for differential expression (DE) which can be set with the test.use parameter in the FindMarkers() function:
  
# ROC test
# t-test
# LRT test based on zero-inflated data
# LRT test based on tobit-censoring models
# Let's compare the four different DE methods for defining cluster 1.

levels(test@ident)
table(test@ident)

my_bimod <- FindMarkers(object = test,
                        ident.1 = 1,
                        thresh.use = 0.25,
                        test.use = "bimod",
                        only.pos = TRUE)

my_roc <- FindMarkers(object = test,
                      ident.1 = 1,
                      thresh.use = 0.25,
                      test.use = "roc",
                      only.pos = TRUE)

my_t <- FindMarkers(object = test,
                    ident.1 = 1,
                    thresh.use = 0.25,
                    test.use = "t",
                    only.pos = TRUE)

my_tobit <- FindMarkers(object = test,
                        ident.1 = 1,
                        thresh.use = 0.25,
                        test.use = "tobit",
                        only.pos = TRUE)
# identical set of genes
dim(my_bimod)
dim(my_roc)
dim(my_t)
dim(my_tobit)

# the rankings of the genes are quite similar between the methods
my_gene <- row.names(my_bimod)
a <- 1:length(my_gene)
b <- match(my_gene, row.names(my_roc))
c <- match(my_gene, row.names(my_t))
d <- match(my_gene, row.names(my_tobit))

# bimod vs. bimod
cor(a, a, method = "spearman")
# bimod vs. roc
cor(a, b, method = "spearman")
# bimod vs. t
cor(a, c, method = "spearman")
# bimod vs. tobit
cor(a, d, method = "spearman")

par(mfrow=c(2,2))
barplot(a, main = 'bimod')
barplot(b, main = 'roc')
barplot(c, main = 't')
barplot(d, main = 'tobit')

# The VlnPlot() and FeaturePlot() functions can be used to visualise marker expression.

VlnPlot(object = test, features.plot = c("Ms4a1", "Cd79a"))

# visualise markers found by the FindMarkers() analysis above
head(my_tobit)

VlnPlot(object = test, features.plot = c("S100a9", "S100a8"))

# you can plot raw UMI counts as well
VlnPlot(object = test,
        features.plot = c("Nkg7", "Pf4"),
        use.raw = TRUE,
        y.log = TRUE)

FeaturePlot(object = test,
            features.plot = c("Nkg7", "Cd3e", "Cd14","Ccr5", "Cd68", "Csf1r", "Cd200", "Foxp3", "Chil3"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

FeaturePlot(object = test,
            features.plot = head(row.names(my_tobit), 9),
            cols.use = c("grey", "blue"))

# The DoHeatmap() function creates a heatmap of genes across all cells. Below, we use plot the top 10 marker genes 
# for the eight clusters.

head(test.markers)
test.markers %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC) -> top10

head(top10)

DoHeatmap(object = test,
          genes.use = top10$gene,
          slim.col.label = TRUE,
          remove.key = TRUE)

# There are canonical markers for each cell type, which can be used to assess the clustering.

FeaturePlot(object = test,
            features.plot = c("Il7r", "Cd14", "Ms4a1", "Cd8a","Ms4a7", "Nkg7", "Cst3", "Ppbp"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

# Relabel Cluster IDs 
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6)
new.cluster.ids <- c("NK Cells",
                     "T-Cells",
                     "Mature Mac",
                     "FoxP3 Tregs",
                     "Monocytes",
                     "DCs",
                     "DC-2")

test@ident <- plyr::mapvalues(x = test@ident,
                              from = current.cluster.ids,
                              to = new.cluster.ids)


TSNEPlot(object = test, do.label = TRUE, pt.size = 0.5)

# Seurat provides the StashIdent() function for keeping cluster IDs; this is useful for testing various parameters and 
# comparing the clusters. For example, adjusting the parameters may lead to the CD4 T cells subdividing into two groups.

# stash cluster identities for later
#test <- StashIdent(object = test, save.name = "ClusterNames_0.6")

#test <- FindClusters(object = test,
#                     reduction.type = "pca",
#                     dims.use = 1:10,
#                     resolution = 0.8,
#                     print.output = FALSE)

# plot two tSNE plots side by side, and colour points based on different criteria
plot1 <- TSNEPlot(object = test,
                  do.return = TRUE,
                  no.legend = TRUE,
                  do.label = TRUE)

#plot2 <- TSNEPlot(object = test,
#                  do.return = TRUE,
#                  group.by = "ClusterNames_0.6",
#                  no.legend = TRUE,
#                  do.label = TRUE)

#plot_grid(plot1, plot2)

# Analyze TLR 7 & 8 levels as well as coreceptors in groups
cluster.averages <- AverageExpression(object = test)
#head(x = cluster.averages[, 1:7])

# Return this information as a Seurat object (enables downstream plotting
# and analysis)
cluster.averages <- AverageExpression(object = test, return.seurat = FALSE, show.progress = FALSE)

# How can I plot the average expression of NK cells vs. T cells?  Pass
# do.hover = T for an interactive plot to identify gene outliers
CellPlot(object = cluster.averages, cell1 = "NK Cells", cell2 = "T-Cells")

# Return the averages for individual genes across the groups
# TLR data
Tlr1 = cluster.averages["Tlr1", 1:7]
Tlr2 = cluster.averages["Tlr2", 1:7]
Tlr3 = cluster.averages["Tlr3", 1:7]
Tlr4 = cluster.averages["Tlr4", 1:7]
Tlr5 = cluster.averages["Tlr5", 1:7]
Tlr6 = cluster.averages["Tlr6", 1:7]
Tlr7 = cluster.averages["Tlr7", 1:7]
Tlr8 = cluster.averages["Tlr8", 1:7]
Tlr11 = cluster.averages["Tlr11", 1:7]
Tlr = rbind(Tlr1, Tlr2, Tlr3, Tlr4, Tlr5, Tlr6, Tlr7, Tlr8, Tlr11)

# Nfkb 
Nfkb1 = cluster.averages["Nfkb1", 1:7]
Nfkb2 = cluster.averages["Nfkb2", 1:7]
Nfkbia = cluster.averages["Nfkbia", 1:7]
Nfkbib = cluster.averages["Nfkbib", 1:7]
Nfkbid = cluster.averages["Nfkbid", 1:7]
Nfkbie = cluster.averages["Nfkbie", 1:7]
Nfkbil1 = cluster.averages["Nfkbil1", 1:7]
Nfkbiz = cluster.averages["Nfkbiz", 1:7]
Nfkb = rbind(Nfkb1, Nfkb2, Nfkbia, Nfkbib, Nfkbid, Nfkbie, Nfkbil1, Nfkbiz)

# Other pathway components
Myd88 = cluster.averages["Myd88", 1:7]
Irf3 = cluster.averages["Irf3", 1:7]
Irf7 = cluster.averages["Irf7", 1:7]
Traf3 = cluster.averages["Traf3", 1:7]
Irak1 = cluster.averages["Irak1", 1:7]
Irak4 = cluster.averages["Irak4", 1:7]
Traf6 = cluster.averages["Traf6", 1:7]
Tollip = cluster.averages["Tollip", 1:7]
Ikk = cluster.averages["Ikk", 1:7]
Tbk1 = cluster.averages["Tbk1", 1:7]
Tlr_pathway = rbind(Myd88, Irf3, Irf7, Traf3, Irak1, Irak4, Traf6, Tollip, Ikk, Tbk1) 

# Plot data
barplot(t(as.matrix(Tlr)), beside=TRUE, ylab = "Normalized Gene Expression", col = brewer.pal(7, "Set1"), 
        legend = new.cluster.ids)
barplot(t(as.matrix(Nfkb)), beside=TRUE, ylab = "Normalized Gene Expression", col = brewer.pal(7, "Set1"), 
        legend = new.cluster.ids)
barplot(t(as.matrix(Tlr_pathway)), beside=TRUE, ylab = "Normalized Gene Expression", col = brewer.pal(7, "Set1"), 
        legend = new.cluster.ids)