# ===============
# Exercise 1: Load data
# ===============
#Load packages
library(zellkonverter)
library(scuttle)
library(scater)
library(scran)
library(ggplot2)

#Load and inspect data (0.5 pt)
gut <- readH5AD("/Users/cmdb/qbb2024-answers/week8/v2_fca_biohub_gut_10x_raw.h5ad")
assayNames(gut) <- "counts"
gut <- logNormCounts(gut)

#Q1: How many genes are quantitated (should be >10,000)?
nrow(gut)
# 13,407 genes

#Q1: How many cells are in the dataset?
ncol(assay(gut))
# 11,788 cells

#Q1: What dimension reduction datasets are present?
reducedDimNames(gut)
# pca, tsne, umap

#Inspect cell metadata (0.5 pt)

#Q2: How many columns are there in colData(gut)?
cell_metadata <- as.data.frame(colData(gut))
#39 columns

#Q2: Which three column names reported by colnames() seem most interesting? Briefly explain why.
#I think percent_mito is the most interesting because I did not think of mitochondrial RNA being a potential contaminant, so it was interesting to learn this.

#Q2: Plot cells according to X_umap using plotReducedDim() and colouring by broad_annotation
plotReducedDim(gut, dim = "X_umap", color_by = "broad_annotation")

# ===============
# Exercise 2: Explore data
# ===============
#Explore gene-level statistics (1 pt)
genecounts <- rowSums(assay(gut))

#Q3: What is the mean and median genecount according to summary()? What might you conclude from these numbers?
genecounts_sum <- summary(genecounts)
#Mean is 3,185 and median is 254. The gene counts vary greatly which is why the median is so different from the mean.

#Q3: What are the three genes with the highest expression after using sort()? What do they share in common?
genecounts_sort <- sort(genecounts, decreasing = TRUE)
#Highest genes are lncRNA:Hsromega, pre-rRNA:CR45845 and lncRNA:roX1. These are all non coding RNAs rather than protein coding mRNAs

#Explore cell-level statistics (1 pt)
#Q4a: Explore the total expression in each cell across all genes (0.5 pt)
cellcounts <- colSums(assay(gut))
cellcounts_hist <- hist(cellcounts)

#Q4a: What is the mean number of counts per cell?
cellcounts_sum <- summary(cellcounts)
#Mean is 3,622 counts per cell. 

#Q4a: How would you interpret the cells with much higher total counts (>10,000)?
#Maybe contamination from other sources (e.g. rRNA)

#Q4b: Explore the number of genes detected in each cell (0.5 pt)
celldetected <- colSums(assay(gut)>0)
celldetected_hist <- hist(celldetected)

#Q4b: What is the mean number of genes detected per cell?
celldetected_sum <- summary(celldetected)
#Mean is 1,059 genes detected per cell

#Q4b: What fraction of the total number of genes does this represent?
#This represents 1,059/13,407 ~ 0.08

#Explore mitochondrial reads (1 pt)
mito <- grep("^mt:", rownames(gut), value = TRUE)
df <- perCellQCMetrics(gut, subsets=list(Mito=mito))
df_sum <- summary(as.data.frame(df))
#sum mean = 3,622 and detected mean = 1,059, which matches previous calculations
colData(gut) <- cbind( colData(gut), df )

#Q5: Visualize percent of reads from mitochondria (1 pt)
plotColData(gut, x = "broad_annotation", y = "subsets_Mito_percent") +
  theme( axis.text.x=element_text( angle=90 ) )

#Q5: Which cell types may have a higher percentage of mitochondrial reads? Why might this be the case?
#Muscle cells have a high percentage of mitochondrial reads. This makes sense because these cells require a lot of energy and will have more mitochondria.

# ===============
# Exercise 3: Identify marker genes
# ===============
#Analyze epithelial cells (3 pt)
#Q6a: Subset cells annotated as “epithelial cell” (1 pt)
coi <- colData(gut)$broad_annotation == "epithelial cell"
epi <- gut[,coi]
plotReducedDim(epi, dim = "X_umap", color_by = "annotation")

marker.info <- scoreMarkers( epi, colData(epi)$annotation )
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

#Q6b: Evaluate top marker genes (2 pt)
#What are the six top marker genes in the anterior midgut? Based on their functions at flybase.org, what macromolecule does this region of the gut appear to specialize in metabolizing?
#Mal-A6, Men-b, vnd, betaTry, Mal-A1, and Nhe2. Based on their functions, the anterior midgut probably specializes in metabolising carbohydrates.

plotExpression(epi, features = c("Mal-A6", "Men-b", "vnd", "betaTry", "Mal-A1", "Nhe2"), x = "annotation") +
  theme( axis.text.x=element_text( angle=90 ) )

#Analyze somatic precursor cells (3 pt)
coi <- colData(gut)$broad_annotation == "somatic precursor cell"
som <- gut[,coi]
marker.info <- scoreMarkers( som, colData(som)$annotation )
chosen <- marker.info[["intestinal stem cell"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
goi <- rownames(ordered)[1:6]
plotExpression(som, features = goi, x = "annotation") +
  theme( axis.text.x=element_text( angle=90 ) )

#Q7: Which two cell types have more similar expression based on these markers?
#Enteroblast and intestinal stem cell

#Q7: Which marker looks most specific for intestinal stem cells?
#DI