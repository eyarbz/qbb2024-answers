# ====================
# Exercise 3
# ====================

# Step 3.1: Loading and filtering the data

BiocManager::install("vsn")
install.packages("ggfortify")
install.packages("hexbin")

library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(ggfortify)

counts <- readr::read_tsv("/Users/cmdb/qbb2024-answers/week5/salmon.merged.gene_counts.tsv")
counts <- column_to_rownames(counts, var = "gene_name")
counts <- dplyr::select(counts, -gene_id)
counts <- dplyr::mutate_if(counts, is.numeric, as.integer)
counts <- counts[rowSums(counts) > 100,]
counts <- dplyr::select(counts, "A1_Rep1": "P2-4_Rep3")

# Step 3.2: Creating DESeq2 model and batch-correction

metadata = tibble(tissue=as.factor(c("A1", "A1", "A1",
                                     "A2-3", "A2-3", "A2-3",
                                     "Cu", "Cu", "Cu", 
                                     "LFC-Fe", "LFC-Fe", "Fe", 
                                     "LFC-Fe", "Fe", "Fe", 
                                     "P1", "P1", "P1", 
                                     "P2-4", "P2-4", "P2-4")), 
                  rep=as.factor(c(1, 2, 3,
                                  1, 2, 3,
                                  1, 2, 3, 
                                  1, 2, 3,
                                  1, 2, 3,
                                  1, 2, 3,
                                  1, 2, 3)))
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), 
                              colData = metadata, 
                              design = ~tissue)
vsd <- vst(dds)
meanSdPlot(assay(vsd))

# Step 3.3: PCA analysis

pca_data = plotPCA(vsd, intgroup=c("rep","tissue"), returnData=TRUE)
PCA_plot <- ggplot(pca_data, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)
ggsave("/Users/cmdb/qbb2024-answers/week5/PCA_plot.png", plot = PCA_plot, dpi = 600)

# Step 3.4: Filtering genes by variance

vsd_matrix <- as.matrix(assay(vsd))
combined = vsd_matrix[,seq(1, 21, 3)]
combined = combined + vsd_matrix[,seq(2, 21, 3)]
combined = combined + vsd_matrix[,seq(3, 21, 3)]
combined = combined / 3
sds = rowSds(combined)

vsd_matrix_filtered <- vsd_matrix[sds>1,]

# Step 3.5: K-means clustering genes

set.seed(42)
kmeans_result <- kmeans(vsd_matrix_filtered, centers = 12)
kmeans_labels <- kmeans_result$cluster

ordering <- order(kmeans_labels)
kmeans_labels_ordered <- kmeans_labels[ordering]
vsd_matrix_filtered_ordered <- vsd_matrix_filtered[ordering,]

heatmap(vsd_matrix_filtered_ordered, Rowv=NA, Colv=NA, RowSideColors=RColorBrewer::brewer.pal(12,"Paired")[kmeans_labels_ordered])

# Step 3.6: Gene ontology enrichment analysis

cluster1_gene_names <- rownames(vsd_matrix_filtered[kmeans_labels == 1,])
write(cluster1_gene_names, file = "/Users/cmdb/qbb2024-answers/week5/cluster1_gene_names.txt")

