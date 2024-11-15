# ==========================
# Exercise 1: Perform principal component analysis
# ==========================

#Step 1.1: Loading data and importing libraries
library(DESeq2)
library(tidyverse)
library(broom)

gene_locations <- read_delim("/Users/cmdb/qbb2024-answers/week7/gene_locations.txt")
gtex_metadata <- read_delim("/Users/cmdb/qbb2024-answers/week7/gtex_metadata_downsample.txt")
blood_counts <- read_delim("/Users/cmdb/qbb2024-answers/week7/gtex_whole_blood_counts_downsample.txt")

blood_counts <- column_to_rownames(blood_counts, var = "GENE_NAME")
gtex_metadata <- column_to_rownames(gtex_metadata, var = "SUBJECT_ID")

#Step 1.2: Create a DESeq2 object
if (all(colnames(blood_counts) == rownames(gtex_metadata)) == TRUE){
  print("Columns match rows!")
}

dds <- DESeqDataSetFromMatrix(countData = blood_counts, 
                              colData = gtex_metadata, 
                              design = ~ SEX + DTHHRDY + AGE)

#Step 1.3: Normalization and PCA
vsd <- vst(dds)

plotPCA(vsd, intgroup = "SEX") +
  scale_color_manual(values = c("limegreen", "hotpink"),
                     labels = c("Male", "Female")) +
  labs(title = "Whole blood gene expression PCA", color = "Biological sex")

plotPCA(vsd, intgroup = "AGE") +
  scale_color_continuous(labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")) +
  labs(title = "Whole blood gene expression PCA", color = "Age group")

plotPCA(vsd, intgroup = "DTHHRDY") +
  scale_color_manual(values = c("turquoise", "salmon"),
                     labels = c("Fast death of natural causes", "Ventilator death")) +
  labs(title = "Whole blood gene expression PCA", color = "Death classification")

#The first PC captures 48% of the variance and the second PC captures 7%. It looks like fast death of natural causes and ventilator case are associated with the principal components since they group on two distinct sides.

# ==========================
# Exercise 2: Perform differential expression analysis
# ==========================

#Step 2.1: 
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()

vsd_df <- bind_cols(gtex_metadata, vsd_df)

m1 <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

#WASH7P does not show significant evidence of sex-differential expression because the p-value associated with this predictor is very high (~0.2).

m2 <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

#SLC25A47 does  show significant evidence of sex-differential expression because the p-value associated with this predictor is less than 0.05 (~0.02). Males have a higher expression of SLC25A47 than females.

#Step 2.2:
dds <- DESeq(dds)

#Step 2.3:
res <- results(dds, name = "SEX_male_vs_female")  %>%
  as_tibble(rownames = "GENE_NAME")

sum(res$padj < 0.1, na.rm = TRUE)

#262 genes exhibit significant differential expression between males and females at a 10% FDR

res_locs <- left_join(res, gene_locations, by = "GENE_NAME")

#Chromosomes Y and X encode the genes that are most strongly upregulated in males versus females. There are more male upregulated genes near the top of the list. This is because female is set as the reference level, so genes that are more upregulated in males will appear.

WASH7P_padj <- res_locs[which(res_locs$GENE_NAME == "WASH7P"), "padj"]
#WASH7P has a padj value of 0.8986810, greater than 0.1, meaning it is not differentially expressed. This matches result from earlier where p value was not low enough to be significant.

SLC25A47_padj <- res_locs[which(res_locs$GENE_NAME == "SLC25A47"), "padj"]
#SLC25A47 has a padj value of 8.322385e-07, which is much less than 0.1, meaning it is differentially expressed. This matches result from earlier where p value was less than 0.5 and considered significant.

#Step 2.4:
res2 <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes")  %>%
  as_tibble(rownames = "GENE_NAME")
sum(res2$padj < 0.1, na.rm = TRUE)

#16069 genes exhibit significant differential expression according to death classification at a 10% FDR
#The PCA showed that PC1 captured 48% of the variance, and death classification separated well on that, so it makes sense that death classification is responsible for genes being more differentially expressed

# ==========================
# Exercise 3: Visualization
# ==========================

res_v_plot <- res %>% dplyr::filter(!is.na(padj)) %>%
  mutate(trans_padj = -log10(padj))

v_plot_sex_DE <- ggplot(res_v_plot, aes(x = log2FoldChange, y = trans_padj)) +
  geom_point(aes(color = padj < 0.1 & abs(log2FoldChange) > 1), size = 3) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Sex differential expression", 
       x = "Fold change (log2)",
       y = "P-adjusted (-log10)",
       color = "P-adjusted < 0.1 and log2 fold change > 1")