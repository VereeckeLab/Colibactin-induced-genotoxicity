################################################################################
########### Load Packages
################################################################################
library(readxl)
library(writexl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gtools)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
library(ggpubr)
library(ggforce)
library(ggrepel)
library(DESeq2)
library(edgeR)
library(apeglm)
library(ashr)
library(org.Mm.eg.db)
library(fgsea)

source("2.Script_Functions.R")

################################################################################
########### Load DESeq2 Object and needed arguments
################################################################################
# Load dds object created in 1.Script_DESeq2_Pre_processing.R
filtered_dds <- readRDS("Data/DESeq2_object.rds")

# Load annotation file
anno_file <- read_csv("Data/GRCm39_full_annotation.txt")

# Set cutoffs
prefilter_cutoff <- 1
sign_cutoff <- 0.05
fc_cutoff <- 0

# define factors
factor_levels <- c("condition","group","genotype","tissue","type","eartag")



################################################################################
########### DE analysis 11G5 vs Nissle 1917 Immune (group B vs F)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("B","F")] 
# Drop unused levels to avoid errors
for (factor in factor_levels){
  subset_dds[[factor]] <- droplevels(subset_dds[[factor]])
}
# Filter genes on counts to avoid low baseMean genes in analysis
num_samples <- min(summary(subset_dds[["group"]]))
keep <- rowSums(counts(subset_dds) >= prefilter_cutoff) >= num_samples
subset_dds <- subset_dds[keep,]
# quick check
head(subset_dds@colData)
vsd <- vst(subset_dds, blind = F)
# QC
plotPCA(vsd, intgroup = "condition") + theme_classic()

# Find differentially expressed genes
title = "11G5 TG vs Nissle 1917 TG (Immune)"
# Set reference group
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "F")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","B","F"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "BvsF_11G5_CD45_TGvsNISSLE_CD45_TG")

################################################################################
########### PCA (Supplementary Figure 1C)
################################################################################
colData(subset_dds)$Condition <- colData(subset_dds)$type
# Change level names
levels(subset_dds$Condition) <- c("11G5","Nissle 1917")
vsd <- vst(subset_dds, blind = F)
colData(subset_dds)

# pca data
pca.data <- plotPCA(vsd,"Condition",returnData = T)

# Get variance percentage for axis labels
percentVar <- round(100 * attr(pca.data, "percentVar"))

# Figure 2A
pl1 <- ggplot(pca.data, aes(PC1, PC2, color=Condition)) +
  geom_point(size=5) +
  geom_mark_ellipse(n = 1000,expand = unit(20,"mm"),show.legend = FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlim(c(-20,20)) +
  ylim(c(-25,25)) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(title) +
  scale_color_manual(values = c("Nissle 1917" = "#6B92CB", "11G5" = "#ED6E68")) +
  coord_fixed() + 
  theme_classic() +
  theme(text = element_text(family="Calibri"),
        axis.title = element_text(face = "bold", color = "black",size = 22,hjust = 0.5),
        axis.text = element_text(color = "black", size = 16),
        plot.title = element_text(face = "bold", color = "black",size = 24,hjust = 0.5),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 22, face = "bold"))
print(pl1)
# Save figure
ggsave("Data/Figures/B_vs_F_PCA.png",pl1, units = "px", width = 1000, height = 1000, dpi = 100, bg = "white")


################################################################################
########### VolcanoPlot (Supplementary Figure 1D)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/BvsF_11G5_CD45_TGvsNISSLE_CD45_TG_filtered.xlsx")
# Make figure
DEG_names <- deframe(results[results$DE_gene == TRUE,c("Gene.name")])
pl2 <- get_paper_volcano_plot(result_table = results,
                              important_genes = DEG_names,
                              crop_value_cutoff = 1e-10,
                              padj_cutoff = 0.05,
                              logFC_cutoff = 0,
                              gene_num = 10000,
                              point_size = 4,
                              colors = c("#6B92CB","lightgrey","#ED6E68"),
                              title = "",
                              figure_path = "",
                              figure_name = paste0(title,".png"))

ggarrange(pl1,pl2 + theme(legend.background = element_rect(size=0.5, linetype="solid",colour ="black")))
ggsave("Data/Figures/Immune_11G5_TG_vs_Nissle_TG.png",
      units = "px",
      width = 1920, height = 1043,
      dpi = 100,
      bg = "white")


# Save Figure
ggsave("Data/Figures/B_vs_F_Volcano.png",pl2, units = "px", width = 1644/3.25, height = 1291/3.1, dpi = 100, bg = "white")





################################################################################
########### DE analysis 11G5 vs Nissle 1917 Epithelial (group C vs G)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("C","G")] 
# Drop unused levels to avoid errors
for (factor in factor_levels){
  subset_dds[[factor]] <- droplevels(subset_dds[[factor]])
}
# Filter genes on counts to avoid low baseMean genes in analysis
num_samples <- min(summary(subset_dds[["group"]]))
keep <- rowSums(counts(subset_dds) >= prefilter_cutoff) >= num_samples
subset_dds <- subset_dds[keep,]
# quick check
head(subset_dds@colData)
vsd <- vst(subset_dds, blind = F)
# QC
plotPCA(vsd, intgroup = "condition") + theme_classic()

# Find differentially expressed genes
title = "11G5 WT vs Nissle 1917 WT (Epithelial)"
# Set reference group
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "G")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","C","G"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "CvsG_11G5_EPI_WTvsNISSLE_EPI_WT")

################################################################################
########### PCA (Supplementary Figure 1C)
################################################################################
colData(subset_dds)$Condition <- colData(subset_dds)$type
# Change level names
levels(subset_dds$Condition) <- c("11G5","Nissle 1917")
vsd <- vst(subset_dds, blind = F)
colData(subset_dds)

# pca data
pca.data <- plotPCA(vsd,"Condition",returnData = T)

# Get variance percentage for axis labels
percentVar <- round(100 * attr(pca.data, "percentVar"))

# Figure 2A
pl1 <- ggplot(pca.data, aes(PC1, PC2, color=Condition)) +
  geom_mark_ellipse(n = 1000,expand = unit(20,"mm"),show.legend = FALSE) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlim(c(-20,20)) +
  ylim(c(-25,25)) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(title) +
  scale_color_manual(values = c("Nissle 1917" = "#6B92CB", "11G5" = "#ED6E68")) +
  coord_fixed() + 
  theme_classic() +
  theme(text = element_text(family="Calibri"),
        axis.title = element_text(face = "bold", color = "black",size = 22,hjust = 0.5),
        axis.text = element_text(color = "black", size = 16),
        plot.title = element_text(face = "bold", color = "black",size = 24,hjust = 0.5),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 22, face = "bold"))
print(pl1)

# Save figure
ggsave("Data/Figures/C_vs_G_PCA.png",pl1, units = "px", width = 1000, height = 1000, dpi = 100, bg = "white")


################################################################################
########### VolcanoPlot (Supplementary Figure 1D)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/CvsG_11G5_EPI_WTvsNISSLE_EPI_WT_filtered.xlsx")
# Make figure
DEG_names <- deframe(results[results$DE_gene == TRUE,c("Gene.name")])
pl2 <- get_paper_volcano_plot(result_table = results,
                              important_genes = DEG_names,
                              crop_value_cutoff = 1e-10,
                              padj_cutoff = 0.05,
                              logFC_cutoff = 0,
                              gene_num = 10000,
                              point_size = 4,
                              colors = c("#6B92CB","lightgrey","#ED6E68"),
                              title = title,
                              figure_path = "Data/Results/",
                              figure_name = paste0(title,".png"))

ggarrange(pl1,pl2+ theme(legend.background = element_rect(size=0.5, linetype="solid",colour ="black")))

# Save Figure
ggsave("Data/Figures/Epithelial_11G5_WT_vs_NISSLE_WT.png",
      units = "px",
      width = 1920, height = 1043,
      dpi = 100,
      bg='white')


################################################################################
########### DE analysis 11G5 vs Nissle 1917 Immune (group D vs H)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("D","H")] 
# Drop unused levels to avoid errors
for (factor in factor_levels){
  subset_dds[[factor]] <- droplevels(subset_dds[[factor]])
}
# Filter genes on counts to avoid low baseMean genes in analysis
num_samples <- min(summary(subset_dds[["group"]]))
keep <- rowSums(counts(subset_dds) >= prefilter_cutoff) >= num_samples
subset_dds <- subset_dds[keep,]
# quick check
head(subset_dds@colData)
vsd <- vst(subset_dds, blind = F)
# QC
plotPCA(vsd, intgroup = "condition") + theme_classic()

# Find differentially expressed genes
title = "11G5 WT vs Nissle 1917 WT (Immune)"
# Set reference group
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "H")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","D","H"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "DvsH_11G5_CD45_WTvsNISSLE_CD45_WT")

################################################################################
########### VolcanoPlot (Supplementary Figure 1K)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/DvsH_11G5_CD45_WTvsNISSLE_CD45_WT_filtered.xlsx")
# Make figure
pl2 <- get_paper_volcano_plot(result_table = results,
                              important_genes = c("Mt1","Mt2"),
                              crop_value_cutoff = 1e-10,
                              padj_cutoff = 0.05,
                              logFC_cutoff = 0,
                              gene_num = 10000,
                              point_size = 4,
                              colors = c("#6B92CB","lightgrey","#ED6E68"),
                              title = title,
                              figure_path = "Data/Results",
                              figure_name = paste0(title,".png"))
pl2 = pl2 + theme(legend.background = element_rect(size=0.5, linetype="solid",colour ="black"))
# Save Figure
ggsave("Data/Figures/D_vs_H_Volcano.png",pl2, units = "px", width = 1000, height = 1200, dpi = 100, bg = "white")


