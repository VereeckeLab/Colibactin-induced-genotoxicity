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
library(DESeq2)
library(edgeR)
library(apeglm)
library(ashr)
library(org.Mm.eg.db)
library(fgsea)
library(ggrepel)
library(genekitr)
library(gplots)
library(GSVA)
library(stringr)


source("2.Script_Functions.R")

################################################################################
########### Load DESeq2 Object and needed arguments
################################################################################
# Load dds object created in 1.Script_DESeq2_Pre_processing.R
filtered_dds <- readRDS("Data/DESeq2_object.rds")

# Load annotation file
anno_file <- read_csv("Data/GRCm39_full_annotation.txt")

# Set cutoffs
prefilter_cutoff <- 1 # Impoves gsea 10 original was 1
sign_cutoff <- 0.05
fc_cutoff <- 0

# define factors
factor_levels <- c("condition","group","genotype","tissue","type","eartag")


################################################################################
########### DE analysis 11G5 vs Nissle 1917 Epithelial (group A vs E)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("A","E")] 
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
title = "11G5 TG vs Nissle 1917 TG (Epithelial)"
# Set reference group
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "E")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
resultsNames(f_dds_DE)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","A","E"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "AvsE_11G5_EPI_TGvsNISSLE_EPI_TG")


################################################################################
########### PCA (Figure 2A)
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
  geom_point(size=6) +
  geom_mark_ellipse(n = 1000,expand = unit(10,"mm"),show.legend = FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlim(c(-15,15)) +
  ylim(c(-15,15)) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA") +
  scale_color_manual(values = c("Nissle 1917" = "#6B92CB", "11G5" = "#ED6E68")) +
  coord_fixed() + 
  theme_classic() +
  theme(text = element_text(family="Calibri"),
        axis.title = element_text(face = "bold", color = "black",size = 23,hjust = 0.5),
        axis.text = element_text(color = "black", size = 17),
        plot.title = element_text(face = "bold", color = "black",size = 25,hjust = 0.5),
        legend.text = element_text(size = 21),
        legend.title = element_text(size = 27, face = "bold"))
print(pl1)
# Save figure
ggsave("Data/Figures/A_vs_E_PCA.png",pl1, units = "px", width = 1000, height = 1000, dpi = 100, bg = "white")


################################################################################
########### VolcanoPlot (Figure 2B)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/AvsE_11G5_EPI_TGvsNISSLE_EPI_TG_filtered.xlsx")
# Make figure
pl2 <- get_paper_volcano_plot(result_table = results,
                              important_genes = c("Adgrg2","Defb37","Ifi44","Trim3","Snta1","Slc37a2"),
                              crop_value_cutoff = 1e-15,
                              padj_cutoff = 0.05,
                              logFC_cutoff = 0,
                              gene_num = 10000,
                              point_size = 2,
                              colors = c("#6B92CB","lightgrey","#ED6E68"),
                              title = "",
                              figure_path = "Data/Results/",
                              figure_name = paste0(title,".png"))

# Extra tweaks
pl2 <- pl2 + theme(text = element_text(family="Calibri"),
          axis.title = element_text(face = "bold", color = "black",size = 16,hjust = 0.5),
          axis.text = element_text(color = "black", size = 12),
          plot.title = element_text(face = "bold", color = "black",size = 18,hjust = 0.5),
          legend.text = element_text(size = 13),
          legend.position = c(0.84, 0.90),
          legend.title = element_text(size = 13, face = "bold"))

# Save Figure
ggsave("Data/Figures/A_vs_E_Volcano.png",pl2, units = "px", width = 1644/3.25, height = 1291/3.1, dpi = 100, bg = "white")

################################################################################
########### GSEA barplot (Figure 2C)
################################################################################
# DISCLAIMER: Keep in mind scores can change slightly when rerunning the analysis!
# (https://github.com/ctlab/fgsea/issues/12)
# Fetch shrunk results using ashr
shrunk_DE_results <- lfcShrink(dds = f_dds_DE,
                            contrast=c("group","A","E"),
                            type="ashr",
                            res=DE_result)
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = shrunk_DE_results,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "AvsE_11G5_EPI_TGvsNISSLE_EPI_TG_GSEA")



# Load shrunk results
shrunk_results <- read_xlsx("Data/Results/AvsE_11G5_EPI_TGvsNISSLE_EPI_TG_GSEA_filtered.xlsx")
# Generate ranked gene file
#ranks_lfc <- shrunk_results[order(shrunk_results$log2FoldChange,decreasing = T) & (shrunk_results$baseMean > 30) ,c("Entrez_GSEA","log2FoldChange")] # clean up useless results
ranks_lfc <- shrunk_results[order(shrunk_results$log2FoldChange,decreasing = T) ,c("Entrez_GSEA","log2FoldChange")]
# If duplicate gene names present, average the values
if( sum(duplicated(ranks_lfc$Entrez_GSEA)) > 0) {
  ranks_lfc <- aggregate(.~Entrez_GSEA, FUN = mean, data = ranks_lfc)
  ranks_lfc <- ranks_lfc[order(ranks_lfc$log2FoldChange, decreasing = T),]
}
# Turn the dataframe into a named vector for fgsea()
ranks_lfc_deframed <- deframe(ranks_lfc)

# load msigdb data
Mm.c2.all.v7.1.entrez <- readRDS("Data/Mm.c2.all.v7.1.entrez.rds")
Mm.c5.bp.v7.1.entrez <- readRDS("Data/Mm.c5.bp.v7.1.entrez.rds")
length(Mm.c2.all.v7.1.entrez)
length(Mm.c5.bp.v7.1.entrez)

# Perform GSEA
fgsea_results <- fgsea(c(Mm.c5.bp.v7.1.entrez,Mm.c2.all.v7.1.entrez),
                       ranks_lfc_deframed, nPermSimple = 1000) #nPermSimple = 1000 (DEFAULT)


# Add column showing significance based on msigdb FAQ (padj < 0.25)
fgsea_results$DE = FALSE
fgsea_results[(fgsea_results$padj < 0.25) & (fgsea_results$NES > 1) ,"DE"] = TRUE
# Convert and add leading edge as gene symbols
for (row_idx in 1:nrow(fgsea_results)){
    # fetch row
    row = fgsea_results[row_idx,] 
    # unlist leadingEdge Convert entrezID to symbols (mouse)
    symbols = na.omit(select(org.Mm.eg.db, keys = unlist(row$leadingEdge), keytype = "ENTREZID", column = "SYMBOL"))
    # Set to results
    fgsea_results[row_idx,"leadingEdge"] = paste0(symbols$SYMBOL, collapse = ',')
}
fgsea_results$leadingEdge = as.character(fgsea_results$leadingEdge)

# Save results
#write_xlsx(fgsea_results,"Data/Results/A_E_GSEA.xlsx") # reload rsults for consistency figures
fgsea_results = read_xlsx("/home/maartenc/Documents/GEO/Maude/Extra/important/A_E_GSEA_DE_MJ.xlsx")

# Select Significant pathways to show
main_selection <- c("GO_negative_regulation_of_double-strand_break_repair",
                    "GO_negative_regulation_of_DNA_repair",
                    "PID_P38_GAMMA_DELTA_PATHWAY",
                    "GO_regulation_of_double-strand_break_repair_via_homologous_recombination",
                    "GO_DNA_damage_response,_signal_transduction_resulting_in_transcription",
                    "REACTOME_GLUCONEOGENESIS",
                    "PID_E2F_PATHWAY",
                    "GO_positive_regulation_of_double-strand_break_repair_via_nonhomologous_end_joining",
                    "KANNAN_TP53_TARGETS_UP",
                    "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE",
                    "FARDIN_HYPOXIA_9",
                    "KOBAYASHI_EGFR_SIGNALING_24HR_UP",
                    "REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX",
                    "REACTOME_O_LINKED_GLYCOSYLATION",
                    "REACTOME_SIGNALING_BY_NOTCH",
                    "GO_regulation_of_cell_growth_by_extracellular_stimulus", 
                    "GO_actin_cytoskeleton_reorganization",
                    "GO_lipopolysaccharide-mediated_signaling_pathway",
                    "PID_IL23_PATHWAY",
                    "GO_response_to_bacterium",
                    "GO_regulation_of_phagocytosis",
                    "GO_interleukin-10_production",
                    "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                    "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
                    "GO_interleukin-17_production",
                    "GO_interleukin-4_production",
                    "RUAN_RESPONSE_TO_TNF_UP",
                    "GO_positive_regulation_of_cell-cell_adhesion",
                    "GO_negative_regulation_of_toll-like_receptor_4_signaling_pathway",
                    "BIERIE_INFLAMMATORY_RESPONSE_TGFB1"
                    )

fgsea_results <- fgsea_results[fgsea_results$pathway %in% main_selection,]
view(fgsea_results)
# function for plotting
g1 <- GSEA_plot(fgsea_results,
          title = "",
          NES_cutoff = 0,
          npathw = 100,
          color_down = "#6B92CB",
          color_up = "#ED6E68",
          reverse_order = FALSE,
          y_text_size = 10) #12
g1 <- g1 + theme(legend.position = "None")
g1
# Save Figure
ggsave("Data/Figures/Epithelial_11G5_TG_vs_Nissle_TG_GSEA.png", units = "px", width = (915) * 1.7, height = 1400, dpi = 100, bg = "white", plot = g1)


###############################################################################
########### GSEA enrichment plot using genekitr (Figure 2D)
################################################################################
# 1st step: prepare pre-ranked gene list and genekitr files
data(geneList, package = "genekitr")
head(geneList)

# 2nd step: prepare gene set
gs <- geneset::getMsigdb(org = "mouse",category = "C2-CP-REACTOME")
gs$geneset[gs$geneset$gs_name == "REACTOME_DNA_REPAIR",]
gs$geneset[gs$geneset$gs_name == "REACTOME_DNA_DAMAGE_BYPASS",]
gs$geneset[gs$geneset$gs_name == "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",]
gs$geneset[gs$geneset$gs_name == "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",]
gs$geneset <- gs$geneset[gs$geneset$gs_name %in% pathways,]

# 3rd step: GSEA analysis
gse <- genGSEA(genelist = deframe(ranks_lfc), geneset = gs, max_gset_size = 1000000, q_cutoff = 1 , min_gset_size = 0, p_cutoff = 1)
pathways <- c("REACTOME_DNA_REPAIR","REACTOME_DNA_DAMAGE_BYPASS","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_NUCLEOTIDE_EXCISION_REPAIR")
genes <- c("Tp53", "Atm", "Atr","Rad9","Rad1","Parp1","Brca1")
plot <- plotGSEA(gse, plot_type = "classic",
                 show_pathway = pathways,
                 #show_gene = genes,
                 main_text_size = 13,
                 #base.size = 15
)  #+ theme(legend.text = 100)
plot + theme(legend.text=element_text(size=30))
plot + theme(plot.title = element_text(size = 30, face = "bold"),
             legend.title=element_text(size=30), 
             legend.text=element_text(size=30))
print(plot)

################################################################################
########### GSVA (Supplementary Figure 1G)
################################################################################
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
dds_norm <- vst(f_dds_DE,blind = T)
vst_df <- as.data.frame(assay(dds_norm))
vst_df <- rownames_to_column(vst_df,"ensembl_id")

# First let's create a mapped data frame we can join to the gene expression values
keytypes(org.Mm.eg.db)
mapped_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Mm.eg.db,
    keys = vst_df$ensembl_id,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Ensembl" = "ensembl_id"))

# check for duplicates
sum(duplicated(mapped_df$entrez_id))
# First let's determine the gene means
gene_means <- rowMeans(dplyr::select(mapped_df,-Ensembl, -entrez_id))
# Add gene_means as a column called gene_means
mapped_df <- mutate(mapped_df,gene_means)
# Reorder the columns so `gene_means` column is upfront
mapped_df <- dplyr::select(mapped_df,Ensembl, entrez_id, gene_means, dplyr::everything())
# Sort so that the highest mean expression values are at the top
mapped_df <- arrange(mapped_df,desc(gene_means))
# Filter out the duplicated rows using `dplyr::distinct()`
filtered_mapped_df <- distinct(mapped_df,entrez_id, .keep_all = TRUE)
sum(duplicated(filtered_mapped_df$entrez_id))

# GSVA can't use the Ensembl IDs so we should drop this column as well as the means
filtered_mapped_df <- dplyr::select(filtered_mapped_df,-Ensembl, -gene_means)
# We need to store our gene identifiers as row names
filtered_mapped_df <- column_to_rownames(filtered_mapped_df,"entrez_id")
# Now we can convert our object into a matrix
filtered_mapped_matrix <- as.matrix(filtered_mapped_df)

#GSVA
gsva_results <- gsva(
  filtered_mapped_matrix,
  Mm.c2.all.v7.1.entrez,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 3,
  # Maximum gene set size
  max.sz = 2000,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

# Print 6 rows results
head(gsva_results[, 1:7])

# Show pathways of intrest
res <- gsva_results[c("AIGNER_ZEB1_TARGETS",
                      "REACTOME_CELL_CELL_COMMUNICATION",
                      "REACTOME_PI3K_AKT_ACTIVATION",
                      "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                      "REACTOME_GAP_JUNCTION_TRAFFICKING_AND_REGULATION",
                      "REACTOME_GAP_JUNCTION_DEGRADATION",
                      "REACTOME_APOPTOTIC_CLEAVAGE_OF_CELL_ADHESION_PROTEINS",
                      "REACTOME_SIGNALING_BY_WNT"),]

# add annotation column
Tissue.names <- c("11G5_1","11G5_2","11G5_3","11G5_4","Nissle_1","Nissle_2","Nissle_3")
colnames(res) <- Tissue.names 
Tissue <- c("11G5","11G5","11G5","11G5","Nissle 1917","Nissle 1917","Nissle 1917") 
anno_col <- data.frame(Tissue)
anno_col$Tissue <- factor(anno_col$Tissue)
rownames(anno_col) <- Tissue.names

# change rownames
names <- rownames(res)
names <- str_replace(names,"REACTOME_","")
names <- str_replace(names,"KEGG_","")
names <- str_replace(names,"AIGNER_","")
names <- str_replace_all(names,"_"," ")
rownames(res) <- names

breaks_list = seq(-0.3,0.3,by = 0.001)
pal <- colorpanel(length(breaks_list),"#6B92CB","white","#ED6E68")
pathway_heatmap <- pheatmap(res,
                            color = pal,
                            breaks = breaks_list,
                            border_color = "black",
                            cellwidth = 30,
                            cellheight = 30,
                            #scale="column",
                            cutree_cols = 2,
                            cluster_cols = T,
                            #annotation_col = anno_col,
                            #annotation_colors = list(Tissue = c(11G5 = "#6ef88a", Nissle = "#d357fe")),
                            #annotation_legend = F,
                            #annotation_names_col = F,
                            show_colnames = F, # Don't show sample labels
                            fontsize = 14, # Shrink the pathway labels a tad
                            fontfamily="Calibri"
)
ggsave("Data/Figures/Epithelial_11G5_TG_vs_Nissle_TG_GSVA.png", units = "px", width = 1920, height = 1043, dpi = 100, bg = "white", plot = pathway_heatmap)
