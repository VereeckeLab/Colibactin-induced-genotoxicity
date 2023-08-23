# Store DE results
save_all_DE_results <- function(dds_DE , #= filtered_subset_dds_DE
                                dds_res , #= subset_results
                                save_dir , #= subset_dir
                                padj_cutoff = 0.05,
                                annotation , #= annotation
                                name_file = name ){
  # Convert to table and add rownames to column with name "ensembl_id"
  results_table <- as_tibble(dds_res, rownames = "ensembl_id")
  # add annotation
  results_table <- left_join(results_table,annotation, by = c("ensembl_id" = "Gene.stable.ID"))
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds_DE, normalized = F)
  cpms_table <- enframe(rowMeans(cpm(raw_counts)))
  colnames(cpms_table) <- c("ensembl_id", "avg_raw_cpm")
  # add avg cpms to table
  results_table <- left_join(results_table,cpms_table, by = c("ensembl_id" = "ensembl_id"))
  
  # Remove all NA from gene.name because this can give error when fetching entrez_id
  results_table <- results_table[!is.na(results_table$Gene.name),]
  # Add tag if gene is a DE gene
  results_table <- mutate(results_table, DE_gene = case_when(padj <= padj_cutoff ~ TRUE, T ~ FALSE))
  # add entrez id from GSEA
  results_table[(results_table$Gene.name == ""),"Gene.name"] <- "Not_found"
  # Get entrez ids from database
  entrez_ids <- results_table$Gene.name
  names(entrez_ids) <- mget(results_table$Gene.name,revmap(org.Mm.egSYMBOL),ifnotfound = NA)
  results_table[["Entrez_GSEA"]] <- names(entrez_ids) 
  # check for NA columns
  colSums(is.na(results_table))
  # Store all unfiltered data
  all_data <- results_table
  
  # Remove empty string gene names ("")
  results_table <- results_table[(results_table$Gene.name != "Not_found"),]
  # Remove duplicate gene names
  results_table <- results_table[!duplicated(results_table$Gene.name), ]
  # Remove non protein coding genes
  results_table <- results_table[(results_table$Gene.type == "protein_coding"),]
  # Remove all NA from padj 
  results_table <- results_table[!is.na(results_table$padj),]
  colSums(is.na(results_table))
  # Remove genes with different entrez id's
  dim(results_table)
  results_table <- results_table[(results_table$Entrez.id == results_table$Entrez_GSEA),]
  dim(results_table)
  # Store all filtered data
  filtered_all_data <- na.omit(results_table)
  colSums(is.na(filtered_all_data))
  
  
  dir.create(save_dir, showWarnings = FALSE)
  write.csv (all_data, file = paste0(save_dir,"/",name_file,".csv"), row.names = F)
  write.csv (filtered_all_data, file = paste0(save_dir,"/",name_file,"_filtered.csv"), row.names = F)
  write_xlsx(all_data, path = paste0(save_dir,"/",name_file,".xlsx"))
  write_xlsx(filtered_all_data, path = paste0(save_dir,"/",name_file,"_filtered.xlsx"))
}

# Volcano Plot
get_paper_volcano_plot <- function(result_table = results,
                                   important_genes,
                                   title="test",
                                   crop_value_cutoff = 1e-10,
                                   padj_cutoff = 0.05,
                                   logFC_cutoff = 1,
                                   gene_num = 10,
                                   colors = c("#00FFFF", "gray50", "#ff1493"),
                                   point_size = 2.5,
                                   text_title_size = 20,
                                   figure_path = "",
                                   figure_name = ""){
  # Adding upregulated or downregulated tag
  library(ggrepel)
  
  n_up_genes <- nrow(result_table[(result_table$padj<=padj_cutoff) & (result_table$log2FoldChange >= logFC_cutoff),])
  n_down_genes <- nrow(result_table[(result_table$padj<=padj_cutoff) & (result_table$log2FoldChange <= -(logFC_cutoff)),])
  fulldata <- result_table
  result_table <- mutate(result_table,Expression = case_when(log2FoldChange >= logFC_cutoff & padj <= padj_cutoff ~ paste0("Up-regulated"),
                                                             log2FoldChange <= -(logFC_cutoff) & padj <= padj_cutoff ~ paste0("Down-regulated"),
                                                             T ~ "Unchanged"))
  
  result_table <- mutate(result_table,Newpvalue = case_when(pvalue > crop_value_cutoff  ~ pvalue,
                                                            T ~ crop_value_cutoff))
  
  result_table <- mutate(result_table,View = case_when(pvalue > crop_value_cutoff  ~ "Outside",
                                                       T ~ "Inside"))
  
  fulldata <- mutate(fulldata,Expression = case_when(log2FoldChange >= logFC_cutoff & padj <= padj_cutoff ~ colors[3],
                                                     log2FoldChange <= -(logFC_cutoff) & padj <= padj_cutoff ~ colors[1],
                                                     T ~ colors[2]))
  
  fulldata <- mutate(fulldata,Newpvalue = case_when(pvalue > crop_value_cutoff  ~ pvalue,
                                                    T ~ crop_value_cutoff))
  
  fulldata <- mutate(fulldata,View = case_when(pvalue > crop_value_cutoff  ~ "Outside",
                                               T ~ "Inside"))
  
  # Filter on padj before doing this then order and take top X
  significant_genes <- result_table[result_table$padj < padj_cutoff,]
  significant_genes_ordered <- significant_genes[order(significant_genes$log2FoldChange,decreasing = T),]
  up <- head(significant_genes_ordered,round(gene_num/2,0))
  down <- tail(significant_genes_ordered,round(gene_num/2,0))
  top_genes <- rbind(up,down)
  top_genes_names <- top_genes$Gene.name
  
  # Get important genes
  found_genes <- c()
  for (gene in important_genes){
    if (gene %in% fulldata$Gene.name){
      found_genes <- append(found_genes,gene)
    }
  }
  data <- result_table[result_table$Gene.name == "nothing",]
  for (gene in found_genes){
    gene_info <- fulldata[fulldata$Gene.name == gene,]
    data <- rbind(data,gene_info)
  }
  top_genes_names <- data$Gene.name
  expr <- data$Expression
  
  volcano_plot <- ggplot(result_table, aes(log2FoldChange, -log(Newpvalue,10))) + # -log10 conversion  
    geom_point(aes(color = Expression, shape = View) , size = point_size) +
    scale_color_manual(breaks = c("Down-regulated","Up-regulated") , values = c("Down-regulated" = colors[1], "Up-regulated" = colors[3], "Unchanged" = colors[2]) ) +
    scale_shape_manual(values=c(20, 16), labels = c("Outside","Inside"))+ # Causes bug in frame!
    geom_label_repel(data = data,
                     mapping = aes(log2FoldChange, -log(Newpvalue,10), label = top_genes_names ,),
                     max.overlaps = 50, nudge_x = 0.5,
                     size = point_size + 2) +
    ggtitle(title) +
    ylim(0,-log(crop_value_cutoff,10)) +
    xlab(expression(bold("log"[2]*"FC"))) + 
    ylab(expression(bold("-log"[10]*"p-value"))) +
    guides(shape = "none") +  # dont show shape in legend
    theme_classic() +
    theme(text = element_text(family="Calibri"),
          axis.title = element_text(face = "bold", color = "black",size = 23,hjust = 0.5),
          axis.text = element_text(color = "black", size = 17),
          plot.title = element_text(face = "bold", color = "black",size = 25,hjust = 0.5),
          legend.text = element_text(size = 21),
          legend.position = c(0.8, 0.9), # OPTIONAL
          #legend.background = element_rect(size=0.5, linetype="solid",colour ="black"), # OPTIONAL
          legend.title = element_text(size = 27, face = "bold"))
  
  
  print(volcano_plot)
  ggsave(filename=figure_name,path = figure_path,height=2036,width=1855,units="px",dpi=200,limitsize=FALSE)
  return(volcano_plot)
}

# GSEA plots
GSEA_plot <- function (fgsea_results,
                       NES_cutoff = 1.5,
                       npathw = 20,
                       color_down = "#00FFFF",
                       color_up = "#ff1493",
                       title = "No title",
                       reverse_order = T,
                       y_text_size = 7 ) {
  n <- npathw/2
  fgsea_results <- mutate(fgsea_results, short_name = str_split_fixed(fgsea_results$pathway,"_",2)[,2])
  fgsea_results$short_name <- str_replace_all(fgsea_results$short_name,"_"," ")
  print(fgsea_results$short_name)
  fgsea_results$short_name <- str_replace_all(fgsea_results$short_name," UP","")
  fgsea_results$short_name <- str_replace_all(fgsea_results$short_name," DN","")
  fgsea_results <- mutate(fgsea_results, Regulated = case_when(NES > 0 ~ "Up-regulated",
                                                               NES < 0 ~ "Down-regulated",
                                                               T ~ "Unchanged"))
  if (reverse_order == T){
    up <- fgsea_results[fgsea_results$Regulated == "Up-regulated",]
    up <- up[order(up$NES, decreasing = F),]
    
    down <- fgsea_results[fgsea_results$Regulated == "Down-regulated",]
    down <- down[order(down$NES, decreasing = F),]
    
    top20 <- rbind(tail(up,n),head(down,n))
  } else {
    up <- fgsea_results[fgsea_results$Regulated == "Up-regulated",]
    up <- up[order(up$NES, decreasing = F),]
    
    down <- fgsea_results[fgsea_results$Regulated == "Down-regulated",]
    down <- down[order(down$NES, decreasing = F),]
    
    top20 <- rbind(head(down,n),tail(up,n))
  }
  
  
  # Check for short names already present and remove them for now
  top20 <- top20[duplicated(top20$short_name) == F,]
  
  # lock in factor level order
  top20$short_name <- factor(top20$short_name, levels = top20$short_name)
  
  p <- ggplot(top20,aes(short_name, NES)) +
    geom_bar(stat= "identity", aes(fill = Regulated))+
    scale_fill_manual(values=c(color_down, color_up)) +
    coord_flip() +
    #scale_y_continuous(limits = c(-2, 2)) +
    geom_hline(yintercept = 0) +
    labs(x = "", y = "Normalized Enrichment Score")+
    theme(text = element_text(family = "Calibri",color = "black"),
          #axis.text.y = element_text(size = y_text_size), 
          #plot.title = element_text(hjust = 1),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_rect(fill = "white"),
          #axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black",size = 18,hjust = 0.5,face = "bold"),
          axis.text.y = element_text(color = "black", size = 21),
          axis.text.x = element_text(color = "black", size = 16),
          plot.title = element_text(face = "bold", color = "black",size = 23,hjust = 0.5), #0.5
          #legend.position = "None",
          legend.position = "right",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 21, face = "bold")
    ) +
    ggtitle(title)
  return(p)
}