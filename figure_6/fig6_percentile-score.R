#=====================================================================================
#Heatmap of gene expression from extracted module
#=====================================================================================
library(DESeq2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(org.Mm.eg.db)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(stats)
library(ggrepel)
library(viridis)

set.seed(78974564)

# Load geneInfo from WGCNA

#Load data from csv
kid_gene_up <- read.csv("geneInfo_from_maroon_paleturquoise_merged.csv")
kid_gene_down <- read.csv("geneInfo_from_darkseagreen4_merged.csv")
kid_gene_up$X <- mapIds(org.Mm.eg.db, kid_gene_up$ensembl.id, keytype = "ENSEMBL", column = "SYMBOL")
kid_gene_down$X <- mapIds(org.Mm.eg.db, kid_gene_down$ensembl.id, keytype = "ENSEMBL", column = "SYMBOL")


liv_gene_up <- read.csv("geneInfo_from_darkred_salmon4_merged.csv")
liv_gene_down <- read.csv("geneInfo_from_darkorange_indianred4_navajowhite2_merged.csv")
liv_gene_up$X <- mapIds(org.Mm.eg.db, liv_gene_up$ensembl.id, keytype = "ENSEMBL", column = "SYMBOL")
liv_gene_down$X <- mapIds(org.Mm.eg.db, liv_gene_down$ensembl.id, keytype = "ENSEMBL", column = "SYMBOL")


### Load MAtrix count from liver and Kidney log Normalized data
liv_datExpr <- read.csv("counts_logCPM_liver.csv",row.names=1)
kid_datExpr <- read.csv("counts_logCPM_kidney.csv",row.names=1)

liv_datExpr <- subset(liv_datExpr, select = -c(title, rho, SD, AVG))
kid_datExpr <- subset(kid_datExpr, select = -c(title, rho, SD, AVG))

## Kidney
##########
#Kidney

# Function to compute percentile transformation
percentile_transformation <- function(logCPM_matrix) {
  
  # Calculate the percentile rank for each value in the matrix
  percentile_matrix <- t(apply(logCPM_matrix, 1, function(x) {
    ranks <- rank(x, na.last = "keep", ties.method = "average")
    percentiles <- (ranks - 1) / (length(x) - 1)
    return(percentiles)
  }))
  
  return(percentile_matrix)
}

kid_datExpr$gene <- rownames(kid_datExpr)

# Filter top gene kidney
  
# Arrange the data by moduleColor and descending order of GS.fgf21.serum
kid_gene_up <- kid_gene_up[order(kid_gene_up$moduleColor, -kid_gene_up$GS.fgf21.serum), ]
kid_gene_down <- kid_gene_down[order(kid_gene_down$moduleColor, -kid_gene_down$GS.fgf21.serum), ]

# Group the data by moduleColor and select the top 10 rows for each group
top_ensembl_ids_up <- kid_gene_up %>%
  na.omit() %>%
  group_by(moduleColor) %>%
  dplyr::filter(GS.fgf21.serum > 0.6) %>%
  slice_head(n = 15) %>%
  dplyr::select(ensembl.id, Symbol, moduleColor, GS.fgf21.serum)

top_ensembl_ids_down <- kid_gene_down %>%
  na.omit() %>%
  group_by(moduleColor) %>%
  dplyr::filter(GS.fgf21.serum < -0.6) %>%
  slice_head(n = 15) %>%
  dplyr::select(ensembl.id, Symbol, moduleColor, GS.fgf21.serum)


# write csv file for Regulon analysis
top_ensembl_ids_ireg <- kid_gene_up %>%
  na.omit() %>%
  group_by(moduleColor) %>%
  dplyr::filter(GS.fgf21.serum > 0.6) %>%
  dplyr::select(ensembl.id, Symbol, moduleColor, GS.fgf21.serum, MM.maroon, MM.paleturquoise)

write.csv(top_ensembl_ids_ireg, "kid_top_genes.csv")


#merge up and down genes
top_ensembl_ids <- rbind(top_ensembl_ids_up,top_ensembl_ids_down)

kid_modGenes <- kid_datExpr %>% dplyr::filter(gene %in% top_ensembl_ids$Symbol)

kid_modGenes <- subset(kid_modGenes, select = - gene)

all.group <- as.data.frame(colnames(kid_modGenes))
all.group$group[grepl("Ctrl", all.group[,1], fixed = TRUE)] = "Ctrl"
all.group$group[grepl("HC", all.group[,1], fixed = TRUE)] = "HC"
all.group$group[grepl("LP", all.group[,1], fixed = TRUE)] = "LP"
all.group$group[grepl("LPHC", all.group[,1], fixed = TRUE)] = "LPHC"
all.group <- as.vector(all.group$group)

# Perform percentile transformation
percentile_matrix <- percentile_transformation(kid_modGenes)
mean_percentile_matrix <- apply(percentile_matrix, 1, function(x)by(x, all.group, mean))
mean_kid_modGenes <- mean_percentile_matrix

####

# heatmap of gene expression from the selected module
module <- data.frame(top_ensembl_ids$moduleColor[top_ensembl_ids$Symbol %in% colnames(mean_kid_modGenes)])

row.names(module) <- colnames(mean_kid_modGenes)
colnames(module)[1] <- "module"
my_colour = list("module" = 
                   c("maroon" = "maroon", 
                     "paleturquoise" = "paleturquoise", 
                     "darkseagreen4" = "darkseagreen4"))

pdf("Heatmap_module_gene_kid_top15_each module_ percentile_transform.pdf",onefile=T,width=13,height=6)
mybreaks <- c(
  seq(0, 1,length.out=100)
) 
color <- viridis(99) 
pheatmap(mean_kid_modGenes,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = color,
         annotation_col = module, 
         annotation_colors = my_colour, 
         annotation_legend = TRUE,
         border_color = NA,
         cellheight = 15,
         cellwidth = 12,
         fontsize_col = 12
)
dev.off()


## Liver
##########
# Liver

#Declare z-score function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# Function to compute percentile transformation
percentile_transformation <- function(logCPM_matrix) {
  
  # Calculate the percentile rank for each value in the matrix
  percentile_matrix <- t(apply(logCPM_matrix, 1, function(x) {
    ranks <- rank(x, na.last = "keep", ties.method = "average")
    percentiles <- (ranks - 1) / (length(x) - 1)
    return(percentiles)
  }))
  
  return(percentile_matrix)
}


liv_datExpr$gene <- rownames(liv_datExpr)

# Filter top gene kidney

# Assuming your dataset is stored in a variable called 'data'
# Arrange the data by moduleColor and descending order of GS.fgf21.serum
liv_gene_up <- liv_gene_up[order(liv_gene_up$moduleColor, -liv_gene_up$GS.FGF21), ]
liv_gene_down <- liv_gene_down[order(liv_gene_down$moduleColor, -liv_gene_down$GS.FGF21), ]

# Group the data by moduleColor and select the top 10 rows for each group
top_ensembl_ids_up <- liv_gene_up %>%
  na.omit() %>%
  group_by(moduleColor) %>%
  #slice_head(n = 10) %>%
  dplyr::filter(GS.FGF21 > 0.8) %>%
  slice_head(n = 15) %>%
  dplyr::select(ensembl.id, Symbol, moduleColor, GS.FGF21)


top_ensembl_ids_down <- liv_gene_down %>%
  na.omit() %>%
  group_by(moduleColor) %>%
  #slice_head(n = 10) %>%
  dplyr::filter(GS.FGF21 < -0.75) %>%
  slice_head(n = 15) %>%
  dplyr::select(ensembl.id, Symbol, moduleColor, GS.FGF21)


# write csv file for Regulon analysis
top_ensembl_ids_ireg <- liv_gene_up %>%
  na.omit() %>%
  group_by(moduleColor) %>%
  #slice_head(n = 10) %>%
  dplyr::filter(GS.FGF21 > 0.8) %>%
  #slice_head(n = 25) %>%
  dplyr::select(ensembl.id, Symbol, moduleColor, GS.FGF21, MM.darkred, MM.salmon4)

write.csv(top_ensembl_ids_ireg, "liver_top_genes.csv")

#merge up and down genes
top_ensembl_ids <- rbind(top_ensembl_ids_up,top_ensembl_ids_down)


liv_modGenes <- liv_datExpr %>% dplyr::filter(gene %in% top_ensembl_ids$Symbol)

liv_modGenes <- subset(liv_modGenes, select = - gene)

all.group <- as.data.frame(colnames(liv_modGenes))
all.group$group[grepl("PER_0", all.group[,1], fixed = TRUE)] = "pct 0%"
all.group$group[grepl("PER_2", all.group[,1], fixed = TRUE)] = "pct 2%"
all.group$group[grepl("PER_6", all.group[,1], fixed = TRUE)] = "pct 6%"
all.group$group[grepl("PER_10", all.group[,1], fixed = TRUE)] = "pct 10%"
all.group$group[grepl("PER_14", all.group[,1], fixed = TRUE)] = "pct 14%"
all.group$group[grepl("PER_18", all.group[,1], fixed = TRUE)] = "pct 18%"
all.group <- as.vector(all.group$group)

# Perform percentile transformation
percentile_matrix <- percentile_transformation(liv_modGenes)
mean_percentile_matrix <- apply(percentile_matrix, 1, function(x)by(x, all.group, mean))
mean_liv_modGenes <- mean_percentile_matrix

# Sort the rows by row names (pct %)
mean_liv_modGenes <- mean_liv_modGenes[order(as.numeric(gsub("%","", gsub("pct ", "", rownames(mean_liv_modGenes)))), decreasing = T),]

# heatmap of gene expression from the selected module
module <- data.frame(top_ensembl_ids$moduleColor[top_ensembl_ids$Symbol %in% colnames(mean_liv_modGenes)])

row.names(module) <- colnames(mean_liv_modGenes)
colnames(module)[1] <- "module"
my_colour = list("module" = 
                   c("darkred" = "darkred", 
                     "salmon4" = "salmon4", 
                     "darkorange" = "darkorange", 
                     "indianred4" = "indianred4", 
                     "navajowhite2" = "navajowhite2"))


pdf("Heatmap_module_gene_liv_top15_each module_2_percentile transform.pdf",onefile=T,width=14,height=6)
mybreaks <- c(
  #seq(0, -0.01, length.out=50),
  seq(0,1,length.out=100)
) 
color <- viridis(99)
pheatmap(mean_liv_modGenes,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = color,
         breaks = mybreaks,
         annotation_col = module, 
         annotation_colors = my_colour, 
         annotation_legend = TRUE,
         border_color = NA,
         cellheight = 15,
         cellwidth = 12,
         fontsize_col = 12
)
dev.off()
