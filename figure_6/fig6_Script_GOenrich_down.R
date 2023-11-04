## GO Gene Set Enrichment Analysis
library("tidyverse")
library(ggplot2)
library(dplyr)
library("clusterProfiler")
library("org.Mm.eg.db")
library(DOSE)
library(enrichplot)
library(stats)
library(ggrepel)
library(grid)

set.seed(45849489)

## Load file 
## need csv file with 2 columns, list of genes, foldchange
## 2 file needed
## gene file with yours DEG or any list of relevant genes
## univers file with all your geneset used for background ajustment 
## (by default it come from the OrgDb you used but better to use counts matrix file from your analysis)

liver <- read.csv("geneInfo_from_darkorange_indianred4_navajowhite2_merged.csv",)
kidney <- read.csv("geneInfo_from_darkseagreen4_merged.csv")

#filter pval fgf21 < .05 + mmbership
filt_liver <- liver #%>% dplyr::filter(GS.FGF21 < -0.4 & (MM.darkorange > .2 | MM.indianred4 > .2 | MM.navajowhite2 > .2))
filt_kidney <- kidney #%>% dplyr::filter(GS.fgf21.serum < -0.4 & MM.darkseagreen4 > .2)

moduleColors.genes <- c(names(table(liver$moduleColor)), names(table(kidney$moduleColor)))

gene <- list()
for (color in moduleColors.genes) {
  gene[[color]] <- rbind(filt_liver %>% filter(moduleColor == color), filt_kidney %>% filter(moduleColor == color))
  gene[[color]] <- gene[[color]]$Symbol
}

univers_liver <- na.omit(read.csv("backgroundgene_liver_all.csv"))
univers_kidney <- na.omit(read.csv("backgroundgene_kidney_all.csv"))

names <- as.vector(names(gene))

###########################################-
## ClusterProfiler Analysis

### 1- GO over-representation analysis

for (i in 1:length(gene)) {
  if (names(gene)[i] %in% names(table(kidney$moduleColor))) {
    assign(paste0("goresult_", names(gene)[i]),
           enrichGO(gene = gene[[i]],
                    universe      = univers_kidney$x,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    keyType = "SYMBOL",
                    readable      = TRUE))
  } else {
    assign(paste0("goresult_", names(gene)[i]),
           enrichGO(gene = gene[[i]],
                    universe      = univers_liver$x,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    keyType = "SYMBOL",
                    readable      = TRUE))
  }
  
}

## simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms
results <- list()
for (i in seq_along(gene)) {
  result <- simplify(get(paste0("goresult_", names(gene)[i])),
                     cutoff = 0.7,
                     by = "p.adjust",
                     select_fun = min,
                     measure = "Wang",
                     semData = NULL)
  results[[names(gene)[i]]] <- result
}


for (i in seq_along(results)) {
  color <- names(results)[i]
  pdf(paste0("dotplot__", color, ".pdf"), width = 7, height = 5)
  print(dotplot(results[[color]], title = paste("GO Term M", color), font.size = 16, showCategory = 10, label_format=30))
  dev.off()
}

for (i in seq_along(results)) {
  color <- names(results)[i]
  pdf(paste0("upsetplot__", color, ".pdf"), width = 8, height = 4)
  print(upsetplot(results[[color]],n=10, font.size=14))
  dev.off()
}

for (i in seq_along(results)) {
  color <- names(results)[i]
  filename <- paste0("enrichGO_filtered_result_", color, ".csv")
  write.csv(results[[color]], file = filename)
}



## Meta analysis with euleR venn diagram with common genes between different approaches
library(VennDiagram)

# function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

venndata <- list(
  A = na.omit(results[[4]]@result$ID),
  B = na.omit(union(union(results[[1]]@result$ID,results[[2]]@result$ID),results[[3]]@result$ID))
)

display_venn(venndata,
             category.names = c("GOTerm Kidney" , "GOTerm Liver"),
             #Circles
             lwd = 2,
             fill = c("#ffffff", "#858585"),
             # Numbers
             cex = 1,
             fontfamily = "arial",
             fontface = "italic",
             print.mode = c("raw", "percent"),
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.fontfamily = "arial",
             cat.default.pos = "outer",
             # Other
             euler.d = TRUE,
             scaled = TRUE,
             ext.text = TRUE
)


overlap <- calculate.overlap(venndata)
write.csv(go2term(overlap$a3), "goterm_overlap_down_list_livervskidney_filt.csv")

#recover Info from pathways from liver
# Convert each list item to a data frame
df1 <- as.data.frame(results[[1]]@result)
df2 <- as.data.frame(results[[2]]@result)
df3 <- as.data.frame(results[[3]]@result)


# Remove any rows with missing values (NA) using na.omit()
df1_clean <- na.omit(df1)
df2_clean <- na.omit(df2)
df3_clean <- na.omit(df3)


# Combine the cleaned data frames using rbind()
liver_ovlp <- rbind(df1_clean, df2_clean, df3_clean)
liver_ovlp <- liver_ovlp[liver_ovlp$ID %in% overlap$a3,]

#recover Info from pathways from kidney
df4 <- as.data.frame(results[[4]]@result)
df4_clean <- na.omit(df4)
kidney_ovlp <- df4_clean
kidney_ovlp <- kidney_ovlp[kidney_ovlp$ID %in% overlap$a3,]



## Scatterplot GOTerm overlap kidney liver with pvalue
goterm.liver <- rbind(results[[1]]@result, results[[2]]@result, results[[3]]@result)
goterm.kidney <- results[[4]]@result

goterm.liver.filt <- dplyr::filter(goterm.liver[!duplicated(goterm.liver$ID),], ID %in% as.vector(overlap$a3))
goterm.liver.filt$type <- "liver"
goterm.kidney.filt <- dplyr::filter(goterm.kidney[!duplicated(goterm.kidney$ID),], ID %in% as.vector(overlap$a3))
goterm.kidney.filt$type <- "kidney"

goterm.merge <- merge(goterm.liver.filt, goterm.kidney.filt, by = "ID" ,all = T)

write.csv(goterm.merge, "goterm_common_down_liver_kidney.csv")

ggplot(goterm.merge, aes(x = -log10(p.adjust.x), y = -log10(p.adjust.y))) +
  geom_point(aes(size=Count.x, colour = p.adjust.y),
             show.legend = T) +
  geom_text_repel(data = subset(goterm.merge, (-log10(p.adjust.x) > 2.6 | -log10(p.adjust.y) > 2.6)),
    aes(label = Description.x)) +
  geom_hline(yintercept = 2, linetype = "dashed", size = .5) + 
  geom_vline(xintercept = 2, linetype = "dashed", size = .5) +
  theme_bw()
