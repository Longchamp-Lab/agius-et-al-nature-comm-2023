## GO Gene Set Enrichment Analysis
library("tidyverse")
library(openxlsx)
library(ggplot2)
library(dplyr)
library(stringr)
library("clusterProfiler")
library("org.Mm.eg.db")
library(DOSE)
library(enrichplot)
library(grid)
library(stats)
library(ggrepel)

set.seed(45849489)

## Load file 
## need csv file with 2 columns, list of genes, foldchange
## 2 file needed
## gene file with yours DEG or any list of relevant genes
## univers file with all your geneset used for background ajustment 
## (by default it come from the OrgDb you used but better to use counts matrix file from your analysis)

liver <- read.csv("geneInfo_from_darkred_salmon4_merged.csv",)
kidney <- read.csv("geneInfo_from_maroon_paleturquoise_merged.csv")

gene <- list(paleturquoise = kidney %>% filter(moduleColor == "paleturquoise"),
             maroon = kidney %>% filter(moduleColor == "maroon"),
             darkred = liver %>% filter(moduleColor == "darkred"),
             salmon4 =liver %>% filter(moduleColor == "salmon4"))

gene[["paleturquoise"]]<- gene[["paleturquoise"]]$Symbol
gene[["maroon"]]<- gene[["maroon"]]$Symbol
gene[["darkred"]]<- gene[["darkred"]]$Symbol
gene[["salmon4"]]<- gene[["salmon4"]]$Symbol

univers_liver <- na.omit(read.csv("backgroundgene_liver_all.csv"))
univers_kidney <- na.omit(read.csv("backgroundgene_kidney_all.csv"))

names <- as.vector(names(gene))

###########################################-
## ClusterProfiler Analysis

### 1- GO over-representation analysis

for (i in 1:length(gene)) {
  if (names(gene)[i] == "paleturquoise" | names(gene)[i] == "maroon") {
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
go_paleturquoise <- simplify(goresult_paleturquoise,
                             cutoff = 0.7,
                             by = "p.adjust",
                             select_fun = min,
                             measure = "Wang",
                             semData = NULL)
go_maroon <- simplify(goresult_maroon,
                      cutoff = 0.7,
                      by = "p.adjust",
                      select_fun = min,
                      measure = "Wang",
                      semData = NULL)
go_salmon4 <- simplify(goresult_salmon4,
                       cutoff = 0.7,
                       by = "p.adjust",
                       select_fun = min,
                       measure = "Wang",
                       semData = NULL)
go_darkred <- simplify(goresult_darkred,
                       cutoff = 0.7,
                       by = "p.adjust",
                       select_fun = min,
                       measure = "Wang",
                       semData = NULL)


# Upsetplot showing over representation analysis

# paleturquoise
pdf("upsetplot_kidup_paleturquoise.pdf", width = 8, height = 4)
upsetplot(go_paleturquoise, n=10, font.size=14)
dev.off()

pdf("dotplot_kidup_paleturquoise.pdf", width = 6, height = 6)
dotplot(go_paleturquoise, title ="GOTerm MPaleturquoise - Kidney Up", font.size = 16, showCategory = 10, label_format = 30)
dev.off()

# maroon
pdf("upsetplot_kidup_maroon.pdf", width = 8, height = 4)
upsetplot(go_maroon, n=10, font.size=14)
dev.off()

pdf("dotplot_kidup_maroon.pdf", width = 6, height = 6)
dotplot(go_maroon, title ="GOTerm Mmaroon - Kidney Up", font.size = 16, showCategory = 10, label_format = 30)
dev.off()

# salmon4
pdf("upsetplot_livup_salmon4.pdf", width = 11, height = 4)
upsetplot(go_salmon4, n=8, font.size=14)
dev.off()

pdf("dotplot_livup_salmon4.pdf", width = 7.5, height = 6)
dotplot(go_salmon4, title ="GOTerm Msalmon4 - Liver Up", font.size = 16, showCategory = 10, label_format = 30)
dev.off()

# darkred
pdf("upsetplot_livup_darkred.pdf", width = 8, height = 4)
upsetplot(go_darkred, n=10, font.size=14)
dev.off()

pdf("dotplot_livup_darkred.pdf", width = 6, height = 6)
dotplot(go_darkred, title ="GOTerm Mdarkred - Liver Up", font.size = 16, showCategory = 10, label_format = 30)
dev.off()


write.csv(ego.filt,"GO Analysis/enrichGO_filtered_result_ORA.csv") # print out result
write.csv(ego,"enrichGO_result.csv") # print out result

## Meta analysis with euleR venn diagram with common genes between differents approaches
library(VennDiagram)



# function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

venndata <- list(
  A = na.omit(union(go_paleturquoise@result$ID,go_maroon@result$ID)),
  B = na.omit(union(go_darkred@result$ID,go_salmon4@result$ID))
)

venndata <- list(
  A = na.omit(go_darkred@result$ID),
  B = na.omit(go_salmon4@result$ID)
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
write.csv(go2term(overlap$a3), "goterm_overlap_list_livervskidney.csv")


## Scatterplot GOTerm overlap kidney liver with pvalue

goterm.liver <- merge(go_darkred@result, go_salmon4@result, all = T)
goterm.kidney <- merge(go_maroon@result, go_paleturquoise@result, all = T)

goterm.liver.filt <- dplyr::filter(goterm.liver[!duplicated(goterm.liver$ID),], ID %in% as.vector(overlap$a3))
goterm.liver.filt$type <- "liver"
goterm.kidney.filt <- dplyr::filter(goterm.kidney[!duplicated(goterm.kidney$ID),], ID %in% as.vector(overlap$a3))
goterm.kidney.filt$type <- "kidney"

goterm.merge <- merge(goterm.liver.filt, goterm.kidney.filt, by = "ID" ,all = T)

write.csv(goterm.merge, "goterm_commonliver_kidney.csv")

ggplot(goterm.merge, aes(x = -log10(p.adjust.x), y = -log10(p.adjust.y))) +
  geom_point(aes(size=Count.x, colour = p.adjust.y),
             show.legend = FALSE) +
  geom_text_repel(data = subset(goterm.merge, (-log10(p.adjust.x) > 9 | -log10(p.adjust.y) > 9)),
    aes(label = Description.x)) +
  geom_hline(yintercept = 8, linetype = "dashed", size = .5) + 
  geom_vline(xintercept = 8, linetype = "dashed", size = .5) +
  theme_bw()
