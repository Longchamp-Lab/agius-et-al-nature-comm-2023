library(WGCNA)
library(DESeq2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(genefilter)
library(org.Mm.eg.db)
library(flashClust)

#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================

# Increase memory to facilitate integration:
options(future.globals.maxSize = 8000 * 1024^2)
allowWGCNAThreads()

set.seed(78974564)

## Load data with edgeR
# Read the gene counts table 
data0 <- read.csv("../fig6_counts.csv",row.names=1) #data0 <- read.csv("counts_logCPM_ungrouped.csv",row.names=1)

# 4 diet only
data0 <- data0[,!grepl("IGF1", colnames(data0), fixed = TRUE)]
data0 <- data0[,!grepl("Ctrl3", colnames(data0), fixed = TRUE)]
data0 <- data0[,!grepl("Ctrl5", colnames(data0), fixed = TRUE)]

## Normalize data
# Only keep rows that have total counts above the cutoff
de_input = data0 %>% dplyr::filter(rowSums(.) >= 50)

de_input = as.matrix(de_input)
  

meta_df <- data.frame( Sample = names(data0)) %>%
  mutate(
    Group = gsub("[0-9]","", Sample) %>% gsub("[_].*","", .)
  )

dds <- DESeqDataSetFromMatrix(round(de_input),
                              meta_df,
                              design = ~1)

## remove all genes with counts < 15 in more than 75% of samples (24*0.75)
## suggested by WGCNA on RNAseq FAQ
## 75% - maintains the TPS gene in dataset. balances low counts/false positives
dds75 <- dds[rowSums(counts(dds) >= 15) >= length(colnames(de_input))*0.75,]
nrow(dds75) # 14849 genes

## Perform the DESeq2 normalization, required before WGCNA analysis
dds2 <- DESeq(dds75)
## Perform a variance-stabilizing transformation
vsd <- getVarianceStabilizedData(dds2)
geneNames = rownames(vsd)
## Many functions expect the matrix to be transposed
datExpr <- t(vsd) 

gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");
# plot sample tree
pdf(file = "1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=50, by=2))
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers,
                        networkType = "signed hybrid",
                        verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate #power=7

#TOM = TOMsimilarityFromExpr(datExpr, power = power)
#dissTOM = 1-TOM 

#adjacency=adjacency(datExpr, power=softPower, type="signed hybrid")
adjacency = adjacency(datExpr, power = power, #nThreads = 22,
                      type = "signed hybrid");
                      #corFnc = "bicor", corOptions = "use = 'p', maxPOutliers = 0.1");

save(sft, adjacency, file = "adjacency.signedhybrid.softpwer.RData")

## Topological Overlap Matrix (TOM)
# Turn adjacency into topological overlap, i.e. translate the adjacency into 
# topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType = "signed", verbose = 5);
dissTOM = 1-TOM;

#===============================================================================
#
#  Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
## generate a clustered gene tree
geneTree <- flashClust(as.dist(dissTOM), method="average")

sizeGrWindow(8, 6)
pdf(file = "3-gene_cluster_alt.pdf", width = 12, height = 9);
plot(geneTree, xlab = "", sub = "", main = "Gene Clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# This sets the minimum number of genes to cluster into a module
minModuleSize <- 30 
# generating modules and assigning them colors
# Module identification using dynamic tree cut: 
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)
table(dynamicMods)

## Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
## Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf(file = "4-module_tree_alt.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================

## Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs) #or bicor
# Cluster module eigengenes
METree <- flashClust(as.dist(MEDiss), method= "average")


#plots tree showing how the eigengenes cluster together
sizeGrWindow(15,10)
MEDissThres <- 0.4
plot(METree, main= "Clustering of module eigengenes", xlab = "", sub = "") %>%
  abline(h = MEDissThres, col = "red")

pdf(file = "5-Module_Dendro.pdf", wi = 15, he = 10)
par(mar = c(.75, 2.75, 0.6, 1) + 0.1)            # The default is ‘c(5, 4, 4, 2) + 0.1’ c(bottom, left, top, right)’.
plot(METree, main= "Clustering of module eigengenes", xlab = "", sub = "") %>%
  abline(h = MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, 
                           corFnc = "bicor", verbose = 5)
# The merged module colors
mergedColors <- merge$colors
length(unique(mergedColors))
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
length(mergedMEs)

# To see the merging upon module colors, plot the gene dendrogram again, with
# the original and merged colors underneath
pdf(file = "5-merged_geneDendro.pdf", wi = 6, he = 4)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.rowText = 1.3)
dev.off()

write.table(merge$oldMEs,file="oldMEs.txt");
write.table(merge$newMEs,file="newMEs.txt");

## In subsequent analysis, the merged module colors in mergedColors will be used
# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresonding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
dim(MEs)

#===============================================================================
#
#  Plot the heatmap of module eigen-genes and samples
#
#===============================================================================

# Heatmap of old module eigen-genes and samples
pdf(file="oldMEs.pdf",heigh=10,width=10)
row.names(merge$oldMEs)=names(data0)
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=12)
dev.off()
# Heatmap of new module eigen-genes and samples
pdf(file="newMEs.pdf",heigh=10,width=10)
row.names(merge$newMEs)=names(data0)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=12)
dev.off()

#=====================================================================================
#
#  Correlation between gene modules and traits
#
#=====================================================================================

## Correlate eigengenes with external traits and look for the most significant associations
# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# Recalculate MEs with color labels
MEsO <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEsO)
# Read data as traits
traits = read.table("fig6_traits.csv", header = T, sep = ",")
rownames(traits) = traits[, 1]
traits = traits[, -c(1)]
traits = traits[, -c(1)]

# sample names should be consistent in eigen genes and traits !!!!
traits = traits[match(rownames(MEs), rownames(traits)), ]
table(rownames(MEs) == rownames(traits))
moduleTraitCor <- bicor(MEs, traits, use = "p") #bicor, cor
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");

#=====================================================================================
#
#  Plot heatmap of module-traits relationship
#
#=====================================================================================

# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(datExpr, traits$fgf21.serum , use = 'p')
gene.signf.corr.pvals <- data.frame(corPvalueStudent(gene.signf.corr, nSamples))

gene.signf.corr.pvals$symbol <- mapIds(org.Mm.eg.db, rownames(gene.signf.corr.pvals),
                                       keytype = "ENSEMBL", column = "SYMBOL")

write.csv(gene.signf.corr.pvals %>% 
            as.data.frame() %>% 
            arrange(corPvalueStudent.gene.signf.corr..nSamples.), "gene.signf.corr.pval.csv")

# Will display correlations and their p-values
pdf(file = "mod-trait.relations.pdf", width = 12, height = 9)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10.1, 3, 1.5) + 0.1);            # The default is ‘c(5, 4, 4, 2) + 0.1’.
# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = moduleTraitCor,
                         xLabels = names(traits),
                         yLabels = names(MEs),
                         ySymbols = names(MEs),
                         colorLabels = FALSE,
                         colors = blueWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 0.7,
                         textAdj = c(0.5, 0.5),
                         zlim = c(-1,1),
                         maxColsPerPage = 15,
                         maxRowsPerPage = 25,
                         main = paste("Module-trait relationships"))
dev.off()

# FGF21 only
pdf(file = "mod-trait.relations_fgf21_format.pdf", width = 5, height = 7)
textMatrix <- paste(signif(moduleTraitCor[,4], 2), "\n(",
                    signif(moduleTraitPvalue[,4], 1), ")", sep = "")
dim(textMatrix) <- dim(data.frame(fgf21.serum = moduleTraitCor[,4]))
par(mar = c(5, 10.1, 3, 1.5) + 0.1);            # The default is ‘c(5, 4, 4, 2) + 0.1’.
# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = data.frame(fgf21.serum = moduleTraitCor[,4]),
                         xLabels = names(traits)[4],
                         yLabels = names(MEs),
                         ySymbols = names(MEs),
                         colorLabels = FALSE,
                         colors = blueWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 1,
                         pointsize = 16,
                         textAdj = c(0.5, 0.5),
                         zlim = c(-1,1),
                         maxColsPerPage = 15,
                         maxRowsPerPage = 25,
                         main = paste("Module-trait relationships"))
dev.off()

# all
pdf(file = "mod-trait.relationships.pdf", width = 12, height = 14)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10.1, 3, 1.5) + 0.1);            # The default is ‘c(5, 4, 4, 2) + 0.1’.
# Display the correlation values within a heatmap plot
labeledHeatmap(moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               textAdj = c(0.5, 0.5),
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignificance
#
#=====================================================================================

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, traits))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


geneTraitSignificance = as.data.frame(cor(datExpr, traits, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(traits), sep="");
names(GSPvalue) = paste("p.GS.", names(traits), sep="");

# Plot the dendrogram
pdf(file = "6-Eigengene_dendrogram.pdf", width = 12, height = 7)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf(file = "6-Eigengene adjacency heatmap.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

#======================================
## fgf21
fgf21 = as.data.frame(traits$fgf21.serum)
# rename colname
names(fgf21) = "fgf21.serum"
# Add the weight to existing module eigengenes 
MEsfgf21 <- orderMEs(cbind(MEs,fgf21))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "6-eigengene.relationships.fgf21.pdf", wi = 6, he = 6) 
plotEigengeneNetworks(MEsfgf21,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()

modules = "darkseagreen4" #up: #maroon #salmon4 down: #darkseagreen4
# Rename to moduleColors
moduleColors = mergedColors
column = match(modules, modNames);
moduleGenes = moduleColors==modules;

pdf(file = "7-Module membership vs gene significance darkseagreen4.pdf", width = 7, height = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 4]),
                   lmFnc = lm,
                   #abline = T,
                   abline.color = "red",
                   xlab = paste("Module Membership in", modules, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")
dev.off()

#=====================================================================================
#
#  Extract gene list from module
#
#=====================================================================================

names(merge$newMEs)
modules = c(substring(names(merge$newMEs)[2], 3)); #maroon 5, paleturquoise 12, darkseagreen4 2
genes = colnames(datExpr)
inModule = is.finite(match(mergedColors,modules))
modGenes = data.frame(ensembl.id = genes[inModule],
                      Symbol = mapIds(org.Mm.eg.db, genes[inModule],
                             keytype = "ENSEMBL", column = "SYMBOL"),
                      moduleColor = mergedColors[inModule],
                      geneTraitSignificance[genes[inModule],],
                      GSPvalue[genes[inModule],])

# Order modules by their significance for fgf21
modOrder = order(-abs(cor(MEs, traits$fgf21.serum, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(modGenes)
  modGenes = data.frame(modGenes, geneModuleMembership[genes[inModule], modOrder[mod]], 
                         MMPvalue[genes[inModule], modOrder[mod]]);
  names(modGenes) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the modGenes variable first by module color, then by geneTraitSignificance
geneOrder = order(modGenes$moduleColor, -abs(modGenes$GS.fgf21.serum));
geneInfo = modGenes[geneOrder, ]

write.csv(geneInfo, paste("geneInfo_from_",paste(modules,collapse = "_"), "_merged.csv", sep = ""))