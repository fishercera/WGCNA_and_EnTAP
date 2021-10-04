### 2018 C R Fisher
### Fifth module of WGCNA process
## This part does need parallelization 
# Main purpose of this module - visualize Topological Overlap Map for a selection of genes
library(WGCNA)
#setwd("D:/Cera Fisher/Google Drive/Treehoppers/ResearchFiles/SRAProject/Hemiptera/Hemiptera2/")
Config <- read.table("WGCNA_Config.txt", header=FALSE, sep="\t")

infile <- Config$V1
traitfile <- Config$V2
basename <- Config$V3
annotations <- Config$V4

options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
allowWGCNAThreads(nThreads = 16)

lnames = load(file = paste(basename, "WGCNA.RData", sep="_"));
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = paste(basename, "modules.RData", sep="_"));
lnames

lnames = load(file = paste(basename, "-block.1.RData", sep=""))
lnames 

datExpr0 <- datExpr0
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

#=====================================================================================
#
#  Code chunk 1 - Plot heatmap for selection of module genes
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
### If you need to recalculate TOM, uncomment the next line and comment out the subsequent. 
#dissTOM = 1-TOMsimilarityFromExpr(datExpr0, power = 6);
### If you saved TOM, you can do this: 
dissTOM <- as.matrix(1-TOM)


nSelect = 2000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
# This power is something you might want to set as a variable
plotDiss = selectTOM^12
diag(plotDiss) = NA

# Call the plot function
png(paste(basename, as.character(nSelect), "selectedGenes_TOMplot.png", sep="_"), width=1000, height=1000)
#sizeGrWindow(9,9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()
# 
# plotTOM <- dissTOM
# diag(plotTOM) = NA
# png(paste(basename, "allGenes_TOMplot.png", sep="_"), width=1000, height=1000)
# par(cex=2)
# TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
# dev.off()

##########

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr0, moduleColors)$eigengenes
# Isolate weight from the clinical traits
trait = as.data.frame(datTraits$Wings);

# To be investigated - what was this about?
  
pdf(paste(basename, "allTraits_allMods", "dendro.pdf", sep="_"), width=10, height=10)
METall <- orderMEs(cbind(MEs,datTraits))
plotEigengeneNetworks(METall, "", marDendro = c(0,8,1,2), 
                      marHeatmap = c(5,8,1,2),
                      cex.lab = 0.8, 
                      xLabelsAngle = 90)

dev.off()
names(trait) = "Wings"
# Add the trait to existing module eigengenes
MET = orderMEs(cbind(MEs, trait))
# Plot the relationships among the eigengenes and the trait
pdf(file=paste(basename, "EigengeneNetworks_dendro", names(trait), sep="_"))

par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2),
                      cex.lab = 0.8, 
                      xLabelsAngle = 90)
dev.off()


#save(dissTOM, plotTOM, file=paste(basename, "TOMplot.RData", sep="_"))
