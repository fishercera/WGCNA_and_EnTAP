### 2018 C R Fisher
### First module of WGCNA process
## This part does not need parallelization. 

# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# Load the WGCNA package
getwd()
setwd(".")
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)
Config <- read.delim("WGCNA_Config.txt", header=FALSE, sep="\t")
infile <- Config$V1[2]
traitfile <- Config$V2[2]
basename <- Config$V3[2]
annots <- Config$V4[2]

#Load the data from the first part
lnames = load(file = paste(basename, "WGCNA.RData", sep="_"))
#The variable lnames contains the names of loaded variables.
lnames

####Step 1: Choose soft-thresholding powers.
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
# Scale-free topology fit index as a function of the soft-thresholding power
par(mfrow = c(1,2));
cex1 = 0.9;

pdf(paste(basename, "SoftThreshold_plots.pdf", sep="_"))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) 
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

print("SFT values")
sft

####Step 2: Construct network

## This part should probably be done AFTER the soft-threshold level is selected. 

# net = blockwiseModules(datExpr0, power = sft$powerEstimate,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = basename, 
#                        maxBlockSize = 20000, 
#                        verbose = 3)

TOM <- adjacency(datExpr0, power=6)
TOM <- TOMsimilarity(adj)
dissTOM <- 1-TOM
rm(TOM)
geneTree <- hclust(as.dist(dissTOM), method="average")
sizeGrWindow(12,9)
png(paste(basename, "TOM_geneTree.png", sep="_"))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

minModuleSize <- 15
dynamicMods <- cutreeDynamic(dendro=geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)

table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

##
save(dynamicMods, dynamicColors, TOM, geneTree, file=paste(basename, "unmerged_and_TOM.RData", sep="_"))
##
sizeGrWindow(8,6)
png(paste(basename, "ModuleColors_and_GeneTree.png", sep="_"))
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
png(paste(basename, "ME_tree", "WithCutHeight.png", sep="_"))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
dev.off()
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

png(paste(basename, "Merged_Modules.png", sep="_"))
#Convert labels to colors for plotting

#Plot the dendrogram and the module colors underneath
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                   c("Dynamic Tree Cut", "Merged dynamic"),
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


## Save modules and TOM to file
lnames = load(file = paste("Ofas_WGCNA_WGCNA.RData", sep="_"))
lnames = load(file = "Ofas_WGCNA_modules.RData")

save(MEs, moduleLabels, moduleColors, geneTree, file = paste(basename, "modules.RData", sep="_"))
save(TOM, dissTOM, file=paste(basename, "-block.1.RData", sep=""))


