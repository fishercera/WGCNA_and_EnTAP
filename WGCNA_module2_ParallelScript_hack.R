### 2018 C R Fisher
### First module of WGCNA process
## This part does need parallelization, the more cpus the better

# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# Load the WGCNA package
getwd()
setwd(".")
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 16)                                                          # Edit this to add functionality for multi-threading
Config <- read.delim("WGCNA_Config.txt", header=FALSE, sep="\t")                         # Read config file and set variables
infile <- Config$V1                                                                     
traitfile <- Config$V2
basename <- Config$V3
annots <- Config$V4

                                                                                          #Load the data from the first part
lnames = load(file = paste(basename, "WGCNA.RData", sep="_"))
                                                                                          #The variable lnames contains the names of loaded variables.
lnames

####Step 1: Choose soft-thresholding powers.
powers = c(c(1:10), seq(from = 12, to=40, by=2))                                         # set a range of powers to try for the soft threshold 
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)                     # Call the network topology analysis function 
                                                                                         # Plot the results:
               # Scale-free topology fit index as a function of the soft-thresholding power
par(mfrow = c(1,2));
cex1 = 0.9;

pdf(paste(basename, "SoftThreshold_plots.pdf", sep="_"))                                  # Print the plot to a PDF
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) 
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.80,col="red")                                                                  # this red-line line corresponds to using an R^2 cut-off of h
                # Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

print("SFT values")                                                                       # Print SFT results to stdout
sft

####Step 2: Construct network

                                                                                          ## This part must be done AFTER the soft-threshold level is selected. 

# net = blockwiseModules(datExpr0, power = sft$powerEstimate,                             ## making the modules blockwise is a little more computer-friendly, 
#                        TOMType = "unsigned", minModuleSize = 30,                        # but seems to cost accuracy, I don't recommend it. 
#                        reassignThreshold = 0, mergeCutHeight = 0.25,                    # use the code below instead to generate the full TOM (topological overlap map)
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,                 # and the dissimilarity matrix for clustering. 
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = basename, 
#                        maxBlockSize = 20000, 
#                        verbose = 3)

#                                               This is the really computationally intensive stuff. For best results, run this on a beefy machine (>12 cores)
#adj <- adjacency(datExpr0, power=6)                                                       # Compute adjacency
#TOM <- TOMsimilarity(adj)                                                                 # Compute similarity matrix

lnames = load(file = "NewTestRunTOM.RData")
str(dissTOM)

#dissTOM <- 1-TOM                                                                          # derive dissimilarity matrix (AKA a distance matrix)
#rm(TOM)                                                                                  # TOM and dissTOM are huge, and the same size. To free up RAM, delete TOM. 
geneTree <- hclust(as.dist(dissTOM), method="average")                                    # Making a gene tree with hclust
sizeGrWindow(12,9)                                                                        # Let's print out the gene tree/hclust dendrogram!
png(paste(basename, "TOM_geneTree.png", sep="_"))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()
## Set this as a variable up top so we can play around with it
minModuleSize <- 15                                                                       # Set a minimum module size. 
                                                                                          # Now we want to cut the tree dynamically to generate modules of similar genes
dynamicMods <- cutreeDynamic(dendro=geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)                              # These are default options, I believe.

table(dynamicMods)                                                                        # dynamicMods holds the initial groups of genes that seem to be most similar

dynamicColors = labels2colors(dynamicMods)                                                # Let's give them color names.
table(dynamicColors)                                                                      # how many genes in each color-module? 

TOM <- 1-dissTOM
##
save(dynamicMods, dynamicColors, TOM, geneTree, file=paste(basename, "unmerged_and_TOM.RData", sep="_")) # We want to use these things again later. 
##
sizeGrWindow(8,6)
png(paste(basename, "ModuleColors_and_GeneTree.png", sep="_"))                            # Print out the gene-tree with the modules each gene belongs to. 
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
### Step XX: Make Module eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)                               # Compute module eigengenes -- MEs are generalizations 
MEs = MEList$eigengenes                                                                   #   of the expression profile of all the genes in a module.
MEDiss = 1-cor(MEs);                                                                      # Let's make a distance matrix of the MEs' correlation with each other
METree = hclust(as.dist(MEDiss), method = "average");                                     # Make a dendrogram of the module eigengenes (to see if we can merge similar modules)
png(paste(basename, "ME_tree", "WithCutHeight.png", sep="_"))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
## Set this as a variable up top so we can play around with it 
MEDissThres = 0.25                                                                        # We're going to merge genes using a dissimilarity threshold of 0.25 (tree heigh)
abline(h=MEDissThres, col = "red")                                                        # This cut-line shows which modules are going to get merged. 
dev.off()
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)  # Now let's merge the close modules. 
mergedColors = merge$colors;                                                              # Get the remaining colors
mergedMEs = merge$newMEs;                                                                 # And get the remaining MEs

png(paste(basename, "Merged_Modules.png", sep="_"))                                       # Let's show what these merged modules look like 
#Convert labels to colors for plotting

#Plot the dendrogram and the module colors underneath
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),                         # Plotting the dendrogram again with our gene tree and module colors
                   c("Dynamic Tree Cut", "Merged dynamic"),
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


## Save modules and TOM to file
lnames = load(file = paste(basename, "WGCNA.RData", sep="_"))                           # Let's make sure we save this data! 
#lnames = load(file = paste(basename, "modules.RData", sep="_"))

save(MEs, moduleLabels, moduleColors, geneTree, file = paste(basename, "modules.RData", sep="_"))
save(TOM, dissTOM, file=paste(basename, "TOM.RData", sep=""))


