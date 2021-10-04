### 2018 C R Fisher
### Fourth module of WGCNA process
## This part does not need parallelization. 
getwd()
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "."
setwd(workingDir) 
# Load the WGCNA package
library(WGCNA)
options(stringsAsFactors = FALSE)
#allowWGCNAThreads(nThreads = 16)

print("Read the config file, pretty please")
Config <- read.table("WGCNA_Config.txt", header=FALSE, sep="\t")
print("I read the table")
infile <- Config$V1
traitfile <- Config$V2
basename <- Config$V3
annotations <- Config$V4

print("These are the variables")
print(paste(infile, traitfile, basename, annotations, sep=" "))
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the expression and trait data saved in the first part
lnames = load(file = paste(basename, "WGCNA.RData", sep="_"));
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = paste(basename, "modules.RData", sep="_"));
lnames

datExpr0 <- data.frame(datExpr0)
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEsList = moduleEigengenes(datExpr0, moduleColors)$eigengenes
## Solved the bug with the above code - needed to be pulling only the $eigengenes part of moduleEigengenes
MEs = orderMEs(MEsList)
### Debug here -- it appears to be a problem that can be solved with an as.numeric coercion
### https://www.programmingr.com/r-error-messages/list-object-cannot-be-coerced-to-type-double/
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

pdf(paste(basename, "mod-trait-relationships.pdf", sep="_"), width=20, height=10)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", formatC(moduleTraitPvalue, format="e", digits=2), ")",
                    sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 15, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# This is also repetition from the last module - doing verboseScatterplots as one-offs
# Define variable weight containing the weight column of datTrait
tissue = as.data.frame(datTraits$Pro);
names(tissue) = "Pro"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
g <- moduleTraitCor[,which(colnames(moduleTraitCor)=="Pro")]

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr0, tissue, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(tissue), sep="");
names(GSPvalue) = paste("p.GS.", names(tissue), sep="");


#=====================================================================================
#
#  Code chunk 5 -- makes a scatter plot, this isn't necessary?
#
#=====================================================================================

#Get the module that has the most significance for the variable we're looking at
#module = moduleColors[which(GSPvalue==min(GSPvalue))]
module = modNames[which(g==max(g))] # this gets the module with the highest correlation for the trait we're interested in 
#module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#build a scatterplot for it
pdf(paste(basename, names(tissue), module, "scatterplot_significant-module.pdf", sep="_"))
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module", sep=" "),
                   ylab = paste("Gene significance for trait:", names(tissue), sep=" "),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()

#=====================================================================================
#
#  Code chunk 7 - get all the gene ids for a particular module and write them to a file
#
#=====================================================================================
# This can be modified for any desired module
## This may be something you would want to define at the top
module="plum2"
MemberGenes <- names(datExpr0)[moduleColors==module]
write.csv(MemberGenes, paste(basename, module, "memberGenes.csv", sep="_"))


#=====================================================================================
#
#  Code chunk 8 - attach annotations to ids
#
#=====================================================================================
# These annotations come from Entap. 
# This is bespoke code. 
#### Need to get Annotations for my OrthoGroups!!
print("Now please read the annotations file please")
print(annotations)
annot = read.delim(annotations, sep="\t")
dim(annot)
names(annot)
names(annot)[1] <- "ECid"
head(names(datExpr0))
head(annot$ECid)
#annot$HVid <- gsub("-R", "-P", annot$HVid)
probes = names(datExpr0)

probes2annot = match(probes, annot$ECid)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0. 

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

names(annot)
# Create the starting data frame
names(annot) <- c("ECid", 
                  "AccessionBlastHit", 
                  "BLASTp.Hit", 
                  "taxon", 
                  "OrthoSeed", 
                  "GeneSymbol", 
                  "EggNOG.Orthogroups", 
                  "Descrip1", 
                  "Descrip2", 
                  "ProteinDomains", 
                  "GO.BP", 
                  "GO.CC", 
                  "GO.MF", 
                  "IPRO.fam", 
                  "Pfam", 
                  "")
geneInfo0 = data.frame(ECid = probes,
                       geneName = annot$BLASTp.Hit[probes2annot],
                       Descrip1 = annot$Descrip1[probes2annot],
                       Descrip2 = annot$Descrip2[probes2annot],
                       ConservedDomains = annot$ProteinDomains[probes2annot], 
                       OrthoSeed = annot$OrthoSeed[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for our given trait of interest
# We had set that trait earlier as "tissue"
modOrder = order(-abs(cor(MEs, tissue, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Pro));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

save(geneInfo, moduleTraitCor, moduleTraitPvalue, MEsList, file=paste(basename, "geneInfo.RData", sep="_"))
write.csv(geneInfo, file = paste(basename, "geneInfo_relevanceForPro.csv", sep="_"))

#save(geneInfo0, geneModuleMembership, geneTraitSignificance, file="Hemiptera_Pro_GeneInfo_WingsTrait.RData")

library(dplyr)

module <- "plum2"

moduleGeneInfo <- filter(geneInfo, moduleColor==module)
write.table(moduleGeneInfo, file=paste(basename, module, names(tissue), "geneInfo.tab", sep="_"), sep="\t")

## This is in development - 
## select hub genes for a given module
## Figure out what we wanted to do with this? 
moduleInfo <- moduleGeneInfo %>% arrange(MM.darkmagenta) %>%
  select(HVid, geneName, Descrip1, Descrip2, MM.darkmagenta, GS.Wings)
hubGenes <- moduleInfo %>% top_n(20, MM.darkmagenta)
