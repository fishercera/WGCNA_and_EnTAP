### 2018 C R Fisher
### First module of WGCNA process
## This part does not need parallelization.

# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# Load the WGCNA package
getwd()
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#allowWGCNAThreads(nThreads = 16)

#setwd("D:/Cera Fisher/Google Drive/")
#Step 0 read the config file and get variable names
Config <- read.delim("WGCNA_Config.txt", header=FALSE, sep="\t")
infile <- Config$V1
traitfile <- Config$V2
basename <- Config$V3
annots <- Config$V4

lnames = load(file = paste(basename, "WGCNA.RData", sep="_"));
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = paste(basename, "modules.RData", sep="_"))
lnames

####Step 1: Setting up module eigengenes with color-names

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

library(data.table)
### Make a key/value data.table of each trait and the highest correlation module name
dict <- data.frame(x=colnames(moduleTraitCor), y=1:length(colnames(moduleTraitCor)))

for (i in seq(1,length(colnames(moduleTraitCor)),1)) {
  moduleName <- names(which(moduleTraitCor[,i]==max(moduleTraitCor[,i])))
  dict[i,] <- c(colnames(moduleTraitCor)[i],moduleName)
}
dict <- data.table(dict)
write.table(dict, paste(basename, "Traits_and_HighestCorModule.tab", sep="_"), sep="\t")

dict2 <- data.frame(x=rownames(moduleTraitCor), y=1:length(rownames(moduleTraitCor)))

for (i in seq(1,length(rownames(moduleTraitCor)),1)) {
  print(rownames(moduleTraitCor)[i])
  traitName <- names(which(moduleTraitCor[i,]==max(moduleTraitCor[i,])))
  dict2[i,] <- c(rownames(moduleTraitCor)[i],traitName)
}
dict2 <- data.table(dict2)
write.table(dict2, paste(basename, "Modules_and_HighestCorTrait.tab", sep="_"), sep="\t")
### Written to a file. 

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))

# Display the correlation values within a heatmap plot
pdf("Eigengene_Trait_Correlation_heatmap.pdf", width=10, height=8)
par(mar = c(6, 13, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()



for (i in seq(1,length(colnames(datTraits)),1)) {
  colnum <- which(colnames(datTraits)==as.character(dict[i,1]))
  trait = as.data.frame(datTraits[,colnum])
  names(trait) = colnames(datTraits)[colnum]
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(cor(datExpr0, trait, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(trait), sep="")
  names(GSPvalue) = paste("p.GS.", names(trait), sep="")
  g <- moduleTraitCor[,which(colnames(moduleTraitCor)==colnames(datTraits)[colnum])]
  module = modNames[which(g==max(g))]
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  numGenes <- as.numeric(table(moduleGenes)[2])
  par(mfrow = c(1,1))
  pdf(paste(basename, module, names(trait), "Scatterplot.pdf", sep="_"))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module\n#Genes: ", numGenes, sep=" "),
                     ylab = paste("Gene significance for trait:", names(trait), sep=" "),
                     main = paste("Module membership vs. gene significance\n "),
                     displayAsZero = 0,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
  dev.off()
}


####

tissue = as.data.frame(datTraits$Pro)
names(tissue) = "Pro" 
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr0, tissue, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(tissue), sep="")
names(GSPvalue) = paste("p.GS.", names(tissue), sep="")

g <- moduleTraitCor[,which(colnames(moduleTraitCor)=="Pro")]
module = modNames[which(g==max(g))]




column = match(module, modNames)
moduleGenes = moduleColors==module
table(moduleGenes)

par(mfrow = c(1,1))
pdf(paste(basename, module, names(tissue), "Scatterplot.pdf", sep="_"))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module", sep=" "),
                   ylab = paste("Gene significance for trait:", names(tissue), sep=" "),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
dev.off()



save(MEs,moduleTraitCor,dict,dict2,moduleTraitPvalue,geneModuleMembership,MMPvalue,GSPvalue, geneTraitSignificance,
     file = paste(basename, "module3.RData", sep="_"))
