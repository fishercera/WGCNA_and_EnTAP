### 2018 C R Fisher
### First module of WGCNA process
## This part does not need parallelization. 

getwd()
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "."
setwd(workingDir) 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)
disableWGCNAThreads()

#Step 0 read the config file and get variable names
Config <- read.delim("WGCNA_Config.txt", header=FALSE, sep="\t")
infile <- Config$V1
traitfile <- Config$V2
basename <- Config$V3
annots <- Config$V4

#Step1 read in the data
Data = read.delim(infile, sep="\t", header=TRUE)
##take a quick look at the data set:
dim(Data)

colnames(Data)
rownames(Data)
#Step2 rearrange and manipulate the data -- only want columns with counts
datExpr0 = as.data.frame(t(Data[,-c(1:3)]))
storage.mode(datExpr0[77,]) <- "integer"

colnames(datExpr0) = Data$ECid


####Step 3: Run builtin data conformity check. 

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

## This part is optional -- in case the function checking to make sure the data is good returns anything bad

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


####Step 4: cluster samples to check for outliers



sampleTree = hclust(dist(datExpr0), method = "ward.D")

png(paste(basename, "Scaled_TPM_SampleClustering_toDetectOutliers.png", sep="_"))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

##HV7 samples are WAY different from everything else. Going to drop them. 
rownames(datExpr0)
#datExpr0 <- datExpr0[c(1:16,22:36),] #Dropping five samples from HV7
rownames(datExpr0)
####Step 5: Read in trait data for samples. 


getwd()
traitfile <- "ECEF-withpronwings-expanded_samplesntraits.txt"
traitData = read.delim(traitfile, sep="\t", header=TRUE)
dim(traitData)
colnames(traitData)
traitData$Sample


# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

Samples = rownames(datExpr0)
traitRows = match(Samples, allTraits$Sample)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()


####Step 6: Re-cluster samples

# We're also going to add a simple heatmap below the dendrogram, showing how
# each of the traits we specified maps against the dendrogram.


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "ward.D")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
png(paste(basename, "sampleDendro_traitHeatmap.png", sep="_"))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
#Save workspace to load in next module
save(datExpr0, datTraits, file = paste(basename, "WGCNA.RData", sep="_"))

