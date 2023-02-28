#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

##data input
##input expression level (e.g. TPM) of each gene (row)

# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored
# "." means current directory. On Windows use a forward slash / instead of the usual \
workingDir = "."
setwd(workingDir)
# Load the WGCNA package
library(tidyverse)
library(WGCNA)

# The following setting is important, do not omit
options(stringsAsFactors = FALSE)
#Read in the female liver data set
InputData = read.csv("IN_TPM.csv")
# Take a quick look at what is in the data set
dim(InputData)
names(InputData)

# transfer dataset
# only TPM input
datExpr0 = as.data.frame(t(InputData[,c(4:length(InputData[1,]))]))
# rename the col name of new dataset based on the col "id" of the original dataset
names(datExpr0) = InputData$ko_id
# give the original treatment name back to the new dataset
rownames(datExpr0) = names(InputData)[c(4:length(InputData[1,]))]
datExpr0$sra <- row.names(datExpr0) 

# input trait data (descriptions or environmental factors of each sample)
traitData = read.csv("IN_Taraocean_Trait.csv")
# filtering data
traitData <- traitData %>%
  filter(temperature >= 15, temperature <= 31, depth <= 5, Chlorophyll > 0.3)

datExpr_filted_full <- inner_join(traitData, datExpr0, by="sra")
datExpr_filted_sra <- datExpr_filted_full %>%
  select(sra)

#only keep selected sra
datExpr0 <- inner_join(datExpr_filted_sra, datExpr0, by="sra")
row.names(datExpr0) <- datExpr0$sra
datExpr0 <- datExpr0 %>%
  select(-sra)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

##check if all genes and samples are good

gsg = goodSamplesGenes(datExpr0, verbose = 3)
# if TRUE then all genes are OK
gsg$allOK
# check summary
summary(gsg)
# remove bad genes if not all OK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

##check outlier samples

sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

##remove samples which higher than "CutH"
CutH = 300000

# Plot a line to show the cut
abline(h = CutH, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = CutH, minSize = 10)
# clust 1 contains the samples we want to keep
table(clust)
keepSamples = (clust==1)
# 'datExpr' will be used in further steps
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# make a clean traits
names(traitData)
# select columns that hold information we need
# only keep numeric data
allTraits <- traitData %>%
  select(nitrate, Chlorophyll, oxy, Salinity, temperature)

#set all negative value as zero
i = 1
while(i <= length(allTraits[,1])){
  j = 1
  while(j <= length(allTraits[1,])){
    if(allTraits[i,j] < 0){
      allTraits[i,j] <- 0
    }
    j = j+1
  }
  print(i/length(allTraits[,1]))
  i= i+1
}

#set all extremely high value as na
i = 1
while(i <= length(allTraits[,1])){
  j = 1
  while(j <= length(allTraits[1,])){
  if(allTraits[i,j] > 9999){
    allTraits[i,j] <- NA
  }
  j = j+1
  }
  print(i/length(allTraits[,1]))
  i= i+1
}

allTraits$sra <- traitData$sra
  
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the traits
datTraits <- inner_join(datExpr_filted_sra, allTraits, by="sra")
rownames(datTraits) <- datTraits$sra
datTraits <- datTraits %>%
  select(-sra)
collectGarbage()

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

##have a visual check of all samples and their traits

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

##output current data, will be used in next step

save(datExpr, datTraits, file = "01-dataInput.RData")


