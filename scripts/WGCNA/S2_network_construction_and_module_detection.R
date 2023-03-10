#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

##input data proceed in last step
##will allow multi-threading within WGCNA

# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored
# "." means current directory. On Windows use a forward slash / instead of the usual \
workingDir = "."
setwd(workingDir)
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. At present this call is necessary
# Any error here may be ignored but you may want to update WGCNA if you see one
# Caution: skip this line if you run RStudio or other third-party R environments
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "01-dataInput.RData")
#The variable lnames contains the names of loaded variables
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

##choose the first number meet the 'abline' as softPower, abline default = 0.9 but 0.8 also accepted

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# set the softpower
softPower = 12
adjacency = adjacency(datExpr, power = softPower, type = "signed")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

##calculate TOM (Topological overlap matrix)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

##set the minimum module size, this will merge small module to a larger module

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize)
table(dynamicMods)

##visualize clusters and modules

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

##merge similar module into larger module

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# set threshold to merge similar modules
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plot final module profiles
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

pdf(file = "wgcna_plot_geneDendro.pdf", wi = 9, he = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

##output eigengenes and related modules

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "02-networkConstruction-stepByStep.RData")
