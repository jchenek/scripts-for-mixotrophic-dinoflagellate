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
library(tidyverse)
library(WGCNA)
library(ggplot2)
# The following setting is important, do not omit
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "01-dataInput.RData")
#The variable lnames contains the names of loaded variables
lnames
# Load network data saved in the second part.
lnames = load(file = "02-networkConstruction-stepByStep.RData")
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

##define module-trait relationships

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# plot module-trait relationships
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
mycolors = colorRampPalette(c("#6699FF","white","#FF3300"))(20)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = mycolors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# pdf(file = "wgcna_plot_Mod-tra_relation.pdf", wi = 8, he = 6)
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = mycolors,
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 1,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# dev.off()

# output Module-trait relationships
Gene_module_info = as.data.frame(moduleColors)
rownames(Gene_module_info) = colnames(datExpr)
write.csv(Gene_module_info,file="wgcna-output-module-trait-relationships.csv")
write.csv(moduleTraitCor,file="wgcna-output-moduleTraitCor.csv")
write.csv(textMatrix,file="wgcna-output-textMatrix.csv")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

##check geneModuleMembership and geneTraitSignificance of defined module and trait

# Define Trait
names(datTraits)
targetTrait = as.data.frame(datTraits$temperature)

# names (colors) of the modules
modNames = substring(names(MEs), 3)

# get geneModuleMembership and geneTraitSignificance
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, targetTrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(targetTrait), sep="")
names(GSPvalue) = paste("p.GS.", names(targetTrait), sep="")

# Define a module
substring(names(MEs), 3)
module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors==module

# quick check Module membership vs. gene significance
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# set thresholds and choose hub genes
# make dataset
hub_data = as.data.frame(abs(geneModuleMembership[moduleGenes, column]))
hub_data$'abs(geneTraitSignificance[moduleGenes, 1])' = abs(geneTraitSignificance[moduleGenes, 1])
rownames(hub_data) = rownames(geneModuleMembership[which(moduleGenes== TRUE),])
names(hub_data) = c("Membership","Significance")

#select high value genes
T_geneModuleMembership = 0.8
T_geneTraitSignificance = 0.2
hub_data <- hub_data %>%
  mutate(Color = case_when(Membership > T_geneModuleMembership & Significance > T_geneTraitSignificance ~ "A",
                           Membership <= T_geneModuleMembership | Significance <= T_geneTraitSignificance ~ "B"))

bub_plot <- ggplot(data = hub_data, mapping = aes(x = Membership, y = Significance, fill = Color)) +
  scale_color_manual(values=c("#3399FF", "#FF6600"))+
  geom_point(shape = 21, colour = "grey", stroke = 0.5,  size = 4.5, alpha = 4/10)+
  geom_hline(yintercept = T_geneTraitSignificance, linetype = "dashed", color = "#33333380", size = 0.6)+
  geom_vline(xintercept = T_geneModuleMembership, linetype = "dashed", color = "#33333380", size = 0.6) +
  theme_light()+
  theme(aspect.ratio = 1/1)
bub_plot
ggsave("wgcna_plot_gene-trait-sig-memb.pdf", plot = bub_plot, width = 5, height = 7, dpi=600)

# genes with color "A" are hub genes 
write.csv(hub_data,file="wgcna-output-Membership-Significance.csv")

