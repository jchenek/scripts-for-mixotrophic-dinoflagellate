#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored
# "." means current directory. On Windows use a forward slash / instead of the usual \
workingDir = "."
setwd(workingDir)
# Load the WGCNA package
library(tidyverse)
library(igraph)
library(WGCNA)
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

# Recalculate topological overlap if needed
# set power according to S2
TOM = TOMsimilarityFromExpr(datExpr, power = 12, networkType = "signed")

# Select modules
substring(names(MEs), 3)
modules = c("brown")

# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# set threshold and remove low value connection
ThreadHold = 0.25
filted_modTOM = modTOM
filted_modTOM[filted_modTOM < ThreadHold] = 0
filted_modTOM <- filted_modTOM[which(rowSums(filted_modTOM)>1),]
filted_modTOM <- filted_modTOM[,which(colSums(filted_modTOM)>1)] 

# # read "wgcna-output-Membership-Significance.csv" and only keep hub genes
# hub_data = read.csv("wgcna-output-Membership-Significance.csv", row.names = 1, header = TRUE)
# filted_modTOM <- filted_modTOM[rownames(hub_data[which(hub_data$Color== "A"),]),]
# filted_modTOM <- filted_modTOM[,rownames(hub_data[which(hub_data$Color== "A"),])] 

# quick check the network with igraph
igraphTOM = graph_from_adjacency_matrix(filted_modTOM, mode = c("undirected"), weighted = TRUE, diag = FALSE)
set.seed(520)
plot(igraphTOM, main="quick check network",
     vertex.color="orange", vertex.frame.color="grey", vertex.size=5, vertex.label=NA, 
     edge.color="grey", edge.width=1, edge.lty=1, edge.curved=TRUE, margin=c(0,0,0,0))

# draw nicer graph
nicer_igraphTOM = igraphTOM
# give degree to nodes
#V(nicer_igraphTOM)$size=igraph::degree(nicer_igraphTOM)%>%log()*2
V(nicer_igraphTOM)$size=7
# give width to edges
#E(nicer_igraphTOM)$width=E(nicer_igraphTOM)$weight%>%log()/2
E(nicer_igraphTOM)$width=E(nicer_igraphTOM)$weight%>%log()/10
# check graph
set.seed(520)
plot(nicer_igraphTOM, main="network",
     vertex.color="orange", vertex.frame.color="grey", vertex.label=NA, 
     edge.color="grey", edge.lty=1, edge.curved=TRUE, margin=c(0,0,0,0))#,layout=layout_with_drl

# color hub gene
hub_data <- read.csv("wgcna-output-Membership-Significance.csv", row.names = 1, header = TRUE)
hub_data <- hub_data %>%
  mutate(gene_id = row.names(hub_data))

# add annotation
gene_list <-  read.csv('IN_target_gene_list.csv', header = FALSE, row.names = 1)
names(gene_list) <- c("name","desc")
gene_list <- gene_list %>%
  mutate(gene_id = row.names(gene_list))

anno_igraphTOM = as.data.frame(V(igraphTOM)$name)
names(anno_igraphTOM) <- "gene_id"
anno_igraphTOM <- left_join(anno_igraphTOM, hub_data,by = "gene_id")
anno_igraphTOM <- left_join(anno_igraphTOM, gene_list,by = "gene_id")
anno_igraphTOM <- anno_igraphTOM %>%
  mutate(c1 = name) %>%
  mutate(c2 = case_when(is.na(name) == FALSE ~ "#FF33FF",
                        Color == "A" ~ paste(adjustcolor("#FF3300", alpha.f = 0.38)),
                        Color == "B" ~ paste(adjustcolor("#00CCCC", alpha.f = 0.55)),
  ))

# annotate each node using excel
# write.csv(anno_igraphTOM,file="anno_igraphTOM.csv")
# edit target gene label and color in excel, and read them back
# e.g.:
# 	gene_id	c1	c2
# 1	gene1	NA	orange
# 2	gene2	gene_id	red
# 3 gene3 NA orange
# 4 gene4 NA orange
# or
# set color with alpha
# the higher the darker
# check the value of "adjustcolor("orange", alpha.f = 0.5)", which is #FFA5004D, and edit excel
# adjustcolor("#3399FF", alpha.f = 0.3) #3399FF4D
# anno_igraphTOM = read.csv("anno_igraphTOM.csv",row.names = NULL)

# give label to nodes
V(nicer_igraphTOM)$name=anno_igraphTOM$c1
# give color to nodes
V(nicer_igraphTOM)$color=anno_igraphTOM$c2
set.seed(520)
plot(nicer_igraphTOM, main="network",
     vertex.label.color="black", vertex.label.font=1, vertex.label.cex=1,
     vertex.frame.color="grey",
     edge.color="grey", edge.lty=1, edge.curved=TRUE, margin=c(0,0,0,0))#,layout=layout_with_drl

#plot one with label
pdf(file = "wgcna_plot_network.pdf", wi = 9, he = 6)
set.seed(520)
plot(nicer_igraphTOM, main="network",
     vertex.label.color="black", vertex.label.font=1, vertex.label.cex=1,
     vertex.frame.color="grey",
     edge.color=adjustcolor("grey", alpha.f = 0.38), edge.lty=1, edge.curved=TRUE, margin=c(0,0,0,0))
dev.off()

#plot on without label
anno_igraphTOM$c1 = NA
V(nicer_igraphTOM)$name=anno_igraphTOM$c1
pdf(file = "wgcna_plot_nolabel_network.pdf", wi = 9, he = 6)
set.seed(520)
plot(nicer_igraphTOM, main="network",
     vertex.label.color="black", vertex.label.font=1, vertex.label.cex=1,
     vertex.frame.color="grey",
     edge.color=adjustcolor("grey", alpha.f = 0.38), edge.lty=1, edge.curved=TRUE, margin=c(0,0,0,0))
dev.off()

# Export filted TOM matrix in csv format  
write.csv(filted_modTOM,file="wgcna-output-correlation.csv")
