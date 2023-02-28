####################################################################################################################
#reference: edgeR: a Bioconductor package for differential expression analysis of digital gene expression data
library(tidyverse)
library(edgeR)

#input unified data
raw_readcount_data <-read.csv('raw_counts.csv', header=T, check.names = FALSE, row.names = 1)

####################################################################################################################
#GLM method (recommend)
#The glm approach to multiple groups is similar to the classic approach, but permits more general comparisons to be made.
#if available, batch effect, as well as other factors, can be added in the design

#only keep counts data
rawcounts <- raw_readcount_data[,2:7]
#define all groups as factor with level
group <- names(rawcounts) %>% str_replace(.,"[0-9]$","") %>% str_replace(.,"^2","S2") %>% as.factor()
#create  DGEList
#the name "counts=" and "group=" are fixed
y = DGEList(counts = rawcounts,group = group)
#TMM normalization
y = calcNormFactors(y,method = "TMM")
#check if DGEList correct
y$samples
levels(y$samples$group)
#visualize replications and remove outlayers
plotMDS(y, col=rep(1:2, each=3))
#create the design
#"0+" here is needed, in this case, all groups can be used for comparisons
design <- model.matrix(~0+group, data = y$samples)
colnames(design) <- levels(y$samples$group)
#check design
design

#estimating GLM dispersions
#estimate common dispersion, trended dispersions and tagwise dispersions in one run (recommended)
y <- estimateDisp(y, design)
# #Alternatively, step-by-step
# y <- estimateGLMCommonDisp(y,design)
# y <- estimateGLMTrendedDisp(y,design)
# y <- estimateGLMTagwiseDisp(y,design)

#Dispersion means biological coeffient of variation (BCV) squared.
#Genes expression typically differs from replicate to replicate by 10% (BCV is 0.1), and its dispersion is 0.01
#value around 0.01 is acceptable
y$common.dispersion
plotBCV(y)

#GLM approach
fit <- glmQLFit(y, design)
#check groups
design
#set comparisons
#makeContrasts(A-B, levels=design): A is experimental group, B is control group
comparisons <- makeContrasts(S23m-S23a, levels=design) #<------ change here to set comparison
#define DEGs
qlf <- glmQLFTest(fit, contrast=comparisons)
plotQLDisp(qlf)
DE_glm <- qlf$table
DE_glm <- DE_glm %>%
  mutate(adj_pvalue = p.adjust(PValue,method = "BH"))

#set threshold and output DEGs
upDEGs_glm <- DE_glm %>%
  filter(logFC > 1, adj_pvalue < 0.01)
downDEGs_glm <- DE_glm %>%
  filter(logFC < -1, adj_pvalue < 0.01)
allDEGs_glm <- rbind(upDEGs_glm,downDEGs_glm)

#output data
comparision_id <- colnames(comparisons) %>% str_replace(.," - ","vs")
write.csv(upDEGs_glm,file= paste(comparision_id,dim(upDEGs_glm)[1],"upDEGs.csv",sep = "_"))
write.csv(downDEGs_glm,file= paste(comparision_id,dim(downDEGs_glm)[1],"downDEGs.csv",sep = "_"))
write.csv(allDEGs_glm,file= paste(comparision_id,dim(allDEGs_glm)[1],"allDEGs.csv",sep = "_"))
write.csv(DE_glm,file= "allgenes_edgeR.csv")

