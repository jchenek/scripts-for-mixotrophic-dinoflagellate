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

#input ko list
gene_list <-  read.csv('IN_target_gene_list.csv', header = FALSE, row.names = 1)
ko_id <- as.character(row.names(gene_list))

target_datExpr <- datExpr %>%
  select(all_of(ko_id)) %>%
  mutate(gene_id = row.names(datExpr))

datTraits <- datTraits %>%
  mutate(gene_id = row.names(datTraits))

target_datExpr_trait <- inner_join(datTraits,target_datExpr, by="gene_id")

#multi-plot and do linear regression
pdf(file = "wgcna_plot_target_gene_profile_compare.pdf", wi = 24, he = 2)
par(mfrow = c(1,13),mar=c(1,0.5,1,0.5),oma=c(1,1,1,1))
#DO NOT change i
i = 7
while (i <= length(target_datExpr_trait[1,])) {
  #linear regression
  fit <- lm( target_datExpr_trait[,i] ~temperature, data = target_datExpr_trait)
  #get R squared and p value
  print(summary(fit))
  #cancel scientific notation
  options(scipen=200)
  r_sqa <- round(summary(fit)$r.squared,2)
  p_val <- round(summary(fit)$coefficients[2,4],4)
  #get gene name
  gene_name <- gene_list[names(target_datExpr_trait)[i],1]
  #plot
  if(p_val < 0.05){
  plot(x=target_datExpr_trait$temperature, y=target_datExpr_trait[,i],
       col = "#6699FF80", pch= 16, cex=3,
       main = NULL,
       xlab = "", ylab= "",
       xaxt = "n", yaxt = "n"
  )
  abline(fit, col = "#FF3333", lty = 1, lwd = 10)
  }
  if(p_val >= 0.05){
    plot(x=target_datExpr_trait$temperature, y=target_datExpr_trait[,i],
         col = "#6699FF80", pch= 16, cex=3,
         main = NULL,
         xlab = "", ylab= "",
         xaxt = "n", yaxt = "n"
    )
    abline(fit, col = "#CCCCCC", lty = 1, lwd = 10)
  }
  i = i+1
}
dev.off()

