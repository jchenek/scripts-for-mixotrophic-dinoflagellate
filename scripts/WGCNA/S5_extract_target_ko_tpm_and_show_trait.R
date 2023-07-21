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
library(psych)
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

write.csv(target_datExpr_trait,file="wgcna-output-target_datExpr_trait.csv")

#multi-plot and do linear regression
pdf(file = "wgcna_plot_target_gene_profile.pdf", wi = 12, he = 8)
par(mfrow = c(3,5))
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
  
  #spearman correlation
  corr_matrix <- corr.test(target_datExpr_trait[,i], target_datExpr_trait$temperature, method = 'spearman')
  my_cor_r <- round(corr_matrix$r,2)
  my_cor_p <- round(corr_matrix$p,2)
  
  #get gene name
  gene_name <- gene_list[names(target_datExpr_trait)[i],1]
  #plot
  plot(x=target_datExpr_trait$temperature, y=target_datExpr_trait[,i],
       col = "#6699FF80", pch= 16, cex=1.2, mgp = c(1.5, 0.5, 0),
       main = gene_name, col.main='#333333',font.main=4,
       xlab="Temperature ¡ãC",ylab="TPM", cex.lab=1, font.lab=2, col.lab='#333333',
       font.axis = 2,col.axis='#333333'
  )
  if(my_cor_p <= 0.05 & p_val <= 0.05){
  abline(fit, col = "#FF3333", lty = 1, lwd = 3)}
  if(my_cor_p > 0.05 & p_val < 0.05){
    abline(fit, col = "yellow", lty = 1, lwd = 3)}
  mtext(bquote(paste("R"^2," = ",.(r_sqa),", ",italic(p),"-value = ",.(p_val), 
                      " cor_r = ",.(my_cor_r)," cor_p = ",.(my_cor_p))),
        side = 3, adj = 0, font = 2, cex=0.5)
  i = i+1
}
dev.off()
