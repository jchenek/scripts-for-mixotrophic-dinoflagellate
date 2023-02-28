####################################################################################################################
library(tidyverse)
library(clusterProfiler)

####################################################################################################################
#prepare background gene info for enricher
#input prior full pathway and module
prior_pathway <- read.csv('./IN_enricher_prior/prekegg2gene_pathway.csv', header=F, check.names = FALSE)
colnames(prior_pathway) <- c("map_id","ko_id")

prior_module <- read.csv('./IN_enricher_prior/prekegg2gene_module.csv', header=F, check.names = FALSE)
colnames(prior_module) <- c("module_id","ko_id")

prior_pathway_info <- read.delim('./IN_enricher_prior/pathway_list_wb.txt', sep = '\t', header=F,
                                 stringsAsFactors = FALSE, check.names = FALSE)
colnames(prior_pathway_info) <- c("map_id1","map_info")
prior_pathway_info <- prior_pathway_info %>%
  mutate(map_id = str_replace(map_id1, "path:",""))

prior_module_info <- read.delim('./IN_enricher_prior/module_list_wb.txt', sep = '\t', header=F,
                                 stringsAsFactors = FALSE, check.names = FALSE)
colnames(prior_module_info) <- c("module_id1","module_info")
prior_module_info <- prior_module_info %>%
  mutate(module_id = str_replace(module_id1, "md:",""))

#input bg data and unify ko_id
bg_enricher <- read.csv('./IN_enr_bg/bg_enricher.csv', header=T, check.names = FALSE)
uni_ko_id <- bg_enricher %>%
  select(ko_id) %>%
  distinct(ko_id, .keep_all = TRUE)

#making pathway2gene and module2gene
pathway2gene <- prior_pathway %>%
  filter(ko_id %in% uni_ko_id[,1])
module2gene <- prior_module %>%
  filter(ko_id %in% uni_ko_id[,1])
####################################################################################################################
#input target gene set
enr_files <- c(list.files("./IN_enr_bg")) 
enr_files <- str_subset(enr_files,"^enr_")

#apply enricher one by one
dir.create("OU_enricher_pathway_output")
dir.create("OU_enricher_module_output") 
i = 1
while(i <= length(enr_files)){
  #path to file
  path <- paste("./IN_enr_bg/",enr_files[i],sep = "")
  
  #create input file
  d_enricher <- read.csv(path,header=T,stringsAsFactors=FALSE)
  gene_enricher <- d_enricher[,2]
  names(gene_enricher) <- as.character(d_enricher[,2])
  
  #enricher function for pathway
  enr_kegg_pathway <- enricher(gene_enricher, TERM2GENE = pathway2gene)
  enr_kegg_pathway_res <- enr_kegg_pathway@result %>%
    mutate(map_id = rownames(enr_kegg_pathway@result))
  enr_kegg_pathway_res <- left_join(enr_kegg_pathway_res,prior_pathway_info,by="map_id") 
  enr_kegg_pathway_res <- enr_kegg_pathway_res %>%
    select(-map_id1)
  rownames(enr_kegg_pathway_res) <- rownames(enr_kegg_pathway@result)
  
  path_out <- paste("./OU_enricher_pathway_output/",enr_files[i],sep = "")
  write.csv(enr_kegg_pathway_res,file= path_out) 
  
  #enricher function for module
  enr_kegg_module <- enricher(gene_enricher, TERM2GENE = module2gene)
  enr_kegg_module_res <- enr_kegg_module@result %>%
    mutate(module_id = rownames(enr_kegg_module@result))
  enr_kegg_module_res <- left_join(enr_kegg_module_res,prior_module_info,by="module_id")
  enr_kegg_module_res <- enr_kegg_module_res %>%
    select(-module_id1)
  rownames(enr_kegg_module_res) <- rownames(enr_kegg_module@result)
  
  path_out <- paste("./OU_enricher_module_output/",enr_files[i],sep = "")
  write.csv(enr_kegg_module_res,file= path_out)
  
  i = i+1
}



