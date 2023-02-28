library(vegan)
library(factoextra)
library(pheatmap)

getwd()
set.seed(520)


#input results from 'get_gene_enri_term_relation.pl'
data1 <- read.csv('gene_enri_term_relation.csv',header=T,row.names = 1)
#check data
rowSums(data1)
colSums(data1)

#kmeans clustering
#set an appropriate cluster number according to 'sum of squared error (SSE)' (wss method)
#k.max value can not be higher than nrow()
nrow(data1)
fviz_nbclust(data1,kmeans,k.max = 15, method = "wss",
             diss = vegdist(data1,method="euclidean"))+ 
  geom_vline(xintercept = 3, linetype = 2)


#set cluster number and clustering
center_num = 3
km_result <- kmeans(data1, centers = center_num, nstart = 100)
#plot1
fviz_cluster(km_result, data = data1,
             geom = c("point", "text"),
             show.clust.cent = T,
             ellipse.type = "euclid",
             star.plot = T,
             repel = TRUE,
             pointsize = 1.5, labelsize = 12, 
             main = "Cluster plot", xlab = NULL,
             ylab = NULL,
             ggtheme = theme_minimal()
)

#hierarchical clustering
x <- vegdist(data1,method="euclidean")
data2 <- hclust(x,"ward.D2")
#plot2
fviz_dend(data2, k = center_num, cex = 1, lwd = 1,
          main ="main",xlab = "xlab"
)

#Anosim analysis
Group <- km_result[["cluster"]] #set group based on kmeans result
data.anosim = anosim(x, grouping = Group, permutations = 999)
summary(data.anosim)

#output clustering results if the three analyses consist with each other
group_out <-as.data.frame(Group)
group_out[,2] <- row.names(group_out)
colnames(group_out) <- c("cluster","map_id")
group_out <- group_out[order(km_result[["cluster"]], decreasing = F),]
write.csv(group_out,file= "OU_km_cluster.csv",row.names = T)

#preparing heatmap
data2 <- data1
data2 <- data2[order(km_result[["cluster"]], decreasing = F),]
data2 <- data2[,order(colSums(data1),decreasing = T)]
#plot
anno_group <- as.data.frame(Group)
anno_group[,1] <- as.character(anno_group[,1])
pheatmap(data2, annotation_row = anno_group,
         cellwidth = 4, cellheight =12, 
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("#FFFFFF","#CC3333"))(10),
         legend=TRUE,  border="#CCCCCC", angle_col = c("315"),#fontface="italic",
         fontsize_row=9, fontsize_col=8,main=""
)
#output the dataframe of heatmap
write.csv(data2,file= "OU_gene_enri_term_frame.csv")

#output the dataframe of each cluster in a directory
dir.create("OU_cluster_ko")
for(i in 1:center_num){
  tar_cluster <- group_out[group_out$cluster == i,]
  data_cluster <- data2[match(rownames(tar_cluster), rownames(data2)),]
  data_cluster <- data_cluster[,which(colSums(data_cluster)>0)]
  data_cluster <- t(data_cluster)
  write.csv(data_cluster,file = paste("./OU_cluster_ko/data_cluster_", i,".csv", sep = ""))
}
