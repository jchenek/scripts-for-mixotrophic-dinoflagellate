library(tidyverse)
library(pheatmap)

#get file list
getwd()
list.files()
chl_files <- list.files("chl")

summary_data <- data.frame(matrix(nrow = 0,ncol = 5))
names(summary_data) <- c("Year","Jan_chl_cov","Apr_chl_cov",
                         "Jul_chl_cov","Oct_chl_cov")


year <- 2003
i <- 1

while (year < 2022) {
  summary_data[i, 1] <- year
  for (month in c(1,4,7,10)) {
    if(month < 10){
      chl <- str_subset(chl_files, paste(".*",year,"-0",month,".*",sep = ""))
    }
    if(month >= 10){
      chl <- str_subset(chl_files, paste(".*",year,"-",month,".*",sep = ""))
    }
    chl_data <- read.csv(paste("chl/",chl,sep = ""),header=F,row.names = NULL)
    chl_data_nas <- chl_data
    
    chl_data_nas[chl_data > 9999] = NA
    
    chl_thread <- 0.3
    chl_high <- sum(chl_data_nas > chl_thread, na.rm = T)
    chl_all <- sum(chl_data_nas > 0, na.rm = T)
    chl_ratio <- round(chl_high/chl_all,2)
    
    j <- ((month - 1)/3) + 2
    summary_data[i,j] <- chl_ratio
  }
  i = i + 1
  year = year + 1
}
summary_data <- summary_data %>%
  mutate(ave_chl_cov = (Jan_chl_cov+Apr_chl_cov+Jul_chl_cov+Oct_chl_cov)/4)

#output data
write.csv(summary_data,file = "nasa_summary_chl_data.csv")

#brief plot and analysis
plot(summary_data$Year,summary_data$ave_chl_cov,type = "o")
fit <- lm(ave_chl_cov ~ Year, data = summary_data)
abline(fit)
summary(fit)

#plot heatmap
year <- 2021

for (month in c(1,4,7,10)) {
  if(month < 10){
    chl <- str_subset(chl_files, paste(".*",year,"-0",month,".*",sep = ""))
  }
  if(month >= 10){
    chl <- str_subset(chl_files, paste(".*",year,"-",month,".*",sep = ""))
  }
  chl_data <- read.csv(paste("chl/",chl,sep = ""),header=F,row.names = NULL)
  chl_data_map <- chl_data
  
  chl_data_map[chl_data > 9999] = NA
  chl_data_map[chl_data > 0 & chl_data <= 0.1] = 0.05
  chl_data_map[chl_data > 0.1 & chl_data <= 0.2] = 0.15
  chl_data_map[chl_data > 0.2 & chl_data <= 0.3] = 0.25
  chl_data_map[chl_data > 0.3 & chl_data <= 999] = 0.35
  
  map_breaks = c(0, 0.1, 0.2, 0.3, 0.4)
  
  pheatmap(chl_data_map,cellwidth = 0.5, cellheight = 0.5, 
           cluster_rows = F,cluster_cols = F,
           na_col = "#CCCCCC",
           color = colorRampPalette(c("white", "#CCFFCC", "#0066FF", "#FF0000"))(4),
           breaks = map_breaks,
           show_rownames = F, show_colnames = F,
           legend= F,
           main = paste(year,"-",month,sep = "")
  )
}

