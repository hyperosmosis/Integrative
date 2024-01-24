setwd("~/Documents/Wu/Weighted_Kernel/Real Data/")

library(tidyverse)
library(KEGGREST)

subject <- read.csv("SubjectData.csv", header = T, stringsAsFactors = F)
subject$Subject_ID <- sub("^", "X", subject$Subject_ID)

klist <- read.csv("K_List.csv", header = F, stringsAsFactors = F)
KO_data <- read.csv("KO_Data.csv", header = T, stringsAsFactors = F)
pathways <- list()
range <- (1:500)
for(i in 1:14){
  path <- keggLink("pathway", klist$V1[range])
  pathways[[i]] <- data.frame(KO_Terms=names(path), value=path, row.names=NULL)
  range <- range + 500
}

path <- keggLink("pathway", klist$V1[7001:7121])
pathways[[15]] <- data.frame(KO_Terms=names(path), value=path, row.names=NULL)
KO_paths <- do.call(rbind.data.frame, pathways)
KO_paths$KO_Terms <- sub("ko:", "", KO_paths$KO_Terms)
KO_paths$value <- sub("path:", "", KO_paths$value)
KO_paths <- KO_paths[- grep("ko", KO_paths$value),]
final_paths <- KO_paths[!duplicated(KO_paths$KO_Terms),]

KO_data <- merge(final_paths, KO_data, by = "KO_Terms")[,-1]
pathway_values <- KO_data %>%
  group_by(value) %>%
  summarise_all(funs(sum))

metab <- read.csv("Metabolite.csv", header = T, stringsAsFactors = F)
metab$Metabolite <- gsub("(.*)_.*", "\\1", metab$Metabolite)

mapList <- read.table("Pathways.txt", header = F, sep = "\t")
colnames(mapList) <- c("Pathway", "Metabolite")
metab.path <- merge(mapList, metab, by = "Metabolite")

# Subset to get same subjects
pathway_values <- pathway_values[,c(TRUE, colnames(pathway_values)[2:617]%in%colnames(metab.path)[3:408])]
metab.path <- metab.path[,c(TRUE, TRUE, colnames(pathway_values)[2:617]%in%colnames(metab.path)[3:408])]
subject <- subject[subject$Subject_ID%in%colnames(pathway_values)[2:348],]

## Remove rows with lots of zeros
miss.value <- sapply(pathway_values[,-1], function(x) min(x[x != 0])/sqrt(2))
for (i in 1:140){
  pathway_values[i, pathway_values[i,] == 0] <- miss.value[i] 
}

metab.path <- metab.path[rowSums(metab.path[,-c(1:2)] == 0) <= 261,] # 75%
miss.value <- apply(metab.path[,-c(1,2)], 1, function(x) min(x[x != 0])/sqrt(2))
for (i in 1:1262){
  metab.path[i, metab.path[i,] == 0] <- miss.value[i] 
}

metab.log <- metab.path
metab.log[,-c(1:2)] <- log(metab.log[,-c(1:2)])
save(metab.log, pathway_values, subject, file = "RealData.rData")

