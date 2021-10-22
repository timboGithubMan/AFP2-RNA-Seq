library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ggrepel)
library(cowplot)
library(BiocManager)
library(devtools)
#library(ComplexHeatmap)
#library(InteractiveComplexHeatmap)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationDbi)
library(org.At.tair.db)
library(enrichplot)
library(SPIA)
library(WGCNA)
library(dplyr)
library(openxlsx)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(ggpubr)
library(gprofiler2)
library(dendsort)
library(dendextend)
library(gplots)
register(SerialParam())
heat_colors <- rev(brewer.pal(11, "RdBu"))
write.xlsx <- function(data, filename) { 
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  conditionalFormatting(wb, "sheet", cols=9:19, rows = 1:nrow(data),
                        style = c("blue", "white", "red"),
                        rule = c(-5,0,5),
                        type = "colourScale")
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}
write.xlsx_noformatting <- function(data, filename) { 
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}
write.xlsx_GO <- function(data, filename) { 
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  conditionalFormatting(wb, "sheet", cols=4:5, rows = 1:nrow(data),
                        style = createStyle(bgFill = "#00B050"),
                        rule = "<=0.05",
                        type = "expression")
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}
write.xlsx_GO_GSEA <- function(data, filename) { 
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  conditionalFormatting(wb, "sheet", cols=3:4, rows = 1:nrow(data),
                        style = createStyle(bgFill = "#00B050"),
                        rule = "<=0.05",
                        type = "expression")
  conditionalFormatting(wb, "sheet", cols=5:6, rows = 1:nrow(data),
                        style = c("blue", "white", "red"),
                        rule = c(-3,0,3),
                        type = "colourScale")
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

genbankToTair <- read.delim(file="data/genbank to tair.txt")
count_data = read.csv(file= "data/RawReadsCY.csv") %>%
  left_join(genbankToTair, by=c("Transcript" = "From"))
count_data = count_data[,c(21,3:20)]
count_data2 <- rename(count_data, geneID=To)
count_data2 <- count_data2[!duplicated(count_data2$geneID),]
count_data2 <- na.omit(count_data2)

count_data <- data.frame(count_data2[,-1], row.names=count_data2[,1])

col_data = read.csv(file= "data/col_data.csv")

annot=read.csv(file="data/annot.csv")[,c(2,3)]
entrez=read.csv(file="data/entrez.csv")[,c(2,3,4)]

k12_5_entrezlfc <- all_LFCs$LFC_YSA_YvsCol
names(k12_5_entrezlfc) <- all_LFCs$Entrez

entrez_col_TSAvsGM <- all_LFCs$LFC_Col_TSAvsGM
names(entrez_col_TSAvsGM) <- all_LFCs$Entrez
entrez_col_ABAvsGM <- all_LFCs$LFC_Col_ABAvsGM
names(entrez_col_ABAvsGM) <- all_LFCs$Entrez
entrez_tsa_YvsCOL <- all_LFCs$LFC_YSA_YvsCol
names(entrez_tsa_YvsCOL) <- all_LFCs$Entrez
entrez_gm_YvsCOL <- all_LFCs$LFC_GM_YvsCol
names(entrez_gm_YvsCOL) <- all_LFCs$Entrez

pathview(gene.data = entrez_tsa_YvsCOL,
               out.suffix = "tsa_YvsCOL",
               pathway.id = "ath03020",
               species = "ath",
               low = list(gene = "red", cpd = "yellow"),
               high = list(gene = "green", cpd = "blue"),
               limit = list(gene = 2, # value gives the max/min limit for foldchanges
                            cpd = 1))


##create data sets
ddsColABA = DESeqDataSetFromMatrix(countData = count_data,
                             colData = col_data,
                             design = ~Type)
ddsColGM = DESeqDataSetFromMatrix(countData = count_data,
                             colData = col_data,
                             design = ~Type)
ddsColTSA <- DESeqDataSetFromMatrix(countData = count_data,
                             colData = col_data,
                             design = ~Type)
ddsYGM <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = col_data,
                                design = ~Type)
ddsYTSA <- DESeqDataSetFromMatrix(countData = count_data,
                                 colData = col_data,
                                 design = ~Type)
ddsColGM$Type <- relevel(ddsColGM$Type, ref="ColGM")
ddsColTSA$Type <- relevel(ddsColTSA$Type, ref="ColTSA")
ddsYGM$Type <-relevel(ddsYGM$Type, ref="YGM")
ddsYTSA$Type <-relevel(ddsYTSA$Type, ref="YTSA")

ddsColABA <- DESeq(ddsColABA)
ddsColGM <- DESeq(ddsColGM)
ddsColTSA <- DESeq(ddsColTSA)
ddsYGM <- DESeq(ddsYGM)
ddsYTSA <- DESeq(ddsYTSA)

normalized_counts <- counts(ddsColGM, normalized=T) %>% #create tibbles of normalized counts, annotated by gene symbol
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

testmin <- apply(normalized_counts[2:19],1,min)
testmax <- apply(normalized_counts[2:19],1,max)
testscounts <- normalized_counts[order(testmax-testmin, decreasing=TRUE),][1:100,]

scaled_counts <- scale_rows(normalized_counts[2:19]) %>% cbind(normalized_counts$gene) 
scaled_counts <- scaled_counts[,c(19,1:18)]
scaled_countsMEANS <- scale_rows(count_means[2:7]) %>% cbind(count_means$gene) 
scaled_countsMEANS <- scaled_countsMEANS[,c(7,1:6)] 
scaled_countsMEANS <- rename(scaled_countsMEANS, gene = 'count_means$gene')
scaled_countsMEDIANS <- scale_rows(count_medians[2:7]) %>% cbind(count_medians$gene) 
scaled_countsMEDIANS <- scaled_countsMEDIANS[,c(7,1:6)] 
scaled_countsMEDIANS <- rename(scaled_countsMEDIANS, gene = 'count_medians$gene')

count_means_var <- apply(count_means[2:7], 1, var)
count_means_mean <- apply(count_means[2:7], 1, mean)
count_means_filtered <- count_means %>%
  filter(gene %in% sigresColABAColGM$gene 
         | gene %in% sigresColABAColTSA$gene 
         | gene %in% sigresColTSAColGM$gene 
         | gene %in% sigresYABAColABA$gene 
         | gene %in% sigresYABAYGM$gene 
         | gene %in% sigresYABAYTSA$gene 
         | gene %in% sigresYGMColGM$gene 
         | gene %in% sigresYTSAColTSA$gene 
         | gene %in% sigresYTSAYGM$gene)

count_medians_filtered <- count_medians %>%
  filter(gene %in% sigresColABAColGM$gene 
         | gene %in% sigresColABAColTSA$gene 
         | gene %in% sigresColTSAColGM$gene 
         | gene %in% sigresYABAColABA$gene 
         | gene %in% sigresYABAYGM$gene 
         | gene %in% sigresYABAYTSA$gene 
         | gene %in% sigresYGMColGM$gene 
         | gene %in% sigresYTSAColTSA$gene 
         | gene %in% sigresYTSAYGM$gene)

scaled_countsMEANS_filtered <- scaled_countsMEANS %>% filter(gene %in% count_means_filtered$gene)
scaled_countsMEDIANS_filtered <- scaled_countsMEDIANS %>% filter(gene %in% count_medians_filtered$gene)

corrAVGHCmedians_filtered<- hclust(as.dist(1- cor(t(count_medians_filtered[2:7]))), method="average")
corrAVGHCmeans_filtered<- hclust(as.dist(1- cor(t(count_means_filtered[2:7]))), method="average")
euclAVGhcscaledmeans_filtered <- hclust(dist(as.matrix(scaled_countsMEANS_filtered[2:7], method="euclidean")), method="average")
corrAVGHCmeansSpearman_filtered<- hclust(as.dist(1- cor(t(count_means_filtered[2:7]), method="spearman")), method="average")
corrAVGHC<- hclust(as.dist(1- cor(t(normalized_counts[2:19]))), method="average")
corrAVGHCmeans<- hclust(as.dist(1- cor(t(count_means[2:7]))), method="average")
corrAVGHCmeansComplete<- hclust(as.dist(1- cor(t(count_means[2:7]), method="pearson")), method="complete")
euclAVGhcscaled <- hclust(dist(as.matrix(scaled_counts[2:19], method="euclidean")), method="average")
euclAVGhcscaledmeans <- hclust(dist(as.matrix(scaled_countsMEANS[2:7], method="euclidean")), method="average")
corrAVGHCmeansSpearman<- hclust(as.dist(1- cor(t(count_means[2:7]), method="spearman")), method="average")
corrAVGHCmeansSpearmanComplete<- hclust(as.dist(1- cor(t(count_means[2:7]), method="spearman")), method="complete")

pearsonAVGmedianssorted <- dendsort(as.dendrogram(corrAVGHCmedians_filtered), type="average")
order_corrAVGHCmedians_filtered_sorted <- scaled_countsMEDIANS_filtered[as.hclust(pearsonAVGmedianssorted)$order,]
order_corrAVGHCmedians_filtered <- scaled_countsMEDIANS_filtered[corrAVGHCmedians_filtered$order,]
colnames(order_corrAVGHCmedians_filtered)[2] <- "ColGM"
colnames(order_corrAVGHCmedians_filtered)[3] <- "ColABA"
colnames(order_corrAVGHCmedians_filtered)[4] <- "ColTSA"
colnames(order_corrAVGHCmedians_filtered)[5] <- "YGM"
colnames(order_corrAVGHCmedians_filtered)[6] <- "YABA"
colnames(order_corrAVGHCmedians_filtered)[7] <- "YTSA"
order_corrAVGHCmeans_filtered <- scaled_countsMEANS_filtered[corrAVGHCmeans_filtered$order,]
order_euclAVGhcscaledmeans_filtered <- scaled_countsMEANS_filtered[euclAVGhcscaledmeans_filtered$order,]
order_corrAVGHCmeansSpearman_filtered <- scaled_countsMEANS_filtered[corrAVGHCmeansSpearman_filtered$order,]
order_corrAVGHC <- scaled_counts[corrAVGHC$order,]
order_corrAVGHCmeans <- scaled_countsMEANS[corrAVGHCmeans$order,]
order_corrAVGHCmeansComplete <- scaled_countsMEANS[corrAVGHCmeansComplete$order,]
order_euclAVGhcscaled <- scaled_counts[euclAVGhcscaled$order,]
order_corrAVGHCmeansSpearman <- scaled_countsMEANS[corrAVGHCmeansSpearman$order,]
order_corrAVGHCmeansSpearmanComplete <- scaled_countsMEANS[corrAVGHCmeansSpearmanComplete$order,]

zerod_scaled_countsMEANS_filtered <- scaled_countsMEANS_filtered
zerod_order_corrAVGHCmeans_filtered <- order_corrAVGHCmeans_filtered
zerod_order_euclAVGHCscaledmeans_filtered  <- order_euclAVGhcscaledmeans_filtered 
zerod_order_corrAVGHCmeansSpearman_filtered <- order_corrAVGHCmeansSpearman_filtered
#zerod_order_corrAVGHCmeans <- order_corrAVGHCmeans
#zerod_order_euclAVGHCscaledmeans <- order_euclAVGhcscaledmeans
#zerod_order_corrAVGHCmeansSpearman <- order_corrAVGHCmeansSpearman
#zerod_order_corrAVGHCmeansSpearmanComplete <- order_corrAVGHCmeansSpearmanComplete

zerod_scaled_countsMEANS_filtered$photosynthesis <- zerod_scaled_countsMEANS_filtered$YTSA * (scaled_countsMEANS_filtered$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
zerod_order_corrAVGHCmeans_filtered$photosynthesis <- zerod_order_corrAVGHCmeans_filtered$YTSA * (order_corrAVGHCmeans_filtered$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
zerod_order_euclAVGHCscaledmeans_filtered$photosynthesis <- zerod_order_euclAVGHCscaledmeans_filtered$YTSA * (order_euclAVGhcscaledmeans_filtered$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
zerod_order_corrAVGHCmeansSpearman_filtered$photosynthesis <- zerod_order_corrAVGHCmeansSpearman_filtered$YTSA * (order_corrAVGHCmeansSpearman_filtered$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
#zerod_order_corrAVGHCmeans$photosynthesis <- zerod_order_corrAVGHCmeans$YTSA * (order_corrAVGHCmeans$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
#zerod_order_euclAVGHCscaledmeans$photosynthesis <- zerod_order_euclAVGHCscaledmeans$YTSA * (order_euclAVGhcscaledmeans$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
#zerod_order_corrAVGHCmeansSpearman$photosynthesis <- zerod_order_corrAVGHCmeansSpearman$YTSA * (order_corrAVGHCmeansSpearman$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000
#zerod_order_corrAVGHCmeansSpearmanComplete$photosynthesis <- zerod_order_corrAVGHCmeansSpearmanComplete$YTSA * (order_corrAVGHCmeansSpearmanComplete$gene %in% pathways_GO[["PHOTOSYNTHESIS"]]) * 10000

zerod_scaled_countsMEANS_filtered$SEED_MATURATION <- zerod_scaled_countsMEANS_filtered$YTSA * (scaled_countsMEANS_filtered$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
zerod_order_corrAVGHCmeans_filtered$SEED_MATURATION <- zerod_order_corrAVGHCmeans_filtered$YTSA * (order_corrAVGHCmeans_filtered$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
zerod_order_euclAVGHCscaledmeans_filtered$SEED_MATURATION <- zerod_order_euclAVGHCscaledmeans_filtered$YTSA * (order_euclAVGhcscaledmeans_filtered$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
zerod_order_corrAVGHCmeansSpearman_filtered$SEED_MATURATION <- zerod_order_corrAVGHCmeansSpearman_filtered$YTSA * (order_corrAVGHCmeansSpearman_filtered$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
#zerod_order_corrAVGHCmeans$SEED_MATURATION <- zerod_order_corrAVGHCmeans$YTSA * (order_corrAVGHCmeans$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
#zerod_order_euclAVGHCscaledmeans$SEED_MATURATION <- zerod_order_euclAVGHCscaledmeans$YTSA * (order_euclAVGhcscaledmeans$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
#zerod_order_corrAVGHCmeansSpearman$SEED_MATURATION <- zerod_order_corrAVGHCmeansSpearman$YTSA * (order_corrAVGHCmeansSpearman$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000
#zerod_order_corrAVGHCmeansSpearmanComplete$SEED_MATURATION <- zerod_order_corrAVGHCmeansSpearmanComplete$YTSA * (order_corrAVGHCmeansSpearmanComplete$gene %in% pathways_GO[["SEED_MATURATION"]]) * 10000

zerod_scaled_countsMEANS_filtered$chromatin_modification <- zerod_scaled_countsMEANS_filtered$YTSA * (scaled_countsMEANS_filtered$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
zerod_order_corrAVGHCmeans_filtered$chromatin_modification <- zerod_order_corrAVGHCmeans_filtered$YTSA * (order_corrAVGHCmeans_filtered$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
zerod_order_euclAVGHCscaledmeans_filtered$chromatin_modification <- zerod_order_euclAVGHCscaledmeans_filtered$YTSA * (order_euclAVGhcscaledmeans_filtered$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
zerod_order_corrAVGHCmeansSpearman_filtered$chromatin_modification <- zerod_order_corrAVGHCmeansSpearman_filtered$YTSA * (order_corrAVGHCmeansSpearman_filtered$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
#zerod_order_corrAVGHCmeans$chromatin_modification <- zerod_order_corrAVGHCmeans$YTSA * (order_corrAVGHCmeans$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
#zerod_order_euclAVGHCscaledmeans$chromatin_modification <- zerod_order_euclAVGHCscaledmeans$YTSA * (order_euclAVGhcscaledmeans$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
#zerod_order_corrAVGHCmeansSpearman$chromatin_modification <- zerod_order_corrAVGHCmeansSpearman$YTSA * (order_corrAVGHCmeansSpearman$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000
#zerod_order_corrAVGHCmeansSpearmanComplete$chromatin_modification <- zerod_order_corrAVGHCmeansSpearmanComplete$YTSA * (order_corrAVGHCmeansSpearmanComplete$gene %in% pathways_GO[["CHROMATIN_MODIFICATION"]]) * 10000

kmeans_countmeans <- kmeans(count_means[2:7], centers = 720)
image(t(count_means[2:7])[, order(kmeans_countmeans$cluster)], main = "Clustered Data")

kmeans_GO <- kmeans(GO_k12[5:16], centers = 720)
image(t(GO_k12[5:16])[, order(kmeans_GO$cluster)], main = "Clustered Data")

pheatmap(order_corrAVGHCmedians_filtered_sorted[,c(2,5,3,6,4,7)],
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-2, 2, length.out = 101),
         show_rownames = F,
         fontsize = 10,
         scale = "none",
         fontsize_row = 10,
         filename = "order_corrAVGHCmedians_filtered_sorted.png",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         height = 15)

pheatmap(zerod_order_corrAVGHCmeans_filtered[,c(2,5,3,6,4,7,8,9,10)],
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-2, 2, length.out = 101),
         show_rownames = F,
         fontsize = 10,
         scale = "none",
         fontsize_row = 10,
         filename = "pearson_average_means_filtered.png",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         height = 15)
Heatmap(zerod_scaled_countsMEANS_filtered[,c(2,5,3,6,4,7,8,9,10)], name = "z-score", cluster_rows = corrAVGHCmeans_filtered,
        cluster_columns = FALSE, column_title = "Pearson Average")
pheatmap(zerod_order_euclAVGHCscaledmeans[,c(2,5,3,6,4,7,8,9,10)],
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-3, 3, length.out = 101),
         show_rownames = F,
         fontsize = 10,
         scale = "none",
         fontsize_row = 10,
         filename = "euclideanscaled_average_means_filtered.png",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         height = 15)
pheatmap(zerod_order_corrAVGHCmeansSpearman_filtered[,c(2,5,3,6,4,7,8,9,10)],
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-3, 3, length.out = 101),
         show_rownames = F,
         fontsize = 10,
         scale = "none",
         fontsize_row = 10,
         filename = "spearman_average_means_filtered.png",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         height = 15)

pearsonAVGmeanssorted <- dendsort(as.dendrogram(corrAVGHCmeans_filtered), type="average")
dend <- color_branches(pearsonAVGmeanssorted, k=20)

dendcut0.35 = cutree(pearsonAVGmedianssorted, h=0.35) 
clusDyn <- cutreeDynamic(corrAVGHCmedians_filtered, distM = as.matrix(as.dist(1- cor(t(count_medians_filtered[2:7]))), method = "hybrid"))

clusters <- cutreeDynamic(corrAVGHCmeans_filtered, distM = as.matrix(as.dist(1- cor(t(count_means_filtered[2:7]))), method = "tree"))
clusters <- clusters[order.dendrogram(as.dendrogram(corrAVGHCmeans_filtered))]
clusters_numbers <- unique(clusters) - (0 %in% clusters)
n_clusters <- length(clusters_numbers)
dend2 <- corrAVGHCmeans_filtered %>% as.dendrogram %>% 
  branches_attr_by_clusters(clusters, values = cols) 

SSE <- (nrow(order_corrAVGHCmedians_filtered[2:7])-1)*sum(apply(order_corrAVGHCmedians_filtered[2:7],2,var))
for (i in 2:20) SSE[i] <- sum(kmeans(order_corrAVGHCmedians_filtered[2:7],
                                     centers=i)$withinss)
SSE <- (nrow(GO_k12[5:16])-1)*sum(apply(GO_k12[5:16],2,var))
for (i in 2:20) SSE[i] <- sum(kmeans(GO_k12[5:16],
                                     centers=i)$withinss)

png("etstset.png", width = 400, height =200)

plot(1:20, SSE, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

library(cluster)
sil <- rep(0, 20)
for(i in 2:20){
  k1to20 <- kmeans(order_corrAVGHCmedians_filtered[2:7], centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(order_corrAVGHCmedians_filtered[2:7]))
  sil[i] <- mean(ss[, 3])
}
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

library(vegan)
fit <- cascadeKM(order_corrAVGHCmedians_filtered[2:7], 1, 20, iter = 100)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)

library(cluster)
gap <- clusGap(order_corrAVGHCmedians_filtered[2:7], kmeans, 20, B = 100, verbose = interactive())
plot(gap, main = "Gap statistic")
abline(v=which.max(gap$Tab[,3]), lty = 2)


colnames(order_corrAVGHCmedians_filtered)[2] <- "1ColGM"
colnames(order_corrAVGHCmedians_filtered)[3] <- "3ColABA"
colnames(order_corrAVGHCmedians_filtered)[4] <- "5ColTSA"
colnames(order_corrAVGHCmedians_filtered)[5] <- "2YGM"
colnames(order_corrAVGHCmedians_filtered)[6] <- "4YABA"
colnames(order_corrAVGHCmedians_filtered)[7] <- "6YTSA"
kClust <- kmeans(order_corrAVGHCmedians_filtered[2:7], centers=6, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster

clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, order_corrAVGHCmedians_filtered[2:7], kClusters)

library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

#plot
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Condition") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Condition",color = "Cluster")
p1

kClust <- kmeans(order_corrAVGHCmedians_filtered[2:7], centers=6, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
order_corrAVGHCmedians_filtered_k6<- data.frame(order_corrAVGHCmedians_filtered, kClusters) %>% arrange(kClusters)
pheatmap(scaledata_k12_arranged[,c(2,5,3,6,4,7)],
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-2, 2, length.out = 101),
         show_rownames = F,
         fontsize = 10,
         scale = "none",
         fontsize_row = 10,
         filename = "medians_filtered_k12_arranged.png",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         height = 15)
pheatmap(scaledata_k12_arranged[8],
         show_rownames = F,
         fontsize = 10,
         scale = "none",
         fontsize_row = 10,
         filename = "medians_filtered_k12_arranged_bar.png",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         height = 15)

kClustGO <- kmeans(GO_k12[5:16], centers=8, nstart = 1000, iter.max = 20)
kClustersGO <- kClustGO$cluster
GO_k12_k8 <- data.frame(GO_k12, kClustersGO) %>% arrange(kClustersGO)

#write functions to save xlsx files (without and with formatting(colored boxes))
write.xlsx_noformatting <- function(data, filename) { 
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}
write.xlsx_GO_kmeans <- function(data, filename) { 
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  conditionalFormatting(wb, "sheet", cols=5:16, rows = 1:nrow(data),
                        style = c("green", "white"),
                        rule = c(0,0.10),
                        type = "colourScale")
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}

scaledata_k15 <- order_corrAVGHCmedians_filtered_k12
torgenes <- read.csv("torgenes.csv")[2:3]
scaledata_k15_TOR <- scaledata_k15 %>% filter(gene %in% torgenes$gene)
numClusters <- scaledata_k15[which.max(scaledata_k15$kClusters),c("kClusters")]
percentTOR <- array(1:numClusters, dim=c(1,numClusters))
for (i in 1:numClusters) {
  percentTOR[i] <- sum(scaledata_k15_TOR$kClusters == i)/sum(scaledata_k15$kClusters == i)
}

scaledata_k12_arranged <- arrange(scaledata_k12, factor(kClusters, levels = c(11,10,9,5,4,7,1,2,12,3,6,8)))
scaledata_k12_arranged$kClusters <- scaledata_k12_arranged$kClusters*10
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 110] <- 1.1
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 100] <- 1.2
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 90] <- 2.0
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 50] <- 3.1
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 40] <- 3.2
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 70] <- 3.3
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 10] <- 4.0
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 20] <- 5.0
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 120] <- 6.0
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 30] <- 7.0
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 60] <- 8.0
scaledata_k12_arranged$kClusters[scaledata_k12_arranged$kClusters == 80] <- 9.0

#perform GO analysis with g:profiler gost function and save spreadsheets
gost_k12_4 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==1, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_4$result, "gost_k12_4.xlsx")
gost_k12_5 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==2, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_5$result, "gost_k12_5.xlsx")
gost_k12_7 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==3, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_7$result, "gost_k12_7.xlsx")
gost_k12_3b <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==4, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_3b$result, "gost_k12_3b.xlsx")
gost_k12_3 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==5, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_3$result, "gost_k12_3.xlsx")
gost_k12_8 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==6, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_8$result, "gost_k12_8.xlsx")
gost_k12_3c <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==7, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_3c$result, "gost_k12_3c.xlsx")
gost_k12_9 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==8, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_9$result, "gost_k12_9.xlsx")
gost_k12_2 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==9, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_2$result, "gost_k12_2.xlsx")
gost_k12_1b <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==10, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_1b$result, "gost_k12_1b.xlsx")
gost_k12_1 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==11, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_1$result, "gost_k12_1.xlsx")
gost_k12_6 <- gost(order_corrAVGHCmedians_filtered_k12[which (order_corrAVGHCmedians_filtered_k12$kClusters==12, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k12_6$result, "gost_k12_6.xlsx")


gost_k12_4[["result"]][["intersection"]]
#make sheet combining p values for all clusters

GO_k12_avgLFC <- gost_k12_4[["result"]][,c(9,16)] %>%
  full_join(gost_k12_5[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_7[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_3b[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_3[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_8[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_3c[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_9[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_2[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_1b[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_1[["result"]][,c(9,16)], by=c('term_id')) %>%
  full_join(gost_k12_6[["result"]][,c(9,16)], by=c('term_id'))
colnames(GO_k12_avgLFC)[2] <- "4"
colnames(GO_k12_avgLFC)[3] <- "5"
colnames(GO_k12_avgLFC)[4] <- "7"
colnames(GO_k12_avgLFC)[5] <- "3b"
colnames(GO_k12_avgLFC)[6] <- "3"
colnames(GO_k12_avgLFC)[7] <- "8"
colnames(GO_k12_avgLFC)[8] <- "3c"
colnames(GO_k12_avgLFC)[9] <- "9"
colnames(GO_k12_avgLFC)[10] <- "2"
colnames(GO_k12_avgLFC)[11] <- "1b"
colnames(GO_k12_avgLFC)[12] <- "1"
colnames(GO_k12_avgLFC)[13] <- "6"

GO_k12 <- gost_k12_4[["result"]][,c(9,10,11,4,3)] %>%
  full_join(gost_k12_5[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_7[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_3b[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_3[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_8[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_3c[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_9[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_2[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_1b[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_1[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k12_6[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size'))
colnames(GO_k12)[5] <- "4_pval"
colnames(GO_k12)[6] <- "5_pval"
colnames(GO_k12)[7] <- "7_pval"
colnames(GO_k12)[8] <- "3b_pval"
colnames(GO_k12)[9] <- "3_pval"
colnames(GO_k12)[10] <- "8_pval"
colnames(GO_k12)[11] <- "3c_pval"
colnames(GO_k12)[12] <- "9_pval"
colnames(GO_k12)[13] <- "2_pval"
colnames(GO_k12)[14] <- "1b_pval"
colnames(GO_k12)[15] <- "1_pval"
colnames(GO_k12)[16] <- "6_pval"
GO_k12[is.na(GO_k12)] <- 1
GO_k12$sign <- as.numeric(as.character(unlist(GO_k12[5])))*as.numeric(as.character(unlist(GO_k12[6])))*as.numeric(as.character(unlist(GO_k12[7])))*as.numeric(as.character(unlist(GO_k12[8])))*as.numeric(as.character(unlist(GO_k12[9])))*as.numeric(as.character(unlist(GO_k12[10])))*as.numeric(as.character(unlist(GO_k12[11])))*as.numeric(as.character(unlist(GO_k12[12])))*as.numeric(as.character(unlist(GO_k12[13])))*as.numeric(as.character(unlist(GO_k12[14])))*as.numeric(as.character(unlist(GO_k12[15])))*as.numeric(as.character(unlist(GO_k12[16])))
GO_k12 <- arrange(GO_k12, sign)
write.xlsx_GO_kmeans(GO_k12_k8,"GO_k12_k8.xlsx")

png(file = "pearsonavgmeanssorted_dendsort.png", bg="transparent", width=600,height=4455)
#plot(rev(dend), leaflab="none", horiz=TRUE)
plot(pearsonAVGmedianssorted,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(clusDyn, pearsonAVGmedianssorted, sort_by_labels_order = T, y_shift=-0.1)

plot(as.dendrogram(corrAVGHCmedians_filtered),
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(clusDyn, corrAVGHCmedians_filtered, sort_by_labels_order = T, y_shift=-0.1)

dev.off()

pheatmap::pheatmap(zerod_scaled_countsMEANS_filtered[,c(2,5,3,6,4,7,8,9,10)],
                   color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                   breaks = seq(-2, 2, length.out = 101),
                   show_rownames = F,
                   fontsize = 10,
                   scale = "none",
                   fontsize_row = 10,
                   cluster_cols=FALSE,
                   cluster_rows=as.hclust(dend),
                   treeheight_row = 0,
                   height = 15,
                   width=8,
                   filename = "pearson.png")


png(file = "heatmap3.png")
heatmap.2(as.matrix(zerod_scaled_countsMEANS_filtered[,c(2,5,3,6,4,7,8,9,10)]),
          dendrogram="row",
          Rowv = dend,
          Colv = "NA",
          col=bluered(100),
          trace="none",
          scale="none",
)
dev.off()

count_medians <- normalized_counts[,c(1)]
count_medians$ColGM = normalized_counts[,c("col_g1","col_g2","col_g3")] %>% 
  rowwise() %>% 
  mutate(median = median(c(col_g1, col_g2, col_g3), na.rm = TRUE)) %>%
  select(median)
count_medians$ColABA = normalized_counts[,c("col_a1","col_a2","col_a3")] %>% 
  rowwise() %>% 
  mutate(median = median(c(col_a1, col_a2, col_a3), na.rm = TRUE)) %>%
  select(median)
count_medians$ColTSA =  normalized_counts[,c("col_t1","col_t2","col_t3")] %>% 
  rowwise() %>% 
  mutate(median = median(c(col_t1, col_t2, col_t3), na.rm = TRUE)) %>%
  select(median)
count_medians$YGM =  normalized_counts[,c("y_g1","y_g2","y_g3")] %>% 
  rowwise() %>% 
  mutate(median = median(c(y_g1, y_g2, y_g3), na.rm = TRUE)) %>%
  select(median)
count_medians$YABA =  normalized_counts[,c("y_a1","y_a2","y_a3")] %>% 
  rowwise() %>% 
  mutate(median = median(c(y_a1, y_a2, y_a3), na.rm = TRUE)) %>%
  select(median)
count_medians$YTSA =  normalized_counts[,c("y_t1","y_t2","y_t3")] %>% 
  rowwise() %>% 
  mutate(median = median(c(y_t1, y_t2, y_t3), na.rm = TRUE)) %>%
  select(median)

count_means <- normalized_counts[,c(1)]
count_means$ColGM = (normalized_counts$col_g1+normalized_counts$col_g2+normalized_counts$col_g3)/3
count_means$ColABA = (normalized_counts$col_a1+normalized_counts$col_a2+normalized_counts$col_a3)/3
count_means$ColTSA = (normalized_counts$col_t1+normalized_counts$col_t2+normalized_counts$col_t3)/3
count_means$YGM = (normalized_counts$y_g1+normalized_counts$y_g2+normalized_counts$y_g3)/3
count_means$YABA = (normalized_counts$y_a1+normalized_counts$y_a2+normalized_counts$y_a3)/3
count_means$YTSA = (normalized_counts$y_t1+normalized_counts$y_t2+normalized_counts$y_t3)/3

## Y-AFP vs Col
resYABAColABA <- results(ddsColABA, contrast=c("Type", "YABA", "ColABA"), alpha = 0.05, tidy=TRUE)
resYGMColGM <- results(ddsColGM, contrast=c("Type", "YGM", "ColGM"), alpha = 0.05, tidy=TRUE)
resYTSAColTSA <- results(ddsColTSA, contrast=c("Type", "YTSA", "ColTSA"), alpha = 0.05, tidy=TRUE)

resYABAColABAlfc <- lfcShrink(ddsColABA, coef="Type_YABA_vs_ColABA", res=resYABAColABA)
resYGMColGMlfc <- lfcShrink(ddsColGM, coef="Type_YGM_vs_ColGM", res=resYGMColGM)
resYTSAColTSAlfc <- lfcShrink(ddsColTSA, coef="Type_YTSA_vs_ColTSA", res=resYTSAColTSA)

## ABAvsGM and TSAvsGM
resColABAColGM <- results(ddsColGM, contrast=c("Type", "ColABA", "ColGM"), alpha = 0.05, tidy=TRUE)
resYABAYGM <- results(ddsYGM, contrast=c("Type", "YABA", "YGM"), alpha = 0.05)
resColTSAColGM <- results(ddsColABA, contrast=c("Type", "ColTSA", "ColGM"), alpha = 0.05)
resYTSAYGM <- results(ddsYGM, contrast=c("Type", "YTSA", "YGM"), alpha = 0.05)

resColABAColGMlfc <- lfcShrink(ddsColGM, coef="Type_ColABA_vs_ColGM", res=resColABAColGM)
resYABAYGMlfc <- lfcShrink(ddsYGM, coef="Type_YABA_vs_YGM", res=resYABAYGM)
resColTSAColGMlfc <- lfcShrink(ddsColGM, coef="Type_ColTSA_vs_ColGM", res=resColTSAColGM)
resYTSAYGMlfc <- lfcShrink(ddsYGM, coef="Type_YTSA_vs_YGM", res=resYTSAYGM)

## ABAvsTSA
resColABAColTSA <- results(ddsColTSA, contrast=c("Type", "ColABA", "ColTSA"), alpha = 0.05, lfcThreshold = 0.41)
resYABAYTSA <- results(ddsYTSA, contrast=c("Type", "YABA", "YTSA"), alpha = 0.05, lfcThreshold = 0.41)

resColABAColTSAlfc <- lfcShrink(ddsColTSA, coef="Type_ColABA_vs_ColTSA", res=resColABAColTSA)
resYABAYTSAlfc <- lfcShrink(ddsYTSA, coef="Type_YABA_vs_YTSA", res=resYABAYTSA)

resYABAColABA_tb <- resYABAColABAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
  sigresYABAColABA <- resYABAColABA_tb %>%
    filter(padj < 0.05)
resYGMColGM_tb <- resYGMColGMlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresYGMColGM <- resYGMColGM_tb %>%
    filter(padj < 0.05)
resYTSAColTSA_tb <- resYTSAColTSAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresYTSAColTSA <- resYTSAColTSA_tb %>%
    filter(padj < 0.05)

resColABAColGM_tb <- resColABAColGMlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresColABAColGM <- resColABAColGM_tb %>%
    filter(padj < 0.05)
resYABAYGM_tb <- resYABAYGMlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresYABAYGM <- resYABAYGM_tb %>%
    filter(padj < 0.05)
resColTSAColGM_tb <- resColTSAColGMlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresColTSAColGM <- resColTSAColGM_tb %>%
    filter(padj < 0.05)
resYTSAYGM_tb <- resYTSAYGMlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresYTSAYGM <- resYTSAYGM_tb %>%
    filter(padj < 0.05)

resColABAColTSA_tb <- resColABAColTSAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresColABAColTSA <- resColABAColTSA_tb %>%
    filter(padj < 0.05)
resYABAYTSA_tb <- resYABAYTSAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From")) 
  sigresYABAYTSA <- resYABAYTSA_tb %>%
    filter(padj < 0.05)

master_sheet <- count_means %>%
  left_join(resColABAColGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYABAYGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resColTSAColGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYTSAYGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYGMColGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYABAColABA_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYTSAColTSA_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resColABAColTSA_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYABAYTSA_tb[,c(1,3)], by=c("gene"="gene")) %>%
  rename(
    LFC_ColABAColGM = log2FoldChange.x,
    LFC_YABAYGM = log2FoldChange.y,
    LFC_ColTSAColGM = log2FoldChange.x.x,
    LFC_YTSAYGM = log2FoldChange.y.y,
    LFC_YGMColGM = log2FoldChange.x.x.x,
    LFC_YABAColABA = log2FoldChange.y.y.y,
    LFC_YTSAColTSA = log2FoldChange.x.x.x.x,
    LFC_ColTSAColABA = log2FoldChange.y.y.y.y,
    LFC_YTSAYABA = log2FoldChange
  )
master_sheet$LFC_YTSAYABA = master_sheet$LFC_YTSAYABA*-1
master_sheet$LFC_ColTSAColABA = master_sheet$LFC_ColTSAColABA*-1
master_sheet$LFC_diff_ABA = master_sheet$LFC_ColABAColGM-master_sheet$LFC_YABAYGM
master_sheet$LFC_diff_TSA = master_sheet$LFC_ColTSAColGM-master_sheet$LFC_YTSAYGM
master_sheet <- left_join(master_sheet, annot, by=c("gene" = "ensgene")) %>% left_join(entrez[,c(2:4)], by=c("gene"="From")) 
master_sheet <- master_sheet[,c(1,19,2:9,17,10:11,18,12:16,20:ncol(master_sheet))]
    


#c(8:10,2:4,14:16,11:13,5:7,17:19)
pheatmap::pheatmap(TSAregulated[,c(3:8)],
                   color = heat_colors,
                   cluster_rows = F,
                   show_rownames = F,
                   border_color = NA,
                   fontsize = 10,
                   scale = "row",
                   fontsize_row = 10,
                   filename = "TSAregulated.png",
                   cluster_cols=FALSE,
                   height = 15)

pheatmap::pheatmap(TSAregulated[,c(18,19)],
                   color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                   cluster_rows = F,
                   show_rownames = F,
                   border_color = NA,
                   fontsize = 10,
                   scale = "none",
                   fontsize_row = 10,
                   filename = "TSAregulatedLFC.png",
                   cluster_cols=FALSE,
                   height = 15,
                   breaks = seq(-5, 5, length.out = 101))


ABATSA = resColABAColGM_tb[,c(1,3)] %>%
  left_join(resYABAYGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resColTSAColGM_tb[,c(1,3)], by=c("gene"="gene")) %>%
  left_join(resYTSAYGM_tb[,c(1,3,2)], by=c("gene"="gene")) %>%
  rename(
    LFC_ColABA = log2FoldChange.x,
    LFC_YABA = log2FoldChange.y,
    LFC_ColTSA = log2FoldChange.x.x,
    LFC_YTSA = log2FoldChange.y.y
  )
ABATSA$abaSign = ABATSA$LFC_ColABA*ABATSA$LFC_YABA
ABATSA$tsaSign = ABATSA$LFC_ColTSA*ABATSA$LFC_YTSA
ABATSA$colSign = ABATSA$LFC_ColTSA*ABATSA$LFC_ColABA
ABATSA$order[ABATSA$LFC_ColABA>0 & ABATSA$LFC_ColTSA>0] <- 1
ABATSA$order[ABATSA$LFC_ColABA<0 & ABATSA$LFC_ColTSA<0] <- 2
ABATSA$order[ABATSA$LFC_ColABA>0 & ABATSA$LFC_ColTSA<0] <- 3
ABATSA$order[ABATSA$LFC_ColABA<0 & ABATSA$LFC_ColTSA>0] <- 4
ABATSAsigs <- ABATSA %>% arrange(order, desc(colSign)) %>%
  filter(gene %in% sigresColABAColGM$gene) %>%
  filter(gene %in% sigresColTSAColGM$gene)
ABATSAsigsmap <- ABATSAsigs %>%
  left_join(count_means, by=c("gene"="gene"))
ABATSAsigsmap <- ABATSAsigsmap[,c(1,11:16)]

TSAvsABA_YnotCol <- master_sheet %>%
  filter(gene %in% sigresYABAYTSA$gene) %>%
  dplyr::filter(!gene %in% sigresColABAColTSA$gene) 
TSAvsABA_YnotCol$sign = TSAvsABA_YnotCol$LFC_YTSAYABA*TSAvsABA_YnotCol$LFC_ColTSAColABA
TSAvsABA_YnotCol$order[TSAvsABA_YnotCol$LFC_YTSAYABA>0] <- 1
TSAvsABA_YnotCol$order[TSAvsABA_YnotCol$LFC_YTSAYABA<0] <- 2
TSAvsABA_YnotCol$delta = TSAvsABA_YnotCol$LFC_YTSAYABA - TSAvsABA_YnotCol$LFC_ColTSAColABA
TSAvsABA_YnotCol <- arrange(TSAvsABA_YnotCol, order, desc(delta))
write.xlsx(TSAvsABA_YnotCol, "TSAvsABA_YnotCol.xlsx")

TSAspecific_AFP2repress <- master_sheet %>%
  filter(gene %in% sigresColTSAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresColABAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresYTSAYGM$gene)
    TSAspecific_AFP2repress$sign = TSAspecific_AFP2repress$LFC_ColABAColGM*TSAspecific_AFP2repress$LFC_ColTSAColGM
    TSAspecific_AFP2repress$order[TSAspecific_AFP2repress$LFC_ColTSAColGM>0] <- 1
    TSAspecific_AFP2repress$order[TSAspecific_AFP2repress$LFC_ColTSAColGM<0] <- 2
    TSAspecific_AFP2repress <- arrange(TSAspecific_AFP2repress, order, desc(abs(LFC_YTSAColTSA)))
    write.xlsx(TSAspecific_AFP2repress, "TSAspecific_AFP2repress.xlsx")

TSAspecific <- master_sheet %>%
  filter(gene %in% sigresColTSAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresColABAColGM$gene)
TSAspecific$sign = TSAspecific$LFC_ColABAColGM*TSAspecific$LFC_ColTSAColGM
TSAspecific$order[TSAspecific$LFC_ColTSAColGM>0] <- 1
TSAspecific$order[TSAspecific$LFC_ColTSAColGM<0] <- 2
TSAspecific <- arrange(TSAspecific, order, desc(abs(LFC_YTSAColTSA)))
write.xlsx(TSAspecific, "TSAspecific.xlsx")

ABAspecific_AFP2repress <- master_sheet %>%
  filter(gene %in% sigresColABAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresColTSAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresYABAYGM$gene)
    ABAspecific_AFP2repress$sign = ABAspecific_AFP2repress$LFC_ColABAColGM*ABAspecific_AFP2repress$LFC_ColTSAColGM
    ABAspecific_AFP2repress$order[ABAspecific_AFP2repress$LFC_ColABAColGM>0] <- 1
    ABAspecific_AFP2repress$order[ABAspecific_AFP2repress$LFC_ColABAColGM<0] <- 2
    ABAspecific_AFP2repress <- arrange(ABAspecific_AFP2repress, order, desc(abs(LFC_YABAColABA)))
    write.xlsx(ABAspecific_AFP2repress, "ABAspecific_AFP2repress.xlsx")
  
ABAspecific <- master_sheet %>%
  filter(gene %in% sigresColABAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresColTSAColGM$gene)
  ABAspecific$sign = ABAspecific$LFC_ColABAColGM*ABAspecific$LFC_ColTSAColGM
  ABAspecific$order[ABAspecific$LFC_ColABAColGM>0] <- 1
  ABAspecific$order[ABAspecific$LFC_ColABAColGM<0] <- 2
  ABAspecific <- arrange(ABAspecific, order, desc(abs(LFC_YABAColABA)))
  write.xlsx(ABAspecific, "ABAspecific.xlsx")
    
  
ABA_AFP2repress <- master_sheet %>%
  filter(gene %in% sigresColABAColGM$gene) %>%
  dplyr::filter(!gene %in% sigresYABAYGM$gene)
    ABA_AFP2repress$order[ABA_AFP2repress$LFC_ColABAColGM>0] <- 1
    ABA_AFP2repress$order[ABA_AFP2repress$LFC_ColABAColGM<0] <- 2
    ABA_AFP2repress <- arrange(ABA_AFP2repress, order, LFC_YTSAYABA)
    write.xlsx(ABA_AFP2repress, "ABA_AFP2repress.xlsx")

TSA_AFP2repress <- master_sheet %>%
      filter(gene %in% sigresColTSAColGM$gene) %>%
      dplyr::filter(!gene %in% sigresYTSAYGM$gene)
    TSA_AFP2repress$order[TSA_AFP2repress$LFC_ColTSAColGM>0] <- 1
    TSA_AFP2repress$order[TSA_AFP2repress$LFC_ColTSAColGM<0] <- 2
    TSA_AFP2repress <- arrange(TSA_AFP2repress, order, LFC_YTSAYABA)
    write.xlsx(TSA_AFP2repress, "TSA_AFP2repress.xlsx")    

ABAregulated <- master_sheet %>%
  filter(gene %in% sigresColABAColGM$gene)
    ABAregulated$order[ABAregulated$LFC_ColABAColGM>0] <- 1
    ABAregulated$order[ABAregulated$LFC_ColABAColGM<0] <- 2
    ABAregulated <- arrange(ABAregulated, order, LFC_diff_ABA)
    write.xlsx(ABAregulated, "ABAregulated.xlsx")

TSAregulated <- master_sheet %>%
    filter(gene %in% sigresColTSAColGM$gene)
    TSAregulated$order[TSAregulated$LFC_ColTSAColGM>0] <- 1
    TSAregulated$order[TSAregulated$LFC_ColTSAColGM<0] <- 2
    TSAregulated <- arrange(TSAregulated, order, LFC_diff_TSA)
    write.xlsx(TSAregulated, "TSAregulated.xlsx")
    
ABA_AFP2repress_byLFC <- slice_max(ABAregulated, abs(LFC_diff_ABA), n = 412) %>% arrange(order, LFC_YTSAYABA)
    write.xlsx(ABA_AFP2repress_byLFC, "ABA_AFP2repress_byLFC.xlsx")


YGMnotYTSA <- master_sheet %>%
  filter(gene %in% sigresYGMColGM$gene) %>%
  dplyr::filter(!gene %in% sigresYTSAColTSA$gene) 
    YGMnotYTSA$order[YGMnotYTSA$LFC_YGMColGM>0] <- 1
    YGMnotYTSA$order[YGMnotYTSA$LFC_YGMColGM<0] <- 2
    YGMnotYTSA <- arrange(YGMnotYTSA, order, LFC_YTSAColTSA)
    write.xlsx(YGMnotYTSA, "YGMnotYTSA.xlsx")
    
ABA_YvsCol <- master_sheet %>%
  filter(gene %in% sigresYABAColABA$gene)
ABA_YvsCol$order[ABA_YvsCol$LFC_YABAColABA>0] <- 1
ABA_YvsCol$order[ABA_YvsCol$LFC_YABAColABA<0] <- 2
ABA_YvsCol <- arrange(ABA_YvsCol, order, LFC_YABAColABA)
write.xlsx(ABA_YvsCol, "ABA_YvsCol.xlsx")

TSA_YvsCol <- master_sheet %>%
  filter(gene %in% sigresYTSAColTSA$gene)
TSA_YvsCol$order[TSA_YvsCol$LFC_YTSAColTSA>0] <- 1
TSA_YvsCol$order[TSA_YvsCol$LFC_YTSAColTSA<0] <- 2
TSA_YvsCol <- arrange(TSA_YvsCol, order, LFC_YTSAColTSA)
write.xlsx(TSA_YvsCol, "TSA_YvsCol.xlsx")

coexpress_AFP2_TSA <- master_sheet %>%
  filter(gene %in% sigresColTSAColGM$gene) %>%
  filter(gene %in% sigresYTSAColTSA$gene)
coexpress_AFP2_TSA$sign = coexpress_AFP2_TSA$LFC_ColTSAColGM * coexpress_AFP2_TSA$LFC_YTSAColTSA

coexpress_ABATSA_down <- master_sheet %>%
  filter(gene %in% sigresColABAColGM$gene[sigresColABAColGM$log2FoldChange<0]) %>%
  filter(gene %in% sigresColTSAColGM$gene[sigresColTSAColGM$log2FoldChange<0])

coexpress_ABATSA_up <- master_sheet %>%
  filter(gene %in% sigresColABAColGM$gene[sigresColABAColGM$log2FoldChange>0]) %>%
  filter(gene %in% sigresColTSAColGM$gene[sigresColTSAColGM$log2FoldChange>0])

NCEDgenes <- c("AT1G30100", "AT1G78390", "AT2G44990", "AT3G14440", "AT3G24220", "AT3G63520", "AT4G18350")
NCEDlfcs <- all_LFCs %>% filter(gene %in% NCEDgenes)
NCEDcounts <- count_means %>% filter(gene %in% NCEDgenes) %>% left_join(annot, by=c("gene" = "ensgene")) 

GMTfile <- read_tsv("data/araslim2.gmt")
GMTlabels <- GMTfile[,c(1)]
GMTdupl <- GMTlabels %>% filter(duplicated(.[["HOEWYK_SE_ROOT_DN"]]))

ranksYABAColABA <- master_sheet[c("gene", "LFC_YABAColABA")] %>% arrange(desc(LFC_YABAColABA))
ranksYABAColABA <- setNames(ranksYABAColABA$LFC_YABAColABA, ranksYABAColABA$gene)
ranksYTSAColTSA <- master_sheet[c("gene", "LFC_YTSAColTSA")] %>% arrange(desc(LFC_YTSAColTSA))
ranksYTSAColTSA <- setNames(ranksYTSAColTSA$LFC_YTSAColTSA, ranksYTSAColTSA$gene)
ranksYGMColGM <- master_sheet[c("gene", "LFC_YGMColGM")] %>% arrange(desc(LFC_YGMColGM))
ranksYGMColGM <- setNames(ranksYGMColGM$LFC_YGMColGM, ranksYGMColGM$gene)

ranksColABAColGM <- master_sheet[c("gene", "LFC_ColABAColGM")] %>% arrange(desc(LFC_ColABAColGM))
ranksColABAColGM <- setNames(ranksColABAColGM$LFC_ColABAColGM, ranksColABAColGM$gene)
ranksYABAYGM <- master_sheet[c("gene", "LFC_YABAYGM")] %>% arrange(desc(LFC_YABAYGM))
ranksYABAYGM <- setNames(ranksYABAYGM$LFC_YABAYGM, ranksYABAYGM$gene)

ranksYTSAYGM <- master_sheet[c("gene", "LFC_YTSAYGM")] %>% arrange(desc(LFC_YTSAYGM))
ranksYTSAYGM <- setNames(ranksYTSAYGM$LFC_YTSAYGM, ranksYTSAYGM$gene)
ranksColTSAColGM <- master_sheet[c("gene", "LFC_ColTSAColGM")] %>% arrange(desc(LFC_ColTSAColGM))
ranksColTSAColGM <- setNames(ranksColTSAColGM$LFC_ColTSAColGM, ranksColTSAColGM$gene)

ranksYTSAYABA <- master_sheet[c("gene", "LFC_YTSAYABA")] %>% arrange(desc(LFC_YTSAYABA))
ranksYTSAYABA <- setNames(ranksYTSAYABA$LFC_YTSAYABA, ranksYTSAYABA$gene)
ranksColTSAColABA <- master_sheet[c("gene", "LFC_ColTSAColABA")] %>% arrange(desc(LFC_ColTSAColABA))
ranksColTSAColABA <- setNames(ranksColTSAColABA$LFC_ColTSAColABA, ranksColTSAColABA$gene)
###
###
ranksABSYABAColABA <- master_sheet[c("gene", "LFC_YABAColABA")] 
ranksABSYABAColABA$LFC_YABAColABA <- abs(ranksABSYABAColABA$LFC_YABAColABA)
ranksABSYABAColABA <- arrange(ranksABSYABAColABA, desc(LFC_YABAColABA))
ranksABSYABAColABA <- setNames(ranksABSYABAColABA$LFC_YABAColABA, ranksABSYABAColABA$gene)

ranksABSYTSAColTSA <- master_sheet[c("gene", "LFC_YTSAColTSA")] 
ranksABSYTSAColTSA$LFC_YTSAColTSA <- abs(ranksABSYTSAColTSA$LFC_YTSAColTSA)
ranksABSYTSAColTSA <- arrange(ranksABSYTSAColTSA, desc(LFC_YTSAColTSA))
ranksABSYTSAColTSA <- setNames(ranksABSYTSAColTSA$LFC_YTSAColTSA, ranksABSYTSAColTSA$gene)

ranksABSYGMColGM <- master_sheet[c("gene", "LFC_YGMColGM")] 
ranksABSYGMColGM$LFC_YGMColGM <- abs(ranksABSYGMColGM$LFC_YGMColGM)
ranksABSYGMColGM <- arrange(ranksABSYGMColGM, desc(LFC_YGMColGM))
ranksABSYGMColGM <- setNames(ranksABSYGMColGM$LFC_YGMColGM, ranksABSYGMColGM$gene)
##
ranksABSColABAColGM <- master_sheet[c("gene", "LFC_ColABAColGM")] 
ranksABSColABAColGM$LFC_ColABAColGM <- abs(ranksABSColABAColGM$LFC_ColABAColGM)
ranksABSColABAColGM <- arrange(ranksABSColABAColGM, desc(LFC_ColABAColGM))
ranksABSColABAColGM <- setNames(ranksABSColABAColGM$LFC_ColABAColGM, ranksABSColABAColGM$gene)

ranksABSYABAYGM <- master_sheet[c("gene", "LFC_YABAYGM")] 
ranksABSYABAYGM$LFC_YABAYGM <- abs(ranksABSYABAYGM$LFC_YABAYGM)
ranksABSYABAYGM <- arrange(ranksABSYABAYGM, desc(LFC_YABAYGM))
ranksABSYABAYGM <- setNames(ranksABSYABAYGM$LFC_YABAYGM, ranksABSYABAYGM$gene)
##
ranksABSYTSAYGM <- master_sheet[c("gene", "LFC_YTSAYGM")] 
ranksABSYTSAYGM$LFC_YTSAYGM <- abs(ranksABSYTSAYGM$LFC_YTSAYGM)
ranksABSYTSAYGM <- arrange(ranksABSYTSAYGM, desc(LFC_YTSAYGM))
ranksABSYTSAYGM <- setNames(ranksABSYTSAYGM$LFC_YTSAYGM, ranksABSYTSAYGM$gene)

ranksABSColTSAColGM <- master_sheet[c("gene", "LFC_ColTSAColGM")] 
ranksABSColTSAColGM$LFC_ColTSAColGM <- abs(ranksABSColTSAColGM$LFC_ColTSAColGM)
ranksABSColTSAColGM <- arrange(ranksABSColTSAColGM, desc(LFC_ColTSAColGM))
ranksABSColTSAColGM <- setNames(ranksABSColTSAColGM$LFC_ColTSAColGM, ranksABSColTSAColGM$gene)
##
ranksABSYTSAYABA <- master_sheet[c("gene", "LFC_YTSAYABA")] 
ranksABSYTSAYABA$LFC_YTSAYABA <- abs(ranksABSYTSAYABA$LFC_YTSAYABA)
ranksABSYTSAYABA <- arrange(ranksABSYTSAYABA, desc(LFC_YTSAYABA))
ranksABSYTSAYABA <- setNames(ranksABSYTSAYABA$LFC_YTSAYABA, ranksABSYTSAYABA$gene)

ranksABSColTSAColABA <- master_sheet[c("gene", "LFC_ColTSAColABA")] 
ranksABSColTSAColABA$LFC_ColTSAColABA <- abs(ranksABSColTSAColABA$LFC_ColTSAColABA)
ranksABSColTSAColABA <- arrange(ranksABSColTSAColABA, desc(LFC_ColTSAColABA))
ranksABSColTSAColABA <- setNames(ranksABSColTSAColABA$LFC_ColTSAColABA, ranksABSColTSAColABA$gene)
##

len <- sapply(pathways_GO, length)
testtail <- tail(pathways_GO[order(len)], 100)

splitcommas <- function(input) {
  output = unlist(strsplit(input, split=","))
  return(output)
}
splitcommas2 <- function(input) {
  output = unlist(strsplit(input, split=", "))
  return(output)
}

slim_bp <- read.delim("data/mapped_slim_bp.txt")
slim_mf <- read.delim("data/mapped_slim_mf.txt")
slims <- full_join(slim_bp, slim_mf)
slims <- slims[,c(2,1,9)]
write_tsv(slims, "slimsgost.tsv")

pathwayLines <- strsplit(readLines("data/slims.gmt"), "\t")
pathways_GOslim <- lapply(pathwayLines, tail, -1)
names(pathways_GOslim) <- sapply(pathwayLines, head, 1)
pathways_GOslim <- lapply(pathways_GOslim, splitcommas2)

pathways_GO <- gmtPathways("data/Ara_GO.gmt") %>% lapply(splitcommas)
pathways_GO <- pathways_GO[order(sapply(pathways_GO, length), decreasing=TRUE)]
pathways_lit <- gmtPathways("data/ara_literature.gmt") %>% lapply(splitcommas)
pathways_lit <- pathways_lit[order(sapply(pathways_lit, length), decreasing=TRUE)]

GSEA_GO_YABAColABA <- fgsea(pathways_GO, ranksYABAColABA, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
GSEA_GO_YTSAColTSA <- fgsea(pathways_GO, ranksYTSAColTSA, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
GSEA_GO_YGMColGM <- fgsea(pathways_GO, ranksYGMColGM, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
GSEA_GO_ColABAColGM <- fgsea(pathways_GO, ranksColABAColGM, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
collapsedPathways <- collapsePathways(GSEA_GO_ColABAColGM[order(pval)][padj < 0.01],
                                      pathways_GO, ranksColABAColGM)
mainPathways <- GSEA_GO_ColABAColGM[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
GSEA_GO_ColABAColGMcollapsed <- GSEA_GO_ColABAColGM[pathway %in% mainPathways]
GSEA_GO_YABAYGM <- fgsea(pathways_GO, ranksYABAYGM, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
GSEA_GO_YTSAYGM <- fgsea(pathways_GO, ranksYTSAYGM, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
GSEA_GO_ColTSAColGM <- fgsea(pathways_GO, ranksColTSAColGM, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
collapsedPathways <- collapsePathways(GSEA_GO_ColTSAColGM[order(pval)][padj < 0.01],
                                      pathways_GO, ranksColTSAColGM)
mainPathways <- GSEA_GO_ColTSAColGM[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
GSEA_GO_ColTSAColGMcollapsed <- GSEA_GO_ColTSAColGM[pathway %in% mainPathways]
GSEA_GO_YTSAYABA <- fgsea(pathways_GO, ranksYTSAYABA, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)
GSEA_GO_ColTSAColABA <- fgsea(pathways_GO, ranksColTSAColABA, eps=0.0, minSize=8, maxSize=1660, nPermSimple = 100000)

GSEA_LIT_YABAColABA <- fgsea(pathways_lit, ranksYABAColABA, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)
GSEA_LIT_YTSAColTSA <- fgsea(pathways_lit, ranksYTSAColTSA, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)
GSEA_LIT_YGMColGM <- fgsea(pathways_lit, ranksYGMColGM, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)
GSEA_LIT_ColABAColGM <- fgsea(pathways_lit, ranksColABAColGM, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)
GSEA_LIT_YABAYGM <- fgsea(pathways_lit, ranksYABAYGM, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)
GSEA_LIT_YTSAYGM <- fgsea(pathways_lit, ranksYTSAYGM, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)
GSEA_LIT_ColTSAColGM <- fgsea(pathways_lit, ranksColTSAColGM, eps=0.0, minSize=8, maxSize=4000, nPermSimple = 1000)

GSEA_GOslim_YABAColABA <- fgsea(pathways_GOslim, ranksYABAColABA)
GSEA_GOslim_YTSAColTSA <- fgsea(pathways_GOslim, ranksYTSAColTSA)
GSEA_GOslim_YGMColGM <- fgsea(pathways_GOslim, ranksYGMColGM)
GSEA_GOslim_ColABAColGM <- fgsea(pathways_GOslim, ranksColABAColGM)
GSEA_GOslim_YABAYGM <- fgsea(pathways_GOslim, ranksYABAYGM)
GSEA_GOslim_YTSAYGM <- fgsea(pathways_GOslim, ranksYTSAYGM)
GSEA_GOslim_ColTSAColGM <- fgsea(pathways_GOslim, ranksColTSAColGM)
GSEA_GOslim_YTSAYABA <- fgsea(pathways_GOslim, ranksYTSAYABA)
GSEA_GOslim_ColTSAColABA <- fgsea(pathways_GOslim, ranksColTSAColABA)

GSEA_GOslimABS_YABAColABA <- fgsea(pathways_GOslim, ranksABSYABAColABA, scoreType="pos")
GSEA_GOslimABS_YTSAColTSA <- fgsea(pathways_GOslim, ranksABSYTSAColTSA, scoreType="pos")
GSEA_GOslimABS_YGMColGM <- fgsea(pathways_GOslim, ranksABSYGMColGM, scoreType="pos")
GSEA_GOslimABS_ColABAColGM <- fgsea(pathways_GOslim, ranksABSColABAColGM, scoreType="pos")
GSEA_GOslimABS_YABAYGM <- fgsea(pathways_GOslim, ranksABSYABAYGM, scoreType="pos")
GSEA_GOslimABS_YTSAYGM <- fgsea(pathways_GOslim, ranksABSYTSAYGM, scoreType="pos")
GSEA_GOslimABS_ColTSAColGM <- fgsea(pathways_GOslim, ranksABSColTSAColGM, scoreType="pos")
GSEA_GOslimABS_YTSAYABA <- fgsea(pathways_GOslim, ranksABSYTSAYABA, scoreType="pos")
GSEA_GOslimABS_ColTSAColABA <- fgsea(pathways_GOslim, ranksABSColTSAColABA, scoreType="pos")
save.image(file='GSEAgoslim.RData')

fwrite(GSEA_GO_YABAColABA, file="GSEA_GO_YABAColABA.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_YTSAColTSA, file="GSEA_GO_YTSAColTSA.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_YGMColGM, file="GSEA_GO_YGMColGM.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_ColABAColGM, file="GSEA_GO_ColABAColGM.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_YABAYGM, file="GSEA_GO_YABAYGM.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_YTSAYGM, file="GSEA_GO_YTSAYGM.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_ColTSAColGM, file="GSEA_GO_ColTSAColGM.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_YTSAYABA, file="GSEA_GO_YTSAYABA.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(GSEA_GO_ColTSAColABA, file="GSEA_GO_CColTSAColABA.tsv", sep="\t", sep2=c("", " ", ""))

figure_test2 <- plotEnrichment(pathways_lit[["ZHANG_SEED_DIFF"]], ranksColABAColGM)
figure_test2

figure_test3 <- plotEnrichment(pathways_lit[["DIJK_DEHYDRATION-STRESS_DN"]], ranksColTSAColGM)
figure_test3

figure_GSEA_GO_SEED_MATURATION <- ggarrange(plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksYABAColABA) + labs(x = paste("YABAColABA ES=",substr(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))), 
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksYTSAColTSA) + labs(x = paste("YTSAColTSA ES=",substr(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),  
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksYGMColGM) + labs(x = paste("YGMColGM ES=",substr(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),  
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksColABAColGM) + labs(x = paste("ColABAColGM ES=",substr(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),  
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksYABAYGM) + labs(x = paste("YABAYGM ES=",substr(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))), 
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksColTSAColGM) + labs(x = paste("ColTSAColGM ES=",substr(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksYTSAYGM) + labs(x = paste("YTSAYGM ES=",substr(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),  
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksYTSAYABA) + labs(x = paste("YTSAYABA ES=",substr(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),  
                            plotEnrichment(pathways_GO[["SEED_MATURATION"]], ranksColTSAColABA) + labs(x = paste("ColTSAColABA ES=",substr(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="SEED_MATURATION", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="SEED_MATURATION", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="SEED_MATURATION", arr.ind=TRUE),3]),format="e",digits=2))),  
                            ncol = 1, nrow = 9)
pdf("GSEA_GO_SEED_MATURATION.pdf", width = 6, height =11)
figure_GSEA_GO_SEED_MATURATION
dev.off()

figure_GSEA_GO_RESPONSE_TO_ABSCISIC_ACID_STIMULUS <- ggarrange(plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksYABAColABA) + labs(x = paste("YABAColABA ES=",substr(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))), 
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksYTSAColTSA) + labs(x = paste("YTSAColTSA ES=",substr(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),  
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksYGMColGM) + labs(x = paste("YGMColGM ES=",substr(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),  
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksColABAColGM) + labs(x = paste("ColABAColGM ES=",substr(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),  
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksYABAYGM) + labs(x = paste("YABAYGM ES=",substr(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))), 
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksColTSAColGM) + labs(x = paste("ColTSAColGM ES=",substr(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksYTSAYGM) + labs(x = paste("YTSAYGM ES=",substr(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),  
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksYTSAYABA) + labs(x = paste("YTSAYABA ES=",substr(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),  
                           plotEnrichment(pathways_GO[["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"]], ranksColTSAColABA) + labs(x = paste("ColTSAColABA ES=",substr(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="RESPONSE_TO_ABSCISIC_ACID_STIMULUS", arr.ind=TRUE),3]),format="e",digits=2))),                             
                           ncol = 1, nrow = 9)
pdf("GSEA_GO_RESPONSE_TO_ABSCISIC_ACID_STIMULUS.pdf", width = 6, height =11)
figure_GSEA_GO_RESPONSE_TO_ABSCISIC_ACID_STIMULUS
dev.off()

figure_GSEA_GO_ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY <- ggarrange(plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksYABAColABA) + labs(x = paste("YABAColABA ES=",substr(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YABAColABA[ which(GSEA_GO_YABAColABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))), 
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksYTSAColTSA) + labs(x = paste("YTSAColTSA ES=",substr(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAColTSA[ which(GSEA_GO_YTSAColTSA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),  
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksYGMColGM) + labs(x = paste("YGMColGM ES=",substr(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YGMColGM[ which(GSEA_GO_YGMColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),  
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksColABAColGM) + labs(x = paste("ColABAColGM ES=",substr(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColABAColGM[ which(GSEA_GO_ColABAColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),  
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksYABAYGM) + labs(x = paste("YABAYGM ES=",substr(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YABAYGM[ which(GSEA_GO_YABAYGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))), 
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksColTSAColGM) + labs(x = paste("ColTSAColGM ES=",substr(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColTSAColGM[ which(GSEA_GO_ColTSAColGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksYTSAYGM) + labs(x = paste("YTSAYGM ES=",substr(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAYGM[ which(GSEA_GO_YTSAYGM$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),  
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksYTSAYABA) + labs(x = paste("YTSAYABA ES=",substr(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_YTSAYABA[ which(GSEA_GO_YTSAYABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),  
                                                               plotEnrichment(pathways_GO[["ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY"]], ranksColTSAColABA) + labs(x = paste("ColTSAColABA ES=",substr(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),5],1,6)," NES=",substr(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),6],1,6),", padj=", formatC(unlist(GSEA_GO_ColTSAColABA[ which(GSEA_GO_ColTSAColABA$pathway=="ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY", arr.ind=TRUE),3]),format="e",digits=2))),                                                                 
                                                               ncol = 1, nrow = 9)
pdf("GSEA_GO_ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY.pdf", width = 6, height =11)
figure_GSEA_GO_ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY
dev.off()

plotGseaTable(pathways_GO["RESPONSE_TO_ABSCISIC_ACID_STIMULUS"], ranksYABAColABA, GSEA_GO_YABAColABA, gseaParam=0.5)

topPathways <- GSEA_GO_YABAColABA %>% 
  top_n(20, wt=abs(NES)) %>% 
  arrange(-NES) %>% 
  pull(pathway)

pdf("topPathways.pdf", width = 12, height =9)
plotGseaTable(pathways_GO[topPathways], 
              ranksYABAColABA, 
              GSEA_GO_YABAColABA, 
              gseaParam = 0.5)
dev.off()

upload_GMT_file(gmtfile = "slimsgost.gmt")

gost_ABA <- gost(query = ABAregulated$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABAup <- gost(query = ABAregulated$gene[ABAregulated$LFC_ColABAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABAdown <- gost(query = ABAregulated$gene[ABAregulated$LFC_ColABAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_ABA_AFP2repress <- gost(query = ABA_AFP2repress$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABA_AFP2repressup <- gost(query = ABA_AFP2repress$gene[ABA_AFP2repress$LFC_ColABAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABA_AFP2repressdown <- gost(query = ABA_AFP2repress$gene[ABA_AFP2repress$LFC_ColABAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_TSA <- gost(query = TSAregulated$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAup <- gost(query = TSAregulated$gene[TSAregulated$LFC_ColTSAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAdown <- gost(query = TSAregulated$gene[TSAregulated$LFC_ColTSAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_TSA_AFP2repress <- gost(query = TSA_AFP2repress$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSA_AFP2repressup <- gost(query = TSA_AFP2repress$gene[TSA_AFP2repress$LFC_ColTSAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSA_AFP2repressdown <- gost(query = TSA_AFP2repress$gene[TSA_AFP2repress$LFC_ColTSAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_ABAspecific <- gost(query = ABAspecific$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABAspecificup <- gost(query = ABAspecific$gene[ABAspecific$LFC_ColABAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABAspecificdown <- gost(query = ABAspecific$gene[ABAspecific$LFC_ColABAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_ABAspecific_AFP2repress <- gost(query = ABAspecific_AFP2repress$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABAspecific_AFP2repressup <- gost(query = ABAspecific_AFP2repress$gene[ABAspecific_AFP2repress$LFC_ColABAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABAspecific_AFP2repressdown <- gost(query = ABAspecific_AFP2repress$gene[ABAspecific_AFP2repress$LFC_ColABAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_TSAspecific <- gost(query = TSAspecific$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAspecificup <- gost(query = TSAspecific$gene[TSAspecific$LFC_ColTSAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAspecificdown <- gost(query = TSAspecific$gene[TSAspecific$LFC_ColTSAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_TSAspecific_AFP2repress <- gost(query = TSAspecific_AFP2repress$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAspecific_AFP2repressup <- gost(query = TSAspecific_AFP2repress$gene[TSAspecific_AFP2repress$LFC_ColTSAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAspecific_AFP2repressdown <- gost(query = TSAspecific_AFP2repress$gene[TSAspecific_AFP2repress$LFC_ColTSAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

#####

gost_TSAvsABA_YnotCol <- gost(query = TSAvsABA_YnotCol$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAvsABA_YnotColup <- gost(query = TSAvsABA_YnotCol$gene[TSAvsABA_YnotCol$LFC_YTSAYABA>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSAvsABA_YnotColdown <- gost(query = TSAvsABA_YnotCol$gene[TSAvsABA_YnotCol$LFC_YTSAYABA<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_YGMnotYTSA <- gost(query = YGMnotYTSA$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_YGMnotYTSA <- gost(query = YGMnotYTSA$gene[YGMnotYTSA$LFC_YGMColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_YGMnotYTSA <- gost(query = YGMnotYTSA$gene[YGMnotYTSA$LFC_YGMColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

#######

gost_ABA_AFP2repress_byLFC <- gost(query = ABA_AFP2repress_byLFC$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABA_AFP2repress_byLFCup <- gost(query = ABA_AFP2repress_byLFC$gene[ABA_AFP2repress_byLFC$LFC_ColABAColGM>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABA_AFP2repress_byLFCdown <- gost(query = ABA_AFP2repress_byLFC$gene[ABA_AFP2repress_byLFC$LFC_ColABAColGM<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_ABA_YvsCol <- gost(query = ABA_YvsCol$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABA_YvsColup <- gost(query = ABA_YvsCol$gene[ABA_YvsCol$LFC_YABAColABA>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_ABA_YvsColdown <- gost(query = ABA_YvsCol$gene[ABA_YvsCol$LFC_YABAColABA<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

gost_TSA_YvsCol <- gost(query = TSA_YvsCol$gene, organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSA_YvsColup <- gost(query = TSA_YvsCol$gene[TSA_YvsCol$LFC_YTSAColTSA>0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")
gost_TSA_YvsColdown <- gost(query = TSA_YvsCol$gene[TSA_YvsCol$LFC_YTSAColTSA<0], organism = "gp__CE2a_FN9H_GFM", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "bonferroni")

GO_ABAvsTSA_YvsCol <- gost_ABA_YvsCol[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA_YvsCol[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    intersection_ABA = intersection.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
    intersection_TSA = intersection.y
  )
GO_ABAvsTSA_YvsCol$order[GO_ABAvsTSA_YvsCol$significant_ABA & GO_ABAvsTSA_YvsCol$significant_TSA] <- 1
GO_ABAvsTSA_YvsCol$order[GO_ABAvsTSA_YvsCol$significant_ABA & !GO_ABAvsTSA_YvsCol$significant_TSA] <- 2
GO_ABAvsTSA_YvsCol$order[!GO_ABAvsTSA_YvsCol$significant_ABA & GO_ABAvsTSA_YvsCol$significant_TSA] <- 3
GO_ABAvsTSA_YvsCol$sign = GO_ABAvsTSA_YvsCol$p_value_ABA*GO_ABAvsTSA_YvsCol$p_value_TSA
GO_ABAvsTSA_YvsCol <- arrange(GO_ABAvsTSA_YvsCol, order, sign)
GO_ABAvsTSA_YvsCol <- GO_ABAvsTSA_YvsCol[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_YvsCol,"GO_ABAvsTSA_YvsCol.xlsx")

GO_ABAvsTSA_YvsCol_up <- gost_ABA_YvsColup[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA_YvsColup[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    intersection_ABA = intersection.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
    intersection_TSA = intersection.y
  )
GO_ABAvsTSA_YvsCol_up$order[GO_ABAvsTSA_YvsCol_up$significant_ABA & GO_ABAvsTSA_YvsCol_up$significant_TSA] <- 1
GO_ABAvsTSA_YvsCol_up$order[GO_ABAvsTSA_YvsCol_up$significant_ABA & !GO_ABAvsTSA_YvsCol_up$significant_TSA] <- 2
GO_ABAvsTSA_YvsCol_up$order[!GO_ABAvsTSA_YvsCol_up$significant_ABA & GO_ABAvsTSA_YvsCol_up$significant_TSA] <- 3
GO_ABAvsTSA_YvsCol_up$sign = GO_ABAvsTSA_YvsCol_up$p_value_ABA*GO_ABAvsTSA_YvsCol_up$p_value_TSA
GO_ABAvsTSA_YvsCol_up <- arrange(GO_ABAvsTSA_YvsCol_up, order, sign)
GO_ABAvsTSA_YvsCol_up <- GO_ABAvsTSA_YvsCol_up[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_YvsCol_up,"GO_ABAvsTSA_YvsCol_up.xlsx")

GO_ABAvsTSA_YvsCol_down <- gost_ABA_YvsColdown[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA_YvsColdown[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    intersection_ABA = intersection.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
    intersection_TSA = intersection.y
  )
GO_ABAvsTSA_YvsCol_down$order[GO_ABAvsTSA_YvsCol_down$significant_ABA & GO_ABAvsTSA_YvsCol_down$significant_TSA] <- 1
GO_ABAvsTSA_YvsCol_down$order[GO_ABAvsTSA_YvsCol_down$significant_ABA & !GO_ABAvsTSA_YvsCol_down$significant_TSA] <- 2
GO_ABAvsTSA_YvsCol_down$order[!GO_ABAvsTSA_YvsCol_down$significant_ABA & GO_ABAvsTSA_YvsCol_down$significant_TSA] <- 3
GO_ABAvsTSA_YvsCol_down$sign = GO_ABAvsTSA_YvsCol_down$p_value_ABA*GO_ABAvsTSA_YvsCol_down$p_value_TSA
GO_ABAvsTSA_YvsCol_down <- arrange(GO_ABAvsTSA_YvsCol_down, order, sign)
GO_ABAvsTSA_YvsCol_down <- GO_ABAvsTSA_YvsCol_down[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_YvsCol_down,"GO_ABAvsTSA_YvsCol_down.xlsx")

GO_ABAvsTSA <- gost_ABA[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    intersection_ABA = intersection.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
    intersection_TSA = intersection.y
  )
GO_ABAvsTSA$order[GO_ABAvsTSA$significant_ABA & GO_ABAvsTSA$significant_TSA] <- 1
GO_ABAvsTSA$order[GO_ABAvsTSA$significant_ABA & !GO_ABAvsTSA$significant_TSA] <- 2
GO_ABAvsTSA$order[!GO_ABAvsTSA$significant_ABA & GO_ABAvsTSA$significant_TSA] <- 3
GO_ABAvsTSA$sign = GO_ABAvsTSA$p_value_ABA*GO_ABAvsTSA$p_value_TSA
GO_ABAvsTSA <- arrange(GO_ABAvsTSA, order, sign)
GO_ABAvsTSA <- GO_ABAvsTSA[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA,"GO_ABAvsTSA.xlsx")

GO_ABAvsTSA_down <- gost_ABAdown[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSAdown[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_down$order[GO_ABAvsTSA_down$significant_ABA & GO_ABAvsTSA_down$significant_TSA] <- 1
GO_ABAvsTSA_down$order[GO_ABAvsTSA_down$significant_ABA & !GO_ABAvsTSA_down$significant_TSA] <- 2
GO_ABAvsTSA_down$order[!GO_ABAvsTSA_down$significant_ABA & GO_ABAvsTSA_down$significant_TSA] <- 3
GO_ABAvsTSA_down$sign = GO_ABAvsTSA_down$p_value_ABA*GO_ABAvsTSA_down$p_value_TSA
GO_ABAvsTSA_down <- arrange(GO_ABAvsTSA_down, order, sign)
GO_ABAvsTSA_down <- GO_ABAvsTSA_down[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_down, "GO_ABAvsTSAdown.xlsx")

GO_ABAvsTSA_up <- gost_ABAup[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSAup[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_up$order[GO_ABAvsTSA_up$significant_ABA & GO_ABAvsTSA_up$significant_TSA] <- 1
GO_ABAvsTSA_up$order[GO_ABAvsTSA_up$significant_ABA & !GO_ABAvsTSA_up$significant_TSA] <- 2
GO_ABAvsTSA_up$order[!GO_ABAvsTSA_up$significant_ABA & GO_ABAvsTSA_up$significant_TSA] <- 3
GO_ABAvsTSA_up$sign = GO_ABAvsTSA_up$p_value_ABA*GO_ABAvsTSA_up$p_value_TSA
GO_ABAvsTSA_up <- arrange(GO_ABAvsTSA_up, order, sign)
GO_ABAvsTSA_up <- GO_ABAvsTSA_up[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_up, "GO_ABAvsTSAup.xlsx")

GO_ABAvsTSA_specific <- gost_ABAspecific[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSAspecific[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_specific$order[GO_ABAvsTSA_specific$significant_ABA & GO_ABAvsTSA_specific$significant_TSA] <- 1
GO_ABAvsTSA_specific$order[GO_ABAvsTSA_specific$significant_ABA & !GO_ABAvsTSA_specific$significant_TSA] <- 2
GO_ABAvsTSA_specific$order[!GO_ABAvsTSA_specific$significant_ABA & GO_ABAvsTSA_specific$significant_TSA] <- 3
GO_ABAvsTSA_specific$sign = GO_ABAvsTSA_specific$p_value_ABA*GO_ABAvsTSA_specific$p_value_TSA
GO_ABAvsTSA_specific <- arrange(GO_ABAvsTSA_specific, order, sign)
GO_ABAvsTSA_specific <- GO_ABAvsTSA_specific[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]

GO_ABAvsTSA_specific_up <- gost_ABAspecificup[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSAspecificup[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_specific_up$order[GO_ABAvsTSA_specific_up$significant_ABA & GO_ABAvsTSA_specific_up$significant_TSA] <- 1
GO_ABAvsTSA_specific_up$order[GO_ABAvsTSA_specific_up$significant_ABA & !GO_ABAvsTSA_specific_up$significant_TSA] <- 2
GO_ABAvsTSA_specific_up$order[!GO_ABAvsTSA_specific_up$significant_ABA & GO_ABAvsTSA_specific_up$significant_TSA] <- 3
GO_ABAvsTSA_specific_up$sign = GO_ABAvsTSA_specific_up$p_value_ABA*GO_ABAvsTSA_specific_up$p_value_TSA
GO_ABAvsTSA_specific_up <- arrange(GO_ABAvsTSA_specific_up, order, sign)
GO_ABAvsTSA_specific_up <- GO_ABAvsTSA_specific_up[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_specific_up,"GO_ABAvsTSA_specific_up.xlsx")

GO_ABAvsTSA_AFP2repress <- gost_ABA_AFP2repress[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA_AFP2repress[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_AFP2repress$order[GO_ABAvsTSA_AFP2repress$significant_ABA & GO_ABAvsTSA_AFP2repress$significant_TSA] <- 1
GO_ABAvsTSA_AFP2repress$order[GO_ABAvsTSA_AFP2repress$significant_ABA & !GO_ABAvsTSA_AFP2repress$significant_TSA] <- 2
GO_ABAvsTSA_AFP2repress$order[!GO_ABAvsTSA_AFP2repress$significant_ABA & GO_ABAvsTSA_AFP2repress$significant_TSA] <- 3
GO_ABAvsTSA_AFP2repress$sign = GO_ABAvsTSA_AFP2repress$p_value_ABA*GO_ABAvsTSA_AFP2repress$p_value_TSA
GO_ABAvsTSA_AFP2repress <- arrange(GO_ABAvsTSA_AFP2repress, order, sign)
GO_ABAvsTSA_AFP2repress <- GO_ABAvsTSA_AFP2repress[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_AFP2repress,"GO_ABAvsTSA_AFP2repress.xlsx")


GO_ABAvsTSA_AFP2repress_down <- gost_ABA_AFP2repressdown[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA_AFP2repressdown[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_AFP2repress_down$order[GO_ABAvsTSA_AFP2repress_down$significant_ABA & GO_ABAvsTSA_AFP2repress_down$significant_TSA] <- 1
GO_ABAvsTSA_AFP2repress_down$order[GO_ABAvsTSA_AFP2repress_down$significant_ABA & !GO_ABAvsTSA_AFP2repress_down$significant_TSA] <- 2
GO_ABAvsTSA_AFP2repress_down$order[!GO_ABAvsTSA_AFP2repress_down$significant_ABA & GO_ABAvsTSA_AFP2repress_down$significant_TSA] <- 3
GO_ABAvsTSA_AFP2repress_down$sign = GO_ABAvsTSA_AFP2repress_down$p_value_ABA*GO_ABAvsTSA_AFP2repress_down$p_value_TSA
GO_ABAvsTSA_AFP2repress_down <- arrange(GO_ABAvsTSA_AFP2repress_down, order, sign)
GO_ABAvsTSA_AFP2repress_down <- GO_ABAvsTSA_AFP2repress_down[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_AFP2repress_down,"GO_ABAvsTSA_AFP2repress_down.xlsx")

GO_ABAvsTSA_AFP2repress_up <- gost_ABA_AFP2repressup[["result"]][,c(9,10,11,4,5,6,16,3,2)] %>%
  left_join(gost_TSA_AFP2repressup[["result"]][,c(9,5,6,16,3,2)], by=c("term_id" = "term_id")) %>%
  rename(
    query_size_ABA = query_size.x,
    intersection_size_ABA = intersection_size.x,
    p_value_ABA = p_value.x,
    significant_ABA = significant.x,
    query_size_TSA = query_size.y,
    intersection_size_TSA = intersection_size.y,
    p_value_TSA = p_value.y,
    significant_TSA = significant.y,
  )
GO_ABAvsTSA_AFP2repress_up$order[GO_ABAvsTSA_AFP2repress_up$significant_ABA & GO_ABAvsTSA_AFP2repress_up$significant_TSA] <- 1
GO_ABAvsTSA_AFP2repress_up$order[GO_ABAvsTSA_AFP2repress_up$significant_ABA & !GO_ABAvsTSA_AFP2repress_up$significant_TSA] <- 2
GO_ABAvsTSA_AFP2repress_up$order[!GO_ABAvsTSA_AFP2repress_up$significant_ABA & GO_ABAvsTSA_AFP2repress_up$significant_TSA] <- 3
GO_ABAvsTSA_AFP2repress_up$sign = GO_ABAvsTSA_AFP2repress_up$p_value_ABA*GO_ABAvsTSA_AFP2repress_up$p_value_TSA
GO_ABAvsTSA_AFP2repress_up <- arrange(GO_ABAvsTSA_AFP2repress_up, order, sign)
GO_ABAvsTSA_AFP2repress_up <- GO_ABAvsTSA_AFP2repress_up[,c(1,2,3,8,13,15,4,5,6,7,10,11,12,9,14)]
write.xlsx_GO(GO_ABAvsTSA_AFP2repress_up,"GO_ABAvsTSA_AFP2repress_up.xlsx")

GO_GSEA_YTSAYABA_vs_ColTSAColABA <- GSEA_GO_ColTSAColABA[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GO_YTSAYABA[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_Col = padj.x,
    padj_YAFP2 = padj.y,
    ES_Col = ES.x,
    ES_YAFP2 = ES.y,
    NES_Col = NES.x,
    NES_YAFP2 = NES.y,
    edge_Col = leadingEdge.x,
    edge_YAFP2 = leadingEdge.y
  )
GO_GSEA_YTSAYABA_vs_ColTSAColABA$order[GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 1
GO_GSEA_YTSAYABA_vs_ColTSAColABA$order[GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & !GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 2
GO_GSEA_YTSAYABA_vs_ColTSAColABA$order[!GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 3
GO_GSEA_YTSAYABA_vs_ColTSAColABA$sign = GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col*GO_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2
GO_GSEA_YTSAYABA_vs_ColTSAColABA <- arrange(GO_GSEA_YTSAYABA_vs_ColTSAColABA, order, sign)
GO_GSEA_YTSAYABA_vs_ColTSAColABA <- GO_GSEA_YTSAYABA_vs_ColTSAColABA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GO_GSEA_YTSAYABA_vs_ColTSAColABA,"GO_GSEA_YTSAYABA_vs_ColTSAColABA.xlsx")

GO_GSEA_YABAColABA_vs_YTSAColTSA <- GSEA_GO_YABAColABA[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GO_YTSAColTSA[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_ABA = padj.x,
    padj_TSA = padj.y,
    ES_ABA = ES.x,
    ES_TSA = ES.y,
    NES_ABA = NES.x,
    NES_TSA = NES.y,
    edge_ABA = leadingEdge.x,
    edge_TSA = leadingEdge.y
  )
GO_GSEA_YABAColABA_vs_YTSAColTSA$order[GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 1
GO_GSEA_YABAColABA_vs_YTSAColTSA$order[GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & !GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 2
GO_GSEA_YABAColABA_vs_YTSAColTSA$order[!GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 3
GO_GSEA_YABAColABA_vs_YTSAColTSA$sign = GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA*GO_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA
GO_GSEA_YABAColABA_vs_YTSAColTSA <- arrange(GO_GSEA_YABAColABA_vs_YTSAColTSA, order, sign)
GO_GSEA_YABAColABA_vs_YTSAColTSA <- GO_GSEA_YABAColABA_vs_YTSAColTSA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GO_GSEA_YABAColABA_vs_YTSAColTSA,"GO_GSEA_YABAColABA_vs_YTSAColTSA.xlsx")


GO_GSEA_ABAvsTSA <- GSEA_GO_ColABAColGM[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GO_ColTSAColGM[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_ABA = padj.x,
    padj_TSA = padj.y,
    ES_ABA = ES.x,
    ES_TSA = ES.y,
    NES_ABA = NES.x,
    NES_TSA = NES.y,
    edge_ABA = leadingEdge.x,
    edge_TSA = leadingEdge.y
  )
GO_GSEA_ABAvsTSA$order[GO_GSEA_ABAvsTSA$padj_ABA<0.05 & GO_GSEA_ABAvsTSA$padj_TSA<0.05] <- 1
GO_GSEA_ABAvsTSA$order[GO_GSEA_ABAvsTSA$padj_ABA<0.05 & !GO_GSEA_ABAvsTSA$padj_TSA<0.05] <- 2
GO_GSEA_ABAvsTSA$order[!GO_GSEA_ABAvsTSA$padj_ABA<0.05 & GO_GSEA_ABAvsTSA$padj_TSA<0.05] <- 3
GO_GSEA_ABAvsTSA$sign = GO_GSEA_ABAvsTSA$padj_ABA*GO_GSEA_ABAvsTSA$padj_TSA
GO_GSEA_ABAvsTSA <- arrange(GO_GSEA_ABAvsTSA, order, sign)
GO_GSEA_ABAvsTSA <- GO_GSEA_ABAvsTSA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GO_GSEA_ABAvsTSA,"GO_GSEA_ABAvsTSA.xlsx")

##################
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA <- GSEA_GOslimABS_ColTSAColABA[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GOslimABS_YTSAYABA[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_Col = padj.x,
    padj_YAFP2 = padj.y,
    ES_Col = ES.x,
    ES_YAFP2 = ES.y,
    NES_Col = NES.x,
    NES_YAFP2 = NES.y,
    edge_Col = leadingEdge.x,
    edge_YAFP2 = leadingEdge.y
  )
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$order[GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 1
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$order[GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & !GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 2
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$order[!GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 3
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$sign = GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col*GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA <- arrange(GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA, order, sign)
GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA <- GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA,"GOslimABS_GSEA_YTSAYABA_vs_ColTSAColABA.xlsx")

GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA <- GSEA_GOslimABS_YABAColABA[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GOslimABS_YTSAColTSA[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_ABA = padj.x,
    padj_TSA = padj.y,
    ES_ABA = ES.x,
    ES_TSA = ES.y,
    NES_ABA = NES.x,
    NES_TSA = NES.y,
    edge_ABA = leadingEdge.x,
    edge_TSA = leadingEdge.y
  )
GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$order[GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 1
GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$order[GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & !GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 2
GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$order[!GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 3
GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$sign = GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA*GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA
GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA <- arrange(GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA, order, sign)
GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA <- GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA,"GOslimABS_GSEA_YABAColABA_vs_YTSAColTSA.xlsx")

GOslimABS_GSEA_ABAvsTSA <- GSEA_GOslimABS_ColABAColGM[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GOslimABS_ColTSAColGM[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_ABA = padj.x,
    padj_TSA = padj.y,
    ES_ABA = ES.x,
    ES_TSA = ES.y,
    NES_ABA = NES.x,
    NES_TSA = NES.y,
    edge_ABA = leadingEdge.x,
    edge_TSA = leadingEdge.y
  )
GOslimABS_GSEA_ABAvsTSA$order[GOslimABS_GSEA_ABAvsTSA$padj_ABA<0.05 & GOslimABS_GSEA_ABAvsTSA$padj_TSA<0.05] <- 1
GOslimABS_GSEA_ABAvsTSA$order[GOslimABS_GSEA_ABAvsTSA$padj_ABA<0.05 & !GOslimABS_GSEA_ABAvsTSA$padj_TSA<0.05] <- 2
GOslimABS_GSEA_ABAvsTSA$order[!GOslimABS_GSEA_ABAvsTSA$padj_ABA<0.05 & GOslimABS_GSEA_ABAvsTSA$padj_TSA<0.05] <- 3
GOslimABS_GSEA_ABAvsTSA$sign = GOslimABS_GSEA_ABAvsTSA$padj_ABA*GOslimABS_GSEA_ABAvsTSA$padj_TSA
GOslimABS_GSEA_ABAvsTSA <- arrange(GOslimABS_GSEA_ABAvsTSA, order, sign)
GOslimABS_GSEA_ABAvsTSA <- GOslimABS_GSEA_ABAvsTSA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GOslimABS_GSEA_ABAvsTSA,"GOslimABS_GSEA_ABAvsTSA.xlsx")

###

GOslim_GSEA_YTSAYABA_vs_ColTSAColABA <- GSEA_GOslim_ColTSAColABA[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GOslim_YTSAYABA[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_Col = padj.x,
    padj_YAFP2 = padj.y,
    ES_Col = ES.x,
    ES_YAFP2 = ES.y,
    NES_Col = NES.x,
    NES_YAFP2 = NES.y,
    edge_Col = leadingEdge.x,
    edge_YAFP2 = leadingEdge.y
  )
GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$order[GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 1
GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$order[GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & !GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 2
GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$order[!GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col<0.05 & GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2<0.05] <- 3
GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$sign = GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_Col*GOslim_GSEA_YTSAYABA_vs_ColTSAColABA$padj_YAFP2
GOslim_GSEA_YTSAYABA_vs_ColTSAColABA <- arrange(GOslim_GSEA_YTSAYABA_vs_ColTSAColABA, order, sign)
GOslim_GSEA_YTSAYABA_vs_ColTSAColABA <- GOslim_GSEA_YTSAYABA_vs_ColTSAColABA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GOslim_GSEA_YTSAYABA_vs_ColTSAColABA,"GOslim_GSEA_YTSAYABA_vs_ColTSAColABA.xlsx")

GOslim_GSEA_YABAColABA_vs_YTSAColTSA <- GSEA_GOslim_YABAColABA[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GOslim_YTSAColTSA[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_ABA = padj.x,
    padj_TSA = padj.y,
    ES_ABA = ES.x,
    ES_TSA = ES.y,
    NES_ABA = NES.x,
    NES_TSA = NES.y,
    edge_ABA = leadingEdge.x,
    edge_TSA = leadingEdge.y
  )
GOslim_GSEA_YABAColABA_vs_YTSAColTSA$order[GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 1
GOslim_GSEA_YABAColABA_vs_YTSAColTSA$order[GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & !GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 2
GOslim_GSEA_YABAColABA_vs_YTSAColTSA$order[!GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA<0.05 & GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA<0.05] <- 3
GOslim_GSEA_YABAColABA_vs_YTSAColTSA$sign = GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_ABA*GOslim_GSEA_YABAColABA_vs_YTSAColTSA$padj_TSA
GOslim_GSEA_YABAColABA_vs_YTSAColTSA <- arrange(GOslim_GSEA_YABAColABA_vs_YTSAColTSA, order, sign)
GOslim_GSEA_YABAColABA_vs_YTSAColTSA <- GOslim_GSEA_YABAColABA_vs_YTSAColTSA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GOslim_GSEA_YABAColABA_vs_YTSAColTSA,"GOslim_GSEA_YABAColABA_vs_YTSAColTSA.xlsx")

GOslim_GSEA_ABAvsTSA <- GSEA_GOslim_ColABAColGM[,c(1,7,3,5,6,8)] %>%
  left_join(GSEA_GOslim_ColTSAColGM[,c(1,3,5,6,8)], by=c("pathway" = "pathway")) %>%
  rename(
    padj_ABA = padj.x,
    padj_TSA = padj.y,
    ES_ABA = ES.x,
    ES_TSA = ES.y,
    NES_ABA = NES.x,
    NES_TSA = NES.y,
    edge_ABA = leadingEdge.x,
    edge_TSA = leadingEdge.y
  )

GOslim_GSEA_ABAvsTSA$order[GOslim_GSEA_ABAvsTSA$padj_ABA<0.05 & GOslim_GSEA_ABAvsTSA$padj_TSA<0.05] <- 1
GOslim_GSEA_ABAvsTSA$order[GOslim_GSEA_ABAvsTSA$padj_ABA<0.05 & !GOslim_GSEA_ABAvsTSA$padj_TSA<0.05] <- 2
GOslim_GSEA_ABAvsTSA$order[!GOslim_GSEA_ABAvsTSA$padj_ABA<0.05 & GOslim_GSEA_ABAvsTSA$padj_TSA<0.05] <- 3
GOslim_GSEA_ABAvsTSA$sign = GOslim_GSEA_ABAvsTSA$padj_ABA*GOslim_GSEA_ABAvsTSA$padj_TSA
GOslim_GSEA_ABAvsTSA <- arrange(GOslim_GSEA_ABAvsTSA, order, sign)
GOslim_GSEA_ABAvsTSA <- GOslim_GSEA_ABAvsTSA[,c(1,2,3,7,5,9,4,6,8,10)]
write.xlsx_GO_GSEA(GOslim_GSEA_ABAvsTSA,"GOslim_GSEA_ABAvsTSA.xlsx")


outputDirectory<-getwd()
enrichResult_YvsCol_ABA<-WebGestaltR(enrichMethod="ORA",organism="athaliana",
                          enrichDatabase="geneontology_Biological_Process_noRedundant",interestGene=sigresYABAColABA$Entrez, 
                          interestGeneType="entrezgene",referenceGene=master_sheet$Entrez, 
                          referenceGeneType="entrezgene",isOutput=TRUE,
                          outputDirectory=outputDirectory,projectName="YvsCol_ABA")
enrichResult_YvsCol_TSA<-WebGestaltR(enrichMethod="ORA",organism="athaliana",
                                     enrichDatabase="geneontology_Biological_Process_noRedundant",interestGene=sigresYTSAColTSA$Entrez, 
                                     interestGeneType="entrezgene",referenceGene=master_sheet$Entrez, 
                                     referenceGeneType="entrezgene",isOutput=TRUE,
                                     outputDirectory=outputDirectory,projectName="YvsCol_TSA")

enrichResult_YvsCol_both <- enrichResult_YvsCol_ABA %>%
  filter(description %in% enrichResult_YvsCol_TSA$description)
enrichResult_YvsCol_TSAonly <- enrichResult_YvsCol_TSA %>%
  filter(!description %in% enrichResult_YvsCol_ABA$description)

#write functions to save xlsx files (without and with formatting(colored boxes))
write.xlsx_noformatting <- function(data, filename) {
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}
write.xlsx_GO_kmeans <- function(data, filename) {
  wb <- createWorkbook()
  addWorksheet(wb, "sheet")
  writeData(wb, "sheet", data)
  conditionalFormatting(wb, "sheet", cols=5:19, rows = 1:nrow(data),
                        style = c("green", "white"),
                        rule = c(0,0.10),
                        type = "colourScale")
  saveWorkbook(wb, file=filename, overwrite=TRUE)
}

#perform GO analysis with g:profiler gost function and save spreadsheets
gost_k15_4 <- gost(scaledata_k15[which (scaledata_k15$kClusters==1, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_4$result, "gost_k15_4.xlsx")
gost_k15_5 <- gost(scaledata_k15[which (scaledata_k15$kClusters==2, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_5$result, "gost_k15_5.xlsx")
gost_k15_7 <- gost(scaledata_k15[which (scaledata_k15$kClusters==3, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_7$result, "gost_k15_7.xlsx")
gost_k15_3b <- gost(scaledata_k15[which (scaledata_k15$kClusters==4, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_3b$result, "gost_k15_3b.xlsx")
gost_k15_3 <- gost(scaledata_k15[which (scaledata_k15$kClusters==5, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_3$result, "gost_k15_3.xlsx")
gost_k15_8 <- gost(scaledata_k15[which (scaledata_k15$kClusters==6, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_8$result, "gost_k15_8.xlsx")
gost_k15_3c <- gost(scaledata_k15[which (scaledata_k15$kClusters==7, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_3c$result, "gost_k15_3c.xlsx")
gost_k15_9 <- gost(scaledata_k15[which (scaledata_k15$kClusters==8, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_9$result, "gost_k15_9.xlsx")
gost_k15_2 <- gost(scaledata_k15[which (scaledata_k15$kClusters==9, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_2$result, "gost_k15_2.xlsx")
gost_k15_1b <- gost(scaledata_k15[which (scaledata_k15$kClusters==10, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_1b$result, "gost_k15_1b.xlsx")
gost_k15_1 <- gost(scaledata_k15[which (scaledata_k15$kClusters==11, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_1$result, "gost_k15_1.xlsx")
gost_k15_6 <- gost(scaledata_k15[which (scaledata_k15$kClusters==12, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_6$result, "gost_k15_6.xlsx")
gost_k15_6 <- gost(scaledata_k15[which (scaledata_k15$kClusters==13, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_6$result, "gost_k15_6.xlsx")
gost_k15_6 <- gost(scaledata_k15[which (scaledata_k15$kClusters==14, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_6$result, "gost_k15_6.xlsx")
gost_k15_6 <- gost(scaledata_k15[which (scaledata_k15$kClusters==15, arr.ind=TRUE),1], organism = "athaliana", significant = FALSE, numeric_ns = "AFFY_ARAGENE", evcodes = TRUE, correction_method = "gSCS")
write.xlsx_noformatting(gost_k15_6$result, "gost_k15_6.xlsx")


#make sheet combining p values for all clusters
GO_k15 <- gost_k15_4[["result"]][,c(9,10,11,4,3)] %>%
  full_join(gost_k15_5[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_7[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_3b[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_3[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_8[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_3c[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_9[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_2[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_1b[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_1[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_6[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_6[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_6[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size')) %>%
  full_join(gost_k15_6[["result"]][,c(9,10,11,4,3)], by=c('term_id', 'source', 'term_name', 'term_size'))
colnames(GO_k15)[5] <- "4_pval"
colnames(GO_k15)[6] <- "5_pval"
colnames(GO_k15)[7] <- "7_pval"
colnames(GO_k15)[8] <- "3b_pval"
colnames(GO_k15)[9] <- "3_pval"
colnames(GO_k15)[10] <- "8_pval"
colnames(GO_k15)[11] <- "3c_pval"
colnames(GO_k15)[12] <- "9_pval"
colnames(GO_k15)[13] <- "2_pval"
colnames(GO_k15)[14] <- "1b_pval"
colnames(GO_k15)[15] <- "1_pval"
colnames(GO_k15)[16] <- "6_pval"
colnames(GO_k15)[17] <- "6_pval"
colnames(GO_k15)[18] <- "6_pval"
colnames(GO_k15)[19] <- "6_pval"
GO_k15[is.na(GO_k15)] <- 1
GO_k15$sign <- as.numeric(as.character(unlist(GO_k15[5])))*as.numeric(as.character(unlist(GO_k15[6])))*as.numeric(as.character(unlist(GO_k15[7])))*as.numeric(as.character(unlist(GO_k15[8])))*as.numeric(as.character(unlist(GO_k15[9])))*as.numeric(as.character(unlist(GO_k15[10])))*as.numeric(as.character(unlist(GO_k15[11])))*as.numeric(as.character(unlist(GO_k15[12])))*as.numeric(as.character(unlist(GO_k15[13])))*as.numeric(as.character(unlist(GO_k15[14])))*as.numeric(as.character(unlist(GO_k15[15])))*as.numeric(as.character(unlist(GO_k15[16])))*as.numeric(as.character(unlist(GO_k15[17])))*as.numeric(as.character(unlist(GO_k15[18])))*as.numeric(as.character(unlist(GO_k15[19])))
GO_k15 <- arrange(GO_k15, sign)
write.xlsx_GO_kmeans(GO_k15,"GO_k15.xlsx")
