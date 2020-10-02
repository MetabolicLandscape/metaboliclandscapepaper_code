library(dplyr)
library(gplots)
library(RColorBrewer)

##colors for heatmap
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-10,-3,length=100),  seq(-2.99,3,length=100),  seq(3.01,10,length=100))

##absdistance, wardD2 clustering for each gene expression dataset seperately
GEO_GSEA<-read.delim("GEO_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(GEO_GSEA)<-gsub("_minus_log_pvalue","",colnames(GEO_GSEA))
GEO_selection<-as.matrix(read.delim("TC_Selection_Top3.txt",header=TRUE,sep="\t"))
GEO_GSEA<-GEO_GSEA[ , !(names(GEO_GSEA) %in% GEO_selection)]

order<-row.names(read.delim("GEO_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(GEO_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

GEO_Reorder_GSEA <- GEO_GSEA[match(order, row.names(GEO_GSEA)),]
GEO_Reorder_GSEA_TCs <- t(GEO_Reorder_GSEA)
GEO_Reorder_GSEA_TCs <- GEO_Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(GEO_Reorder_GSEA_TCs)),]
GEO_GSEA_Ordered <- t(GEO_Reorder_GSEA_TCs)
GEO_GSEA_Ordered[1:5,1:5]

pdf("GEO_GSEA_wardD2_absdist_heatmap_NOTMETABOLIC.pdf", width=((30/132)*723), height=70, pointsize=14)
heatmap.2(as.matrix(GEO_GSEA_Ordered), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none",  Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()


##for TCGA
TCGA_GSEA<-read.delim("TCGA_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(TCGA_GSEA)<-gsub("_minus_log_pvalue","",colnames(TCGA_GSEA))
TCGA_selection<-as.matrix(read.delim("TC_Selection_Top3_TCGA.txt",header=TRUE,sep="\t"))
TCGA_GSEA<-TCGA_GSEA[ , !(names(TCGA_GSEA) %in% TCGA_selection)]

order<-row.names(read.delim("TCGA_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(TCGA_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

Reorder_GSEA <- TCGA_GSEA[match(order, row.names(TCGA_GSEA)),]
Reorder_GSEA_TCs <- t(Reorder_GSEA)
Reorder_GSEA_TCs <- Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(Reorder_GSEA_TCs)),]
TCGA_GSEA_Ordered <- t(Reorder_GSEA_TCs)
TCGA_GSEA_Ordered[1:5,1:5]

pdf("TCGA_GSEA_wardD2_absdist_heatmap_NOTMETABOLIC.pdf", width=((30/151)*723), height=70, pointsize=14)
heatmap.2(as.matrix(TCGA_GSEA_Ordered), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none",  Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()


##for CCLE
CCLE_GSEA<-read.delim("CCLE_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(CCLE_GSEA)<-gsub("_minus_log_pvalue","",colnames(CCLE_GSEA))
CCLE_selection<-as.matrix(read.delim("TC_Selection_Top3_CCLE.txt",header=TRUE,sep="\t"))
CCLE_GSEA<-CCLE_GSEA[ , !(names(CCLE_GSEA) %in% CCLE_selection)]

order<-row.names(read.delim("CCLE_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(CCLE_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

Reorder_GSEA <- CCLE_GSEA[match(order, row.names(CCLE_GSEA)),]
Reorder_GSEA_TCs <- t(Reorder_GSEA)
Reorder_GSEA_TCs <- Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(Reorder_GSEA_TCs)),]
CCLE_GSEA_Ordered <- t(Reorder_GSEA_TCs)

CCLE_GSEA_Ordered[1:5,1:5]

pdf("CCLE_GSEA_wardD2_absdist_heatmap_NOTMETABOLIC.pdf", width=((30/136)*330), height=70, pointsize=14)
heatmap.2(as.matrix(CCLE_GSEA_Ordered), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()


##for GDSC
GDSC_GSEA<-read.delim("GDSC_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(GDSC_GSEA)<-gsub("_minus_log_pvalue","",colnames(GDSC_GSEA))
GDSC_selection<-as.matrix(read.delim("TC_Selection_Top3_GDSC.txt",header=TRUE,sep="\t"))
GDSC_GSEA<-GDSC_GSEA[ , !(names(GDSC_GSEA) %in% GDSC_selection)]

order<-row.names(read.delim("GDSC_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(GDSC_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

Reorder_GSEA <- GDSC_GSEA[match(order, row.names(GDSC_GSEA)),]
Reorder_GSEA_TCs <- t(Reorder_GSEA)
Reorder_GSEA_TCs <- Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(Reorder_GSEA_TCs)),]
GDSC_GSEA_Ordered <- t(Reorder_GSEA_TCs)
GDSC_GSEA_Ordered[1:5,1:5]

pdf("GDSC_GSEA_wardD2_absdist_heatmap_NOTMETABOLIC.pdf", width=((30/137)*330), height=70, pointsize=14)
heatmap.2(as.matrix(GDSC_GSEA_Ordered), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()






##Reorder all on GEO-geneset-order to be able to align y-axes in heatmaps

GEO_GSEA<-read.delim("GEO_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t")
order<-row.names(GEO_GSEA)
glimpse(order)


TCGA_GSEA<-read.delim("TCGA_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t")

Reorder_GSEA <- TCGA_GSEA[match(order, row.names(TCGA_GSEA)),]
Reorder_GSEA[1:5,1:5]

pdf("TCGA_GSEA_wardD2_absdist_heatmap_GEOorder.pdf", width=30, height=70, pointsize=14)
heatmap.2(as.matrix(Reorder_GSEA), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none",  Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()

CCLE_GSEA<-read.delim("CCLE_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t")
Reorder_GSEA <- CCLE_GSEA[match(order, row.names(CCLE_GSEA)),]
Reorder_GSEA[1:5,1:5]
pdf("CCLE_GSEA_wardD2_absdist_heatmap_GEOorder.pdf", width=30, height=70, pointsize=14)
heatmap.2(as.matrix(Reorder_GSEA), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none",  Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()

GDSC_GSEA<-read.delim("GDSC_GSEA_wardD2_absdist.txt",header=TRUE,sep="\t")
Reorder_GSEA <- GDSC_GSEA[match(order, row.names(TCGA_GSEA)),]
Reorder_GSEA[1:5,1:5]
pdf("GDSC_GSEA_wardD2_absdist_heatmap_GEOorder.pdf", width=30, height=70, pointsize=14)
heatmap.2(as.matrix(Reorder_GSEA), breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none",  Colv=FALSE, Rowv=FALSE, key=FALSE)
dev.off()





##Idem, but for NON METABOLIC TCS
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-10,-3,length=100),  seq(-2.99,3,length=100),  seq(3.01,10,length=100))

GEO_GSEA<-read.delim("GEO_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(GEO_GSEA)<-gsub("_minus_log_pvalue","",colnames(GEO_GSEA))
GEO_selection<-as.matrix(read.delim("TC_Selection_Top3.txt",header=TRUE,sep="\t"))
GEO_GSEA<-GEO_GSEA[ , !(names(GEO_GSEA) %in% GEO_selection)]

order<-row.names(read.delim("GEO_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(GEO_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

GEO_Reorder_GSEA <- GEO_GSEA[match(order, row.names(GEO_GSEA)),]
GEO_Reorder_GSEA_TCs <- t(GEO_Reorder_GSEA)
GEO_Reorder_GSEA_TCs <- GEO_Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(GEO_Reorder_GSEA_TCs)),]
GEO_GSEA_Ordered <- t(GEO_Reorder_GSEA_TCs)
GEO_GSEA_Ordered[1:5,1:5]

pdf("GEO_GSEA_wardD2_absdist_heatmap_ALLDATASETS_GENESET_CLUSTERED_NOTMETABOLIC.pdf", width=15, height=30, pointsize=14)
heatmap.2(as.matrix(GEO_GSEA_Ordered),main = "GEO GSEA data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,15),   col=my_palette,  dendrogram="none", Colv=FALSE, cexRow = 0.4, cexCol = 0.4, Rowv=FALSE, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()

TCGA_GSEA<-read.delim("TCGA_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(TCGA_GSEA)<-gsub("_minus_log_pvalue","",colnames(TCGA_GSEA))
TCGA_selection<-as.matrix(read.delim("TC_Selection_Top3_TCGA.txt",header=TRUE,sep="\t"))
TCGA_GSEA<-TCGA_GSEA[ , !(names(TCGA_GSEA) %in% TCGA_selection)]

order<-row.names(read.delim("TCGA_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(TCGA_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

Reorder_GSEA <- TCGA_GSEA[match(order, row.names(TCGA_GSEA)),]
Reorder_GSEA_TCs <- t(Reorder_GSEA)
Reorder_GSEA_TCs <- Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(Reorder_GSEA_TCs)),]
TCGA_GSEA_Ordered <- t(Reorder_GSEA_TCs)
TCGA_GSEA_Ordered[1:5,1:5]

pdf("TCGA_GSEA_wardD2_absdist_heatmap_ALLDATASETS_GENESET_CLUSTERED_NOTMETABOLIC.pdf", width=15, height=30, pointsize=14)
heatmap.2(as.matrix(TCGA_GSEA_Ordered),main = "TCGA GSEA data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,15),   col=my_palette,  dendrogram="none", Colv=FALSE, cexRow = 0.4, cexCol = 0.4, Rowv=FALSE, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()



CCLE_GSEA<-read.delim("CCLE_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(CCLE_GSEA)<-gsub("_minus_log_pvalue","",colnames(CCLE_GSEA))
CCLE_selection<-as.matrix(read.delim("TC_Selection_Top3_CCLE.txt",header=TRUE,sep="\t"))
CCLE_GSEA<-CCLE_GSEA[ , !(names(CCLE_GSEA) %in% CCLE_selection)]

order<-row.names(read.delim("CCLE_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(CCLE_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

Reorder_GSEA <- CCLE_GSEA[match(order, row.names(CCLE_GSEA)),]
Reorder_GSEA_TCs <- t(Reorder_GSEA)
Reorder_GSEA_TCs <- Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(Reorder_GSEA_TCs)),]
CCLE_GSEA_Ordered <- t(Reorder_GSEA_TCs)

CCLE_GSEA_Ordered[1:5,1:5]

pdf("CCLE_GSEA_wardD2_absdist_heatmap_ALLDATASETS_GENESET_CLUSTERED_NOTMETABOLIC.pdf", width=15, height=30, pointsize=14)
heatmap.2(as.matrix(CCLE_GSEA_Ordered),main = "CCLE GSEA data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,15),   col=my_palette,  dendrogram="none", Colv=FALSE, cexRow = 0.4, cexCol = 0.4, Rowv=FALSE, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()



GDSC_GSEA<-read.delim("GDSC_Combined_databases_Metabolic_Selection_GSEA_Data.txt",header=TRUE,sep="\t")
colnames(GDSC_GSEA)<-gsub("_minus_log_pvalue","",colnames(GDSC_GSEA))
GDSC_selection<-as.matrix(read.delim("TC_Selection_Top3_GDSC.txt",header=TRUE,sep="\t"))
GDSC_GSEA<-GDSC_GSEA[ , !(names(GDSC_GSEA) %in% GDSC_selection)]

order<-row.names(read.delim("GDSC_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t"))
glimpse(order)

ith_TC_dist <- as.dist(1-cor(abs(GDSC_GSEA)))
ith_TC_hc <- hclust(ith_TC_dist, method="ward.D2")
ith_TC_dend <- as.dendrogram(ith_TC_hc)
ith_TC_order <- as.data.frame(ith_TC_hc$order)
ith_TC_order[,2] <- ith_TC_hc$labels[ith_TC_order[,1]]
ith_TC_order[1:5,]

Reorder_GSEA <- GDSC_GSEA[match(order, row.names(GDSC_GSEA)),]
Reorder_GSEA_TCs <- t(Reorder_GSEA)
Reorder_GSEA_TCs <- Reorder_GSEA_TCs[match(ith_TC_order[,2], row.names(Reorder_GSEA_TCs)),]
GDSC_GSEA_Ordered <- t(Reorder_GSEA_TCs)
GDSC_GSEA_Ordered[1:5,1:5]

pdf("GDSC_GSEA_wardD2_absdist_heatmap_ALLDATASETS_GENESET_CLUSTERED_NOTMETABOLIC.pdf", width=15, height=30, pointsize=14)
heatmap.2(as.matrix(GDSC_GSEA_Ordered),main = "GDSC GSEA data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,15),   col=my_palette,  dendrogram="none", Colv=FALSE, cexRow = 0.4, cexCol = 0.4, Rowv=FALSE, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()



