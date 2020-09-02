library(dplyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(BBmisc)


genelevel<-read.delim("Genelevel_using_jetscore_Consensus_Independent_Components_20173006_GEO_with_normal_.txt", header=TRUE, sep=" ")
genelevel[1:5,1:5]
colnames(genelevel)<-paste(colnames(genelevel),"_GEO",sep="")
mTC_selection<-as.matrix(read.delim("TC_Selection_Top3_GEO_REDUCED.txt",header=FALSE,sep="\t"))

data<-select(genelevel,mTC_selection)
write.table(data,"GEO_Genelevel_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")


#clustering all
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("GEO_Genelevel_REVISED_mTC_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("GEO_Genelevel_REVISED_mTC_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "GEO_Genelevel_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
Genelevel_Data_Subset_Ordered<-as.data.frame(Genelevel_Data_Subset_Ordered)

#clustering just above 10

N<-10
out1 <- paste("GEO_Genelevel_REVISED_mTCs_Geneselection_SD",N,".txt",sep="")
out2 <- paste("GEO_Genelevel_REVISED_mTCs_subset_SD",N,"_genes.txt",sep="")
Significant_Genes <- vector()
for (i in 1:length(mTC_selection)){
  oneTC <- select(Genelevel_Data_Subset_Ordered,mTC_selection[i])
  oneTC[,2]<-rownames(oneTC)
  colnames(oneTC)<-c("GeneWeight","EntrezID")
  oneTC_SORTED<-arrange(oneTC,desc(GeneWeight))
  oneTC_SORTED$EntrezID<-as.integer(oneTC_SORTED$EntrezID)
  oneTC_Significant <- filter(oneTC_SORTED, abs(GeneWeight) >= N)
  oneTC_Significant_EntrezIDOnly <- select(oneTC_Significant, EntrezID)
  Significant_Genes <- rbind(Significant_Genes,oneTC_Significant_EntrezIDOnly)
}
Significant_Genes<-distinct(Significant_Genes)
write.table(Significant_Genes,out1,sep="\t")
Genelevel_Data_Significant_Subset <- subset(Genelevel_Data_Subset_Ordered,rownames(Genelevel_Data_Subset_Ordered) %in% Significant_Genes$EntrezID)
write.table(Genelevel_Data_Significant_Subset,out2,row.names=TRUE,col.names=TRUE,sep="\t")

data<-Genelevel_Data_Significant_Subset
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("GEO_Genelevel_REVISED_mTC_SD10_GENES_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("GEO_Genelevel_REVISED_mTC_SD10_GENES_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "GEO_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")

outfile<-paste("GEO_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-3,-1,length=100),  seq(-0.99,1,length=100),  seq(1.01,3,length=100))
pdf(outfile, width=15, height=45, pointsize=14)
heatmap.2(as.matrix(Genelevel_Data_Subset_Ordered),main = "GEO revisedmTC sample level data clustered SD10 genes",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,5),   col=my_palette,  dendrogram="none", Colv=FALSE, Rowv=FALSE, cexRow = 0.4, cexCol = 0.4, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()


genelevel<-read.delim("Consensus_Independent_Components_20170906_Duplicate_removed_TCGA_data_.txt", header=TRUE, sep="\t")
genelevel[1:5,1:5]
row.names(genelevel)<-genelevel$X
genelevel$X<-NULL
colnames(genelevel)<-paste(colnames(genelevel),"_TCGA",sep="")
mTC_selection<-as.matrix(read.delim("TC_Selection_Top3_TCGA_REDUCED.txt",header=FALSE,sep="\t"))
data<-select(genelevel,mTC_selection)
data[1:5,1:5]
write.table(data,"TCGA_Genelevel_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("TCGA_Genelevel_REVISED_mTC_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("TCGA_Genelevel_REVISED_mTC_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "TCGA_Genelevel_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
Genelevel_Data_Subset_Ordered<-as.data.frame(Genelevel_Data_Subset_Ordered)
N<-10
out1 <- paste("TCGA_Genelevel_REVISED_mTCs_Geneselection_SD",N,".txt",sep="")
out2 <- paste("TCGA_Genelevel_REVISED_mTCs_subset_SD",N,"_genes.txt",sep="")
Significant_Genes <- vector()
for (i in 1:length(mTC_selection)){
  oneTC <- select(Genelevel_Data_Subset_Ordered,mTC_selection[i])
  oneTC[,2]<-rownames(oneTC)
  colnames(oneTC)<-c("GeneWeight","EntrezID")
  oneTC_SORTED<-arrange(oneTC,desc(GeneWeight))
  oneTC_SORTED$EntrezID<-as.integer(oneTC_SORTED$EntrezID)
  oneTC_Significant <- filter(oneTC_SORTED, abs(GeneWeight) >= N)
  oneTC_Significant_EntrezIDOnly <- select(oneTC_Significant, EntrezID)
  Significant_Genes <- rbind(Significant_Genes,oneTC_Significant_EntrezIDOnly)
}
Significant_Genes<-distinct(Significant_Genes)
write.table(Significant_Genes,out1,sep="\t")
Genelevel_Data_Significant_Subset <- subset(Genelevel_Data_Subset_Ordered,rownames(Genelevel_Data_Subset_Ordered) %in% Significant_Genes$EntrezID)
write.table(Genelevel_Data_Significant_Subset,out2,row.names=TRUE,col.names=TRUE,sep="\t")
data<-Genelevel_Data_Significant_Subset
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("TCGA_Genelevel_REVISED_mTC_SD10_GENES_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("TCGA_Genelevel_REVISED_mTC_SD10_GENES_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "TCGA_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
outfile<-paste("TCGA_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-3,-1,length=100),  seq(-0.99,1,length=100),  seq(1.01,3,length=100))
pdf(outfile, width=15, height=45, pointsize=14)
heatmap.2(as.matrix(Genelevel_Data_Subset_Ordered),main = "TCGA revisedmTC sample level data clustered SD10 genes",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,5),   col=my_palette,  dendrogram="none", Colv=FALSE, Rowv=FALSE, cexRow = 0.4, cexCol = 0.4, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()




genelevel<-read.delim("Genelevel_using_jetscore_Consensus_Independent_Components_20170209_CCLE_data_.txt", header=TRUE, sep=" ")
genelevel[1:5,1:5]
colnames(genelevel)<-paste(colnames(genelevel),"_CCLE",sep="")
mTC_selection<-as.matrix(read.delim("TC_Selection_Top3_CCLE_REDUCED.txt",header=FALSE,sep="\t"))
data<-select(genelevel,mTC_selection)
write.table(data,"CCLE_Genelevel_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("CCLE_Genelevel_REVISED_mTC_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("CCLE_Genelevel_REVISED_mTC_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "CCLE_Genelevel_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
Genelevel_Data_Subset_Ordered<-as.data.frame(Genelevel_Data_Subset_Ordered)
N<-10
out1 <- paste("CCLE_Genelevel_REVISED_mTCs_Geneselection_SD",N,".txt",sep="")
out2 <- paste("CCLE_Genelevel_REVISED_mTCs_subset_SD",N,"_genes.txt",sep="")
Significant_Genes <- vector()
for (i in 1:length(mTC_selection)){
  oneTC <- select(Genelevel_Data_Subset_Ordered,mTC_selection[i])
  oneTC[,2]<-rownames(oneTC)
  colnames(oneTC)<-c("GeneWeight","EntrezID")
  oneTC_SORTED<-arrange(oneTC,desc(GeneWeight))
  oneTC_SORTED$EntrezID<-as.integer(oneTC_SORTED$EntrezID)
  oneTC_Significant <- filter(oneTC_SORTED, abs(GeneWeight) >= N)
  oneTC_Significant_EntrezIDOnly <- select(oneTC_Significant, EntrezID)
  Significant_Genes <- rbind(Significant_Genes,oneTC_Significant_EntrezIDOnly)
}
Significant_Genes<-distinct(Significant_Genes)
write.table(Significant_Genes,out1,sep="\t")
Genelevel_Data_Significant_Subset <- subset(Genelevel_Data_Subset_Ordered,rownames(Genelevel_Data_Subset_Ordered) %in% Significant_Genes$EntrezID)
write.table(Genelevel_Data_Significant_Subset,out2,row.names=TRUE,col.names=TRUE,sep="\t")
data<-Genelevel_Data_Significant_Subset
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("CCLE_Genelevel_REVISED_mTC_SD10_GENES_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("CCLE_Genelevel_REVISED_mTC_SD10_GENES_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "CCLE_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
outfile<-paste("CCLE_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-3,-1,length=100),  seq(-0.99,1,length=100),  seq(1.01,3,length=100))
pdf(outfile, width=15, height=45, pointsize=14)
heatmap.2(as.matrix(Genelevel_Data_Subset_Ordered),main = "CCLE revisedmTC sample level data clustered SD10 genes",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,5),   col=my_palette,  dendrogram="none", Colv=FALSE, Rowv=FALSE, cexRow = 0.4, cexCol = 0.4, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()









genelevel<-read.delim("Genelevel_Consensus_Independent_Components_20170209_GDSC_data_.txt", header=TRUE, sep="\t")
genelevel[1:5,1:5]
colnames(genelevel)<-paste(colnames(genelevel),"_GDSC",sep="")
mTC_selection<-as.matrix(read.delim("TC_Selection_Top3_GDSC_REDUCED.txt",header=FALSE,sep="\t"))
data<-select(genelevel,mTC_selection)
write.table(data,"GDSC_Genelevel_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("GDSC_Genelevel_REVISED_mTC_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("GDSC_Genelevel_REVISED_mTC_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "GDSC_Genelevel_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
Genelevel_Data_Subset_Ordered<-as.data.frame(Genelevel_Data_Subset_Ordered)
N<-10
out1 <- paste("GDSC_Genelevel_REVISED_mTCs_Geneselection_SD",N,".txt",sep="")
out2 <- paste("GDSC_Genelevel_REVISED_mTCs_subset_SD",N,"_genes.txt",sep="")
Significant_Genes <- vector()
for (i in 1:length(mTC_selection)){
  oneTC <- select(Genelevel_Data_Subset_Ordered,mTC_selection[i])
  oneTC[,2]<-rownames(oneTC)
  colnames(oneTC)<-c("GeneWeight","EntrezID")
  oneTC_SORTED<-arrange(oneTC,desc(GeneWeight))
  oneTC_SORTED$EntrezID<-as.integer(oneTC_SORTED$EntrezID)
  oneTC_Significant <- filter(oneTC_SORTED, abs(GeneWeight) >= N)
  oneTC_Significant_EntrezIDOnly <- select(oneTC_Significant, EntrezID)
  Significant_Genes <- rbind(Significant_Genes,oneTC_Significant_EntrezIDOnly)
}
Significant_Genes<-distinct(Significant_Genes)
write.table(Significant_Genes,out1,sep="\t")
Genelevel_Data_Significant_Subset <- subset(Genelevel_Data_Subset_Ordered,rownames(Genelevel_Data_Subset_Ordered) %in% Significant_Genes$EntrezID)
write.table(Genelevel_Data_Significant_Subset,out2,row.names=TRUE,col.names=TRUE,sep="\t")
data<-Genelevel_Data_Significant_Subset
genedist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(genedist, method="ward.D2")
dend <- as.dendrogram(hc)
order <- as.data.frame(hc$order)
order$LABEL <- hc$labels[order[,1]]
pdf("GDSC_Genelevel_REVISED_mTC_SD10_GENES_Gene_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(dend)
dev.off()
TCdist <- as.dist(1-cor(as.matrix(data)))
TChc <- hclust(TCdist, method="ward.D2")
TCdend <- as.dendrogram(TChc)
TCorder <- as.data.frame(TChc$order)
TCorder$LABEL <- TChc$labels[TCorder[,1]]
pdf("GDSC_Genelevel_REVISED_mTC_SD10_GENES_TC_Dend_ward.pdf", width=20,height=15,pointsize=15)
plot(TCdend)
dev.off()
Genelevel_Data_Subset_Reorder_Genes <- data[match(order[,2], row.names(data)),]
Genelevel_Data_Subset_Reorder_Genes_TCs <- t(Genelevel_Data_Subset_Reorder_Genes)
Genelevel_Data_Subset_Reorder_Genes_TCs <- Genelevel_Data_Subset_Reorder_Genes_TCs[match(TCorder[,2], row.names(Genelevel_Data_Subset_Reorder_Genes_TCs)),]
Genelevel_Data_Subset_Ordered <- t(Genelevel_Data_Subset_Reorder_Genes_TCs)
write.table(Genelevel_Data_Subset_Ordered, "GDSC_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")
outfile<-paste("GDSC_Genelevel_REVISED_Metabolic_Sel_SD10_GENES_1-cor_ward_clustered.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-3,-1,length=100),  seq(-0.99,1,length=100),  seq(1.01,3,length=100))
pdf(outfile, width=15, height=45, pointsize=14)
heatmap.2(as.matrix(Genelevel_Data_Subset_Ordered),main = "GDSC revisedmTC sample level data clustered SD10 genes",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(3,5),   col=my_palette,  dendrogram="none", Colv=FALSE, Rowv=FALSE, cexRow = 0.4, cexCol = 0.4, keysize=0.4,lhei=c(1,10),lwid=c(1,10))
dev.off()

