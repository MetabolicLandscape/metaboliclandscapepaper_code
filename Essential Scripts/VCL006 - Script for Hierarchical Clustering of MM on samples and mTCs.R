library(dplyr)

MMGEO<-read.delim("Consensus_Mix_Matrix_20173006_GEO_with_normal_.txt", header=TRUE, sep="\t")
MMGEO[1:5,1:5]
GEO_mTCs<-as.matrix(read.delim("TC_Selection_Top3_GEO_REDUCED.txt", header=FALSE, sep="\t"))
MMTCGA<-read.delim("Consensus_Mix_Matrix_20170906_Duplicate_removed_TCGA_data_.txt", header=TRUE, sep="\t")
MMTCGA[1:5,1:5]
TCGA_mTCs<-as.matrix(read.delim("TC_Selection_Top3_TCGA_REDUCED.txt", header=FALSE, sep="\t"))
MMGDSC<-read.delim("Consensus_Mix_Matrix_20170209_GDSC_data_.txt", header=TRUE, sep="\t")
MMGDSC[1:5,1:5]
GDSC_mTCs<-as.matrix(read.delim("TC_Selection_Top3_GDSC_REDUCED.txt", header=FALSE, sep="\t"))
MMCCLE<-read.delim("Consensus_Mix_Matrix_20170209_CCLE_data_.txt", header=TRUE, sep="\t")
MMCCLE[1:5,1:5]
CCLE_mTCs<-as.matrix(read.delim("TC_Selection_Top3_CCLE_REDUCED.txt", header=FALSE, sep="\t"))

MMGEO<-as.data.frame(t(MMGEO))
colnames(MMGEO)<-paste(colnames(MMGEO),"_GEO",sep="")
MMGEO<-select(MMGEO,GEO_mTCs)
write.table(MMGEO,"GEO_MM_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")
MMTCGA<-as.data.frame(t(MMTCGA))
colnames(MMTCGA)<-paste(colnames(MMTCGA),"_TCGA",sep="")
MMTCGA<-select(MMTCGA,TCGA_mTCs)
write.table(MMTCGA,"TCGA_MM_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")
MMGDSC<-as.data.frame(t(MMGDSC))
colnames(MMGDSC)<-paste(colnames(MMGDSC),"_GDSC",sep="")
MMGDSC<-select(MMGDSC,GDSC_mTCs)
write.table(MMGDSC,"GDSC_MM_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")
MMCCLE<-as.data.frame(t(MMCCLE))
colnames(MMCCLE)<-paste(colnames(MMCCLE),"_CCLE",sep="")
MMCCLE<-select(MMCCLE,CCLE_mTCs)
write.table(MMCCLE,"CCLE_MM_REVISED_Metabolic_Sel.txt",col.names=TRUE,row.names=TRUE,sep="\t")


MM<-MMGEO
sample_dist <- as.dist(1-cor(t(MM)))
sample_hc <- hclust(sample_dist, method="ward.D2")
sample_dend <- as.dendrogram(sample_hc)
sample_order <- as.data.frame(sample_hc$order)
sample_order[,2] <- sample_hc$labels[sample_order[,1]]
glimpse(sample_order)
TC_dist <- as.dist(1-cor(MM))
TC_hc <- hclust(TC_dist, method="ward.D2")
TC_dend <- as.dendrogram(TC_hc)
TC_order <- as.data.frame(TC_hc$order)
TC_order[,2] <- TC_hc$labels[TC_order[,1]]
glimpse(TC_order)
nrow(MM)
Reorder_One <- MM[match(sample_order[,2], row.names(MM)),]
Reorder_One[1:5,1:5]
nrow(Reorder_One)
Reorder_Two <- t(Reorder_One)
Reorder_Two <- Reorder_Two[match(TC_order[,2], row.names(Reorder_Two)),]
Reordered <- t(Reorder_Two)
Reordered[1:5,1:5]
nrow(Reordered)
write.table(Reordered,"GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")

MM<-MMTCGA
sample_dist <- as.dist(1-cor(t(MM)))
sample_hc <- hclust(sample_dist, method="ward.D2")
sample_dend <- as.dendrogram(sample_hc)
sample_order <- as.data.frame(sample_hc$order)
sample_order[,2] <- sample_hc$labels[sample_order[,1]]
TC_dist <- as.dist(1-cor(MM))
TC_hc <- hclust(TC_dist, method="ward.D2")
TC_dend <- as.dendrogram(TC_hc)
TC_order <- as.data.frame(TC_hc$order)
TC_order[,2] <- TC_hc$labels[TC_order[,1]]
Reorder_One <- MM[match(sample_order[,2], row.names(MM)),]
Reorder_Two <- t(Reorder_One)
Reorder_Two <- Reorder_Two[match(TC_order[,2], row.names(Reorder_Two)),]
Reordered <- t(Reorder_Two)
write.table(Reordered,"TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")


MM<-MMGDSC
sample_dist <- as.dist(1-cor(t(MM)))
sample_hc <- hclust(sample_dist, method="ward.D2")
sample_dend <- as.dendrogram(sample_hc)
sample_order <- as.data.frame(sample_hc$order)
sample_order[,2] <- sample_hc$labels[sample_order[,1]]
TC_dist <- as.dist(1-cor(MM))
TC_hc <- hclust(TC_dist, method="ward.D2")
TC_dend <- as.dendrogram(TC_hc)
TC_order <- as.data.frame(TC_hc$order)
TC_order[,2] <- TC_hc$labels[TC_order[,1]]
Reorder_One <- MM[match(sample_order[,2], row.names(MM)),]
Reorder_Two <- t(Reorder_One)
Reorder_Two <- Reorder_Two[match(TC_order[,2], row.names(Reorder_Two)),]
Reordered <- t(Reorder_Two)
write.table(Reordered,"GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")


MM<-MMCCLE
sample_dist <- as.dist(1-cor(t(MM)))
sample_hc <- hclust(sample_dist, method="ward.D2")
sample_dend <- as.dendrogram(sample_hc)
sample_order <- as.data.frame(sample_hc$order)
sample_order[,2] <- sample_hc$labels[sample_order[,1]]
TC_dist <- as.dist(1-cor(MM))
TC_hc <- hclust(TC_dist, method="ward.D2")
TC_dend <- as.dendrogram(TC_hc)
TC_order <- as.data.frame(TC_hc$order)
TC_order[,2] <- TC_hc$labels[TC_order[,1]]
Reorder_One <- MM[match(sample_order[,2], row.names(MM)),]
Reorder_Two <- t(Reorder_One)
Reorder_Two <- Reorder_Two[match(TC_order[,2], row.names(Reorder_Two)),]
Reordered <- t(Reorder_Two)
write.table(Reordered,"CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", col.names=TRUE, row.names=TRUE, sep="\t")

