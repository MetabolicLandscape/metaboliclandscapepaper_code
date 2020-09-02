library(dplyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(BBmisc)

TYPE_ORDER<-c("Breast cancer","CNS cancer","Endocrine cancer","Gastrointestinal cancer","HNSCC","Leukemia","Lung cancer","Lymphoma","Melanoma","Myeloma","Sarcoma","Thyroid cancer","Urogenital cancer","Mesothelioma","Thymoma","Normal")
NUMBER_TYPES<-length(TYPE_ORDER)
colors<-c("Breast cancer"="yellowgreen","CNS cancer"="yellow3","Endocrine cancer"="wheat4","Gastrointestinal cancer"="violetred3","HNSCC"="turquoise4","Leukemia"="tomato4","Lung cancer"="thistle4","Lymphoma"="tan4","Melanoma"="steelblue4","Myeloma"="springgreen4","Sarcoma"="snow4","Thyroid cancer"="slateblue4","Urogenital cancer"="skyblue4","Mesothelioma"="sienna4","Thymoma"="seagreen4","Normal"="salmon4")
annotation<-read.delim("GEO_Novel_Sample_annotation.txt",header=TRUE,sep="\t")
annotation$CLASS<-NULL
annotation$TYPE<-annotation$TYPE_NOVEL
annotation$TYPE_NOVEL<-NULL

data<-read.delim("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")
cutoffheight<-6.0

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
sample_groups<-cutree(hc,k=NULL,h=cutoffheight)
sample_groups<-as.data.frame(sample_groups)
sample_groups$LABEL<-row.names(sample_groups)
sample_groups[1:5,]
sample_groups$V2<-sample_groups$LABEL
alt_groups<-left_join(order,sample_groups,by="V2")
alt_groups<-select(alt_groups,sample_groups,LABEL)
colnames(alt_groups)<-c("GROUP","LABEL")
alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
alt_groups[1:5,]
outfile<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_LABELS.txt",sep="")
write.table(alt_groups,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
k1<-max(alt_groups$GROUP)
k1
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col1a<-sample(color,k1,replace=TRUE)
col1b<-c(1:k1)
col1<-as.data.frame(cbind(col1a,col1b))
colnames(col1)<-c("col1a","GROUP")
col1$GROUP<-as.numeric(col1$GROUP)
merge<-left_join(alt_groups,col1,by="GROUP")
merge[1:5,]
outfile<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_COLORS.txt",sep="")
write.table(merge,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
RowColors<-as.vector(merge$col1a)
outfile<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"NOT_NORMALIZED.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-0.3,-0.05,length=100),  seq(-0.049,0.05,length=100),  seq(0.051,0.3,length=100))
pdf(outfile, width=5, height=15, pointsize=14)
heatmap.2(as.matrix(data),main = "GEO sample level data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=FALSE, Rowv=FALSE, keysize=1)
dev.off()
Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(data)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes
write.table(Data_Ordered_Classes, "GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.txt", sep="\t")
ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=TYPE_ORDER) + scale_shape_manual(values=c(rep(95, NUMBER_TYPES))) + scale_size_manual(values=c(rep((20), NUMBER_TYPES))) + scale_color_manual(values=colors)
ggsave("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.pdf", width = NUMBER_TYPES, height= 170, dpi=75, limitsize=FALSE)










data<-read.delim("TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")
cutoffheight<-5.2
annotation<-read.delim("TCGA_Novel_Sample_annotation.txt",header=TRUE,sep="\t")
annotation$CLASS<-NULL
annotation$TYPE<-annotation$TYPE_NOVEL
annotation$TYPE_NOVEL<-NULL

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
sample_groups<-cutree(hc,k=NULL,h=cutoffheight)
sample_groups<-as.data.frame(sample_groups)
sample_groups$LABEL<-row.names(sample_groups)
sample_groups[1:5,]
sample_groups$V2<-sample_groups$LABEL
alt_groups<-left_join(order,sample_groups,by="V2")
alt_groups<-select(alt_groups,sample_groups,LABEL)
colnames(alt_groups)<-c("GROUP","LABEL")
alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
alt_groups[1:5,]
outfile<-paste("TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_LABELS.txt",sep="")
write.table(alt_groups,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
k1<-max(alt_groups$GROUP)
k1
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col1a<-sample(color,k1,replace=TRUE)
col1b<-c(1:k1)
col1<-as.data.frame(cbind(col1a,col1b))
colnames(col1)<-c("col1a","GROUP")
col1$GROUP<-as.numeric(col1$GROUP)
merge<-left_join(alt_groups,col1,by="GROUP")
merge[1:5,]
outfile<-paste("TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_COLORS.txt",sep="")
write.table(merge,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
RowColors<-as.vector(merge$col1a)
outfile<-paste("TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"NOT_NORMALIZED.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-0.3,-0.05,length=100),  seq(-0.049,0.05,length=100),  seq(0.051,0.3,length=100))
pdf(outfile, width=5, height=15, pointsize=14)
heatmap.2(as.matrix(data),main = "TCGA sample level data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=FALSE, Rowv=FALSE, keysize=1)
dev.off()
Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(data)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes
write.table(Data_Ordered_Classes, "TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.txt", sep="\t")
ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=TYPE_ORDER) + scale_shape_manual(values=c(rep(95, NUMBER_TYPES))) + scale_size_manual(values=c(rep((20), NUMBER_TYPES))) + scale_color_manual(values=colors)
ggsave("TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.pdf", width = NUMBER_TYPES, height= 170, dpi=75, limitsize=FALSE)




data<-read.delim("GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")
cutoffheight<-3.4

annotation<-read.delim("GDSC_Novel_Sample_annotation.txt",header=TRUE,sep="\t")
annotation$CLASS<-NULL
annotation$TYPE<-annotation$TYPE_NOVEL
annotation$TYPE_NOVEL<-NULL

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
sample_groups<-cutree(hc,k=NULL,h=cutoffheight)
sample_groups<-as.data.frame(sample_groups)
sample_groups$LABEL<-row.names(sample_groups)
sample_groups[1:5,]
sample_groups$V2<-sample_groups$LABEL
alt_groups<-left_join(order,sample_groups,by="V2")
alt_groups<-select(alt_groups,sample_groups,LABEL)
colnames(alt_groups)<-c("GROUP","LABEL")
alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
alt_groups[1:5,]
outfile<-paste("GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_LABELS.txt",sep="")
write.table(alt_groups,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
k1<-max(alt_groups$GROUP)
k1
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col1a<-sample(color,k1,replace=TRUE)
col1b<-c(1:k1)
col1<-as.data.frame(cbind(col1a,col1b))
colnames(col1)<-c("col1a","GROUP")
col1$GROUP<-as.numeric(col1$GROUP)
merge<-left_join(alt_groups,col1,by="GROUP")
merge[1:5,]
outfile<-paste("GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_COLORS.txt",sep="")
write.table(merge,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
RowColors<-as.vector(merge$col1a)
outfile<-paste("GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"NOT_NORMALIZED.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-0.3,-0.05,length=100),  seq(-0.049,0.05,length=100),  seq(0.051,0.3,length=100))
pdf(outfile, width=5, height=15, pointsize=14)
heatmap.2(as.matrix(data),main = "GDSC sample level data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=FALSE, Rowv=FALSE, keysize=1)
dev.off()
Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(data)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes
write.table(Data_Ordered_Classes, "GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.txt", sep="\t")
ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=TYPE_ORDER) + scale_shape_manual(values=c(rep(95, NUMBER_TYPES))) + scale_size_manual(values=c(rep((20), NUMBER_TYPES))) + scale_color_manual(values=colors)
ggsave("GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.pdf", width = NUMBER_TYPES, height= 170, dpi=75, limitsize=FALSE)





data<-read.delim("CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")
cutoffheight<-4.2

annotation<-read.delim("CCLE_Novel_Sample_annotation.txt",header=TRUE,sep="\t")
annotation$CLASS<-NULL
annotation$TYPE<-annotation$TYPE_NOVEL
annotation$TYPE_NOVEL<-NULL

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
sample_groups<-cutree(hc,k=NULL,h=cutoffheight)
sample_groups<-as.data.frame(sample_groups)
sample_groups$LABEL<-row.names(sample_groups)
sample_groups[1:5,]
sample_groups$V2<-sample_groups$LABEL
alt_groups<-left_join(order,sample_groups,by="V2")
alt_groups<-select(alt_groups,sample_groups,LABEL)
colnames(alt_groups)<-c("GROUP","LABEL")
alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
alt_groups[1:5,]
outfile<-paste("CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_LABELS.txt",sep="")
write.table(alt_groups,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
k1<-max(alt_groups$GROUP)
k1
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col1a<-sample(color,k1,replace=TRUE)
col1b<-c(1:k1)
col1<-as.data.frame(cbind(col1a,col1b))
colnames(col1)<-c("col1a","GROUP")
col1$GROUP<-as.numeric(col1$GROUP)
merge<-left_join(alt_groups,col1,by="GROUP")
merge[1:5,]
outfile<-paste("CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"_SAMPLE_COLORS.txt",sep="")
write.table(merge,outfile, col.names=TRUE, row.names=TRUE, sep="\t")
RowColors<-as.vector(merge$col1a)
outfile<-paste("CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h",cutoffheight,"NOT_NORMALIZED.pdf",sep="")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
col_breaks <- c(seq(-0.3,-0.05,length=100),  seq(-0.049,0.05,length=100),  seq(0.051,0.3,length=100))
pdf(outfile, width=5, height=15, pointsize=14)
heatmap.2(as.matrix(data),main = "CCLE sample level data clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=FALSE, Rowv=FALSE, keysize=1)
dev.off()
Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(data)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes
write.table(Data_Ordered_Classes, "CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.txt", sep="\t")
ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=TYPE_ORDER) + scale_shape_manual(values=c(rep(95, NUMBER_TYPES))) + scale_size_manual(values=c(rep((20), NUMBER_TYPES))) + scale_color_manual(values=colors)
ggsave("CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_SIDEBAR_SIMPLIFIED.pdf", width = NUMBER_TYPES, height= 170, dpi=75, limitsize=FALSE)

