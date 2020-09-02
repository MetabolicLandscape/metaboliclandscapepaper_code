library(dplyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(BBmisc)


GEO_Annotation <- read.delim("GEO_Novel_Sample_annotation.txt",header=TRUE,sep="\t")
GEO_Annotation <- select(GEO_Annotation,"GSM_IDENTIFIER","TYPE")
colnames(GEO_Annotation)<-c("GSM_IDENTIFIER","TYPE")
tissue_classes<-c("Normal brain","Brain cancer - anaplastic astrocytoma","Brain cancer - anaplastic oligoastrocytoma","Brain cancer - anaplastic oligodendroglioma","Brain cancer - astrocytoma","Brain cancer - atypical teratoid / rhabdoid tumor","Brain cancer - ependymoma","Brain cancer - glioblastoma","Brain cancer - medulloblastoma","Brain cancer - meningioma","Brain cancer - oligoastrocytoma","Brain cancer - oligodendroglioma","Brain cancer - pilocytic astrocytoma","Primary central nervous system lymphoma","primitive neuroectodermal tumor","B lymphocytes","Multiple myeloma","Diffuse large B-Cell lymphoma","Follicular lymphoma","Mantle cell lymphoma","Acute lymphoblastic leukemia","Chronic lymphoblastic leukemia","Acute myeloid leukemia","Burkitt lymphoma","T-cell lymphoma","Normal breast","Breast cancer - TNBC","Breast cancer - ER-neg / HER2-pos","Breast cancer - ER-pos / HER2-pos","Breast cancer - ER-pos / HER2-neg","Normal colon_rectum","Colorectal cancer","Normal esophagus","Esophageal cancer - adenocarcinoma","Esophageal cancer - squamous cell carcinoma","Normal stomach","Gastric cancer","Normal liver","HCC","Cholangiocarcinoma","Periampullary cancer","Normal spleen","Normal lung","Lung cancer - Adenocarcinoma","Lung cancer - Neuroendocrine","Lung cancer - Squamous cell carcinoma","Normal pancreas","Pancreas cancer","Normal oral cavity","Normal parotid gland","Normal pharynx","Normal sinuses_nasal cavity","HNSCC","HNSCC - nasopharynx","Adrenal cancer - adrenocortical","Adrenal cancer - neuroblastoma","Adrenal cancer - pheochromocytoma","Renal cancer - chromophobe","Renal cancer - clear cell","Renal cancer - clear cell sarcoma","Renal cancer - nephroblastoma","Renal cancer - papillary","Normal skin","Melanoma - cutaneous","Melanoma - uveal","Normal heart","Normal muscle","Normal thyroid","Thyroid cancer - anaplastic","Thyroid cancer - follicular","Thyroid cancer - papillary","Bladder cancer","Normal prostate","Prostate cancer","Cervical cancer","Ovarian cancer","Vulva cancer","Ewing's sarcoma","Leiomyosarcoma","Liposarcoma","Myelodysplastic syndrome","Osteosarcoma","Synovial sarcoma","Sarcoma NOS","Undifferentiated sarcoma")
melanomas<-c("Melanoma - cutaneous","Melanoma - uveal","Na1","na2","na3","na4","na5","na6","na7","na8","na9","na10","na11","na12","na13")
brains<-c("Normal brain","Brain cancer - astrocytoma","Brain cancer - anaplastic astrocytoma","Brain cancer - pilocytic astrocytoma","Brain cancer - oligoastrocytoma","Brain cancer - anaplastic oligoastrocytoma","Brain cancer - oligodendroglioma","Brain cancer - anaplastic oligodendroglioma","Brain cancer - atypical teratoid / rhabdoid tumor","Brain cancer - ependymoma","Brain cancer - glioblastoma","Brain cancer - medulloblastoma","Brain cancer - meningioma","Na1","Na2")
breasts<-c("Normal breast","Breast cancer - TNBC","Breast cancer - ER-neg / HER2-pos","Breast cancer - ER-pos / HER2-pos","Breast cancer - ER-pos / HER2-neg","Na1","na2","na3","na4","na5","na6","na7","na8","na9","na10")
my_palette <- colorRampPalette(c("dodgerblue3", "white", "red2"))(n = 299)
colors<-c("yellowgreen","yellow3","wheat4","violetred3","turquoise4","tomato4","thistle4","tan4","steelblue4","springgreen4","snow4","slateblue4","skyblue4","sienna4","seagreen4","salmon4","rosybrown4","red4","purple4","plum4")

##clusterdivisions
row_labels<-read.delim("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_h6_SAMPLE_COLORS.txt",header=TRUE,sep="\t")
row_labels<-select(row_labels,LABEL,col1a)


infile2<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt",sep="")
outfile2<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_Samples_h6_MELANOMA.pdf",sep="")
outfile1<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_Samples_h6_BRAIN.pdf",sep="")
outfile3<-paste("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_Samples_h6_BREAST.pdf",sep="")
Data_Ordered<-read.delim(infile2,header=TRUE,sep="\t")

Data_Ordered_Classes<-as.data.frame(row.names(Data_Ordered))
colnames(Data_Ordered_Classes)<-c("GSM_IDENTIFIER")
Data_Ordered_Classes<-left_join(Data_Ordered_Classes,GEO_Annotation, by= "GSM_IDENTIFIER")
Data_Ordered_Classes$TYPE<-as.character(Data_Ordered_Classes$TYPE)
Data_Ordered_Classes$GSM_IDENTIFIER<-factor(Data_Ordered_Classes$GSM_IDENTIFIER, levels = unique(Data_Ordered_Classes$GSM_IDENTIFIER))

Filtering<-filter(GEO_Annotation,TYPE %in% melanomas)
Filtering_brain<-filter(GEO_Annotation,TYPE %in% brains)
Filtering_breast<-filter(GEO_Annotation,TYPE %in% breasts)
nrow(Filtering)
nrow(Filtering_brain)
nrow(Filtering_breast)


#filtering step



Data_Ordered$Sample<-row.names(Data_Ordered)
Data_Needed<-filter(Data_Ordered,Sample %in% Filtering$GSM_IDENTIFIER)
row.names(Data_Needed)<-Data_Needed$Sample
Samples<-as.data.frame(Data_Needed$Sample)
colnames(Samples)<-"LABEL"
Data_Needed$Sample<-NULL


##rowcolors
row_labels_needed<-left_join(Samples,row_labels,by="LABEL")
RowColors<-as.vector(row_labels_needed$col1a)
glimpse(RowColors)

col_breaks <- c(seq(-0.3,-0.05,length=100),  seq(-0.049,0.05,length=100),  seq(0.051,0.3,length=100))
pdf(outfile2, width=33.33, height=10, pointsize=14)
heatmap.2(as.matrix(Data_Needed),main = "GEO mixing matrix SOTA clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=TRUE, Rowv=FALSE, keysize=1)
dev.off()

Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(Data_Needed)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,GEO_Annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes


ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=melanomas) + scale_shape_manual(values=c(rep(95, length(melanomas)))) + scale_size_manual(values=c(rep((20), length(melanomas)))) + scale_color_manual(values=colors)

ggsave("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_Samples_h6_MELANOMA-sidebar.pdf", width = 15, height= 30, dpi=75, limitsize=FALSE)













Data_Ordered$Sample<-row.names(Data_Ordered)
Data_Needed<-filter(Data_Ordered,Sample %in% Filtering_brain$GSM_IDENTIFIER)
row.names(Data_Needed)<-Data_Needed$Sample
Samples<-as.data.frame(Data_Needed$Sample)
colnames(Samples)<-"LABEL"
Data_Needed$Sample<-NULL


##rowcolors
row_labels_needed<-left_join(Samples,row_labels,by="LABEL")
RowColors<-as.vector(row_labels_needed$col1a)



pdf(outfile1, width=33.33, height=15, pointsize=14)
heatmap.2(as.matrix(Data_Needed),main = "GEO mixing matrix SOTA clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=TRUE, Rowv=FALSE, keysize=1)
dev.off()

Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(Data_Needed)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,GEO_Annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes


ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=brains) + scale_shape_manual(values=c(rep(95, length(brains)))) + scale_size_manual(values=c(rep((20), length(brains)))) + scale_color_manual(values=colors)

ggsave("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_Samples_h6_BRAINS-sidebar.pdf", width = 15, height= 170, dpi=75, limitsize=FALSE)









Data_Ordered$Sample<-row.names(Data_Ordered)
Data_Needed<-filter(Data_Ordered,Sample %in% Filtering_breast$GSM_IDENTIFIER)
row.names(Data_Needed)<-Data_Needed$Sample
Samples<-as.data.frame(Data_Needed$Sample)
colnames(Samples)<-"LABEL"
Data_Needed$Sample<-NULL


##rowcolors
row_labels_needed<-left_join(Samples,row_labels,by="LABEL")
RowColors<-as.vector(row_labels_needed$col1a)

pdf(outfile3, width=33.33, height=15, pointsize=14)
heatmap.2(as.matrix(Data_Needed),main = "GEO mixing matrix SOTA clustered",  breaks=col_breaks, notecol="black",  density.info="none",  trace="none",   margins =c(5,5),   col=my_palette,  dendrogram="none", RowSideColors=RowColors, Colv=TRUE, Rowv=FALSE, keysize=1)
dev.off()

Genelevel_Ordered_Rownames_Classes<-NULL
Genelevel_Ordered_Rownames_Classes$Genes<-row.names(Data_Needed)
Genelevel_Ordered_Rownames_Classes<-as.data.frame(Genelevel_Ordered_Rownames_Classes)
colnames(Genelevel_Ordered_Rownames_Classes)<-"GSM_IDENTIFIER"
Genelevel_Ordered_Rownames_Classes<-left_join(Genelevel_Ordered_Rownames_Classes,GEO_Annotation, by= "GSM_IDENTIFIER")
Genelevel_Ordered_Rownames_Classes$TYPE<-as.character(Genelevel_Ordered_Rownames_Classes$TYPE)
Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER<-factor(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER, levels = unique(Genelevel_Ordered_Rownames_Classes$GSM_IDENTIFIER))
Data_Ordered_Classes<-Genelevel_Ordered_Rownames_Classes


ggplot(data = Data_Ordered_Classes) + geom_point(mapping = aes(y = GSM_IDENTIFIER, x = TYPE, shape=TYPE, color=TYPE, size=TYPE, fill=TYPE)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size=30)) + theme(axis.text.y = element_text(size=4)) + theme(legend.position="none") + scale_y_discrete(limits = rev(levels(Data_Ordered_Classes$GSM_IDENTIFIER))) + scale_x_discrete(limits=breasts) + scale_shape_manual(values=c(rep(95, length(breasts)))) + scale_size_manual(values=c(rep((20), length(breasts)))) + scale_color_manual(values=colors)

ggsave("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered_Samples_h6_BREASTS-sidebar.pdf", width = 15, height= 200, dpi=75, limitsize=FALSE)

