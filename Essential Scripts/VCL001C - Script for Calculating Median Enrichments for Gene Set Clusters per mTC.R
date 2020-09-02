library(dplyr)

GSEAdata<-read.delim("GEO_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t")
ncolGSEA<-ncol(GSEAdata)
groups<-as.data.frame(read.delim("ALL_DATASETS_GENESET_CLUSTER_GROUPS.txt",header=TRUE,sep="\t"))
groups$GENESET<-row.names(groups)
colnames(groups)<-c("GROUP","GENESET")
GSEAdata$GENESET<-row.names(GSEAdata)
combined<-left_join(GSEAdata,groups,by="GENESET")
groups_vector<-unique(combined$GROUP)

medians<-matrix(nrow=length(groups_vector),ncol=ncolGSEA)
colnames(medians)<-colnames(GSEAdata)[1:ncolGSEA]
row.names(medians)<-1:nrow(medians)
for(j in 1:length(groups_vector)){
	row.names(medians)[j]<-paste("GROUP_",groups_vector[j],"_Median",sep="")
}

for(i in 1:ncolGSEA){
	for (j in 1:length(groups_vector)){
		GSEA_filter_groupj<-filter(combined, GROUP == groups_vector[j])
		medians[j,i]<-median(GSEA_filter_groupj[,i])
	}
}
medians<-t(medians)
write.table(medians,"GEO_GSEA_Geneset_Clusters_MEDIANS.txt",col.names=TRUE,row.names=TRUE,sep="\t")



GSEAdata<-read.delim("TCGA_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t")
ncolGSEA<-ncol(GSEAdata)
groups<-read.delim("ALL_DATASETS_GENESET_CLUSTER_GROUPS.txt",header=TRUE,sep="\t")
groups$GENESET<-row.names(groups)
colnames(groups)<-c("GROUP","GENESET")
GSEAdata$GENESET<-row.names(GSEAdata)
combined<-left_join(GSEAdata,groups,by="GENESET")
groups_vector<-unique(combined$GROUP)
medians<-matrix(nrow=length(groups_vector),ncol=ncolGSEA)
colnames(medians)<-colnames(GSEAdata)[1:ncolGSEA]
row.names(medians)<-1:nrow(medians)
for(j in 1:length(groups_vector)){
	row.names(medians)[j]<-paste("GROUP_",groups_vector[j],"_Median",sep="")
}

for(i in 1:ncolGSEA){
	for (j in 1:length(groups_vector)){
		GSEA_filter_groupj<-filter(combined, GROUP == groups_vector[j])
		medians[j,i]<-median(GSEA_filter_groupj[,i])
	}
}
medians<-t(medians)
write.table(medians,"TCGA_GSEA_Geneset_Clusters_MEDIANS.txt",col.names=TRUE,row.names=TRUE,sep="\t")


GSEAdata<-read.delim("CCLE_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t")
ncolGSEA<-ncol(GSEAdata)
groups<-read.delim("ALL_DATASETS_GENESET_CLUSTER_GROUPS.txt",header=TRUE,sep="\t")
groups$GENESET<-row.names(groups)
colnames(groups)<-c("GROUP","GENESET")
GSEAdata$GENESET<-row.names(GSEAdata)
combined<-left_join(GSEAdata,groups,by="GENESET")
groups_vector<-unique(combined$GROUP)
medians<-matrix(nrow=length(groups_vector),ncol=ncolGSEA)
colnames(medians)<-colnames(GSEAdata)[1:ncolGSEA]
row.names(medians)<-1:nrow(medians)
for(j in 1:length(groups_vector)){
	row.names(medians)[j]<-paste("GROUP_",groups_vector[j],"_Median",sep="")
}

for(i in 1:ncolGSEA){
	for (j in 1:length(groups_vector)){
		GSEA_filter_groupj<-filter(combined, GROUP == groups_vector[j])
		medians[j,i]<-median(GSEA_filter_groupj[,i])
	}
}
medians<-t(medians)
write.table(medians,"CCLE_GSEA_Geneset_Clusters_MEDIANS.txt",col.names=TRUE,row.names=TRUE,sep="\t")



GSEAdata<-read.delim("GDSC_GSEA_wardD2_absdist_ALLDATASETS_GENESET_CLUSTERED.txt",header=TRUE,sep="\t")
ncolGSEA<-ncol(GSEAdata)
groups<-read.delim("ALL_DATASETS_GENESET_CLUSTER_GROUPS.txt",header=TRUE,sep="\t")
groups$GENESET<-row.names(groups)
colnames(groups)<-c("GROUP","GENESET")
GSEAdata$GENESET<-row.names(GSEAdata)
combined<-left_join(GSEAdata,groups,by="GENESET")
groups_vector<-unique(combined$GROUP)
medians<-matrix(nrow=length(groups_vector),ncol=ncolGSEA)
colnames(medians)<-colnames(GSEAdata)[1:ncolGSEA]
row.names(medians)<-1:nrow(medians)
for(j in 1:length(groups_vector)){
	row.names(medians)[j]<-paste("GROUP_",groups_vector[j],"_Median",sep="")
}

for(i in 1:ncolGSEA){
	for (j in 1:length(groups_vector)){
		GSEA_filter_groupj<-filter(combined, GROUP == groups_vector[j])
		medians[j,i]<-median(GSEA_filter_groupj[,i])
	}
}
medians<-t(medians)
write.table(medians,"GDSC_GSEA_Geneset_Clusters_MEDIANS.txt",col.names=TRUE,row.names=TRUE,sep="\t")
