library(dplyr)

h<-5.6

annotation<-read.delim("GEO_Novel_Sample_annotation.txt", header=TRUE, sep="\t")
annotation<-select(annotation,GSM_IDENTIFIER,TYPE)
labelfile<-paste("GEO_MixingMatrix_Metabolic_Sel_1-cor_ward_clustered_h",h,"_SAMPLE_LABELS.txt",sep="")
labels<-read.delim(labelfile, header=TRUE, sep="\t")
colnames(labels)<-c("GROUP","GSM_IDENTIFIER")

annotatedgroups<-left_join(annotation,labels,by="GSM_IDENTIFIER")
annotatedgroups$TYPE<-as.factor(annotatedgroups$TYPE)

clustersummary<-matrix(nrow=length(unique(annotation$TYPE)),ncol=max(annotatedgroups$GROUP))
colnames(clustersummary)<-1:ncol(clustersummary)
rownames(clustersummary)<-sort(unique(annotation$TYPE))

for (i in 1:ncol(clustersummary)){
	clustercontent<-filter(annotatedgroups, GROUP == i)
	clustertable<-as.data.frame(table(clustercontent$TYPE))
	clustertable<-arrange(clustertable,Var1)
	clustersummary[,i]<-clustertable$Freq
}

outfile<-paste("GEO_MixingMatrix_Metabolic_Sel_1-cor_ward_clustered_h",h,"_CLUSTER_COMPOSITION_SUMMARY.txt",sep="")
write.table(clustersummary,outfile,col.names=TRUE,row.names=TRUE,sep="\t")

