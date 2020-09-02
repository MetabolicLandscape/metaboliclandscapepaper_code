library(dplyr)

data<-read.delim("GEO_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
for(h in seq(from=0.2, to=10, by=0.2)){
	groups<-cutree(hc,k=NULL,h=h)
	outfile<-paste("GEO_MM_REVISED_samples_h",h,"_histogram_groups.txt",sep="")
	max(groups)
	groups<-as.data.frame(groups)
	groups$V2<-row.names(groups)
	alt_groups<-left_join(order,groups,by="V2")
	alt_groups<-select(alt_groups,groups,V2)
	colnames(alt_groups)<-c("GROUP","LABEL")
	alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
	quant<-as.data.frame(table(alt_groups$GROUP))
	quant[1:5,]
	is.numeric(quant$Freq)
	write.table(quant,outfile,sep="\t")
}


data<-read.delim("TCGA_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")
sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
for(h in seq(from=0.2, to=10, by=0.2)){
	groups<-cutree(hc,k=NULL,h=h)
	outfile<-paste("TCGA_MM_REVISED_samples_h",h,"_histogram_groups.txt",sep="")
	max(groups)
	groups<-as.data.frame(groups)
	groups$V2<-row.names(groups)
	alt_groups<-left_join(order,groups,by="V2")
	alt_groups<-select(alt_groups,groups,V2)
	colnames(alt_groups)<-c("GROUP","LABEL")
	alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
	quant<-as.data.frame(table(alt_groups$GROUP))
	quant[1:5,]
	is.numeric(quant$Freq)
	write.table(quant,outfile,sep="\t")
}


data<-read.delim("GDSC_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
for(h in seq(from=0.2, to=10, by=0.2)){
	groups<-cutree(hc,k=NULL,h=h)
	outfile<-paste("GDSC_MM_REVISED_samples_h",h,"_histogram_groups.txt",sep="")
	max(groups)
	groups<-as.data.frame(groups)
	groups$V2<-row.names(groups)
	alt_groups<-left_join(order,groups,by="V2")
	alt_groups<-select(alt_groups,groups,V2)
	colnames(alt_groups)<-c("GROUP","LABEL")
	alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
	quant<-as.data.frame(table(alt_groups$GROUP))
	quant[1:5,]
	is.numeric(quant$Freq)
	write.table(quant,outfile,sep="\t")
}


data<-read.delim("CCLE_MM_REVISED_Metabolic_Sel_1-cor_ward_clustered.txt", header=TRUE, sep="\t")

sampledist <- as.dist(1-cor(t(as.matrix(data))))
hc <- hclust(sampledist, method="ward.D2")
order<-as.data.frame(hc$order)
order[,2]<-hc$labels[order[,1]]
order[1:5,]
for(h in seq(from=0.2, to=10, by=0.2)){
	groups<-cutree(hc,k=NULL,h=h)
	outfile<-paste("CCLE_MM_REVISED_samples_h",h,"_histogram_groups.txt",sep="")
	max(groups)
	groups<-as.data.frame(groups)
	groups$V2<-row.names(groups)
	alt_groups<-left_join(order,groups,by="V2")
	alt_groups<-select(alt_groups,groups,V2)
	colnames(alt_groups)<-c("GROUP","LABEL")
	alt_groups$GROUP<-as.numeric(alt_groups$GROUP)
	quant<-as.data.frame(table(alt_groups$GROUP))
	quant[1:5,]
	is.numeric(quant$Freq)
	write.table(quant,outfile,sep="\t")
}

