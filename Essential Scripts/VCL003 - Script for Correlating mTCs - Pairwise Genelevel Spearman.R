library(dplyr)

GEO_data<-read.delim("GEO_Genelevel_REVISED_Metabolic_Sel.txt",header=TRUE,sep="\t")
GEO_data[1:3,1:3]
TCGA_data<-read.delim("TCGA_Genelevel_REVISED_Metabolic_Sel.txt",header=TRUE,sep="\t")
TCGA_data[1:3,1:3]
CCLE_data<-read.delim("CCLE_Genelevel_REVISED_Metabolic_Sel.txt",header=TRUE,sep="\t")
CCLE_data[1:3,1:3]
GDSC_data<-read.delim("GDSC_Genelevel_REVISED_Metabolic_Sel.txt",header=TRUE,sep="\t")
GDSC_data[1:3,1:3]

correlations<-data.frame(matrix(ncol = 7, nrow = 0))
x <- c("TC_GEO", "TC_TCGA", "ngenesA","ngenesB","overlap","cor", "pvalue")
colnames(correlations) <- x
for (i in 1:ncol(GEO_data)){
	correlations_set<-data.frame(matrix(ncol = 7, nrow = 0))
	x <- c("TC_GEO", "TC_TCGA", "ngenesA","ngenesB","overlap","cor", "pvalue")
	colnames(correlations_set) <- x
	GEO_ith_TC<-GEO_data[i]
	GEO_ith_TC<-subset(GEO_ith_TC,abs(GEO_ith_TC)>3)
	ngenesA<-nrow(GEO_ith_TC)
	genesA<-row.names(GEO_ith_TC)
	correlations_set[1:ncol(TCGA_data),1]<-colnames(GEO_data)[i]
	correlations_set[1:ncol(TCGA_data),3]<-ngenesA
        for (j in 1:ncol(TCGA_data)){
		TCGA_jth_TC<-TCGA_data[j]
		TCGA_jth_TC<-subset(TCGA_jth_TC,abs(TCGA_jth_TC)>3)
		ngenesB<-nrow(TCGA_jth_TC)
		genesB<-row.names(TCGA_jth_TC)
		genes_overlap<-intersect(genesB,genesA)
		correlations_set[j,2]<-colnames(TCGA_data)[j]
		correlations_set[j,4]<-ngenesB
		GEO_ith_TC$gene_names<-row.names(GEO_ith_TC)
		TCGA_jth_TC$gene_names<-row.names(TCGA_jth_TC)
		compareset<-inner_join(TCGA_jth_TC,GEO_ith_TC,by="gene_names")
		overlap<-nrow(compareset)
		correlations_set[j,5]<-overlap
		if (overlap>2){
			testresults<-cor.test(as.matrix(compareset[1]),as.matrix(compareset[3]),method="spearman")
         		correlations_set[j,6]<-as.numeric(testresults["estimate"])
			correlations_set[j,7]<-as.numeric(testresults["p.value"])
		}
                
        }
	correlations<-rbind(correlations,correlations_set)
	tail(correlations)
}
write.table(correlations,"GEO_TCGA_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",row.names=FALSE,col.names=TRUE,sep="\t")


correlations<-data.frame(matrix(ncol = 7, nrow = 0))
x <- c("TC_GEO", "TC_CCLE", "ngenesA","ngenesB","overlap","cor", "pvalue")
colnames(correlations) <- x
for (i in 1:ncol(GEO_data)){
	correlations_set<-data.frame(matrix(ncol = 7, nrow = 0))
	x <- c("TC_GEO", "TC_CCLE", "ngenesA","ngenesB","overlap","cor", "pvalue")
	colnames(correlations_set) <- x
	GEO_ith_TC<-GEO_data[i]
	GEO_ith_TC<-subset(GEO_ith_TC,abs(GEO_ith_TC)>3)
	ngenesA<-nrow(GEO_ith_TC)
	genesA<-row.names(GEO_ith_TC)
	correlations_set[1:ncol(CCLE_data),1]<-colnames(GEO_data)[i]
	correlations_set[1:ncol(CCLE_data),3]<-ngenesA
        for (j in 1:ncol(CCLE_data)){
		CCLE_jth_TC<-CCLE_data[j]
		CCLE_jth_TC<-subset(CCLE_jth_TC,abs(CCLE_jth_TC)>3)
		ngenesB<-nrow(CCLE_jth_TC)
		genesB<-row.names(CCLE_jth_TC)
		genes_overlap<-intersect(genesB,genesA)
		correlations_set[j,2]<-colnames(CCLE_data)[j]
		correlations_set[j,4]<-ngenesB
		GEO_ith_TC$gene_names<-row.names(GEO_ith_TC)
		CCLE_jth_TC$gene_names<-row.names(CCLE_jth_TC)
		compareset<-inner_join(CCLE_jth_TC,GEO_ith_TC,by="gene_names")
		overlap<-nrow(compareset)
		correlations_set[j,5]<-overlap
		if (overlap>2){
			testresults<-cor.test(as.matrix(compareset[1]),as.matrix(compareset[3]),method="spearman")
         		correlations_set[j,6]<-as.numeric(testresults["estimate"])
			correlations_set[j,7]<-as.numeric(testresults["p.value"])
		}
                
        }
	correlations<-rbind(correlations,correlations_set)
	tail(correlations)
}
write.table(correlations,"GEO_CCLE_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",row.names=FALSE,col.names=TRUE,sep="\t")




correlations<-data.frame(matrix(ncol = 7, nrow = 0))
x <- c("TC_GEO", "TC_GDSC", "ngenesA","ngenesB","overlap","cor", "pvalue")
colnames(correlations) <- x
for (i in 1:ncol(GEO_data)){
	correlations_set<-data.frame(matrix(ncol = 7, nrow = 0))
	x <- c("TC_GEO", "TC_GDSC", "ngenesA","ngenesB","overlap","cor", "pvalue")
	colnames(correlations_set) <- x
	GEO_ith_TC<-GEO_data[i]
	GEO_ith_TC<-subset(GEO_ith_TC,abs(GEO_ith_TC)>3)
	ngenesA<-nrow(GEO_ith_TC)
	genesA<-row.names(GEO_ith_TC)
	correlations_set[1:ncol(GDSC_data),1]<-colnames(GEO_data)[i]
	correlations_set[1:ncol(GDSC_data),3]<-ngenesA
        for (j in 1:ncol(GDSC_data)){
		GDSC_jth_TC<-GDSC_data[j]
		GDSC_jth_TC<-subset(GDSC_jth_TC,abs(GDSC_jth_TC)>3)
		ngenesB<-nrow(GDSC_jth_TC)
		genesB<-row.names(GDSC_jth_TC)
		genes_overlap<-intersect(genesB,genesA)
		correlations_set[j,2]<-colnames(GDSC_data)[j]
		correlations_set[j,4]<-ngenesB
		GEO_ith_TC$gene_names<-row.names(GEO_ith_TC)
		GDSC_jth_TC$gene_names<-row.names(GDSC_jth_TC)
		compareset<-inner_join(GDSC_jth_TC,GEO_ith_TC,by="gene_names")
		overlap<-nrow(compareset)
		correlations_set[j,5]<-overlap
		if (overlap>3){
			compareset[1:2,1:3]
			testresults<-cor.test(as.matrix(compareset[1]),as.matrix(compareset[3]),method="spearman")
            		correlations_set[j,6]<-as.numeric(testresults["estimate"])
			correlations_set[j,7]<-as.numeric(testresults["p.value"])
		}
                
        }
	correlations<-rbind(correlations,correlations_set)
	tail(correlations)
}
write.table(correlations,"GEO_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",row.names=FALSE,col.names=TRUE,sep="\t")



correlations<-data.frame(matrix(ncol = 7, nrow = 0))
x <- c("TC_TCGA", "TC_GDSC", "ngenesA","ngenesB","overlap","cor", "pvalue")
colnames(correlations) <- x
for (i in 1:ncol(TCGA_data)){
	correlations_set<-data.frame(matrix(ncol = 7, nrow = 0))
	x <- c("TC_TCGA", "TC_GDSC", "ngenesA","ngenesB","overlap","cor", "pvalue")
	colnames(correlations_set) <- x
	TCGA_ith_TC<-TCGA_data[i]
	TCGA_ith_TC<-subset(TCGA_ith_TC,abs(TCGA_ith_TC)>3)
	ngenesA<-nrow(TCGA_ith_TC)
	genesA<-row.names(TCGA_ith_TC)
	correlations_set[1:ncol(GDSC_data),1]<-colnames(TCGA_data)[i]
	correlations_set[1:ncol(GDSC_data),3]<-ngenesA
        for (j in 1:ncol(GDSC_data)){
		GDSC_jth_TC<-GDSC_data[j]
		GDSC_jth_TC<-subset(GDSC_jth_TC,abs(GDSC_jth_TC)>3)
		ngenesB<-nrow(GDSC_jth_TC)
		genesB<-row.names(GDSC_jth_TC)
		genes_overlap<-intersect(genesB,genesA)
		correlations_set[j,2]<-colnames(GDSC_data)[j]
		correlations_set[j,4]<-ngenesB
		TCGA_ith_TC$gene_names<-row.names(TCGA_ith_TC)
		GDSC_jth_TC$gene_names<-row.names(GDSC_jth_TC)
		compareset<-inner_join(GDSC_jth_TC,TCGA_ith_TC,by="gene_names")
		overlap<-nrow(compareset)
		correlations_set[j,5]<-overlap
		if (overlap>3){
			compareset[1:2,1:3]
			testresults<-cor.test(as.matrix(compareset[1]),as.matrix(compareset[3]),method="spearman")
            		correlations_set[j,6]<-as.numeric(testresults["estimate"])
			correlations_set[j,7]<-as.numeric(testresults["p.value"])
		}
                
        }
	correlations<-rbind(correlations,correlations_set)
	tail(correlations)
}
write.table(correlations,"TCGA_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",row.names=FALSE,col.names=TRUE,sep="\t")



correlations<-data.frame(matrix(ncol = 7, nrow = 0))
x <- c("TC_TCGA", "TC_CCLE", "ngenesA","ngenesB","overlap","cor", "pvalue")
colnames(correlations) <- x
for (i in 1:ncol(TCGA_data)){
	correlations_set<-data.frame(matrix(ncol = 7, nrow = 0))
	x <- c("TC_TCGA", "TC_CCLE", "ngenesA","ngenesB","overlap","cor", "pvalue")
	colnames(correlations_set) <- x
	TCGA_ith_TC<-TCGA_data[i]
	TCGA_ith_TC<-subset(TCGA_ith_TC,abs(TCGA_ith_TC)>3)
	ngenesA<-nrow(TCGA_ith_TC)
	genesA<-row.names(TCGA_ith_TC)
	correlations_set[1:ncol(CCLE_data),1]<-colnames(TCGA_data)[i]
	correlations_set[1:ncol(CCLE_data),3]<-ngenesA
        for (j in 1:ncol(CCLE_data)){
		CCLE_jth_TC<-CCLE_data[j]
		CCLE_jth_TC<-subset(CCLE_jth_TC,abs(CCLE_jth_TC)>3)
		ngenesB<-nrow(CCLE_jth_TC)
		genesB<-row.names(CCLE_jth_TC)
		genes_overlap<-intersect(genesB,genesA)
		correlations_set[j,2]<-colnames(CCLE_data)[j]
		correlations_set[j,4]<-ngenesB
		TCGA_ith_TC$gene_names<-row.names(TCGA_ith_TC)
		CCLE_jth_TC$gene_names<-row.names(CCLE_jth_TC)
		compareset<-inner_join(CCLE_jth_TC,TCGA_ith_TC,by="gene_names")
		overlap<-nrow(compareset)
		correlations_set[j,5]<-overlap
		if (overlap>3){
			compareset[1:2,1:3]
			testresults<-cor.test(as.matrix(compareset[1]),as.matrix(compareset[3]),method="spearman")
            		correlations_set[j,6]<-as.numeric(testresults["estimate"])
			correlations_set[j,7]<-as.numeric(testresults["p.value"])
		}
                
        }
	correlations<-rbind(correlations,correlations_set)
	tail(correlations)
}
write.table(correlations,"TCGA_CCLE_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",row.names=FALSE,col.names=TRUE,sep="\t")



correlations<-data.frame(matrix(ncol = 7, nrow = 0))
x <- c("TC_CCLE", "TC_GDSC", "ngenesA","ngenesB","overlap","cor", "pvalue")
colnames(correlations) <- x
for (i in 1:ncol(CCLE_data)){
	correlations_set<-data.frame(matrix(ncol = 7, nrow = 0))
	x <- c("TC_CCLE", "TC_GDSC", "ngenesA","ngenesB","overlap","cor", "pvalue")
	colnames(correlations_set) <- x
	CCLE_ith_TC<-CCLE_data[i]
	CCLE_ith_TC<-subset(TCGA_ith_TC,abs(CCLE_ith_TC)>3)
	ngenesA<-nrow(CCLE_ith_TC)
	genesA<-row.names(CCLE_ith_TC)
	correlations_set[1:ncol(GDSC_data),1]<-colnames(CCLE_data)[i]
	correlations_set[1:ncol(GDSC_data),3]<-ngenesA
        for (j in 1:ncol(GDSC_data)){
		GDSC_jth_TC<-GDSC_data[j]
		GDSC_jth_TC<-subset(GDSC_jth_TC,abs(GDSC_jth_TC)>3)
		ngenesB<-nrow(GDSC_jth_TC)
		genesB<-row.names(GDSC_jth_TC)
		genes_overlap<-intersect(genesB,genesA)
		correlations_set[j,2]<-colnames(GDSC_data)[j]
		correlations_set[j,4]<-ngenesB
		CCLE_ith_TC$gene_names<-row.names(CCLE_ith_TC)
		GDSC_jth_TC$gene_names<-row.names(GDSC_jth_TC)
		compareset<-inner_join(GDSC_jth_TC,CCLE_ith_TC,by="gene_names")
		overlap<-nrow(compareset)
		correlations_set[j,5]<-overlap
		if (overlap>3){
			compareset[1:2,1:3]
			testresults<-cor.test(as.matrix(compareset[1]),as.matrix(compareset[3]),method="spearman")
            		correlations_set[j,6]<-as.numeric(testresults["estimate"])
			correlations_set[j,7]<-as.numeric(testresults["p.value"])
		}
                
        }
	correlations<-rbind(correlations,correlations_set)
	tail(correlations)
}
write.table(correlations,"CCLE_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",row.names=FALSE,col.names=TRUE,sep="\t")


GEO<-GEO_data
TCGA<-TCGA_data
CCLE<-CCLE_data
GDSC<-GDSC_data

GEO_IDs<-row.names(GEO)
length(GEO_IDs)
is.numeric(GEO_IDs)
TCGA_IDs<-row.names(TCGA)
length(TCGA_IDs)
GDSC_IDs<-row.names(GDSC)
length(GDSC_IDs)
is.numeric(GDSC_IDs)
CCLE_IDs<-row.names(CCLE)
length(CCLE_IDs)

A<-GEO_IDs
B<-TCGA_IDs
C<-GDSC_IDs
D<-CCLE_IDs
repeats<-10000

check<-read.delim("GEO_TCGA_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",header=TRUE,sep="\t")
check$overlap_significance<-NA
check$overlap_percentage<-NA
for (i in 1:nrow(check)){
	counter<-0
	ngenesA<-check[i,3]
	ngenesB<-check[i,4]
	for (j in 1:repeats){		
		genesA<-sample(A,ngenesA, replace=FALSE)
		genesB<-sample(B,ngenesB, replace=FALSE)
		genes_overlap<-intersect(genesB,genesA)
		overlap<-length(genes_overlap)
		if(length(genes_overlap)>check[i,5]){
			counter<-counter+1
		}
	}
	signif<-(counter/repeats)
	check[i,8]<-signif
	check[i,9]<-((check[i,5]/min(check[i,3:4]))*100)
}

write.table(check,"GEO_TCGA_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations_CHECKED.txt",col.names=TRUE,sep="\t")



check<-read.delim("GEO_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",header=TRUE,sep="\t")
check$overlap_significance<-NA
check$overlap_percentage<-NA
for (i in 1:nrow(check)){
	counter<-0
	ngenesA<-check[i,3]
	ngenesB<-check[i,4]
	for (j in 1:repeats){		
		genesA<-sample(A,ngenesA, replace=FALSE)
		genesB<-sample(C,ngenesB, replace=FALSE)
		genes_overlap<-intersect(genesB,genesA)
		overlap<-length(genes_overlap)
		if(length(genes_overlap)>check[i,5]){
			counter<-counter+1
		}
	}
	signif<-(counter/repeats)
	check[i,8]<-signif
	check[i,9]<-((check[i,5]/min(check[i,3:4]))*100)
}

write.table(check,"GEO_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations_CHECKED.txt",col.names=TRUE,sep="\t")


check<-read.delim("GEO_CCLE_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",header=TRUE,sep="\t")
check$overlap_significance<-NA
check$overlap_percentage<-NA
for (i in 1:nrow(check)){
	counter<-0
	ngenesA<-check[i,3]
	ngenesB<-check[i,4]
	for (j in 1:repeats){		
		genesA<-sample(A,ngenesA, replace=FALSE)
		genesB<-sample(D,ngenesB, replace=FALSE)
		genes_overlap<-intersect(genesB,genesA)
		overlap<-length(genes_overlap)
		if(length(genes_overlap)>check[i,5]){
			counter<-counter+1
		}
	}
	signif<-(counter/repeats)
	check[i,8]<-signif
	check[i,9]<-((check[i,5]/min(check[i,3:4]))*100)
}

write.table(check,"GEO_CCLE_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations_CHECKED.txt",col.names=TRUE,sep="\t")


check<-read.delim("TCGA_CCLE_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",header=TRUE,sep="\t")
check$overlap_significance<-NA
check$overlap_percentage<-NA
for (i in 1:nrow(check)){
	counter<-0
	ngenesA<-check[i,3]
	ngenesB<-check[i,4]
	for (j in 1:repeats){		
		genesA<-sample(B,ngenesA, replace=FALSE)
		genesB<-sample(D,ngenesB, replace=FALSE)
		genes_overlap<-intersect(genesB,genesA)
		overlap<-length(genes_overlap)
		if(length(genes_overlap)>check[i,5]){
			counter<-counter+1
		}
	}
	signif<-(counter/repeats)
	check[i,8]<-signif
	check[i,9]<-((check[i,5]/min(check[i,3:4]))*100)
}

write.table(check,"TCGA_CCLE_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations_CHECKED.txt",col.names=TRUE,sep="\t")


check<-read.delim("TCGA_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",header=TRUE,sep="\t")
check$overlap_significance<-NA
check$overlap_percentage<-NA
for (i in 1:nrow(check)){
	counter<-0
	ngenesA<-check[i,3]
	ngenesB<-check[i,4]
	for (j in 1:repeats){		
		genesA<-sample(B,ngenesA, replace=FALSE)
		genesB<-sample(C,ngenesB, replace=FALSE)
		genes_overlap<-intersect(genesB,genesA)
		overlap<-length(genes_overlap)
		if(length(genes_overlap)>check[i,5]){
			counter<-counter+1
		}
	}
	signif<-(counter/repeats)
	check[i,8]<-signif
	check[i,9]<-((check[i,5]/min(check[i,3:4]))*100)
}

write.table(check,"TCGA_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations_CHECKED.txt",col.names=TRUE,sep="\t")


check<-read.delim("CCLE_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations.txt",header=TRUE,sep="\t")
check$overlap_significance<-NA
check$overlap_percentage<-NA
for (i in 1:nrow(check)){
	counter<-0
	ngenesA<-check[i,3]
	ngenesB<-check[i,4]
	for (j in 1:repeats){		
		genesA<-sample(D,ngenesA, replace=FALSE)
		genesB<-sample(C,ngenesB, replace=FALSE)
		genes_overlap<-intersect(genesB,genesA)
		overlap<-length(genes_overlap)
		if(length(genes_overlap)>check[i,5]){
			counter<-counter+1
		}
	}
	signif<-(counter/repeats)
	check[i,8]<-signif
	check[i,9]<-((check[i,5]/min(check[i,3:4]))*100)
}

write.table(check,"CCLE_GDSC_REVISED_mTCs_genelevel_filtered_SD3_spearman_correlations_CHECKED.txt",col.names=TRUE,sep="\t")
		