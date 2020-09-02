library(dplyr)
##Combining data
Cibersort<-read.csv("tcga_output_cibersort.csv",header=TRUE,sep=",")
glimpse(Cibersort)
Cibersort$sample_id<-gsub("-",".",Cibersort$sample_id)
Mixmatrix<-read.delim("Consensus_Mix_Matrix_20170906_Duplicate_removed_TCGA_data_.txt",header=TRUE,sep="\t")
Mixmatrix<-as.data.frame(t(Mixmatrix))
Mixmatrix$sample_id<-row.names(Mixmatrix)
Combined<-left_join(Cibersort,Mixmatrix,by="sample_id")
glimpse(Combined)

ImmunoFractions<-colnames(Combined[2:23])
ImmunoFractions
TCs<-colnames(Combined[28:length(colnames(Combined))])
outfilepearson<-paste("Immuno_Focused_TCGA_TCs_pearson_correlations_ABSOLUTE.txt",sep="")
outfilespearman<-paste("Immuno_Focused_TCGA_TCs_spearman_correlations_ABSOLUTE.txt",sep="")

cortest<-c("IMMUNE_CELLTYPE","TC","Correlation","PValue")
for (i in 1:length(ImmunoFractions)){
  ith_immu_data<-select(Combined,ImmunoFractions[i],TCs)
  for (j in 1:length(TCs)){
    test_data<-select(ith_immu_data,ImmunoFractions[i],TCs[j])
    colnames(test_data)<-c("IMMUNE_CELLTYPE","TC")
    testresults<-cor.test(test_data$IMMUNE_CELLTYPE,test_data$TC,method="pearson")
    info<-c(ImmunoFractions[i],TCs[j],as.numeric(testresults["estimate"]),as.numeric(testresults["p.value"]))
    cortest<-rbind(cortest,info)
  }
}
write.table(cortest,outfilepearson,row.names=TRUE,col.names=TRUE,sep="\t")


cortest2<-c("IMMUNE_CELLTYPE","TC","Correlation","PValue")
for (i in 1:length(ImmunoFractions)){
  ith_immu_data<-select(Combined,ImmunoFractions[i],TCs)
  for (j in 1:length(TCs)){
    test_data<-select(ith_immu_data,ImmunoFractions[i],TCs[j])
    colnames(test_data)<-c("IMMUNE_CELLTYPE","TC")
    testresults<-cor.test(test_data$IMMUNE_CELLTYPE,test_data$TC,method="spearman")
    info<-c(ImmunoFractions[i],TCs[j],as.numeric(testresults["estimate"]),as.numeric(testresults["p.value"]))
    cortest2<-rbind(cortest2,info)
  }
}
write.table(cortest2,outfilespearman,row.names=TRUE,col.names=TRUE,sep="\t")




##Combining data
Cibersort<-read.csv("geo100k_output_cibersort.csv",header=TRUE,sep=",")
Cibersort[1:5,1:5]
Mixmatrix<-read.delim("Consensus_Mix_Matrix_20173006_GEO_with_normal_.txt",header=TRUE,sep="\t")
Mixmatrix<-as.data.frame(t(Mixmatrix))
Mixmatrix$sample_id<-row.names(Mixmatrix)
Combined<-left_join(Cibersort,Mixmatrix,by="sample_id")


ImmunoFractions<-colnames(Combined[2:23])
ImmunoFractions
TCs<-colnames(Combined[28:length(colnames(Combined))])
outfilepearson<-paste("Immuno_Focused_GEO_TCs_pearson_correlations_ABSOLUTE.txt",sep="")
outfilespearman<-paste("Immuno_Focused_GEO_TCs_spearman_correlations_ABSOLUTE.txt",sep="")

cortest<-c("IMMUNE_CELLTYPE","TC","Correlation","PValue")
for (i in 1:length(ImmunoFractions)){
  ith_immu_data<-select(Combined,ImmunoFractions[i],TCs)
  for (j in 1:length(TCs)){
    test_data<-select(ith_immu_data,ImmunoFractions[i],TCs[j])
    colnames(test_data)<-c("IMMUNE_CELLTYPE","TC")
    testresults<-cor.test(test_data$IMMUNE_CELLTYPE,test_data$TC,method="pearson")
    info<-c(ImmunoFractions[i],TCs[j],as.numeric(testresults["estimate"]),as.numeric(testresults["p.value"]))
    cortest<-rbind(cortest,info)
  }
}
write.table(cortest,outfilepearson,row.names=TRUE,col.names=TRUE,sep="\t")


cortest2<-c("IMMUNE_CELLTYPE","TC","Correlation","PValue")
for (i in 1:length(ImmunoFractions)){
  ith_immu_data<-select(Combined,ImmunoFractions[i],TCs)
  for (j in 1:length(TCs)){
    test_data<-select(ith_immu_data,ImmunoFractions[i],TCs[j])
    colnames(test_data)<-c("IMMUNE_CELLTYPE","TC")
    testresults<-cor.test(test_data$IMMUNE_CELLTYPE,test_data$TC,method="spearman")
    info<-c(ImmunoFractions[i],TCs[j],as.numeric(testresults["estimate"]),as.numeric(testresults["p.value"]))
    cortest2<-rbind(cortest2,info)
  }
}
write.table(cortest2,outfilespearman,row.names=TRUE,col.names=TRUE,sep="\t")


