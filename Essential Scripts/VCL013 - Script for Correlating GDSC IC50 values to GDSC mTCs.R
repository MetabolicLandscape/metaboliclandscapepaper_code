library(dplyr)
library(pspearman)
IC50data<-read.delim("GDSC_MixingMatrix_IC50_combined.txt",header=TRUE,sep="\t")
DRUG_IDENTIFIERS<-unique(IC50data$DRUG_ID)
DRUG_NAMES<-unique(IC50data$DRUG_NAME)
DRUG_NAMES<-as.character(DRUG_NAMES)
TCs<-colnames(IC50data[8:length(colnames(IC50data))])
outfilepearson<-paste("Drug_Focused_GDSC_TCs_pearson_correlations.txt",sep="")
outfilespearman<-paste("Drug_Focused_GDSC_TCs_spearman_correlations.txt",sep="")

cortest<-c("DRUG_ID","DRUG_NAME","TC","Correlation","PValue")
for (i in 1:length(DRUG_IDENTIFIERS)){
  ith_drug_IC50data<-filter(IC50data,DRUG_ID == DRUG_IDENTIFIERS[i])
  for (j in 1:length(TCs)){
    test_IC50data<-select(ith_drug_IC50data,IC50_VALUE,TCs[j])
    colnames(test_IC50data)<-c("IC50_VALUE","TC")
    testresults<-cor.test(test_IC50data$IC50_VALUE,test_IC50data$TC,method="pearson")
    info<-c(DRUG_IDENTIFIERS[i],DRUG_NAMES[i],TCs[j],as.numeric(testresults["estimate"]),as.numeric(testresults["p.value"]))
    cortest<-rbind(cortest,info)
  }
}
write.table(cortest,outfilepearson,row.names=TRUE,col.names=TRUE,sep="\t")


cortest2<-c("DRUG_ID","DRUG_NAME","TC","Correlation","PValue")
for (i in 1:length(DRUG_IDENTIFIERS)){
  ith_drug_IC50data<-filter(IC50data,DRUG_ID == DRUG_IDENTIFIERS[i])
  for (j in 1:length(TCs)){
    test_IC50data<-select(ith_drug_IC50data,IC50_VALUE,TCs[j])
    colnames(test_IC50data)<-c("IC50_VALUE","TC")
    testresults<-spearman.test(test_IC50data$IC50_VALUE,test_IC50data$TC,approximation="t-distribution")
    info<-c(DRUG_IDENTIFIERS[i],DRUG_NAMES[i],TCs[j],as.numeric(testresults["estimate"]),as.numeric(testresults["p.value"]))
    cortest2<-rbind(cortest2,info)
  }
}
write.table(cortest2,outfilespearman,row.names=TRUE,col.names=TRUE,sep="\t")