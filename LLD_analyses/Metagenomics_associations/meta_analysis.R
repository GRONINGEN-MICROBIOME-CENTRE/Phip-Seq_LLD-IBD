library (meta)

LLD <- read.delim("LLD.probe.Microbiome.Nov2021.txt")
IBD <- read.delim("results_AB_microbiome_IBD_filtered_1_year.txt")

colnames(IBD)=paste("IBD",colnames(IBD), sep = "_")
colnames(LLD)=paste("LLD",colnames(LLD), sep = "_")

LLD$factor=paste(LLD$LLD_taxon,LLD$LLD_probe,sep = "_")
IBD$factor=paste(IBD$IBD_taxon,IBD$IBD_probe,sep = "_")

MicroMeta=merge(LLD,IBD,by="factor", all = T)
#25
MicroMeta$I2=NA
#26
MicroMeta$Q=NA
#27
MicroMeta$pval_Q=NA
#28
MicroMeta$z_fixed=NA
#29
MicroMeta$pval_fixed=NA
#30
MicroMeta$z_random=NA
#31
MicroMeta$pval_random=NA

for (i in 1:nrow(MicroMeta)){
  if(!is.na(MicroMeta[i,8]) & !is.na(MicroMeta[i,19])){
    mytest=metagen(TE=c(MicroMeta[i,8],MicroMeta[i,19]), seTE =c(MicroMeta[i,9],MicroMeta[i,20]), comb.fixed=TRUE, comb.random=TRUE) 
    MicroMeta[i,25]=mytest$I2
    MicroMeta[i,26]=mytest$Q
    MicroMeta[i,27]=mytest$pval.Q
    MicroMeta[i,28]=mytest$zval.fixed
    MicroMeta[i,29]=mytest$pval.fixed
    MicroMeta[i,30]=mytest$zval.random
    MicroMeta[i,31]=mytest$pval.random
  }
}

