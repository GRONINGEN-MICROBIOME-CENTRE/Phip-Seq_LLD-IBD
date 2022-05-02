library(reshape2)
library(vegan)
library(ape)
library(ggplot2)


Metadata <- read.delim("Metadata_IBD_cohort.txt")
final_selection= read.delim("final_list.txt")
Antibody <- read.csv("abprevalencetable_2815IBD.csv", row.names = 1)
taxa <- read.delim("metaphlan3_lld_base_fup_500fg_300ob_ibd_2713samples_merged_taxa_abundance.txt", row.names = 1)
ID <- read.delim("IDs.txt", header = F)
ID_IBD <- read.delim("ids_ibd.txt", header = T)
ID_IBD2=ID_IBD[!grepl("LLDeep",ID_IBD$PID),]

taxa2=as.data.frame(t(taxa))
taxa3=merge(ID_IBD2,taxa2, by.x="IDs_MGS", by.y="row.names") 
taxa3$IDs_MGS=NULL
taxa3$PF_RD=NULL
row.names(taxa3)=taxa3$PID
taxa3$PID=NULL
taxa3=as.data.frame(t(taxa3))
taxa4 <- as.data.frame(sapply(taxa3, as.numeric))
row.names(taxa4)=row.names(taxa3)
metaphlan = taxa4[,colSums(taxa4,na.rm = T)>0]
species = metaphlan[,grep("s__",colnames(metaphlan))]
species = cbind(species,metaphlan[,"UNKNOWN",drop =F])

do_clr_externalWeighting = function(interest_matrix, core_matrix){
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
	if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

metaphlan_transformed = do_clr_externalWeighting(metaphlan,species)
metaphlan_filt = metaphlan_transformed[,colSums(metaphlan>0)> 0.1 *nrow(metaphlan)]
metaphlan_filt = metaphlan_filt[,colnames(metaphlan_filt)!="UNKNOWN"]

Metadata_ids=merge(ID,Metadata, by.y = "Participant.ID", by.x = "V1")
IDs2=Metadata_ids[,c("V2","Sample_id_WIS")]

metaphlan_filt2=merge(IDs2,metaphlan_filt,by.x="V2", by.y="row.names")
metaphlan_filt2$V2=NULL
row.names(metaphlan_filt2)=metaphlan_filt2$Sample_id_WIS
metaphlan_filt2$Sample_id_WIS=NULL

mgs_ids=metaphlan_filt2[,1,drop=F]

#row.names(Antibody)=Antibody$Sample
#Antibody$Sample=NULL
Antibody2=merge(Antibody,mgs_ids,by="row.names")
row.names(Antibody2)=Antibody2$Row.names
Antibody2$Row.names=NULL
Antibody2$k__Archaea=NULL

mgs_ids=Antibody2[,1,drop=F]
metaphlan_filt3=merge(metaphlan_filt2,mgs_ids,by="row.names")
row.names(metaphlan_filt3)=metaphlan_filt3$Row.names
metaphlan_filt3$Row.names=NULL
metaphlan_filt3$agilent_121687=NULL


meta=Metadata_ids[,c("Sample_id_WIS","Age","Gender", "BMI", "Plate_id.1","Time.difference.serum.feces..days.")]
for (cc in colnames(meta)) {if (sum(is.na(meta[[cc]])) > 0) {meta[[cc]][is.na(meta[[cc]])] <- median(meta[[cc]],na.rm = T)}}
Antibody3=merge(meta,Antibody2,by.x="Sample_id_WIS", by.y="row.names")
row.names(Antibody3)=Antibody3$Sample_id_WIS
Antibody3$Sample_id_WIS=NULL

covars = Antibody3[,c(1:5)]
Antibody4 = Antibody3[,-c(1:5)]

info = read.table("./Desktop/antibody/info_4analysis.txt",sep="\t")
my_order=read.table("./Desktop/antibody/my_order.txt",sep="\t", header=T)

library(plyr)
info3=join(my_order,info, by="ID")



covars = covars
colnames(covars)=c("Age", "Sex", "BMI", "plate_id", "time")
covars$BMI=NULL
covars2=covars
covars=subset(covars,covars$time<365)

dataMat_T = subset(Antibody4, row.names(Antibody4) %in% row.names(covars))
pheno= subset(metaphlan_filt3, row.names(metaphlan_filt3) %in% row.names(covars))
info = info3


library(foreach)

result_1 = foreach(i = 1:100,.combine = rbind)%:%
  foreach(j = 1:ncol(dataMat_T),.combine = rbind)%do%{
    print(colnames(pheno)[i])
    c1 = summary(glm(dataMat_T[,j] ~ pheno[,i] + covars$Sex + covars$Age + covars$plate_id,family = "binomial"))
    if (rownames(c1$coef)[2] == "pheno[, i]"){
    data.frame(taxon = colnames(pheno)[i],
               probe = colnames(dataMat_T)[j],
               probe.protein = info[colnames(dataMat_T)[j],"prot"],
               probe.NA = sum(is.na(dataMat_T[,j])),
               probe.positive = sum(dataMat_T[,j],na.rm = T),
               probe.negative = sum(dataMat_T[,j]==0,na.rm = T),
               beta = c1$coef[2,1],se=c1$coef[2,2],F=c1$coef[2,3],P=c1$coef[2,4])
    
  }}
  
result_2 = foreach(i = 101:343,.combine = rbind)%:%
    foreach(j = 1:ncol(dataMat_T),.combine = rbind)%do%{
        print(colnames(pheno)[i])
        c1 = summary(glm(dataMat_T[,j] ~ pheno[,i] + covars$Sex + covars$Age + covars$plate_id,family = "binomial"))
        if (rownames(c1$coef)[2] == "pheno[, i]"){
            data.frame(taxon = colnames(pheno)[i],
                       probe = colnames(dataMat_T)[j],
                       probe.protein = info[colnames(dataMat_T)[j],"prot"],
                       probe.NA = sum(is.na(dataMat_T[,j])),
                       probe.positive = sum(dataMat_T[,j],na.rm = T),
                       probe.negative = sum(dataMat_T[,j]==0,na.rm = T),
                       beta = c1$coef[2,1],se=c1$coef[2,2],F=c1$coef[2,3],P=c1$coef[2,4])
            
        }}
        
        

results_3=rbind(results_1, results_2)

