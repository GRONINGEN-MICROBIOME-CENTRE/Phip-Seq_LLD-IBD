args = commandArgs(trailingOnly=TRUE)
args = as.integer(args)

pheno = read.table("pheno_data/pheno_4analysis2.txt")
dataMat_T = read.table("pheno_data/probes_4analysis_b1.txt")
covars = read.table("pheno_data/covars_4analysis2.txt")
info = read.table("info_4analysis.txt",sep="\t")

library(foreach)

result = foreach(i = args[2]:args[3],.combine = rbind)%:%
  foreach(j = 1:ncol(dataMat_T),.combine = rbind)%do%{
    print(i)
    c1 = summary(glm(dataMat_T[,j] ~ pheno[,i] + covars$Sex + covars$Age + covars$plate_id,family = "binomial"))
    if (rownames(c1$coef)[2] == "pheno[, i]"){
    data.frame(pheno = colnames(pheno)[i],
               probe = colnames(dataMat_T)[j],
               probe.protein = info[colnames(dataMat_T)[j],"prot"],
               probe.NA = sum(is.na(dataMat_T[,j])),
               probe.positive = sum(dataMat_T[,j],na.rm = T),
               probe.negative = sum(dataMat_T[,j]==0,na.rm = T),
               beta = c1$coef[2,1],se=c1$coef[2,2],F=c1$coef[2,3],P=c1$coef[2,4])
    
  }}
write.table(result,paste0("results.probes.Pheno/result.batch",args[1],".txt"),sep="\t",quote = F,row.names=F)
