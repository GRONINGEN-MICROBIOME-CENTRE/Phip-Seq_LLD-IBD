library(tidyverse)
library(foreach)
print("Reading data")
Quant = read.table("Result.tsv") %>% as_tibble()
Cov = read_tsv("../../Data/Covariates_LLD&IBD.tsv")
AB = readRDS("../../Data/Immuno_matrix_postSelection.rds")
Bugs = read_tsv("Data/bacteria_associated_phenos.txt", col_names=F)
print("Formatting Shortbred results")
colnames(Quant) = Quant[1,]
Quant %>% filter(! ID == "ID") -> Quant
Quant %>%  filter(Family %in% Bugs$X1) -> Quant
Quant %>% select(ID, Family, Count) %>% spread(Family, Count) -> Quant_spred
print("Formatting AB matrix")
AB %>% filter(grepl("32_", ID)) -> ID
AB$ID = sapply(AB$ID, function(x){ str_split(x, "_")[[1]][2] } )

left_join(AB, select(Cov, c(ID, Sample_name))) %>% mutate(ID = Sample_name) %>% select(-Sample_name) -> AB


AB %>% select(one_of(colnames(Quant_spred))) -> AB

print("Formatting Covariates")
Cov %>% mutate(ID = Sample_name) %>% select(-Sample_name) -> Cov
print(dim(AB))
print(dim(Quant_spred))
print("running associations")
foreach( ab = colnames(select(AB, -ID)), .combine='rbind') %do% {
	AB %>% select(c("ID",ab)) -> I1
	Quant_spred  %>% select(c("ID",ab))  -> I2
	left_join(I2, I1, by="ID", suffix=c("_shortbred", "_ab")) %>% drop_na() -> I3
	left_join(I3, Cov, by="ID") -> I4
	
	Prevalence = select(I4, paste0(ab, "_ab")) %>% sum()
	Prevalence = Prevalence/dim(I4)[1]
	Mean_abundance = select(I4, paste0(ab, "_shortbred")) %>% as_vector() %>% as.numeric() %>% sum()	
	Mean_abundance = Mean_abundance/dim(I4)[1]

	Model = as.formula( paste0( "as.numeric(", ab,"_shortbred) ~", ab, "_ab", " + Age + Sex + plate_id") )
	summary(lm(Model, I4)) -> Result
	as.data.frame(Result$coefficients)[2,]  %>% as_tibble() %>% mutate(Peptide = ab, Prevalence = Prevalence, Mean_abundance =Mean_abundance) -> Result
	return(Result)
} -> Associations

Associations %>% arrange(`Pr(>|t|)`) -> Associations
Associations %>% drop_na() %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) -> Associations 
Associations %>% print()
write_tsv(Associations, "Summary_Stats.tsv")

