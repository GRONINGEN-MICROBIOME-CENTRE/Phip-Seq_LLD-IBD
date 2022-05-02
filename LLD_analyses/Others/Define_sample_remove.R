library(tidyverse)

#Some samples have really low antibody counts. Identify which samples ( Q_25 - IQR ). Rerun phenotypic assocaitions without those samples.

read_tsv("Data/Covariates_LLD&IBD.tsv") -> D_all
D_all %>% filter(grepl("LL", Sample_name)) -> D

summary(D$Probe_n) -> Probe_counts
IQR = Probe_counts[5] - Probe_counts[2] #639
Threshold = Probe_counts[2] - IQR
print(paste0("Minimal number of peptides to have ", as.character(Threshold)))
D %>% mutate(Remove = ifelse( Probe_n < Threshold, T, F)) -> D
D %>% group_by(Remove) %>% summarise( n() ) %>% print()
filter(D, Remove == T) %>% arrange(Probe_n) %>% print()
write_tsv(D, "Data/LLD_covariates.tsv") 

