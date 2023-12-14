########################################################################
##Multivariate and univariate associations between PhIP-Seq and OMICS###
########################################################################
#v3 Sergio Andreu-Sanchez
#This script does the multi-omics analyses. This include: HAllA, sCCA and lasso prediction. 
#Check other scripts for: CMV analysis, rhinoviral analysis, scRNA-Seq

####Package needed

library(tidyverse) #data handling
library(PMA) #sparse CCA ; Not in cluster
library(caret) #split data methods
library(readxl) #open excel files: Annotation
library(patchwork) #put plots together
library(MetBrewer) #colors for plots
library(ggforce) #extension of ggplot, mainly for geom_sina
library(qgraph) #generate and plot correlation network ; Not in cluster
library(pheatmap) #heatmap plots
library(glmnet) #lasso
#library(xgboost) #if prediction is to be done using boosting trees, this package needs to be loaded
library(ppcor) #partial correlations to reproduce halla results


Location_dir = "~/Documents/GitHub/CMV_Associations/"


#Commonly used functions
source(paste0(Location_dir, "scripts/Functions.R"))
source(paste0(Location_dir, "scripts/Functions_multivariate.R"))

#Make sure commonly used functions are from dplyr
select = dplyr::select
rename = dplyr::rename

#Color palette
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
cp <- c(`Animal` = "dodgerblue2", `Bacteria` = "#E31A1C", `fungus` = "green4", `Phage`="#6A3D9A", `Plant`= "#FF7F00", `Protozoan`="black",`Virus`= "gold1")

#Reproducibility seed
set.seed(999)

########################
###1. READ DATA#########
########################

#Data needed:
#   - Immuno_matrix_postSelection: ~2,800 prevalent antibody responses. This is a binary table. Includes different cohorts. So it is subset to use only LLD 
#   - Quantifiaction_PhIPseq_TPM_twist_raw.tsv / Quantifiaction_PhIPseq_TPM_agilent_raw.tsv ; these are continuous tables of TPM computed with salmon
#   - Merged_Phenotypes_LLD.csv: LLD phenotypes. We add estrogen information and CMV_prediciton (from CMV script).
#   - CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt: Raw OLINK CVDIII
#   - Combined_data_uncorrected.tsv: Aging biomarkers from Andreu-Sanchez Comm.Biology 2022
#   - NMR.rds: NMR lipidomics measurements. MiMIR package is needed to use the tables of coefficients to compute metabolic bioage
#   - SupplementaryTable1.xlsx: This is the supp. table from Andreu-Sanchez Immunity 2023. Will be used to provide annotation to peptides
#   - BIOS_cell_types_DeconCell_2019-03-08.LLD_subset.txt: Decon deconvolution of cell types basde on bulk RNA-seq



################
################

#Antibodies PhIP-Seq
readRDS(paste0(Location_dir,"Data/Immuno_matrix_postSelection.rds")) -> ImmunoMatrix
Prepare_AB = function(ImmunoMatrix, N = "32"){
  'Keep samples from a specific cohort: 32=Lifelines Deep
  Merge PhIPSeq profiels with covariates. Change IDs to LLD IDs and remove NAs (mainly follow-up samples)
  Remove antibodies with missing data'
  if (! N == "*"){
    N = paste0(N, "_")
    ImmunoMatrix %>% filter(grepl(N, ID)) -> ImmunoMatrix
    ImmunoMatrix$ID = str_replace(ImmunoMatrix$ID , N, "")
  }
  read_tsv("~/Resilio Sync/Antibodies_WIS (1)/Plate records/Covariates.tsv") -> Phip_cov
  ImmunoMatrix %>% rename(sample_id = ID) %>% left_join(. , select(Phip_cov, c(ID, sample_id))) %>% select(-sample_id) %>% filter(! is.na(ID)) -> ImmunoMatrix
  ImmunoMatrix[ , colSums(is.na(ImmunoMatrix)) == 0] -> ImmunoMatrix #Removed antibodies with at least 1 NA
  return(ImmunoMatrix)
}
Prepare_AB(ImmunoMatrix, "33") -> ImmunoMatrix_IBD
Prepare_AB(ImmunoMatrix, "32") -> ImmunoMatrix

#Split matrix per chemestry: For exploration of differences between libraries. Not really important.
colnames(ImmunoMatrix)[grepl("agilent", colnames(ImmunoMatrix))] -> Agilent_peptides
colnames(ImmunoMatrix)[grepl("twist", colnames(ImmunoMatrix))] -> Twist_peptides
ImmunoMatrix %>% select(! Agilent_peptides) -> ImmunoMatrix_twist
ImmunoMatrix %>% select(! Twist_peptides) -> ImmunoMatrix_agilent

#For continuous data. TPMs of PhIP-Seq obtained using SAlmon on pair-end data. Sequences to map to were provided by WIS
Format_IP = F #This take quite some time to execute, thus, I saved the data as an rds after processing, and can be directly read without executing this code.
if ( Format_IP == T){
  'Open TPM counts and format them. I tested different log-ratio methods. End up suing median of ratios, but results are similar'
  read_tsv(paste0(Location_dir,"Data/Quantifiaction_PhIPseq_TPM_twist_raw.tsv")) %>% 
    as.data.frame() %>% column_to_rownames("Name") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% as_tibble() %>% write_rds(., "~/Documents/GitHub/CMV_Associations/Data/Quantifiaction_PhIPseq_TPM_twist_raw.rds")
  read_tsv(paste0(Location_dir,"Data/Quantifiaction_PhIPseq_TPM_agilent_raw.tsv")) %>% 
    as.data.frame() %>% column_to_rownames("Name") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% as_tibble() %>% write_rds(., "~/Documents/GitHub/CMV_Associations/Data/Quantifiaction_PhIPseq_TPM_agilent_raw.rds")
  
  #TPMs are accounting for differences in peptide length (agilent vs twist)
  
  #ALR transformation
  Agilent_tpm %>% select("agilent_221096", "agilent_133222") %>% apply(1, Geom_mean) -> Denom1
  Twist_tpm %>% select("twist_54788", "twist_54435","twist_47588" ) %>% apply(1, Geom_mean) -> Denom2
  #Median of ratios normalization (Ala DEseq2)
  Agilent_tpm %>% select("agilent_221096", "agilent_133222") %>% apply(2, Geom_mean) -> pseudoreference
  Agilent_tpm[1,] %>% select("agilent_221096", "agilent_133222") %>% apply(1, function(x){ median(x/pseudoreference) } ) -> Denom3
  #Transform
  Agilent_tpm %>% select(colnames(ImmunoMatrix_agilent)) %>% select(-ID) %>% apply(2, function(X){ log10((X+1)/Denom3) } ) %>% as_tibble() %>% mutate(ID=Agilent_tpm$ID, .before=1) -> Agilent_alr
  Twist_tpm %>% select(colnames(ImmunoMatrix_twist)) %>% select(-ID) %>% apply(2, function(X){ log10((X+1)/Denom3) } ) %>% as_tibble() %>% mutate(ID=Twist_tpm$ID, .before=1) -> Twist_alr
  #Merge
  Agilent_alr$ID %>% sapply(function(x){ str_split(x, "_")[[1]][4] } ) -> Agilent_alr$ID
  Twist_alr$ID %>% sapply(function(x){ str_split(x, "_")[[1]][4] } ) -> Twist_alr$ID
  
  left_join(Agilent_alr, Twist_alr) -> All_alr
  All_alr %>% filter(! ID %in% c("3260408512", "SF00115588", "SF00374694", "119410350", "218-32780" ) ) -> All_alr
  Prepare_AB(All_alr, N = "*") -> All_alr
  All_alr %>% write_rds(., paste0(Location_dir,"Data/Quantifiaction_PhIPseq_alr.rds"))
}
All_alr = read_rds( paste0(Location_dir,"Data/Quantifiaction_PhIPseq_alr.rds"))
All_alr %>% filter(grepl("LLD", ID)) %>% filter(! grepl("_F", ID) ) -> All_alr


######################
#Covariates LLD#######
######################
Covariates <- read_delim(paste0(Location_dir,"Data/Merged_Phenotypes_LLD.csv"), delim="|") #this is a comprenhensive list of phenotypes. We will only use a few to control for covariates
#Add estrogens medication to phenotypes
read_tsv(paste0(Location_dir,"Data/lld_estrogens.tsv")) -> Estrogens #Add strogens, which have been shown to affect Olink CMV panel
left_join(Covariates, Estrogens) -> Covariates
#Add prediciton of CMV. Check CMV prediciton script. This will be used to control for CMV effect, since many other peptide ab response tend to co-occur
read_tsv(paste0(Location_dir,"/Results/LLD_CMV_profile.tsv")) %>% filter(TimePoint == "Baseline") %>% select(ID, CMV_prediction) %>% left_join(Covariates, . ) -> Covariates


##################
##Other OMICS#####
##################


#Olink proteomics: CVD3 panel
read_tsv(paste0(Location_dir,"Data/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt")) -> Proteomics
Proteomics %>% t() %>% as.data.frame() %>%  rownames_to_column("ID") %>% as_tibble() %>% filter(!ID == "Protein") %>% `colnames<-`(c("ID", Proteomics$Protein))  -> Proteomics
Proteomics %>% select(-CCL22) -> Proteomics #Unfortunately, according to Olink this protein was removed from the panel due to technical issues. See: https://olink.com/bdnf-info/


#Biological aging biomarkers
Age_markers <- read_tsv(paste0(Location_dir,"Data/Combined_data_uncorrected.tsv"))
Age_markers %>% select(- colnames(Age_markers)[grepl("INT", colnames(Age_markers))] ) %>% rename(ID = Sample) -> Age_markers
Telomere_markers <- colnames(Age_markers)[grepl("MTL",colnames(Age_markers))]
Telomere_markers = Age_markers %>% dplyr::select(c("ID",Telomere_markers)) %>% mutate(MTL_gran = as.numeric(MTL_gran), `MTL_CD20+`=as.numeric(`MTL_CD20+`), `MTL_CD45+_CD20-`=as.numeric(`MTL_CD45+_CD20-`),`MTL_CD57+`=as.numeric(`MTL_CD57+`)) %>%
  `colnames<-`(c("ID", "TL_Lymphocytes", "TL_Granulocytes", "TL_Naive_T-cells", "TL_Memory_T-cells","TL_B-cells","TL_NK-cells")) #Mind that TL-NK cells contain also fully differentiated T-cells
Methylation_markers = Age_markers %>% select(ID, `Methyl.Age.Hannum.`, `Methyl.Age.Weidner.`, `Methyl.Age.Horvath.`) %>% drop_na()
sjTRECs = Age_markers %>% mutate(sjTREC = as.numeric(dCT))  %>% select(ID, sjTREC) %>% drop_na() -> sjTRECs

#Cytokines
Covariates %>% dplyr::select(ID,Citrullin, IL.1.Beta, IL.6, IL.8, IL.10, IL.12P70, `IL-18`, `IL-18bpx`, TNF.Alpha, `Leptin `, Resistin, `Adiponectin `, `AAT `, `hs-CRP`) %>% drop_na() %>% rename ( "IL_18" = `IL-18` , "AAT"=`AAT `, "Adiponectin"=`Adiponectin `, "Leptin"=`Leptin `, "hsCRP"=`hs-CRP`, "IL_18bpx" =`IL-18bpx`) -> Cytokines

#NMR metabolic age
read_rds(paste0(Location_dir,"Data/NMR.rds")) -> NMR
#Calculation of Metabolic age using MiMIR data
load("~/Documents/GitHub/MiMIR/data/metabo_names_translator.rda") #metabo_names_translator
load("~/Documents/GitHub/MiMIR/data/mort_betas.rda") #mort_betas #Deelen et al MetaboHealth
load("~/Documents/GitHub/MiMIR/data/CVD_score_betas.rda")  #"CVD_score_betas"
load("~/Documents/GitHub/MiMIR/data/PARAM_metaboAge.rda")

#NA are 0s. Pseudocunt of 1. Log10 and standarization (mean = 0, sd=1)
MetaboHealth = NMR %>% select(-ID) %>% apply(2, function(x){ x[is.na(x)] = 0 ; x = x+ 1 ; log10(x) -> y ; return(scale(y)) }) %>% as_tibble() %>% 
  mutate(ID = NMR$ID, MetaboHealth = XXL_VLDL_L*log(0.8) + S_HDL_L*log(0.87) + VLDL_D*log(0.85) + PUFA_FA_ratio*log(0.78) +
           Glc*log(1.16) + Lac*log(1.06) + His*log(0.93) + Ile*log(1.23) + Leu*log(0.82) + Gp*log(1.32) +
           Val*log(0.87) + Phe*log(1.13) +  AcAce*log(1.08) + Alb*log(0.89) ) %>% select(ID, MetaboHealth)
MetaboAge  = NMR %>% `colnames<-`( tolower(colnames(NMR)) ) %>% rename(faw3_fa=faw3_fa_ratio, faw6_fa=  faw6_fa_ratio, mufa_fa=mufa_fa_ratio, pufa_fa=pufa_fa_ratio, sfa_fa=sfa_fa_ratio) %>% select(-id)  %>% select(PARAM_metaboAge$MET) %>%
  apply(1, function(x){ x[is.na(x)] = 0 ; return( (x-PARAM_metaboAge$MEAN)/PARAM_metaboAge$SD ) }) %>% as_tibble()  %>%
  apply(2, function(x){ c(1,as_vector(x)) %*% PARAM_metaboAge$FIT_COEF } ) %>% as_tibble() %>% mutate(ID = NMR$ID, .before=1) %>% rename(MetaboAge= value)

#Merge all biological aging markers
Aging_health_markers = full_join(full_join(full_join(full_join(Methylation_markers, MetaboAge), MetaboHealth), sjTRECs), Telomere_markers)
write_tsv(Aging_health_markers, paste0(Location_dir,"Data/Aging_and_health_markers.tsv"))


##############################################
#Get information about phip-seq databases#####
##############################################

readxl::read_excel("~/Resilio Sync/Antibodies_WIS (1)/Results/Supplemenatry_Tables/SupplementaryTable1.xlsx", sheet="SupTable1.1") %>% dplyr::rename(Feature = Peptide)  -> Annotation

##This code is just necessary if enrichment analysis is performed from the output of sCCA

#Source database from peptide
DB_phipseq = list() #list of annotations
DB_phipseq_DF = tibble() #data frame of annotations
N = 1
#For enrichment of database analysis
for (i in unique(Annotation$Source_database)){
  Annotation %>% filter(Source_database == i) -> DB
  append(DB_phipseq, list(DB$Feature)) -> DB_phipseq
  names(DB_phipseq)[[N]] = i
  DB_phipseq_DF = rbind(DB_phipseq_DF, tibble(DataBase = i, Peptide=DB$Feature) )
  N = N + 1
}
#Taxonomic origin of peptide
Taxonomy_phipseq = list()
Taxonomy_phipseq_DF = tibble()
N=1
#For taxonomy enrichment
for (i in unique(Annotation$High_taxonomy)){
  Annotation %>% filter(High_taxonomy == i) -> DB
  append(Taxonomy_phipseq, list(DB$Feature)) -> Taxonomy_phipseq
  names(Taxonomy_phipseq)[[N]] = i
  Taxonomy_phipseq_DF = rbind(Taxonomy_phipseq_DF, tibble(Taxonomy = i, Peptide=DB$Feature) )
  
  N = N + 1
}


#######################
#Prepare Cell counts###
#######################

#Cell predictions
all_cell_types = read_tsv(paste0(Location_dir,"Data/BIOS_cell_types_DeconCell_2019-03-08.LLD_subset.txt")) -> all_cell_types
#format
all_cell_types %>% gather(key = var_name, value = value, 2:ncol(all_cell_types)) %>% spread_(key = names(all_cell_types)[1],value = 'value') %>% rename(ID = var_name) -> all_cell_types
# negative values to 0, recalculate to 100 and ALR transformation using monocytes as denominator
Clean_and_ALR(all_cell_types) -> cell_matrix_alr
cell_matrix_xlr = cell_matrix_alr #We choose to use ALR as log-ratio transformation

# logs of 0 give infinite, make NA and remove samples
cell_matrix_xlr[sapply(cell_matrix_xlr, is.infinite)] <- NA
cell_matrix_xlr %>% drop_na() -> cell_matrix_xlr
#Make a PCA to control for cell composition
res.pca <- prcomp(dplyr::select(cell_matrix_xlr,-ID), scale = TRUE)
Variability = 100* (res.pca$sdev/sum(res.pca$sdev))
PCs = as_tibble(res.pca$x) %>% mutate(ID = cell_matrix_xlr$ID)
PCs %>% dplyr::select(c(ID, PC1, PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)) -> Top_PCs
#The PCs will be used if there is the need to reduce dimensionality of the whole cell type composition
PCs_formula = "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10" 
Covariates %>% left_join(., Top_PCs) -> Covariates


###############################
###2. MATCH AND SPLIT##########
###############################
IDs_SC = read_tsv("Downloads/donor_md.txt") %>% rename(ID = donor)
Do_upset = function(All_IDs){
  tibble(ID = All_IDs) %>% mutate(PhIPSeq = ifelse(ID %in% ImmunoMatrix$ID, 1, 0), 
                                  Telomeres = ifelse(ID %in% Telomere_markers$ID, 1, 0),
                                  #Proteomics = ifelse(ID %in% Proteomics$ID, 1, 0),
                                  Cytokines = ifelse(ID %in% Cytokines$ID, 1, 0),
                                  Methylation_clocks = ifelse(ID %in% Methylation_markers$ID, 1, 0),
                                  sjTRECs = ifelse(ID %in% sjTRECs$ID, 1, 0),
                                  bulkRNAseq = ifelse(ID %in% cell_matrix_xlr$ID, 1, 0 ),
                                  scRNAseq = ifelse(ID %in% IDs_SC$ID, 1, 0 ),
                                  Metabolic_clocks = ifelse(ID %in% left_join(MetaboAge, MetaboHealth)$ID, 1, 0) ) %>%
    as.data.frame() %>% column_to_rownames("ID") %>% UpSetR::upset(nsets = 10, nintersects = 15, order.by = "freq", text.scale=2 ) -> DataSet_overview
  return(DataSet_overview)
  
}

ImmunoMatrix %>% filter(! grepl("_F", ID)) -> ImmunoMatrix_BL #get only baseline samples 
unique(c(Cytokines$ID, Telomere_markers$ID, Proteomics$ID, ImmunoMatrix_BL$ID)) -> All_IDs
DataSet_overview = Do_upset(All_IDs)

print("Selection of Test data")
set.seed(9854)
#Split into train and validation. Validtation data will be used to replicate assocations in train
sample( ImmunoMatrix_BL$ID, size = 0.2*length(ImmunoMatrix_BL$ID) ) -> Test_data

print("Upset of test")
Do_upset(Test_data) -> Test_overview

print("Upset of train")
Do_upset(tibble(ID = All_IDs) %>% filter(! ID %in% Test_data ) %>% select(ID) %>% as_vector() ) -> Train_overview

DataSet_overview ; Test_overview ; Train_overview


##########################################################
###3. Identify covariates of interest per dataset#########
##########################################################
#Covariates
KN =  c("ID", "Age", "Sex", "smk_now", "Lymph", "Eryth", "Mono", "Neutro", "Thrombo", "Eosino", "Blood_baso", "estrogens")

#We will ran permanova to identify which covariates are significantly associated witht he inter-sample dissimilarity
#Permanova
Telomere_markers %>% Run_permanova_covariates(. , Covariates, KN, 500 ) %>% mutate(Datset = "Telomeres") -> permanova_cov
Cytokines %>% Run_permanova_covariates(. , Covariates, KN, 500 ) %>% mutate(Datset = "Cytokines") %>% rbind(permanova_cov , .) -> permanova_cov
Proteomics %>% Run_permanova_covariates(. , Covariates, KN, 500 ) %>% mutate(Datset = "Proteomics") %>% rbind(permanova_cov ,. ) -> permanova_cov
Metabolites_clr %>% Run_permanova_covariates(. , Covariates, KN, 500 ) %>% mutate(Datset = "Metabolomics") %>% rbind(permanova_cov ,. ) -> permanova_cov
Methylation_markers %>% Run_permanova_covariates(. , Covariates, KN, 500 ) %>% mutate(Datset = "Methylation_age") %>% rbind(permanova_cov ,. ) -> permanova_cov

#Plot PERMANOVA
Naming = c("Anthro_Age", "Anthro_Sex", "Life_smk_now", "Cell_Lymph", "Cell_Eryth", "Cell_Mono", "Cell_Neutro", "Cell_Thrombo", "Cell_Eosino", "Cell_Blood_baso", "Drugs_estrogens")

permanova_cov %>% mutate(Feature = rep(Naming, 5)) %>% ggplot(aes(x=R2, y=Feature, fill= `Pr(>F)`<0.05 )) + geom_bar(stat="identity")  + theme_bw() + facet_wrap(~Datset,scales = "free") +  scale_fill_manual(values = c("red", "blue")) -> Variability_per_layer
Variability_per_layer %>% ggsave(paste0(Location_dir,"Results/Figures/Sup_PERMANOVA_confounders.png", .))

permanova_cov %>% filter(`Pr(>F)`<0.05)  #Proteomics evrything , Cytokines Sex/Ery; telomeres: Age/Sex/Neutro
#Univariate testing for markers that are not multivariable
lm( paste0("sjTREC ~ ", paste(KN[2:length(KN)], collapse = "+")), left_join(sjTRECs, Covariates) ) %>% summary() #use KN_methylation
lm( paste0("MetaboHealth ~ ", paste(KN[2:length(KN)], collapse = "+")), left_join(MetaboHealth, Covariates) ) %>% summary() #use KN_methylation
lm( paste0("MetaboAge ~ ", paste(KN[2:length(KN)], collapse = "+")), left_join(MetaboAge, Covariates) ) %>% summary() 
#Define covariates per dataset
KN_cytokines = c("ID", "Age", "Sex",  "Lymph", "Eryth", "Mono", "Neutro", "Thrombo", "Eosino", "Blood_baso","estrogens")
KN_methylation = c("ID", "Age", "Sex",  "Lymph", "Mono", "Neutro")
KN_health = c("ID", "Sex",  "Lymph", "Mono", "Eryth")

################################
######Antibody exposure score####
#################################
ImmunoMatrix %>% select(-ID) %>% apply(1, function(x){ x= x[!is.na(x)] ; sum(x) }  ) -> N_antibodies
left_join( Covariates, ImmunoMatrix %>% mutate(Antibody_N = N_antibodies) %>% select(ID, Antibody_N) ) -> Covariates

Run_associations_exposureScore = function(Covariates, Antibody_column ){
  cor.test( as_vector(Covariates[Antibody_column]), Covariates$Age, method="spearman") %>% print()
  left_join(Aging_health_markers, Covariates, by="ID") -> to_test
  results_score = tibble()
  for (marker in colnames(Aging_health_markers)  ){
    if (marker == "ID"){ next }
    to_test %>% select("Age", "Sex", marker, Antibody_column ) %>% drop_na() -> for_model
    pcor.test(x = as.vector(for_model[marker]) , y = as.vector(for_model[Antibody_column]) , z =for_model %>% select(c("Age","Sex")) , method = "spearman") -> Correlation_t
    Correlation_t %>% as_tibble() %>% mutate(bio_age = marker) %>% rbind(results_score, . ) ->  results_score
  }
  print(results_score %>% arrange(p.value ) %>%  mutate(FDR= p.adjust(p.value, "fdr")) )
}

Run_associations_exposureScore(Covariates,"Antibody_N")

#This can be stratified in different pathogen groups
#Pathogens - any
read_excel("~/Downloads/pathogens_final.xlsx") -> Pathogens
Annotation %>% filter(Taxa %in% filter(Pathogens, GROUP=="T")$TAXA) ->Pathogen_peptides
ImmunoMatrix %>% select(one_of(Pathogen_peptides$Feature)) %>% apply(1, function(x){ x= x[!is.na(x)] ; sum(x) }  ) -> N_antibodies_patho
left_join( Covariates, ImmunoMatrix %>% mutate(AntibodyPathogen_N = N_antibodies_patho) %>% select(ID, AntibodyPathogen_N) ) -> Covariates

Run_associations_exposureScore(Covariates,"AntibodyPathogen_N")

read_excel("~/Resilio Sync/Antibodies_WIS (1)/Results/Supplemenatry_Tables/SuplemmentaryTable2.xlsx", sheet = 8) -> DF
DF %>% filter(pheno == "Age") %>% select(peptide, Taxa, FDR, beta) %>% filter(beta > 0) -> Pos 
DF %>% filter(pheno == "Age") %>% select(peptide, Taxa, FDR, beta) %>% filter(beta < 0) -> Neg
ImmunoMatrix %>% select(Pos$peptide) %>% apply(1, function(x){ x= x[!is.na(x)] ; sum(x) }  ) -> N_antibodies_aging
ImmunoMatrix %>% select(Neg$peptide) %>% apply(1, function(x){ x= x[!is.na(x)] ; sum(x) }  ) -> N_antibodies_youth
left_join( Covariates, ImmunoMatrix %>% mutate(Young_ab = N_antibodies_youth, Old_ab = N_antibodies_aging ) %>% select(ID, Young_ab, Old_ab) ) -> Covariates

Run_associations_exposureScore(Covariates,"Young_ab")
Covariates %>% left_join(Aging_health_markers, . ) %>% ggplot(aes(y=MetaboHealth, x=Young_ab, col=as.factor(Sex)) ) + geom_point() + geom_smooth(method = "lm") + theme_bw()

Run_associations_exposureScore(Covariates,"Old_ab")




########################
###3. RUN HALLA#########
########################
###Prepare data for halla
###Residualizing omics layers####

#Proteomics
Remove_variability( filter(Proteomics, ! ID %in% Test_data ), Covariates, Covariates_n = KN[2:length(KN)] ) -> Proteomics_residuals
#Cytokines
Remove_variability( filter(Cytokines, ! ID %in% Test_data ), Covariates, Covariates_n = KN_cytokines[2:length(KN_cytokines)] ) -> Cytokines_residuals
#Telomeres
Remove_variability( filter(Telomere_markers,! ID %in% Test_data ), Covariates, Covariates_n =c("Age","Sex") ) -> Telomeres_residuals
#Methyl clock
Remove_variability( filter(Methylation_markers,! ID %in% Test_data ), Covariates, Covariates_n = KN_methylation[2:length(KN_methylation)] ) -> Methylation_residuals
#Biological aging. TL + methyl age, etc
Remove_variability( filter(Aging_health_markers,! ID %in% Test_data ), Covariates, Covariates_n =c("Age","Sex") ) -> Aging_health_residuals


write_tsv(For_halla(Proteomics_residuals), paste0(Location_dir,"Data/Residuals/Olink_residuals.tsv"))
write_tsv(For_halla(Cytokines_residuals), paste0(Location_dir,"Data/Residuals/Cytokines_residuals.tsv"))
write_tsv(For_halla(Telomeres_residuals), paste0(Location_dir,"Data/Residuals/Telomeres_residuals.tsv"))
write_tsv(For_halla(Methylation_residuals), paste0(Location_dir,"Data/Residuals/Methylation_residuals.tsv"))
write_tsv(For_halla(Aging_health_residuals), paste0(Location_dir,"Data/Residuals/Aging_residuals.tsv"))


#Filter based on prevalence?
write_tsv(For_halla(ImmunoMatrix %>% filter(!ID %in% Test_data)), paste0(Location_dir,"Data/Residuals/ImmunoMatrix_train.tsv"))
write_tsv(For_halla(ImmunoMatrix_agilent %>% filter(!ID %in% Test_data)), paste0(Location_dir,"Data/Residuals/ImmunoMatrix_train_agilent.tsv"))
write_tsv(For_halla(ImmunoMatrix_twist %>% filter(!ID %in% Test_data)), paste0(Location_dir,"Data/Residuals/ImmunoMatrix_train_twist.tsv"))
#Antibodies continuous
Remove_variability( filter(All_alr,! ID %in% Test_data ), Covariates, Covariates_n = c("Age", "Sex") ) -> Antibodies_residuals
write_tsv(For_halla(Antibodies_residuals %>% filter(!ID %in% Test_data)), paste0(Location_dir,"Data/Residuals/ImmunoMatrix_train_continuous.tsv"))


###Run halla on the shell
halla_To_run_in_shell = function(){
  Command = "  
  conda activate Microbiome
  #Important: I changed the `preprocess` function so that it should not discretize 1/0 traits (they are already discret). File: /Users/sergio/miniconda3/envs/Microbiome/lib/python3.10/site-packages/halla/utils/data.py
  #According to the simulations in the paper, increases in FNR between 0.1-0.4 dont result in important increasing of FDR. True FDR is always below the FDR threshold set.
  
  #Run with a less stringent FDR?
  #Stimulation
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Stimulation_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Stimulation_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.3 --rank_cluster best
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Stimulation_residuals_gonl.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Stimulation_train_gonl --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.3 --rank_cluster best
  #All aging
   halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Aging_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Aging_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.3 --rank_cluster best

  #Cytokines
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Cytokines_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Cytokines_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.3 --rank_cluster best
  #Telomeres
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Telomeres_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Telomeres_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.2 --rank_cluster best
  #Proteomics
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Olink_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Proteomics_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.2 --rank_cluster best
  #Metabolomics
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Metabolites_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Metabolites_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.2 --rank_cluster best
  #Methylation clock
  halla -x ~/Documents/GitHub/CMV_Associations/Data/Residuals/ImmunoMatrix_train.tsv -y ~/Documents/GitHub/CMV_Associations/Data/Residuals/Methylation_residuals.tsv --out_dir ~/Documents/GitHub/CMV_Associations/Results/HALLA/Clock_train --seed 95 --num_threads 4 -m spearman  --diagnostic_plot --fdr_alpha 0.05  --fnr_thresh 0.2 --rank_cluster best
  "
}


################################
###4. REPRODUCE RESULTS#########
################################
read_tsv(paste0(Location_dir,"Results/HALLA/Aging_train/sig_clusters.txt")) -> Halla_Aging
read_tsv(paste0(Location_dir,"Results/HALLA/Cytokines_train/sig_clusters.txt")) -> Halla_cyt #binary data
read_tsv(paste0(Location_dir,"Results/HALLA/Cytokines_train_cont/sig_clusters.txt")) -> Halla_cyt2 #Cont. data

#Aging biomarkers
#1. validation controlling for age and sex
Check_cluster(C = Halla_Aging, Original_X = ImmunoMatrix, Original_Y = Aging_health_markers, Make_Plot = F,Confounder_names = c("Age", "Sex") ) -> Aging_halla_check2
#2. validation controlling also for CMV
Check_cluster(C = Halla_Aging, Original_X = ImmunoMatrix, Original_Y = Aging_health_markers, Make_Plot = F, Confounder_names = c("Age", "Sex", "CMV_prediction") ) -> Aging_halla_check3
rbind(Aging_halla_check2 %>% mutate(Controlled="Age+Sex"), Aging_halla_check3 %>% mutate(Controlled="Age+Sex+CMV")) %>% write_tsv(. , paste0(Location_dir,"Results/halla_aging_test.tsv") )

#3. Generate HEATMAP
#All trained-singificant hits
str_split(Halla_Aging$cluster_X, ";") %>% unlist() -> all_check2
#Read all associations from HAllA. Keep only the Peptides that had at least one association in train
read_tsv(paste0(Location_dir,"Results/HALLA/Aging_train/all_associations.txt")) %>% filter(X_features %in% all_check2) -> For_aging_plot
#Make wide format: Telomeres as columns, each value a correlation with Features (row)
For_aging_plot %>% select( X_features, Y_features, association) %>% spread(Y_features, association) -> For_aging_plot_wid
#For annotation, rename 
For_aging_plot_wid %>% left_join(Annotation %>% rename(X_features = Feature ) %>% select(X_features, Taxa, High_taxonomy, Source_database)  ) %>% mutate(Taxa = ifelse(Taxa %in% c("Enterovirus", "Human rhinovirus", "root") , "Rhinovirus", Taxa)  ) %>%
  mutate(Taxa = ifelse(Taxa %in% c("Human betaherpesvirus 5 (CMV)") , "CMV", Taxa)  ) %>%
  mutate(Taxa = ifelse(Taxa %in% c("Human gammaherpesvirus 4 (EBV)") , "EBV", Taxa)  ) %>%
  mutate(Taxa = ifelse(High_taxonomy == "Phage" , "Phage", Taxa)  ) %>%
  mutate(Taxa = ifelse(Source_database == "is_allergens" , "Allergen", Taxa)  ) %>%
  mutate(Taxa = ifelse(Source_database == "is_gut_microbiome" , "microbiota", Taxa)  ) %>% 
  mutate(Taxa = ifelse(Source_database == "is_VFDB" , "Virulence factor", Taxa)) -> For_annotation
For_annotation %>% select(X_features, Taxa) %>% as.data.frame() %>% column_to_rownames("X_features") -> For_annotation2
#Generate a matrix of * where the annotation should go in the heatmap
Annotation_wide_aging = tibble()
for (peptide_row in unique(For_aging_plot_wid$X_features)){
  #Go per row. Is this peptide significant in test?
  Aging_halla_check2 %>% filter(X == peptide_row ) %>% filter(P_test < 0.05, P_train<0.05) %>% filter(sign(rho_train) == sign(rho_test) ) -> To_add
  Aging_halla_check3 %>% filter(X == peptide_row ) %>% filter(P_test < 0.05, P_train<0.05) %>% filter(sign(rho_train) == sign(rho_test) ) -> To_add3
  
  Make_wide = tibble()
  for (TL in c("TL_NK-cells","TL_Lymphocytes","TL_Granulocytes", "TL_Naive_T-cells","TL_Memory_T-cells","MetaboAge"  ) ){
    ifelse(dim(filter(To_add3, Y==TL))[1] == 0, "" , "*") -> R1
    #ifelse(dim(filter(To_add2, Y==TL))[1] == 0, "" , "+") ->R2
    ifelse(dim(filter(To_add, Y==TL))[1] == 0, "" , "*") -> R3
    Anno = paste0(R1, R3)
    tibble(X_features=peptide_row, Telomere=TL, A=Anno ) %>% rbind(Make_wide, . ) -> Make_wide
  }
  Make_wide %>% spread(Telomere,A) %>% rbind(Annotation_wide_aging, .) ->Annotation_wide_aging 
  Annotation_wide_aging %>% select(X_features, `TL_NK-cells`, TL_Lymphocytes, TL_Granulocytes, `TL_Naive_T-cells`,  `TL_Memory_T-cells`, MetaboAge)->Annotation_wide_aging 
  
}
Annotation_wide_aging %>% as.data.frame() %>% column_to_rownames("X_features") -> Annotation_wide_aging
#Make heatmap and save
pdf(paste0(Location_dir,"Results/Figures/Heatmap_halla_Aging.pdf"))
For_aging_plot_wid %>% select(X_features, `TL_NK-cells`, TL_Lymphocytes, TL_Granulocytes, `TL_Naive_T-cells`,  `TL_Memory_T-cells`, MetaboAge) %>% as.data.frame() %>% column_to_rownames("X_features") %>%
  pheatmap::pheatmap(. , annotation_row = For_annotation2 , fontsize_row = 4, display_numbers = Annotation_wide_aging,
                     number_color = "black", fontsize_number=8 )
dev.off()                     


#Cytokines
#Validation. Using Inv rank trasnformed cytokine data. All covariates.
Check_cluster(Halla_cyt  , ImmunoMatrix, Cytokines%>% select(-ID) %>% apply(2, InvRank) %>% as_tibble() %>% mutate(ID = Cytokines$ID), Make_Plot=F) -> Cytokines_halla_check
Cytokines_halla_check %>% select(X, Y, rho_train) %>% spread(Y,rho_train)
write_tsv(Cytokines_halla_check , paste0(Location_dir,"Results/halla_cytokines_test.tsv") )


#For plots
#Get all associations, make them into a matrix
read_tsv(paste0(Location_dir,"/Results/HALLA/Cytokines_train/all_associations.txt")) %>% filter(X_features %in% Halla_cyt$cluster_X, Y_features %in% Halla_cyt$cluster_Y ) -> For_cytokine_plot
For_cytokine_plot %>% select( X_features, Y_features, association) %>% spread(Y_features, association) -> For_cytokine_plot_wid
For_cytokine_plot_wid %>% as.data.frame() %>% column_to_rownames("X_features") %>% pheatmap::pheatmap(. , annotation_row = Annotation %>% filter(Feature %in%  For_cytokine_plot_wid$X_features ) %>% select(Feature, Taxa) %>% as.data.frame() %>% column_to_rownames("Feature") )
#Adenovirus C strong association with IL_18bpx --> this is also reproduced in Component 1 of sCCA (top load)
#Firmicutes: Not supported associaton with TNF.alpha
#EBV nuclear antigen (x3) associations with IL_18bpx: Also reproduced in component 1 of sCCA (twist_34374 & agilent_3387 & agilent_7263)

#Cytokine analysis but using continuous PhIP-Seq data
Check_cluster(Halla_cyt2, All_alr, Cytokines%>% select(-ID) %>% apply(2, InvRank) %>% as_tibble() %>% mutate(ID = Cytokines$ID), Make_Plot = F ) -> Cytokines_halla_check2
read_tsv(paste0(Location_dir,"Results/HALLA/Cytokines_train_cont/all_associations.txt")) %>% #filter(X_features %in% Halla_cyt2$cluster_X, Y_features %in% Halla_cyt2$cluster_Y ) -> For_cytokine_plot
  select( X_features, Y_features, association) %>% filter(Y_features %in% c("TNF.alpha", "IL_18", "IL_18bpx", "Resistin") ) %>% filter(X_features %in% filter(Cytokines_halla_check2,P_test<0.05 )$X )  %>% spread(Y_features, association) -> For_cytokine_plot_wid
Annotation_wide= tibble()
for (peptide_row in unique(For_cytokine_plot_wid$X_features)){
  Cytokines_halla_check2 %>% filter(X == peptide_row ) %>% filter(P_test < 0.05) -> To_add
  tibble(X_features=peptide_row , IL_18=ifelse(dim(filter(To_add, Y=="IL_18"))[1] == 0, "" , "*"),
         IL_18bpx=ifelse(dim(filter(To_add, Y=="IL_18bpx"))[1] == 0, "" , "*"), Resistin=ifelse(dim(filter(To_add, Y=="Resistin"))[1] == 0, "" , "*")) %>% rbind(Annotation_wide, .) ->Annotation_wide 
}
Annotation_wide %>% as.data.frame() %>% column_to_rownames("X_features") -> Annotation_wide
For_cytokine_plot_wid %>% as.data.frame() %>% column_to_rownames("X_features") %>% pheatmap::pheatmap(. , annotation_row = Annotation %>% filter(Feature %in%  For_cytokine_plot_wid$X_features ) %>% select(Feature, High_taxonomy) %>% as.data.frame() %>% column_to_rownames("Feature"), fontsize_row = 6, display_numbers = Annotation_wide )



#######################################
###5. SPARSE CCA ######################
#######################################

#Subset peptides to perform sCCA on
ImmunoMatrix %>% select(-ID) %>% apply(2, function(x){ sum(x)/length(x) } ) %>% as.data.frame() %>% rownames_to_column("Peptide") %>% as_tibble() %>% filter(`.`<0.1) -> Keep_peptide


############
#Proteomics#
############
filter(Proteomics, ! ID %in% Test_data) %>% select(-ID) %>% apply(., 2, InvRank) %>% as_tibble() %>% mutate(ID=filter(Proteomics, ! ID %in% Test_data)$ID, .before=1) %>% Remove_variability(. , Covariates, Covariates_n = KN[2:length(KN) ] ) -> Proteomics_for_CCA
filter(Proteomics,  ID %in% Test_data) %>% select(-ID) %>% apply(., 2, InvRank) %>% as_tibble() %>% mutate(ID=filter(Proteomics,  ID %in% Test_data)$ID, .before=1) %>% Remove_variability(. , Covariates, Covariates_n = KN[2:length(KN) ] ) -> Proteomics_for_CCA_test

Results_proteomics = Analysis(Proteomics_for_CCA, ImmunoMatrix, test_omics = Proteomics_for_CCA_test)

Plot_loads(Results_proteomics, Component_n= 1,Component_x = "Component_Proteomics") -> PlotProteomicsCCA
R_ora_proteomics = Run_ora(Results_proteomics,pathways = Taxonomy_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 1,Data_type = "Antibody")  
R_ora_proteomics2 = Run_ora(Results_proteomics,pathways = DB_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 1,Data_type = "Antibody")  

write_rds(Results_proteomics, paste0(Location_dir,"Results/sCCA/Proteomics.rds"))


############
#Cytokines#
############
#Inv rank transform train and test independently
filter(Cytokines, ! ID %in% Test_data) %>% select(-ID) %>% apply(., 2, InvRank) %>% as_tibble() %>% mutate(ID=filter(Cytokines, ! ID %in% Test_data)$ID, .before=1) %>% Remove_variability(. , Covariates, Covariates_n = KN_cytokines[2:length(KN_cytokines) ] ) -> Cytokines_for_CCA
filter(Cytokines,  ID %in% Test_data) %>% select(-ID) %>% apply(., 2, InvRank) %>% as_tibble() %>% mutate(ID=filter(Cytokines,  ID %in% Test_data)$ID, .before=1) %>% Remove_variability(. , Covariates, Covariates_n = KN_cytokines[2:length(KN_cytokines) ] ) -> Cytokines_for_CCA_test
Results_cytokines = Analysis(Cytokines_for_CCA, ImmunoMatrix, test_omics = Cytokines_for_CCA_test)

R1_ora = Run_ora(Results_cytokines,pathways = DB_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 1,Data_type = "Antibody")
R2_ora = Run_ora(Results_cytokines,pathways = Taxonomy_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 1,Data_type = "Antibody")  
R2_ora = Run_ora(Results_cytokines,pathways = Taxonomy_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 4,Data_type = "Antibody")  
Results_cytokines[[1]] %>% filter(Component == 4) %>% filter(Data=="Antibody") %>% arrange(desc(abs(Coefficient))) %>% View()

write_rds(Results_cytokines, paste0(Location_dir,"Results/sCCA/Cytokines.rds"))

#Plots cytokine top loads
##Component 1
Results_cytokines[[1]] %>% filter(Data=="Antibody", Component == 1 ) %>% arrange(desc(abs(Coefficient))) %>% head(10) %>% mutate(Taxa = ifelse(Taxa == "root", "Human adenovirus C", Taxa) ) %>% ggplot(aes(x=Feature, y=Coefficient, col=Taxa )) + theme_bw() + coord_flip() + geom_segment(aes(x=Feature, xend=Feature, y=0, yend=Coefficient), color="black") +
  geom_point(size=4) + scale_color_manual(values = c(met.brewer(name="Klimt", n=4,type="discrete")) ) -> PlotLoads
ggsave(paste0(Location_dir,"Results/Figures/top_loads_sCCA_cytokines.pdf"),PlotLoads, height = 2)
##Component 4
Results_cytokines[[1]] %>% filter(Data=="Antibody", Component == 4 ) %>% arrange(desc(abs(Coefficient))) %>% head(10)  %>% ggplot(aes(x=Feature, y=Coefficient, col=Taxa )) + theme_bw() + coord_flip() + geom_segment(aes(x=reorder(Feature, Coefficient), xend=reorder(Feature, Coefficient), y=0, yend=Coefficient), color="black") +
  geom_point(size=4)  + scale_color_manual(values = c("#df9ed4",  "#469d76", "#c93f55", "#eacc62", "#3c4b99", "#924099") ) -> PlotLoads
ggsave(paste0(Location_dir,"Results/Figures/top_loads_sCCA_cytokines_4.pdf"),PlotLoads, height = 2)



############  
#Telomeres##
############
Telomere_markers %>% filter(! ID %in% Test_data) %>% Impute_telomere_data(.) -> List_telomere_imputation_train
Telomere_markers_train = List_telomere_imputation_train[[1]]
Remove_variability(Telomere_markers_train, Covariates, Covariates_n = c("Age", "Sex") ) -> Telomeres_corrected_train

Results_telomeres = Analysis( Telomeres_corrected_train, ImmunoMatrix,test_omics =  Remove_variability(filter(Telomere_markers, ID %in% Test_data) , Covariates, Covariates_n = c("Age", "Sex") ),  penaltyX = seq(0.05,0.4,length=10) )

#Check Results
Plot_loads(Results_telomeres, Component_n= 1,Component_x = "Component_Telomere") -> PlotTelomereCCA
Plot_loads(Results_telomeres, Component_n= 2,Component_x = "Component_Telomere") -> PlotTelomereCCA2

R1_ora = Run_ora(Results_telomeres,pathways = DB_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 1,Data_type = "Antibody")
R2_ora = Run_ora(Results_telomeres,pathways = Taxonomy_phipseq_DF , Universe = colnames(select(ImmunoMatrix,-ID)) , Component_n = 1,Data_type = "Antibody")  
#Enrichment of viruses both as positive and negative loads
write_rds(Results_telomeres, paste0(Location_dir,"Results/sCCA/Telomeres.rds"))



###############################
#####Prediction models#########
###############################
set.seed(9878)

#Some extra functions for ML
Predict_omic = function(Cyto, X_omic, Y_omic,test_ids,  Model = c("xgboost", "lasso"), score = "rho", Covariates_n=c("Age", "Sex")  ){
  'Train a model to predict a Feature using antibody-bound peptide profiles from PhIP-Seq. It accepts IDs that belong to a test set. Computes R2 on the predicted values'
  #Score can be rho or r2
  #Functions to calculate r2
  Pseudo_Rsqr = function(x1, x2){ R2 = 100*(1 - (mse(x1,x2)/var(x1))) ; if(R2<0){ print("Negative R2, prediction is worse than the mean") ; R2 = 0} ;  return(paste0(R2, "%"))   }
  mse = function(x1, x2){  mean((x1 - x2)^2)  }
  
  print("X_omic should contain all variables to test, so remove variables BEFORE calling the function")
  print(paste0("Please make sure that ", Cyto, " is in your Y_omic"))
  
  
  print("Merging data and splitting in train and test sets")
  Covariates %>% select(c("ID",Covariates_n)) %>% left_join(. , Y_omic, by="ID")  %>% left_join(. , X_omic, by="ID") %>% drop_na() -> For_model
  print(paste0("Complete data. Samples:",dim(For_model)[1], " Features:",dim(For_model)[2] ))
  #Define train and test
  For_model %>% filter(! ID %in%  test_ids) %>% select(-ID) -> For_model_train
  print(paste0("Train data. Samples:",dim(For_model_train)[1], " Features:",dim(For_model_train)[2] ))
  For_model %>% filter( ID %in%  test_ids) %>% select(-ID) -> For_model_test
  print(paste0("Test data. Samples:",dim(For_model_test)[1], " Features:",dim(For_model_test)[2] ))
  
  Results = list()
  Entry = 1
  
  Xs = c( colnames(X_omic %>% select(-ID)) , Covariates_n )
  Ys =  as_vector(select(For_model_train, Cyto))
  for ( model_name in Model){
    print( paste0("Running ", model_name) )
    if (model_name == "xgboost"){
      library(xgboost)
      bst <- xgboost(data = For_model_train %>% select(Xs) %>% as.matrix(), label = Ys, 
                     max_depth = 6, eta = 0.0025, nrounds = 4000, subsample = 0.6, min_child_weight = 20,
                     nthread = 2, objective = "reg:squarederror", verbose=0) #verbose=1 prints train-rmse per iteration
      print( paste0(model_name, " is finished. Obtaining scores") )
      predict(bst, newdata =  as.matrix(For_model_test %>% select( Xs ) ) ) -> y_hat
      
      if (score == "rho"){
        score_r = cor(as_vector(select(For_model_test, Cyto)), y_hat)
      } else if (score == "r2") {
        score_r = Pseudo_Rsqr( as_vector(select(For_model_test, Cyto)), y_hat )
      }
      Results[[Entry]] = score_r
      Entry = Entry+1
      
    } else if( model_name == "lasso"){
      library(glmnet)
      bst <- glmnet::cv.glmnet( y = Ys , x= For_model_train %>% select(Xs) %>% as.matrix(),
                                alpha =1 , nfolds = 10, type.measure="mse", standardize = F, standardize.response = T )#, lambda=c(0.1, 0.5 ,1, 10) )
      print( paste0(model_name, " is finished. Obtaining scores") )
      plot(bst)
      
      Cof = coef(bst, s = "lambda.min") %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("Feature") %>% as_tibble() %>% arrange(desc(abs(`1`)))
      predict(bst, newx =  select( For_model_test, Xs )  %>% as.matrix(), s= "lambda.min"   ) -> y_hat
      
      if (score == "rho"){
        score_r = cor(as_vector(select(For_model_test, Cyto)), y_hat )
      } else if (score == "r2") {
        score_r = Pseudo_Rsqr( as_vector(select(For_model_test, Cyto)), y_hat )
      }
      Results[[Entry]] = score_r
      Entry = Entry+1
      Results[[Entry]] = Cof
      Entry = Entry+1
      
    }
  }
  #Base model
  print("Running base model")
  bst3 = lm( as.formula(paste0("`", Cyto, "` ~ ", paste0(Covariates_n, collapse="+"))), For_model_train )
  print("Base model is done. Obtaining scores")
  
  predict(bst3, newdata = For_model_test ) -> y_hat_cov
  
  if (score == "rho"){
    Base_score = cor(For_model_test %>% select(Cyto) %>% as_vector() ,y_hat_cov)
  }else if(score == "r2"){
    Base_score = Pseudo_Rsqr( as_vector(select(For_model_test, Cyto)),  y_hat_cov ) 
  }
  Results[[Entry]] = Base_score
  return(Results)
}
Iterate_features = function(Omic_layer, folds, Covariates_names=c("Age", "Sex")){
  'For a given omic layer, iterates through features. Per feature, iterates per fold and applies Predict_omic to obtain an R2'
  
  Fold_results = tibble()
  Score_results = tibble()
  for (Cytokine in c(colnames(Omic_layer)) ){
    fold_n = 1
    if (Cytokine == "ID"){ next }
    print(paste0("Running ", Cytokine))
    for (fold in folds){
      print(paste0("Fold ", fold_n))
      #  Test_ids = Covariates  %>% filter( ID %in%  For_model_test$ID) %>% select(ID) %>% as_vector()
      results_prediction = Predict_omic(Cytokine, X_omic =  All_alr %>% select(colnames(ImmunoMatrix_agilent)) , Y_omic = Omic_layer, test_ids = Covariates$ID[fold],
                                        Model = c("lasso"), score= "rho", Covariates_n=Covariates_names )
      results_prediction2 = Predict_omic(Cytokine, X_omic =  ImmunoMatrix %>% select(colnames(ImmunoMatrix_agilent)) , Y_omic = Omic_layer, test_ids = Covariates$ID[fold],
                                         Model = c("lasso"), score= "rho", Covariates_n=Covariates_names )
      
      Fold_results = rbind(Fold_results, tibble(Fold = fold_n, Rho_lasso =  results_prediction[[1]][1], Rho_null = results_prediction[[3]], Mode="Continuous", Y= Cytokine) )
      Fold_results = rbind(Fold_results, tibble(Fold = fold_n, Rho_lasso =  results_prediction2[[1]][1], Rho_null = results_prediction2[[3]], Mode="Binary", Y=Cytokine ) )
      
      Score_results %>% rbind(., results_prediction[[2]] %>% mutate(Fold=fold_n, Mode="Continuous",Y=Cytokine )  ) %>% rbind(. , results_prediction2[[2]] %>% mutate(Fold=fold_n, Mode="Binary",Y=Cytokine )) -> Score_results
      
      fold_n = fold_n + 1
    }
  }  
  
  return(list(Fold_results, Score_results))
  
}
Only_one_feature = function(Omic_layer,folds,Covariates_names, Name_pheno, Prevalence_threshold ){
  'Applies predict_omic to get an R2 from a feature, using train/test split. This function does not iterate through all features. Use for validation with a different set of covariates as M0'
  Fold_results = tibble()
  Score_results = tibble()
  fold_n = 1
  #Xs to run
  ImmunoMatrix_agilent %>% select(-ID) %>% apply(2, function(x){x=x[!is.na(x)] ;  sum(x)/length(x) }  ) %>% as.data.frame() %>% rownames_to_column("Peptide") %>% as_tibble() %>%
    filter(`.` > Prevalence_threshold) %>% select(Peptide) %>% as_vector() -> To_run
  To_run = c("ID", To_run)
  for (fold in folds){
    print(paste0("Fold ", fold_n))
    #  Test_ids = Covariates  %>% filter( ID %in%  For_model_test$ID) %>% select(ID) %>% as_vector()
    results_prediction = Predict_omic(Name_pheno, X_omic =  All_alr %>% select(To_run) , Y_omic = Omic_layer, test_ids = Covariates$ID[fold],
                                      Model = c("lasso"), score= "rho", Covariates_n=Covariates_names )
    results_prediction2 = Predict_omic(Name_pheno, X_omic =  ImmunoMatrix %>% select(To_run) , Y_omic = Omic_layer, test_ids = Covariates$ID[fold],
                                       Model = c("lasso"), score= "rho", Covariates_n=Covariates_names )
    
    Fold_results = rbind(Fold_results, tibble(Fold = fold_n, Rho_lasso =  results_prediction[[1]][1], Rho_null = results_prediction[[3]], Mode="Continuous", Y= Name_pheno) )
    Fold_results = rbind(Fold_results, tibble(Fold = fold_n, Rho_lasso =  results_prediction2[[1]][1], Rho_null = results_prediction2[[3]], Mode="Binary", Y=Name_pheno ) )
    
    Score_results %>% rbind(., results_prediction[[2]] %>% mutate(Fold=fold_n, Mode="Continuous",Y=Name_pheno )  ) %>% rbind(. , results_prediction2[[2]] %>% mutate(Fold=fold_n, Mode="Binary",Y=Name_pheno )) -> Score_results
    
    fold_n = fold_n + 1
  }
  
  return(list(Fold_results, Score_results))
  
  
}


###Apply Iterate_features to get continuous and binary predictions of a given omci layer. Applied to Cytokines, Telomeres and OLINK
Repeat_cv_result_cyt = tibble()
Repeat_cv_result_tl = tibble()
Repeat_cv_result_proteomics = tibble()
for (CV_repeat in seq(10)){
  #We can repeat this step in a loop, so that we have X CV cases.
  print(paste0("CV repetition", CV_repeat))
  folds <- createFolds(Covariates$ID, k = 10)
  
  Fold_results = Iterate_features(Cytokines, folds)
  Fold_results_telomere =  Iterate_features(Telomere_markers, folds)
  Fold_results_proteomics =  Iterate_features(Proteomics %>% select(-ID) %>% apply(2, as.numeric) %>% as_tibble() %>% mutate(ID = Proteomics$ID, .before=1) , folds, Covaraites_names=KN[2:length(KN)] )

  Fold_results[[2]] -> Scores
  Fold_results[[1]] -> Fold_results
  
  
  Fold_results_telomere[[2]] -> Scores_telomere
  Fold_results_telomere[[1]] -> Fold_results_telomere
  
  
  Fold_results_proteomics[[2]] -> Scores_preomics
  Fold_results_proteomics[[1]] -> Fold_results_proteomics
  
 
  #Formatting
  Fold_results %>% filter(Mode=="Continuous") %>% select(Fold, Y, Rho_null) %>% rename(Null=Rho_null)  -> Null_model
  Fold_results %>% select(-Rho_null) %>% spread(Mode, Rho_lasso) %>% left_join(Null_model, by=c("Fold", "Y")) %>% gather(Model, Rho, c(Binary, Continuous, Null)) %>% mutate(Rho = ifelse(is.na(Rho), 0, Rho )) -> Fold_results
  
  Fold_results_telomere %>% filter(Mode=="Continuous") %>% select(Fold, Y, Rho_null) %>% rename(Null=Rho_null)  -> Null_model
  Fold_results_telomere %>% select(-Rho_null) %>% spread(Mode, Rho_lasso) %>% left_join(Null_model,  by=c("Fold", "Y")) %>% gather(Model, Rho, c(Binary, Continuous, Null)) %>% mutate(Rho = ifelse(is.na(Rho), 0, Rho )) -> Fold_results_telomere
  
  Fold_results_proteomics %>% filter(Mode=="Continuous") %>% select(Fold, Y, Rho_null) %>% rename(Null=Rho_null)  -> Null_model
  Fold_results_proteomics %>% select(-Rho_null) %>% spread(Mode, Rho_lasso) %>% left_join(Null_model,  by=c("Fold", "Y")) %>% gather(Model, Rho, c(Binary, Continuous, Null)) %>% mutate(Rho = ifelse(is.na(Rho), 0, Rho )) -> Fold_results_proteomics
  
  
  Repeat_cv_result_cyt %>% rbind(. , Fold_results %>% mutate(Iteration = CV_repeat) ) -> Repeat_cv_result_cyt
  Repeat_cv_result_tl %>% rbind(. , Fold_results_telomere %>% mutate(Iteration = CV_repeat) ) -> Repeat_cv_result_tl
  Repeat_cv_result_proteomics %>% rbind(. , Fold_results_proteomics %>% mutate(Iteration = CV_repeat) ) -> Repeat_cv_result_proteomics

}  
#Visualize results and save
write_tsv(Repeat_cv_result_cyt, paste0(Location_dir,"~Results/Lasso_predictions_cytokines.tsv"))
Repeat_cv_result_cyt %>% ggplot(aes(x=Rho, fill=Model)) + 
  geom_density(alpha=0.4) + facet_wrap(~Y, scales = "free") +  theme_bw() + geom_vline(xintercept = 0)
write_tsv(Repeat_cv_result_tl, paste0(Location_dir,"Results/Lasso_predictions_telomeres.tsv"))
Repeat_cv_result_tl %>% ggplot(aes(x=Rho, fill=Model)) + 
  geom_density(alpha=0.4) + facet_wrap(~Y) +  theme_bw() + geom_vline(xintercept = 0)
write_tsv(Repeat_cv_result_proteomics, paste0(Location_dir,"Results/Lasso_predictions_proteomics.tsv"))
#Boxplot of results
Repeat_cv_result_cyt %>% filter(Y  ==  "IL_18bpx")  %>% ggplot(aes(x=Model, y= Rho, fill=Model)) + geom_boxplot() + ggforce::geom_sina()  +  theme_bw() + coord_flip() + geom_hline(yintercept = 0) + scale_fill_manual(values = c(met.brewer(name="Klimt", n=3,type="discrete")) ) + ylab("Pearson's correlation\n Predicted IL.18bpx vs Measured IL.18bpx")
Repeat_cv_result_cyt %>% filter(Y  ==  "Resistin")  %>% ggplot(aes(x=Model, y= Rho, fill=Model)) + geom_boxplot() + ggforce::geom_sina()  +  theme_bw() + coord_flip() + geom_hline(yintercept = 0) + scale_fill_manual(values = c(met.brewer(name="Klimt", n=3,type="discrete")) ) + ylab("Pearson's correlation\n Predicted Resistin vs Resistin")
Repeat_cv_result_tl %>% filter(Y  ==  "TL_NK-cells")  %>% ggplot(aes(x=Model, y= Rho, fill=Model)) + geom_boxplot() + ggforce::geom_sina()  +  theme_bw() + coord_flip() + geom_hline(yintercept = 0) + scale_fill_manual(values = c(met.brewer(name="Klimt", n=3,type="discrete")) ) + ylab("Pearson's correlation\n Predicted NK-cell TL vs Measured NK-cell TL")

##Statistical Comparison between models
#Cyto
Comparison_ml = tibble()
for (Cyt in colnames(Cytokines)){
  if (Cyt == "ID"){ next }
  Repeat_cv_result_cyt %>% mutate(Model = factor(Model, levels = c("Null", "Continuous","Binary") )) %>% filter(Y == Cyt) %>%  lm(Rho ~ Model, .) %>% summary() -> R
  as.data.frame(R$coefficients)  %>% rownames_to_column("Term") %>% as_tibble() %>% filter(! Term == "(Intercept)" ) %>% mutate(Feature = Cyt) %>% rbind(Comparison_ml, . ) -> Comparison_ml
}
Comparison_ml %>% filter(Estimate > 0)
#TL
for (TL in colnames(Telomere_markers)){
  if (TL == "ID"){ next }
  Repeat_cv_result_tl %>% mutate(Model = factor(Model, levels = c("Null", "Continuous", "Binary") )) %>% filter(Y == TL) %>%  lm(Rho ~ Model, .) %>% summary() -> R
  as.data.frame(R$coefficients)  %>% rownames_to_column("Term") %>% as_tibble() %>% filter(! Term == "(Intercept)" ) %>% mutate(Feature = TL) %>% rbind(Comparison_ml, . ) -> Comparison_ml
}
Comparison_ml %>% filter(Estimate > 0)
#OLINK
for (Prot in unique(Repeat_cv_result_proteomics$Y)){
  if (Prot == "ID"){ next }
  Repeat_cv_result_proteomics %>% mutate(Model = factor(Model, levels = c("Null", "Continuous", "Binary") )) %>% filter(Y == Prot) %>%  lm(Rho ~ Model, .) %>% summary() -> R
  as.data.frame(R$coefficients)  %>% rownames_to_column("Term") %>% as_tibble() %>% filter(! Term == "(Intercept)" ) %>% mutate(Feature = Prot) %>% rbind(Comparison_ml, . ) -> Comparison_ml
}
Comparison_ml %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) %>% filter(Estimate > 0) %>% arrange(`Pr(>|t|)`)
Repeat_cv_result_proteomics %>% filter(Y %in% filter(Comparison_ml, Estimate>0)$Feature  ) %>%  ggplot(aes(x=Model, y= Rho, fill=Model) ) + 
  geom_boxplot() + ggforce::geom_sina(alpha=0.2)  +  theme_bw() + coord_flip() + geom_hline(yintercept = 0) + scale_fill_manual(values = c(met.brewer(name="Klimt", n=3,type="discrete")) ) + facet_wrap(~Y)
Repeat_cv_result_proteomics %>% filter(Y == "vWF") %>% ggplot(aes(x=Model, y= Rho, fill=Model) ) + 
  geom_boxplot() + ggforce::geom_sina(alpha=0.2)  +  theme_bw() + coord_flip() + geom_hline(yintercept = 0) + scale_fill_manual(values = c(met.brewer(name="Klimt", n=3,type="discrete")) ) + facet_wrap(~Y)


