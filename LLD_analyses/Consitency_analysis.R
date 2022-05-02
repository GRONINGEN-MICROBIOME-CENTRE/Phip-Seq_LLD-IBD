args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(viridis)
library(ggforce)
if (length(args) == 1){ 
	if (args[[1]] == "Repeat"){
		Prepare = T
	} else { Prepare = F }
} else { Prepare = F}

if (Prepare == T){
	#Get data with followup
	
        Cov = read_tsv("Data/Covariates.tsv") #ID , sample_id
        Data = readRDS("Data/Immuno_matrix_postSelection.rds") #read_tsv("../Data/Immuno_matrix_Alex.tsv") # ID "32_1039563132"

        Data %>% filter(grepl("32_", ID)) -> Data
        Data$ID = str_remove(Data$ID, '"')
        Data$ID = as.vector(sapply(Data$ID, FUN= function(x){ str_split(x, "_")[[1]][2] } ))
        Data$sample_id = Data$ID
        select(Data, -ID) -> Data
        print(Data)
        #apply(select(Data, -sample_id) , 2, FUN = function(x){ x = as.numeric(x) ; sum(x, na.rm = T )/length(x)  } )  %>% as.data.frame() %>% rownames_to_column("Probe") -> D
        #write_tsv(D, "Prevalence_antibodies_total.tsv")
        Cov %>% filter(sample_id %in% Data$sample_id) -> Cov
        left_join(Cov, Data, by = "sample_id") -> DD

        DD %>% filter(grepl("LL", ID)) -> DD
        DD %>% mutate(Time_point = ifelse(grepl("_F", ID), "FA", "BL")) -> DD
        DD$ID =  as.vector(sapply(DD$ID, FUN= function(x){ str_split(x, "_F")[[1]][1] } ))
        DD_FA = filter(DD, Time_point == "FA")
        DD_BL = filter(DD, Time_point == "BL")


        DD_FA %>% filter(ID %in% DD_BL$ID) -> DD_FA
        DD %>% filter( ID %in% DD_FA$ID) -> DD

        saveRDS(DD, file = "Consitency_LLD/Data/Data_analysis.rds")
}else{
        # Restore the object
        DD = readRDS(file = "Consitency_LLD/Data/Data_analysis.rds")
}
##Add CMV info as covariate
read_tsv("Data/Covariates_LLD&IBD.tsv") -> Cov
read_tsv("Results/Ordination/Clusters.tsv") -> Clusters
left_join(Cov, Clusters) %>% mutate(ID = Sample_name) %>% select(ID, Cluster) -> Clusters
left_join(DD, Clusters) -> DD
##Add EBV info as covariate
read_tsv("Data/EBV_clusters.txt") -> Clusters_EBV
Clusters_EBV$ID = as.vector(sapply(Clusters_EBV$ID, FUN= function(x){ str_split(x, "_")[[1]][2] } ))
left_join(Cov, Clusters_EBV) %>% mutate(ID = Sample_name) %>% select(ID, EBV_cluster) %>% mutate(EBV_cluster=1*EBV_cluster) -> Clusters_EBV
left_join(DD, Clusters_EBV) -> DD


#Check that all samples included have two timepoints
DD %>% group_by(Time_point) %>% summarise(n())


Names = c("sample_id", "Sex", "Age", "plate_id", "Time_point","Disease_status", "Cluster", "EBV_cluster")

apply(select(DD, -c(Names, "ID")), 1, function(x){ x[!is.na(x)] ; sum(x) } ) ->  N_ab
DD  %>% mutate(N_ab = N_ab) %>% select(ID, plate_id, Disease_status, N_ab, Time_point) %>% arrange(ID) -> DD_n
print(DD_n)
#Compare Number of antibodies enriched between time points
wilcox.test( filter(DD_n, Time_point == "BL")$N_ab , filter(DD_n, Time_point == "FA")$N_ab, paired=T ) -> test
summary(lm(N_ab ~ Time_point + plate_id, DD_n)) %>% print()
print(test)
DD_n %>% drop_na() %>%  group_by(Time_point) %>% summarise(M = mean(N_ab)) %>% print()
DD_n %>% group_by(Time_point, plate_id) %>% summarise(n()) %>% arrange(plate_id) %>% print()


#Participant's distances

f <- function (i, j, dist_obj) {
  #Get location of the distance between two objects, given a distance matrix
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
  k <- (2 * n - j) * (j - 1) / 2 + (i - j)
  k[!valid] <- NA_real_
  k
}


Prevalence_f= function(x, Thresh = 0.1){
        V = sum(x, na.rm= T)/length(x)
        return( V >= Thresh)
}

AB_profile = select(DD, -c("ID", Names)) -> AB_profile
#Keep_ab = apply(AB_profile, 2, FUN = Prevalence_f) #Appplied on samples with two time points
#Keep_ab = colnames(AB_profile)[Keep_ab]

#AB_profile %>% select(Keep_ab) -> AB_profile

Get_distances = function(DD, Distance_matrix){
	#Get distance between time points per sample pair belonging to same person
	Distances_common_tibble = tibble()
	for (People in unique(DD$ID)){
        	Look_up = which(DD$ID == People)
	        k = f(Look_up[1], Look_up[2], Distance_matrix)
        	if (is.na(k)) { k = f(Look_up[2], Look_up[1], Distance_matrix) }
	        if (is.na(k)) { next }
        	Distance = Distance_matrix[k]
	        Distances_common_tibble = rbind(Distances_common_tibble, tibble(ID = People, k= k, Distance= Distance_matrix[k] ))
	}

	Distances_common = Distances_common_tibble$k
	Common = Distance_matrix[Distances_common]
	Not_common =  Distance_matrix[-Distances_common]
	return( list(Common, Not_common, mean(Common)  ) )
}

Distance_Consistency = function(Distance_matrix,DD, suffix=""){
	#1. Get distances between longitudinal samples
	Get_distances(DD, Distance_matrix) -> Not_permuted
	Common = Not_permuted[[1]] ;  Not_common = Not_permuted[[2]]
	
	t = Not_permuted[[3]] #Mean distance, statistic to compare to the null distribution of permutation means
	#2. Use permutations to estimate significance
	Permutations = 2000
	Null_distribution = c()
	#Loop where permutations happen and distances are obtaine. We use mean distance between permuted longitudinal samples to produce a null distribution
	for (p in Permutations){
		DD_p = DD
		DD_p$ID = DD_p$ID[sample(nrow(DD_p)) ]
		Get_distances(DD_p, Distance_matrix) -> permuted
		Null_distribution = c(Null_distribution, permuted[[3]] ) #[[3]] is mean of distances between longitudinal samples
	}
	#Stat of a comparison without using permutations
	wilcox.test(Common, Not_common) -> Stat_w
	P_diff = Stat_w$p.value
	print( paste(c("Wilcox distance related-unrealated: ", as.character(P_diff)), collapse="") )
	#P-value when compared to the null distribution of permuted stats 	
	p = length(Null_distribution[Null_distribution<=t])/Permutations
	print( paste(c("Pvalue to null distribution of means from permutations: ", as.character(p)), collapse="") )
	#Save distance between samples
	Distance_tibble = rbind( tibble(Distance = Common, Longitudinal = T), tibble(Distance = Not_common, Longitudinal = F) )
	saveRDS(Distance_tibble, paste0("Input_figures/Fig1D",suffix,".rds"))
	#Plot distance distributions
	Distance_tibble %>% ggplot(aes(x=Longitudinal, y = Distance)) + ggforce::geom_sina(aes(col=Longitudinal )) +theme_bw() + scale_color_manual(values = c("#6a3d9a", "#fb9a99")) + ylab("Distance") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) -> Distance_box
	Distance_tibble %>% ggplot(aes(x=Distance, fill = Longitudinal)) + geom_density(alpha=0.5) +theme_bw() + scale_colour_viridis_d("viridis") -> Distance_distribution
	ggsave(paste0("Consitency_LLD/Results/Plots/Distance_box",suffix,".pdf"), Distance_box) ; ggsave(paste0("Consitency_LLD/Results/Plots/Distance_distribution",suffix,".pdf"), Distance_distribution)
	
	#3. What factors determine change?
	print("Identify factors that predict change")
	Distance_tibble %>% filter(Longitudinal == T) -> Distances_common_tibble
	Distances_common_tibble %>% arrange(desc(Distance))
	left_join( mutate(Distances_common_tibble, ID = unique(DD$ID)) , select(DD, c(ID, Sex, Age, Cluster,EBV_cluster) ) , by="ID") %>% drop_na() -> Distances_common_tibble
	Distances_common_tibble %>% mutate(Age_n = Age-min(Age) ) -> Distances_common_tibble
	lm(Distance ~ Sex+ Age_n + Cluster + EBV_cluster, Distances_common_tibble) -> M
	summary(M) %>% print()
	write_tsv(Distances_common_tibble, paste0("Consitency_LLD/Distance_per_subject",suffix,".tsv") )
}

#Run analysis using Jaccard distance
vegan::vegdist(AB_profile, method="jaccard", na.rm = T) -> Distance_matrix_jaccard
Distance_Consistency(Distance_matrix_jaccard, DD)

q()
print("Boostrap analysis")
set.seed(54)
Boostrap_distances = list()
for (i in c(0.2, 0.4, 0.6, 0.8)){
	print( paste0("Boostrap ", as.character(i)))
	sample( seq(1,dim(AB_profile)[2]), i*dim(AB_profile)[2], replace=F ) -> Indexes 
	AB_profile[,Indexes] -> AB_boostrap
	print(dim(AB_boostrap))
	vegan::vegdist(AB_boostrap, method="jaccard", na.rm = T) -> Distance_matrix_boostrap
	Distance_Consistency(Distance_matrix_boostrap, DD, paste0("_jaccardBoostrap_", as.character(i*100) ) )
	Name = paste0("Permutation_", as.character(i))
	Boostrap_distances[[Name]] = Distance_matrix_boostrap
}

for (Entry in Boostrap_distances){
	for (Entry2 in Boostrap_distances){
		vegan::mantel(Entry, Entry2) -> R
		R$statistic %>% print()
	}
}
#Stratify analysis based on Cluster (CMV status)
AB_profile[DD$Cluster == 2,] -> AB_CMV
vegan::vegdist(AB_CMV, method="jaccard", na.rm = T) -> Distance_matrix_jaccard_cmv
Distance_Consistency(Distance_matrix_jaccard_cmv, DD[DD$Cluster == 2,], "_jaccard_CMV")
AB_profile[DD$Cluster == 1,] -> AB_noCMV
vegan::vegdist(AB_noCMV, method="jaccard", na.rm = T) -> Distance_matrix_jaccard_nocmv
Distance_Consistency(Distance_matrix_jaccard_nocmv, DD[DD$Cluster == 1,], "_jaccard_noCMV")

#Run analysis using Manhattan distance
vegan::vegdist(AB_profile, method="manhattan", na.rm = T) -> Distance_matrix_manhattan
#Matrix-matrix comparison
vegan::mantel(Distance_matrix_jaccard, Distance_matrix_manhattan) %>% print()
Distance_Consistency(Distance_matrix_manhattan, DD, "_manhattan")

###Analysis per antibody

Keep_ab = colnames(AB_profile)
Prepare = T
if (Prepare  == T){
        #By antibody consistency
        Antibody_table = tibble()
        for ( Antibody in Keep_ab  ){
                if (Antibody %in% c(Names, "ID")){ next }
                DD %>% select(c(Antibody)) %>% as.vector() %>% as_vector() -> AB

                DD %>% select(Names) %>% mutate(Anb = AB) -> Antibody_data
        
                Antibody_data1 =   filter(Antibody_data, Time_point == "BL") %>% arrange(sample_id)
                Antibody_data2 = filter(Antibody_data, Time_point == "FA") %>% arrange(sample_id)

                Antibody_data1 %>% mutate( BL = Antibody_data1$Anb , FA = Antibody_data2$Anb  ) %>% select(-c("Time_point", Anb)) ->Antibody_data_wide
                #Antibody_data %>% spread("Time_point", Antibody) -> Antibody_data_wide
        
        
                NA_rate_BL = sum(is.na(Antibody_data_wide$BL)) / dim(Antibody_data_wide)[1]
                NA_rate_FA = sum(is.na(Antibody_data_wide$FA)) / dim(Antibody_data_wide)[1]

                Antibody_data_wide %>% drop_na() -> Antibody_data_wide
        
                Prevalence_BL = sum(Antibody_data_wide$BL , na.rm= T) / dim(Antibody_data_wide)[1]
                Prevalence_FA = sum( Antibody_data_wide$FA , na.rm= T) / dim(Antibody_data_wide)[1]

                #C = cor(Antibody_data_wide$BL, Antibody_data_wide$FA, method = "pearson" ) #"p"
                Consistency = sum( Antibody_data_wide$BL == Antibody_data_wide$FA ) / dim(Antibody_data_wide)[1] 
		Gains = sum(Antibody_data_wide$FA > Antibody_data_wide$BL)
		Loss = sum(Antibody_data_wide$FA < Antibody_data_wide$BL)
		Gains = Gains/(Gains + Loss)		

                N_T = tibble("Antibody" = Antibody, "Consistency" = Consistency, "Gain(%)" = Gains, "Prevalence_BL" = Prevalence_BL, "Prevalence_FA" = Prevalence_FA, "NA_rate_BL" =NA_rate_BL, "NA_rate_FA"=NA_rate_FA )
		        
                Antibody_table = rbind(Antibody_table, N_T)     
}
        saveRDS(Antibody_table, file = "Consitency_LLD/Results/Antibody_consistency.rds")

} else { Antibody_table =  readRDS(file ="Consitency_LLD/Results/Antibody_consistency.rds")  }

Antibody_table %>% ggplot(aes(x=Consistency )) + geom_density() + theme_bw() -> Plot_consistency
Antibody_table %>% group_by(Mean_change >= 0 ) %>% summarise(n()) %>% print()
ggsave("Consitency_LLD/Results/Plots/Consistency_distribution.png", Plot_consistency)
write_tsv(Antibody_table, "Consitency_LLD/Results/Consistency_by_ab.tsv") 



