library(tidyverse)
library(vegan)
library(ape)
library(patchwork)
library(stats)
library(viridis)


#Data and covariates reading
Data = readRDS(file = "Data/Immuno_matrix_postSelection.rds")
#Removal of samples that belong to CeD cohort
Cohort = sapply(Data$ID, function(x){ str_split(x, "_")[[1]][1] } ) %>% as_vector() -> cohorts
Data %>% mutate(Cohort = cohorts) %>% filter(!Cohort == "34") -> Data
set.seed(99)

Get_clusters = function(Data){
        Data %>% select_if(~ !any(is.na(.)))  -> Data
        kmeans(Data, 2)$cluster -> Clusters
        return(Clusters)
}


#Create function to do PCA and extract eigen values / eigen vectors
DO_PCA = function(Data, permanova=F, Cluster=NULL){
        Data %>% select(-Cohort) -> Data
        sapply(Data$ID, function(x){ str_split(x, "_")[[1]][2] }) -> New_ID
        Data$ID = New_ID
        Covariates = read_tsv("Data/Covariates_LLD&IBD.tsv")
	#Remove Followup samples
	Covariates %>% filter(!grepl("F", Sample_name)) -> Covariates
        Covariates %>% filter(ID %in% Data$ID ) -> Covariates
        Covariates %>% drop_na() -> Covariates

        Data %>% filter(ID %in% Covariates$ID) -> Data
        Data  %>% select(-ID) -> Data2

        print( paste(c("Data size: ", as.character(dim(Data)[1]) ), collapse = "" ) )
        print( paste(c("Data with Metadata: ", as.character(dim(Data2)[1]) ), collapse = "" ) )

        #Get_clusters(Data2) -> Cluster_belonging
        #PCA
        Data2 %>% select_if(~ !any(is.na(.))) -> Data2
        PCA <- rda(Data2, scale = FALSE)
        scores(PCA,seq(100))$sites %>% as_tibble() %>% mutate(ID = Data$ID) -> First_100_PCs
        left_join(First_100_PCs, Covariates) -> First_100_PCs
        #Get eigenvalues, divided by its sum and multiply by 100 to get the % of variance explained
        #Variance_per_axis = (as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) * 100
        summary(PCA) -> PCA_s
        #Get normalized Eigenvalues, and multiply by 100 to get the % of variance explained
        PCA_s$cont$importance["Proportion Explained",] * 100 -> Variance_per_axis
        #Get PCs
        PCs = PCA_s$sites
	PCs %>% as_tibble() %>% mutate(ID = Data$ID) -> PCs
        left_join( PCs, Covariates) -> PCs    
    
        #Scree plot
        PCA_s$cont$importance[c("Cumulative Proportion", "Proportion Explained"),] %>% t() %>%  as.data.frame() %>% rownames_to_column("PC")  -> Scree_input
        Scree_input %>% mutate(PC = factor(PC, levels= Scree_input$PC) ) -> Scree_input
        Scree_input %>% ggplot(aes(x= PC, y = `Proportion Explained`) )   + geom_bar(stat="identity") + theme_bw() + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank()) -> Scree_plot
        #PCs needed to reach 90%
        Scree_input %>% filter(`Cumulative Proportion` >= 0.9) %>% head(1) -> Min_PC

        #Check major contributors of PCs
        scores(PCA,seq(100))$species %>% as.data.frame()  %>% rownames_to_column("order") %>% arrange(desc(abs(PC2))) -> Top_PC
        #order,pos,len_seq,aa_seq,full.name
        read_csv("Data/df_info_AT_v2_2_SZ_New.csv") %>% select(order, `full.name`) -> Annotation
        left_join(Top_PC, Annotation) -> Top_PC
        print("Top probes PC2")
        print(as_tibble(Top_PC))
        select(Top_PC, c(order, `full.name`)) %>% as_tibble() %>% head(n=10) %>% print()


        #Do PCA
	if (is.null(Cluster)){
		Get_clusters(as_tibble(First_100_PCs[,1:10])) -> Cluster_belonging
        } else { Cluster_belonging = Cluster ; print(Cluster_belonging) }
	PCs %>% as_tibble() %>%  mutate(Cluster = as.factor(Cluster_belonging)) %>%
        ggplot() + geom_point(aes(x=PC1, y=PC2, col=Cluster)) + theme_bw() + xlab(paste(c("PC1(", as.character(round(Variance_per_axis[1], 2)), "%)"), collapse="")) + ylab(paste(c("PC2(", as.character(round(Variance_per_axis[2], 2)), "%)"), collapse="")) +  scale_color_manual(values = c("#440154FF","#21908CFF") )  -> PCA_plot
	Cluster_belonging = tibble(ID = First_100_PCs$ID, Cluster = Cluster_belonging)

	#Do PERMANOVA
	if ( permanova == T){
	PCs %>% as_tibble() %>%  mutate(Cluster = as.factor(Cluster_belonging)) -> Input_PCA
	saveRDS(Input_PCA, "Input_figures/Figure_1C.rds")
	print("computing distance")
	Distance = vegdist( Data2 , method="jaccard")
	print("permanova")
	if (length( unique(Covariates$Disease_status) ) >1 ){
		adonis2(formula = Distance ~ plate_id + Age + Sex + Disease_status , data= left_join(Data,Covariates) , na.action = na.exclude, permutations = 2000 ) -> permanova 
	}else{
		adonis2(formula = Distance ~ plate_id + Age + Sex, data= left_join(Data,Covariates) , na.action = na.exclude, permutations = 2000 ) -> permanova
	}
	print(permanova)
	}
	Component_importance =  as.data.frame(PCA_s$cont$importance) %>% rownames_to_column("Legend") %>% as_tibble()
        return( list(PCA_plot, Min_PC, Scree_plot, Top_PC, Component_importance, First_100_PCs, Cluster_belonging) )

}


print("All data")
DO_PCA(Data) -> Overall_results

print("Variance explained by first two components")
Overall_results[[5]][,1:2] %>% print()

print("LLD")
DO_PCA( filter(Data, Cohort == "32"), permanova=F) -> LLD_results
PC_100 = LLD_results[[6]]
Clusters = LLD_results[[7]]
write_tsv(PC_100,"Results/Ordination/Top_100_PCs.tsv")

write_tsv(Clusters,"Results/Ordination/Clusters.tsv")
q()
#LLD without CMV
CMV = read_csv("Data/OLD/CMV_data.csv")
CMV = colnames(CMV)[3:dim(CMV)[2]]
Subset = select(Data, -one_of(CMV))
print(Clusters)
DO_PCA( filter(Subset, Cohort == "32"),  permanova=F, Cluster = Clusters$Cluster) -> LLD_results_noCMV
ggsave("Results/Ordination/PCA_LLD_noCMV.pdf", LLD_results_noCMV[[1]])
print("IBD")
DO_PCA( filter(Data, Cohort == "33") ) -> IBD_results

#Save PCA, Min_number of PCs, scree plot and Top probes for PC1
SAVE = T
if (SAVE == T){
#1. PCA
ggsave("Results/Ordination/PCA_LLD.pdf", LLD_results[[1]])
ggsave("Results/Ordination/PCA_IBD.pdf", IBD_results[[1]])
ggsave("Results/Ordination/PCA_LLD&IBD.pdf", Overall_results[[1]])
ggsave("Results/Ordination/PCA_LLD_noCMV.pdf", LLD_results_noCMV[[1]])
#2. Min number of PCs
write_tsv(Overall_results[[2]],"Results/Ordination/Number_PCs_90perc_var_LLD&IBD.tsv")
write_tsv(LLD_results[[2]], "Results/Ordination/Number_PCs_90perc_var_LLD.tsv")
write_tsv(IBD_results[[2]], "Results/Ordination/Number_PCs_90perc_var_IBD.tsv")
#3. Scree plot
ggsave("Results/Ordination/Scree_LLD.pdf", LLD_results[[3]])
ggsave("Results/Ordination/Scree_IBD.pdf", IBD_results[[3]])
ggsave("Results/Ordination/Scree_LLD&IBD.pdf", Overall_results[[3]])
#Top probes PC1
write_tsv(Overall_results[[4]],"Results/Ordination/Top_loads_LLD&IBD.tsv")
write_tsv(LLD_results[[4]],"Results/Ordination/Top_loads_LLD.tsv")
write_tsv(IBD_results[[4]],"Results/Ordination/Top_loads_IBD.tsv")
write_tsv(LLD_results_noCMV[[4]], "Results/Ordination/Top_loads_LLD_noCMV.tsv")
}
Overall_results[[6]] -> First_100_PCs
write_tsv( as_tibble(First_100_PCs), "Results/Ordination/LLD_PCA")
colnames(First_100_PCs)[grepl("PC", colnames(First_100_PCs))] -> PC_names
PC_association = tibble()
for (PC in PC_names){
        formula =  paste(c(PC,"~Age+Sex+Disease_status" ), collapse="")
        summary(lm(formula, First_100_PCs)) -> model_r
        as.data.frame(model_r$coefficients)["Disease_statusIBD", "Pr(>|t|)"] -> P_IBD
        formula2 = paste(c(PC,"~Age+Sex+Probe_n" ), collapse="")
	formula3 = paste(c(PC,"~Probe_n" ), collapse="")
        summary(lm(formula2, First_100_PCs)) -> model_r2
        as.data.frame(model_r2$coefficients)["Probe_n", "Pr(>|t|)"] -> P_Nprobes
        rbind(PC_association, tibble(PC=PC, P_IBD = P_IBD, P_N_probes=P_Nprobes) ) -> PC_association
}

print("R2 of probe number in PC1")
print(summary(lm("PC1 ~ Probe_n", First_100_PCs))$r.squared)
write_tsv(PC_association, "Results/Ordination/PCs_vs_IBD.tsv")


First_100_PCs %>% ggplot() + geom_point(aes(x=PC3, y=PC20, col=Disease_status)) +  scale_color_manual(values = c("#440154FF","#21908CFF") ) + theme_bw() -> PCA_IBD
ggsave("Results/Ordination/PCA_LLD_vs_IBD.pdf", PCA_IBD)
write_tsv( as_tibble(Overall_results[[5]]), "Results/Ordination/Variance_explained.tsv")
write_tsv( as_tibble(LLD_results[[5]]), "Results/Ordination/Variance_explained_LLD.tsv")
write_tsv( as_tibble(IBD_results[[5]]), "Results/Ordination/Variance_explained_IBD.tsv")


