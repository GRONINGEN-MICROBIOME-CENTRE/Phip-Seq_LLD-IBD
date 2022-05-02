#S. Andreu-Sanchez
#Uses WGCNA in a presence/absence profile to identify modules of highly correlated peptides
#Also: Module correlation. Module composition plotting. Network plotting.
#Sequence similarity clustring

library(WGCNA) # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/OverviewWGCNA.pdf
library(tidyverse)
library(readxl)
library(igraph)
library(vegan)
library(ggtree)

set.seed(1299)

#Colors, 25 different colors easily identified by human eye
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
        "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
#Functions for WGCNA part
process =  function(Data){
        #cleanup of the data before other functions can be applied
        Data %>% select(-ID) -> Data
        Data[ , apply(Data, 2, function(x) !any(is.na(x)))] -> Data
        return(Data)
}
choose_power = function(Data, Distance_m = F, Distance= NULL){
        #Choosing ideal power function.
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        if (Distance_m == F){
                sft = pickSoftThreshold(Data, powerVector = powers, verbose = 5, corFnc = cor )
        } else{
                sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(Distance), powerVector = c(c(1:10), seq(from = 12, to=20, by=2)), verbose=5 )
                
        }
        sizeGrWindow(9, 5)
        par(mfrow = c(1,2));
        cex1 = 0.9;

        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
                main = paste("Scale independence"));
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col="red");
        abline(h=0.90,col="red")

        plot(sft$fitIndices[,1], sft$fitIndices[,5],
                xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
                main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

Do_analysis2 =  function(Data, k=7, cutoff=30, minCorrMerge=0.5){
        #It seems that the automatic block methods is ?always? using unsigned adjacency metrics. So for reproducing paper methos, type shoudl be unsigned
        #While a correlation-based method accounts for co-absence, using Jaccard instead to produce an adjancecy metric does not
        #Jaccard-based clustering shows a certain cluster all over the place. While coventional cutting methods of hierarchical clustering cutting
        #use a unique cutoff (below which we call cluster), WGCNA is using a dynamic tree cut ( https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/BranchCutting/ )
        #which should be less affected by outliers. In this case, we are using dynamic hybrid, which makes clusters appear all over the place.
        # https://academic.oup.com/bioinformatics/article/24/5/719/200751?login=true paper discribing the cutting
        
        #1. Adjacency
        #Used in V1 paper
        adjacency_metric_cor = adjacency(Data, power=k, type = "unsigned")
        #Using a binary metric: jaccard.
        adjacency_metric =  (1 - as.matrix(vegdist(t(Data), method = "jaccard")) )^k
        
        #2. Getting TOM
        TOM = TOMsimilarity(adjacency_metric)
        dissTOM = 1-TOM
        TOM_cor = TOMsimilarity(adjacency_metric_cor)
        dissTOM_cor = 1-TOM_cor
        #Compare TOMs
        vegan::mantel(dissTOM, dissTOM_cor) #Stat r: 0.23, P<0.001 (999 perm)
        
        #3. Hierarchical clustering
        geneTree = hclust(as.dist(dissTOM), method = "average") 
        geneTree_cor = hclust(as.dist(dissTOM_cor), method = "average") 
        #Plotting of similarity of the different features
        plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
        #4. Dynamic tree cut  
        # Module identification using dynamic tree cut: Add a minimum number of peptides per module: cutoff
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                    deepSplit = 2, pamRespectsDendro = FALSE,
                                    minClusterSize = cutoff  )
        Modules_corr = cutreeDynamic(dendro = geneTree_cor, distM = dissTOM_cor,
                                     deepSplit = 2, pamRespectsDendro = FALSE,
                                     minClusterSize = cutoff  )
        
        #Some summaries, number of samples per module and give colors to those samples
        table(dynamicMods)
        dynamicColors = labels2colors(dynamicMods)
        table(dynamicColors)
        # Plot the dendrogram and colors underneath
        sizeGrWindow(8,6)
        plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Gene dendrogram and module colors")
        #5. Eigengenes and merging of modules which eigengene is closer than a  certain (correlation) threshold with the others
        #Eigengenes to compute similarity between modules
        MEList = moduleEigengenes(Data, colors = dynamicColors)
        MEs = MEList$eigengenes
        # Calculate dissimilarity of module eigengenes
        MEDiss = 1-cor(MEs)
        # Cluster module eigengenes
        METree = hclust(as.dist(MEDiss), method = "average")
        # Plot the result
        sizeGrWindow(7, 6)
        plot(METree, main = "Clustering of module eigengenes",
             xlab = "", sub = "")
        #Set a cutoff for merging modules into one (minCorrMerge), display and merge
        abline(h=minCorrMerge, col = "red")
        # Call an automatic merging function
        merge = mergeCloseModules(Data, dynamicColors, cutHeight = minCorrMerge, verbose = 3)
        mergedColors = merge$colors
        # Eigengenes of the new merged modules:
        mergedMEs = merge$newMEs
        #Comparison between modules before and after merging
        sizeGrWindow(12, 9)
        plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
        #Some renaming
        moduleColors = mergedColors
        colorOrder = c("grey", standardColors(50));
        moduleLabels = match(moduleColors, colorOrder)-1;
        MEs = mergedMEs
        
        
        table(moduleLabels) -> overview
        #Visualize similarities between the eigen representation of the modules
        Correlation_eigen =  cor(MEs)
        Correlation_eigen[Correlation_eigen>0.999] = NA
        pheatmap::pheatmap( Correlation_eigen ) -> Plot_heatmap
        #Plotting distance between all probes,and coloring by module
        plotTOM = as.matrix(TOM)
        diag(plotTOM) = NA;
        sizeGrWindow(9,9)
        TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all probes") -> Plot_comparison
        
        tibble(Probe = colnames(Data), Cluster = moduleLabels, Cluster_cor =as.vector(Modules_corr) ) -> Probes_and_cluster

        
        return(list(moduleLabels, moduleColors, MEs, geneTree, Plot_tree, Plot_comparison,Probes_and_cluster))
}
Do_analysis = function(Data, k=7, cutoff=30){
        #Build network in one command. 
        
        
        net = blockwiseModules(Data, power = k,corType = "pearson",
                       TOMType = "unsigned", minModuleSize = cutoff,
                       reassignThreshold = 0, mergeCutHeight = 0.5,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F, randomSeed = 1234, saveTOMFileBase = "unsignedTOM",
                       verbose = 3)

        #Distance matrix is saved in BlockwiseTOM-block.1.RData ; can be loaded by load(lockwiseTOM-block.1.RData)
        #load("BlockwiseTOM-block.1.RData") #It is called TOM

        #Data of interest
        moduleLabels = net$colors  #classification in labels
        moduleColors = labels2colors(net$colors) #color for the labels in figure
        MEs = net$MEs #Eigengenes
        geneTree = net$dendrograms[[1]] #dendogram
        #Overview of clusters
        table(moduleLabels) -> overview
        #Plot tree
        #plotDendroAndColors(geneTree, moduleColors, "Module colors" ,dendroLabels = F, addGuide = T) -> Plot_tree
        Plot_tree = NULL
        #Visualize similarities between the eigen representation of the modules
        Correlation_eigen =  cor(MEs)
        Correlation_eigen[Correlation_eigen>0.999] = NA
        pheatmap::pheatmap( Correlation_eigen ) -> Plot_heatmap
        #Plotting distance between all probes,and coloring by module
        #plotTOM = as.matrix(TOM)
        #diag(plotTOM) = NA;
        #sizeGrWindow(9,9)
        #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all probes") -> Plot_comparison
        Plot_comparison = NULL #Making heatmap takes really long time
        tibble(Probe = colnames(Data), Cluster = moduleLabels) -> Probes_and_cluster
        
        return(list(moduleLabels, moduleColors, MEs, geneTree, Plot_tree, Plot_comparison,Probes_and_cluster, net))
}       


#Read data and remove probes with NA
Data = readRDS("~/Desktop/Immuno_matrix_postSelection.rds") #Complete matrix of 0/1s + IDs, script 01
process(Data) -> Data_all
process( filter(Data, grepl("32_", ID))) -> Data_LLD #Removal of samples that are not LLD. Might consider removing longitudinal samples.
read_csv("~/Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset//2821probes_annotation.csv") -> Annotation #Peptide annotation
colnames(Annotation)[1] = "Probe"

Check_power_different_distances = function(Data_LLD){
        ##LLD##
        choose_power(Data_LLD)
        #Jaccard
        choose_power(Data_LLD, Distance_m = T, Distance = vegan::vegdist(Data_LLD, "jaccard") )
        #Kulczynski, this one is really slow. With a power of 20 it is around 0.6, although it seems that the R2 increases linealy with the power, is not recommended to go really high
        choose_power(Data_LLD, Distance_m = T, Distance = prabclus::kulczynski(Data_LLD) )
        
}
Other_analyses = function(Data, Data_LLD){
        #Analyses not included in the paper. Different number of min peptides per module and including IBD samples
        
        ###OVERALL (LLD+IBD)###
        #Check power
        choose_power(Data_all) #It is recommended to use a number aboce R2 of 0.9, the highest possible before saturation, which in this case is 7
        Do_analysis(Data_all, 7, 10) -> results_all
        ME_all = results_all[[3]]
        Probes_and_cluster = results_all[[7]]
        #Check which are the probes in the same groups
        Probes_and_cluster %>% filter(Cluster != 0) -> Groups
        left_join(Groups, Annotation) -> Groups
        #1 --> CMV , 2 ---> Flagellin, 3 --> Bacteria (allergen?)
        Groups %>% select(Cluster, Taxa, Description,`comments, synonim, blast, details`) %>% filter(Cluster == 1) %>% print(n=100) #Cluster CMV, turquoise
        Groups %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 2) %>% print(n=100) #Cluster CMV, turquoise
        Groups %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 3) %>% print(n=100) #Cluster CMV, turquoise

        Do_analysis(Data_LLD, 7) -> results_LLD
        Probes_and_cluster_LLD = results_LLD[[7]]
        Probes_and_cluster_LLD %>% filter(Cluster != 0) -> Groups_LLD ;left_join(Groups_LLD, Annotation) -> Groups_LLD

        Groups_LLD %>% select(Cluster, Taxa, Description,`comments, synonim, blast, details`) %>% filter(Cluster == 1) %>% print(n=100) #Cluster CMV, turquoise
        Groups_LLD %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 2) %>% print(n=100) #Cluster CMV, turquoise
        Groups_LLD %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 3) %>% print(n=100) #Cluster CMV, turquoise

        results_LLD[[2]][results_LLD[[1]] == 1] #light blue
        results_LLD[[2]][results_LLD[[1]] == 2] #dark blue
        results_LLD[[2]][results_LLD[[1]] == 3] #brown\

        Groups_LLD %>% select(Cluster, Probe,aa_seq) %>% filter(Cluster ==3) -> Cluster_heterogenous
        Cluster_heterogenous %>% write_tsv("Desktop/Cluster_heterogenous.tsv")
        write_tsv(Groups_LLD, path = "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Module_belonging_LLD.tsv")

        #Repeat on different cluster number
        #With 20 there is one additional module
        Do_analysis(Data_LLD, 7, 20)  -> results_LLD2
        Probes_and_cluster_LLD2 = results_LLD2[[7]] ; Probes_and_cluster_LLD2 %>% filter(Cluster != 0) -> Groups_LLD2 ;left_join(Groups_LLD2, Annotation) -> Groups_LLD2
        Groups_LLD2 %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 4) %>% print(n=100)
        #New cluster consists of 26 probes. No clear taxa or domain, mainly bugs
        results_LLD2[[2]][results_LLD2[[1]] == 4] #yellow cluster
}

##1. Identify modules with at least 10 probes per module
#Check_power_different_distances(Data_LLD)
#With minimum 10 pepties per module
Do_analysis(Data_LLD, 7, 10)  -> results_LLD3
Net =  results_LLD3[[8]]
Probes_and_cluster_LLD3 = results_LLD3[[7]] ; Probes_and_cluster_LLD3 %>% filter(Cluster != 0) -> Groups_LLD3 ;left_join(Groups_LLD3, Annotation) -> Groups_LLD3
table(results_LLD3[[2]]) #22 clusters
results_LLD3[[3]] %>% as_tibble() %>% mutate(ID = filter(Data,grepl("32_", ID))$ID, .before=1) %>% select(-ME0) -> EigenGenes


#Need to run overall analysis to find a network common with IBD and LLD
ME_all %>% as_tibble() %>% mutate(ID = Data$ID, .before=1) %>% select(-ME0)  -> EigenGenesAll
left_join(EigenGenes, EigenGenesAll, by="ID", suffix=c("", "_all")) %>% select(-ID)  %>% cor() -> Cor_all #%>% heatmap()
hclust(as.dist( 1- Cor_all), method="average") -> CLUSTER
cutree(CLUSTER, h = 0.05) -> CLUSTER
for (C in CLUSTER){
        CLUSTER[CLUSTER==C] -> Pair
        if (! length(Pair) == 2 ){ next }
        names(Pair)[grepl("_all", names(Pair))] -> All_cluster
        Number_change = names(Pair)[! grepl("_all", names(Pair))]
        str_replace(All_cluster, "_all", "") -> All_cluster
        EigenGenesAll[ paste0(Number_change, "_complete") ] = EigenGenesAll[,All_cluster]
}
select(EigenGenesAll, c("ID", colnames(EigenGenesAll)[grepl("_complete", colnames(EigenGenesAll))] )) -> EigenGenesAll
colnames(EigenGenesAll) = sapply(colnames(EigenGenesAll), function(x){ str_split(x, "_")[[1]][1] } )
write_tsv(EigenGenes, "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/EigenGenes.tsv")
write_tsv(EigenGenesAll, "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/EigenGenes_all.tsv")

#Save
write_tsv(Groups_LLD3, "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Cluster_AtLeast10.tsv")

########################
#Boostrap analysis######
########################
Do_boostrap_analysis = function(Data_LLD, Probes_and_cluster_LLD3, N_iterations=10){
        set.seed(9899)
        ##1. Run boostraps and save clusters
        Boostrap_clusters = tibble()
        #Number of permutations are the number of iterations to do. In each iteration 4 subsets (20%, 40%, 60% and 80%) of data is subsampled
        for (i in seq(N_iterations) ){
                print (paste0("Round number", as.character(i) ) )
                for (a in c(0.2, 0.4, 0.6, 0.8)){
                      print (paste0("Data subset", as.character(a) ) )
                      sample(seq(dim(Data_LLD)[1]),a*dim(Data_LLD)[1],replace=T)  -> Boost
                      Data_LLD[Boost,] -> Boost
                      Do_analysis(Boost, 7, 10) -> Results_boost
                      Clusters =  Results_boost[[7]]
                      Clusters %>% mutate(Percentage = a, Round = i) -> Clusters
                      rbind(Boostrap_clusters, Clusters) -> Boostrap_clusters
                }
        }
        ##2. Find out which clusters are related.
        Boostrap_clusters %>% mutate( ID = paste0(Percentage,"-", Round) ) -> Boostrap_clusters
        Similarity_clusters = tibble()
        #For each Cluster in the whole dataset, check which are the most similar clusters in each subset of data
        #For that, we compute Jaccard similarity (1-jaccard) and pick the top closest module. If no module is over 50% similarity it is assumed no matching cluster is present
        for (Cluster_n in unique(Probes_and_cluster_LLD3$Cluster)){
                Probes_and_cluster_LLD3 %>% mutate(Presence = (Cluster == Cluster_n))  -> Real
                for (Round_n in unique(Boostrap_clusters$ID) ){
                        Similarity_round = tibble()
                        Boostrap_clusters %>% filter(ID  == Round_n) -> Perm
                        for (Cluster_n2 in unique(Perm$Cluster)){
                                Perm %>% mutate(Presence = (Cluster == Cluster_n2)) -> Perm2
                                Similarity = 1 - vegdist( t(data.frame(as.numeric(Real$Presence), as.numeric(Perm2$Presence))), "jaccard")[1]
                                Similarity_round %>% rbind( tibble(Sim = Similarity, Cluster_real = Cluster_n, Cluster_perm = Cluster_n2 ))  -> Similarity_round       
                        }
                        which.max(Similarity_round$Sim) -> Max_d
                        which.max(Similarity_round$Sim[-Max_d] ) -> Max_d2
                        Similarity_round[Max_d,] -> Best_match
                        Similarity_round[Max_d2,] -> Second_best
                        Diff = Best_match$Sim - Second_best$Sim
                        Result_round = tibble(Cluster = Cluster_n, Match =  Best_match$Cluster_perm, Second_match = Second_best$Cluster_perm, Similarity = Best_match$Sim, Difference = Diff, Round = Round_n)
                        Similarity_clusters %>% rbind(Result_round) -> Similarity_clusters
                }
        }
        Similarity_clusters$Round %>% sapply(function(x){ str_split(x, "-")[[1]][1]  } ) -> Subset_perm
        Similarity_clusters %>% mutate(Subset = Subset_perm ) -> Similarity_clusters
        #Plot cluster similarity
        Similarity_clusters %>% ggplot(aes(x = as.factor(Cluster), y = Similarity)) + geom_boxplot(outlier.shape = NA) + theme_bw()+
                 geom_hline(yintercept = 0.5)  + coord_flip() + geom_jitter(aes(col=Subset)) -> Plot1
        Similarity_clusters %>% filter(Similarity>0.5) %>% ggplot(aes(x = as.factor(Cluster), y = Similarity)) + geom_boxplot(outlier.shape = NA) + 
                geom_jitter(aes(col=Subset)) + theme_bw() + coord_flip() -> Plot2
        ggsave("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Bootstrap/Cluster_cluster_top_similarity.pdf", Plot1)
        
        
        ##2.5 Change the cluster value to the ones where they match best.
        #Relatively long execution time ~10min
        Boostrap_clusters2 =  Boostrap_clusters
        library(doParallel)
        detectCores()
        registerDoParallel(6)
        start.time <- Sys.time()
        foreach( row_n = seq(nrow(Boostrap_clusters)), .combine = "rbind" ) %dopar% {
                Subset2 = Boostrap_clusters
                Subset2[row_n,] -> Boostrap_row
                Similarity_clusters %>% filter(Round == Boostrap_row$ID, Match == Boostrap_row$Cluster) %>% filter(Similarity>0.5) -> Match
                if(! dim(Match)[1] == 1 ){
                        C = NA
                }else{ C = Match$Cluster }
                Subset2[row_n,2] = C 
                return(Subset2[row_n,])
        } -> Boostrap_clusters2
        end.time <- Sys.time()
        end.time - start.time
        
        #Go row by row, check the cluster and the matching cluster for that boostrap and change the value of the clsuter to the one they match with
        #for (row_n in  seq(nrow(Boostrap_clusters))){
        #        Boostrap_clusters2[row_n,] -> Boostrap_row
        #        Similarity_clusters %>% filter(Round == Boostrap_row$ID, Match == Boostrap_row$Cluster) %>% filter(Similarity>0.5) -> Match
        #        if(! dim(Match)[1] == 1 ){
        #                C = NA
        #        }else{ C = Match$Cluster }
        #        Boostrap_clusters2[row_n,2] = C
        #}
        Boostrap_clusters2 %>% drop_na() -> Boostrap_clusters2
        
        ##3. Get consistency of clusters that are meant to represent the same
        # Done by counting the % of times the cluster in the complete dataset matches the homologous cluster
        Boostrap_results = tibble()
        for ( Peptide in Probes_and_cluster_LLD3$Probe ){
                Probes_and_cluster_LLD3 %>% filter(Probe == Peptide) -> Peptide_choice
                Boostrap_clusters2 %>% filter(Probe == Peptide) %>% group_by(Percentage, Cluster == Peptide_choice$Cluster) %>% summarise(N = n()) %>%
                        filter(`Cluster == Peptide_choice$Cluster` == T) %>%  mutate(Stability_perc = 100*N/N_iterations ) %>% dplyr::select(-c(`Cluster == Peptide_choice$Cluster`, N)) %>%
                         mutate(Peptide = Peptide_choice$Probe, cluster=Peptide_choice$Cluster ,.before=1) -> Res_boost
                rbind(Boostrap_results,Res_boost) -> Boostrap_results
        }
        Boostrap_results %>% spread(Percentage, Stability_perc) -> Boostrap_results
        Boostrap_results %>% gather(Subset, Consistency, 3:6, factor_key=TRUE) %>% #filter(cluster == 8) %>% drop_na() %>% print(n=20)
              drop_na() %>%  ggplot(aes(x=Subset, y=Consistency) ) + facet_wrap(~cluster) + geom_boxplot() + theme_bw() -> Plot3
        ggsave("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Bootstrap/Peptide_cluster_consistency.pdf", Plot3)
        
        
        Boostrap_results %>% gather(Subset, Consistency, 3:6, factor_key=TRUE) %>% group_by(cluster, Subset) %>% drop_na() %>% summarise(mean(Consistency), median(Consistency))
        
        #Save Bootstrap data
        write_tsv(Boostrap_results, path = "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Bootstrap/Boostrap_results.tsv")
        write_tsv(Boostrap_clusters2, path = "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Bootstrap/Boostrap_clusters.tsv")
        write_tsv(Boostrap_clusters, path = "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Bootstrap/Boostrap_raw_clusters.tsv")
        write_tsv(Similarity_clusters , path = "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Bootstrap/Cluster_similarity.tsv")

}

Do_boostrap_analysis(Data_LLD, Probes_and_cluster_LLD3, N_iterations=50)


#############Correlation between modules using eigengenes
Net$MEs %>% as_tibble() -> Eigengenes
cor(Eigengenes) -> Cor_modules
Cor_eig = tibble()
for (Ei in colnames(Eigengenes)){
        for (Ei2 in colnames(Eigengenes)){
                if (Ei == Ei2){ next }
                ID =paste(sort(c(Ei, Ei2)), collapse="-")
                if ( ID %in% Cor_eig$ID ){ next }
                cor.test(as_vector(select(Eigengenes, Ei)), as_vector(select(Eigengenes, Ei2))) -> test
                rbind(Cor_eig, tibble(ID=ID, P=test$p.value, Estimate=test$estimate)) -> Cor_eig 
                
        }
}
Cor_eig %>% mutate(P_b =p.adjust(P, "bonferroni")) %>%arrange(P_b) %>% filter(P_b < 0.05 ) %>% filter(!grepl("ME0", ID))
cor.test(Eigengenes$ME18, Eigengenes$ME17) ; cor.test(Eigengenes$ME18, Eigengenes$ME14)
Cor_eig %>% mutate(P_b =p.adjust(P, "bonferroni")) %>%arrange(P_b) %>% filter(! grepl("ME0", ID)) %>% filter(P_b < 0.05)

pheatmap::pheatmap(Cor_modules) #Modules 18, 17,14 and ?0? are related
#######################################

#####################################################
###Plot piechart of composition per cluster##########
#####################################################
#Load cluster with annotation
readxl::read_xlsx("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Correlation_vs_Similarity/Cluster_AtLeast10.xlsx") -> Annotation_groups
Annotation_groups %>% group_by(Cluster, High_taxonomy) %>% summarise(N = n()) -> Composition_cluster
New_cluster = tibble()
Colors_do = c25[2: length(c25)] ; Colors_do[6] = c25[11]
annoCol<-list(High_taxonomy=c(Aves=Colors_do[1], Bacteria = Colors_do[2], Fungi=Colors_do[3], 
                         Invertebrate = Colors_do[4], Mammal=Colors_do[5], Phage=Colors_do[6],
                         Plant = Colors_do[7], Virus = Colors_do[8]))
tibble_colors = tibble(High_taxonomy = names(annoCol$High_taxonomy), Value = annoCol$High_taxonomy )
#Make fractions of each category per cluster
for (C in unique(Composition_cluster$Cluster)){
        Composition_cluster %>% filter(Cluster == C) -> subset_composition
        sum(subset_composition$N) -> Total 
        subset_composition %>% mutate(Fraction = N/Total ) -> subset_composition
        rbind(New_cluster, subset_composition) -> New_cluster
        
}
#Plot
ggplot(New_cluster, aes(x="", y=Fraction, fill=High_taxonomy)) + geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) + theme_void() + facet_wrap(~Cluster) + scale_fill_manual(values=Colors_do) -> Piecharts
ggsave("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Correlation_vs_Similarity/Piecharts_moduleComposition.pdf")

##################################
##Network visualization###########
##################################

Data_LLD %>% select(Groups_LLD3$Probe) -> Module_info
Network_labels = tibble(Name = names(Net$colors), Cluster = Net$colors,  Color =labels2colors(Net$colors) )
Network_labels %>% filter(Name %in% colnames(Module_info)) -> Network_labels

dev.off()

Make_graph = function(Module_info, Network_labels, Annotation_groups,tibble_colors, Scale_by_prev = T, Title="Co-occurrence modules"){
  set.seed(8794)
  g <- graph.adjacency( as.matrix(as.dist(cor(Module_info, method="pearson"))),
        mode="undirected", weighted=TRUE,diag=FALSE)
  #https://stdworkflow.com/252/the-igraph-package-visualizes-gene-regulatory-networks
  g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
  E(g)[which(E(g)$weight<0)]$color <- "darkblue"
  E(g)[which(E(g)$weight>0)]$color <- "darkred"
  E(g)$weight <- abs(E(g)$weight)
  g <- delete_edges(g, E(g)[which(E(g)$weight<0.3)])
  V(g)$name <- Network_labels$Cluster
  V(g)$shape <- "circle" ; V(g)$color <- "skyblue" ; V(g)$vertex.frame.color <- "white"
  V(g)$label.cex = Network_labels$Cluster

  if (Scale_by_prev == T){
    scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
    vSizes <- (scale01(apply(Module_info, 2, mean)) ) * 10 #Scaled by prevalnece
  } else { vSizes = 10 }
  edgeweights <- E(g)$weight * 2.0
  mst <- mst(g, algorithm="prim")
  
  mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
  mst.clustering <- make_clusters(mst, membership=as.vector(Network_labels$Cluster ))

  #Color choosing so that they are consistent with the ones used for piecharts
  left_join(tibble(Probe = colnames(Module_info)) ,  select(Annotation_groups, c(Probe, High_taxonomy))) -> ForColors
  #as.factor(ForColors$High_taxonomy) -> ForColors2
  #Colors_do[1: length(levels(ForColors2))] -> Colors_do
  #tibble( High_taxonomy = levels(ForColors2), Color = Colors_do ) -> Colors
  left_join(select(Annotation_groups, c(Probe, High_taxonomy)), tibble_colors  ) -> Colors
  left_join(ForColors , Colors) -> Colors
  V(mst)$color <- Colors$Value
  #####
  plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE, vertex.size=vSizes,
     vertex.label.dist=0, vertex.label.color="white", asp=FALSE,
     vertex.label.cex=0.4, edge.width=edgeweights, edge.arrow.mode=0,
    main=Title, vertex.color=V(mst)$color)
}
Make_graph(Module_info, Network_labels, Annotation_groups)

##############################################
###Generate logos and conservation plots######
##############################################
Tree_with_Alignment = function(Cluster, Module_info, Labels_cluster,Annotation_groups, Distance_to_use = "Similarity"){
  #Plot tree with sequence similarity
  #Color scheme from MEME
  Colors_aa = c("dark blue","dark blue","dark blue","dark blue","dark blue","dark blue","dark blue","dark blue",
                "#C8A2C8","#C8A2C8",
                "orange","pink", "yellow", "light blue",
                "red", "red",
                "dark green", "dark green", "dark green", "dark green", "grey")
  names(Colors_aa) = tolower(c("A", "C", "F","I","L","M","V","W",
                               "D","E",
                               "G", "H", "P","Y",
                                "K", "R", 
                                "N","Q", "S","T", "-"))
  #Color scheme from internet
  #Colors_aa = c("light green","light green","red","red","green","blue","blue","blue","pink","#C8A2C8","#C8A2C8","blue","#C8A2C8","dark green", "dark green", "dark blue", "dark green","dark green","orange","orange",NA)
  #names(Colors_aa) = tolower(c("G", "A", "S","T","C","V","I","L","P","F","Y", "M", "W", "N", "Q", "H","D", "E","K", "R", "-"))
                
  if (Distance_to_use == "Similarity"){
    Distance = paste0("/Users/sergio/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Clusters/Distance_",as.character(Cluster),".mat")
    read.table(Distance, row.names=1, skip=1) -> Distance
    lapply(rownames(Distance), function(x){ str_split(x,":")[[1]][1] } ) %>% unlist() -> Names
    rownames(Distance) = Names
    colnames(Distance) = rownames(Distance)
  }else{
    select(Module_info, Labels_cluster$Name) %>% cor() -> Distance
    abs(Distance-1) %>% as.dist() %>% hclust(method = "average") -> Tree
  }
  Distance %>% as.dist() %>% hclust(method = "average") -> Tree
  ape::as.phylo(Tree) -> Tree
  filter(Annotation_groups, Probe %in% Labels_cluster$Name) -> Anno
  ggtree(Tree)  %<+% Anno + geom_tippoint(aes(color = High_taxonomy), size=3) + geom_tiplab(size=2,offset=0.05) + xlim(-.1, 2) + scale_color_manual(values=Cols_plot$Value) + theme(legend.position="none") -> p
  msaplot(p, paste0("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Alignments/",as.character(Cluster),".fasta"), offset=0.55, width=3, color = Colors_aa) + theme(legend.position="none")  -> MSA_plot
  H = 74 #74
  ggsave(filename = paste0("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Alignments/Tree_and_MSA_",as.character(Cluster),".pdf") , MSA_plot,width =120 ,height =H ,units = "mm")
  
}
Logo_from_alignment = function(Cluster){
  library(Biostrings)
  library(ggseqlogo)
  #Colors_aa = c("light green","light green","red","red","green","blue","blue","blue","pink","#C8A2C8","#C8A2C8","blue","#C8A2C8","dark green", "dark green", "dark blue", "dark green","dark green","orange","orange",NA)
  #names(Colors_aa) = c("G", "A", "S","T","C","V","I","L","P","F","Y", "M", "W", "N", "Q", "H","D", "E","K", "R", "-")
  Colors_aa = c("dark blue","dark blue","dark blue","dark blue","dark blue","dark blue","dark blue","dark blue",
                "#C8A2C8","#C8A2C8",
                "orange","pink", "yellow", "light blue",
                "red", "red",
                "dark green", "dark green", "dark green", "dark green", "grey")
  names(Colors_aa) = c("A", "C", "F","I","L","M","V","W",
                               "D","E",
                               "G", "H", "P","Y",
                               "K", "R", 
                               "N","Q", "S","T", "-")
  
  FASTA = paste0("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Alignments/",as.character(Cluster),".fasta")
  s = readAAStringSet(FASTA, format = "fasta")
  
  #Compute and save Inofrmation content plot
  consensusMatrix(s)  -> Count_s
  Count_s %>% apply(2, function(Ex){ 
          #Ex + 1 -> Ex # We add a pseudocount
          sum(Ex) -> Denominator_freq
          Fex = Ex / Denominator_freq #Probability of each AA
          #Fex[2:length(Fex)] -> Fex #We do not take the gap
          -1 * sum(Fex[Fex != 0] * log2(Fex[Fex != 0])) -> Entropy
          Information = log2(21) - Entropy
          return(Information)
  }) -> Information_per_aa
  Count_s %>% apply(2, function(Ex){ 1 - (Ex[1]/sum(Ex)) }) -> Prevalence
  Count_s %>% apply(2, function(Ex){rownames(Count_s)[ which.max(Ex)] } ) -> Majority_voting
  tibble(Seq = seq(length(Information_per_aa)), H =Information_per_aa, Prevalence= as.numeric(Prevalence), Consensus=as.factor(Majority_voting) ) %>% ggplot(aes(x=Seq, y = H, alpha=Prevalence)) + geom_bar(stat="identity", aes(fill=Consensus)) + theme_bw() +
  geom_hline(yintercept = log2(21),linetype = 'dotted', alpha =0.5)  + theme(panel.grid = element_blank(),
          axis.title = element_blank(), axis.text.x = element_blank()) + scale_fill_manual(values = Colors_aa ) + theme(legend.position="none")  -> Information_plot
  ggsave(filename = paste0("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Alignments/Information_content_",as.character(Cluster),".pdf") , Information_plot, width = 207.151, height = 16.826, units = "mm" )
  
  #Generate sequence logo
  
  #Clip extremes with many gaps. Many Gaps >= 3/4 total sequences with gaps
  consensusMatrix(s)[1,] >= 3*length(s)/4 -> Many_gaps
  cumsum(Many_gaps) -> Many_gaps_cum
  cumsum(rev(Many_gaps)) -> Many_gaps_cum_rev
  sapply(seq(Many_gaps_cum)[2:length(Many_gaps_cum)], function(x){ Many_gaps_cum[x-1] ==   Many_gaps_cum[x]  }  ) -> Identical
  sapply(seq(Many_gaps_cum_rev)[2:length(Many_gaps_cum_rev)], function(x){ Many_gaps_cum_rev[x-1] ==   Many_gaps_cum_rev[x]  }  ) -> Identical_rev
  
  match(T, Identical) -> Clip_position_b
  match(T, Identical_rev) -> Clip_position_e

  subseq(s, Clip_position_b + 1,  length(Many_gaps) - (Clip_position_e - 1)  ) -> s_subset
  #plot
  ggplot()  + geom_logo( as.character(s_subset)) + theme_logo() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())   -> Plot
  ggsave(filename = paste0("~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Alignments/Logo_",as.character(Cluster),".pdf") , Plot)
}
  
#Make plots for some specific clusters
for (Cluster_i in unique(Network_labels$Cluster) ){
  Network_labels %>% filter(Cluster == Cluster_i) -> Labels_cluster
  New_cluster %>% mutate(High_taxonomy = factor(High_taxonomy)) %>% filter(Cluster==Cluster_i) -> To_plot
  Cols_plot = filter(tibble_colors, High_taxonomy %in% To_plot$High_taxonomy)
  #Graph  
  Make_graph(select(Module_info, Labels_cluster$Name), Labels_cluster, filter(Annotation_groups, Probe %in% Labels_cluster$Name), Scale_by_prev = F, Title=NA, tibble_colors)
  #Pie chart
  To_plot %>% ggplot(aes(x="", y=Fraction, fill=High_taxonomy)) + geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) + theme_void() +  theme(legend.position="none")  + scale_fill_manual(values=Cols_plot$Value) %>% print()
  #Tree and alignment plot
  Tree_with_Alignment(Cluster_i, Module_info, Labels_cluster,Annotation_groups, Distance_to_use = "Similarity")
  #Make Logo and information content plots
  Logo_from_alignment(Cluster_i)
}

#######################################
####Analisis all distance matrices#####
#######################################
#For this distance matrices between peptides in each module should be available

Compute_heatmap = function(Distance, Name, Subset_info, correlations){
        Out =   "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Correlation_vs_Similarity/" #Output path
        breaksList = seq(0, 1, by = 0.1)
        read.table(Distance, row.names=1, skip=1) -> Distance
        lapply(rownames(Distance), function(x){ str_split(x,":")[[1]][1] } ) %>% unlist() -> Names
        rownames(Distance) = Names
        colnames(Distance) = rownames(Distance)
        Annotation_groups %>% select(Probe, High_taxonomy) %>% as.data.frame() %>% column_to_rownames("Probe") -> Annotation
        png(file= paste(Out, Name, "_similarity.png" ) )
        pheatmap::pheatmap(Distance, color=viridis::viridis(10), breaks = breaksList, main = Name, annotation_row = Annotation, annotation_colors =  annoCol )
        dev.off()
        
        Distance2 = Distance[upper.tri(Distance)]
        Similarity_score = mean(Distance2)
        Similarity_deviation = sd(Distance2)
        
        Data_LLD %>% select(colnames(Distance)) %>% cor() -> Corr_matrix
        test = hclust(dist(Corr_matrix))
        ORDER = rownames(Corr_matrix)[test$order]
        as.data.frame(Corr_matrix)[ORDER,ORDER] -> Corr_matrix ; Distance[ORDER, ORDER] -> Distance

        New_matrix = 1- Distance
        New_matrix[lower.tri(New_matrix)] <- Corr_matrix[lower.tri(Corr_matrix)]
        #New_matrix %>% drop_na() -> New_matrix
        pdf(file= paste(Out, Name, "_simVScorr.pdf" ) )
        corrplot::corrplot(as.matrix(New_matrix), diag=FALSE, tl.col="black")
        dev.off()
        
        New_matrix = as.matrix(New_matrix)
        diag(New_matrix) = NA
        png(file= paste(Out, Name, "_simVScorr_AnnoHeatmap.png" ) )
        pheatmap::pheatmap(New_matrix, color=RColorBrewer::brewer.pal(name = "Blues", n = 9), breaks = breaksList, annotation_row = Annotation, annotation_colors =  annoCol, cluster_rows = F, cluster_cols  = F, annotation_legend = F)
        dev.off()
        #Get P-value
        vegan::mantel( as.dist(Corr_matrix) , as.dist(1-Distance), method = "spearman", permutations = 2000, na.rm = TRUE) -> Correlation_result
        r = Correlation_result$statistic
        p = Correlation_result$signif
        Result = tibble(Mean_similarity = Similarity_score, Sd_similarity=Similarity_deviation, r_mantel=r, P_mantel=p)
        return(Result)
        
}

files <- list.files(path="/Users/sergio/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Clusters", pattern="Distance_[0-9]+", full.names=TRUE, recursive=FALSE)
read_tsv("/Users/sergio/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Cluster_AtLeast10.tsv") -> Groups #Cluster belonging

Similarity_cluster = tibble()
for (Fi in files){
        basename(tools::file_path_sans_ext(Fi) ) -> N
        Annotation_groups %>% filter(Cluster == as.numeric(str_split(N, "_")[[1]][2]) ) ->Subset_info
        Compute_heatmap(Distance =Fi, Name = N, Subset_info, correlations) -> D
        print(paste(N, D))
        rbind(Similarity_cluster, mutate(D, Cluster = str_split(N, "_")[[1]][2]) ) -> Similarity_cluster
}

Groups %>% filter(Cluster == 1) %>% select(Probe,  Taxa, Description) %>% print(n=99)
Similarity_cluster %>% filter(P_mantel>0.05)
write_tsv(Similarity_cluster,path = "Similarity_scores.tsv")



###############################################
#################Analysis Binary Network#######
###############################################
library(IsingFit)
#IsingFit_parl is the IsingFit function, but including 6 core parallelization from the foreach package
IsingFit_parl <-function(x, family='binomial', AND = TRUE, gamma = 0.25, plot = TRUE, lowerbound.lambda = NA,...){
  t0 <- Sys.time()
  xx <- x
  
  ## Check to prevent error of lognet() in package glmnet
  checklognet <- function(y){
    res <- c() # 0: too little variance, 1: good to go
    y=as.factor(y)
    ntab=table(y)
    minclass=min(ntab)
    if(minclass<=1) res=0 else res=1
    return(res)
  }
  NodesToAnalyze <- apply(x,2,checklognet) !=0
  names(NodesToAnalyze) <- colnames(x)
  if (!any(NodesToAnalyze)) stop("No variance in dataset")
  if (any(!NodesToAnalyze))
  { warning(paste("Nodes with too little variance (not allowed):",paste(colnames(x)[!NodesToAnalyze],collapse = ", "))) }
  ##
  
  x <- as.matrix(x)
  allthemeans <- colMeans(x)
  x <- x[,NodesToAnalyze,drop=FALSE]
  nvar <- ncol(x)
  p <- nvar - 1
  
  ##Parallel computation of lasso. GLMNET is used. Default alpha=1 (lasso), it tries several different lambda parameters.
  nlambdas <- rep(0,nvar)
  intercepts <- betas <- lambdas <- list(vector,nvar)
  doParallel::registerDoParallel(6)
  print("Lasso computation")
  foreach( i = seq(nvar) ) %dopar% {
    a <- glmnet::glmnet(x = x[,-i], y = x[,i], family = family)
    return( list(i, a$a0, a$beta, a$lambda) )
  } -> Output
  for (a in Output){
    i = a[[1]]
    intercepts[[i]] <- a[[2]] ; betas[[i]] <- a[[3]] ;  lambdas[[i]] <- a[[4]] ; nlambdas[i] <- length( a[[4]]) 
  }
  print("Lasso completed")
  
  P <- logl <- sumlogl <- J <- matrix(0, max(nlambdas), nvar)
  #Count how many non 0 parameters were picked by lasso
  for (i in 1:nvar){
    J[1:ncol(betas[[i]]),i] <- colSums(as.matrix(betas[[i]]!=0))
  }
  logl_M <- P_M <- array(0, dim=c(nrow(x),max(nlambdas), nvar) )
  N <- nrow(x)
  print("P-matrix computation")
  #For each Node (feature) compute per each value of lambda a P_matrix and its log (likelihood?)
  foreach( i = 1: nvar ) %dopar% {
    betas.ii <- as.matrix( betas[[i]] )
    int.ii <- intercepts[[i]]
    y <- matrix( 0 , nrow=N , ncol= ncol(betas.ii) ) 
    xi <- x[,-i]
    NB <- nrow( betas.ii) # number of rows in beta
    #For each peptide, get the prediction of y for all lambda
    for (bb in 1:NB){   # bb <- 1
      y <- y + betas.ii[rep(bb,N),] * xi[,bb]
    }
    y <- matrix( int.ii , nrow=N , ncol=ncol(y) , byrow=TRUE ) + y
    # number of NAs
    n_NA <- max(nlambdas)-ncol(y)
    if (n_NA > 0 ){ 
      for ( vv in 1:n_NA){ 
        y <- cbind( y , NA ) 
      } 
    }
    # calculate P matrix
    P = exp(y*x[,i])/(1+exp(y))
    M = log(P) 
    return(list(i, P, M))
  } -> P_matrix_info
  for (a in P_matrix_info){
    i = a[[1]]
    P_M[,,i] <- a[[2]]
    logl_M[,,i] <- a[[3]]
  }
  print("P-matrix computation finished")
  
  logl_Msum <- colSums( logl_M , 1, na.rm=FALSE )
  sumlogl <- logl_Msum 
  sumlogl[sumlogl==0]=NA
  penalty <- J * log(nrow(x)) + 2 * gamma * J * log(p)
  EBIC <- -2 * sumlogl + penalty
  
  lambda.mat <- matrix(NA,nrow(EBIC),ncol(EBIC))
  for (i in 1:nvar){
    lambda.mat[,i] <- c(lambdas[[i]],rep(NA,nrow(EBIC)-length(lambdas[[i]])))
  }
  
  if(!is.na(lowerbound.lambda)){
    EBIC <- EBIC/(lambda.mat>=lowerbound.lambda)*1
  }
  
  lambda.opt <- apply(EBIC,2,which.min)
  lambda.val <- rep(NA,nvar)
  thresholds <- 0
  for(i in 1:length(lambda.opt)){
    lambda.val[i] <- lambda.mat[lambda.opt[i],i]
    thresholds[i] <- intercepts[[i]][lambda.opt[i]]
  }
  weights.opt <- matrix(,nvar,nvar)
  for (i in 1:nvar){
    weights.opt[i,-i] <- betas[[i]][,lambda.opt[i]]
  }
  asymm.weights <- weights.opt
  diag(asymm.weights)=0
  if (AND==TRUE) {
    adj <- weights.opt
    adj <- (adj!=0)*1
    EN.weights <- adj * t(adj)
    EN.weights <- EN.weights * weights.opt
    meanweights.opt <- (EN.weights+t(EN.weights))/2
    diag(meanweights.opt) <- 0 
  } else {
    meanweights.opt <- (weights.opt+t(weights.opt))/2
    diag(meanweights.opt) <- 0
  }
  graphNew <- matrix(0,length(NodesToAnalyze),length(NodesToAnalyze))
  graphNew[NodesToAnalyze,NodesToAnalyze] <- meanweights.opt
  colnames(graphNew) <- rownames(graphNew) <- colnames(xx)
  threshNew <- ifelse(allthemeans > 0.5, -Inf, Inf)
  threshNew[NodesToAnalyze] <- thresholds
  if (plot==TRUE) notplot=FALSE else notplot=TRUE
  Res <- list(weiadj = graphNew, thresholds = threshNew, q = q, gamma = gamma, 
              AND = AND, time = Sys.time() - t0, asymm.weights = asymm.weights,
              lambda.values = lambda.val)
  class(Res) <- "IsingFit"
  return(Res)
}
Get_clusters = function(Distance, Data, minCorrMerge=0.5, cutoff=10){
  geneTree = hclust(as.dist(Distance), method = "average")
  
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = Distance,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = cutoff  )

  #Some summaries, number of samples per module and give colors to those samples
  table(dynamicMods)
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  # Eigengenes and merging of modules which eigengene is closer than a  certain (correlation) threshold with the others
  #Eigengenes to compute similarity between modules
  MEList = moduleEigengenes(Data, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  # Call an automatic merging function
  merge = mergeCloseModules(Data, dynamicColors, cutHeight = minCorrMerge, verbose = 3)
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  #Comparison between modules before and after merging
  #Some renaming
  moduleColors = mergedColors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs
  
  return(list(moduleColors, MEs))
  
}
Compare_clusters = function(Eigen_reference, Eigen){
  Eigen %>% as_tibble() %>% mutate(ID = filter(Data,grepl("32_", ID))$ID, .before=1) %>% select(-MEgrey) -> EigenGenesBinary

  full_join(EigenGenesBinary, Eigen_reference, by="ID", suffix=c("_binary", "_corr")) %>% select(-ID)  %>% cor() -> Cor_all #%>% heatmap()
  hclust(as.dist( 1- Cor_all), method="average") -> CLUSTER
  cutree(CLUSTER, h = 0.05) -> CLUSTER
  Matches = tibble()
  for (C in CLUSTER){
    CLUSTER[CLUSTER==C] -> Pair
    if (! length(Pair) == 2 ){ next }
    names(Pair)[grepl("ME[[:digit:]]+$", names(Pair))] -> Corr_cluster
    Number_change = names(Pair)[! grepl("ME[[:digit:]]+$", names(Pair)) ]
    Matches = rbind(Matches, tibble(New = names(Pair)[1] , Reference = names(Pair)[2] ) )
    EigenGenesBinary[ Corr_cluster ] = EigenGenesBinary[,Number_change]
  }
  return(list(EigenGenesBinary, Matches))
}
set.seed(9897)
#This takes some time to work, so load the RDS file with the output
IsingFit_parl(Data_LLD, AND=TRUE) -> res_bin #this is a weighted network
saveRDS(res_bin, "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/BinaryNetwork/Binary_network.rds")

##Normalization of the adjacency matrix so that it goes from 0-1
res_bin$weiadj -> Adj # 0 - inf
Adj_norm = Adj/max(Adj)
Adj_norm[Adj_norm<0] = 0

#Create TOM from the adjacency matrix?
TOM = TOMsimilarity(Adj_norm) #topological overlap of two nodes for unweighted networks ; not sure this would be the right call
#TOM are meant to identify better clusters in RNAseq data. We can anyway use a dissimilarity measure like 1-Corr (if no negative values I assume). We can also transform the negative values

#Define distance matrices
Dist_tom = 1-TOM
Dist = 1 - Adj_norm

Get_clusters(Dist, Data_LLD) -> Results1
Get_clusters(Dist_tom, Data_LLD) -> Results2

read_excel("~/Resilio Sync/Antibodies_WIS (1)/Results/Supplemenatry_Tables/SupplementaryTable1.xlsx", sheet = 2) -> AB_annotation

#Do any of the cluster belongings resamble the correlation based?
#8 of 12 have homologous
#method returns the dataframe of the eigengens in Eigen + The homologous egengenes from reference
Compare_clusters(Eigen = Results1[[2]], Eigen_reference = EigenGenes) -> NewClusters
NewClusters[[2]] %>% distinct() -> Matches
NewClusters[[1]] -> NewClusters
tibble( Peptide = colnames(Data_LLD), Cluster = paste0("ME", Results1[[1]]) ) %>% filter(! Cluster == "MEgrey") -> ClustersPeptides
colnames(Matches) = c("Cluster", "Homologous_WGCNA")
left_join(ClustersPeptides, Matches) %>% arrange(Cluster) %>% left_join( select(AB_annotation, c(Peptide, Taxa, Description))) -> ClustersBinary
write_tsv(ClustersBinary, "Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/BinaryNetwork/ClustersBinary.tsv")
#9 of  have hologous 21
Compare_clusters(Eigen = Results2[[2]], Eigen_reference = EigenGenes)

#Compare both methods
Results1[[2]] %>% as_tibble() %>% mutate(ID = filter(Data,grepl("32_", ID))$ID, .before=1) %>% select(-MEgrey) -> EigenGenesBinary
colnames(EigenGenesBinary) = c("ID", paste0("ME", seq(1:(dim(EigenGenesBinary)[2]-1) ) )  )
#11 homologous refs out of 23 found in TOM out of 12 possible homologous in the non-TOM: missing 8
Compare_clusters(Eigen = Results2[[2]], Eigen_reference = EigenGenesBinary)

#Check which clusters were reproduced and which are new
NewClusters %>% select(-ID) %>% cor() -> ClustersCor
ClustersCor %>% as.data.frame() %>% rownames_to_column("Cluster") %>% as_tibble() %>% filter( grepl( "ME[[:digit:]]+$",Cluster )) %>% select(- rownames(ClustersCor)[grepl( "ME[[:digit:]]+$",rownames(ClustersCor))]) -> ClustersCor
ClustersCor %>% apply(1, function(x){  names(x)[as.numeric(x) == 1] } ) -> Reproduced
Reproduced[! is.na(Reproduced)] -> Reproduced
#New 4 clusters
NewClusters %>% select(one_of(colnames(ClustersCor))) %>% select(- Reproduced ) -> NewClusters
tibble( Peptide = colnames(Data_LLD), Cluster = paste0("ME", Results1[[1]]) ) %>% filter( Cluster %in% colnames(NewClusters)) -> NewClustersPeptides
#Add Annotation
left_join(NewClustersPeptides, select(AB_annotation, c(Peptide, Taxa, Description))) -> NewClustersPeptides
#CLusters with several identical bacteria: 
#Pink (4 Bacteroides dorei and 3 Bacteroides), Red (root 3, B. fragilis 2, Lachnospiraceae 2 ), Brown (Ligonella 2), Green (Clostridiales 2)
NewClustersPeptides %>% group_by(Cluster, Taxa) %>% summarise(N = n()) %>% arrange(desc(N))
#Brown: 2 allergens + phages + Bugs -> No high identity
#Green: Mosqito + Bugs  --> No high identity
#Pink: Bacteroides: TOnB receptor, STN domain-containing (foundin TOnB) --> Homologous
#Red: Bacteria (D-galactose/ D-glucose-binding protein ,  Peripla_BP_4 domain-containing protein, Peptidase_M48 domain-containing protein)  + EBV --> Common motif in bacterial
write_tsv(NewClustersPeptides,"Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/BinaryNetwork/NewClustersBinary.tsv")

###Seqeuence Similarity of new modules 
files <- list.files(path="/Users/sergio/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Clusters", pattern="Distance_[^0-9]+", full.names=TRUE, recursive=FALSE)
Similarity_cluster = tibble()
for (Fi in files){
  basename(tools::file_path_sans_ext(Fi) ) -> N
  NewClustersPeptides %>% filter(Cluster ==str_split(N, "_")[[1]][2])  ->Subset_info
  Compute_heatmap(Fi, N, Subset_info, correlations) -> D
  print(paste(N, D))
  rbind(Similarity_cluster, mutate(D, Cluster = str_split(N, "_")[[1]][2]) ) -> Similarity_cluster
}

Groups %>% filter(Cluster == 1) %>% select(Probe,  Taxa, Description) %>% print(n=99)
Similarity_cluster %>% filter(P_mantel>0.05)
write_tsv(Similarity_cluster,path = "Similarity_scores.tsv")


#New figure 2
Data %>% select( c("ID", Groups_LLD3$Probe) ) %>% select( -ID ) %>% cor() -> Corr_Matrix_Modules
left_join( Network_labels, select(mutate(Annotation_groups, Name = Probe), c(Name, Cluster, High_taxonomy) )  ) -> Network_labels
Network_labels %>% mutate(Cluster = as.factor(Cluster)) %>% as.data.frame() %>% column_to_rownames("Name") %>% select(-Color) -> Annotation_heatmap
#Color
annotation_colors = unique(Network_labels$Color ) 
names(annotation_colors) = unique(Network_labels$Cluster)

annotation_colors2 = Colors_do[1: length(unique(Network_labels$High_taxonomy)) ]
names(annotation_colors2) = sort(unique(Network_labels$High_taxonomy))


list(  Cluster = annotation_colors, High_taxonomy= annotation_colors2  ) -> annotation_colors
pheatmap::pheatmap(Corr_Matrix_Modules, annotation_row = Annotation_heatmap, annotation_colors =  annotation_colors, labels_row = rep("", nrow(Corr_Matrix_Modules) ) , labels_col = rep("", nrow(Corr_Matrix_Modules) ) )
Annotation_groups %>% group_by(Cluster, High_taxonomy) %>% summarise(N = n()) -> Composition_cluster


