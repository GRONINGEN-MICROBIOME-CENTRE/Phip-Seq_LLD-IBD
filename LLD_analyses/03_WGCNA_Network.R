#S. Andreu-Sanchez
#Uses WGCNA in a presence/absence profile to identify modules of highly correlated peptides
#Also: Module correlation. Module composition plotting. Network plotting.
#Sequence similarity clustring

library(WGCNA) # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/OverviewWGCNA.pdf
library(tidyverse)
library(readxl)
library(igraph)
library(vegan)

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
        Do_analysis(Data_all, 7) -> results_all
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
                      sample(seq(dim(Data_LLD)[1]),a*dim(Data_LLD)[1],replace=F)  -> Boost
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
readxl::read_xlsx("Cluster_AtLeast10.xlsx") -> Annotation_groups
Annotation_groups %>% group_by(Cluster, High_taxonomy) %>% summarise(N = n()) -> Composition_cluster
New_cluster = tibble()
Colors_do = c25[2: length(c25)] ; Colors_do[6] = c25[11]
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
ggsave("Piecharts_moduleComposition.pdf")

##################################
##Network visualization###########
##################################

Data_LLD %>% select(Groups_LLD3$Probe) -> Module_info
Network_labels = tibble(Name = names(Net$colors), Cluster = Net$colors,  Color =labels2colors(Net$colors) )
Network_labels %>% filter(Name %in% colnames(Module_info)) -> Network_labels

dev.off()


g <- graph.adjacency(
        as.matrix(as.dist(cor(Module_info, method="pearson"))),
        mode="undirected",
        weighted=TRUE,
        diag=FALSE
)
#https://stdworkflow.com/252/the-igraph-package-visualizes-gene-regulatory-networks
#g = graph.adjacency(
#        as.matrix( vegdist(Data_LLD, method = "jaccard") ),
#        mode="undirected",
#        weighted=TRUE,
#        diag=FALSE
#)

g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
E(g)[which(E(g)$weight<0)]$color <- "darkblue"
E(g)[which(E(g)$weight>0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight<0.3)])
V(g)$name <- Network_labels$Cluster
V(g)$shape <- "circle" ; V(g)$color <- "skyblue" ; V(g)$vertex.frame.color <- "white"
V(g)$label.cex = Network_labels$Cluster

scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(Module_info, 1, mean)) ) * 10 #Scaled by prevalnece
edgeweights <- E(g)$weight * 2.0
mst <- mst(g, algorithm="prim")

mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=as.vector(Network_labels$Cluster ))

#Color choosing so that they are consistent with the ones used for piecharts
left_join(tibble(Probe = colnames(Module_info)) ,  select(Annotation_groups, c(Probe, High_taxonomy))) -> ForColors
as.factor(ForColors$High_taxonomy) -> ForColors2
Colors_do[1: length(levels(ForColors2))] -> Colors_do
tibble( High_taxonomy = levels(ForColors2), Color = Colors_do ) -> Colors
left_join(ForColors , Colors) -> Colors
V(mst)$color <- Colors$Color
#####
plot(mst,
     layout=layout.fruchterman.reingold,
     edge.curved=TRUE,
     vertex.size=vSizes,
     vertex.label.dist=0,
     vertex.label.color="white",
     asp=FALSE,
     vertex.label.cex=0.4,
     edge.width=edgeweights,
     edge.arrow.mode=0,
     main="Co-occurrence modules", vertex.color=V(mst)$color)




#######################################
####Analisis all distance matrices#####
#######################################
#For this distance matrices between peptides in each module should be available

Compute_heatmap = function(Distance, Name, Subset_info, correlations){
        Out =   "" #Output path
        breaksList = seq(0, 1, by = 0.1)
        read.table(Distance, row.names=1, skip=1) -> Distance
        lapply(rownames(Distance), function(x){ str_split(x,":")[[1]][1] } ) %>% unlist() -> Names
        rownames(Distance) = Names
        colnames(Distance) = rownames(Distance)
        png(file= paste(Out, Name, "_similarity.png" ) )
        pheatmap::pheatmap(Distance, color=viridis::viridis(10), breaks = breaksList, main = Name)
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
        
        #Get P-value
        vegan::mantel( as.dist(Corr_matrix) , as.dist(1-Distance), method = "spearman", permutations = 2000, na.rm = TRUE) -> Correlation_result
        r = Correlation_result$statistic
        p = Correlation_result$signif
        Result = tibble(Mean_similarity = Similarity_score, Sd_similarity=Similarity_deviation, r_mantel=r, P_mantel=p)
        return(Result)
        
}

files <- list.files(path="~/Desktop/Clusters", pattern="Distance_*", full.names=TRUE, recursive=FALSE)
read_tsv("Cluster_AtLeast10.tsv") -> Groups #Cluster belonging

Similarity_cluster = tibble()
for (Fi in files){
        basename(tools::file_path_sans_ext(Fi) ) -> N
        Groups %>% filter(Cluster == as.numeric(str_split(N, "_")[[1]][2]) ) ->Subset_info
        Compute_heatmap(Fi, N, Subset_info, correlations) -> D
        print(paste(N, D))
        rbind(Similarity_cluster, mutate(D, Cluster = str_split(N, "_")[[1]][2]) ) -> Similarity_cluster
}

Groups %>% filter(Cluster == 1) %>% select(Probe,  Taxa, Description) %>% print(n=99)
Similarity_cluster %>% filter(P_mantel>0.05)
write_tsv(Similarity_cluster,path = "Similarity_scores.tsv")





