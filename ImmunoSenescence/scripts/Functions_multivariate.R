library(tidyverse)
library(caret)
library(PMA)
library(patchwork)
library(ppcor)
library(fgsea)
library(clusterProfiler)



  Leave_one_out_CV = function(X, Y, penaltyX, penaltyY){
    #Leave-one-out cross-validation. different combinations of i/j are tried. A CCA model is fit with all data -1. The remaining is then used to match X and Y when multiplied by the coefficients estimated by the training
    scoreXcv <- c()
    scoreYcv <- c()
    corr_CRC <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
    num_samples <- nrow(X)
    
    
    for( i in 1:length(penaltyX)){
      for(j in 1:length(penaltyY)){
        for(k in 1:num_samples){
          print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
          #compute weights with sample k held out:
          res <- CCA(X[-k,] %>% apply(2, as.numeric) ,Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
          ## Compute scores for k'th sample for first pair of canonical variables
          ## Take weight of features (res$u and res$v) computed using all except 
          ## the kth sample and multiple it by values for the kth sample in the 
          ## feature matrix X and Y. 
          ## %*% implies matrix multiplication. 
          scoreXcv[k] <- X[k,]%>% apply(2, as.numeric) %*%res$u ## single value
          scoreYcv[k] <- Y[k,] %>% apply(2, as.numeric) %*%res$v ## single value
        }
        ## correlation between scores for X and Y for all held out samples.
        corr_CRC[i,j] = cor(scoreXcv,scoreYcv) 
      }
    }
    return(corr_CRC)
  }
  KFold_CV = function(X, Y, penaltyX, penaltyY, K=5){
    #Leave-one-out cross-validation. different combinations of i/j are tried. A CCA model is fit with all data -1. The remaining is then used to match X and Y when multiplied by the coefficients estimated by the training
    scorecv <- c()
    corr_CRC <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
    num_samples <- nrow(X)
    Folds = caret::createFolds( as_tibble(t(Y)), k = K, list = TRUE, returnTrain = FALSE)
    
    for( i in 1:length(penaltyX)){
      for(j in 1:length(penaltyY)){
        for(n_k in 1:K){
          k = Folds[[n_k]]
          print(paste0("Index: i = ",i,", j =", j," k = ",n_k)); flush.console()
          #compute weights with sample k held out:
	  X_k = X[-k,] %>% apply(2, as.numeric)
	  Y_k = Y[-k,]
          Y[-k,] %>% apply(2, sd) -> S
          Remove = names(S[S==0])
	  Y_k %>% select(-Remove) -> Y_k
          res <- CCA(X_k ,Y_k, penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
          ## Compute scores for k'th sample for first pair of canonical variables
          ## Take weight of features (res$u and res$v) computed using all except 
          ## the kth sample and multiple it by values for the kth sample in the 
          ## feature matrix X and Y. 
          ## %*% implies matrix multiplication. 
          ScoreX <- X_k%>% apply(2, as.numeric) %*%res$u ## single value
          ScoreY <- Y_k %>% apply(2, as.numeric) %*%res$v ## single value
          scorecv[n_k] = cor(ScoreX, ScoreY)
        }
        ## correlation between scores for X and Y for all held out samples.
        corr_CRC[i,j] = mean(scorecv[!is.na(scorecv)]) 
      }
    }
    return(corr_CRC)
  }
  Held_out_gridSearch = function(X, Y, penaltyX, penaltyY){
    #Leave-one-out cross-validation. different combinations of i/j are tried. A CCA model is fit with all data -1. The remaining is then used to match X and Y when multiplied by the coefficients estimated by the training
    scoreXcv <- c()
    scoreYcv <- c()
    corr_CRC <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
    Held_out = sample(nrow(X), nrow(X)*0.2 )
    for( i in 1:length(penaltyX)){
      for(j in 1:length(penaltyY)){
        print(paste0("Index: i = ",i,", j =", j)); flush.console()
        #compute weights with sample k held out:
        res <- CCA(X[-Held_out,] %>% apply(2, as.numeric) ,Y[-Held_out,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
        ## Compute scores for k'th sample for first pair of canonical variables
        ## Take weight of features (res$u and res$v) computed using all except 
        ## the kth sample and multiple it by values for the kth sample in the 
        ## feature matrix X and Y. 
        ## %*% implies matrix multiplication. 
        scoreXcv <- X[Held_out,]%>% apply(2, as.numeric) %*%res$u ## single value
        scoreYcv <- Y[Held_out,] %>% apply(2, as.numeric) %*%res$v ## single value
        ## correlation between scores for X and Y for all held out samples.
        corr_CRC[i,j] = cor(scoreXcv,scoreYcv) 
      }
    }
    return(corr_CRC)
  }
  #Fit CCA
  run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
    CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                    penaltyx=penaltyX,penaltyz=penaltyZ,
                    v=vInit) ## standardize=T by default
    
    ## add rownames to output factors
    rownames(CCA.out$u) <- colnames(X)
    rownames(CCA.out$v) <- colnames(Z)
    ## compute contribution of selected features to each of the samples.
    CCA_var_proteomics <- X %*% CCA.out$u ## canonical variance for genes 
    CCA_var_AB <- Z %*% CCA.out$v ## canonical variance for microbes
    
    return(list(CCA.out, CCA_var_proteomics, CCA_var_AB))
    
  }
  #Get average of non-0 components among all components
  get_avg_features <- function(cca_cov, CCA.K){
    num_features <- 0
    for(k in 1:CCA.K){
      num_features <- num_features + length(which(cca_cov[,k]!=0))
    }
    avg_features <- num_features/CCA.K
  }
  
  #Associate covariates to the canonical variables
  Check_confounders = function( Results ){
    #Do the Canonical variables correlate with confounders?
    #Proteomics
    as_tibble(Results[[3]][[2]]) %>% mutate(ID = Results[[4]], .before=1) %>% left_join(. ,select(Covariates, c("Age", "Sex","smk_now", "ID") )) -> Proteomics_components
    #AB
    as_tibble(Results[[3]][[3]]) %>% mutate(ID = Results[[4]], .before=1) %>% left_join(. ,select(Covariates, c("Age", "Sex","smk_now", "ID") )) -> AB_components
    Confounding_effects = tibble()
    for (i in 1:10){
      Formula = as.formula( paste0("V",i," ~ Sex + Age + smk_now" ) )
      lm(Formula, Proteomics_components) %>% summary() -> D1
      lm(Formula, AB_components) %>% summary() -> D2
      
      D1$coefficients %>% as.data.frame() %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Component = i, Data="Proteomics") -> Confounding_effects_prot
      D2$coefficients %>% as.data.frame() %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Component = i, Data="PhIP-Seq") -> Confounding_effects_phip
      rbind(Confounding_effects_prot, Confounding_effects_phip) %>% rbind(Confounding_effects, .) -> Confounding_effects
      
    }
    Confounding_effects %>% arrange(`Pr(>|t|)` ) %>% filter(! Feature == "(Intercept)") %>% filter(Component %in% filter(Results[[2]], P<0.05)$Component  ) %>% print()
    return(Confounding_effects)
  }
  
  #Build network of significant canonical variables
  Make_network = function( Results, DF = Proteomics2, P_t = 0.05, Feature_name="Proteins_CVD", Min_corr=0.2  ){
    for ( N in filter(Results[[2]], P < P_t)$Component ){
      Results[[1]] %>% dplyr::filter(Component == N ) -> Component_network
      left_join(DF, ImmunoMatrix) %>% select(Component_network$Feature) %>% drop_na() %>% cor() -> Correlations_network
      Groups = list(Feature_name = which( rownames(Correlations_network) %in% filter(Component_network, Data=="Protein")$Feature) , Ab_bound_peptide= which( rownames(Correlations_network) %in% filter(Component_network, Data=="Antibody")$Feature ) )
      qgraph(Correlations_network,minimum=Min_corr,cut=0.4,vsize=2,legend=TRUE,borders=FALSE,  groups= Groups)
    }
  }
  
  Plot_loads = function(Results, Component_n= 5, Color = "Taxa", Component_x = "Component_protein", Labels_x_size= 10 ){
    #Color can be either "Taxa" or "DB"
    #Get pepetide annotation
    readxl::read_excel("~/Resilio Sync/Antibodies_WIS (1)/Results/Supplemenatry_Tables/SupplementaryTable1.xlsx", sheet="SupTable1.1") %>% dplyr::rename(Feature = Peptide)  -> Annotation
    
    if (Color == "Taxa"){
      Results[[1]] %>% filter(Component == Component_n) %>% filter(Data=="Antibody") %>% left_join(. ,Annotation) %>% filter(!is.na(Coefficient) ) %>%
        arrange(desc(abs(Coefficient))) %>% ggplot( aes(x=reorder(Feature, Coefficient),y= Coefficient, fill=High_taxonomy ) ) + geom_bar(stat="identity") + coord_flip() + 
        theme(axis.title.y=element_blank(), axis.text.y= element_blank()  ,axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank()) + scale_fill_manual(values = cp ) +  guides(fill = guide_legend(title.position = "left", ncol=1,keywidth=0.5, label.theme= element_text(size=5),title=NULL )) -> Loads_AB
    } else {
      Results[[1]] %>% filter(Component == Component_n) %>% filter(Data=="Antibody") %>% left_join(. ,Annotation) %>%
        arrange(desc(abs(Coefficient))) %>% ggplot( aes(x=reorder(Feature, Coefficient),y= Coefficient, fill=Source_database ) ) + geom_bar(stat="identity") + coord_flip() + 
        theme(axis.title.y=element_blank(), axis.text.y=element_text(size=2)  ,axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank()) + scale_fill_discrete(type = c25 ) +  guides(fill = guide_legend(title.position = "left", ncol=1) )-> Loads_AB
    }   
    
    Results[[1]] %>% filter(Component == Component_n) %>% filter(Data=="Protein") %>%
      arrange(desc(abs(Coefficient))) -> FP
    Stats = Results[[2]] %>% filter(Component == Component_n)
    FP %>% ggplot( aes(x=reorder(Feature, Coefficient),y= Coefficient ) ) +
      geom_segment( aes(x=reorder(Feature, Coefficient), xend=reorder(Feature, Coefficient), y=0, yend=Coefficient), color="grey") + geom_point( color=ifelse(FP$Coefficient>0, "salmon", "steel blue" ), size=1.5) + #geom_bar(stat="identity")
      theme(axis.title.x=element_blank(), axis.text.x=element_text(size=4, angle = 90, vjust = 0.5, hjust=1),axis.ticks.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank() ) + guides(color="none") + ylim( ifelse(min(FP$Coefficient)>0, -0.1,min(FP$Coefficient)) , ifelse( max(FP$Coefficient)<0, 0.1,max(FP$Coefficient) )) + geom_hline(yintercept = 0,linetype = "dashed") +  theme(axis.text.x=element_text(size=Labels_x_size))  -> Loads_protein
    Results[[6]] %>% filter(Component == Component_n) %>% ggplot(aes(x=Component_protein, y=Component_phipseq)) + geom_point() + theme_bw() + geom_smooth(method = "lm") +  annotate("text",  x=Inf, y = Inf, label = paste0("Pearson_rho=",as.character(round(Stats$Cor_coefficient_p,2)) ,"\nP=",as.character(round(Stats$P_p,3))), vjust=1, hjust=1) + xlab(Component_x) -> Correlation_Components
    
    layout <- "
BAAAA
BAAAA
BAAAA
DCCCC
"
    Correlation_Components+ Loads_AB+ Loads_protein +  guide_area() +  plot_layout(design = layout, guides = 'collect') -> Plot
    return(Plot)
    
  }
  
  #Run CCA analysis
  Analysis = function(omics, ImmunoMatrix_y, test_omics ,min_Prevalence_AB = 0.1, penaltyX = seq(0.1,0.4,length=10) , penaltyY=seq(0.2,0.8,length=15), Correct=F ){
    set.seed(7863)
    
    #Prepare data together
    print("Matching datasets")
    left_join(omics, ImmunoMatrix_y, by="ID") %>% drop_na()  -> All_together
    All_together %>% select(colnames(omics)) %>% select(-ID) -> omics_match
    All_together %>% select(colnames(ImmunoMatrix_y)) %>% select(-ID) -> ImmunoMatrix_match
    ImmunoMatrix_match %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalence_ab
    tibble(Ab = names(Prevalence_ab), Prev =  Prevalence_ab) %>% arrange(Prev) -> Prevalence_ab
    ImmunoMatrix_match %>% dplyr::select(- filter(Prevalence_ab, Prev < min_Prevalence_AB)$Ab ) -> Filtered_ImmunoMatrix_match
    
    X <- omics_match
    Y <- Filtered_ImmunoMatrix_match
    print( paste0("Number of Ab: ", dim(Y) ) )
    print( paste0("Number of other layer: ", dim(X) ) )
    
    #Split Train/Test
    #Held_out = sample(nrow(X), nrow(X)*0.2 )
    X_tr = X#[-Held_out,] ; X_te = X[Held_out,]
    X_te = test_omics %>% filter(ID %in% ImmunoMatrix_y$ID) %>% arrange(ID) %>% select(-ID)
    if (Correct == F){
    Y_tr = Y#[-Held_out,] ; Y_te = Y[Held_out,]
    Y_te = ImmunoMatrix_y %>% filter(ID %in% test_omics$ID) %>% arrange(ID) %>% select(colnames(Y_tr))	
    } else {
    Y_tr = Remove_variability(Y %>% mutate(ID = All_together$ID) , Covariates, Covariates_n = c("Age", "Sex")  ) %>% select(-ID)	
    Y_te = test_omics %>% filter(ID %in% ImmunoMatrix_y$ID) %>% arrange(ID) %>% Remove_variability(. , Covariates, Covariates_n = c("Age", "Sex")  ) %>% select(-ID)
    }	
    X_te = test_omics %>% filter(ID %in% ImmunoMatrix_y$ID) %>% arrange(ID) %>% select(-ID)
    Y_te = ImmunoMatrix_y %>% filter(ID %in% test_omics$ID) %>% arrange(ID) %>% select(colnames(Y_tr)) 
    print( paste0("Held out samples: ", dim(X_te)[1]) )
    ################################################
    ## select tuning parameters using grid-search###
    ################################################
    #https://github.com/blekhmanlab/host_gene_microbiome_interactions/blob/main/sparseCCA/grid_search_sparseCCA.R
    
    
    
    #Gid search
    #penaltyX <- seq(0.1,0.4,length=10)
    #penaltyY <- seq(0.2,0.8,length=15)
    #Find best parameters using one of the grid search functions. From fastest to slowest: Held_out_gridSearch > KFold_CV > Leave_one_out_CV
    print("Running grid-search")
    #Res = Held_out_gridSearch(X,Y,penaltyX, penaltyY) #Fast
    Res = KFold_CV(X_tr,Y_tr,penaltyX, penaltyY) #Relatively fast
    #Check the best combination of hyperparam
    Coordinates = which(Res==max(Res), arr.ind=T)
    best_pX = penaltyX[Coordinates[1]]
    best_pY = penaltyY[Coordinates[2]]
    print( paste0("Best hyperparameters: omics_penalty - ", best_pX, "   PhIP-Seq_penalty - ", best_pY ) ) 
    ###RUN CCA####
    print("Running CCA")
    #Run with the whole dataset
    if ( dim(X)[2]< 10 ){ N_component = dim(X)[2] 
    } else {    N_component = 10 }
    CCA_res = run_sparseCCA(X%>% apply(2, as.numeric), Y%>% apply(2, as.numeric), N_component, best_pX, best_pY)
    CCA_res[[1]]$cors #canonical correlations
    #Test in test set, per component
    print("Test CCA")
    Stat_test = tibble()
    Test_scores = tibble()
    for (i in 1:N_component){
      X_te %>% apply(2, as.numeric) %*% CCA_res[[1]]$u[,i] -> X_score
      Y_te %>% apply(2, as.numeric) %*% CCA_res[[1]]$v[,i] -> Y_score
      Test_scores %>% rbind(. , tibble(Component = i, Component_protein = X_score, Component_phipseq = Y_score) ) ->  Test_scores
      cor.test(X_score, Y_score, method = "pearson") -> Correlation_ho
      cor.test(X_score, Y_score, method = "spearman") -> Correlation_ho2
      tibble(Component = i, Cor_coefficient_p = Correlation_ho$estimate, Cor_97.5_p=Correlation_ho$conf.int[2], Cor_2.5_p =Correlation_ho$conf.int[1] ,P_p = Correlation_ho$p.value, Cor_coefficient_s = Correlation_ho2$estimate, Cor_97.5_s=Correlation_ho2$conf.int[2], Cor_2.5_s =Correlation_ho2$conf.int[1] ,P_s = Correlation_ho2$p.value) %>% rbind(Stat_test, .) -> Stat_test
    }  
    print("Components stats (based on held-out data)")
    print(Stat_test)
    #Quantify how many non-0 data we have
    avg_proteins <- get_avg_features(CCA_res[[1]]$u, N_component)
    avg_AB <- get_avg_features(CCA_res[[1]]$v, N_component)
    print( paste0("Average number of non-0 proteins per component:", avg_proteins  ) )
    print( paste0("Average number of non-0 peptides per component:", avg_AB  ) )
    
    #Put the non-0 data in a long table
    print("Preparing final table")
    Table_correspondence = tibble()
    for (i in seq(1, N_component )){
      I = CCA_res[[1]]$u[,i] 
      I = I[!I==0]
      I2 = CCA_res[[1]]$v[,i] 
      I2 = I2[!I2==0]
      tibble(Component = i, Data = "Antibody", Feature = names(I2), Coefficient = I2 ) -> AB_tibble
      tibble(Component = i, Data = "Protein", Feature = names(I), Coefficient = I ) -> Protein_tibble
      rbind(AB_tibble, Protein_tibble) %>% rbind(Table_correspondence, . ) -> Table_correspondence
      
    }
    
    #Add peptide annotation
    readxl::read_excel("~/Resilio Sync/Antibodies_WIS (1)/Results/Supplemenatry_Tables/SupplementaryTable1.xlsx", sheet="SupTable1.1") %>% rename(Feature = Peptide) %>% select(Feature, full.name, Taxa) -> Annotation
    left_join(Table_correspondence, Annotation) -> Table_correspondence
    
    
    CCA_res[[2]] %>% as_tibble() %>% gather(Component, Component_protein) -> I1
    as_tibble(CCA_res[[3]]) %>% gather(Component, Component_phipseq ) -> I2 
    cbind(I1 , select(I2, -Component) ) %>% as_tibble() -> Components1
    Components1 %>% as_tibble() %>% mutate(Component = Components1$Component %>% sapply(function(x){ str_replace(x, "V", "Component_")} ) ) %>%
      ggplot(aes(x=Component_protein, y=Component_phipseq)) + geom_point() + theme_bw()  + facet_wrap(~Component, scales = "free") + geom_smooth(method = "lm") + ggtitle("Train correlations") -> Plot_train
    
    Test_scores  %>% mutate(Component = paste0("Component_", Component) ) %>%  ggplot(aes(x=Component_protein, y=Component_phipseq)) + geom_point() + theme_bw()  + facet_wrap(~Component, scales = "free") + geom_smooth(method = "lm") + ggtitle("Test correlations") -> Plot_test
    
    
    return( list(Table_correspondence, Stat_test, CCA_res, All_together$ID, list(Plot_train, Plot_test), Test_scores, list(X,Y) ) )
    
  }
  
  Remove_variability = function(DF, Covariates_df = Covariates, Covariates_n = c("Age" , "Sex")){
    
    DF %>% left_join(. , select(Covariates_df, c(Covariates_n, "ID") ))  -> To_correct
    DF_corrected = NULL  
    
    for (Feature in colnames(DF)){
      if (Feature == "ID"){ next }
      To_correct %>% select(c("ID", Covariates_n, Feature) ) %>% drop_na() -> To_correct_m 
      lm( paste0("`", Feature, "`~ ", paste(Covariates_n, collapse="+" ) ), To_correct_m) -> M
      Corrected = M$resid
      Corrected = tibble(Corrected)
      colnames(Corrected) = Feature
      Corrected %>% mutate(ID = To_correct_m$ID) -> Corrected
      if (is.null(DF_corrected) == TRUE){ DF_corrected = Corrected
      } else { full_join(DF_corrected, Corrected, by="ID") -> DF_corrected }
    }
    DF_corrected %>% as_tibble()  -> DF_corrected
    return(DF_corrected)
  }
  
  Run_gsea = function(RNA_results,pathways, Component_n=1, Data_type = "Protein"){
    RNA_results[[1]] %>% filter(Component == Component_n) %>% arrange(Coefficient) %>% filter(Data == Data_type) -> R1
    #Prepare ranks
    if (Data_type == "Protein"){
      genes_f %>% left_join(R1) %>% mutate(Coefficient = ifelse(is.na(Coefficient), 0, Coefficient) )  %>% mutate(Feature = ensembl_gene_id) -> For_rank
      ranks = For_rank$Coefficient
      names(ranks) = as.character(For_rank$Feature)
    } else {
      ranks = R1$Coefficient
    }
    #Run GSEA
    fgseaRes <- fgsea(pathways, ranks, maxSize=1000, minSize=2,nproc=1,nperm = 1000)
    fgseaRes %>% as_tibble() %>% arrange(pval) %>% print()
    return(fgseaRes)
  }
  
  Run_ora = function(RNA_results,pathways,Universe, Component_n=1, Data_type = "Protein"){
    RNA_results[[1]] %>% filter(Component == Component_n) %>% arrange(Coefficient) %>% filter(Data == Data_type) -> R1
    #Prepare ranks
    if (Data_type == "Protein"){
      genes_f %>% left_join(R1) %>% mutate(Coefficient = ifelse(is.na(Coefficient), 0, Coefficient) ) %>% mutate(Feature = ensembl_gene_id) -> For_rank
    } else { 
      For_rank  = R1
    }
    Positive = filter(For_rank, Coefficient > 0)$Feature 
    Negative = filter(For_rank, Coefficient < 0)$Feature 
    #Run ORA
    if (length(Positive) > 0){
      ora_results1 <- enricher(  gene =  Positive ,  pvalueCutoff = 0.1,  pAdjustMethod = "BH", universe = Universe, minGSSize = 2, TERM2GENE = pathways,maxGSSize = 5000)
      ora_results1@result %>% as_tibble() %>% mutate(Direction = "Positive") %>% print()
    } else { ora_results1 = NA }
    if (length(Negative) > 0 ){
      ora_results2 <- enricher(  gene =  Negative ,  pvalueCutoff = 0.1,  pAdjustMethod = "BH", universe = Universe, minGSSize = 2, TERM2GENE = pathways, maxGSSize = 5000)
      ora_results2@result %>% as_tibble() %>% mutate(Direction = "Negative") %>% print()
      
    } else { ora_results2 = NA}
    
    
    
    return(list(ora_results1, ora_results2))
  }
  
  InvRank =  function(x){
    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  }
  index_to_mean <- function(my_index, my_mean){ return(my_mean[my_index]) }
  Quantile_normalization = function(X, means=NA){
    X2 =  t(X)
    X_rank = apply(X2,2,rank,ties.method="min")
    if (class(Means) == "logical" ){
      Means = X2 %>% apply(., 2, sort) %>% apply(df_sorted, 1, mean)
    }
    X_final <- apply(X_rank, 2, index_to_mean, my_mean=Means)
    return(X_final, Means)
  }
  
    Run_permanova_covariates = function(X, Covariates, Cov_names, Permutation_n=500){
    X  %>% select(-ID) %>% apply(2, as.numeric) %>% as_tibble() %>% mutate(ID = X$ID)  %>% left_join(., select(Covariates, c(Cov_names)) ) %>% drop_na() -> P_temp
    P_temp %>% select( - Cov_names ) %>% apply(2, scale)  %>% vegan::vegdist(. , method = "euclidean" ) -> Distance
    Formula = as.formula( paste0("Distance ~ ", paste(Cov_names[2:length(Cov_names)] , collapse="+")))
    vegan::adonis2(Formula  , P_temp, by = "margin", permutations = Permutation_n  ) -> Results_permanova
    Results_permanova %>% as.data.frame() %>% rownames_to_column("Feature") %>% as_tibble() %>% filter(Feature %in% Cov_names) %>% mutate(Permutation_n = Permutation_n) %>% return()
  }
    For_halla = function(X){
    X %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(ID = colnames(select(X, -ID)), .before=1) %>% `colnames<-`( c("ID", X$ID )) -> X2
    return(X2)
  }

    Check_cluster = function(C, Original_X = ImmunoMatrix, Original_Y = Cytokines , Test_IDs = Test_data, Cov = Covariates, Confounder_names= c("Age", "Sex", "smk_now", "Lymph", "Eryth", "Mono", "Neutro", "Thrombo", "Eosino", "Blood_baso", "estrogens"), Make_Plot=T ){
    #Check halla results
    Results = tibble()
    Plots = list()
    for (Cluster in C$cluster_rank){
      C %>% filter(cluster_rank == Cluster) -> R
      for (A in str_split(R$cluster_X, ";")[[1]] ){
        for (B in str_split(R$cluster_Y, ";")[[1]] ){
          left_join(left_join(Original_X, Original_Y, by="ID"), Covariates, by="ID", suffix=c("", ".rep" ))  %>% select(c("ID", A,B, Confounder_names )) %>% drop_na() -> D
          D %>% filter(ID %in% Test_data) %>% as.data.frame()   -> to_Test
          D %>% filter(! ID %in% Test_data) %>% as.data.frame()  -> to_Train
          pcor.test(x = to_Train[[2]] , y = to_Train[[3]] , z =to_Train[4:dim(to_Train)[2]] , method = "spearman") -> Correlation_train
          pcor.test(x = to_Test[[2]] , y = to_Test[[3]] , z =to_Test[4:dim(to_Test)[2]] , method = "spearman") -> Correlation_test
          Prev_train = sum(to_Train[[2]] != 0) / length(to_Train[[2]])
          Prev_test = sum(to_Test[[2]] != 0) / length(to_Test[[2]])
          tibble(Cluster = Cluster, X=A, Y=B, rho_train = Correlation_train$estimate, P_train=Correlation_train$p.value ,N_train=dim(to_Train)[1], Prevalence_train=Prev_train, rho_test = Correlation_test$estimate, P_test=Correlation_test$p.value  , N_test=dim(to_Test)[1], Prevalence_test=Prev_test) -> Result_temp
          rbind(Results, Result_temp) -> Results
	 if (Make_Plot == F){ next } 
          #plotting
          rbind( to_Train %>% mutate(Dataset = "Train"), to_Test %>% mutate(Dataset = "Test") ) -> To_plot
          ggplot() + geom_boxplot(aes(x=as.factor(as_vector(select(To_plot, A))), y =as_vector(select(To_plot, B)), fill=To_plot$Dataset ), outlier.shape = NA) +  ggforce::geom_sina(aes(x=as.factor(as_vector(select(To_plot, A))), y =as_vector(select(To_plot, B)), col=To_plot$Dataset ), alpha=0.5)   + theme_bw() + xlab(A) + ylab(B) -> Plot
          Plots = append(Plots, print(Plot))
        }  
      }    
      
    }
    if (Make_Plot == F){ return(Results) }
    return(list(Results, Plots))
  }

   Format_rnaseq = function(normalized_counts){
    normalized_counts %>% select(-gene) %>% t() %>% as_tibble() -> normalized_counts_t
    colnames(normalized_counts_t) = normalized_counts$gene
    normalized_counts_t %>% mutate(ID = colnames(select(normalized_counts, -gene)), .before=1) -> normalized_counts_t
    return(normalized_counts_t)
  }
