library(tidyverse)
library(caret)


Impute_telomere_data = function(Telomere_markers){
  set.seed(100)
  Tels = c("TL_Lymphocytes","TL_Granulocytes","`TL_Naive_T-cells`","`TL_Memory_T-cells`", "`TL_B-cells`","`TL_NK-cells`" )#colnames(Telomere_markers)[2: dim(Telomere_markers)[2]]
  Telomere_markers %>% drop_na() -> Complete_telomere
  N = dim(Complete_telomere)[1]
  sample(x= N, replace=F, size= N*0.2) -> Test_samples
  Complete_telomere[Test_samples, 1:(length(Tels)+1)] -> TEST
  Complete_telomere %>% filter(!ID %in% TEST$ID) -> TRAIN
  Model_list = list() ; n =1 ; Acc_tibble = tibble()
  for (i in Tels){
    Formula = as.formula( paste( c( i, "~", paste(Tels[! Tels == i], collapse="+")), collapse= " " )  )
    lm(Formula, TRAIN) -> Model
    predict.lm(Model, TEST) -> Predicted
    Real = as_vector(TEST[colnames(TEST) == str_replace_all(i,"`","") ])
    caret::R2(obs = Real, pred = Predicted) -> pred_r2
    rbind(Acc_tibble, tibble(Telomere= i, R2_test = pred_r2 )) -> Acc_tibble
    Model_list[[n]] = Model
    n = n+1
    
    Predicted = predict.lm( Model , Telomere_markers[is.na(Telomere_markers[colnames(Telomere_markers) == str_replace_all(i,"`","") ]), 1:dim(Telomere_markers)[2]]   )
    Telomere_markers[is.na(Telomere_markers[colnames(Telomere_markers) == str_replace_all(i,"`","") ]),  colnames(Telomere_markers) == str_replace_all(i,"`","")]   = as.numeric(Predicted)
  }
  
  return(list(Telomere_markers, Acc_tibble))
  
}


Clean_and_ALR = function(all_cell_types){
  cell_matrix <- as.matrix(subset(all_cell_types, select=-ID))
  #Some predictions go below 0, make them 0 and go back to 100%
  cell_matrix[cell_matrix<0] = 0
  cell_matrix %>% apply(1, function(x){ 100*(x/sum(x)) } ) %>% t() %>% as_tibble() %>% `colnames<-`(colnames(cell_matrix))  -> cell_matrix 
  #ALR using `Monocytes (CD14+)` as denominator
  Denominator = cell_matrix$`Monocytes (CD14+)`
  cell_matrix %>% dplyr::select(-`Monocytes (CD14+)`) %>% apply(2, function(x){ ALR = log10( x/Denominator )  ; ALR[ abs(ALR) == Inf] = NA  ; return(ALR) }  ) %>% as_tibble() %>% mutate(ID =all_cell_types$ID ) -> cell_matrix_alr
  return(cell_matrix_alr)
}

Clean_and_CLR = function(all_cell_types){
	cell_matrix <- as.matrix(subset(all_cell_types, select=-ID))
	cell_matrix[cell_matrix<0] = 0
	cell_matrix %>% apply(1, function(x){ 100*(x/sum(x)) } ) %>% t() %>% as_tibble() %>% `colnames<-`(colnames(cell_matrix))  -> cell_matrix
	#Compute geometric mean and use as denominador
	#Many cell proportions are overlapping: I took a look and seems that cells can roughly be devided in lymphocytes (T-cells, B-cells and NK-cells), Granulocytes (no subdivision in the data) and Monocyes (CD14+). All the other cell types belong to one of them. Thus, I will use the geometric mean among them as the denominaotr
	c("Granulocytes", "Lymphocytes", "Monocytes (CD14+)" ) -> Non_overlapping
	Denominator = apply( select(cell_matrix, Non_overlapping)   ,1, function(x){ exp(mean(log(x)))   }  ) 
	cell_matrix %>% apply(2, function(x){ CLR = log10( x/Denominator ) ; CLR[ abs(CLR) == Inf] = NA ; return(CLR) } )  %>% as_tibble() %>% mutate(ID =all_cell_types$ID ) -> cell_matrix_clr
	return(cell_matrix_clr)

}

LDA_function = function(Train, Test, Print=T){
  #### LDA ####
  #1, Fit
  lda(CMV_prediction~., Train) -> LDA
  #2. Compute LDA component and plot
  dplyr::select(Train, -CMV_prediction) ->  M
  scores = as_tibble(as.matrix(M) %*% LDA$scaling) %>% mutate(CMV_prediction = Train$CMV_prediction)
  scores  %>% ggplot(aes(x=as.factor(CMV_prediction), y=LD1)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw()
  #3. Get accuracy
  p1 <- predict(LDA, Train)$class
  tab <- table(Predicted = p1, Actual = Train$CMV_prediction)
  sum(diag(tab))/sum(tab)
  if (Print == T ){
  #4. statistical test
  glm(as.factor(CMV_prediction) ~ LD1, scores, family=binomial() ) %>% summary() %>% print()
  }
  #5. Fit in test, compute accuracy and test
  p2 <- predict(LDA, Test)$class
  tab <- table(Predicted = p2, Actual = Test$CMV_prediction)
  sum(diag(tab))/sum(tab)
  scores_test = as_tibble(as.matrix(dplyr::select(Test, -CMV_prediction)) %*% LDA$scaling) %>% mutate(CMV_prediction = Test$CMV_prediction) 
  #plot
  scores_test %>% ggplot(aes(x=as.factor(CMV_prediction), y=LD1)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw()
  #stat test
  if (Print == T){
  glm(as.factor(CMV_prediction) ~ LD1, scores_test, family=binomial() ) %>% summary() %>% print()
  }
  return( list(LDA, scores, scores_test  ) )

}

RF_function =  function(Train, Test){
  Train %>% drop_na() -> Train
  Test %>% drop_na() -> Test
  set.seed(999)
  ntree <- 500
  tunegrid <- expand.grid(.mtry=  4:(sqrt(ncol(Train)-1)+10) )
  control <- trainControl(method="repeatedcv", number=5, repeats=3, search="grid")

  model <- train(as.factor(CMV_prediction)~., data = Train, method ="rf", metric = "Accuracy",tuneGrid=tunegrid, trControl=control, ntree=ntree ) #ranger is a faser implementation of randomforst


  prediction_for_roc_curve <- predict(model, dplyr::select(Test, -CMV_prediction), type ="prob")

  roc(Test$CMV_prediction , prediction_for_roc_curve[,2]) -> ROC
  plot.roc(ROC)
  print(ROC)
  
}

Permutations_2 = function( DF ,Covariates, N_Permutations, Real_P, Interaction=F ){
  #We will calculate FDR values through permutation
  Null_distribution = c()
  for (Permutation_n in seq(N_Permutations) ){
    Covariates_permuted = Covariates
    Covariates_permuted$CMV_prediction = sample(Covariates_permuted$CMV_prediction, size = dim(Covariates)[1], replace=F )
    if (Interaction == T){
      Interaction_association(DF, Covariates_permuted) %>% filter(grepl(":", Regressor)) -> Permuted_results
    } else { Basic_association(DF, Covariates_permuted)  -> Permuted_results  }
    Null_distribution = c(Null_distribution, Permuted_results$`Pr(>|t|)`)
  }
  FDRs = c()
  for (n  in Real_P){
  FDR = sum(Null_distribution <= n)/ N_Permutations
  if( FDR > 1){ FDR = 1 }
  FDRs = c(FDRs, FDR)
  }
  return(FDRs)
}


make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

Geom_mean = function(x){
    exp(mean(log(x)))
  }

 CLR = function(D){
    log2(D / Geom_mean(D) )
    
  }

