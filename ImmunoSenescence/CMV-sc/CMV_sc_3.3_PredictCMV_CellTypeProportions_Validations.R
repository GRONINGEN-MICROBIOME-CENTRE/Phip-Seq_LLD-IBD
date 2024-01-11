#!/usr/bin/env Rscript

# # setting working directory (cluster or local)
# path_cluster <- '/gpfs/projects/bsc83/Projects/scRNAseq/aripol1/CMV-sc/'
# path_local <- '/home/aripol1/Desktop/bsc/Projects/scRNAseq/aripol1/CMV-sc/'
# 
# if(file.exists(path_cluster)){
#   setwd(paste(path_cluster))
# }else if(file.exists(path_local)){
#   setwd(paste(path_local))
# }

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--cell_level"), action="store", default="cell_type", type='character',
              help="predicted.celltype.l1 or cell_type"),
  make_option(c("--phenotype"), action="store", default="CMV_status", type='character',
              help="predicted.celltype.l1 or cell_type"),
  make_option(c("--train_dataset"), action="store", default="v2,v3", type='character',
              help="Train dataset."),
  make_option(c("--validation_dataset"), action="store", default="pilot3", type='character',
              help="Test dataset."),
  make_option(c("--cmv_metadata"), action="store", default="data/Metadata_matchLLD_CMV.tsv", type='character',
              help="CMV metadata"),
  make_option(c("--in_dir"), action="store", default='data/CMV_sc_1_CODA/Azimuth_celltype_proportions', type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_3.3_PredictCMV_CellTypeProportions_Validations', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(data.table))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(tibble))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(broom))
shhh(library(broom.mixed))
shhh(library(ggplot2))
shhh(library(caret))
shhh(library(glmnet))
shhh(library(glmnetUtils))
shhh(library(pROC))
shhh(library(yardstick))
shhh(library(psych))

#################### Define functions #################### 
# Transform data (CLR)
## accessory functions
Geom_mean <- function(x){
  exp(mean(log(x)))
}
CLR <- function(D){
  log2(D / Geom_mean(D))
}

## main function
clr_func <- function(df){
  # Filter out cell types based on minimum nCells/donor: Long to wide --> rows = cell_type, cols = donor, value = n
  ## Check nCells/donor
  DF_n <- reshape2::dcast(df, donor ~ celltype, value.var = "n")
  rownames(DF_n) <- DF_n[,1]
  DF_n <- DF_n[,-1]
  DF_n <- t(DF_n)
  DF_n[is.na(DF_n)] <- 0
  
  ## Define filter: min 5 donors with at least 5 cells/donor
  nDonors_filt.byCT <- apply(DF_n, 1, function(x){sum(x>5)})
  CT_filt <- names(nDonors_filt.byCT[nDonors_filt.byCT>=5])
  
  # Long to wide --> rows = cell_type, cols = donor, value = freq
  DF_proportions <- reshape2::dcast(df, donor ~ celltype, value.var = "freq")
  rownames(DF_proportions) <- DF_proportions[,1]
  DF_proportions <- DF_proportions[,-1]
  DF_proportions <- t(DF_proportions)
  DF_proportions[is.na(DF_proportions)] <- 0
  
  # Add pseudocount
  pseudocount <- 1/3*min(DF_proportions[!DF_proportions==0])
  
  # CLR
  DF_proportions %>% 
    apply(2, function(x){CLR(x+pseudocount)}) -> DF_proportions.clr
  
  # # Check
  # apply(DF_proportions, 1, summary)
  # apply(DF_proportions.clr, 1, summary)
  
  ## Apply previously defined filter
  DF_proportions.clr.filt <- DF_proportions.clr[rownames(DF_proportions.clr)%in%CT_filt,]
  
  return(DF_proportions.clr.filt)
}

# LMER function (regress out date and dataset)
## main function
lmer_by_ct <- function(cell_type, regress_out, df, md){
  print(cell_type)
  # Data
  vec_i <- df[rownames(df)==cell_type,]
  df_i <- as.data.frame(vec_i)
  colnames(df_i) <- 'freq'
  df_i$donor <- rownames(df_i)
  df_i <- merge(df_i, md, by = 'donor')
  rownames(df_i) <- df_i$donor
  df_i <- df_i[,-1]
  Group_order.CMV_status <- c('0','1')
  Group_order.Sex <- c('M','F')
  if('Sex'%in%colnames(df_i)){
    df_i[['Sex']] <- factor(df_i[['Sex']],
                            levels = Group_order.Sex)
  }
  if('CMV_status'%in%colnames(df_i)){
    df_i[['CMV_status']] <- factor(df_i[['CMV_status']],
                                   levels = Group_order.CMV_status)
  }
  
  # Formula
  # covs_fixed <- covs_df[covs_df$type=='fixed',]$covariate
  # covs_random <- covs_df[covs_df$type=='random',]$covariate
  # if(length(covs_fixed)>0){
  #   fixed_fmla <- paste(covs_fixed,collapse='+')
  # }
  # random_fmla <- NULL
  # if(length(covs_random)>0){
  #   random_fmla <- paste(paste0('(1|',covs_random,')'),collapse='+')
  # }
  # fmla <- paste(c(fixed_fmla, random_fmla), collapse = '+')
  # if(interaction){
  #   vars_interaction <- covs_df[!is.na(covs_df$interaction) & covs_df$interaction=='interaction',]$covariate
  #   interaction_fmla <- paste(vars_interaction, collapse = ':')
  #   fmla <- paste0(c(fmla, interaction_fmla), collapse = '+')
  # }
  
  covs_random <- regress_out
  random_fmla <- paste(paste0('(1|',covs_random,')'),collapse='+')
  fmla <- paste0('freq ~ ',random_fmla)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer and extracting residuals: ',fmla))
  mod <-  lmerTest::lmer(form, data = df_i)
  mod.residuals <- residuals(mod)
  return(mod.residuals)
}

# Predict CMV status
## accessory functions
### report function
report_model <- function(train, test, phe){
  # Select train/test data
  train_lasso <- train
  test_lasso <- test
  
  # Report
  train$type <- 'train'
  test$type <- 'test'
  df_to_report <- rbind(train, test)
  group_vars <- c('type', phe)
  df_to_report %>% 
    group_by_at(group_vars, .drop = FALSE) %>%
    summarise(n = n()) %>%
    mutate(prop = n / sum(n)) %>% as.data.frame() -> df.c
  
  return(df.c)
}

### fit function
# Lasso_n = 4 #Number of folds to split the training data  ( (K-1)/K from total data ) for cross-validation to pick up the best hyperparameter lambda (regularization strength)
fit_model <- function(train, test, phe, Lasso_n = 4){
  print("Training Lasso...")
  
  # Model 2: Logistic regression with Lasso normalization, using all covariates ##family= binomial => logistic regression, alpha=1 => lasso
  # http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/153-penalized-regression-essentials-ridge-lasso-elastic-net/#lasso-regression
  train_lasso <- train
  test_lasso <- test
  
  # Find the best lambda using cross-validation
  ## define formula
  fmla <- paste0(phe, ' ~.')
  form <- as.formula(fmla)
  print(paste0('Fitting Logistic regression with Lasso normalization: ',fmla))
  
  ## alpha is an hyperparameter that lets you choose between l1 and l2 penalty. An alpha of 1 corresponds to lasso regression
  model_lasso = glmnetUtils::cv.glmnet(form, 
                                       train_lasso, 
                                       alpha = 1, 
                                       nfolds = Lasso_n, 
                                       family = "binomial", 
                                       type.measure="class")
  
  ## best lambda --> features' importance
  best_lambda <- model_lasso$lambda.min
  coef_matrix <- coef(model_lasso, s = best_lambda)
  important_features.idx <- which(coef_matrix != 0)
  important_features.vec <- coef_matrix[important_features.idx,]
  if(length(important_features.idx)==1){
    names(important_features.vec) <- rownames(coef_matrix)[important_features.idx]
  }
  important_features.df <- data.frame(feature = names(important_features.vec),
                                      value = unname(important_features.vec))
  important_features.df <- important_features.df[order(-important_features.df$value),]
  
  # Fit the final model on the training data --> Predict on newdata
  print("Testing Lasso...")
  test_lasso.newdata <- test_lasso %>% select(-all_of(phe)) %>% as.data.frame()
  lasso_probabilities = predict(model_lasso, 
                                newdata = test_lasso.newdata, 
                                type="response", 
                                s="lambda.min" )
  
  lasso_predictions = rep(0, dim(test_lasso)[1])
  lasso_predictions[lasso_probabilities>.5] = 1
  lasso_probabilities.vec <- lasso_probabilities[,1]
  
  # Performance metrics (https://www.shiksha.com/online-courses/articles/roc-auc-vs-accuracy/#Difference-between-ROC---AUC-and-Accuracy)
  ## ROC-AUC = TPR(True Positive Rate) = TP / (TP + FN)FPR( False Positive Rate) = FP / (FP + TN)
  ## Accuracy = (TP + TN) / (TP + TN + FP + FN)
  ## ROC AUC compares the relation between True Positive Rate and False Positive Rate, while Accuracy is simply the percentage of correct predictions.
  ## But in other cases (for extremely Imbalanced data), ROC AUC may be more important.
  
  ## Prepare table
  performances.df <- as.data.frame(lasso_probabilities)
  performances.df$lambda.min_diff <- 1-performances.df$lambda.min
  class_cnames <- paste(phe,rev(levels(train_lasso[[phe]])),sep='')
  colnames(performances.df) <- class_cnames
  truth_vec <- as.character(test_lasso[[phe]])
  predicted_vec <- as.character(lasso_predictions)
  performances.df$truth <- truth_vec
  performances.df$truth <- factor(performances.df$truth,
                                  levels = c('1','0'))
  performances.df$predicted <- predicted_vec
  performances.df$predicted <- factor(performances.df$predicted,
                                      levels = c('1','0'))
  
  ## Metrics
  # kappa.out <- cohen.kappa(performances.df[,c('truth','predicted')])
  # kappa_df <- data.frame(metric = 'kappa', 
  #                        estimate = kappa.out$kappa)
  conf_mat.out <- conf_mat(performances.df, truth = truth, estimate = predicted)
  classification_metrics <- metric_set(accuracy, mcc, f_meas)
  classification_metrics.out <- classification_metrics(performances.df, truth = truth, estimate = predicted)
  class_name <- class_cnames[1]
  roc_auc.out <- roc_auc(performances.df, truth, class_name)
  metrics_df <- rbind(classification_metrics.out, roc_auc.out) %>% as.data.frame()
  metrics_df <- metrics_df[,-2]
  colnames(metrics_df) <- gsub('\\.','',colnames(metrics_df))
  # metrics_df <- rbind(metrics_df, kappa_df)
  metrics_df$Model <- "Logistic_Lasso"
  metrics_wide <- dcast(metrics_df, Model ~ metric, value.var='estimate')
  
  # Curves
  curve_list <- list(roc_curve = roc_curve(performances.df, truth, class_name),
                     gain_curve = gain_curve(performances.df, truth, class_name),
                     lift_curve <- lift_curve(performances.df, truth, class_name),
                     pr_curve <- pr_curve(performances.df, truth, class_name))
  
  # Output
  out <- list(performance = metrics_wide,
              features = important_features.df,
              curves = curve_list,
              confusion_matrix = conf_mat.out)
  return(out)
}

## main function
validation.func <- function(train_ds, validation_ds, N, l, phenotype){
  print(paste0('# Iteration (N): ', as.character(N)))
  print(paste0('# Training dataset: ', train_ds))
  
  # Data
  train_test.df <- droplevels(l[[train_ds]])
  
  # Split training data in train/test
  trainIndex <- createDataPartition(train_test.df$CMV=="1", p = .8, 
                                    list = FALSE, 
                                    times = 1)
  train_df = train_test.df[trainIndex,]
  test_df = train_test.df[-trainIndex,]
  
  # Validation: external data (validation) or 20% of the training dataset (test)
  validation.df <- test_df
  out.sdir <- paste0(out.dir, 'test/')
  if(!is.null(validation_ds)){
    print(paste0('# Validation dataset: ', validation_ds))
    out.sdir <- paste0(out.dir, 'validation/')
    validation.df <- droplevels(l[[validation_ds]])
  }
  if(!dir.exists(out.sdir)){dir.create(out.sdir, recursive = T)}
  table(train_df$CMV)
  table(test_df$CMV)
  
  # Apply report on the different folds
  report_model.df <- report_model(train = train_df, 
                                  test = validation.df, 
                                  phe = phenotype)
  report_model.df$Iteration <- paste0('Iteration',as.character(N))
  cat('\n')
  
  # Apply fit model on the different folds (validation)
  fit_model.res <- fit_model(train = train_df,
                             test = validation.df,
                             phe = phenotype)
  
  out_names <- c('performance', 'features', 'curves', 'confusion_matrix')
  fit_model.out <- sapply(out_names, function(i) fit_model.res[[i]], simplify = FALSE)
  fit_model.out.p_f <- fit_model.out[names(fit_model.out)%in%c('performance', 'features')]
  fit_model.list <- lapply(fit_model.out.p_f, function(x){
    x$Iteration <- paste0('Iteration',as.character(N))
    return(x)
  })
  
  # Save curves per iteration
  fit_model.out.c <- fit_model.out[names(fit_model.out)%in%c('curves')]
  it_tag <- paste0('Iteration',as.character(N))
  curves_model.list <- list()
  curves_model.list[[it_tag]] <- fit_model.out.c
  curves_model.fn <- paste0(out.sdir, it_tag, '.curves.rds')
  saveRDS(curves_model.list, curves_model.fn)
  
  # Save conf matrix per iteration
  fit_model.out.cm <- fit_model.out[names(fit_model.out)%in%c('confusion_matrix')]
  it_tag <- paste0('Iteration',as.character(N))
  cm_model.list <- list()
  cm_model.list[[it_tag]] <- fit_model.out.cm
  cm_model.fn <- paste0(out.sdir, it_tag, '.confusion_matrix.rds')
  saveRDS(cm_model.list, cm_model.fn)
  
  # Output
  report_model.list <- list(report = report_model.df)
  out <- c(fit_model.list, report_model.list)
  
  cat('\n')
  cat('\n')
  return(out)
}

# Get output from prediction
## main function
get_output <- function(i){
  print(i)
  logistic_lasso.res <- logistic_lasso.list[[i]]
  out_names <- c('performance', 'features', 'report')
  logistic_lasso.out <- sapply(out_names, function(i) 
    lapply(logistic_lasso.res, function(iteration) iteration[[i]]), simplify = FALSE)
  logistic_lasso <- lapply(logistic_lasso.out, function(out_name) do.call("rbind",out_name))
  logistic_lasso.fn <- paste0(out.dir, i, '/logistic_lasso.rds')
  print(paste0('Saving logistic lasso results: ', logistic_lasso.fn))
  saveRDS(logistic_lasso, logistic_lasso.fn)
  return(logistic_lasso)
}

## Boxplot (among datasets)
## main function
bp_performance <- function(p_metrics, cols_vec, df, width_var = 6){
  df_melt <- melt(df, measure.vars = names(p_metrics))
  df_melt <- droplevels(df_melt[df_melt$variable%in%names(p_metrics),])
  df_melt %>% 
    group_by(Dataset, variable) %>%
    summarise(mean=mean(value),
              median=median(value)) %>% as.data.frame() -> df_melt.stats
  df_stats <- df_melt.stats[,c('Dataset', 'variable', 'median')]
  colnames(df_stats)[colnames(df_stats)=='median'] <- 'value' 
  df_tmp <- df_melt[,c('Dataset', 'variable', 'Iteration', 'value')]
  df_all <- df_tmp
  df_all$Iteration <- factor(df_all$Iteration,
                             levels = names(cols_vec))
  df_all$Metric <- unname(p_metrics[df_all$variable])
  df_all$Metric <- factor(df_all$Metric,
                          levels = unname(p_metrics))
  df_all$Dataset <- factor(df_all$Dataset,
                           levels = c('test', 'validation'))
  p <- ggplot(df_all, aes(x=Dataset, y=value)) +
    geom_boxplot(outlier.alpha = 0, fatten = 1) +
    geom_jitter(stroke=NA, width = 0.2, aes(color = Iteration)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    theme_bw() +
    ylab('Performance metric value') + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) + 
    scale_color_manual(values = cols_vec) + 
    facet_grid(. ~ Metric,
               scales = 'free')
  metrics_tag <- paste(names(p_metrics),collapse='.')
  p.fn <- paste0(out.dir, metrics_tag, '.png')
  ggsave(p.fn, p, width = width_var, height = 4)
  return(p)
}

#################### Set Variables and load Data #################### 
# Datasets
# opt$train_dataset <- 'v2'
# opt$train_dataset <- 'v3'
# opt$validation_dataset <- 'pilot3'
train_dataset.vec <- unlist(str_split(opt$train_dataset, ','))
datasets <- c(train_dataset.vec, opt$validation_dataset)
datasets_tag <- paste(datasets, collapse='_')

# Report
print('############################')
print(paste0('Datasets: ', datasets_tag))
print(paste0('Cell level: ', opt$cell_level))
print('############################')
cat('\n')

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', datasets_tag, '/', opt$cell_level, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main results directory: ',out.dir))

# Proportions
## read data
props_fn <- paste0(in.dir, 'proportions.rds')
print(paste0('Reading proportions file in: ',props_fn))
props.df <- readRDS(props_fn)
props.df <- droplevels(props.df[props.df$dataset%in%datasets,])

## change cell types names
cell_types <- unique(props.df$cell_type)
cell_type.vec <- gsub(' ', '_', cell_types)
names(cell_type.vec) <- cell_types
props.df$celltype <- unname(cell_type.vec[props.df$cell_type])
props.df <- droplevels(props.df)

## donor metadata
donor_md <- unique(props.df[,c('donor', 'Gender', 'Age', 'date', 'lane', 'dataset')])
colnames(donor_md)[colnames(donor_md)=='Gender'] <- 'Sex'
donor_md$lane <- as.factor(donor_md$lane)
donor_md$date <- as.factor(donor_md$date)

# Summarise info by celltype
props.df %>%
  group_by(celltype, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count.df
count.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats.df
celltype.stats.df$celltype.label <- paste0(celltype.stats.df$celltype,
                                           ' (n=', celltype.stats.df$n, ';freq=', as.character(round(celltype.stats.df$freq,2)), ')')

# Prepare variables
contrast_vec <- c('SexF', 'Age', 'CMV_status1')
names(contrast_vec) <- c('Sex', 'Age', 'CMV_status')

# Cell metadata
cell.md_fn <- paste0(in.dir, 'metadata.rds')
print(paste0('Reading proportions file in: ',cell.md_fn))
cell.md <- readRDS(cell.md_fn)
cell.md <- droplevels(cell.md[cell.md$dataset%in%datasets,])
cell_level.vec <- c('predicted.celltype.l1', 'cell_type') 
names(cell_level.vec) <- c('l1.pred', 'l2.pred')
cnames <- c('barcode', 'donor', 'date', names(cell_level.vec))
colnames(cell.md)[colnames(cell.md)=='Gender'] <- 'Sex'
cell.md <- cell.md[,cnames]
colnames(cell.md)[which(colnames(cell.md)%in%names(cell_level.vec))] <- unname(cell_level.vec)
cnames <- c('barcode', 'donor', opt$cell_level)
cell.md <- cell.md[,cnames]
colnames(cell.md)[which(colnames(cell.md)==opt$cell_level)] <- 'cell_type'

# CMV metadata
## read data
cmv_df <- read.delim(opt$cmv_metadata, check.names = FALSE)

## subset by sc-donors
sc_donors <- unique(props.df$donor)
print(paste0('# of sc donors: ', length(sc_donors)))
cmv_df.sc <- droplevels(cmv_df[cmv_df$LLDEEP_ID%in%sc_donors,])
print(paste0('# of match donors: ', nrow(cmv_df.sc)))
cmv_df.sc <- cmv_df.sc[!is.na(cmv_df.sc$CMV_Baseline) | !is.na(cmv_df.sc$CMV_Followup),]
print(paste0('# of match donors (without NAs in CMV_Baseline and CMV_Followup): ', nrow(cmv_df.sc)))

## CMV status
cmv_df.sc$CMV_status <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$CMV_Followup, cmv_df.sc$CMV_Baseline)
cmv_df.sc$Age_diff_SCvsLLD <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$Age_diff_SCvsLLD2, cmv_df.sc$Age_diff_SCvsLLD1)
cmv_df.sc$CMV_status <- as.factor(as.character(cmv_df.sc$CMV_status))
table(cmv_df.sc$CMV_status)
summary(cmv_df.sc$Age_diff_SCvsLLD)

## add CMV metadata to seurat metadata
props.df <- merge(props.df, cmv_df.sc, by.x = 'donor', by.y = 'LLDEEP_ID')
length(unique(props.df$donor))

# nDonors per cell type and CMV status --> to later filter datasets - cell type pairs in the pseudobulk-DEA meta-analysis
props.df %>%
  group_by(celltype, CMV_status, .drop = FALSE) %>%
  dplyr::count() -> donor_count.df

donor_count.df %>%
  group_by(celltype) %>%
  mutate(freq = n / sum(n)) %>% 
  as.data.frame() -> donor.stats.df

celltype.stats.df_link <- celltype.stats.df
colnames(celltype.stats.df_link)[c(2,3)] <- c('nCells', 'freqCells')
colnames(donor.stats.df)[c(3,4)] <- c('nDonors', 'freqDonors')
donor_celltype.stats.df <- merge(celltype.stats.df_link, donor.stats.df, by = 'celltype')
donor_celltype.stats.df <- donor_celltype.stats.df[order(-donor_celltype.stats.df$nCells),]
donor_celltype.stats.fn <- paste0(out.dir, 'donor_celltype.stats.txt')
write.table(donor_celltype.stats.df, donor_celltype.stats.fn, row.names = FALSE, quote = FALSE, sep = '\t')

# Pick donor metadata
donor_md <- unique(props.df[,c('donor', 'Sex.SC', 'Age.SC', 'CMV_status', 'dataset', 'date')])
colnames(donor_md)[colnames(donor_md)%in%c('Sex.SC', 'Age.SC')] <- c('Sex','Age')
donor_md$Sex <- ifelse(donor_md$Sex=='Male', 'M',
                       ifelse(donor_md$Sex=='Female', 'F', donor_md$Sex))
donor_md %>% mutate_if(is.character, as.factor) -> donor_md
donor_md$Sex <- factor(donor_md$Sex,
                       c('M','F'))
donor_md$CMV_status <- factor(donor_md$CMV_status,
                              c('0','1'))

#################### Transform data (CLR) ####################
# Apply function
props_clr.df <- clr_func(props.df)

#################### Data pre-processing ####################
# Apply LMER function to regress out date and dataset
## Prepare variables
celltypes.vec <- sort(apply(props_clr.df, 1, mean),decreasing=T)
celltypes <- names(celltypes.vec)
regress_out.vars <- c('dataset', 'date')

## Get residuals
residuals.list <- sapply(celltypes, function(i) lmer_by_ct(cell_type = i,
                                                           regress_out = regress_out.vars, 
                                                           df = props_clr.df, 
                                                           md = donor_md), simplify = FALSE)
residuals.df <- do.call("rbind", residuals.list)
# residuals.df <- props_clr.df #testing

# Merge with donor metadata
donor_md.fn <- paste0(out.dir, 'donor_md.txt')
write.table(donor_md, donor_md.fn, row.names = FALSE, quote = FALSE, sep = '\t')
residuals_clr.tmp <- as.data.frame(t(residuals.df))
residuals_clr.tmp$donor <- rownames(residuals_clr.tmp)
rownames(residuals_clr.tmp) <- NULL
df_model <- merge(donor_md, residuals_clr.tmp, by = 'donor')
rownames(df_model) <- df_model$donor
df_model <- df_model[,-which(colnames(df_model)=='donor')]
Group_order.CMV_status <- c('0','1')
Group_order.Sex <- c('M','F')
if('Sex'%in%colnames(df_model)){
  df_model[['Sex']] <- factor(df_model[['Sex']],
                              levels = Group_order.Sex)
}
if('CMV_status'%in%colnames(df_model)){
  df_model[['CMV_status']] <- factor(df_model[['CMV_status']],
                                     levels = Group_order.CMV_status)
}

#################### Predict CMV status ####################
# Define type of data (train or validation)
train_dataset <-  paste(train_dataset.vec, collapse='_')
if(length(train_dataset.vec)>1){
  df_model$dataset <- ifelse(df_model$dataset%in%train_dataset.vec, 
                             train_dataset, opt$validation_dataset)
}

# Split data by chemistry (dataset)
df_model.byDS <- split(df_model, df_model$dataset)
df_model.byDS <- lapply(df_model.byDS, function(x) x[,-which(colnames(x)%in%regress_out.vars)])
N.vec <- seq(1:5)

## Apply function
#### Validation: 20% of the initial data (test dataset)
set.seed(123)
logistic_lasso_test.res <- sapply(N.vec, function(iteration) validation.func(train_ds = train_dataset,
                                                                             validation_ds = NULL,
                                                                             N = iteration,
                                                                             l = df_model.byDS, 
                                                                             phenotype = opt$phenotype), simplify = FALSE)

#### Validation: external data (validation)
set.seed(123)
logistic_lasso_val.res <- sapply(N.vec, function(iteration) validation.func(train_ds = train_dataset,
                                                                            validation_ds = opt$validation,
                                                                            N = iteration,
                                                                            l = df_model.byDS, 
                                                                            phenotype = opt$phenotype), simplify = FALSE)

## Retrieve outputs
logistic_lasso.list <- list(test = logistic_lasso_test.res,
                            validation = logistic_lasso_val.res)
logistic_lasso.tmp <- sapply(names(logistic_lasso.list), function(ds) get_output(ds), simplify = FALSE)
items <- c('performance', 'features', 'report')
logistic_lasso.out <- sapply(items, function(i) 
  lapply(logistic_lasso.tmp, function(x) x[[i]]), simplify = FALSE)
logistic_lasso.out <- lapply(logistic_lasso.out, function(x){
  xx <- do.call("rbind", x)
  xx$Dataset <- str_split_fixed(rownames(xx), '\\.', 2)[,1]
  return(xx)
})
logistic_lasso.fn <- paste0(out.dir, 'logistic_lasso.rds')
print(paste0('Saving logistic lasso results: ', logistic_lasso.fn))
saveRDS(logistic_lasso.out, logistic_lasso.fn)

#################### Plots: boxplots (performance metrics) #################### 
# Boxplots
## Variables
# performance_metrics <- c('ROC AUC', 'Cohen\'s kappa', 'F1', 'MCC', 'Accuracy')
# performance_metrics.names <- c('roc_auc', 'kappa', 'f_meas', 'mcc', 'accuracy')
performance_metrics <- c('ROC AUC', 'F1', 'MCC', 'Accuracy')
performance_metrics.names <- c('roc_auc', 'f_meas', 'mcc', 'accuracy')
names(performance_metrics) <- performance_metrics.names
# performance_metrics.names_selected <- c('accuracy', 'roc_auc')
# performance_metrics.selected <- performance_metrics[names(performance_metrics)%in%performance_metrics.names_selected]

cols_it <- c('#264653', '#2a9d8f', '#e9c46a', '#f4a261', '#e76f51')
cols_it.names <- paste('Iteration',as.character(seq(1:5)), sep='')
names(cols_it) <- cols_it.names

## Among datasets
bp_performance.res <- bp_performance(p_metrics = performance_metrics, 
                                     cols_vec = cols_it, 
                                     df = logistic_lasso.out$performance)
