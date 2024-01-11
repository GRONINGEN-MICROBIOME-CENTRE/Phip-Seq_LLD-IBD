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
  make_option(c("--cmv_metadata"), action="store", default="data/Metadata_matchLLD_CMV.tsv", type='character',
              help="CMV metadata"),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--in_dir"), action="store", default='data/CMV_sc_1_CODA/Azimuth_celltype_proportions', type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_3.1_PredictCMV_CellTypeProportions', type='character',
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
shhh(library(viridis))

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

# Predict CMV status
## accessory functions
### report function
report_model <- function(i, l, df, phe){
  print(paste0("Report Fold: ", i))
  
  # Select train/test data
  select_rows <- l[[i]]
  test <- droplevels(df[select_rows,])
  train <- droplevels(df[-select_rows,])
  
  # Report
  train$type <- 'train'
  test$type <- 'test'
  df_to_report <- rbind(train, test)
  group_vars <- c('type', phe)
  df_to_report %>% 
    group_by_at(group_vars, .drop = FALSE) %>%
    summarise(n = n()) %>%
    mutate(prop = n / sum(n)) %>% as.data.frame() -> df.c
  
  df.c$CV <- i
  
  return(df.c)
}

### fit function
# Lasso_n = 4 #Number of folds to split the training data  ( (K-1)/K from total data ) for cross-validation to pick up the best hyperparameter lambda (regularization strength)
fit_model <- function(i, l, df, phe, Lasso_n = 4){
  print(paste0("Fitting Fold: ", i ))
  
  # Select train/test data
  select_rows <- l[[i]]
  test <- droplevels(df[select_rows,])
  train <- droplevels(df[-select_rows,])
  
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
  important_features.df$CV <- i
  important_features.df <- important_features.df[order(-important_features.df$value),]
  
  # Fit the final model on the training data
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
  conf_mat.out <- conf_mat(performances.df, truth = truth, estimate = predicted)
  classification_metrics <- metric_set(accuracy, mcc, f_meas)
  classification_metrics.out <- classification_metrics(performances.df, truth = truth, estimate = predicted)
  class_name <- class_cnames[1]
  roc_auc.out <- roc_auc(performances.df, truth, class_name)
  metrics_df <- rbind(classification_metrics.out, roc_auc.out) %>% as.data.frame()
  metrics_df <- metrics_df[,-2]
  colnames(metrics_df) <- gsub('\\.','',colnames(metrics_df))
  metrics_df$Model <- "Logistic_Lasso"
  metrics_wide <- dcast(metrics_df, Model ~ metric, value.var='estimate')
  metrics_wide$CV <- i
  
  ## Curves
  curve_list <- list(roc_curve = roc_curve(performances.df, truth, class_name),
                     gain_curve = gain_curve(performances.df, truth, class_name),
                     lift_curve = lift_curve(performances.df, truth, class_name),
                     pr_curve = pr_curve(performances.df, truth, class_name))
  
  # Output
  out <- list(performance = metrics_wide,
              features = important_features.df,
              curves = curve_list,
              confusion_matrix = conf_mat.out)
  return(out)
}

## main function
logistic_lasso.func <- function(ds, N, K = 3, l, phenotype){
  print(paste0('######### Dataset: ', ds, ' #########'))
  print(paste0('# Iteration (N): ', as.character(N)))
  print(paste0('# Folds (K): ', as.character(K)))
  
  # Data
  df_i <- droplevels(l[[ds]])
  
  # Split data in folds
  row_folds <- caret::createFolds(1:nrow(df_i), k=K) 
  
  # Apply report on the different folds
  report_model.res <- sapply(names(row_folds), function(fold) report_model(i = fold,
                                                                           l = row_folds,
                                                                           df = df_i,
                                                                           phe = phenotype), simplify = FALSE)
  report_model.df <- do.call("rbind", report_model.res)
  report_model.df$Iteration <- paste0('Iteration',as.character(N))
  cat('\n')
  
  # Apply fit model on the different folds
  fit_model.res <- sapply(names(row_folds), function(fold) fit_model(i = fold,
                                                                     l = row_folds,
                                                                     df = df_i,
                                                                     phe = phenotype), simplify = FALSE)
  out_names <- c('performance', 'features', 'curves')
  fit_model.out <- sapply(out_names, function(i) 
    lapply(fit_model.res, function(fold) fold[[i]]), simplify = FALSE)
  fit_model.out.p_f <- fit_model.out[names(fit_model.out)%in%c('performance', 'features')]
  fit_model.list <- lapply(fit_model.out.p_f, function(x){
    xx <- do.call("rbind", x)
    xx$Iteration <- paste0('Iteration',as.character(N))
    return(xx)
  })
  
  # Save curves per iteration
  fit_model.out.c <- fit_model.out[names(fit_model.out)%in%c('curves')]
  it_tag <- paste0('Iteration',as.character(N))
  curves_model.list <- list()
  curves_model.list[[it_tag]] <- fit_model.out.c
  curves_model.fn <- paste0(out.dir, ds, '.', it_tag, '.curves.rds')
  saveRDS(curves_model.list, curves_model.fn)
  
  # Save conf matrix per iteration
  fit_model.out.cm <- fit_model.out[names(fit_model.out)%in%c('confusion_matrix')]
  it_tag <- paste0('Iteration',as.character(N))
  cm_model.list <- list()
  cm_model.list[[it_tag]] <- fit_model.out.cm
  cm_model.fn <- paste0(out.dir, ds, '.', it_tag, '.confusion_matrix.rds')
  saveRDS(cm_model.list, cm_model.fn)
  
  # Output
  report_model.list <- list(report = report_model.df)
  out <- c(fit_model.list, report_model.list)
  
  cat('\n')
  cat('\n')
  return(out)
}

# Boxplots (performance metrics)
## main function (per dataset)
bp_performance_ds <- function(ds, p_metrics, cols_vec, size_vec, alpha_vec, df, width_var = 6.5){
  df_i <- droplevels(df[df$Dataset==ds,])
  df_melt <- melt(df_i, measure.vars = names(p_metrics))
  df_melt <- droplevels(df_melt[df_melt$variable%in%names(p_metrics),])
  df_melt$Type <- 'Individual'
  df_melt %>% 
    group_by(variable, Iteration) %>%
    summarise(mean=mean(value),
              median=median(value)) %>% as.data.frame() -> df_melt.stats
  df_melt.stats$Type <- 'Average'
  df_stats <- df_melt.stats[,c('variable', 'Iteration', 'median', 'Type')]
  colnames(df_stats)[colnames(df_stats)=='median'] <- 'value' 
  df_tmp <- df_melt[,c('variable', 'Iteration', 'value', 'Type')]
  df_all <- rbind(df_stats, df_tmp)
  df_all$Iteration <- factor(df_all$Iteration,
                             levels = names(cols_vec))
  df_all$Metric <- unname(p_metrics[df_all$variable])
  df_all$Metric <- factor(df_all$Metric,
                          levels = unname(p_metrics))
  p <- ggplot(df_all, aes(x=Metric, y=value, group=Metric)) +
    geom_boxplot(outlier.alpha = 0, fatten = 1) +
    geom_jitter(stroke=NA, width = 0.2, aes(color = Iteration, alpha = Type, size = Type)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    theme_bw() +
    ggtitle(paste0(ds)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) + 
    scale_color_manual(values = cols_vec) + 
    scale_alpha_manual(values = alpha_vec,
                       guide = "none") + 
    scale_size_manual(values = size_vec,
                      guide = "none")
  metrics_tag <- paste(names(p_metrics),collapse='.')
  p.fn <- paste0(out.dir, ds, '_', metrics_tag, '.png')
  ggsave(p.fn, p, width = width_var, height = 4)
  return(p)
}

## main function (among datasets)
bp_performance <- function(p_metrics, cols_vec, size_vec, alpha_vec, df, width_var = 6){
  df_melt <- melt(df, measure.vars = names(p_metrics))
  df_melt <- droplevels(df_melt[df_melt$variable%in%names(p_metrics),])
  df_melt$Type <- 'Individual'
  df_melt %>% 
    group_by(Dataset, variable, Iteration) %>%
    summarise(mean=mean(value),
              median=median(value)) %>% as.data.frame() -> df_melt.stats
  df_melt.stats$Type <- 'Average'
  df_stats <- df_melt.stats[,c('Dataset', 'variable', 'Iteration', 'median', 'Type')]
  colnames(df_stats)[colnames(df_stats)=='median'] <- 'value' 
  df_tmp <- df_melt[,c('Dataset', 'variable', 'Iteration', 'value', 'Type')]
  df_all <- rbind(df_stats, df_tmp)
  df_all$Iteration <- factor(df_all$Iteration,
                             levels = names(cols_vec))
  df_all$Metric <- unname(p_metrics[df_all$variable])
  df_all$Metric <- factor(df_all$Metric,
                          levels = unname(p_metrics))
  df_all$Dataset <- factor(df_all$Dataset,
                           levels = datasets)
  p <- ggplot(df_all, aes(x=Dataset, y=value)) +
    geom_boxplot(outlier.alpha = 0, fatten = 1) +
    geom_jitter(stroke=NA, width = 0.2, aes(color = Iteration, alpha = Type, size = Type)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    theme_bw() +
    ylab('Performance metric value') + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) + 
    scale_color_manual(values = cols_vec) + 
    scale_alpha_manual(values = alpha_vec,
                       guide = "none") + 
    scale_size_manual(values = size_vec,
                      guide = "none") + 
    facet_grid(. ~ Metric,
               scales = 'free')
  metrics_tag <- paste(names(p_metrics),collapse='.')
  p.fn <- paste0(out.dir, metrics_tag, '.png')
  ggsave(p.fn, p, width = width_var, height = 3.5)
  return(p)
}

# Barplots (feature importance)
bp_features <- function(ds, bio_covs, tech_covs, out_covs, df){
  df_i <- droplevels(df[df$Dataset==ds & df$feature!='(Intercept)',])
  bio_tech_out_covs <- c(bio_covs, tech_covs, out_covs)
  ct.covs <- setdiff(colnames(df_model), bio_tech_out_covs)
  bio_covs.regex <- paste(bio_covs,collapse='|')
  tech_covs.regex <- tech_covs
  features <- unique(df_i$feature)
  bio.features <- features[grepl(bio_covs.regex, features)]
  tech.features <- features[grepl(tech_covs.regex, features)]
  bio_tech.features <- c(bio.features, tech.features)
  ct.features <- setdiff(ct.covs, bio_tech.features)
  df_i$Type <- ifelse(df_i$feature%in%bio.features, 'Biological', 
                      ifelse(df_i$feature%in%tech.features, 'Technical', 'CellTypes'))
  it.order <- paste('Iteration',as.character(seq(1:5)), sep='')
  
  # Order
  df_i$Iteration <- factor(df_i$Iteration,
                           levels = it.order)
  df_i$Type <- factor(df_i$Type,
                      levels = c('CellTypes', 'Biological', 'Technical'))
  df_i %>% 
    group_by(feature) %>% 
    count() %>% 
    arrange(desc(n)) %>% as.data.frame() -> df_i.c
  df_i.tmp <- unique(df_i[,c('feature','Type')])
  df_i.c <- merge(df_i.c, df_i.tmp)
  df_i.c <- df_i.c[order(df_i.c$Type, -df_i.c$n),]
  feature.order <- unique(df_i.c$feature)
  df_i$feature <- factor(df_i$feature,
                         levels = feature.order)
  
  p <- ggplot(df_i, aes(x = feature, y = n, fill = median_value)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Iteration ~ Type, scales = "free_x", space = "free") +
    scale_fill_viridis() + 
    theme_bw() +
    ggtitle(paste0(ds)) +
    ylab('N') +
    xlab(NULL) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text = element_text(hjust = 0.5, face="bold", size = 10),
          panel.spacing = unit(1, "lines"),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank())
  p.fn <- paste0(out.dir, ds, '_features.png')
  ggsave(p.fn, p, width = 10, height = 7)
  return(p)
}

#################### Set Variables and load Data #################### 
# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', datasets_tag, '/', opt$cell_level, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Datasets: ', datasets_tag))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

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
cell.md_fn <- paste0(in.dir, '/metadata.rds')
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
props_clr.tmp <- as.data.frame(t(props_clr.df))
props_clr.tmp$donor <- rownames(props_clr.tmp)
rownames(props_clr.tmp) <- NULL
df_model <- merge(donor_md, props_clr.tmp, by = 'donor')
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
# Split data by chemistry (dataset)
df_model.byDS <- split(df_model, df_model$dataset)
df_model.byDS <- lapply(df_model.byDS, function(x) x[,-which(colnames(x)=='dataset')])
N.vec <- seq(1:5)

# Apply function
set.seed(123)
logistic_lasso.res <- sapply(names(df_model.byDS), function(dataset)
  sapply(N.vec, function(iteration) logistic_lasso.func(ds = dataset,
                                                        N = iteration,
                                                        l = df_model.byDS, 
                                                        phenotype = opt$phenotype), simplify = FALSE), simplify = FALSE)
  
# Retrieve outputs
out_names <- c('performance', 'features', 'report')
logistic_lasso.out <- sapply(out_names, function(i) 
  lapply(logistic_lasso.res, function(dataset)
    lapply(dataset, function(iteration) iteration[[i]])), simplify = FALSE)

logistic_lasso.byDS <- lapply(logistic_lasso.out, function(out_name) 
  lapply(out_name, function(dataset) do.call("rbind",dataset)))
logistic_lasso <-  lapply(logistic_lasso.byDS, function(out_name) do.call("rbind", out_name))
logistic_lasso <- lapply(logistic_lasso, function(x){
  x$Dataset <- str_split_fixed(rownames(x),'\\.',2)[,1]
  rownames(x) <- NULL
  return(x)})
logistic_lasso.performance <- logistic_lasso$performance
logistic_lasso.features <- logistic_lasso$features
logistic_lasso.report <- logistic_lasso$report
logistic_lasso.fn <- paste0(out.dir, 'logistic_lasso.rds')
print(paste0('Saving logistic lasso results: ', logistic_lasso.fn))
saveRDS(logistic_lasso, logistic_lasso.fn)

# Outputs' stats
## performance
group_vars <- c('Dataset', 'Iteration')
logistic_lasso.performance %>% 
  group_by_at(group_vars, .drop = FALSE) %>%
  summarize(accuracy.mean = median(accuracy),
            f_meas.mean = median(f_meas),
            mcc.mean = median(mcc),
            roc_auc.mean = median(roc_auc)) %>%
  as.data.frame() -> performance.summary

## features
group_vars <- c('Dataset', 'Iteration', 'feature')
logistic_lasso.features %>% 
  group_by_at(group_vars, .drop = FALSE) %>%
  summarize(n = n(),
            mean_value = mean(value),
            median_value = median(value)) %>%
  arrange(Dataset, Iteration, desc(n), desc(abs(median_value))) %>%
  as.data.frame() -> features.summary

## report
group_vars <- c('Dataset', 'type', 'Iteration', 'CMV_status')
logistic_lasso.report %>% 
  group_by_at(group_vars, .drop = FALSE) %>%
  summarize(N_mean = mean(n),
            N_min = min(n),
            N_max = max(n)) %>% as.data.frame() -> report.summary

summary_list <- list(performance = performance.summary,
                     features = features.summary,
                     report = report.summary)
summary_list.fn <- paste0(out.dir, 'summary_list.rds')
print(paste0('Saving logistic lasso (summary) results: ', summary_list.fn))
saveRDS(summary_list, summary_list.fn)

#################### Plots: boxplots (performance metrics) + barplots (feature importance) ####################
# Boxplots (performance metrics)
## Variables
performance_metrics <- c('ROC AUC', 'F1', 'MCC', 'Accuracy')
performance_metrics.names <- c('roc_auc', 'f_meas', 'mcc', 'accuracy')
names(performance_metrics) <- performance_metrics.names
# performance_metrics.names_selected <- c('accuracy', 'roc_auc')
# performance_metrics.selected <- performance_metrics[names(performance_metrics)%in%performance_metrics.names_selected]

cols_it <- c('#264653', '#2a9d8f', '#e9c46a', '#f4a261', '#e76f51')
cols_it.names <- paste('Iteration',as.character(seq(1:5)), sep='')
names(cols_it) <- cols_it.names
size_type <- c(1.5, 3.5)
names(size_type) <- c('Individual', 'Average')
alpha_type <- c(0.6, 1)
names(alpha_type) <- c('Individual', 'Average')

## Apply functions
### Per dataset
bp_performance_ds.res <- sapply(datasets, function(i) bp_performance_ds(ds = i,
                                                                        p_metrics = performance_metrics, 
                                                                        cols_vec = cols_it, 
                                                                        size_vec = size_type, 
                                                                        alpha_vec = alpha_type, 
                                                                        df = logistic_lasso.performance), simplify = FALSE)

### Among datasets
bp_performance.res <- bp_performance(p_metrics = performance_metrics, 
                                     cols_vec = cols_it, 
                                     size_vec = size_type, 
                                     alpha_vec = alpha_type, 
                                     df = logistic_lasso.performance)


# Barplots (feature importance)
## Variables
bio.covs <- c('Sex','Age')
tech.covs <- 'date'
out.covs <- c('CMV_status', 'dataset')
cols_it <- c('#264653', '#2a9d8f', '#e9c46a', '#f4a261', '#e76f51')
cols_it.names <- paste('Iteration',as.character(seq(1:5)), sep='')

## Apply function
bp_features.res <- sapply(datasets, function(i) bp_features(ds = i,
                                                            bio_covs = bio.covs, 
                                                            tech_covs = tech.covs, 
                                                            out_covs = out.covs, 
                                                            df = features.summary), simplify = FALSE)