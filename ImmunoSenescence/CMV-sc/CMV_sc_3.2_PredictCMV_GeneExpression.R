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
  make_option(c("--cell_type"), action="store", default=NULL, type='character',
              help="predicted.celltype.l1 or cell_type"),
  make_option(c("--cell_level"), action="store", default="cell_type", type='character',
              help="predicted.celltype.l1 or cell_type"),
  make_option(c("--phenotype"), action="store", default="CMV_status", type='character',
              help="predicted.celltype.l1 or cell_type"),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--random"), action="store", default="date", type='character',
              help="Date, lane or any variables."),
  make_option(c("--pb_dir"), action="store", default='output/CMV_sc_2.1_DEA_Pseudobulk', type='character',
              help="Pseudobulk directory"),
  make_option(c("--md_dir"), action="store", default='output/CMV_sc_2.2_DEA_LimmaDream', type='character',
              help="Metadata directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_3.2_PredictCMV_GeneExpression', type='character',
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
shhh(library(edgeR))
shhh(library(yardstick))
shhh(library(viridis))
shhh(library(ggpubr))
shhh(library(ggrepel))

#################### Define functions #################### 
# Read input files
## main function
read_data <- function(dataset, pseudobulk_dir, metadata_dir){
  print(dataset)
  
  # define input fns
  pseudobulk_fn <- paste0(pseudobulk_dir, dataset, '/', opt$cell_level, '/', opt$cell_type, '/pseudobulk_matrices.rds')
  metadata_fn <- grep('vars.rds', list.files(metadata_dir, full.names = TRUE), value = TRUE)
  
  # counts and metadata
  print(paste0('Reading pseudobulk: ', pseudobulk_fn))
  counts <- readRDS(pseudobulk_fn)$geneExpr_LimmaVoom$counts
  print(paste0('Reading phenotype metadata: ', metadata_fn))
  metadata <- readRDS(metadata_fn)$metadata
  metadata <- droplevels(metadata[metadata$dataset==dataset,])
  cnames <- c('LLDEEP_ID', 'Gender', 'Age', 'date')
  if(opt$phenotype=='CMV_status'){cnames <- c(cnames, opt$phenotype)}
  metadata_df <- metadata[,cnames]
  colnames(metadata_df)[colnames(metadata_df)=='Gender'] <- 'Sex'
  colnames(metadata_df)[colnames(metadata_df)=='LLDEEP_ID'] <- 'donor'
  counts <- counts[,colnames(counts)%in%rownames(metadata)]
  if(!identical(colnames(counts), rownames(metadata))){
    counts <- counts[,match(rownames(metadata), colnames(counts))]
  }
  counts_df <- as.data.frame(t(counts))
  counts_df$donor <- rownames(counts_df)
  
  # merge data
  df_model <- merge(metadata_df, counts_df, by = 'donor')
  df_model <- droplevels(df_model)
  df_model$date <- as.factor(df_model$date)
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
  return(df_model)
}

# Predict CMV status
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
  system.time(model_lasso <- glmnetUtils::cv.glmnet(form, 
                                                    train_lasso, 
                                                    alpha = 1,
                                                    nfolds = Lasso_n,
                                                    family = "binomial", 
                                                    type.measure="class"))
  
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
  print('Predicting the phenotype...')
  test_lasso.newdata <- test_lasso %>% select(-all_of(phe)) %>% as.data.frame()
  
  system.time(lasso_probabilities <- predict(model_lasso,
                                             newdata = test_lasso.newdata,
                                             type="response",
                                             s="lambda.min"))
  
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
  out_names <- c('performance', 'features', 'curves', 'confusion_matrix')
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
## main function
bp_performance <- function(ds, p_metrics, cols_vec, size_vec, alpha_vec, df, width_var = 6.5){
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

# Barplots (feature importance)
## main function
bp_features <- function(ds, logFC_th = 0, df){
  # pick DEGs
  phe_tag <- ifelse(opt$phenotype!='Age', 
                    paste0(opt$phenotype, '0_', opt$phenotype, '1'), opt$phenotype)
  dea_fn <- paste0(md.dir, phe_tag, '.rds')
  dea_df <- readRDS(dea_fn)[[phe_tag]]
  degs <- dea_df[dea_df$adj.P.Val<=0.05 & abs(dea_df$logFC)>=logFC_th,]$gene
  logFC_th_tag <- paste0('logFC', as.character(logFC_th))
  
  # annotate features
  df_i <- droplevels(df[df$Dataset==ds & df$feature!='(Intercept)',])
  df_i$Type <- ifelse(df_i$feature%in%degs, 'DEGs', 'non_DEGs')
  it.order <- paste('Iteration',as.character(seq(1:5)), sep='')
  
  # Order
  df_i$Iteration <- factor(df_i$Iteration,
                           levels = it.order)
  df_i$Type <- factor(df_i$Type,
                      levels = c('DEGs', 'non_DEGs'))
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
    ggtitle(paste0(ds, ' (', logFC_th_tag, ')')) +
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
  p.fn <- paste0(out.dir, ds, '_', logFC_th_tag, '_features.png')
  ggsave(p.fn, p, width = 14, height = 7)
  return(p)
}

# Scatter plot (feature importance from prediction vs. logFC from DEA)
## main function
sp_features <- function(ds, metric, logFC_th = 0, df, cols_vec){
  print(ds)
  print(metric)
  
  # pick DEGs
  phe_tag <- ifelse(opt$phenotype!='Age', 
                    paste0(opt$phenotype, '0_', opt$phenotype, '1'), opt$phenotype)
  dea_fn <- paste0(md.dir, phe_tag, '.rds')
  dea_df <- readRDS(dea_fn)[[phe_tag]]
  degs <- dea_df[dea_df$adj.P.Val<=0.05 & abs(dea_df$logFC)>=logFC_th,]$gene
  logFC_th_tag <- paste0('logFC', as.character(logFC_th))
  dea_df$Type <- ifelse(dea_df$gene%in%degs, 'DEGs', 'non_DEGs')
  
  # annotate features
  df_i <- droplevels(df[df$Dataset==ds & df$feature!='(Intercept)',])
  df_i %>% 
    group_by(feature, .drop = FALSE) %>%
    summarize(mean = mean(mean_value),
              median = median(median_value)) %>% as.data.frame() -> df_i
  df_i <- merge(df_i, dea_df, by.x = 'feature', by.y = 'gene')
  it.order <- paste('Iteration',as.character(seq(1:5)), sep='')
  
  # Order
  df_i$Type <- factor(df_i$Type,
                      levels = c('DEGs', 'non_DEGs'))
  df_i$log10pval <- -log10(df_i$adj.P.Val)
  df_i$abslogFC <- abs(df_i$logFC)
  # min_x <- ifelse(metric=='logFC', -1.5, 0)
  min_x <- 0
  title_var <- paste0(ds, ' > ', logFC_th_tag, ' (', metric, ')')
  metric_th <- ifelse(metric=='logFC', -1, 3)
  feature_th <- 1
  
  # logFC metric
  if(metric=='logFC'){
    p <- ggplot(df_i, aes(x = .data[[metric]], y = median, color = Type, fill = Type)) +
      geom_point()+
      geom_smooth(method=lm)+
      geom_hline(yintercept=0, linetype='dashed') +
      geom_vline(xintercept=0, linetype='dashed') +
      stat_cor(method = "spearman", aes(color = Type), label.x = -0.5)+
      scale_color_manual(values = cols_vec) + 
      scale_fill_manual(values = cols_vec) + 
      theme_bw() +
      ylab('Feature importance (median)') + 
      ggtitle(title_var) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust=0.5, face="bold"),
            strip.text = element_text(hjust = 0.5, face="bold", size = 10),
            panel.spacing = unit(1, "lines"),
            axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size = 10),
            axis.text.y = element_text(size = 10)) + 
      geom_text_repel(data = df_i[df_i$median>=feature_th & df_i$logFC<=metric_th,], 
                      aes(x=logFC, y=median, label=feature))
  }
  
  # log10pval
  if(metric=='log10pval'){
    p <- ggplot(df_i, aes(x = .data[[metric]], y = median)) +
      geom_point(aes(color = Type))+
      geom_smooth(method=lm)+
      geom_vline(xintercept=-log10(0.05), linetype='dashed') +
      stat_cor(method = "spearman", label.x = min_x)+
      scale_color_manual(values = cols_vec) + 
      scale_fill_manual(values = cols_vec) + 
      theme_bw() +
      ylab('Feature importance (median)') + 
      ggtitle(title_var) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust=0.5, face="bold"),
            strip.text = element_text(hjust = 0.5, face="bold", size = 10),
            panel.spacing = unit(1, "lines"),
            axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size = 10),
            axis.text.y = element_text(size = 10)) + 
      geom_text_repel(data = df_i[df_i$median>=feature_th & df_i$log10pval>=metric_th,], 
                      aes(x=log10pval, y=median, label=feature))
  }
  p.fn <- paste0(out.dir, ds, '_', logFC_th_tag, '_features.', metric ,'.png')
  ggsave(p.fn, p, width = 4.5, height = 4)
  print(p.fn)
  return(p)
}

#################### Set Variables and load Data #################### 
# Cell type
# opt$cell_type <- 'CD8_TEM'
# opt$cell_type <- 'CD4_CTL'

# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Input/Output directories
pb.dir <- paste0(opt$pb_dir, '/')
md.dir <- paste0(opt$md_dir, '/', datasets_tag, '/',  
                 opt$cell_level, '/', opt$cell_type, '/', 
                 opt$random, '/', opt$phenotype, '/')
out.dir <- paste0(opt$out_dir, '/', datasets_tag, '/',
                  opt$phenotype, '/', opt$cell_level, '/', opt$cell_type, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$type))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Datasets: ', datasets_tag))
print(paste0('Pseudobulk dir: ', pb.dir))
print(paste0('Metadata dir: ', md.dir))
print(paste0('Output dir: ', out.dir))
print('############################')
cat('\n')

#################### Predict CMV status ####################
# Read data by chemistry (dataset)
df_model.byDS <- sapply(datasets, function(i) read_data(dataset = i,
                                                        pseudobulk_dir = pb.dir, 
                                                        metadata_dir = md.dir), simplify = FALSE)
# df_model.byDS <- lapply(df_model.byDS, function(x) x[,1:1000])
N.vec <- seq(1:5)

# Apply function
set.seed(222)
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
  summarize(accuracy.median = median(accuracy),
            f_meas.median = median(f_meas),
            mcc.median = median(mcc),
            roc_auc.median = median(roc_auc)) %>%
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

#################### Plots: boxplots (performance metrics) + barplots (feature importance) + scatter plots (feature importance from prediction vs. logFC from DEA) ####################
# Boxplots (performance metrics)
## Variables
performance_metrics <- c('Accuracy', 'ROC AUC', 'F1', 'MCC')
performance_metrics.names <- c('accuracy', 'roc_auc', 'f_meas', 'mcc')
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

## Apply function
bp_performance.res <- sapply(datasets, function(i) bp_performance(ds = i,
                                                                  p_metrics = performance_metrics, 
                                                                  cols_vec = cols_it, 
                                                                  size_vec = size_type, 
                                                                  alpha_vec = alpha_type, 
                                                                  df = logistic_lasso.performance), simplify = FALSE)

# Barplots (feature importance - DEGs)
## Variables
cols_it <- c('#264653', '#2a9d8f', '#e9c46a', '#f4a261', '#e76f51')
cols_it.names <- paste('Iteration',as.character(seq(1:5)), sep='')

## Apply function
bp_features.res <- sapply(datasets, function(i) bp_features(ds = i,
                                                            df = features.summary), simplify = FALSE)

# Scatter plot (feature importance - DEGs)
## Variables
cols_type <- c('#264653', '#8FBACC')
names(cols_type) <- c('DEGs', 'non_DEGs') 

## Apply function
metrics <- c('logFC', 'log10pval')
sp_features.res <- sapply(datasets, function(i) 
  sapply(metrics, function(j) sp_features(ds = i,
                                          metric = j,
                                          df = features.summary, 
                                          cols_vec = cols_type), simplify = FALSE), simplify = FALSE)