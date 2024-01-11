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
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--cell_types"), action="store", default="data/CMV_sc_3.2_PredictCMV_GeneExpression_Plots/cell_types.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--in_dir"), action="store", default='output/CMV_sc_3.2_PredictCMV_GeneExpression', type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_3.2_PredictCMV_GeneExpression_Plots', type='character',
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

#################### Define functions #################### 
# Read input files
read_data <- function(in_dir){
  cts <- gsub('/.+','',list.files(in_dir, recursive=TRUE, pattern='logistic_lasso.rds'))
  performance.byCT <- sapply(cts, function(i) readRDS(paste0(in_dir, i, '/logistic_lasso.rds'))$performance, simplify = FALSE)
  performance.df <- do.call("rbind", performance.byCT)
  performance.df$cell_type <- str_split_fixed(rownames(performance.df), '\\.', 2)[,1]
  return(performance.df)
}

# Boxplots (performance metrics)
## main function (per dataset)
bp_performance_ds <- function(ds, p_metrics, cols_vec, size_vec, alpha_vec, cts_order, df, width_var = 10){
  df_i <- droplevels(df[df$Dataset==ds,])
  df_melt <- melt(df_i, measure.vars = names(p_metrics))
  df_melt <- droplevels(df_melt[df_melt$variable%in%names(p_metrics),])
  df_melt$Type <- 'Individual'
  df_melt %>% 
    group_by(variable, Iteration, cell_type) %>%
    summarise(mean=mean(value, na.rm=TRUE),
              median=median(value, na.rm=TRUE)) %>% as.data.frame() -> df_melt.stats
  df_melt.stats$Type <- 'Average'
  df_stats <- df_melt.stats[,c('variable', 'Iteration', 'cell_type', 'median', 'Type')]
  colnames(df_stats)[colnames(df_stats)=='median'] <- 'value' 
  df_tmp <- df_melt[,c('variable', 'Iteration', 'cell_type', 'value', 'Type')]
  df_all <- rbind(df_stats, df_tmp)
  df_all <- droplevels(df_all[!is.na(df_all$value),])
  df_all <- merge(df_all, cts_df, by = 'cell_type', all.x = TRUE)
  df_all$direction <- ifelse(is.na(df_all$direction), 'non significant', df_all$direction)
  df_all$Iteration <- factor(df_all$Iteration,
                             levels = names(cols_vec))
  df_all$direction <- factor(df_all$direction,
                             levels = c('increase', 'decrease', 'non significant'))
  df_all$Metric <- unname(p_metrics[df_all$variable])
  df_all$Metric <- factor(df_all$Metric,
                          levels = unname(p_metrics))
  cts_in <- cts_order[cts_order%in%unique(df_all$cell_type)]
  df_all$cell_type <- factor(df_all$cell_type,
                             levels = cts_in)
  p <- ggplot(df_all, aes(x=cell_type, y=value)) +
    geom_boxplot(outlier.alpha = 0, fatten = 1) +
    geom_jitter(stroke=NA, width = 0.2, aes(color = Iteration, alpha = Type, size = Type)) +
    theme_bw() +
    ggtitle(paste0(ds)) +
    facet_grid(Metric ~ direction, scales = "free", space = "free_x") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text.y = element_text(angle=0, hjust=0.5, face="bold", size = 10),
          strip.text.x = element_text(hjust=0.5, face="bold", size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 11),
          axis.text.y = element_text(size = 11),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          panel.spacing = unit(0.5, "lines")) + 
    scale_color_manual(values = cols_vec) + 
    scale_alpha_manual(values = alpha_vec,
                       guide = "none") + 
    scale_size_manual(values = size_vec,
                      guide = "none")
  metrics_tag <- paste(names(p_metrics),collapse='.')
  p.fn <- paste0(out.dir, ds, '_', metrics_tag, '.png')
  ggsave(p.fn, p, width = width_var, height = 7)
  return(p)
}

## main function (among datasets)
bp_performance <- function(cts = c('CD4_CTL', 'CD8_TEM'), p_metrics, cols_vec, size_vec, alpha_vec, cts_order, df, width_var = 6){
  df_i <- droplevels(df[df$cell_type%in%cts,])
  df_melt <- melt(df_i, measure.vars = names(p_metrics))
  df_melt <- droplevels(df_melt[df_melt$variable%in%names(p_metrics),])
  df_melt$Type <- 'Individual'
  df_melt %>% 
    group_by(Dataset, variable, Iteration, cell_type) %>%
    summarise(mean=mean(value, na.rm=TRUE),
              median=median(value, na.rm=TRUE)) %>% as.data.frame() -> df_melt.stats
  df_melt.stats$Type <- 'Average'
  df_stats <- df_melt.stats[,c('Dataset', 'variable', 'Iteration', 'cell_type', 'median', 'Type')]
  colnames(df_stats)[colnames(df_stats)=='median'] <- 'value' 
  df_tmp <- df_melt[,c('Dataset', 'variable', 'Iteration', 'cell_type', 'value', 'Type')]
  df_all <- rbind(df_stats, df_tmp)
  df_all <- droplevels(df_all[!is.na(df_all$value),])
  df_all <- merge(df_all, cts_df, by = 'cell_type', all.x = TRUE)
  df_all$direction <- ifelse(is.na(df_all$direction), 'non significant', df_all$direction)
  df_all$Iteration <- factor(df_all$Iteration,
                             levels = names(cols_vec))
  df_all$direction <- factor(df_all$direction,
                             levels = c('increase', 'decrease', 'non significant'))
  df_all$Metric <- unname(p_metrics[df_all$variable])
  df_all$Metric <- factor(df_all$Metric,
                          levels = unname(p_metrics))
  df_all$cell_type <- factor(df_all$cell_type,
                             levels = cts)
  p <- ggplot(df_all, aes(x=cell_type, y=value)) +
    geom_boxplot(outlier.alpha = 0, fatten = 1) +
    geom_jitter(stroke=NA, width = 0.2, aes(color = Iteration, alpha = Type, size = Type)) +
    theme_bw() +
    ylab('Performance metric value') + 
    facet_grid(Dataset ~ Metric, scales = "free") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text.y = element_text(angle=0, hjust=0.5, face="bold", size = 10),
          strip.text.x = element_text(hjust=0.5, face="bold", size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 11),
          axis.text.y = element_text(size = 11),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          panel.spacing = unit(0.5, "lines")) + 
    scale_color_manual(values = cols_vec) + 
    scale_alpha_manual(values = alpha_vec,
                       guide = "none") + 
    scale_size_manual(values = size_vec,
                      guide = "none")
  metrics_tag <- paste(names(p_metrics),collapse='.')
  p.fn <- paste0(out.dir, metrics_tag, '.png')
  ggsave(p.fn, p, width = width_var, height = 4.5)
  return(p)
}

#################### Set Variables and load Data #################### 
# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Input/Output directories
suffix <- paste0(datasets_tag, '/', opt$phenotype, '/', opt$cell_level, '/')

## Input directory
in.dir <- paste0(opt$in_dir, '/', suffix)

## Output directory
out.dir <- paste0(opt$out_dir, '/', suffix)
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Cell types
cts_fn <- opt$cell_types
cts_df <- read.table(cts_fn, header=TRUE)

# Report
print('############################')
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Datasets: ', datasets_tag))
print(paste0('Input dir: ', in.dir))
print(paste0('Output dir: ', out.dir))
print('############################')
cat('\n')

# Read input files
logistic_lasso.performance <- read_data(in.dir)

#################### Plots: boxplots (performance metrics) + barplots (feature importance) ####################
# Boxplots (performance metrics)
## Variables
### summary (not needed)
group_vars <- c('Dataset', 'Iteration', 'cell_type')
logistic_lasso.performance %>% 
  group_by_at(group_vars, .drop = FALSE) %>%
  summarize(accuracy.median = median(accuracy),
            f_meas.median = median(f_meas),
            mcc.median = median(mcc),
            roc_auc.median = median(roc_auc)) %>%
  as.data.frame() -> performance.summary

performance_metrics <- c('ROC AUC', 'F1', 'MCC', 'Accuracy')
performance_metrics.names <- c('roc_auc', 'f_meas', 'mcc', 'accuracy')
names(performance_metrics) <- performance_metrics.names
# performance_metrics.names_selected <- c('accuracy', 'roc_auc')
# performance_metrics.selected <- performance_metrics[names(performance_metrics)%in%performance_metrics.names_selected]

### features
cols_it <- c('#264653', '#2a9d8f', '#e9c46a', '#f4a261', '#e76f51')
cols_it.names <- paste('Iteration',as.character(seq(1:5)), sep='')
names(cols_it) <- cols_it.names
size_type <- c(0.75, 1.75)
names(size_type) <- c('Individual', 'Average')
alpha_type <- c(0.6, 1)
names(alpha_type) <- c('Individual', 'Average')

ct_fn <- 'data/Azimuth_l2.cell_type.tab'
ct.order <- read.table(ct_fn)$V1

## Per dataset
bp_performance_ds.res <- sapply(datasets, function(i) bp_performance_ds(ds = i,
                                                                        p_metrics = performance_metrics, 
                                                                        cols_vec = cols_it, 
                                                                        size_vec = size_type, 
                                                                        alpha_vec = alpha_type, 
                                                                        cts_order = ct.order, 
                                                                        df = logistic_lasso.performance), simplify = FALSE)

## Among datasets
bp_performance.res <- bp_performance(p_metrics = performance_metrics, 
                                     cols_vec = cols_it, 
                                     size_vec = size_type, 
                                     alpha_vec = alpha_type, 
                                     cts_order = ct.order, 
                                     df = logistic_lasso.performance)
