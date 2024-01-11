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
  make_option(c("--cell_types"), action="store", default='CD4_CTL,CD8_TEM', type='character',
              help="Cell types."),
  make_option(c("--logFC_th"), action="store", default=0.5, type='double',
              help="logFC threshold."),
  make_option(c("--phenotype"), action="store", default='CMV_status', type='character',
              help="Cell level."),
  make_option(c("--cell_level"), action="store", default="cell_type", type='character',
              help="Azimuth's cell level."),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--gene"), action="store", default="B3GAT1", type='character',
              help="Check gene expression in the UMAP."),
  make_option(c("--assay"), action="store", default='SCT', type='character',
              help="RNA/SCT"),
  make_option(c("--dea_dir"), action="store", default="output/CMV_sc_2.2_DEA_LimmaDream", type='character',
              help="DEA main directory"),
  make_option(c("--so_dir"), action="store", default='data/CMV_sc_2.4_DEA_AddModuleScore/CMV_UMAP_Azimuth_by_dataset/', type='character',
              help="Seurat object directory"),
  make_option(c("--out_dir"), action="store", default="output/CMV_sc_2.5_DEA_AddModuleScore", type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages and functions
print('Loading R packages...')
shhh(library(Seurat))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(RColorBrewer))
shhh(library(ggplot2))
shhh(library(ggpubr))
shhh(library(ggrepel))
shhh(library(ggpubr))
shhh(library(tidyverse))
shhh(library(patchwork))
shhh(library(viridis))
shhh(library(Seurat))
shhh(library(scCustomize))
shhh(library(qs))

#################### Define functions #################### 
# AddModuleScore + UMAPs
## main function
AddModuleScore_UMAP.byDS <- function(ct, ds, cts = NULL, phe, assay_so, logFC_th, gene, out_dir){
  print(ct)
  print(ds)
  print(cts)
  
  # Phenotype
  phe_tag <- phe
  if(phe%in%c('Gender', 'CMV_status')){
    if(phe=='Gender'){phe_tag <- 'GenderM_GenderF'}
    if(phe=='CMV_status'){phe_tag <- 'CMV_status0_CMV_status1'}
  }
  
  # Cell type
  ct_i <- gsub('_', ' ', ct)
  
  # Read Seurat object
  so_fn <- paste0(so.dir, ds, '/', assay_so, '.so.rds')
  print(paste0('Reading seurat object: ', so_fn))
  so <- readRDS(so_fn)
  DefaultAssay(so) <- assay_so
  cells <- rownames(so@meta.data[so@meta.data$cell_type==ct_i,])
  pbmc <- so[,cells]
  print(dim(pbmc))
  
  # Read DEA
  if(!is.null(cts)){
    cell_types <- celltypes
    ct_int <- ct
    ct_other <- cell_types[cell_types!=ct_int]
    cts_tag <- paste(c(ct_int, ct_other),collapse='.')
    print(paste0('DEGs from more than one cell type: ', cts_tag))
    out_dir <- paste0(out_dir, cts_tag, '.', cts, '/')
    if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}
    degs_list.byCT <- sapply(cell_types, function(i){
      dea_fn <- paste0(dea.dir, i, '/date/', phe, '/', phe_tag, '.rds')
      print(paste0('Reading DEA object: ', dea_fn))
      dea_df <- readRDS(dea_fn)[[1]]
      degs_df <- dea_df[dea_df$adj.P.Val<=0.05,]
      if(logFC_th!=0){
        degs_df <- degs_df[abs(degs_df$logFC)>logFC_th,]
      }
      genes.up <- degs_df[degs_df$logFC>0,]$gene
      genes.down <- degs_df[degs_df$logFC<0,]$gene
      degs_list <- list(degs_up = genes.up,
                        degs_down = genes.down)
      return(degs_list)
    }, simplify = FALSE)
    vars <- c('degs_up', 'degs_down')
    degs_list.byDIR <- sapply(vars, function(i) lapply(degs_list.byCT, function(x) x[[i]]), simplify = FALSE)
    
    # ModuleScores
    if(cts=='shared'){
      module_up_down.list <- lapply(degs_list.byDIR, function(x) Reduce("intersect",x))
      
    }else if(cts=='unique'){
      module_up_down.list <- lapply(degs_list.byDIR, function(x) setdiff(x[[ct_int]], x[[ct_other]]))
    }
  }else{
    print(paste0('DEGs from one cell type: ', ct))
    dea_fn <- paste0(dea.dir, ct, '/date/', phe, '/', phe_tag, '.rds')
    print(paste0('Reading DEA object: ', dea_fn))
    dea_df <- readRDS(dea_fn)[[1]]
    degs_df <- dea_df[dea_df$adj.P.Val<=0.05,]
    if(logFC_th!=0){
      degs_df <- degs_df[abs(degs_df$logFC)>logFC_th,]
    }
    genes.up <- degs_df[degs_df$logFC>0,]$gene
    genes.down <- degs_df[degs_df$logFC<0,]$gene
    
    # ModuleScores
    module_up_down.list <- list(degs_up = genes.up,
                                degs_down = genes.down)
  }
  
  nbin_var <- 24
  nbin_vars <- nbin_var:1
  for(i in nbin_vars){
    print(i)
    pbmc_check <- try(AddModuleScore(object = pbmc, 
                                     features = module_up_down.list,
                                     nbin = i,
                                     name = 'DEGs'))
    if(class(pbmc_check)!='try-error'){break}
  }
  pbmc <- pbmc_check
  module.vec <- c('Up_DEGs','Down_DEGs')
  colnames(pbmc@meta.data)[colnames(pbmc@meta.data)%in%c('DEGs1','DEGs2')] <- module.vec
  
  ## UMAP DEGs splitted
  degs_splitted.umap <- FeaturePlot(pbmc,
                                    features = module.vec)
  umap.degs_splitted.fn <- paste0(out_dir, ds, '.', ct, '.UMAP.',paste(module.vec, collapse=''),'.png')
  ggsave(umap.degs_splitted.fn, degs_splitted.umap, width = 8, height = 5)
  
  ## UMAP DEGs together --> up_score (DeltaDE score)
  metadata <- pbmc@meta.data
  metadata %>%
    group_by(assignment) %>%
    summarize(mean_up = mean(Up_DEGs, na.rm = TRUE),
              mean_down = mean(Down_DEGs, na.rm = TRUE),
              median_up = median(Up_DEGs, na.rm = TRUE),
              median_down = median(Down_DEGs, na.rm = TRUE)) %>%
    as.data.frame() -> md.stats
  md.stats_melt <- reshape2::melt(md.stats, id.vars='assignment')
  metadata.donor <- droplevels(unique(metadata[,c('assignment', 'CMV_status', 'Age','Gender','lane')]))
  df <- merge(md.stats_melt, metadata.donor, by = 'assignment')
  df$variable <- factor(df$variable, 
                        levels = c('mean_up','mean_down','median_up','median_down'))
  metadata %>%
    group_by(assignment, bare_barcode_lane) %>%
    summarise(ms_total = sum(Up_DEGs,Down_DEGs),
              Up_DEGs = Up_DEGs,
              Down_DEGs = Down_DEGs,
              ms_total.abs = sum(abs(Up_DEGs),abs(Down_DEGs)),
              Up_DEGs.abs = abs(Up_DEGs),
              Down_DEGs.abs = abs(Down_DEGs)) %>%
    mutate(ratio_up = Up_DEGs/ms_total,
           ratio_down = Down_DEGs/ms_total,
           ratio_up.abs = Up_DEGs.abs/ms_total.abs,
           ratio_down.abs = Down_DEGs.abs/ms_total.abs,
           ratio.up_down = Up_DEGs/Down_DEGs,
           up_score = Up_DEGs-Down_DEGs) %>% 
    as.data.frame() -> metadata.ratios
  
  idx <- match(rownames(pbmc@meta.data), metadata.ratios$bare_barcode_lane)
  metadata.ratios <- metadata.ratios[idx,]
  metadata.ratios_fn <- paste0(out_dir, ds, '.', ct, '.metadata.ratios.rds')
  saveRDS(metadata.ratios, metadata.ratios_fn)
  pbmc@meta.data <- cbind(pbmc@meta.data, metadata.ratios[,c('ratio.up_down','up_score')]) 
  var_stat.vec <- c('up_score', 'ratio.up_down')
  
  # Check DeltaDE score
  features_var <- 'up_score'
  title_var <- expression(Delta ~ 'DE score')
  umap.up_score <- FeaturePlot(pbmc,
                               features = features_var) + scale_colour_gradient2() + 
    labs(colour = title_var) + 
    theme(plot.title=element_blank()) + xlab('UMAP 1') + ylab('UMAP 2')
  umap.up_score.fn <- paste0(out_dir,  ds, '.', ct, '.UMAP.', features_var, '.gradient.png')
  ggsave(umap.up_score.fn, umap.up_score, width = 5, height = 4.5)
  
  # UMAPs CMV_status
  group.by_var <- phe
  umap.CMV <- DimPlot(pbmc,
                      reduction = 'umap',
                      group.by = group.by_var) +
    labs(colour = paste0(group.by_var)) + 
    theme(plot.title = element_blank(),
          legend.title = element_text(hjust = 0.5, vjust = 0.5)) + xlab('UMAP 1') + ylab('UMAP 2')
  umap.CMV.fn <- paste0(out_dir, ds, '.', ct, '.UMAP.',group.by_var,'.png')
  print(paste0('Saving UMAP colored by date in: ',umap.CMV.fn))
  ggsave(umap.CMV.fn, umap.CMV, width = 5, height = 5)
  
  # FeaturePlot for a specific gene (opt$gene = B3GAT1+ cells)
  FeaturePlot.gene <- FeaturePlot_scCustom(seurat_object = pbmc, features = gene)
  FeaturePlot.gene.fn <- paste0(out_dir, ds, '.', ct, '.UMAP.', gene,'.png')
  print(paste0('Saving UMAP colored by ', gene, ' in: ', FeaturePlot.gene.fn))
  ggsave(FeaturePlot.gene.fn, FeaturePlot.gene, width = 5, height = 5)
  
  # Check gene-specific (B3GAT1) expression VS. up_score
  ####### Correlation in B3GAT1+ cells ####### 
  ## pick cells with B3GAT1 expression
  g_vec <- as.vector(pbmc[which(rownames(pbmc)==gene),]@assays$RNA@counts)
  expr_idx <- which(g_vec>0)
  g_vec_norm <- as.vector(pbmc[which(rownames(pbmc)==gene),]@assays$RNA@data)
  g_vec_norm.expr <- g_vec_norm[expr_idx]
  
  ## pick up score from cells with B3GAT1 expression
  cell_md <- pbmc@meta.data
  cell_md.expr <- cell_md[expr_idx,]
  
  ## merge cell metadata and gene expression
  cell_md.expr$gene <- g_vec_norm.expr
  
  sp <- ggplot(cell_md.expr, aes(x = gene, y = up_score)) +
    geom_point() +
    geom_smooth(method=lm) +
    stat_cor(method = "spearman") +
    xlab(paste0(gene)) + 
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          strip.text = element_text(hjust = 0.5, face="bold", size = 10),
          panel.spacing = unit(1, "lines"),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size = 10),
          axis.text.y = element_text(size = 10))
  sp_fn <- paste0(out_dir, ds, '.', ct, '.', gene, '.scatter_plot.png')
  ggsave(sp_fn, sp)
  
  ####### Wilcoxon in B3GAT1+/- cells ####### 
  cell_md$gene <- g_vec_norm
  gene_tag <- c(paste0(ct,' ', gene, '+'), paste0(ct,' ', gene, '-'))
  cell_md$gene_discrete <- ifelse(cell_md$gene>0, gene_tag[1], gene_tag[2]) 
  y_lab <- expression(Delta ~ 'DE score')
  cols_vec <- c('#8D90E2', '#FDE787')
  names(cols_vec) <- gene_tag
  bp <- ggboxplot(cell_md, x = "gene_discrete", y = "up_score",
                  fill = "gene_discrete")
  bp <- bp + stat_compare_means() + 
    ylab(y_lab) + xlab(NULL) + 
    scale_fill_manual(values = cols_vec, guide = "none")
  bp_fn <- paste0(out_dir, ds, '.', ct, '.', gene, '.box_plot.png')
  ggsave(bp_fn, bp, width = 4.5, height = 4.5)
}

###################### Set Variables #####################
# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Cell types
celltypes <- unlist(str_split(opt$cell_types, ','))

# Assay
# opt$assay <- 'RNA'

# logFC threshold
# opt$logFC_th <- 0
logFC_th.tag <- paste0('logFC_', as.character(opt$logFC_th))

# Input directories
dea.dir <- paste0(opt$dea_dir, '/', datasets_tag, '/', opt$cell_level, '/')
so.dir <- paste0(opt$so_dir, '/')
out.dir <- paste0(opt$out_dir, '/', datasets_tag, '/', 
                  opt$assay, '/', opt$cell_level, '/', opt$phenotype, '/', logFC_th.tag, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

#################### AddModuleScore + UMAPs #################### 
# Apply function
AddModuleScore_UMAP.out <- sapply(celltypes, function(i)
  sapply(datasets, function(j) AddModuleScore_UMAP.byDS(ct = i,
                                                        ds = j,
                                                        phe = opt$phenotype, 
                                                        assay_so = opt$assay, 
                                                        logFC_th = opt$logFC_th, 
                                                        gene = opt$gene, 
                                                        out_dir = out.dir), simplify = FALSE), simplify = FALSE)