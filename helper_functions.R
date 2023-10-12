library(ggforce)
library(ggplot2)
filter_cells <- function(seurat, count_cutoff, n_genes_cutoff, mitochondria_cutoff){
  qc_df <- data.frame(nCount=seurat$nCount_RNA, nFeature=seurat$nFeature_RNA,
                      perc_mito=seurat$percent_mt,
                      perc_hemo=seurat$percent_hb,
                      perc_ribo=seurat$percent_rb,
                      orig.ident=seurat$orig.ident,
                      barcode=names(seurat$orig.ident))
  selected <- WhichCells(seurat, expression = percent_mt < mitochondria_cutoff)
  to_filter <- qc_df[log10(qc_df$nCount) < count_cutoff | qc_df$nFeature < n_genes_cutoff | 
                       !qc_df$barcode%in%selected,]
  seurat_filt <- seurat[,!colnames(seurat) %in% to_filter$barcode ]
  return(seurat_filt)
}

filter_genes <- function(seurat){
  # Filter Mitocondrial genes
  seurat_noMT <- seurat[!grepl("^MT-", rownames(seurat)), ]
  seurat_noMT_noRibosome <- seurat_noMT[!grepl("^RP[SL]", rownames(seurat_noMT)), ]
  
  if('ADT' %in% names(seurat@assays)){
    seurat_noMT_noRibosome[['ADT']] <- seurat[['ADT']]
  }
  return(seurat_noMT_noRibosome)
}



calculate_bcmvn <- function(seurat){
  cat('Run PCA\n')
  seurat <- RunPCA(seurat, verbose = F, npcs = 20)
  cat('Run UMAP\n')
  suppressWarnings(seurat <- RunUMAP(seurat, dims = 1:10, verbose = F))
  
  # optimize parameters
  cat('Run paramsweep')
  suppressWarnings(sweep.res <- paramSweep_v3(seurat, sct = T))
  sweep.stats <- summarizeSweep(sweep.res,GT = FALSE) 
  
  bcmvn <- find.pK(sweep.stats)
  return(bcmvn)
  
}

add_vdj <- function(vdj_path, seurat_obj){
  vdj <- data.frame()
  tryCatch(vdj <- read.table(vdj_path, header=T,sep=','), 
           error=function(e) vdj <- data.frame())
  if(nrow(vdj)==0){
    return(seurat_obj)
  }
  
  if(length(vdj[vdj$d_gene=='',]$d_gene) > 0){
    vdj[vdj$d_gene=='',]$d_gene <- NA
  }
  if(length(vdj[vdj$c_gene=='',]$c_gene) > 0){
    vdj[vdj$c_gene=='',]$c_gene <- NA
  }
  
  # add each of the vdj columns to the seurat object if the barcodes matchs
  for(col in colnames(vdj)[!grepl('barcode',colnames(vdj))]){
    seurat_obj@meta.data[[col]] <- NA
    seurat_obj@meta.data[[col]] <- vdj[match(rownames(seurat_obj@meta.data),vdj$barcode),][[col]]
  }
  
  return(seurat_obj)
}


filter_doublets <- function(seurat, seurat_doublets){
  # remove doublets based on DoubletFinder
  seurat <- seurat[,colnames(seurat) %in% colnames(seurat_doublets)]
  df_column <- colnames(seurat_doublets@meta.data)[grepl('DF.classifications', colnames(seurat_doublets@meta.data))]
  seurat@meta.data$DF.classifications <- seurat_doublets@meta.data[rownames(seurat@meta.data),][[df_column]]
  seurat <- seurat[, colnames(seurat) %in% rownames(seurat@meta.data[seurat@meta.data$DF.classifications=='Singlet',])]
  
  # remove doublets based on Vireo
  if('doublet' %in% seurat@meta.data$donor){
    seurat <- seurat[,!colnames(seurat) %in% rownames(seurat@meta.data[seurat@meta.data$donor=='doublet',])]
  }
  return(seurat)
}

regress_celltypes <- function(seurat){
  # perform cell cycle analysis (make sure to specify the "assay" parameter
  seurat <- CellCycleScoring(
    seurat,
    s.features = s.genes,
    g2m.features = g2m.genes,
    assay = 'SCT',
    set.ident = TRUE)
  
  # normalise again but this time including also the cell cycle scores
  suppressWarnings(seurat <- SCTransform(
    seurat,
    assay = 'RNA',
    new.assay.name = 'SCT',
    vars.to.regress = c('percent_mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'),
    return.only.var.genes = FALSE
  ))
  return(seurat)
}



calculate_mad <- function(seurat_list, features=c("percent_rb","percent_mt","nCount_RNA","nFeature_RNA")){
  MAD_df_list <- lapply(seurat_list, function(x){
    as.data.frame(matrix(ncol = 2, nrow = length((which(features %in% colnames(x@meta.data))))))
  })
  
  MAD_df_list <- lapply(names(MAD_df_list), function(x){
    rownames(MAD_df_list[[x]]) <- colnames(seurat_list[[x]]@meta.data)[which(colnames(seurat_list[[x]]@meta.data) %in% features)]
    colnames(MAD_df_list[[x]]) <- c("Median","MAD")
    for (QC in rownames(MAD_df_list[[x]])){
      MAD_df_list[[x]][QC, "Median"] <- median(seurat_list[[x]]@meta.data[,QC], na.rm = TRUE)
      MAD_df_list[[x]][QC, "MAD"] <- mad(seurat_list[[x]]@meta.data[,QC], center = MAD_df_list[[x]][QC, "Median"],  
                                         constant = 1.4826, na.rm = TRUE,low = FALSE, high = FALSE)
    }
    MAD_df_list[[x]]$Pool <- x
    MAD_df_list[[x]]$QC_Metric <- rownames(MAD_df_list[[x]])
    return(MAD_df_list[[x]])
  })
  MAD_df <- do.call(rbind,MAD_df_list)
  MAD_df$Pool <- as.factor(MAD_df$Pool)
  return(MAD_df)
}




mad_figures <- function(seurat_list, MAD_df, outdir, mito_ribo_features=c('percent_mt','percent_rb')){
  
  violins_MADperPOOL <- list()
  violins <- list()
  
  merged_qc <- data.frame()
  for(pool in seurat_list){
    merged_qc <- rbind(merged_qc, pool@meta.data)
  }
  
  for (QC in unique(MAD_df$QC_Metric)){
    #for(QC in c('nCount_RNA')){
    if (QC %in% mito_ribo_features){
      violins_MADperPOOL[[QC]] <- ggplot(merged_qc, aes_string( x = 'orig.ident', 
                                                                y = QC)) +
        geom_violin() +
        geom_sina(size = 1, alpha = 0.6) +
        theme_classic() +
        ylim(0, min(max(merged_qc[,QC]),100)) +
        labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+1.5*MAD, yend=Median+1.5*MAD, col = "1.5 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
        
        scale_color_manual("MAD", values=c("Median"="grey","1 MAD"="blue3", "1.5 MAD"="darkviolet", "2 MAD" = "firebrick1", "3 MAD" = "darkorange1", "4 MAD" = "gold1")) +
        theme(text = element_text(size=14),
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
      
      
    } else {
      violins_MADperPOOL[[QC]] <- ggplot(merged_qc, aes_string( x = 'orig.ident', y = QC)) +
        geom_violin() +
        geom_sina(size = 1, alpha = 0.6) +
        theme_classic() +
        ylim(0, NA) +
        labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-1.5*MAD, yend=Median-1.5*MAD, col = "1.5 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
        
        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
        scale_color_manual("MAD", values=c("Median"="grey","1 MAD"="blue3","1.5 MAD"="darkviolet", "2 MAD" = "firebrick1", "3 MAD" = "darkorange1", "4 MAD" = "gold1")) +
        theme(text = element_text(size=14),
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
    }
    violins[[QC]] <- ggplot(merged_qc, aes( x = orig.ident, y = merged_qc[,QC])) +
      geom_violin() +
      geom_sina(size = 1, alpha = 0.6) +
      theme_classic() +
      ylim(0, NA) +
      labs(title = paste0(QC, " Violin Plot"), y = QC, x = "orig.ident") +
      theme(text = element_text(size=14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  for (QC in names(violins_MADperPOOL)){
    ggsave(violins_MADperPOOL[[QC]], filename = paste0(outdir, "/", QC, "_violin_MADper_Pool.png"), width = 29.7, height = 21 ,units = c("cm"))
    ggsave(violins[[QC]], filename = paste0(outdir, "/", QC, "_violin_noMADlines.png"), width = 29.7, height = 21 ,units = c("cm"))
  }
  return(violins_MADperPOOL)
}


### Define a function to identify cells that are outliers based on certain MAD from the median
## This will return a seurat object that contains new columns that indicate whether each cell is an outlier for that QC metric
mad_function <- function(seurat, column, number_mad){
  mad <- mad(seurat@meta.data[,column])
  low <- median(seurat@meta.data[,column]) - number_mad*mad
  print("The lower bound is:")
  print(low)
  seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low),"NotOutlier", "Outlier")
  return(seurat)
}
mad_function_mt <- function(seurat, column, number_mad){
  mad <- mad(seurat@meta.data[,column])
  high <- median(seurat@meta.data[,column]) + number_mad*mad
  print("The upper bound is:")
  print(high)
  seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
  return(seurat)
}


Merge_Seurat_List <- function(
    list_seurat,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "SeuratProject"
) {
  # Check list_seurat is list
  if (!inherits(x = list_seurat, what = "list")) {
    cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
  }
  
  # Check list_seurat is only composed of Seurat objects
  for (i in 1:length(x = list_seurat)) {
    if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
      cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
    }
  }
  
  merged_object <- reduce(list_seurat, function(x, y) {
    merge(x = x, y = y, add.cell.ids = add.cell.ids, merge.data = merge.data, project = project)
  })
}


