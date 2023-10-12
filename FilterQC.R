#########################################################
# Quality Filtering scRNA-seq data
# Created by Tarran Rupal and Niek de Klein
# Adapted by Thais Crippa
# Last data modification: 08.Oct.2023
#########################################################

### Install Packages
packages <- c('Seurat', "stringr", "ggplot2" ,"ggpubr", "plyr", "tidyverse", "grid", "gridExtra", "ggExtra", "patchwork", "RColorBrewer", "ggbee", "warm", "ggsci", "viridis")

if(sum(as.numeric(!packages %in% installed.packages())) != 0){
  instalador <- packages[!packages %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(packages, require, character = T) 
} else {
  sapply(packages, require, character = T) 
}


library('Seurat')
library(stringr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggExtra)
library(patchwork)
library(RColorBrewer)
library(ggbeeswarm)
library(ggsci)
library(viridis)

# Plot Parameters
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

Theme <-  theme_pubr() + theme(
  axis.ticks=element_blank(),
  panel.grid.major = element_line(),
  panel.border = element_blank(),
  axis.title = element_text(size =12),
  plot.title = element_text(size = 18, hjust = .5),
  legend.text = element_text(size = 14))

### Set environment - Change for your own directory
basedir <- '/Users/tc24/Desktop/scRNAseq'
setwd('/Users/tc24/Desktop/scRNAseq')
source('/Users/tc24/Desktop/scRNAseq/helper_functions.R')

#Creating subdirectories
dir.create(outdir,showWarnings = F)
figdir <- paste0(outdir,'/figures/')
dir.create(outdir,showWarnings = F)
figdir <- paste0(outdir,'/figures/mad/')

seurat_list <- readRDS('Jaguar_unclearData.RDS')

#########################################################
# 1s STEP - FILTERING ON COUNTS GENES AND MITO READS
#########################################################

mad_filt <- calculate_mad(seurat_list)

### Inspect elements
violins_filt <- mad_figures(seurat_list, mad_filt, paste0('/figures/mad/'))

# Filtering Manually by MAD (Median Allele Deviation)
### Filter is done by inspecting the lines of MAD and, by pool, 
### delimiting the MAD (1,1.5, 2, 3, 4) that is most close to the bottom 
### of the round part of the violin plot

###########FILTERING ON NCOUNT (reads counts) ###########
violins_filt$nCount_RNA

seurat_list$BCELL_MXP1 <- mad_function(seurat = seurat_list$BCELL_MXP1,
                                       column = "nCount_RNA", number_mad =1.5)
seurat_list$BCELL_MXP2 <- mad_function(seurat = seurat_list$BCELL_MXP2,
                                       column = "nCount_RNA", number_mad = 1.5)
seurat_list$BCELL_MXP3 <- mad_function(seurat = seurat_list$BCELL_MXP3,
                                       column = "nCount_RNA", number_mad = 2)
seurat_list$PBMC_MXP1 <- mad_function(seurat = seurat_list$PBMC_MXP1,
                                      column = "nCount_RNA", number_mad = 2)
seurat_list$PBMC_MXP2 <- mad_function(seurat = seurat_list$PBMC_MXP2,
                                      column = "nCount_RNA", number_mad = 2)
seurat_list$PBMC_MXP3 <- mad_function(seurat = seurat_list$PBMC_MXP3,
                                      column = "nCount_RNA", number_mad = 2)

###########FILTERING ON NFEATURE (RNA) ###########
violins_filt$nFeature_RNA

seurat_list$BCELL_MXP1 <- mad_function(seurat = seurat_list$BCELL_MXP1,
                                       column = "nFeature_RNA", number_mad =1.5)
seurat_list$BCELL_MXP2 <- mad_function(seurat = seurat_list$BCELL_MXP2,
                                       column = "nFeature_RNA", number_mad = 1.5)
seurat_list$BCELL_MXP3 <- mad_function(seurat = seurat_list$BCELL_MXP3,
                                       column = "nFeature_RNA", number_mad = 2)
seurat_list$PBMC_MXP1 <- mad_function(seurat = seurat_list$PBMC_MXP1,
                                      column = "nFeature_RNA", number_mad = 2)
seurat_list$PBMC_MXP2 <- mad_function(seurat = seurat_list$PBMC_MXP2,
                                      column = "nFeature_RNA", number_mad = 2.5)
seurat_list$PBMC_MXP3 <- mad_function(seurat = seurat_list$PBMC_MXP3,
                                      column = "nFeature_RNA", number_mad = 2.5)


###########FILTERING ON PERCENT MITOCHONDRIAL ###########
# We are going to remove only greatest numbers of mit perc.
violins_filt$percent_mt

seurat_list$BCELL_MXP1 <- mad_function_mt(seurat = seurat_list$BCELL_MXP1,
                                          column = "percent_mt", number_mad = 3)
seurat_list$PBMC_MXP2<- mad_function_mt(seurat = seurat_list$PBMC_MXP2,
                                        column = "percent_mt", number_mad = 3)
seurat_list$BCELL_MXP3 <- mad_function_mt(seurat = seurat_list$BCELL_MXP3,
                                          column = "percent_mt", number_mad = 3)
seurat_list$PBMC_MXP1<- mad_function_mt(seurat = seurat_list$PBMC_MXP1,
                                        column = "percent_mt", number_mad = 3)
seurat_list$BCELL_MXP2 <- mad_function_mt(seurat = seurat_list$BCELL_MXP2,
                                          column = "percent_mt", number_mad = 3)
seurat_list$PBMC_MXP3 <- mad_function_mt(seurat = seurat_list$PBMC_MXP3,
                                         column = "percent_mt", number_mad = 3)

seurat_list_filt <- seurat_list 
for(pool in names(seurat_list_filt)){
  seurat_list_filt[[pool]] <- subset(seurat_list_filt[[pool]],
                                     subset = percent_mt_mad == "NotOutlier" &
                                       nCount_RNA_mad == "NotOutlier" &
                                       nFeature_RNA_mad == "NotOutlier")
}

#########################################################
# 2nd STEP - Checking plots
#########################################################

# RNA plot
## No filter
ggplot(seurat_list$PBMC_MXP2@meta.data, aes(sample, nCount_RNA)) + 
  geom_quasirandom(aes(group = sample, colour = nCount_RNA_mad)) + Theme
# Filtered
ggplot(seurat_list_filt$PBMC_MXP2@meta.data, aes(sample, nCount_RNA)) + 
  geom_quasirandom(aes(group = sample, colour = nCount_RNA_mad)) + Theme

# MIT plot
## No filter
ggplot(seurat_list$PBMC_MXP2@meta.data, aes(sample, percent_mt)) + 
  geom_quasirandom(aes(group = sample, colour = percent_mt_mad)) + ylim(0,25) + Theme
# Filtered
ggplot(seurat_list_filt$PBMC_MXP2@meta.data, aes(sample, percent_mt)) + 
  geom_quasirandom(aes(group = sample, colour = percent_mt_mad)) + ylim(0,25) + Theme

#########################################################
# SAVE OBJ
#########################################################

#save filtered seurat list
saveRDS(seurat_list_filt,'jaguar_qcFiltered.RDS')
#save annotated filtered seurat list
saveRDS(seurat_list,'jaguar_qcFiltered_annotations.RDS')


#########################################################
# 3th STEP - QC FIGURES AFTER FILTERING
#########################################################

# merge seurat data
merged_metadata_filt <- data.frame()
for(sobj in names(seurat_list_filt)){
  merged_metadata_filt <- rbind.fill(merged_metadata_filt,seurat_list_filt[[sobj]]@meta.data)
  merged_metadata_filt <- as_tibble(merged_metadata_filt)
}

merged_metadata_annotate <- data.frame()
for(sobj in names(seurat_list)){
  merged_metadata_annotate <- rbind.fill(merged_metadata_annotate,seurat_list[[sobj]]@meta.data)
  merged_metadata_annotate <- as_tibble(merged_metadata_annotate)
}

#### som distribution plots after filtering ####
merged_metadata_filt$sample <- factor(merged_metadata_filt$sample,
                                      levels=c('PBMC_MXP1','PBMC_MXP2','PBMC_MXP3',
                                               'BCELL_MXP1','BCELL_MXP2','BCELL_MXP3')) 

p1 <- ggplot(merged_metadata_filt, aes(log10(nCount_RNA),fill=donor))+
  geom_histogram( bins=100)+
  Theme+ theme(legend.position = 'none')+
  scale_fill_manual(values = mycolors) +
  facet_wrap(~sample, scale='free_y',ncol=1)
p2 <- ggplot(merged_metadata_filt, aes(nFeature_RNA,fill=donor))+
  geom_histogram(bins=100)+
  Theme + theme(legend.position = 'none')+
  facet_wrap(~sample, scale='free_y',ncol=1)+
  scale_fill_manual(values = mycolors) 
p1+p2
ggsave('figures/mad/histo_nFeature-nCount_perDonor_allFilt.pdf', width=10, height=6)

ggplot(merged_metadata_filt, aes(sample, percent_mt, fill = sample)) + 
  geom_violin(alpha = 0.4) + geom_boxplot(width = 0.4) + Theme + 
  scale_fill_manual(values = mycolors) +ylim(0,10)
ggsave('figures/mad/violin_mito_perSample_allFilt.pdf', width=10, height=6)

ggplot(merged_metadata_filt, aes(nCount_RNA, nFeature_RNA))+
  geom_point(data=merged_metadata_filt, 
             aes(colour=percent_mt))+ Theme+
  facet_wrap(~sample, scale='free')
ggsave('figures/mad/mitoHighlight_nFeature-nCount_perSample_allFilt.pdf', width=10, height=6)

ggplot(merged_metadata_filt, aes(internalID_matched, nFeature_RNA)) + geom_boxplot(aes(fill = internalID_matched), width=0.5, outlier.size = 0.2) +
  scale_fill_viridis(discrete = TRUE) + xlab("") + theme(legend.position = "none")
ggsave('figures/mad/nFeatureperSample_allFilt.pdf', width=10, height=6)



