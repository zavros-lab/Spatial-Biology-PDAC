library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(writexl)
################################################################################
############################## Load/Prep Data ##################################
################################################################################
### Get Working Directory
working_directory = getwd()

### Load Seurat
data = readRDS(paste0(working_directory, 'SMI-PDAC-Res_TMA_IIT-Ver2-All_V1-Merged-WithNiche-062624.rds'))


### Load cluster Annotations
annots = read.table('/PDAC_cluster_identifications.csv',
                    sep = ',',
                    header = 1)
colnames(annots) = c('seurat_clusters', 'Cell Type')
annots$seurat_clusters = as.factor(annots$seurat_clusters)

annots[4,2] = 'Immune cells1'
annots[10,2] = "Acinar/CSC"
annots[14,2] = 'CTL'
annots[19,2] = "CD8+ T cells1"
annots[22, 2] = "CD8+/CD4+"
annots[26,2] = "CD8+ T cells2"

### Join Annotations
temp = left_join(data@meta.data, annots, by = 'seurat_clusters')
rownames(temp) = rownames(data@meta.data)

### Add cell annotations to meta
data@meta.data = temp

### Gene List: NOT FOUND: GFAP, AMY1X, CD15, DES
genes = c('CD44', 'TACSTD2', 'EPCAM', 'VIM', 'ACTA2', 'ITGAM',
          'CD163', 'SOX9', 'KRT10', 'CD8A', 'CD8B', 'GZMB', 'SPP1')


### Set Colors
cols = c("CSC" = 'darkolivegreen',
         "Schwann cells" = 'hotpink2',
         "Pancreatic stellate cells (PaSC)" = 'deepskyblue2',
         "Immune cells1" = 'aquamarine4',###
         "Alpha cells/MSC" = 'darkorchid3',
         "Alpha cells/NK" = 'brown4',
         "Monocytes" = 'darkgoldenrod',
         "Delta cells" = 'mediumseagreen',
         "MDSC/Macrophages" = 'olivedrab4',
         "Acinar/CSC" = 'firebrick1', ###
         "Acinar cells" = 'yellow',
         "Fibroblasts" = 'springgreen',
         "Immune cells2" = 'navyblue',
         "CTL" = 'chartreuse2', ### 
         "iCAF" = 'cyan',
         "Beta cells" = 'slateblue3',
         "Hepatic stellate cells" = 'chocolate1',
         "Acinar cells" = 'aquamarine3',  
         "CD8+ T cells1" = 'gold2',      ###               
         "myCAF" = 'blue2',                           
         "Kupffer cells" = 'magenta2',           
         "CD8+/CD4+" = 'royalblue4',      ###              
         "Ductal cells" = 'mediumorchid',
         "CD4+ T cells" = 'tan',                    
         "Tumor cells" = 'paleturquoise4',
         "CD8+ T cells2" = 'mediumpurple1') ###

### Sample
data$tissue %>% unique()
samples = c('TMA_N_D4', 'TMA_C2', 'TMA_C3', 'TMA_C4', 'TMA_C6', "R210211","PDAC_IIT_009","PDAC_IIT_010","PDAC_IIT_011", "PDAC_IIT_015", "PDAC_IIT_004", "R230170")

### Remove spaces from tissue names
data$tissue = gsub(' ', '_', data$tissue)

### Not in
`%nin%` = Negate(`%in%`)

### UMAP as reference
DimPlot(data, reduction = 'umap', group.by = 'Cell Type', cols = cols)

################################################################################
########################### Generate FeaturePlots ##############################
################################################################################
### Parse 
for(i in 1:length(samples)){
### Subset seurat
temp = subset(data, subset = tissue == paste0(samples[i]))

### Add coordinates as reduction
coords = as.matrix(temp@meta.data[c('CenterX_global_px', 'CenterY_global_px')])
colnames(coords) = c('InSitu_1', 'InSitu_2')
InSitu <- CreateDimReducObject(embeddings = coords,
                               assay = 'Nanostring', 
                               key = 'InSitu_')
temp[["InSitu"]] <- InSitu

### Set pt.size
if(ncol(temp)>=15000){
  size = 0.75
}else if(ncol(temp)>=10000 & ncol(temp) < 15000){
  size = 1
}else if(ncol(temp) < 10000){
  size = 2
}

### Plot Cell types
pdf(file = paste0(working_directory, '/gene_expression_analysis/',samples[i],'/', samples[i], '_celltype.pdf'), width = 12, height = 12)
print(DimPlot(temp, group.by = c('Cell Type'), reduction = 'InSitu', pt.size =size, cols = cols) + coord_equal())
dev.off()

### Plot Niches
pdf(file = paste0(working_directory, '/gene_expression_analysis/',samples[i],'/', samples[i], '_niche.pdf'), width = 12, height = 12)
print(DimPlot(temp, group.by = c('nclust'), reduction = 'InSitu', pt.size =size)+coord_equal())
dev.off()

### cycle through genes
for(j in 1:length(genes)){
  pdf(file = paste0(working_directory,'/gene_expression_analysis/',samples[i],'/', samples[i], '_', genes[j],'.pdf'), width = 12, height = 12)
  print(FeaturePlot(temp, reduction = 'InSitu', order = T, features = genes[j], pt.size = size, ncol = 1)+ coord_equal())
  dev.off()
  }
}


################################################################################
##################### GLOBAL CD44 & TROP2 Quantification #######################
################################################################################
### Tumor Site Groupings
data$tumor_site = data$tissue
data$tumor_site = ifelse(
  data$tumor_site == 'R210015',
  'Pancreas Intestinal-type',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == 'R210211' | data$tumor_site == 'R230170' | data$tumor_site == 'TMA_A2' | data$tumor_site == 'TMA_A3' | data$tumor_site == 'TMA_A4' | data$tumor_site == 'TMA_B5' | data$tumor_site == 'TMA_B6' | data$tumor_site == 'TMA_C2' | data$tumor_site == 'TMA_C3' | data$tumor_site == 'TMA_C4' | data$tumor_site == 'TMA_C5' | data$tumor_site == 'TMA_C6',
  'Pancreas',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == 'R220019',
  'PNet',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == 'TMA_B3' | data$tumor_site == 'TMA_B4',
  'Mucinous Adenocarcinoma',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == 'TMA_N_D3' | data$tumor_site == 'TMA_N_D4' | data$tumor_site == 'TMA_N_D5' | data$tumor_site == 'TMA_N_D6',
  'Cancer Adjacent Pancreas',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == "PDAC_IIT_004" | data$tumor_site == "PDAC_IIT_009" | data$tumor_site == "PDAC_IIT_010" | data$tumor_site == "PDAC_IIT_015",
  'Liver Metastatic',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == 'PDAC_IIT_006' | data$tumor_site == 'PDAC_IIT_011',
  'Lung Metastatic',
  data$tumor_site
)

data$tumor_site = ifelse(
  data$tumor_site == 'TMA_A5' | data$tumor_site == 'TMA_A6' | data$tumor_site == 'TMA_B2',
  'Invasive Duodenum',
  data$tumor_site
)



### Data Preprocessing
DefaultAssay(data) = 'Nanostring'
data = NormalizeData(data, normalization.method = 'LogNormalize')
data = ScaleData(data)

### Violin Plot of CD44 Expression per Cell Type
VlnPlot(data, features = 'CD44', group.by = 'Cell Type') + NoLegend()

### Violin Plot of TACSTD2 Expression per Cell Type
VlnPlot(data, features = 'TACSTD2', group.by = 'Cell Type') + NoLegend()


### Violin Plot of CD44 Expression per tumor site
VlnPlot(data, features = 'CD44', group.by = 'tumor_site') + NoLegend()

### Violin Plot of TACSTD2 Expression per tumor site
VlnPlot(data, features = 'TACSTD2', group.by = 'tumor_site') + NoLegend()


### Violin Plot of CD44 Expression per niche
VlnPlot(data, features = 'CD44', group.by = 'nclust') + NoLegend()

### Violin Plot of TACSTD2 Expression per niche
VlnPlot(data, features = 'TACSTD2', group.by = 'nclust') + NoLegend()

################################################################################
######################### SUBSET CD44 Quantification ###########################
################################################################################
### Exclude Samples
data_subset = subset(data, tissue %nin% c('R220019', 'R210015', 'TMA_B3', 'TMA_B4'))

### Violin Plot of CD44 Expression per Cell Type
VlnPlot(data_subset, features = 'CD44', group.by = 'Cell Type') + NoLegend()

### Violin Plot of TACSTD2 Expression per Cell Type
VlnPlot(data_subset, features = 'TACSTD2', group.by = 'Cell Type') + NoLegend()


### Violin Plot of CD44 Expression per tumor site
VlnPlot(data_subset, features = 'CD44', group.by = 'tumor_site') + NoLegend()

### Violin Plot of TACSTD2 Expression per tumor site
VlnPlot(data_subset, features = 'TACSTD2', group.by = 'tumor_site') + NoLegend()


### Violin Plot of CD44 Expression per niche
VlnPlot(data_subset, features = 'CD44', group.by = 'nclust') + NoLegend()

### Violin Plot of TACSTD2 Expression per niche
VlnPlot(data_subset, features = 'TACSTD2', group.by = 'nclust') + NoLegend()



################################################################################
#################### SUBSET Epithelial CD44 Quantification #####################
################################################################################
### Subset
data_subset$celltype = data_subset$`Cell Type`
data_epi = subset(data_subset, subset = celltype %in% c('Acinar cells', 'Acinar/CSC', 'CSC', 'Tumor cells'))

### Violin Plot of CD44 Expression per Cell Type
pdf(file = paste0(working_directory, '/Volumes/LABTAKEHOME/Pancreas_cosmx/gene_expression_analysis/Quant/CD44_epi_cells.pdf'), width = 5, height = 5)
VlnPlot(data_epi, features = 'CD44', group.by = 'Cell Type') + NoLegend()
dev.off()

### Violin Plot of TACSTD2 Expression per Cell Type
pdf(file = paste0(working_directory, '/gene_expression_analysis/Quant/TACSTD2_epi_cells.pdf'), width = 5, height = 5)
VlnPlot(data_epi, features = 'TACSTD2', group.by = 'Cell Type') + NoLegend()
dev.off()

### Violin Plot of CD44 Expression per tumor site
pdf(file = paste0(working_directory, '/gene_expression_analysis/Quant/CD44_tissue_site.pdf'), width = 5, height = 5)
VlnPlot(data_epi, features = 'CD44', group.by = 'tumor_site') + NoLegend()
dev.off()

### Violin Plot of TACSTD2 Expression per tumor site
pdf(file = paste0(working_directory, '/gene_expression_analysis/Quant/TACSTD2_tissue_site.pdf'), width = 5, height = 5)
VlnPlot(data_epi, features = 'TACSTD2', group.by = 'tumor_site') + NoLegend()
dev.off()

### Violin Plot of CD44 Expression per niche
pdf(file = paste0(working_directory, '/gene_expression_analysis/Quant/CD44_niche.pdf'), width = 5, height = 5)
VlnPlot(data_epi, features = 'CD44', group.by = 'nclust') + NoLegend()
dev.off()

### Violin Plot of TACSTD2 Expression per niche
pdf(file = paste0(working_directory, '/gene_expression_analysis/Quant/TACSTD2_niche.pdf'), width = 5, height = 5)
VlnPlot(data_epi, features = 'TACSTD2', group.by = 'nclust') + NoLegend()
dev.off()




################################################################################
########### Differentially Expressed Genes Across All Cell Types ###############
################################################################################
### Set Idents 
Idents(data) = 'Cell Type'

### Find Markers
markers = FindAllMarkers(data, only.pos = T)

### Filter out control probes
markers <- markers[!grepl("Neg|System", markers$gene), ]

### Get top 10
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

png(file = paste0(working_directory, '/gene_expression_analysis/Quant/hm_celltype.png'), width = 1500, height = 1200)
DoHeatmap(data, features = top5$gene) + NoLegend()
dev.off()

### Write markers to xlsx
write_xlsx(markers, path = paste0(working_directory, '/gene_expression_analysis/Quant/final/DEGs.xlsx'))





################################################################################
######### Visualize Gene Expression Differences between Cell types #############
################################################################################
### Subset
data_subset$celltype = data_subset$`Cell Type`
data_epi = subset(data_subset, subset = celltype %in% c('Acinar cells', 'Acinar/CSC', 'CSC', 'Ductal cells'))

### Add cell number 
data_epi$celltypes = data_epi$`Cell Type`
data_epi$celltypes = ifelse(
  data_epi$celltypes == 'Acinar cells',
  '17: Acinar cells',
  data_epi$celltypes
)

data_epi$celltypes = ifelse(
  data_epi$celltypes == 'Acinar/CSC',
  '9/10: Acinar/CSC',
  data_epi$celltypes
)

data_epi$celltypes = ifelse(
  data_epi$celltypes == 'CSC',
  '20: Ductal CSC',
  data_epi$celltypes
)

data_epi$celltypes = ifelse(
  data_epi$celltypes == 'Ductal cells',
  '22: Cancer cells',
  data_epi$celltypes
)


cols = c("20: Ductal CSC" = 'darkolivegreen',
         "9/10: Acinar/CSC" = 'firebrick1', 
         "17: Acinar cells" = 'aquamarine3',
         "22: Cancer cells" = 'mediumorchid')
VlnPlot(data_epi, features = c('CD44','SOX2', 'KRT19'), group.by = 'celltypes',cols = cols) + NoLegend()


##### CD44 #####
# Ensure grouping is a factor
data_epi$celltypes <- factor(data_epi$celltypes)

# Extract expression + group data
cd44_data <- FetchData(data_epi, vars = c("CD44", "celltypes"))

# Generate all pairwise combinations
group_levels <- levels(data_epi$celltypes)
pairwise <- combn(group_levels, 2, simplify = FALSE)

# Run Wilcoxon tests for each pair
results <- lapply(pairwise, function(pair) {
  group1 <- cd44_data %>% filter(celltypes == pair[1]) %>% pull(CD44)
  group2 <- cd44_data %>% filter(celltypes == pair[2]) %>% pull(CD44)
  p <- wilcox.test(group1, group2)$p.value
  list(group1 = pair[1], group2 = pair[2], p = p)
})

# Filter for significant comparisons
sig_results <- Filter(function(x) x$p < 0.05, results)

# Create comparison list for plotting
sig_comparisons <- lapply(sig_results, function(x) c(x$group1, x$group2))

# get pval
print(data.frame(
  comparison = sapply(sig_results, function(x) paste(x$group1, "vs", x$group2)),
  p_value = sapply(sig_results, function(x) signif(x$p, 3))
))

# Define label.y values 
label_y_positions <- seq(7, 9, length.out = length(sig_comparisons))

#Plot 
pdf(file = paste0(working_directory, '/gene_expression_analysis/Quant/final/cd44_individual.pdf'), width = 5, height = 6)
cd44 = VlnPlot(data_epi, features = "CD44", group.by = "celltypes", cols = cols) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", 
                     comparisons = sig_comparisons, 
                     label.y = label_y_positions) +
  NoLegend() +
  theme(axis.title.x = element_blank()) +
  coord_cartesian(ylim = c(0, 9.5))+
  ylim(-4,30)
cd44
dev.off()




##### SOX2 #####
# Ensure grouping is a factor
data_epi$celltypes <- factor(data_epi$celltypes)

# Extract expression + group data
SOX2_data <- FetchData(data_epi, vars = c("SOX2", "celltypes"))

# Generate all pairwise combinations
group_levels <- levels(data_epi$celltypes)
pairwise <- combn(group_levels, 2, simplify = FALSE)

# Run Wilcoxon tests for each pair
results <- lapply(pairwise, function(pair) {
  group1 <- SOX2_data %>% filter(celltypes == pair[1]) %>% pull(SOX2)
  group2 <- SOX2_data %>% filter(celltypes == pair[2]) %>% pull(SOX2)
  p <- wilcox.test(group1, group2)$p.value
  list(group1 = pair[1], group2 = pair[2], p = p)
})

# Filter for significant comparisons
sig_results <- Filter(function(x) x$p < 0.05, results)

# Create comparison list for plotting
sig_comparisons <- lapply(sig_results, function(x) c(x$group1, x$group2))

# Print p-values
print(data.frame(
  comparison = sapply(sig_results, function(x) paste(x$group1, "vs", x$group2)),
  p_value = sapply(sig_results, function(x) signif(x$p, 3))
))

# Define label.y values 
label_y_positions <- seq(7, 9, length.out = length(sig_comparisons))

# plot with manual y-label positions
pdf(file = paste0('/gene_expression_analysis/Quant/final/SOX2_individual.pdf'), width = 5, height = 6)
SOX2 = VlnPlot(data_epi, features = "SOX2", group.by = "celltypes", cols = cols) +
  NoLegend() +
  theme(axis.title.y = element_blank()) +
  coord_cartesian(ylim = c(0, 9.5))+
  ylim(-4,30)
SOX2
dev.off()



##### KRT19 #####
# Ensure grouping is a factor
data_epi$celltypes <- factor(data_epi$celltypes)

# Extract expression + group data
KRT19_data <- FetchData(data_epi, vars = c("KRT19", "celltypes"))

# Generate all pairwise combinations
group_levels <- levels(data_epi$celltypes)
pairwise <- combn(group_levels, 2, simplify = FALSE)

# Run Wilcoxon tests for each pair
results <- lapply(pairwise, function(pair) {
  group1 <- KRT19_data %>% filter(celltypes == pair[1]) %>% pull(KRT19)
  group2 <- KRT19_data %>% filter(celltypes == pair[2]) %>% pull(KRT19)
  p <- wilcox.test(group1, group2)$p.value
  list(group1 = pair[1], group2 = pair[2], p = p)
})

# Filter for significant comparisons
sig_results <- Filter(function(x) x$p < 0.05, results)

# Create comparison list for plotting
sig_comparisons <- lapply(sig_results, function(x) c(x$group1, x$group2))

# Print p-values
print(data.frame(
  comparison = sapply(sig_results, function(x) paste(x$group1, "vs", x$group2)),
  p_value = sapply(sig_results, function(x) signif(x$p, 3))
))

# Define label.y values 
label_y_positions <- seq(7, 9, length.out = length(sig_comparisons))

#  plot with manual y-label positions
pdf(file = paste0('/gene_expression_analysis/Quant/final/KRT19_individual.pdf'), width = 5, height = 6)
KRT19 = VlnPlot(data_epi, features = "KRT19", group.by = "celltypes", cols = cols) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", 
                     comparisons = sig_comparisons, 
                     label.y = label_y_positions) +
  NoLegend() +
  theme(axis.title = element_blank(),
        title = element_blank()) +
  coord_cartesian(ylim = c(0, 9.5))+
  ylim(-4,30)
KRT19
dev.off()



##### Composite #####
pdf(file = paste0('/gene_expression_analysis/Quant/final/Composite.pdf'), width = 12, height = 6)
(cd44|SOX2|KRT19)
dev.off()



