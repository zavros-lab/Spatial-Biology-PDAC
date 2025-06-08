# Required Libraries
library(Seurat) # Seurat v5
library(ggplot2)


# Data Input
data.combined <- readRDS("SMI-PDAC-Res_TMA_IIT-Ver2-All_V1-Merged-WithNiche-062624.rds")
data.combined$custom.ident <- data.combined@active.ident

# Subsetting
# All
data <- data.combined
# R210015
data <- subset(data.combined, subset = tissue == "R210015")
# R210211
data <- subset(data.combined, subset = tissue == "R210211")
# R230170
data <- subset(data.combined, subset = tissue == "R230170")
# R220019
data <- subset(data.combined, subset = tissue == "R220019")
# IIT 004, 006, 009, 010, 011, 015
data <- subset(data.combined, subset = tissue == "PDAC IIT 006")
# TMA - PDAC
data <- subset(data.combined, subset = tissue == "TMA_C3")
# TMA - Normal
data <- subset(data.combined, subset = tissue == "TMA_N_D4")
# All TMA
data.combined$is_tma <- startsWith(data.combined$tissue, "TMA")
data <- subset(data.combined, subset = is_tma == T)

# Colors used for plotting
ccolss <- c(
  `0` = '#556b2f', #0
  `1` = '#ed6fc8', #1
  `2` = '#45a9f8', #2
  `3` = '#53bcc2', #3
  `4` = '#bc83f8', #4
  `5` = '#741536', #5
  `6` = '#c69732', #6
  `7` = '#56ba6f', #7
  `8` = '#929c57', #8
  `9` ='#ff4500', #9
  `10` ='#ffff00', #10
  `11` ='#00ff7f', #11
  `12` ='#000e7a', #12
  `13` = '#a1fa4e', #13
  `14` = '#72fbfd', #14
  `15` = '#742ef5', #15
  `16` = '#f08533', #16
  `17` ='#a9d9d2', #17
  `18` ='#f5c747', #18
  `19` ='#0029f5', #19
  `20` = '#fb00ff', #20
  `21` = '#415c96', #21
  `22` = '#be2ffa', #22
  `23` = '#f0d2b9', #23
  `24` = '#839f97', #24
  `25` = '#b1ade0', #25
  `26` = '#abcdef', #26
  `27` = '#ababab', #27
  `28` = '#f58682', #28
  `29` = '#29677f', #29
  `30` = '#f4c6e5' #30
)

ccolss_patients <- c(
  `TMA_N_D6` = '#556b2f',
  `TMA_N_D5` = '#ed6fc8',
  `TMA_N_D4` = '#0029f5',
  `TMA_N_D3` = '#6ebac1',
  `TMA_A2` = '#bc83f8',
  `TMA_A3` = '#ecc0ec',
  `TMA_A4` = '#c69732',
  `TMA_A5` = '#71b776',
  `TMA_A6` = 'maroon',
  `TMA_B2` ='#ff4500', 
  `TMA_B3` ='#67d85a',
  `TMA_B4` ='#00ff7f',
  `TMA_B5` ='#000e7a',
  `TMA_B6` = '#a1fa4e',
  `TMA_C2` = '#72fbfd',
  `TMA_C3` = '#742ef5',
  `TMA_C4` = '#f08533',
  `TMA_C5` ='#a9d9d2', 
  `TMA_C6` ='#f4c048', 
  `R210015` ='#45a9f8', 
  `R210211` = '#fb00ff',
  `R230170` = '#415c96',
  `R220019` = '#7f7433',
  'PDAC IIT 004'='#ffff00',
  'PDAC IIT 006'='#f4c6e5',
  'PDAC IIT 009'='#701e45',
  'PDAC IIT 010'='#5a6a37',
  'PDAC IIT 011'='#e78c86',
  'PDAC IIT 015'='#bebebe'
)

#### UMAP Plots
# UMAP Plots grouped by Seurat cluster
DimPlot(data, label = T, repel = T, raster = F, cols = ccolss) # Labeled
DimPlot(data, label = T, label.box = T, repel = T, raster = F, cols = ccolss) # Labeled with boxes

# UMAP Plots grouped by patient
DimPlot(data, label = F,  raster = F, repel = T, shuffle=T, group.by = "tissue", cols = ccolss_patients)
DimPlot(data, label = T, label.box=T, raster = F, repel = T, shuffle=T, group.by = "tissue", cols = ccolss_patients)

# UMAP plots of single FOVs
DimPlot(subset(data, subset = fov == 9), split.by = "tissue")

#### Spatial UMAP Plots
# Overall
ImageDimPlot(data, fov = NULL, border.size = 0.03, border.color = "white",cols = ccolss)

# Specific FOVs
sub <- subset(data, subset = fov %in% c(16,28))
ImageDimPlot(sub, cols = ccolss, border.size = 0.03, border.color = "black")

# Specific FOV including molecules
fov_target <- 9
cells_of_interest <- data$id[(data$fov == fov_target) & (data$Run_Tissue_name == 'TMA_N_D4')]
zoom_fov <- apply(data@images$Run5713.PIT$centroids@coords[data@images$Run5713.PIT$centroids@cells %in% cells_of_interest,], 2, range)
ImageDimPlot(subset(data, subset = fov == fov_target), fov = NULL, border.size = 0.3, cols = cellcolors, border.color = "white", mols.size = 1.5) + xlim(zoom_fov[, 2]) + ylim(zoom_fov[, 1])

#### Violin Plots
# Compare patients
Idents(data) <- data$tissue
VlnPlot(data, features = c("CCND1"), slot = "scale.data", pt.size = 0) + scale_y_continuous(limits = c(0.05, 40))
Idents(data) <- data$custom.ident

# Compare clusters
Idents(data) <- data$custom.ident
VlnPlot(data, features = c("MALAT1"), slot = "scale.data", pt.size = 0)

#### Dot Plots
# Compare gene expression across a selected list of clusters
gene.list <- c("EPCAM","CD44","TACSTD2", "CD274","CD163", "CD68","ACTA2","IL6","LIFR","VIM")
DotPlot(data, features = gene.list, idents = c(1,2,7,9,12)) + coord_flip() + RotatedAxis() + theme(axis.text.x = element_text(size = 9))

# Compare FOVs within the same Seurat cluster
gene.list <- c("EPCAM", "CD44")
cluster <- 1
DotPlot(subset(data, idents = cluster), features = gene.list, group.by = "fov")

# Compare gene expression across all clusters
gene.list <- c("EPCAM","CD44","SOX2","TACSTD2")
DotPlot(data, features = gene.list) + coord_flip() + RotatedAxis() + theme(axis.text.x = element_text(size = 9))

# Compare gene expression across different patients
gene.list <- c("EPCAM", "CD44", "SOX2", "TACSTD2")
DotPlot(data, features = gene.list, group.by = "tissue") + coord_flip() + RotatedAxis() + theme(axis.text.x = element_text(size = 9))

#### Find top differentially expressed markers in a specific cluster
markers <- FindMarkers(data, min.pct = 0.0001, logfc.treshold = 0.0, ident.1 = 10)
markers <- markers[!grepl("SystemControl", rownames(markers)),]
markers <- markers[!grepl("Negative", rownames(markers)),]

#### Find top differentially expressed markers in a specific FOV
d <- data
Idents(d) <- d$fov
markers <- FindMarkers(d, min.pct = 0.0001, logfc.treshold = 0.0, ident.1 = 10)
markers <- markers[!grepl("SystemControl", rownames(markers)),]
markers <- markers[!grepl("Negative", rownames(markers)),]

#### Find top differentially expressed markers in every patient and cluster combination
data$unique.ident <- paste0(data$tissue, "-", data$custom.ident)
Idents(data) <- data$unique.ident
list <- unique(data$unique.ident)

for(i in list) {
  print(paste("Starting", i))
  
  cell_count <- length(data$id[which(data$unique.ident == i)])
  if(cell_count >= 3)
  {
    markers <- FindMarkers(data, min.pct = 0.0001, logfc.threshold = 0.0, ident.1 = i)
    markers <- markers[!grepl("SystemControl", rownames(markers)),]
    markers <- markers[!grepl("Negative", rownames(markers)),]
    if(nrow(markers) != 0)
    {
      markers$cluster <- i
    }
    file.name <- paste0(i, "-markers.csv")
    write.csv(markers, file.name)
    print(paste("Finished", i))
  }
  else
  {
    print(paste("Skipping due to low cell count:", i))
  }
}


# Niche Colors
niche_cols <- c(
  `niche1` = '#0000ff',
  `niche2` = '#ff4500',
  `niche3` = '#ffff00',
  `niche4` = '#c71585',
  `niche5` = '#1e90ff',
  `niche6` = '#ffdab9',
  `niche7` = '#00fa9a'
)

# UMAP Grouped by Niches
DimPlot(data, split.by = "nclust", label = T, cols = ccolss, raster=F, label.size = 7, repel=T)


# Spatial UMAP grouped by niches
ImageDimPlot(data, group.by = "nclust", fov = NULL, border.size = 0.03, cols = niche_cols, border.color = "black", mols.size = 1.5)

# Niche Bar Plot grouped by clusters
ggplot(data@meta.data, aes(x = nclust, fill = custom.ident)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = ccolss) +
  labs(x = '', y = 'Percent cells', fill = 'Cell type', title = 'Niche compositions across all data') + theme(axis.text.x = element_text(size = 11))+ theme(axis.text.y = element_text(size = 11))


# Niche Bar Plot grouped by FOVs
data$custom.ident <- data@active.ident
md <- data@meta.data
ggplot(md, aes(x = nclust, fill = custom.ident)) +
  geom_bar(position = 'fill') +
  facet_wrap(tissue~fov, scale = 'free_x', nrow = 2) +
  scale_fill_manual(values = ccolss) +
  labs(x = '', y = 'Percent cells', fill = 'Cell type', title = 'Niche compositions within each FOV') +
  theme(axis.text.x = element_text(angle = 90)) + theme(axis.text.x = element_text(size = 11))+ theme(axis.text.y = element_text(size = 9))

#### scType (adapted from https://github.com/IanevskiAleksandr/sc-type)
lapply(c("dplyr","HGNChelper","openxlsx"), library, character.only = T)

db_ <- "/home/u8/chakraj/SMI-PDAC/ScTypeDB_PDAC.xlsx"
tissue <- "Pancreas"

# scType Load
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = as.matrix(data[["SCT"]]@data), scaled = F,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
data$custom.ident <- data@active.ident
cL_resutls = do.call("rbind", lapply(unique(data@meta.data$custom.ident), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$custom.ident==cl, , drop = F]), drop = F]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$custom.ident==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay onto UMAP
data@meta.data$sctype_classes <- ""
for(j in unique(sctype_scores$cluster)){
  cl_type <- sctype_scores[sctype_scores$cluster==j,];
  data@meta.data$sctype_classes[data@meta.data$custom.ident == j] <- as.character(cl_type$type[1])
}

DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classes', label.size = 4)
write.xlsx(cL_resutls, "scType-Scores-SMI-PDAC.xlsx")

# Visualize Bubbles (alternative visualization for scType)
# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();

for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

# Remove duplicates
nodes <- distinct(nodes)

mygraph <- graph_from_data_frame(edges, vertices=nodes)

#Make the graph (with cluster labels)
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
gggr

# Make the graph (without cluster labels, with edited font)
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.25)))
gggr

# Multiplot with UMAP
scater::multiplot(DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)

# BubblePlot for specific clusters
# Keep 1 specific cluster and plot
cl <- 0
nodes_regex <- paste0("_", cl, "$|cluster ", cl, "$")
nodes_subset <- nodes[grepl(nodes_regex, nodes$cluster),]
edges_regex <- paste0("cluster ", cl)
edges_subset <- edges[edges$from %in% edges_regex,]
mygraph <- graph_from_data_frame(edges_subset, vertices=nodes_subset)

gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*5)))
gggr

# For Loop: Plot every cluster individually
for(cl in levels(data$custom.ident)) {
  nodes_regex <- paste0("_", cl, "$|cluster ", cl, "$")
  nodes_subset <- nodes[grepl(nodes_regex, nodes$cluster),]
  edges_regex <- paste0("cluster ", cl)
  edges_subset <- edges[edges$from %in% edges_regex,]
  mygraph <- graph_from_data_frame(edges_subset, vertices=nodes_subset)
  
  gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
    geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
    theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*5)))
  print(gggr)
}








