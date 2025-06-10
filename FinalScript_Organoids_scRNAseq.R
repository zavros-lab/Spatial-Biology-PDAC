#Zavros scRNAseq Organoid - Clean script
#For GitHub repository

# 1. Data loading and merging ---------------------------------------------

##Libraries
install.packages("remotes")
remotes::install_github("LTLA/scuttle")
library("Seurat")
library("dplyr")
library("EnhancedVolcano")
library("clusterProfiler")
library("org.Hs.eg.db")
library("scCustomize")
library("dplyr")
library("scuttle")
library("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install('glmGamPoi')
library("scater")
install.packages("RColorBrewer")
library("RColorBrewer")
BiocManager::install('multtest')
install.packages('metap')
devtools::install_github('immunogenomics/presto')

# Define the folder containing the .rds files
folder_path <- "~/Datasets/Zavros/Organoids"

# Get a list of all .rds files in the folder
rds_files <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)

# Function to read all Seurat objects
read_seurat_objects <- function(file_list) {
  seurat_list <- lapply(file_list, readRDS)
  names(seurat_list) <- basename(file_list)  # Assign file names as list names
  return(seurat_list)
}

# Load all Seurat objects into a list
seurat_objects <- read_seurat_objects(rds_files)

# Check the name of each object
print(names(seurat_objects))

# Extract the first Seurat object
first_seurat <- seurat_objects[[1]]
remaining_seurats <- seurat_objects[-1]

# Merge all objects while preserving sample IDs
merged_seurat <- merge(
  x = first_seurat, 
  y = remaining_seurats, 
  add.cell.ids = names(seurat_objects), 
  project = "Zavros_Organoids"
)

# Check if cell names have been correctly prefixed
head(colnames(merged_seurat))
tail(colnames(merged_seurat))

# Confirm unique prefixes
unique(sapply(X = strsplit(colnames(merged_seurat), split = "_"), FUN = "[", 1))

# Check sample distribution
table(merged_seurat$orig.ident)

# Check the merged object
print(merged_seurat)

rm(seurat_objects)
gc()

# 2. Doublet detection ----------------------------------------------------

library(DropletUtils)

#Conversion Seurat to SCE
org.sce <-as.SingleCellExperiment(merged_seurat)

#Check if samples names are present, and create a metadata column called "Sample"
head(colData(org.sce))
org.sce$Sample <- org.sce$ident

# Doublet detection
## Preprocessing
### QC
#Discard only very low quality cells. According to `scDblFinder` vignette discarding cells with < 200 UMIs is reasonable as an initial qc (more [here](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html#usage)). 

qc.dbl.df <- perCellQCMetrics(org.sce)

qc.dbl.sum <- qc.dbl.df$sum < 200



### Discard low quality cells.

fltrd.org.sce <- org.sce[, !qc.dbl.sum]

table(nexprs(fltrd.org.sce, byrow = T) > 0)

set.seed(1234)
non.zero.features <- nexprs(fltrd.org.sce, byrow = T) > 0
table(non.zero.features)

fltrd.org.sce <- fltrd.org.sce[non.zero.features, ]

### Normalisation
set.seed(1234)
dbl.clusters.samples <- scran::quickCluster(fltrd.org.sce,
                                            block=fltrd.org.sce$Sample
)                                                                             

library(scuttle)
library(scran)
library(DropletUtils)
fltrd.org.sce <-computeSumFactors(fltrd.org.sce,
                                  clusters = dbl.clusters.samples)

print(sizeFactors(fltrd.org.sce) %>% head(n = 10))

fltrd.org.sce <-
  logNormCounts(fltrd.org.sce)

### Variance modelling

set.seed(1234)
dec.dbl <-
  modelGeneVarByPoisson(fltrd.org.sce)
top.dec.dbl <-
  getTopHVGs(dec.dbl, n = 5000)

### Dimensionality reduction.

set.seed(1234)
fltrd.org.sce <-
  denoisePCA(fltrd.org.sce,
             technical = dec.dbl,
             subset.row = top.dec.dbl) 

reducedDimNames(fltrd.org.sce)

sample.no <- length(unique(fltrd.org.sce$Sample))
sample.colours  <- colorRampPalette(brewer.pal(11, "Spectral"))

plotReducedDim(fltrd.org.sce, colour_by = "Sample", dimred = "PCA") + scale_color_manual(values = sample.colours(sample.no))

set.seed(1234)
fltrd.org.sce <- runUMAP(
  fltrd.org.sce,
  dimred = "PCA"#,
)

reducedDimNames(fltrd.org.sce)

plotReducedDim(fltrd.org.sce,
               colour_by = "Sample",
               dimred = "UMAP"#,
) + scale_color_manual(values = sample.colours(sample.no))

plotReducedDim(fltrd.org.sce,
               colour_by = "MUC4",
               dimred = "UMAP",
               text_by = "Sample",
               # text_col = "black"
)


### Clustering
require('scDblFinder')
set.seed(1234)

fltrd.org.sce$fastClusters <-
  fastcluster(fltrd.org.sce, rdname = "PCA") %>% as.factor()

colLabels(fltrd.org.sce) <- fltrd.org.sce$fastClusters

table(colLabels(fltrd.org.sce))

fltrd.org.sce$Sample <- factor(fltrd.org.sce$Sample, levels = fltrd.org.sce$Sample %>% unique())

plotUMAP(fltrd.org.sce, colour_by = "label") +
  plotUMAP(fltrd.org.sce, colour_by = "Sample")  + scale_color_manual(values = sample.colours(sample.no))

## scDblFinder
### Run scDblFinder

fltrd.org.sce <-
  scDblFinder(
    fltrd.org.sce,
    samples = "Sample",
    clusters = "fastClusters",
  )

table(fltrd.org.sce$scDblFinder.class)

plotUMAP(fltrd.org.sce, colour_by="scDblFinder.class") +
  plotUMAP(fltrd.org.sce, colour_by="Sample") #+ scale_color_manual(values = sample.colours(sample.no))

saveRDS(fltrd.org.sce, file = "~/Datasets/Zavros/fltrd.org.sce.RDS")


### Doublet detection metrics.

dbl.class.table.list <-
  lapply(fltrd.org.sce$Sample %>% unique(), function(x)
    fltrd.org.sce[, fltrd.org.sce$Sample == x, ]$scDblFinder.class %>% tibble::as_tibble() %>%
      group_by(value) %>% tally()
  )

names(dbl.class.table.list) <- fltrd.org.sce$Sample %>% unique()

library(purrr)
dbl.tibble <- map_dfr(dbl.class.table.list, bind_rows, .id = "Sample")
dbl.tibble <- dbl.tibble %>% group_by(Sample) %>% mutate( percentage = 100 * n /sum(n)) %>% ungroup()

dbl.tibble %>%
  filter(value == "doublet") %>%
  arrange(desc(percentage)) %>% 
  print(n=12)

print(as.data.frame(colData(fltrd.org.sce)) %>%
        dplyr::count(Sample) %>%
        dplyr::arrange(desc(n)))

# Discard doublets

dbl.fltrd.org.sce <-
  fltrd.org.sce[, fltrd.org.sce$scDblFinder.class == "singlet"]

print(as.data.frame(colData(dbl.fltrd.org.sce)) %>%
        dplyr::count(scDblFinder.class) %>%
        dplyr::arrange(desc(n)))

saveRDS(dbl.fltrd.org.sce,
        file = "~/Datasets/Zavros/dbl.fltrd.org.sce.RDS")

# 3. QC -------------------------------------------------------------------

# Load SCE object with doublets removed.
dbl.fltrd.org.sce <- readRDS(file = "~/Datasets/Zavros/dbl.fltrd.org.sce.RDS")
rowData(dbl.fltrd.org.sce) <- DataFrame(Symbol = rownames(dbl.fltrd.org.sce))
head(rowData(dbl.fltrd.org.sce))

colnames(rowData(dbl.fltrd.org.sce))
colnames(dbl.fltrd.org.sce)
diet.dbl.fltrd.org.sce <-
  SingleCellExperiment(
    assays = list(counts = counts(dbl.fltrd.org.sce)),
    colData = colData(dbl.fltrd.org.sce),
    rowData = rowData(dbl.fltrd.org.sce)
  )

identical(counts(dbl.fltrd.org.sce),
          counts(diet.dbl.fltrd.org.sce))

identical(colData(dbl.fltrd.org.sce),
          colData(diet.dbl.fltrd.org.sce))

identical(rowData(dbl.fltrd.org.sce),
          rowData(diet.dbl.fltrd.org.sce))

rm(dbl.fltrd.org.sce)
gc()

diet.dbl.fltrd.org.sce$sizeFactor <- NULL
diet.dbl.fltrd.org.sce$fastClusters <- NULL
diet.dbl.fltrd.org.sce$label <- NULL

# QC

mito.genes <- grep("^MT-", rowData(diet.dbl.fltrd.org.sce)$Symbol)
ribo.genes <- grep("^RP[SL]", rowData(diet.dbl.fltrd.org.sce)$Symbol)
mito.genes.names <- grep("^MT-", rowData(diet.dbl.fltrd.org.sce)$Symbol, value = T)
ribo.genes.names <- grep("^RP[SL]", rowData(diet.dbl.fltrd.org.sce)$Symbol, value = T)


cell.qc.df <-
  perCellQCMetrics(diet.dbl.fltrd.org.sce,
                   subsets = list(mito = mito.genes,
                                  ribo = ribo.genes))

colData(diet.dbl.fltrd.org.sce) <- cbind(colData(diet.dbl.fltrd.org.sce), cell.qc.df)
colData(diet.dbl.fltrd.org.sce)$aggregation_order <- gsub("_.*", "", rownames(colData(diet.dbl.fltrd.org.sce)))


## Fixed thresholds

fixed.mito <- diet.dbl.fltrd.org.sce$subsets_mito_percent > 30
table(fixed.mito)

fixed.ribo <-  diet.dbl.fltrd.org.sce$subsets_ribo_percent > 20
table(fixed.ribo)

plotColData(diet.dbl.fltrd.org.sce, x="Sample", y="subsets_mito_percent",
            colour_by=I(fixed.mito)) +
  plotColData(diet.dbl.fltrd.org.sce, x="Sample", y="subsets_ribo_percent",
              colour_by=I(fixed.ribo))

fixed.detected <- diet.dbl.fltrd.org.sce$detected < 250 
table(fixed.detected)

plotColData(diet.dbl.fltrd.org.sce, x="Sample", y="sum")

fixed.sum <- diet.dbl.fltrd.org.sce$sum < 500 
table(fixed.sum)

plotColData(diet.dbl.fltrd.org.sce, x="Sample", y="sum",
            colour_by=I(fixed.sum)) +
  plotColData(diet.dbl.fltrd.org.sce, x="Sample", y="detected",
              colour_by=I(fixed.detected))

table(fixed.mito | fixed.ribo |fixed.detected | fixed.sum)


# Discard low quality cells.

cells.discard <- fixed.mito | fixed.ribo | fixed.detected | fixed.sum 
table(cells.discard)

qc.dbl.fltrd.org.sce <- diet.dbl.fltrd.org.sce[, !cells.discard]

plotColData(
  qc.dbl.fltrd.org.sce,
  x = "Sample",
  y = "detected"
)

colData(qc.dbl.fltrd.org.sce) %>% as.data.frame() %>% group_by(Sample) %>%
  summarise(median = median(sum, na.rm = TRUE)) %>% arrange(median)

summary(qc.dbl.fltrd.org.sce$subsets_ribo_percent)


# Discard  mitochondrial and ribosomal genes.

mito.genes.names <- grep("^MT-", rowData(qc.dbl.fltrd.org.sce)$Symbol, value = T)
ribo.genes.names <- grep("^RP[SL]", rowData(qc.dbl.fltrd.org.sce)$Symbol, value = T)

no.mito.ribo.genes <- !rowData(qc.dbl.fltrd.org.sce)$Symbol %in% c(mito.genes.names, ribo.genes.names)

table(no.mito.ribo.genes)

qc.dbl.fltrd.org.sce <-
  qc.dbl.fltrd.org.sce[no.mito.ribo.genes 
  ]

dim(qc.dbl.fltrd.org.sce)

saveRDS(qc.dbl.fltrd.org.sce,
        file ="~/Datasets/Zavros/final_qc.dbl.fltrd.org.sce.RDS")


# 4. Seurat Pre-Processing -------------------------------------------------------

# Load QCed SCE object.

qc.dbl.fltrd.org.sce <- readRDS("~/Datasets/Zavros/final_qc.dbl.fltrd.org.sce.RDS")
library("RColorBrewer")

sample.no <- length(unique(qc.dbl.fltrd.org.sce$Sample))
sample.colours  <- colorRampPalette(brewer.pal(11, "Spectral"))

# Create a Seurat object (QCed with Scater)
rownames(qc.dbl.fltrd.org.sce) <-
  uniquifyFeatureNames(rownames(qc.dbl.fltrd.org.sce),
                       rowData(qc.dbl.fltrd.org.sce)$Symbol)

rowData(qc.dbl.fltrd.org.sce) %>% as_tibble (rownames = "rownames") %>% dplyr::filter(duplicated(Symbol) | duplicated(Symbol, fromLast = TRUE))

org.seurat <-
  CreateSeuratObject(
    counts = counts(qc.dbl.fltrd.org.sce),
    meta.data = colData(qc.dbl.fltrd.org.sce) %>% as.data.frame()
  )
rm(qc.dbl.fltrd.org.sce)
gc()

#Removal of 1 sample from the clinical trial
org.seurat <- subset(org.seurat, subset = Sample != "PDAC-008-V2B")
table(org.seurat$Sample)
org.seurat$Sample <- droplevels(org.seurat$Sample)
table(org.seurat$Sample)

# 5. Statistics Seurat object ---------------------------------------------


#STATISTICS::

ncol(org.seurat)
meta_data <- org.seurat@meta.data

# Add nFeature_RNA (gene count per cell) if not already present
if (!"nFeature_RNA" %in% colnames(meta_data)) {
  meta_data$nFeature_RNA <- org.seurat[["RNA"]]@meta.features$nFeature_RNA
}

# Summary stats
summary_df <- meta_data %>%
  group_by(Sample) %>%
  summarise(
    num_cells = n(),
    mean_genes_per_cell = mean(nFeature_RNA),
    median_genes_per_cell = median(nFeature_RNA)
  )

# Overall stats
median_cells_per_sample <- median(summary_df$num_cells)
mean_cells_per_sample <- mean(summary_df$num_cells)

cat("Median number of cells per organoid (sample):", median_cells_per_sample, "\n")
cat("Mean number of cells per organoid (sample):", mean_cells_per_sample, "\n")

mean_genes_all <- mean(org.seurat@meta.data$nFeature_RNA)
median_genes_all <- median(org.seurat@meta.data$nFeature_RNA)

cat("Mean number of genes per cell (entire dataset):", mean_genes_all, "\n")
cat("Median number of genes per cell (entire dataset):", median_genes_all, "\n")

ggplot(summary_df, aes(x = Sample, y = num_cells)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Number of Cells per Sample",
    x = "Sample (Organoid)",
    y = "Number of Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(org.seurat@meta.data, aes(x = Sample, y = nFeature_RNA)) +
  geom_boxplot(fill = "orange", outlier.size = 0.5) +
  theme_minimal() +
  labs(
    title = "Distribution of Genes per Cell per Sample",
    x = "Sample (Organoid)",
    y = "Number of Genes per Cell"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 6. Normalisation SCTransform and RPCA integration (Integrate Layers v5) ------------------------------------------


#NORMALISATION WITH SCTRANSFORM AND INTEGRATION AS MENTIONED HERE: https://satijalab.org/seurat/articles/integration_introduction.html#perform-integration-with-sctransform-normalized-datasets
# SCTransform 

#Normalisation using SCTransform

# split datasets and process without integration
org.seurat[["RNA"]] <- split(org.seurat[["RNA"]], f = org.seurat$Sample)
oopts <- options(future.globals.maxSize = 6.0 * 1e9)  ## 6.0 GB
org.seurat <- SCTransform(org.seurat, vst.flavor = "v2",  method = "glmGamPoi", vars.to.regress="subsets_mito_percent", ncells=ncol(org.seurat),
                          verbose = TRUE)
org.seurat <- RunPCA(org.seurat)
org.seurat <- RunUMAP(org.seurat, dims = 1:30)
DimPlot(org.seurat, reduction = "umap", group.by = c("Sample"))

# integrate datasets
options(future.globals.maxSize = 40.0 * 1e9) 
org.seurat <- IntegrateLayers(
  object = org.seurat,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = FALSE
)

ElbowPlot(org.seurat, ndims = 50)

org.seurat <- FindNeighbors(org.seurat, reduction = "integrated.dr", dims = 1:30)
org.seurat <- FindClusters(org.seurat, resolution = 0.7)

org.seurat <- RunUMAP(org.seurat, dims = 1:30, reduction = "integrated.dr")

DimPlot(org.seurat, reduction = "umap", group.by = "Sample")
DimPlot(org.seurat, reduction = "umap", group.by = c("seurat_clusters"))
FeaturePlot(org.seurat,features="GAFP")
saveRDS(org.seurat, file = "~/Datasets/Zavros/SCTransform_RPCAInt_Org.seurat.rds")


# 7. Cell type annotation -------------------------------------------------

# Set identities
Idents(org.seurat) <- "seurat_clusters"

# SCT assay prep
org.seurat <- PrepSCTFindMarkers(org.seurat)

#2022 publication +in-house analysis
VlnPlot(org.seurat,features=c("KRT19","MUC1","KRT17","ACTN4","LMO7"),ncol=1) #Tumor cells in-house + 2022
VlnPlot(org.seurat,features=c("COL1A1","CALD1","COL6A3","COL1A2","FN1","POSTN","DCN")) #CAF
VlnPlot(org.seurat,features=c("PDGFRB","FAP")) #CAF
VlnPlot(org.seurat,features=c("CHAT","TAC1","CHGA","SST")) #neurons

VlnPlot(org.seurat,features=c("CPB1","CELA3A","CTRC","PNLIP","PDIA2")) #Acinar
VlnPlot(org.seurat,features=c("CFTR")) #ADM

VlnPlot(org.seurat,features=c("KRT5")) #Malignant?

FeaturePlot(org.seurat,features=c("C3","CALD1","MMP11"))

VlnPlot(org.seurat,features=c("MMP11","CCN2","POSTN","ACTA2")) #myCAF from this publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC9327514/, CCN2=CTGF=Connective tissue growth factor 
VlnPlot(org.seurat,features=c("IL6","CXCL1","CXCL12","CCL2","PDGFRA","HAS1")) #iCAF from this publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC9327514/
VlnPlot(org.seurat,features=c("H2AB1","CD74","SAA3")) #apCAF from this publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC9327514/

#2023 Publication
VlnPlot(org.seurat,features=c("KRT19","EPCAM"),ncol=1) #Tumor cells 2023
VlnPlot(org.seurat,features=c("MKI67"),ncol=1) #2023
VlnPlot(org.seurat,features="DCN") #used 2023 to define Mes. cells
VlnPlot(org.seurat,features=c("ACTA2","MMP11","COL1OA1")) #myCAF
VlnPlot(org.seurat,features=c("C3","C7","CFD","PTGDS")) #iCAF
VlnPlot(org.seurat,features=c("ACTA2","RGS5","CSPG4","NOTCH3")) #pericytes
VlnPlot(org.seurat,features=c("SOX10","S100B","PLP1")) #schwann

#From the manuscript
VlnPlot(org.seurat,features=c("SOX9","CD44","TACSTD2"))
VlnPlot(org.seurat,features=c("EPCAM", "VCAM1", "ITGB1","CD24")) 
FeaturePlot(org.seurat,features="KRT19") 

#Other tests

VlnPlot(org.seurat,features=c("PROM1", "CD24", "CXCR4","ABCG2","ESA")) #Marker of stem in PDAC
VlnPlot(org.seurat,features=c("PTF1A", "BHLHA15", "NR5A2")) #Acinar : From manuscript
VlnPlot(org.seurat,features=c("SOX2", "HNF1B", "ONECUT1","PDX1","CA2","PROM1","SPP1")) #Ductal : From manuscript
VlnPlot(org.seurat,features=c("FAP", "SERPINA2", "SPINK7","VIM","MMP19","COL1A2","COL3A1")) #CAF : From manuscript

VlnPlot(org.seurat,features=c("PTF1A", "BHLHA15", "NR5A2","CPE",ncol=2)) #Acinar : From Yana paper

VlnPlot(org.seurat,features=c("SOX2", "HNF1B", "ONECUT1","PDX1","CA2","PROM1","SPP1")) #Ductal : From Yana paper

DotPlot(org.seurat,features=c("TP53","SMAD4","KRAS"))

#Cluster Annotation (manual analysis based on markers)
org.seurat$cell_type <- "Tumor cells"
org.seurat$cell_type[org.seurat$seurat_clusters %in% c(4, 9, 5, 13)] <- "Proliferating tumor cells"
org.seurat$cell_type[org.seurat$seurat_clusters == 19] <- "CAFs"
org.seurat$cell_type[org.seurat$seurat_clusters == 18] <- "CD44+ TROP2+ CSCs"
org.seurat$cell_type[org.seurat$seurat_clusters %in% c(3, 11, 14, 16, 21, 22)] <- "CSCs"

DimPlot(org.seurat,group.by="cell_type")

DimPlot(org.seurat, group.by = "cell_type", repel = TRUE) +
  ggtitle("Cell Type Annotation") +                        # Add title               # Custom colors
  labs(color = "Cell Type Annotation") +                  # Change legend title
  theme_classic()                                         # Clean theme (optional)

saveRDS(org.seurat, file = "~/Datasets/Zavros/Annotated_SCTransform_RPCAInt_Org.seurat.rds")

# 9.Cluster proportion analysis per sample ----------------------------------------------------------------------

# Tabulate number of cells per cluster per sample
df <- org.seurat@meta.data[, c("seurat_clusters", "Sample")]
cluster_counts <- df %>%
  group_by(Sample, seurat_clusters) %>%
  summarise(n = n()) %>%
  ungroup()

# Tile plot
ggplot(cluster_counts, aes(x = Sample, y = seurat_clusters, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Cluster Abundance per Organoid",
       x = "Organoid (Sample)",
       y = "Seurat Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(cluster_counts, aes(x = Sample, y = seurat_clusters, size = n)) +
  geom_point(color = "darkred", alpha = 0.7) +
  scale_size_continuous(name = "Cell Count") +
  labs(title = "Cluster Representation per Organoid (Dot Size = Cell Count)",
       x = "Organoid (Sample)",
       y = "Seurat Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compute proportions per Sample
cluster_props <- cluster_counts %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n))

# Stacked barplot of proportions
ggplot(cluster_props, aes(x = Sample, y = prop, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cluster Composition per Organoid (Proportions)",
       x = "Organoid (Sample)",
       y = "Proportion of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Differential expression--------------------------------------------------------------------

# Set identities
Idents(org.seurat) <- "seurat_clusters"

# SCT assay prep
org.seurat <- PrepSCTFindMarkers(org.seurat)

# Get unique clusters
clusters <- levels(org.seurat)

# Loop through each cluster and find markers
all_markers <- lapply(clusters, function(clust) {
  markers <- FindMarkers(org.seurat, assay = "SCT", ident.1 = clust, verbose = FALSE)
  markers$gene <- rownames(markers)
  top_markers <- head(markers[order(markers$p_val_adj), ], 200)
  return(top_markers)
})

# Name each list element by cluster
names(all_markers) <- paste0("Cluster_", clusters)

#View ex. cluster 0
head(all_markers[["Cluster_0"]], 200)

#Save
for (i in names(all_markers)) {
  write.csv(all_markers[[i]], file = paste0(i, "_top200_markers.csv"), row.names = FALSE)
}


# Function to search one or more genes in top markers
find_gene_in_clusters <- function(gene_list, marker_list) {
  gene_list <- toupper(gene_list)  # make case-insensitive
  result <- lapply(names(marker_list), function(cluster_name) {
    markers <- marker_list[[cluster_name]]
    hits <- markers[toupper(markers$gene) %in% gene_list, ]
    if (nrow(hits) > 0) {
      hits$Cluster <- cluster_name
      return(hits)
    } else {
      return(NULL)
    }
  })
  # Combine results into one data frame
  result_df <- do.call(rbind, result)
  if (is.null(result_df)) {
    message("No genes found in any cluster.")
    return(NULL)
  } else {
    return(result_df)
  }
}

find_gene_in_clusters(c("CD44", "SOX9", "TACSTD2"), all_markers)
find_gene_in_clusters(c("VIM", "ACTA2", "KRT19"), all_markers)

find_gene_in_clusters(c("FN1", "CALD1", "AMY1A"), all_markers)
find_gene_in_clusters(c("EPCAM", "CALD1", "AMY1A"), all_markers)
