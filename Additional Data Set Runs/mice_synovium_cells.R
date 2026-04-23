# This code is adapted from the CellMentor demo and performs
# dimensionality reduction on mice synovium cells from this dataset:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200898

# Set seed for reproducibility
set.seed(100)

# Load in necessary libraries
library(CellMentor)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
library(scRNAseq)

# Load datasets (Control and MMS)
cont_matrix <- Read10X(data.dir = "/projectnb/ds596/students/kbeaty26/cont/")
mms_matrix  <- Read10X(data.dir = "/projectnb/ds596/students/kbeaty26/MMS/")

# Create reference labels with clustering
seu_cont <- CreateSeuratObject(counts = cont_matrix, min.cells = 3, min.features = 200)
seu_cont <- NormalizeData(seu_cont)
seu_cont <- FindVariableFeatures(seu_cont)
seu_cont <- ScaleData(seu_cont)
seu_cont <- RunPCA(seu_cont)
seu_cont <- FindNeighbors(seu_cont, dims = 1:20)
seu_cont <- FindClusters(seu_cont, resolution = 0.3)
seu_cont <- RunUMAP(seu_cont, dims = 1:20)

# Extract cluster labels as reference cell types for CellMentor
reference_celltypes <- as.character(Idents(seu_cont))
names(reference_celltypes) <- colnames(seu_cont)
reference_matrix <- as(GetAssayData(seu_cont, layer = "counts"), "dgCMatrix")
query_matrix <- as(mms_matrix, "dgCMatrix")

# Align genes
shared_genes <- intersect(rownames(reference_matrix), rownames(query_matrix))
cat("Shared genes between cont and MMS:", length(shared_genes), "\n")

reference_matrix <- reference_matrix[shared_genes, ]
query_matrix     <- query_matrix[shared_genes, ]

# Create CellMentor object
csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix    = reference_matrix,
  ref_celltype  = reference_celltypes,
  data_matrix   = query_matrix,
  norm          = TRUE,
  most.variable = TRUE,
  scale         = TRUE,
  scale_by      = "cells",
  verbose       = TRUE,
  num_cores     = 1
)

# Perform parameter search
optimal_params <- CellMentor(
  csfnmf_obj,
  alpha_range = c(1, 5),
  beta_range  = c(1, 5),
  gamma_range = c(0.1),
  delta_range = c(1),
  num_cores   = 5,
  verbose     = TRUE
)

best_model <- optimal_params$best_model
K_VALUE    <- best_model@parameters$rank
print(optimal_params$best_params)

# Project MMS query onto learned space
h_project <- project_data(
  W         = best_model@W,
  X         = best_model@matrices@data,
  num_cores = 5,
  verbose   = TRUE
)

# Predict cell types using CellMentor latent space
H_ref    <- as.matrix(best_model@H)
ref_labs <- best_model@annotation$celltype
names(ref_labs) <- rownames(best_model@annotation)
H_query  <- as.matrix(h_project)

cosine_sim <- function(A, B) {
  A_norm <- sweep(A, 2, sqrt(colSums(A^2)), "/")
  B_norm <- sweep(B, 2, sqrt(colSums(B^2)), "/")
  t(A_norm) %*% B_norm
}

sim_matrix       <- cosine_sim(H_query, H_ref)
best_match_idx   <- apply(sim_matrix, 1, which.max)
predicted_labels <- unname(as.character(ref_labs[best_match_idx]))

# Create Seurat object
rownames(query_matrix) <- make.unique(rownames(query_matrix))
seu_mms <- CreateSeuratObject(counts = query_matrix)
seu_mms$condition          <- "MMS"
seu_mms$predicted_celltype <- predicted_labels

seu_mms[["CellMentor"]] <- CreateDimReducObject(
  embeddings = t(as.matrix(h_project)),
  key        = "CellMentor_",
  assay      = DefaultAssay(seu_mms),
  loadings   = as.matrix(best_model@W)
)

seu_mms <- RunUMAP(seu_mms, reduction = "CellMentor", dims = 1:K_VALUE)

# Visualizations
# Plot 1: MMS predicted cell types
DimPlot(seu_mms,
        group.by = "predicted_celltype",
        label    = TRUE,
        repel    = TRUE) +
  ggtitle("GSE200898 MMS Synovium — Predicted Cell Type\n(Control Synovium as Reference)") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 2: Control reference clusters for comparison
DimPlot(seu_cont,
        label = TRUE,
        repel = TRUE) +
  ggtitle("GSE200898 Control Synovium — Seurat Clusters\n(Used as CellMentor Reference)") +
  theme(plot.title = element_text(hjust = 0.5))

# Find top expressed genes for clusters
all_markers <- FindAllMarkers(
  seu_cont,
  only.pos    = TRUE,   
  min.pct     = 0.25, 
  logfc.threshold = 0.25
)

# Print top 10 genes per cluster ordered by fold change
top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

top10 %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = Inf)
