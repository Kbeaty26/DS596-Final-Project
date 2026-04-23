
# Reproduction of the CellMentor tutorial on the human fibroblast dataset.

library(Matrix)
library(CellMentor)
library(Seurat)

# Subset function

create_subset <- function(matrix, celltypes, cells_per_type = 30) {
  unique_types <- unique(celltypes)
  selected_cells <- c()
  for (cell_type in unique_types) {
    type_cells <- names(celltypes)[celltypes == cell_type]
    n_to_select <- min(cells_per_type, length(type_cells))
    selected <- sample(type_cells, n_to_select)
    selected_cells <- c(selected_cells, selected)
  }
  
  list(
    matrix = matrix[, selected_cells],
    celltypes = celltypes[selected_cells]
  )
}




###

# Loading reference data (lung fibroblast)

lung_matrix <- GetAssayData(LUNG_fibro_tumour, assay = "SCT", layer = "data")
lung_celltypes <- as.character(Idents(LUNG_fibro_tumour))
names(lung_celltypes) <- colnames(LUNG_fibro_tumour)

lung_ref_subset <- create_subset(lung_matrix, lung_celltypes, cells_per_type = 50)

# Loading query data (breast fibroblast)

breast_matrix <- GetAssayData(BREAST_fibro_tumour, assay = "RNA", layer = "data")
breast_celltypes <- as.character(Idents(BREAST_fibro_tumour))
names(breast_celltypes) <- colnames(BREAST_fibro_tumour)
breast_q_subset <- create_subset(breast_matrix, breast_celltypes, cells_per_type = 50)

# Updated names- run after subsetting

reference_matrix <- as.matrix(lung_ref_subset$matrix)
reference_celltypes <- lung_ref_subset$celltypes
query_matrix <- as.matrix(breast_q_subset$matrix)
query_celltypes <- breast_q_subset$celltypes

# CellMentor object creation

csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix = reference_matrix,
  ref_celltype = reference_celltypes,
  data_matrix = query_matrix,
  norm = TRUE,
  most.variable = TRUE,
  scale = TRUE,
  scale_by = "cells",
  verbose = TRUE,
  num_cores = 1
)

# Find optimal parameters
optimal_params <- CellMentor(
  csfnmf_obj,
  alpha_range = c(1, 5),      # Limited alpha range
  beta_range = c(1, 5),       # Limited beta range
  gamma_range = c(0.1),     # use only one gamma for speed
  delta_range = c(1),       # use only one delta for speed
  num_cores = 5,
  verbose = TRUE
)

# Get best model
best_model <- optimal_params$best_model
K_VALUE <- best_model@parameters$rank
print(optimal_params$best_params)

# Project query data onto the learned space
h_project <- project_data(
  W = best_model@W, # Learned gene weights
  X = best_model@matrices@data,  # Query data matrix
  num_cores = 5,
  verbose = TRUE
)

rownames(query_matrix) <- make.unique(rownames(query_matrix))
seu_breast <- CreateSeuratObject(counts = query_matrix)
seu_breast$celltype <- query_celltypes

# Add CellMentor dimensionality reduction to Seurat object
seu_breast$CellMentor <- CreateDimReducObject(
  embeddings = t(as.matrix(h_project)),
  key = "CellMentor_",
  assay = DefaultAssay(seu_breast),
  loadings = as.matrix(best_model@W)
)

# Visualization
seu_breast <- RunUMAP(seu_breast, reduction = 'CellMentor', dims= 1:K_VALUE)
DimPlot(seu_breast, group.by = 'celltype')
