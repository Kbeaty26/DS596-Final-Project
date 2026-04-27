
# Reproduction of the CellMentor tutorial on the human cancer fibroblast dataset.

## Package imports

library(Matrix)
library(CellMentor)
library(Seurat)
library(ggplot2)
library(dplyr)

## Extension import 

.libPaths(c("/projectnb/ds596/students/aliviap/R_libs_4.5.2", .libPaths()))
source("diagnostics.R")
source("plotting.R")

###

## Subset creation helper function

create_subset <- function(matrix, celltypes, cells_per_type = 5000) {
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

## Loading reference data (lung fibroblast)

LUNG_fibro_tumour <- readRDS("/projectnb/ds596/students/mivargas/scRNA-seq/LUNG_fibro_tumour.rds")
lung_matrix <- GetAssayData(LUNG_fibro_tumour, assay = "SCT", layer = "data")
lung_fibro_types <- as.character(Idents(LUNG_fibro_tumour))
names(lung_fibro_types) <- colnames(LUNG_fibro_tumour)
lung_ref_subset <- create_subset(lung_matrix, lung_fibro_types, cells_per_type = 5000)

## Loading query data (breast fibroblast)

BREAST_fibro_tumour <- readRDS("/projectnb/ds596/students/mivargas/scRNA-seq/BREAST_fibro_tumour.rds")
breast_matrix <- GetAssayData(BREAST_fibro_tumour, assay = "RNA", layer = "data")
breast_fibro_types <- as.character(BREAST_fibro_tumour$CAFtype)
names(breast_fibro_types) <- colnames(BREAST_fibro_tumour)
breast_q_subset <- create_subset(breast_matrix, breast_fibro_types, cells_per_type = 5000)

## Updated names- run after subsetting

reference_matrix <- as.matrix(lung_ref_subset$matrix)
reference_fibro_types <- lung_ref_subset$celltypes
query_matrix <- as.matrix(breast_q_subset$matrix)
query_fibro_types <- breast_q_subset$celltypes

## CellMentor object creation

csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix = reference_matrix,
  ref_celltype = reference_fibro_types,
  data_matrix = query_matrix,
  norm = TRUE,
  most.variable = TRUE,
  scale = TRUE,
  scale_by = "cells",
  verbose = TRUE,
  num_cores = 10
)

## Find optimal parameters for CellMentor
optimal_params <- CellMentor(
  csfnmf_obj,
  alpha_range = c(1, 5),
  beta_range = c(3, 4),
  gamma_range = c(0.1),
  delta_range = c(1),
  num_cores = 10,
  verbose = TRUE
)

## Get the best model from the optimal parameters
best_model <- optimal_params$best_model
K_VALUE <- best_model@parameters$rank
print(optimal_params$best_params)

## Project query data onto the learned space
h_project <- project_data(
  W = best_model@W, # Learned gene weights
  X = best_model@matrices@data,  # Query data matrix
  num_cores = 10,
  verbose = TRUE
)

## Extension begins here 

### Extract reference matrix and embedding from the trained model
X_ref  <- best_model@matrices@ref    
H_ref  <- best_model@H               
X_query <- best_model@matrices@data   

### Reference cell type labels

ref_ct <- tryCatch(
  best_model@annotation$celltype,
  error = function(e) reference_fibro_types[colnames(X_ref)]
)

## I ran this before running the extension, then went back and ran extension to compare.

rownames(query_matrix) <- make.unique(rownames(query_matrix))
seu_breast <- CreateSeuratObject(counts = query_matrix)
seu_breast$celltype <- query_fibro_types

### Add CellMentor dimensionality reduction to Seurat object
seu_breast[["CellMentor"]] <- CreateDimReducObject(
  embeddings = t(as.matrix(h_project)),
  key = "CellMentor_",
  assay = DefaultAssay(seu_breast),
  loadings = as.matrix(best_model@W)
)

### Diagnostic run
diag_out <- run_projection_diagnostics(
  W                 = best_model@W,
  H_ref             = H_ref,
  X_ref             = X_ref,
  ref_celltypes     = ref_ct,
  H_query           = h_project,
  X_query           = X_query,
  novelty_threshold = 0.95,    # flag cells above 95th percentile of reference error
  novelty_method    = "per_class",  # threshold calibrated per cell type
  distance_metric   = "cosine",
  verbose           = TRUE
)

diagnostics <- diag_out$diagnostics

### Attaching true labels

diagnostics$true_label <- query_fibro_types[diagnostics$cell]
write.csv(diagnostics, "cancerfibro_diag5000.csv", row.names = FALSE)

## Visualizations (Normal UMAP)
seu_breast <- RunUMAP(seu_breast, reduction = 'CellMentor', dims= 1:K_VALUE)
DimPlot(seu_breast, group.by = 'celltype', label = TRUE, repel = TRUE) + ggtitle("Clustering of Fibroblast Types in Human Breast Cancer")

## Diagnostic visualization

p1 <- plot_error_distributions(
  ref_errors   = diag_out$ref_errors,
  query_errors = diagnostics$reconstruction_error,
  threshold    = quantile(diag_out$ref_errors, 0.95)
)
ggsave("error_distributions_fibroblast5000.png", p1, width = 6, height = 4, dpi = 150)

p2 <- plot_confidence_distribution(diagnostics)
ggsave("confidence_distribution_fibroblast5000.png", p2, width = 6, height = 4, dpi = 150)

## Seurat novelty plot

umap_coords <- Seurat::Embeddings(seu_breast, "umap")
diagnostics_ordered <- diagnostics[match(rownames(umap_coords), diagnostics$cell), ]

panels <- plot_umap_diagnostics(
  umap_coords = umap_coords,
  diagnostics = diagnostics_ordered,
  true_labels = diagnostics_ordered$true_label
)
ggsave("umap_novelty5000.png", panels$novel, width = 6, height = 5, dpi = 150)

## ARI/NMI run

library(aricode)
library(tidyr)

evaluate_clustering <- function(seu, K, resolution = 0.1) {
  seu <- FindNeighbors(seu, reduction = "CellMentor", dims = 1:K, verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  clusters <- as.integer(seu$seurat_clusters)
  truth    <- as.integer(factor(seu$celltype))
  ari <- ARI(clusters, truth)
  nmi <- NMI(clusters, truth)
  n_clusters <- length(unique(clusters))
  n_types    <- length(unique(seu$celltype))
  print(sprintf(
    "resolution=%.2f | clusters=%d (true types=%d) | ARI=%.3f | NMI=%.3f",
    resolution, n_clusters, n_types, ari, nmi
  ))
  list(seu = seu, ari = ari, nmi = nmi,
       n_clusters = n_clusters, resolution = resolution)
}

eval_res <- evaluate_clustering(seu_breast, K_VALUE, resolution = 0.1)

resolutions  <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
sweep_results <- lapply(resolutions, function(r) {
  evaluate_clustering(seu_breast, K_VALUE, resolution = r)
})

sweep_df <- data.frame(
  resolution = sapply(sweep_results, `[[`, "resolution"),
  n_clusters = sapply(sweep_results, `[[`, "n_clusters"),
  ARI        = sapply(sweep_results, `[[`, "ari"),
  NMI        = sapply(sweep_results, `[[`, "nmi")
)

print("Resolution sweep summary:")
print(sweep_df)

## Plot sweep results
sweep_long <- pivot_longer(sweep_df, cols = c("ARI", "NMI"),
                           names_to = "metric", values_to = "score")

ggplot(sweep_long, aes(x = resolution, y = score,
                       colour = metric, group = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.11, y = 0.05, label = "paper default",
           hjust = 0, colour = "grey40", size = 3.5) +
  scale_colour_manual(values = c("ARI" = "#2166ac", "NMI" = "#d6604d")) +
  labs(x = "Clustering resolution", y = "Score", colour = "Metric",
       title = "CellMentor clustering performance vs resolution\n(Breast CAF)") +
  theme_bw(base_size = 12)

## Plot UMAP of predicted labels

ggsave("umap_pred5000.png", panels$predicted, width = 6, height = 5, dpi = 150)

umap_truth <- DimPlot(seu_breast, group.by = 'celltype', label = TRUE, repel = TRUE, ) + ggtitle("Breast Cancer CAF: Ground Truth")
ggsave(filename = "umap_truth5000.png", plot = umap_truth, width = 6, height = 5, dpi = 150)
