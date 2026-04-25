# This code is adapted from the CellMentor demo and performs
#   dimensionality reduction on mice synovium cells from this dataset:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200898
# Strategy: Cluster control+MMS together in Seurat, annotate all clusters, 
#   then do a 70/30 stratified split for CellMentor evaluation

set.seed(100)

library(CellMentor)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
library(aricode)
library(tidyr)

# Load both datasets (control and MMS)
cont_matrix <- Read10X(data.dir = "/projectnb/ds596/students/kbeaty26/cont/")
mms_matrix  <- Read10X(data.dir = "/projectnb/ds596/students/kbeaty26/MMS/")

# Find shared genes before merging
shared_genes <- intersect(rownames(cont_matrix), rownames(mms_matrix))
cat("Shared genes between cont and MMS:", length(shared_genes), "\n")

cont_matrix <- cont_matrix[shared_genes, ]
mms_matrix  <- mms_matrix[shared_genes, ]

# Create individual Seurat objects with condition labels (Control and MMS)
seu_cont <- CreateSeuratObject(counts = cont_matrix, min.cells = 3, min.features = 200)
seu_mms  <- CreateSeuratObject(counts = mms_matrix,  min.cells = 3, min.features = 200)

seu_cont$condition <- "Control"
seu_mms$condition  <- "MMS"

# Merge and run joint Seurat clustering
# To ensure both conditions are clustered with the same labels
seu_merged <- merge(seu_cont, y = seu_mms,
                    add.cell.ids = c("Control", "MMS"),
                    project = "Synovium")
seu_merged <- JoinLayers(seu_merged)

seu_merged <- NormalizeData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- ScaleData(seu_merged)
seu_merged <- RunPCA(seu_merged)
seu_merged <- FindNeighbors(seu_merged, dims = 1:20)
seu_merged <- FindClusters(seu_merged, resolution = 0.2)
seu_merged <- RunUMAP(seu_merged, dims = 1:20)

# Clusters distributions across conditions
cat("\nCluster × Condition table:\n")
print(table(Idents(seu_merged), seu_merged$condition))

# Merged cluster visualization before annotation
# Used this plot to decide on cell type labels and gene markers
DimPlot(seu_merged,
        group.by = "seurat_clusters",
        label    = TRUE,
        repel    = TRUE,
        split.by = "condition") +
  ggtitle("Merged Synovium — Seurat Clusters by Condition\n(Before annotation)") +
  theme(plot.title = element_text(hjust = 0.5))

# Find markers to help with annotation
all_markers <- FindAllMarkers(
  seu_merged,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

cat("\nTop 10 marker genes per cluster:\n")
top10 %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = Inf)

# Annotate clusters
# Cross-reference top10 markers with the paper's Fig. 5a UMAP, Fig. 5b dot plot, and
#   existing control annotations.

celltype_names <- c(
  "0"  = "Penk_universal_fibroblast",
  "1"  = "Vascular_endothelial",
  "2"  = "Cx3cr1_resident_macrophage",
  "3"  = "Capillary_endothelial",
  "4"  = "Lrrc15_myofibroblast",
  "5"  = "Fap_sublining_fibroblast",
  "6"  = "Prg4_lining_fibroblast",
  "7"  = "MHCII_macrophage",
  "8"  = "Smooth_muscle_cell",
  "9"  = "Neutrophil",           # filtered out -- not in paper
  "10" = "Adipose_fibroblast",
  "11" = "Vascular_endothelial_2",
  "12" = "Proliferating_cell",
  "13" = "Lymphocyte",
  "14" = "Erythrocyte",          # filtered out -- not in paper
  "15" = "Schwann_cell"         # filtered out -- not in paper
)

seu_merged$celltype <- unname(celltype_names[as.character(Idents(seu_merged))])

cat("\nCell type counts across full merged dataset:\n")
print(sort(table(seu_merged$celltype), decreasing = TRUE))

# Annotated UMAP split by condition -- check MMS-specific clusters appear
DimPlot(seu_merged,
        group.by = "celltype",
        label    = TRUE,
        repel    = TRUE,
        split.by = "condition") +
  ggtitle("Merged Synovium — Annotated Cell Types by Condition") +
  theme(plot.title  = element_text(hjust = 0.5),
        legend.position = "bottom")


# Extract combined count matrix and labels
clusters_to_remove <- c("9","14","15")

seu_merged_filtered <- subset(seu_merged,
                              idents = clusters_to_remove,
                              invert = TRUE)

combined_counts    <- as(GetAssayData(seu_merged_filtered, layer = "counts"), "dgCMatrix")
combined_celltypes <- seu_merged_filtered$celltype
names(combined_celltypes) <- colnames(seu_merged_filtered)

cat("\nFinal cell type distribution before split:\n")
print(sort(table(combined_celltypes), decreasing = TRUE))


# 70/30 stratified split
create_split <- function(counts, celltypes, ref_frac = 0.70, seed = 100) {
  set.seed(seed)
  all_types <- unique(celltypes)
  ref_idx   <- c()
  qry_idx   <- c()
  
  for (ct in all_types) {
    idx    <- which(celltypes == ct)
    n_ref  <- max(1, round(length(idx) * ref_frac))
    chosen <- sample(idx, n_ref)
    ref_idx <- c(ref_idx, chosen)
    qry_idx <- c(qry_idx, setdiff(idx, chosen))
  }
  
  list(
    ref_counts    = counts[, ref_idx],
    ref_celltypes = setNames(celltypes[ref_idx], colnames(counts)[ref_idx]),
    qry_counts    = counts[, qry_idx],
    qry_celltypes = setNames(celltypes[qry_idx], colnames(counts)[qry_idx])
  )
}

splits <- create_split(combined_counts, combined_celltypes, ref_frac = 0.70)

ref_counts    <- splits$ref_counts
ref_celltypes <- splits$ref_celltypes
qry_counts    <- splits$qry_counts
qry_celltypes <- splits$qry_celltypes

print(sprintf("Reference: %d cells | Query: %d cells",
                ncol(ref_counts), ncol(qry_counts)))
print("Cell type distribution in reference:")
print(table(ref_celltypes))
print("Cell type distribution in query:")
print(table(qry_celltypes))

# Create CellMentor object and run parameter search
csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix    = ref_counts,
  ref_celltype  = ref_celltypes,
  data_matrix   = qry_counts,
  norm          = TRUE,
  most.variable = TRUE,
  scale         = TRUE,
  scale_by      = "cells",
  verbose       = TRUE,
  num_cores     = 10
)

optimal_params <- CellMentor(
  csfnmf_obj,
  alpha_range = c(1, 5),
  beta_range  = c(1, 5),
  gamma_range = c(0.1),
  delta_range = c(1),
  num_cores   = 10,
  verbose     = TRUE
)

best_model <- optimal_params$best_model
K_VALUE    <- best_model@parameters$rank

print("Best parameters:")
print(optimal_params$best_params)
print(sprintf("Optimal rank K = %d", K_VALUE))

# Project query data onto learned latent space
h_project <- project_data(
  W         = best_model@W,
  X         = best_model@matrices@data,
  num_cores = 10,
  verbose   = TRUE
)

# Build Seurat object for query and run UMAP
rownames(qry_counts) <- make.unique(rownames(qry_counts))
seu_qry <- CreateSeuratObject(counts = qry_counts)
seu_qry$celltype  <- qry_celltypes
seu_qry$condition <- seu_merged$condition[match(colnames(qry_counts),
                                                colnames(seu_merged))]

seu_qry[["CellMentor"]] <- CreateDimReducObject(
  embeddings = t(as.matrix(h_project)),
  key        = "CellMentor_",
  assay      = DefaultAssay(seu_qry),
  loadings   = as.matrix(best_model@W)
)

seu_qry <- RunUMAP(seu_qry, reduction = "CellMentor", dims = 1:K_VALUE)


# UMAP plots
# Plot 1: colored by ground truth cell type
DimPlot(seu_qry,
        group.by = "celltype",
        label    = TRUE,
        repel    = TRUE) +
  ggtitle("CellMentor UMAP — Query Cells\nColored by ground truth cell type") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 2: colored by condition to check Control vs MMS separation
DimPlot(seu_qry,
        group.by = "condition",
        cols     = c("Control" = "#2166ac", "MMS" = "#d6604d")) +
  ggtitle("CellMentor UMAP — Query Cells\nColored by condition") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 3: highlight MMS-specific cell types from the paper
mms_specific <- c("Lrrc15_myofibroblast", "Mmp9_inflammatory_macrophage")
seu_qry$is_mms_specific <- seu_qry$celltype %in% mms_specific

DimPlot(seu_qry,
        group.by = "is_mms_specific",
        cols     = c("TRUE" = "#E63946", "FALSE" = "grey80"),
        order    = "TRUE",
        pt.size  = 1.2) +
  ggtitle("MMS-specific subsets highlighted\n(Lrrc15+ myofibroblasts & Mmp9+ macrophages)") +
  theme_bw()

# Quantitative evaluation — ARI and NMI with resolution sweep
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

# Primary evaluation at resolution 0.1 (paper default)
print("\nClustering evaluation at resolution = 0.1")
eval_res <- evaluate_clustering(seu_qry, K_VALUE, resolution = 0.1)

# Resolution sweep
resolutions  <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
sweep_results <- lapply(resolutions, function(r) {
  evaluate_clustering(seu_qry, K_VALUE, resolution = r)
})

sweep_df <- data.frame(
  resolution = sapply(sweep_results, `[[`, "resolution"),
  n_clusters = sapply(sweep_results, `[[`, "n_clusters"),
  ARI        = sapply(sweep_results, `[[`, "ari"),
  NMI        = sapply(sweep_results, `[[`, "nmi")
)

print("\nResolution sweep summary:")
print(sweep_df)

# Plot sweep results
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
       title = "CellMentor clustering performance vs resolution\n(Combined Control+MMS Synovium)") +
  theme_bw(base_size = 12)