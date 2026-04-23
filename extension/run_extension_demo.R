set.seed(100)

# SCC paths
R_LIBS   <- "/projectnb/ds596/students/aliviap/R_libs_4.5.2"
CM_LOCAL <- "/projectnb/ds596/students/aliviap/CellMentor_repro/scripts/CellMentor"
OUTDIR   <- "/projectnb/ds596/students/aliviap/extension/figures"

.libPaths(c(R_LIBS, .libPaths()))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Parse --subset flag
args     <- commandArgs(trailingOnly = TRUE)
SUBSET   <- "--subset" %in% args
NUM_CORES <- if (SUBSET) 1 else 4

cat("Mode:", if (SUBSET) "SUBSET (fast test)" else "FULL", "\n")

# Load packages
suppressPackageStartupMessages({
  if (!requireNamespace("CellMentor", quietly = TRUE)) {
    devtools::load_all(CM_LOCAL)
  } else {
    library(CellMentor)
  }
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# Load our extension
SCRIPT_DIR <- "/projectnb/ds596/students/aliviap/extension"
source(file.path(SCRIPT_DIR, "diagnostics.R"))
source(file.path(SCRIPT_DIR, "plotting.R"))


# 1. Load data
message("\n[1] Loading Baron and Muraro datasets from CellMentor package...")

baron  <- hBaronDataset()
muraro <- muraro_dataset()

reference_matrix    <- baron$data
reference_celltypes <- baron$celltypes
query_matrix        <- muraro$data
query_celltypes     <- muraro$celltypes   # known only for evaluation

message(sprintf("    Baron:  %d genes x %d cells, %d types",
                nrow(reference_matrix), ncol(reference_matrix),
                length(unique(reference_celltypes))))
message(sprintf("    Muraro: %d genes x %d cells, %d types",
                nrow(query_matrix), ncol(query_matrix),
                length(unique(query_celltypes))))


# 2. Optional subsample for fast testing 
if (SUBSET) {
  message("[2] Subsampling (30 cells/type) for fast test run...")
  subset_balanced <- function(mat, ct, n = 30) {
    sel <- unlist(lapply(unique(ct), function(t) {
      pool <- names(ct)[ct == t]
      sample(pool, min(n, length(pool)))
    }))
    list(matrix = mat[, sel], celltypes = ct[sel])
  }
  ref <- subset_balanced(reference_matrix, reference_celltypes, 30)
  qry <- subset_balanced(query_matrix,    query_celltypes,     30)
  reference_matrix    <- ref$matrix;  reference_celltypes <- ref$celltypes
  query_matrix        <- qry$matrix;  query_celltypes     <- qry$celltypes
} else {
  message("[2] Using full dataset...")
}


# 3. Build CellMentor object
message("\n[3] Building CSFNMF object...")

csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix    = reference_matrix,
  ref_celltype  = reference_celltypes,
  data_matrix   = query_matrix,
  norm          = TRUE,
  most.variable = TRUE,
  scale         = TRUE,
  scale_by      = "cells",
  verbose       = TRUE,
  num_cores     = NUM_CORES
)


# 4. Run CellMentor + parameter search
message("\n[4] Running CellMentor parameter search...")

optimal_params <- CellMentor(
  csfnmf_obj,
  alpha_range = c(1, 5),
  beta_range  = c(1, 5),
  gamma_range = c(0.1),
  delta_range = c(1),
  num_cores   = NUM_CORES,
  verbose     = TRUE
)

best_model <- optimal_params$best_model
K_VALUE    <- best_model@parameters$rank
message("Best rank K = ", K_VALUE)
message("Best params:")
print(optimal_params$best_params)


# 5. Project query
message("\n[5] Projecting query data...")

h_project <- project_data(
  W         = best_model@W,
  X         = best_model@matrices@data,
  num_cores = NUM_CORES,
  verbose   = TRUE
)


# 6. Run projection diagnostics (our extension)
message("\n[6] Running projection diagnostics (extension)...")

# Extract pieces from the trained model.
# The tryCatch() blocks handle any slot-name differences across CellMentor versions.
X_ref <- tryCatch(best_model@matrices@ref,
  error = function(e) stop(
    "Slot best_model@matrices@ref not found.\n",
    "Run: slotNames(best_model@matrices)\n",
    "and update this line to match.\n"))

H_ref <- tryCatch(best_model@H,
  error = function(e) stop(
    "Slot best_model@H not found.\n",
    "Run: slotNames(best_model)\n"))

ref_ct <- tryCatch(
  best_model@annotation$celltype,
  error = function(e) reference_celltypes[colnames(X_ref)]
)

X_query <- best_model@matrices@data

diag_out <- run_projection_diagnostics(
  W                 = best_model@W,
  H_ref             = H_ref,
  X_ref             = X_ref,
  ref_celltypes     = ref_ct,
  H_query           = h_project,
  X_query           = X_query,
  novelty_threshold = 0.95,
  novelty_method    = "per_class",
  distance_metric   = "cosine",
  verbose           = TRUE
)

diagnostics            <- diag_out$diagnostics
diagnostics$true_label <- query_celltypes[diagnostics$cell]

# Write CSV
out_csv <- file.path(OUTDIR, "baron_muraro_diagnostics.csv")
write.csv(diagnostics, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)

# Print summary
cat("\n==== Projection diagnostics summary ====\n")
cat("N query cells:       ", nrow(diagnostics), "\n")
cat("Mean recon error:    ", round(mean(diagnostics$reconstruction_error), 4), "\n")
cat("Mean confidence:     ", round(mean(diagnostics$confidence), 3), "\n")
cat("Flagged as novel:    ", sum(diagnostics$is_potential_novel),
    sprintf("(%.1f%%)\n", 100 * mean(diagnostics$is_potential_novel)))

cat("\nNovelty rate by true cell type:\n")
nov_tbl <- diagnostics |>
  group_by(true_label) |>
  summarise(n = n(),
            n_novel  = sum(is_potential_novel),
            pct_novel = round(100 * mean(is_potential_novel), 1),
            mean_error = round(mean(reconstruction_error), 4),
            mean_conf  = round(mean(confidence), 3),
            .groups = "drop") |>
  arrange(desc(pct_novel))
print(as.data.frame(nov_tbl))


# 7. Plots
message("\n[7] Generating plots...")

# Error distributions
p1 <- plot_error_distributions(
  ref_errors   = diag_out$ref_errors,
  query_errors = diagnostics$reconstruction_error,
  threshold    = quantile(diag_out$ref_errors, 0.95)
)
ggsave(file.path(OUTDIR, "01_error_distributions.png"),
       p1, width = 6, height = 4, dpi = 150)
message("Wrote: 01_error_distributions.png")

# Confidence distribution
p2 <- plot_confidence_distribution(diagnostics)
ggsave(file.path(OUTDIR, "02_confidence_distribution.png"),
       p2, width = 6, height = 4, dpi = 150)
message("Wrote: 02_confidence_distribution.png")


# 8. UMAP (requires Seurat)
message("\n[8] Building UMAP...")

rownames(query_matrix) <- make.unique(rownames(query_matrix))
seu <- CreateSeuratObject(counts = query_matrix)
seu$celltype <- query_celltypes

seu[["CellMentor"]] <- CreateDimReducObject(
  embeddings = t(as.matrix(h_project)),
  key        = "CellMentor_",
  assay      = DefaultAssay(seu),
  loadings   = as.matrix(best_model@W)
)
seu <- RunUMAP(seu, reduction = "CellMentor", dims = 1:K_VALUE, verbose = FALSE)

umap_coords <- Embeddings(seu, "umap")
diagnostics_ordered <- diagnostics[match(rownames(umap_coords), diagnostics$cell), ]

panels <- plot_umap_diagnostics(
  umap_coords = umap_coords,
  diagnostics = diagnostics_ordered,
  true_labels = diagnostics_ordered$true_label
)

# Use patchwork if available, otherwise save individually
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- (panels$true | panels$predicted) /
              (panels$confidence | panels$error) /
              (panels$novel | patchwork::plot_spacer())
  ggsave(file.path(OUTDIR, "03_umap_all_panels.png"),
         combined, width = 12, height = 13, dpi = 150)
  message("Wrote: 03_umap_all_panels.png")
} else {
  for (nm in names(panels)) {
    ggsave(file.path(OUTDIR, sprintf("03_umap_%s.png", nm)),
           panels[[nm]], width = 6, height = 5, dpi = 150)
  }
  message("Wrote: 03_umap_*.png (patchwork not available, saved individually)")
}


# -9. Save model + diagnostics for reuse 
message("\n[9] Saving results...")
saveRDS(
  list(best_model  = best_model,
       h_project   = h_project,
       diagnostics = diagnostics,
       ref_errors  = diag_out$ref_errors,
       seu         = seu),
  file = file.path(OUTDIR, "baron_muraro_results.rds")
)

sink(file.path(OUTDIR, "session_info.txt"))
sessionInfo()
sink()

message("\nAll done. Outputs in: ", OUTDIR)
