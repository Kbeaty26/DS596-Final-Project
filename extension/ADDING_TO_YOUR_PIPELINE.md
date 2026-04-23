# Adding the Projection Diagnostics Extension to Your CellMentor Pipeline

**You do not need to re-run CellMentor.** The extension works on the outputs
your pipeline has already produced.

---

Make sure your saved results include these objects. If you saved them in an
`.rds` file, load it and check:

```r
res <- readRDS("your_results.rds")

best_model    # the fitted CSFNMF object from CellMentor()
h_project     # the K x n_query matrix from project_data()
```

---

## Step 1 - Source the extension functions

Add this near the top of your script, after loading CellMentor:

```r
.libPaths(c("/projectnb/ds596/students/aliviap/R_libs_4.5.2", .libPaths()))

source("/projectnb/ds596/students/aliviap/extension/diagnostics.R")
source("/projectnb/ds596/students/aliviap/extension/plotting.R")
```

---

## Step 2 - Extract the pieces the extension needs

Add this block immediately after your `project_data()` call:

```r
# Extract reference matrix and embedding from the trained model
X_ref  <- best_model@matrices@ref     # preprocessed reference expression matrix
H_ref  <- best_model@H                # reference cell embeddings (K x n_ref)
X_query <- best_model@matrices@data   # preprocessed query expression matrix

# Reference cell type labels
# Try the model's annotation slot first; fall back to your original labels
ref_ct <- tryCatch(
  best_model@annotation$celltype,
  error = function(e) your_reference_celltypes[colnames(X_ref)]
  #                   ^^^^^^^^^^^^^^^^^^^^^^^^
  #                   replace with whatever you named your reference labels
)
```

---

## Step 3 - Run the diagnostics

```r
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
```

`diagnostics` is now a data.frame with one row per query cell containing:

| Column | Description |
|---|---|
| `cell` | cell barcode |
| `predicted_class` | nearest reference centroid (by cosine distance) |
| `nearest_distance` | cosine distance to that centroid |
| `margin` | gap between nearest and second-nearest centroid |
| `confidence` | softmax probability of predicted class |
| `reconstruction_error` | normalized `‖Wh - x‖² / ‖x‖²` |
| `novelty_threshold` | per-class threshold used |
| `is_potential_novel` | `TRUE` if error exceeds threshold |

---

## Step 4 - Attach your true labels (if known) and save

```r
# If you know the true query cell types (for evaluation):
diagnostics$true_label <- your_query_celltypes[diagnostics$cell]

write.csv(diagnostics, "my_dataset_diagnostics.csv", row.names = FALSE)
```

---

## Step 5 - Generate the diagnostic plots

```r
library(ggplot2)

# Error distribution: reference vs. query, with novelty threshold marked
p1 <- plot_error_distributions(
  ref_errors   = diag_out$ref_errors,
  query_errors = diagnostics$reconstruction_error,
  threshold    = quantile(diag_out$ref_errors, 0.95)
)
ggsave("error_distributions.png", p1, width = 6, height = 4, dpi = 150)

# Confidence histogram colored by novelty flag
p2 <- plot_confidence_distribution(diagnostics)
ggsave("confidence_distribution.png", p2, width = 6, height = 4, dpi = 150)
```

If you have a UMAP already built (Seurat object with a UMAP reduction):

```r
umap_coords <- Seurat::Embeddings(seu, "umap")
diagnostics_ordered <- diagnostics[match(rownames(umap_coords), diagnostics$cell), ]

panels <- plot_umap_diagnostics(
  umap_coords = umap_coords,
  diagnostics = diagnostics_ordered,
  true_labels = diagnostics_ordered$true_label   # or NULL if unknown
)
# panels is a list: $true, $predicted, $confidence, $error, $novel
# combine with patchwork or save individually:
ggsave("umap_novelty.png", panels$novel, width = 6, height = 5, dpi = 150)
```

---

## Quick summary of novelty rates

Print a per-cell-type novelty summary:

```r
library(dplyr)
diagnostics |>
  group_by(true_label) |>
  summarise(
    n         = n(),
    n_novel   = sum(is_potential_novel),
    pct_novel = round(100 * mean(is_potential_novel), 1),
    mean_error = round(mean(reconstruction_error), 4),
    .groups = "drop"
  ) |>
  arrange(desc(pct_novel)) |>
  print()
```

---

## Troubleshooting

**"Slot not found" error on `@matrices@ref` or `@annotation`**
Run `slotNames(best_model@matrices)` and `slotNames(best_model)` in R to see
the actual slot names in your installed version of CellMentor, then edit the
extraction block in Step 2 accordingly.

**Extension functions not found**
Make sure you ran the `source()` calls in Step 1 and that the paths to
`diagnostics.R` and `plotting.R` are correct.

**Very high novelty rates (>50% flagged)**
This usually means the query and reference are from very different tissues or
species. Try switching to `novelty_method = "global"` which uses a single
population-wide threshold rather than per-class thresholds.
