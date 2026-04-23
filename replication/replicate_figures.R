set.seed(100)

R_LIBS  <- "/projectnb/ds596/students/aliviap/R_libs_4.5.2"
RESULTS <- "/projectnb/ds596/students/aliviap/extension/figures/baron_muraro_results.rds"
OUTDIR  <- "/projectnb/ds596/students/aliviap/extension/figures/replication"

.libPaths(c(R_LIBS, .libPaths()))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(CellMentor)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(aricode)   # ARI + NMI
  library(patchwork)
})

stopifnot(file.exists(RESULTS))
message("Loading saved results from: ", RESULTS)
res <- readRDS(RESULTS)

best_model   <- res$best_model
h_project    <- res$h_project
diagnostics  <- res$diagnostics
seu          <- res$seu

# True Muraro labels aligned to the Seurat cell order
true_labels <- seu$celltype
K_VALUE     <- best_model@parameters$rank


# Figure A: Publication-style UMAP of query projection
# Mirrors Figure A (page 7) of Hevdeli et al. 2025: Muraro query projected
# onto the Baron reference latent space, colored by true cell type.

message("Building Figure A: publication-style projection UMAP...")

umap_df <- as.data.frame(Embeddings(seu, "umap"))
umap_df$celltype <- true_labels

# Use the same 8 Muraro cell types
muraro_types  <- sort(unique(umap_df$celltype))
palette_muraro <- c(
  acinar      = "#E45756",
  alpha       = "#F28E2B",
  beta        = "#76B7B2",
  delta       = "#59A14F",
  ductal      = "#EDC948",
  endothelial = "#B07AA1",
  epsilon     = "#FF9DA7",
  gamma       = "#9C755F"
)
# Fall back for any unexpected labels
fill_colors <- palette_muraro[muraro_types]
fill_colors[is.na(fill_colors)] <- "grey60"
names(fill_colors) <- muraro_types

figA <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_manual(values = fill_colors, name = "Cell type") +
  labs(
    title    = "CellMentor: Muraro query projected onto Baron reference",
    subtitle = sprintf("Baron (reference, n=8,569) → Muraro (query, n=2,042)  |  K=%d", K_VALUE),
    x = "UMAP 1", y = "UMAP 2",
    caption  = "Reproduces Fig. A (p.7), Hevdeli et al. Nat. Commun. 2025"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.caption     = element_text(size = 8, color = "grey50")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(OUTDIR, "figA_umap_projection.png"),
       figA, width = 7, height = 5, dpi = 200)
message("Wrote: figA_umap_projection.png")


# Figure B: ARI / NMI vs clustering resolution
# Sweeps Seurat clustering resolution and reports ARI + NMI against true
# Muraro cell type labels at each resolution. Mirrors the resolution sweep
# approach used in the paper's benchmarking.

message("Building Figure B: ARI/NMI resolution sweep...")

resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

seu <- FindNeighbors(seu, reduction = "CellMentor",
                     dims = 1:K_VALUE, verbose = FALSE)

sweep_results <- lapply(resolutions, function(r) {
  seu_r <- FindClusters(seu, resolution = r, verbose = FALSE)
  cluster_labels <- as.character(Idents(seu_r))

  ari_val <- ARI(cluster_labels, as.character(true_labels))
  nmi_val <- NMI(cluster_labels, as.character(true_labels))

  data.frame(resolution = r, ARI = ari_val, NMI = nmi_val)
})

sweep_df <- do.call(rbind, sweep_results)
print(sweep_df)
write.csv(sweep_df, file.path(OUTDIR, "figB_resolution_sweep.csv"), row.names = FALSE)

# Find paper default (resolution = 0.1)
default_res <- 0.1

figB <- sweep_df |>
  tidyr::pivot_longer(cols = c(ARI, NMI), names_to = "Metric", values_to = "Score") |>
  ggplot(aes(x = resolution, y = Score, color = Metric, group = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = default_res, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = default_res, y = 0.05,
           label = "paper default", hjust = -0.08, size = 3.2, color = "grey40") +
  scale_color_manual(values = c(ARI = "#4C72B0", NMI = "#DD8452")) +
  scale_x_continuous(breaks = resolutions) +
  labs(
    title    = "CellMentor clustering performance vs. resolution",
    subtitle = "Baron → Muraro projection  |  evaluated against true Muraro labels",
    x = "Clustering resolution", y = "Score", color = "Metric"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "figB_ari_nmi_resolution.png"),
       figB, width = 6, height = 4.5, dpi = 200)
message("Wrote: figB_ari_nmi_resolution.png")


# Figure C: ARI / NMI at paper default (numeric comparison
# Reports the single-number ARI + NMI at resolution = 0.1 for direct
# comparison with the values reported in the CellMentor paper (Table 1 / S1).

message("Computing ARI/NMI at default resolution (0.1) for paper comparison...")

seu_default <- FindClusters(seu, resolution = 0.1, verbose = FALSE)
ari_default <- ARI(as.character(Idents(seu_default)), as.character(true_labels))
nmi_default <- NMI(as.character(Idents(seu_default)), as.character(true_labels))

metrics_df <- data.frame(
  dataset       = "Baron (ref) → Muraro (query)",
  method        = "CellMentor",
  resolution    = 0.1,
  ARI           = round(ari_default, 3),
  NMI           = round(nmi_default, 3),
  K             = K_VALUE,
  alpha         = best_model@parameters$alpha,
  beta          = best_model@parameters$beta
)
print(metrics_df)
write.csv(metrics_df, file.path(OUTDIR, "figC_ari_nmi_table.csv"), row.names = FALSE)
message("Wrote: figC_ari_nmi_table.csv")

cat(sprintf("\nARI at resolution 0.1: %.3f\nNMI at resolution 0.1: %.3f\n",
            ari_default, nmi_default))
cat("Compare these values with Table 1 / Supplementary Table S1 of Hevdeli et al. 2025\n")

message("\nAll replication figures written to: ", OUTDIR)
