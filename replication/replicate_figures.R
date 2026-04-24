# =============================================================================
# replicate_figures.R
# -----------------------------------------------------------------------------
# Reproduces publication-adjacent figures using the Baron -> Muraro dataset,
# the only non-deprecated dataset in the current CellMentor package.
#
# Requires: baron_muraro_results.rds produced by run_extension_demo.R
#
# Produces (in extension/figures/replication/):
#   figA_umap_projection.png     -- mirrors Figure A (p.7) of the CellMentor paper
#   figB_ari_nmi_resolution.png  -- ARI/NMI vs clustering resolution sweep
#   figB_resolution_sweep.csv    -- numeric scores at each resolution
#   figC_ari_nmi_table.csv       -- ARI/NMI at default resolution (compare Table 1)
#
# Usage (SCC):
#   qsub replication/submit_replication.sh
# Or interactively:
#   Rscript replication/replicate_figures.R
# =============================================================================

set.seed(100)

BASE    <- "/projectnb/ds596/students/aliviap/DS596-Final-Project"
R_LIBS  <- "/projectnb/ds596/students/aliviap/R_libs_4.5.2"
RESULTS <- file.path(BASE, "extension", "figures", "baron_muraro_results.rds")
OUTDIR  <- file.path(BASE, "extension", "figures", "replication")

.libPaths(c(R_LIBS, .libPaths()))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Validate results file before doing anything
if (!file.exists(RESULTS)) {
  stop("Results file not found at: ", RESULTS,
       "\nRun extension/run_extension_demo.R first.")
}
message("Loading results from: ", RESULTS)

# Install tidyr if missing (needed for pivot_longer)
if (!requireNamespace("tidyr", quietly = TRUE)) {
  message("Installing tidyr...")
  install.packages("tidyr", lib = R_LIBS, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(CellMentor)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(aricode)
  library(patchwork)
})

res <- readRDS(RESULTS)
best_model  <- res$best_model
h_project   <- res$h_project
seu         <- res$seu
true_labels <- seu$celltype
K_VALUE     <- best_model@parameters$rank

message(sprintf("Loaded: %d query cells, K=%d", ncol(h_project), K_VALUE))


# ---- Figure A: Publication-style UMAP ---------------------------------------

message("Building Figure A: publication-style projection UMAP...")

umap_df <- as.data.frame(Embeddings(seu, "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype <- true_labels

muraro_types <- sort(unique(umap_df$celltype))
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
fill_colors <- palette_muraro[muraro_types]
fill_colors[is.na(fill_colors)] <- "grey60"
names(fill_colors) <- muraro_types

figA <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 1.0, alpha = 0.85) +
  scale_color_manual(values = fill_colors, name = "Cell type") +
  labs(
    title    = "CellMentor: Muraro query projected onto Baron reference",
    subtitle = sprintf("Baron (ref, n=8,569) → Muraro (query, n=2,042)  |  K=%d", K_VALUE),
    x = "UMAP 1", y = "UMAP 2",
    caption = "Reproduces Fig. A (p.7), Hevdeli et al. Nat. Commun. 2025"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(size = 8, color = "grey50")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(OUTDIR, "figA_umap_projection.png"),
       figA, width = 7, height = 5, dpi = 200)
message("Wrote: figA_umap_projection.png")


# ---- Figure B: ARI/NMI resolution sweep ------------------------------------

message("Building Figure B: ARI/NMI resolution sweep...")

seu <- FindNeighbors(seu, reduction = "CellMentor",
                     dims = 1:K_VALUE, verbose = FALSE)

resolutions  <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
sweep_results <- lapply(resolutions, function(r) {
  seu_r <- FindClusters(seu, resolution = r, verbose = FALSE)
  data.frame(
    resolution = r,
    ARI = ARI(as.character(Idents(seu_r)), as.character(true_labels)),
    NMI = NMI(as.character(Idents(seu_r)), as.character(true_labels))
  )
})
sweep_df <- do.call(rbind, sweep_results)
print(sweep_df)
write.csv(sweep_df, file.path(OUTDIR, "figB_resolution_sweep.csv"), row.names = FALSE)

figB <- sweep_df |>
  pivot_longer(cols = c(ARI, NMI), names_to = "Metric", values_to = "Score") |>
  ggplot(aes(x = resolution, y = Score, color = Metric, group = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.1, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = 0.1, y = min(sweep_df$ARI, sweep_df$NMI) - 0.03,
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


# ---- Figure C: ARI/NMI at paper default -------------------------------------

message("Computing ARI/NMI at default resolution (0.1)...")

seu_default <- FindClusters(seu, resolution = 0.1, verbose = FALSE)
ari_val <- ARI(as.character(Idents(seu_default)), as.character(true_labels))
nmi_val <- NMI(as.character(Idents(seu_default)), as.character(true_labels))

metrics_df <- data.frame(
  dataset    = "Baron (ref) -> Muraro (query)",
  method     = "CellMentor",
  resolution = 0.1,
  ARI        = round(ari_val, 3),
  NMI        = round(nmi_val, 3),
  K          = K_VALUE,
  alpha      = best_model@parameters$alpha,
  beta       = best_model@parameters$beta
)
print(metrics_df)
write.csv(metrics_df, file.path(OUTDIR, "figC_ari_nmi_table.csv"), row.names = FALSE)
message("Wrote: figC_ari_nmi_table.csv")

cat(sprintf("\nARI at resolution 0.1: %.3f\nNMI at resolution 0.1: %.3f\n",
            ari_val, nmi_val))
cat("Compare with Table 1 / Supplementary Table S1 of Hevdeli et al. 2025\n")

message("\nAll replication figures written to: ", OUTDIR)
