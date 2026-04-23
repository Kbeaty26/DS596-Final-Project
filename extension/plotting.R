suppressPackageStartupMessages({
  library(ggplot2)
})


#' Plot distribution of reference vs. query reconstruction errors
#'
#' @param ref_errors numeric vector of reference cell errors
#' @param query_errors numeric vector of query cell errors
#' @param threshold optional numeric threshold to draw as a vertical line
#' @return ggplot
#' @export
plot_error_distributions <- function(ref_errors, query_errors, threshold = NULL) {
  df <- rbind(
    data.frame(error = ref_errors, dataset = "Reference"),
    data.frame(error = query_errors, dataset = "Query")
  )
  p <- ggplot(df, aes(x = error, fill = dataset)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c(Reference = "#4C72B0", Query = "#DD8452")) +
    labs(x = "Normalized reconstruction error", y = "Density",
         fill = NULL,
         title = "Per-cell reconstruction error") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  if (!is.null(threshold)) {
    p <- p + geom_vline(xintercept = threshold, linetype = "dashed",
                        color = "firebrick", linewidth = 0.6) +
      annotate("text", x = threshold, y = 0, label = " threshold ",
               vjust = -0.5, hjust = 0, color = "firebrick", size = 3.5)
  }
  p
}


#' Plot confidence score distribution, colored by whether the cell was flagged
#' as potentially novel
#'
#' @param diagnostics data.frame from run_projection_diagnostics()
#' @return ggplot
#' @export
plot_confidence_distribution <- function(diagnostics) {
  ggplot(diagnostics, aes(x = confidence, fill = is_potential_novel)) +
    geom_histogram(bins = 40, position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("FALSE" = "#4C72B0", "TRUE" = "#C44E52"),
                      labels = c("FALSE" = "In reference", "TRUE" = "Flagged novel"),
                      name = NULL) +
    labs(x = "Confidence score", y = "Query cells",
         title = "Confidence distribution with novelty flag") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}


#' Overlay diagnostics on a UMAP embedding
#'
#' @param umap_coords 2-column matrix or data.frame (cell x {UMAP_1, UMAP_2})
#' @param diagnostics data.frame from run_projection_diagnostics()
#' @param true_labels optional, true cell type labels (for evaluation)
#' @return list of ggplot objects
#' @export
plot_umap_diagnostics <- function(umap_coords, diagnostics, true_labels = NULL) {
  df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    predicted = diagnostics$predicted_class,
    confidence = diagnostics$confidence,
    error = diagnostics$reconstruction_error,
    is_novel = diagnostics$is_potential_novel
  )
  if (!is.null(true_labels)) df$true_label <- true_labels

  base <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

  p_predicted <- ggplot(df, aes(UMAP_1, UMAP_2, color = predicted)) +
    geom_point(size = 0.5, alpha = 0.8) +
    labs(title = "Predicted cell type", color = NULL) +
    base +
    guides(color = guide_legend(override.aes = list(size = 2.5)))

  p_confidence <- ggplot(df, aes(UMAP_1, UMAP_2, color = confidence)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_viridis_c(option = "C") +
    labs(title = "Confidence score") +
    base

  p_error <- ggplot(df, aes(UMAP_1, UMAP_2, color = error)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_viridis_c(option = "B") +
    labs(title = "Reconstruction error") +
    base

  p_novel <- ggplot(df, aes(UMAP_1, UMAP_2, color = is_novel)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#C44E52"),
                       labels = c("FALSE" = "In reference", "TRUE" = "Novel?"),
                       name = NULL) +
    labs(title = "Novelty flag") +
    base +
    guides(color = guide_legend(override.aes = list(size = 2.5)))

  out <- list(predicted = p_predicted,
              confidence = p_confidence,
              error = p_error,
              novel = p_novel)

  if (!is.null(true_labels)) {
    p_true <- ggplot(df, aes(UMAP_1, UMAP_2, color = true_label)) +
      geom_point(size = 0.5, alpha = 0.8) +
      labs(title = "True cell type", color = NULL) +
      base +
      guides(color = guide_legend(override.aes = list(size = 2.5)))
    out$true <- p_true
  }

  out
}
