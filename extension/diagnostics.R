suppressPackageStartupMessages({
  library(Matrix)
})


# 1. Per-cell reconstruction error

#' Compute per-cell reconstruction error
#'
#' For each cell j, returns \eqn{||W h_j - x_j||_2^2}. Optionally normalizes by
#' \eqn{||x_j||_2^2} to produce a scale-invariant relative error, which is more
#' comparable across cells with different total expression.
#'
#' @param W gene x factor matrix (e.g. best_model@W)
#' @param H factor x cell matrix (e.g. h_project from project_data())
#' @param X gene x cell matrix of observed data, must match preprocessing used
#'   during training (i.e. best_model@matrices@data for the query, or the
#'   preprocessed reference matrix for reference cells)
#' @param normalize logical. If TRUE, returns relative error ||Wh - x||^2 / ||x||^2
#' @return numeric vector of length ncol(X)
#' @export
compute_reconstruction_error <- function(W, H, X, normalize = TRUE) {
  stopifnot(ncol(W) == nrow(H), nrow(W) == nrow(X), ncol(H) == ncol(X))

  W <- as.matrix(W)
  H <- as.matrix(H)
  X <- as.matrix(X)

  X_hat <- W %*% H
  residuals <- X - X_hat
  errors <- colSums(residuals^2)

  if (normalize) {
    norms <- colSums(X^2)
    # avoid divide-by-zero for cells that are all-zero after preprocessing
    errors <- errors / (norms + .Machine$double.eps)
  }

  names(errors) <- colnames(X)
  errors
}


# 2. Centroid-based confidence scores

#' Compute per-cell confidence scores in latent space
#'
#' Uses the reference H matrix and cell type labels to build a centroid for
#' each class, then for each query cell returns:
#'   - predicted_class: nearest centroid by cosine distance
#'   - nearest_distance: distance to that centroid
#'   - margin: (second-nearest distance) - (nearest distance); higher = more
#'     confident the cell belongs to one class rather than being between classes
#'   - confidence: softmax-style probability of the predicted class, computed
#'     from negated distances
#'
#' @param H_query factor x cell matrix for query data
#' @param H_ref factor x cell matrix for reference data
#' @param ref_celltypes character vector of length ncol(H_ref)
#' @param distance_metric "cosine" (default) or "euclidean"
#' @return data.frame with one row per query cell
#' @export
compute_confidence_scores <- function(H_query, H_ref, ref_celltypes,
                                      distance_metric = c("cosine", "euclidean")) {
  distance_metric <- match.arg(distance_metric)
  stopifnot(nrow(H_query) == nrow(H_ref),
            length(ref_celltypes) == ncol(H_ref))

  H_query <- as.matrix(H_query)
  H_ref <- as.matrix(H_ref)

  unique_types <- sort(unique(ref_celltypes))

  # Build class centroids: K x n_classes
  centroids <- sapply(unique_types, function(ct) {
    cells <- which(ref_celltypes == ct)
    if (length(cells) == 1) H_ref[, cells] else rowMeans(H_ref[, cells, drop = FALSE])
  })
  colnames(centroids) <- unique_types

  # Distance matrix: n_classes x n_query
  if (distance_metric == "cosine") {
    q_norms <- sqrt(colSums(H_query^2)) + .Machine$double.eps
    c_norms <- sqrt(colSums(centroids^2)) + .Machine$double.eps
    similarity <- t(centroids) %*% H_query
    similarity <- sweep(similarity, 2, q_norms, "/")
    similarity <- sweep(similarity, 1, c_norms, "/")
    dist_mat <- 1 - similarity
  } else {
    # euclidean
    dist_mat <- sapply(seq_len(ncol(H_query)), function(j) {
      sqrt(colSums((centroids - H_query[, j])^2))
    })
    rownames(dist_mat) <- unique_types
  }

  # Per-cell summaries
  predicted_idx <- apply(dist_mat, 2, which.min)
  predicted_class <- unique_types[predicted_idx]
  nearest_distance <- dist_mat[cbind(predicted_idx, seq_len(ncol(dist_mat)))]

  if (nrow(dist_mat) >= 2) {
    sorted_dist <- apply(dist_mat, 2, sort)
    margin <- sorted_dist[2, ] - sorted_dist[1, ]
  } else {
    margin <- rep(NA_real_, ncol(dist_mat))
  }

  # Softmax-style confidence over classes (lower distance -> higher weight)
  neg_d <- -dist_mat
  probs <- apply(neg_d, 2, function(x) {
    e <- exp(x - max(x))
    e / sum(e)
  })
  confidence <- apply(probs, 2, max)

  data.frame(
    cell = colnames(H_query),
    predicted_class = predicted_class,
    nearest_distance = nearest_distance,
    margin = margin,
    confidence = confidence,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# 3. Novel cell type flagging

#' Flag query cells as potentially novel based on reconstruction error
#'
#' Uses the reference's own reconstruction error distribution to calibrate a
#' threshold, then flags query cells whose reconstruction error exceeds that
#' threshold. In "global" mode a single threshold is used. In "per_class" mode
#' each query cell is compared against the error distribution of its predicted
#' class.
#'
#' @param query_errors numeric vector of reconstruction errors for query cells
#' @param ref_errors numeric vector of reconstruction errors for reference cells
#' @param ref_celltypes character vector, labels for reference cells
#' @param predicted_classes character vector of length(query_errors), the
#'   predicted class for each query cell (from compute_confidence_scores)
#' @param method "global" or "per_class"
#' @param threshold quantile of reference errors above which cells are flagged
#'   (default 0.95, i.e. cells worse than the worst 5% of reference cells)
#' @return data.frame with per-cell novelty flag and the threshold value used
#' @export
flag_novel_cells <- function(query_errors, ref_errors, ref_celltypes = NULL,
                             predicted_classes = NULL,
                             method = c("per_class", "global"),
                             threshold = 0.95) {
  method <- match.arg(method)
  stopifnot(threshold > 0, threshold < 1)

  is_novel <- logical(length(query_errors))
  thresh_used <- numeric(length(query_errors))

  if (method == "global") {
    t_val <- as.numeric(quantile(ref_errors, threshold, na.rm = TRUE))
    is_novel <- query_errors > t_val
    thresh_used[] <- t_val
  } else {
    stopifnot(!is.null(ref_celltypes), !is.null(predicted_classes),
              length(ref_celltypes) == length(ref_errors),
              length(predicted_classes) == length(query_errors))

    class_thresh <- tapply(ref_errors, ref_celltypes,
                           function(x) as.numeric(quantile(x, threshold, na.rm = TRUE)))

    for (i in seq_along(query_errors)) {
      pc <- predicted_classes[i]
      t_val <- if (!is.na(pc) && pc %in% names(class_thresh)) {
        class_thresh[[pc]]
      } else {
        as.numeric(quantile(ref_errors, threshold, na.rm = TRUE))
      }
      thresh_used[i] <- t_val
      is_novel[i] <- query_errors[i] > t_val
    }
  }

  data.frame(
    cell = names(query_errors),
    reconstruction_error = as.numeric(query_errors),
    novelty_threshold = thresh_used,
    is_potential_novel = is_novel,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# 4. Wrapper: run full projection diagnostics

#' Run projection diagnostics on a trained CellMentor model
#'
#' Given a fit CellMentor model and query data, this function calls the
#' existing CellMentor project_data() and adds:
#'   - per-cell reconstruction error (query AND reference)
#'   - centroid-based confidence scores (predicted class, margin, confidence)
#'   - novel cell type flags
#'
#' @param W gene x factor matrix from best_model@W
#' @param H_ref factor x cell reference embedding (best_model@H, or as
#'   returned from the training phase)
#' @param X_ref preprocessed gene x cell reference matrix used during training
#' @param ref_celltypes character vector of reference cell type labels
#' @param H_query factor x cell query embedding from project_data()
#' @param X_query preprocessed gene x cell query matrix
#'   (best_model@matrices@data)
#' @param novelty_threshold quantile for novelty flag (default 0.95)
#' @param novelty_method "per_class" or "global"
#' @param distance_metric "cosine" or "euclidean"
#' @return list with $diagnostics (per-cell data.frame), $ref_errors, and
#'   $class_centroids
#' @export
run_projection_diagnostics <- function(W, H_ref, X_ref, ref_celltypes,
                                       H_query, X_query,
                                       novelty_threshold = 0.95,
                                       novelty_method = c("per_class", "global"),
                                       distance_metric = c("cosine", "euclidean"),
                                       verbose = TRUE) {
  novelty_method <- match.arg(novelty_method)
  distance_metric <- match.arg(distance_metric)

  if (verbose) message("[diag] Computing reference reconstruction errors...")
  ref_errors <- compute_reconstruction_error(W, H_ref, X_ref, normalize = TRUE)

  if (verbose) message("[diag] Computing query reconstruction errors...")
  query_errors <- compute_reconstruction_error(W, H_query, X_query, normalize = TRUE)

  if (verbose) message("[diag] Computing confidence scores...")
  conf_df <- compute_confidence_scores(H_query, H_ref, ref_celltypes,
                                       distance_metric = distance_metric)

  if (verbose) message("[diag] Flagging potentially novel cells...")
  novel_df <- flag_novel_cells(
    query_errors = query_errors,
    ref_errors = ref_errors,
    ref_celltypes = ref_celltypes,
    predicted_classes = conf_df$predicted_class,
    method = novelty_method,
    threshold = novelty_threshold
  )

  diagnostics <- merge(conf_df, novel_df, by = "cell", sort = FALSE)

  if (verbose) {
    n_novel <- sum(diagnostics$is_potential_novel)
    message(sprintf("[diag] Flagged %d/%d cells as potentially novel (%.1f%%)",
                    n_novel, nrow(diagnostics), 100 * n_novel / nrow(diagnostics)))
  }

  list(
    diagnostics = diagnostics,
    ref_errors = ref_errors,
    parameters = list(
      novelty_threshold = novelty_threshold,
      novelty_method = novelty_method,
      distance_metric = distance_metric
    )
  )
}
