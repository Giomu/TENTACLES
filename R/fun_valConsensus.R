#' @title Select Top Gene Combinations Based on Clustering Metrics
#'
#' @description
#' Ranks and selects the top gene combinations based on the average value of a specified metric
#' across multiple clustering methods. Returns a data frame with the top N combinations and all evaluated metrics.
#'
#' @param results A named list containing the clustering results for each gene combination. Each element should be
#' a list of clustering results for the given gene combination, as produced by \code{valConsensus}.
#' @param N An integer specifying the number of top combinations to return.
#' @param metric A character string indicating the evaluation metric to use for ranking.
#' Possible values are \code{"Accuracy"}, \code{"Precision"}, \code{"Recall"}, and \code{"FScore"}.
#'
#' @return A data frame containing the top N gene combinations, the average metric across all clustering methods
#' (\code{Mean_Metric}), and the individual metrics for each method.
#'
#' @details
#' This function aggregates the results from multiple clustering methods for all tested gene combinations.
#' For each combination, it computes the mean of the specified metric across all methods and selects the top N combinations.
#' The resulting data frame includes the gene combination identifier, the average metric, and the individual
#' metrics for each method and metric type.
#'
#' @seealso \code{\link{valConsensus}}
#'
#' @examples
#' \dontrun{
#' # Example structure of 'results' (as produced by valConsensus)
#' results <- list(
#'   GeneA_GeneB = list(
#'     KMeans = list(Accuracy = 0.85, Precision = 0.82, Recall = 0.80, FScore = 0.81),
#'     GMM    = list(Accuracy = 0.86, Precision = 0.83, Recall = 0.82, FScore = 0.82)
#'   ),
#'   GeneC_GeneD = list(
#'     KMeans = list(Accuracy = 0.78, Precision = 0.75, Recall = 0.74, FScore = 0.74),
#'     GMM    = list(Accuracy = 0.80, Precision = 0.76, Recall = 0.76, FScore = 0.75)
#'   )
#' )
#' top_df <- selectTopCombinations(results, N = 1, metric = "Accuracy")
#' }
#' @noRd
selectTopCombinations <- function(results, N, metric) {
  allowed_metrics <- c("Accuracy", "Precision", "Recall", "FScore")
  if (!metric %in% allowed_metrics) {
    cli::cli_abort("Invalid metric selected. Allowed metrics are: {paste(allowed_metrics, collapse = ', ')}.")
  }

  methods <- names(results[[1]])  # assume all combinations have same methods
  metric_names <- paste0(rep(methods, each = length(allowed_metrics)), "_", allowed_metrics)

  results_list <- list()
  for (combo_name in names(results)) {
    combo_results <- results[[combo_name]]

    # Extract metric values from all clustering methods for this gene combination
    metric_values <- unlist(lapply(combo_results, `[[`, metric))
    mean_metric <- mean(metric_values)

    # Extract all metrics for all methods
    row_data <- unlist(lapply(methods, function(m) {
      unlist(combo_results[[m]][allowed_metrics])
    }))

    # Store as row
    results_list[[combo_name]] <- c(combo_name, mean_metric, row_data)
  }

  # Convert list to data frame
  results_df <- as.data.frame(do.call(rbind, results_list), stringsAsFactors = FALSE)
  colnames(results_df) <- c("Gene_Combination", "Mean_Metric", metric_names)

  # Sort and select top N
  results_df$Mean_Metric <- as.numeric(results_df$Mean_Metric)
  results_df[, -1] <- lapply(results_df[, -1], as.numeric)

  # Sort and select top N
  results_df <- results_df[order(-results_df$Mean_Metric), , drop = FALSE]
  results_df <- results_df[seq_len(min(N, nrow(results_df))), , drop = FALSE]

  return(results_df)
}

#' @title Evaluate Clustering Results with Metric Orientation
#'
#' @description
#' Evaluates clustering results against true class labels using standard metrics.
#' Automatically checks both cluster label orientations (i.e., 0/1 or 1/0 assignment) and
#' selects the orientation with the best metric value.
#'
#' @param pred A vector of predicted cluster assignments (typically 0 and 1).
#' @param truth A vector of true class labels (must be coercible to factors with two levels).
#' @param metric A character string specifying the primary evaluation metric to use for orientation selection.
#' Must be one of \code{"Accuracy"}, \code{"Precision"}, \code{"Recall"}, or \code{"FScore"}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{clusters}{A vector of predicted cluster assignments (after optimal orientation, as integers).}
#'   \item{Accuracy}{The accuracy of the clustering result.}
#'   \item{Precision}{The precision of the clustering result.}
#'   \item{Recall}{The recall of the clustering result.}
#'   \item{FScore}{The F1 score of the clustering result.}
#' }
#'
#' @details
#' Since unsupervised clustering labels (e.g., 0/1) may be assigned arbitrarily compared to
#' ground truth, this function evaluates both the predicted and inverted cluster assignments.
#' It then returns the set of metrics corresponding to the orientation that yields the
#' highest value for the selected metric.
#'
#' @seealso \code{\link[yardstick]{accuracy}}, \code{\link[yardstick]{precision}},
#'   \code{\link[yardstick]{recall}}, \code{\link[yardstick]{f_meas}}
#'
#' @examples
#' \dontrun{
#' true_labels <- factor(c(1, 0, 1, 0, 1))
#' predicted_clusters <- c(1, 0, 1, 0, 0)
#' evaluate_one_side(predicted_clusters, true_labels, metric = "Accuracy")
#' }
#' @noRd
# Helper function for assessment based on a choose metric
evaluate_one_side <- function(pred, truth, metric = metric) {
  df <- data.frame(
    labels = factor(truth),
    cluster = factor(pred),
    cluster_inv = factor(1 - as.numeric(as.character(pred)))
  )

  metric_fun <- switch(metric,
                       Accuracy  = yardstick::accuracy,
                       Precision = yardstick::precision,
                       Recall    = yardstick::recall,
                       FScore    = yardstick::f_meas,
                       cli::cli_abort("Unsupported metric: {metric}")
  )

  val_cluster     <- metric_fun(df, truth = labels, estimate = cluster, event_level = "second")[[".estimate"]]
  val_cluster_inv <- metric_fun(df, truth = labels, estimate = cluster_inv, event_level = "second")[[".estimate"]]

  best_col <- if (val_cluster >= val_cluster_inv) "cluster" else "cluster_inv"
  df_pred <- data.frame(labels = df$labels, pred = df[[best_col]])

  return(list(
    clusters  = as.integer(as.character(df_pred$pred)),
    Accuracy  = yardstick::accuracy(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
    Precision = yardstick::precision(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
    Recall    = yardstick::recall(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
    FScore    = yardstick::f_meas(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]]
  ))
}

# Helper function to build the binary gene x combination matrix of the top
# combinations. combos is a list of gene vectors in rank order; columns are
# named C1, C2, ... and genes are sorted by the number of combinations
# containing them, then by name
gene_combination_matrix <- function(combos) {
  genes <- unique(unlist(combos, use.names = FALSE))
  presence <- vapply(combos, function(combo) as.integer(genes %in% combo), integer(length(genes)))
  gene_matrix <- matrix(
    presence,
    nrow = length(genes),
    dimnames = list(genes, paste0("C", seq_along(combos)))
  )
  gene_matrix[order(-rowSums(gene_matrix), genes), , drop = FALSE]
}

# Helper function to count, for each gene of a gene x combination matrix, the
# combinations containing it
gene_frequency <- function(gene_matrix) {
  data.frame(
    Gene = rownames(gene_matrix),
    Count = as.integer(rowSums(gene_matrix)),
    Combinations = apply(gene_matrix, 1, function(x) {
      paste(colnames(gene_matrix)[x == 1], collapse = ", ")
    }),
    row.names = NULL
  )
}

#' @title valConsensus
#' @description This function validates the consensus genes in a provided dataset using clustering methods.
#'
#' @param df.count A data frame with the gene expression data, with samples in rows and genes in
#' columns. A matrix is not supported.
#' @param gene.list A character vector with the consensus genes to evaluate, at least two. All of them
#' must be columns of \code{df.count}.
#' @param class A vector with the true class of each sample (row of \code{df.count}), coded as 0 (control)
#' and 1 (case). It can be numeric or a factor with levels \code{"0"} and \code{"1"}. Class 1 is the positive
#' class for Precision, Recall and FScore.
#' @param N An integer specifying the number of top gene combinations to return. If it is larger than
#' the number of possible combinations, it is reduced to that number with a message.
#' @param metric A character string specifying the metric used to rank the combinations and to orient
#' the clusters (see Details). Possible values are \code{"Accuracy"}, \code{"Precision"}, \code{"Recall"},
#' and \code{"FScore"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{allResults}{A list with the clustering results for each combination of genes and clustering
#'   methods.}
#'   \item{topCombinations}{A data frame with the top \code{N} combinations, one per row, sorted by
#'   decreasing \code{Mean_Metric}. \code{Gene_Combination} holds the genes of the combination joined
#'   by \code{"_"}, \code{Mean_Metric} the mean of \code{metric} across the six methods, and the
#'   remaining columns, named \code{<Method>_<Metric>} (e.g. \code{KMeans_FScore}), the four metrics
#'   of each method.}
#'   \item{geneFrequency}{A data frame with, for each gene in the top combinations, the number of
#'   combinations containing it (\code{Count}) and their identifiers (\code{Combinations}). Combination
#'   \code{Ci} is the i-th row of \code{topCombinations}.}
#'   \item{geneCombinations}{A binary matrix with genes in rows and top combinations (\code{C1},
#'   \code{C2}, ...) in columns, with 1 when the gene is in the combination. Genes are ordered as in
#'   \code{geneFrequency}.}
#' }
#'
#' @details
#' For every combination of 2 to \code{length(gene.list)} genes (\eqn{2^k - k - 1} combinations
#' for \eqn{k} genes), the expression of the genes is scaled (centered and divided by the
#' standard deviation) and the samples are split into two clusters by six methods:
#' \describe{
#'   \item{KMeans}{k-means (MacQueen algorithm) on the scaled data.}
#'   \item{GMM}{Gaussian mixture model with two components (\code{mclust::Mclust()}).}
#'   \item{HC}{Hierarchical clustering with complete linkage on Euclidean distances, cut into two clusters.}
#'   \item{PCA}{k-means on the first two principal components.}
#'   \item{tSNE}{k-means on a two-dimensional t-SNE embedding (perplexity 10).}
#'   \item{UMAP}{k-means on a two-dimensional UMAP embedding (default settings).}
#' }
#'
#' Cluster labels are arbitrary, so each clustering is evaluated with both assignments of
#' its clusters to classes 0 and 1, and the one with the higher \code{metric} is kept. Accuracy,
#' Precision, Recall and FScore are then computed for that assignment.
#'
#' Combinations are ranked by the mean of \code{metric} across the six methods, and the top
#' \code{N} are returned. For these, the function counts how many combinations contain each gene
#' and prints two plots: the metrics heatmap (\code{\link{plotTopMetrics}}) and the gene presence
#' heatmap (\code{\link{plotGeneCombinations}}).
#'
#' The random seed is fixed at 123 inside the function, so results are reproducible but
#' cannot be varied. t-SNE with perplexity 10 requires at least 31 samples. The number of
#' combinations grows exponentially with the number of genes: 15 genes give 32,752
#' combinations.
#'
#' @importFrom mclust Mclust mclustBIC
#' @importFrom yardstick accuracy precision recall f_meas
#'
#' @examples
#' \dontrun{
#' data("acc.count", "acc.clin", package = "TENTACLES")
#' pp <- preProcess(acc.count[, 1:100], acc.clin, is.normalized = FALSE,
#'                  batch = "patient.gender", plot = FALSE)
#' data <- pp@processed$adjusted.data
#' genes <- c("A1BG", "A1CF", "A2BP1", "A2LD1")
#' vc <- valConsensus(
#'   df.count = data[, genes], gene.list = genes,
#'   class = data$class, N = 5, metric = "FScore"
#' )
#' vc$geneFrequency
#' }
#'
#' @export
valConsensus <- function(df.count, gene.list, class, N = 10, metric = "FScore") {

  cli::cli_h1("Validating Consensus Genes")

  # Fail before testing any combination, naming the missing features
  missing_genes <- setdiff(gene.list, colnames(df.count))
  if (length(missing_genes) > 0) {
    cli::cli_abort(c(
      "{length(missing_genes)} feature{?s} of {.arg gene.list} {?is/are} not in {.arg df.count}:",
      "x" = "{.val {missing_genes}}"
    ))
  }

  gene_list <- gene.list
  labels <- class

  # Maximum number of genes to consider
  max_genes <- length(unique(gene_list))

  # List to save results for all combinations
  all_results <- list()
  # Genes of each combination, keyed like all_results. Kept as vectors because
  # gene names may contain the "_" used to build the keys
  combo_genes <- list()

  # Compute total number of combinations
  total_combinations <- as.numeric(sum(sapply(2:max_genes, function(n) choose(length(gene_list), n))))
  cli::cli_alert_info("Applying models to {total_combinations} combination{?s} ...")
  if (N > total_combinations) {
    cli::cli_alert_warning(
      "{.arg N} = {N} is larger than the {total_combinations} possible combination{?s}. Using N = {total_combinations}."
    )
    N <- total_combinations
  }
  if (total_combinations > 1000000) {
    cli::cli_alert_warning("The number of combinations is > 1M. This may take a long time to complete. Consider reducing the number of genes.")
  }

  # Initialize progress bar
  progress_bar <- cli::cli_progress_bar(
    format = "Tested combinations: {cli::pb_bar} {cli::pb_percent} ({cli::pb_current}/{cli::pb_total})",
    total = total_combinations,
    clear = FALSE
  )

  withr::with_seed(
    seed = 123,
    code = {
      # For each number of genes up to the maximum
      for (n in 2:max_genes) {
        gene_combinations <- combn(gene_list, n, simplify = FALSE)

        # For each combination of genes
        for (combo in gene_combinations) {
          # Increment progress bar
          cli::cli_progress_update(inc = 1)
          # Subset the table with the specified genes
          count_subset <- df.count[, combo, drop = FALSE]

          # Ensure data is numeric and scaled
          count_subset <- as.data.frame(lapply(count_subset, as.numeric))
          count_subset <- scale(count_subset)
          rownames(count_subset) <- rownames(df.count)

          clustering_results <- list()

          # K-Means Clustering
          kmeans_result <- stats::kmeans(count_subset, centers = 2, iter.max = 100, algorithm = "MacQueen")
          clustering_results$KMeans <- evaluate_one_side(kmeans_result$cluster - 1, labels, metric)

          # Gaussian Mixture Model (GMM)
          gmm_result <- mclust::Mclust(count_subset, G = 2, verbose = FALSE)
          clustering_results$GMM <- evaluate_one_side(gmm_result$classification - 1, labels, metric)

          # Hierarchical Clustering
          hc_result <- stats::hclust(dist(count_subset), method = "complete")
          hc_clusters <- stats::cutree(hc_result, k = 2) - 1
          clustering_results$HC <- evaluate_one_side(hc_clusters, labels, metric)

          # PCA + KMeans
          pca_result <- stats::prcomp(count_subset, scale = FALSE)
          pca_clusters <- stats::kmeans(pca_result$x[, 1:2], centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1
          clustering_results$PCA <- evaluate_one_side(pca_clusters, labels, metric)

          # t-SNE + KMeans
          tsne_result <- Rtsne::Rtsne(count_subset, dims = 2, perplexity = 10, check_duplicates = FALSE)
          tsne_clusters <- stats::kmeans(tsne_result$Y, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1
          clustering_results$tSNE <- evaluate_one_side(tsne_clusters, labels, metric)

          # UMAP + KMeans
          umap_result <- umap::umap(count_subset)
          umap_clusters <- stats::kmeans(umap_result$layout, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1
          clustering_results$UMAP <- evaluate_one_side(umap_clusters, labels, metric)

          # Save results for this combination
          all_results[[paste(combo, collapse = "_")]] <- clustering_results
          combo_genes[[paste(combo, collapse = "_")]] <- combo
        }
      }
    }
  )

  # End Progress bar
  cli::cli_progress_done()
  cli::cli_alert_success("All combinations tested successfully!")

  # Add a new element to the list to store the top genes combinations by running the selectTopCombinations function
  cli::cli_h3("Top Combinations Selection")
  cli::cli_alert_info("Selecting top {N} combination{?s} based on {metric} ...")
  topCombinations <- selectTopCombinations(all_results, N = N, metric = metric)
  cli::cli_alert_success("Top combinations selected successfully!")

  # Plot the top combinations using the plotTopMetrics function
  cli::cli_h3("Top Combinations Plotting")
  cli::cli_alert_info("Generating plots of top {N} combination{?s} ...")
  p <- plotTopMetrics(topCombinations)
  print(p)
  cli::cli_alert_success("Top combinations plotted successfully!")

  # Count how often each gene appears in the top combinations
  cli::cli_h3("Gene Frequency in Top Combinations")
  geneCombinations <- gene_combination_matrix(combo_genes[topCombinations$Gene_Combination])
  geneFrequency <- gene_frequency(geneCombinations)
  print(plotGeneCombinations(geneCombinations))
  cli::cli_alert_success("Gene frequency plotted successfully!")

  # Store the results in a list
  results <- list(
    allResults = all_results,
    topCombinations = topCombinations,
    geneFrequency = geneFrequency,
    geneCombinations = geneCombinations
  )

  return(results)
}
