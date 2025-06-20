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
selectTopCombinations <- function(results, N, metric) {
  allowed_metrics <- c("Accuracy", "Precision", "Recall", "FScore")
  if (!metric %in% allowed_metrics) stop("Invalid metric selected.")

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
#'
# Helper function for assessment based on a choose metric
evaluate_one_side <- function(pred, truth, metric = metric) {
  df <- data.frame(
    labels = factor(truth),
    cluster = factor(pred),
    cluster_inv = factor(1 - as.numeric(as.character(pred)))
  )

  metric_fun <- switch(metric,
                       Accuracy = yardstick::accuracy,
                       Precision = yardstick::precision,
                       Recall = yardstick::recall,
                       FScore = yardstick::f_meas,
                       stop("Unsupported metric")
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

#' @title valConsensus
#' @description This function validates the consensus genes in a provided dataset using clustering methods.
#'
#' @param df.count A numeric matrix containing the gene expression data. Rows represent samples and columns represent genes.
#' @param gene.list A character vector specifying the consensus genes to be evaluated.
#' @param class A factor vector specifying the true class labels for the samples.
#' @param N An integer specifying the number of top gene combinations to select based on the metric.
#' @param metric A character string specifying the metric to use for selecting the top gene combinations.
#' Possible values are 'Accuracy', 'Precision', 'Recall', and 'FScore'.
#'
#' @return A list containing the clustering results for each combination of genes and clustering methods.
#'
#' @details
#' The function evaluates the clustering performance of different combinations of genes using various clustering methods.
#' The function calculates the accuracy, precision, recall, and F1 score for each clustering method.
#' The function implements the following clustering methods: K-Means, Gaussian Mixture Model (GMM),
#' Hierarchical Clustering, k-Means on PCA dimensions, k-Means on t-SNE dimensions, and k-Means on UMAP dimensions.
#'
#' @importFrom mclust Mclust mclustBIC
#' @importFrom yardstick accuracy precision recall f_meas
#'
#' @examples
#' \dontrun{
#' count_table <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' gene_list <- c("Gene1", "Gene2", "Gene3")
#' class <- factor(sample(1:2, 100, replace = TRUE))
#' vc <- valConsensus(
#'   df.count = count_table, gene.list = gene_list,
#'   class = class, N = 10, metric = "FScore"
#' )
#' }
#'
#' @export
valConsensus <- function(df.count, gene.list, class, N = 10, metric = "FScore") {

  cli::cli_h1("Validating Consensus Genes")

  gene_list <- gene.list
  labels <- class

  # Maximum number of genes to consider
  max_genes <- length(unique(gene_list))

  # List to save results for all combinations
  all_results <- list()

  # Compute total number of combinations
  total_combinations <- as.numeric(sum(sapply(2:max_genes, function(n) choose(length(gene_list), n))))
  cli::cli_alert_info("Applying models to {total_combinations} combinations ...")
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
        }
      }
    }
  )

  # End Progress bar
  cli::cli_progress_done()
  cli::cli_alert_success("All combinations tested successfully!")

  # Add a new element to the list to store the top genes combinations by running the selectTopCombinations function
  cli::cli_h3("Top Combinations Selection")
  cli::cli_alert_info("Selecting top {N} combinations based on {metric} ...")
  topCombinations <- selectTopCombinations(all_results, N = N, metric = metric)
  cli::cli_alert_success("Top combinations selected successfully!")

  # Plot the top combinations using the plotTopMetrics function
  cli::cli_h3("Top Combinations Plotting")
  cli::cli_alert_info("Generating plots of top {N} combinations ...")
  p <- plotTopMetrics(topCombinations)
  print(p)
  cli::cli_alert_success("Top combinations plotted successfully!")

  # Store the results in a list
  results <- list(
    allResults = all_results,
    topCombinations = topCombinations
  )

  return(results)
}
