selectTopCombinations <- function(results, N, metric) {
  allowed_metrics <- c("Accuracy", "Precision", "Recall", "FScore")
  if (!metric %in% allowed_metrics) stop("Invalid metric selected.")

  methods <- names(results[[1]])  # assume all combinations have same methods
  metric_names <- paste0(rep(methods, each = length(allowed_metrics)), "_", allowed_metrics)

  purrr::map_dfr(names(results), function(combo_name) {
    combo_results <- results[[combo_name]]

    # Extract metric values from all clustering methods for this gene combination
    metric_values <- unlist(lapply(combo_results, `[[`, metric))

    # Compute the average of the selected metric across methods
    mean_metric <- mean(metric_values)

    row_data <- unlist(lapply(methods, function(m) {
      unlist(combo_results[[m]][allowed_metrics])
    }))

    tibble::tibble(
      Gene_Combination = combo_name,
      Mean_Metric = mean_metric,
      !!!setNames(as.list(row_data), metric_names)
    )
  }) %>%
    dplyr::arrange(dplyr::desc(Mean_Metric)) %>%
    dplyr::slice_head(n = N)
}

# Helper function for assessment based on a choose metric
evaluate_one_side <- function(pred, truth, metric = metric) {
  tib <- tibble::tibble(
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

  val_cluster     <- metric_fun(tib, truth = labels, estimate = cluster, event_level = "second")[[".estimate"]]
  val_cluster_inv <- metric_fun(tib, truth = labels, estimate = cluster_inv, event_level = "second")[[".estimate"]]

  best_col <- if (val_cluster >= val_cluster_inv) "cluster" else "cluster_inv"
  tib_pred <- tibble::tibble(labels = tib$labels, pred = tib[[best_col]])

  return(list(
    clusters  = as.integer(as.character(tib_pred$pred)),
    Accuracy  = yardstick::accuracy(tib_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
    Precision = yardstick::precision(tib_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
    Recall    = yardstick::recall(tib_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
    FScore    = yardstick::f_meas(tib_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]]
  ))
}

#' @title valConsensus
#' @description This function validates the consensus genes in a provided dataset using clustering methods.
#'
#' @param df.count A numeric matrix containing the gene expression data. Rows represent samples and columns represent genes.
#' @param gene.list A character vector specifying the consensus genes to be evaluated.
#' @param labels A factor vector specifying the true class labels for the samples.
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
#' @importFrom mclust Mclust
#' @importFrom yardstick accuracy precision recall f_meas
#'
#' @examples
#' \dontrun{
#' count_table <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' gene_list <- c("Gene1", "Gene2", "Gene3")
#' labels <- factor(sample(1:2, 100, replace = TRUE))
#' vc <- valConsensus(
#'   df.count = count_table, gene.list = gene_list,
#'   labels = labels, N = 10, metric = "FScore"
#' )
#' }
#'
#' @export
valConsensus <- function(df.count, gene.list, labels, N = 10, metric = "FScore") {

  cli::cli_h1("Validating Consensus Genes")

  gene_list <- gene.list
  labels <- labels

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
