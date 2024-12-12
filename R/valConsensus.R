utils::globalVariables("Mean_Metric")

selectTopCombinations <- function(results, N, metric) {
  # Initialize an empty list to store the flattened results
  flattened_results <- list()

  # Loop through each gene combination
  for (combo_name in names(results)) {
    combo_results <- results[[combo_name]]

    # Calculate the mean of the selected metric across all methods
    total_metric <- mean(sapply(combo_results, function(method_result) method_result[[metric]]))

    # Store the data, including all metrics for each method
    flattened_results[[length(flattened_results) + 1]] <- data.frame(
      Gene_Combination = combo_name,
      Mean_Metric = total_metric,
      KMeans_Accuracy = combo_results$KMeans[["Accuracy"]],
      KMeans_Precision = combo_results$KMeans[["Precision"]],
      KMeans_Recall = combo_results$KMeans[["Recall"]],
      KMeans_FScore = combo_results$KMeans[["FScore"]],
      GMM_Accuracy = combo_results$GMM[["Accuracy"]],
      GMM_Precision = combo_results$GMM[["Precision"]],
      GMM_Recall = combo_results$GMM[["Recall"]],
      GMM_FScore = combo_results$GMM[["FScore"]],
      HC_Accuracy = combo_results$HC[["Accuracy"]],
      HC_Precision = combo_results$HC[["Precision"]],
      HC_Recall = combo_results$HC[["Recall"]],
      HC_FScore = combo_results$HC[["FScore"]],
      PCA_Accuracy = combo_results$PCA[["Accuracy"]],
      PCA_Precision = combo_results$PCA[["Precision"]],
      PCA_Recall = combo_results$PCA[["Recall"]],
      PCA_FScore = combo_results$PCA[["FScore"]],
      tSNE_Accuracy = combo_results$tSNE[["Accuracy"]],
      tSNE_Precision = combo_results$tSNE[["Precision"]],
      tSNE_Recall = combo_results$tSNE[["Recall"]],
      tSNE_FScore = combo_results$tSNE[["FScore"]],
      UMAP_Accuracy = combo_results$UMAP[["Accuracy"]],
      UMAP_Precision = combo_results$UMAP[["Precision"]],
      UMAP_Recall = combo_results$UMAP[["Recall"]],
      UMAP_FScore = combo_results$UMAP[["FScore"]],
      stringsAsFactors = FALSE
    )
  }

  # Combine all rows into a single dataframe
  final_df <- dplyr::bind_rows(flattened_results)

  # Select the top N combinations based on the mean of the selected metric
  top_combinations <- final_df %>%
    dplyr::arrange(desc(Mean_Metric)) %>%
    dplyr::slice_head(n = N)

  return(top_combinations)
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
#' @importFrom MLmetrics Accuracy Precision Recall F1_Score
#'
#' @examples
#' \dontrun{
#' count_table <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' gene_list <- c("Gene1", "Gene2", "Gene3")
#' labels <- factor(sample(1:2, 100, replace = TRUE))
#' vc <- valConsensus(
#'   df.count = count_table, gene.list = gene_list,
#'   class = labels, N = 10, metric = "FScore"
#' )
#' }
#'
#' @export
valConsensus <- function(df.count, gene.list, class, N = 10, metric = "FScore") {
  gene_list <- gene.list

  cli::cli_h1("Validating Consensus Genes")
  # Numero massimo di geni da considerare
  max_genes <- length(unique(gene_list))
  labels <- class

  # Lista per salvare i risultati per tutte le combinazioni
  all_results <- list()

  # Compute total number of combinations
  total_combinations <- as.integer(sum(sapply(2:max_genes, function(n) choose(length(gene_list), n))))
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
          # Subset la tabella con i geni specificati
          count_subset <- df.count[, combo, drop = FALSE]

          # Ensure data is numeric and scaled
          count_subset <- as.data.frame(lapply(count_subset, as.numeric))
          count_subset <- scale(count_subset)
          rownames(count_subset) <- rownames(df.count)

          clustering_results <- list()

          # K-Means Clustering
          kmeans_result <- stats::kmeans(count_subset, centers = 2, iter.max = 100, algorithm = "MacQueen")
          clustering_results$KMeans <- list(
            clusters = kmeans_result$cluster,
            # avg.silwidth = fpc::cluster.stats(stats::dist(count_subset), kmeans_result$cluster)$avg.silwidth,
            Accuracy = max(MLmetrics::Accuracy(labels, kmeans_result$cluster), MLmetrics::Accuracy(labels, 3 - kmeans_result$cluster)),
            Precision = max(MLmetrics::Precision(labels, kmeans_result$cluster, positive = "1"), MLmetrics::Precision(labels, 3 - kmeans_result$cluster, positive = "1")),
            Recall = max(MLmetrics::Recall(labels, kmeans_result$cluster, positive = "1"), MLmetrics::Recall(labels, 3 - kmeans_result$cluster, positive = "1")),
            FScore = max(MLmetrics::F1_Score(labels, kmeans_result$cluster, positive = "1"), MLmetrics::F1_Score(labels, 3 - kmeans_result$cluster, positive = "1"))
          )

          # Gaussian Mixture Model (GMM)
          gmm_result <- mclust::Mclust(count_subset, G = 2, verbose = FALSE)
          clustering_results$GMM <- list(
            clusters = gmm_result$classification,
            # avg.silwidth = fpc::cluster.stats(stats::dist(count_subset), gmm_result$classification)$avg.silwidth,
            Accuracy = max(MLmetrics::Accuracy(labels, gmm_result$classification), MLmetrics::Accuracy(labels, 3 - gmm_result$classification)),
            Precision = max(MLmetrics::Precision(labels, gmm_result$classification, positive = "1"), MLmetrics::Precision(labels, 3 - gmm_result$classification, positive = "1")),
            Recall = max(MLmetrics::Recall(labels, gmm_result$classification, positive = "1"), MLmetrics::Recall(labels, 3 - gmm_result$classification, positive = "1")),
            FScore = max(MLmetrics::F1_Score(labels, gmm_result$classification, positive = "1"), MLmetrics::F1_Score(labels, 3 - gmm_result$classification, positive = "1"))
          )

          # Hierarchical Clustering
          hc_result <- stats::hclust(dist(count_subset), method = "complete")
          hc_clusters <- stats::cutree(hc_result, k = 2)
          clustering_results$HC <- list(
            clusters = hc_clusters,
            # avg.silwidth = fpc::cluster.stats(stats::dist(count_subset), hc_clusters)$avg.silwidth,
            Accuracy = max(MLmetrics::Accuracy(labels, hc_clusters), MLmetrics::Accuracy(labels, 3 - hc_clusters)),
            Precision = max(MLmetrics::Precision(labels, hc_clusters, positive = "1"), MLmetrics::Precision(labels, 3 - hc_clusters, positive = "1")),
            Recall = max(MLmetrics::Recall(labels, hc_clusters, positive = "1"), MLmetrics::Recall(labels, 3 - hc_clusters, positive = "1")),
            FScore = max(MLmetrics::F1_Score(labels, hc_clusters, positive = "1"), MLmetrics::F1_Score(labels, 3 - hc_clusters, positive = "1"))
          )

          # PCA + KMeans
          pca_result <- stats::prcomp(count_subset, scale = FALSE)
          pca_data <- data.frame(pca_result$x)
          pca_clusters <- stats::kmeans(pca_data[, 1:2], centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster
          clustering_results$PCA <- list(
            clusters = pca_clusters,
            # avg.silwidth = fpc::cluster.stats(stats::dist(pca_data[, 1:2]), pca_clusters)$avg.silwidth,
            Accuracy = max(MLmetrics::Accuracy(labels, pca_clusters), MLmetrics::Accuracy(labels, 3 - pca_clusters)),
            Precision = max(MLmetrics::Precision(labels, pca_clusters, positive = "1"), MLmetrics::Precision(labels, 3 - pca_clusters, positive = "1")),
            Recall = max(MLmetrics::Recall(labels, pca_clusters, positive = "1"), MLmetrics::Recall(labels, 3 - pca_clusters, positive = "1")),
            FScore = max(MLmetrics::F1_Score(labels, pca_clusters, positive = "1"), MLmetrics::F1_Score(labels, 3 - pca_clusters, positive = "1"))
          )

          # t-SNE + KMeans
          tsne_result <- Rtsne::Rtsne(count_subset, dims = 2, perplexity = 10)
          tsne_data <- tsne_result$Y
          tsne_clusters <- stats::kmeans(tsne_data, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster
          clustering_results$tSNE <- list(
            clusters = tsne_clusters,
            # avg.silwidth = fpc::cluster.stats(stats::dist(tsne_data), tsne_clusters)$avg.silwidth,
            Accuracy = max(MLmetrics::Accuracy(labels, tsne_clusters), MLmetrics::Accuracy(labels, 3 - tsne_clusters)),
            Precision = max(MLmetrics::Precision(labels, tsne_clusters, positive = "1"), MLmetrics::Precision(labels, 3 - tsne_clusters, positive = "1")),
            Recall = max(MLmetrics::Recall(labels, tsne_clusters, positive = "1"), MLmetrics::Recall(labels, 3 - tsne_clusters, positive = "1")),
            FScore = max(MLmetrics::F1_Score(labels, tsne_clusters, positive = "1"), MLmetrics::F1_Score(labels, 3 - tsne_clusters, positive = "1"))
          )

          # UMAP + KMeans
          umap_result <- umap::umap(count_subset)
          umap_data <- as.data.frame(umap_result$layout)
          umap_clusters <- stats::kmeans(umap_data, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster
          clustering_results$UMAP <- list(
            clusters = umap_clusters,
            # avg.silwidth = fpc::cluster.stats(stats::dist(umap_data), umap_clusters)$avg.silwidth,
            Accuracy = max(MLmetrics::Accuracy(labels, umap_clusters), MLmetrics::Accuracy(labels, 3 - umap_clusters)),
            Precision = max(MLmetrics::Precision(labels, umap_clusters, positive = "1"), MLmetrics::Precision(labels, 3 - umap_clusters, positive = "1")),
            Recall = max(MLmetrics::Recall(labels, umap_clusters, positive = "1"), MLmetrics::Recall(labels, 3 - umap_clusters, positive = "1")),
            FScore = max(MLmetrics::F1_Score(labels, umap_clusters, positive = "1"), MLmetrics::F1_Score(labels, 3 - umap_clusters, positive = "1"))
          )

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
