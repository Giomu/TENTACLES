# Helper function to pick the top N combinations: returns the indices of the N
# highest mean metrics, in decreasing order. Ties keep the order in which the
# combinations were tested, and NA means go last
selectTopCombinations <- function(mean_metric, N) {
  ranked <- order(mean_metric, decreasing = TRUE)
  ranked[seq_len(min(N, length(ranked)))]
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
#' highest value for the selected metric. An orientation whose metric is \code{NA} (e.g. Precision
#' when all samples are predicted as 0) loses to one where it can be computed.
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
  # Fixed levels, so a method that puts all samples in one cluster still gives
  # two-level factors, as yardstick requires
  df <- data.frame(
    labels = factor(truth, levels = c(0, 1)),
    cluster = factor(pred, levels = c(0, 1)),
    cluster_inv = factor(1 - as.numeric(as.character(pred)), levels = c(0, 1))
  )

  metric_fun <- switch(metric,
                       Accuracy  = yardstick::accuracy,
                       Precision = yardstick::precision,
                       Recall    = yardstick::recall,
                       FScore    = yardstick::f_meas,
                       cli::cli_abort("Unsupported metric: {metric}")
  )

  # With a single cluster, the orientation predicting no 1s has undefined
  # Precision and FScore (NA). That is expected, so yardstick's warning is muffled
  # and the NA is kept in the returned metrics
  withCallingHandlers({
    val_cluster     <- metric_fun(df, truth = labels, estimate = cluster, event_level = "second")[[".estimate"]]
    val_cluster_inv <- metric_fun(df, truth = labels, estimate = cluster_inv, event_level = "second")[[".estimate"]]

    # An orientation whose metric is NA loses to one where it can be computed
    best_col <- if (is.na(val_cluster_inv) || (!is.na(val_cluster) && val_cluster >= val_cluster_inv)) {
      "cluster"
    } else {
      "cluster_inv"
    }
    df_pred <- data.frame(labels = df$labels, pred = df[[best_col]])

    result <- list(
      clusters  = as.integer(as.character(df_pred$pred)),
      Accuracy  = yardstick::accuracy(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
      Precision = yardstick::precision(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
      Recall    = yardstick::recall(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]],
      FScore    = yardstick::f_meas(df_pred, truth = labels, estimate = pred, event_level = "second")[[".estimate"]]
    )
  }, yardstick_warning_precision_undefined = function(w) invokeRestart("muffleWarning"))

  return(result)
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

# Helper function to turn a gene x combination matrix into a data frame with
# one row per gene: the binary combination columns, then the number of
# combinations containing the gene and their identifiers
gene_frequency <- function(gene_matrix) {
  freq <- data.frame(Gene = rownames(gene_matrix), gene_matrix, row.names = NULL, check.names = FALSE)
  freq$Count <- as.integer(rowSums(gene_matrix))
  freq$Combinations <- unname(apply(gene_matrix, 1, function(x) {
    paste(colnames(gene_matrix)[x == 1], collapse = ", ")
  }))
  freq
}

#' @title valConsensus
#' @description This function validates the consensus genes in a provided dataset using clustering methods.
#'
#' @param df.count A data frame with the gene expression data, with samples in rows and genes in
#' columns. A matrix is not supported.
#' @param gene.list A character vector with the consensus genes to evaluate, at least two. All of them
#' must be columns of \code{df.count}.
#' @param class A vector with the true class of each sample (row of \code{df.count}), coded as 0 (control)
#' and 1 (case), with both classes present. It can be numeric or a factor with levels \code{"0"} and \code{"1"}.
#' Class 1 is the positive class for Precision, Recall and FScore.
#' @param N An integer specifying the number of top gene combinations to return. If it is larger than
#' the number of possible combinations, it is reduced to that number with a message.
#' @param metric A character string specifying the metric used to rank the combinations and to orient
#' the clusters (see Details). Possible values are \code{"Accuracy"}, \code{"Precision"}, \code{"Recall"},
#' and \code{"FScore"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{allResults}{A list with the results of all tested combinations:
#'   \describe{
#'     \item{trueClass}{An integer vector with the true class (0/1) of each sample, named by the
#'     row names of \code{df.count}.}
#'     \item{clusters}{An integer array with dimensions \code{[sample, combination, method]} holding
#'     the cluster (0/1) predicted for each sample, already assigned to the classes as described in
#'     Details. Combinations are named by their genes joined by \code{"_"} and methods are
#'     \code{"KMeans"}, \code{"GMM"}, \code{"HC"}, \code{"PCA"}, \code{"tSNE"} and \code{"UMAP"}.
#'     For example, \code{cbind(trueClass = trueClass, clusters[, "GENE1_GENE2", ])} gives the true
#'     and predicted classes of the combination of GENE1 and GENE2.}
#'     \item{metrics}{A data frame with one row per combination, in the order of the
#'     \code{combination} dimension of \code{clusters}. \code{Gene_Combination} holds the genes of the
#'     combination joined by \code{"_"}, \code{Mean_Metric} the mean of \code{metric} across the six
#'     methods, and the remaining columns, named \code{<Method>_<Metric>} (e.g.
#'     \code{KMeans_FScore}), the four metrics of each method, computed from \code{trueClass} and
#'     \code{clusters}.}
#'   }}
#'   \item{topCombinations}{The rows of \code{allResults$metrics} for the top \code{N} combinations,
#'   sorted by decreasing \code{Mean_Metric}.}
#'   \item{geneFrequency}{A data frame with one row per gene in the top combinations, sorted by
#'   decreasing \code{Count} and then by name. Columns: \code{Gene}; one binary column per top
#'   combination (\code{C1}, \code{C2}, ...), with 1 when the gene is in the combination; \code{Count},
#'   the number of combinations containing the gene; and \code{Combinations}, their identifiers
#'   (e.g. \code{"C1, C3"}). Combination \code{Ci} is the i-th row of \code{topCombinations}. The
#'   binary matrix can be extracted with
#'   \code{as.matrix(geneFrequency[grep("^C[0-9]+$", names(geneFrequency))])}.}
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
#' Precision, Recall and FScore are then computed for that assignment. If a method puts all samples
#' in a single cluster, the assignment where \code{metric} can be computed is kept; Precision and
#' FScore are \code{NA} when the kept assignment predicts no sample as 1, which can only happen
#' when \code{metric} is \code{"Accuracy"}.
#'
#' Combinations are ranked by the mean of \code{metric} across the six methods, and the top
#' \code{N} are returned. For these, the function counts how many combinations contain each gene
#' and prints two plots: the metrics heatmap (\code{\link{plotTopMetrics}}) and the gene presence
#' heatmap (\code{\link{plotGeneCombinations}}).
#'
#' The random seed is fixed at 123 inside the function, so results are reproducible but
#' cannot be varied. t-SNE with perplexity 10 requires at least 31 samples. The number of
#' combinations grows exponentially with the number of genes: 15 genes give 32,752
#' combinations. The \code{clusters} array takes 4 bytes per sample, combination and method,
#' about 79 MB for 15 genes and 100 samples; a warning is shown when it exceeds 1 GB.
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
  # Metrics need classes coded 0/1; other values would become NA and be dropped
  bad_classes <- setdiff(unique(as.character(class)), c("0", "1"))
  if (length(bad_classes) > 0) {
    cli::cli_abort(c(
      "{.arg class} must be coded as 0 (control) and 1 (case).",
      "x" = "Found {.val {bad_classes}}."
    ))
  }
  if (length(class) != nrow(df.count)) {
    cli::cli_abort(
      "{.arg class} has {length(class)} value{?s}, but {.arg df.count} has {nrow(df.count)} sample{?s}."
    )
  }
  # Clusters are compared with the classes, so both must be present
  if (length(unique(as.character(class))) < 2) {
    cli::cli_abort("{.arg class} must contain both classes, 0 and 1.")
  }

  allowed_metrics <- c("Accuracy", "Precision", "Recall", "FScore")
  if (!metric %in% allowed_metrics) {
    cli::cli_abort("{.arg metric} must be one of {.val {allowed_metrics}}, not {.val {metric}}.")
  }

  gene_list <- unique(gene.list)
  if (length(gene_list) < 2) {
    cli::cli_abort("{.arg gene.list} must contain at least two distinct genes.")
  }
  labels <- class
  methods <- c("KMeans", "GMM", "HC", "PCA", "tSNE", "UMAP")

  # All combinations of 2 or more genes, in the order they are tested. Genes are
  # kept as vectors because gene names may contain the "_" used to build the keys
  combo_genes <- unlist(
    lapply(2:length(gene_list), function(n) combn(gene_list, n, simplify = FALSE)),
    recursive = FALSE
  )
  combo_keys <- vapply(combo_genes, paste, character(1), collapse = "_")
  total_combinations <- length(combo_genes)

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
  # The clusters array holds one 4-byte integer per sample, combination and method
  clusters_gb <- nrow(df.count) * total_combinations * length(methods) * 4 / 1024^3
  if (clusters_gb > 1) {
    cli::cli_alert_warning("Storing the predicted clusters will take about {round(clusters_gb, 1)} GB of memory.")
  }

  # Results of all combinations, filled in place: predicted clusters (already
  # oriented by evaluate_one_side) and the four metrics of each method
  clusters <- array(
    NA_integer_,
    dim = c(nrow(df.count), total_combinations, length(methods)),
    dimnames = list(sample = rownames(df.count), combination = combo_keys, method = methods)
  )
  metric_values <- matrix(
    NA_real_,
    nrow = total_combinations,
    ncol = length(methods) * length(allowed_metrics),
    dimnames = list(NULL, paste0(rep(methods, each = length(allowed_metrics)), "_", allowed_metrics))
  )

  # Initialize progress bar
  progress_bar <- cli::cli_progress_bar(
    format = "Tested combinations: {cli::pb_bar} {cli::pb_percent} ({cli::pb_current}/{cli::pb_total})",
    total = total_combinations,
    clear = FALSE
  )

  withr::with_seed(
    seed = 123,
    code = {
      for (i in seq_len(total_combinations)) {
        combo <- combo_genes[[i]]
        # Increment progress bar
        cli::cli_progress_update(inc = 1)
        # Subset the table with the specified genes
        count_subset <- df.count[, combo, drop = FALSE]

        # Ensure data is numeric and scaled
        count_subset <- as.data.frame(lapply(count_subset, as.numeric))
        count_subset <- scale(count_subset)
        rownames(count_subset) <- rownames(df.count)

        # Predicted clusters (0/1) of each method. The methods run in a fixed
        # order, so the random draws, and the results, are reproducible
        predictions <- list()

        # K-Means Clustering
        kmeans_result <- stats::kmeans(count_subset, centers = 2, iter.max = 100, algorithm = "MacQueen")
        predictions$KMeans <- kmeans_result$cluster - 1

        # Gaussian Mixture Model (GMM)
        gmm_result <- mclust::Mclust(count_subset, G = 2, verbose = FALSE)
        predictions$GMM <- gmm_result$classification - 1

        # Hierarchical Clustering
        hc_result <- stats::hclust(dist(count_subset), method = "complete")
        predictions$HC <- stats::cutree(hc_result, k = 2) - 1

        # PCA + KMeans
        pca_result <- stats::prcomp(count_subset, scale = FALSE)
        predictions$PCA <- stats::kmeans(pca_result$x[, 1:2], centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1

        # t-SNE + KMeans
        tsne_result <- Rtsne::Rtsne(count_subset, dims = 2, perplexity = 10, check_duplicates = FALSE)
        predictions$tSNE <- stats::kmeans(tsne_result$Y, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1

        # UMAP + KMeans
        umap_result <- umap::umap(count_subset)
        predictions$UMAP <- stats::kmeans(umap_result$layout, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1

        # Orient the clusters, compute the metrics and save them for this combination
        for (m in methods) {
          evaluation <- evaluate_one_side(predictions[[m]], labels, metric)
          clusters[, i, m] <- evaluation$clusters
          metric_values[i, paste0(m, "_", allowed_metrics)] <- unlist(evaluation[allowed_metrics])
        }
      }
    }
  )

  # End Progress bar
  cli::cli_progress_done()
  cli::cli_alert_success("All combinations tested successfully!")

  all_metrics <- data.frame(
    Gene_Combination = combo_keys,
    Mean_Metric = rowMeans(metric_values[, paste0(methods, "_", metric), drop = FALSE]),
    metric_values,
    check.names = FALSE
  )

  cli::cli_h3("Top Combinations Selection")
  cli::cli_alert_info("Selecting top {N} combination{?s} based on {metric} ...")
  top_index <- selectTopCombinations(all_metrics$Mean_Metric, N = N)
  topCombinations <- all_metrics[top_index, , drop = FALSE]
  rownames(topCombinations) <- NULL
  cli::cli_alert_success("Top combinations selected successfully!")

  # Plot the top combinations using the plotTopMetrics function
  cli::cli_h3("Top Combinations Plotting")
  cli::cli_alert_info("Generating plots of top {N} combination{?s} ...")
  p <- plotTopMetrics(topCombinations)
  print(p)
  cli::cli_alert_success("Top combinations plotted successfully!")

  # Count how often each gene appears in the top combinations
  cli::cli_h3("Gene Frequency in Top Combinations")
  geneFrequency <- gene_frequency(gene_combination_matrix(combo_genes[top_index]))
  print(plotGeneCombinations(geneFrequency))
  cli::cli_alert_success("Gene frequency plotted successfully!")

  # Store the results in a list
  results <- list(
    allResults = list(
      trueClass = stats::setNames(as.integer(as.character(class)), rownames(df.count)),
      clusters = clusters,
      metrics = all_metrics
    ),
    topCombinations = topCombinations,
    geneFrequency = geneFrequency
  )

  return(results)
}
