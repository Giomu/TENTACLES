# Helper function to pick the top N combinations: returns the indices of the N
# highest mean metrics, in decreasing order. Ties keep the order in which the
# combinations were tested. NA means, from combinations where a method failed,
# are left out, so fewer than N indices may be returned
selectTopCombinations <- function(mean_metric, N) {
  ranked <- order(mean_metric, decreasing = TRUE, na.last = NA)
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

# Helper function to split the samples into two clusters (0/1) with a Gaussian
# mixture model. Mclust fits several covariance models and keeps the best by
# BIC. The inner M-step loop of models such as VEE has no practical iteration
# limit by default, and on degenerate data (e.g. a gene with many tied values)
# it never converges, so it is capped at 1000; such a model then returns NULL
# and is skipped, as the other degenerate models are. When every model fails,
# Mclust returns NULL, which becomes an error recorded by valConsensus
gmm_clusters <- function(x) {
  # With exactly two distinct points, both clusters have zero variance and no
  # model can be fitted; the two points are then the two clusters
  points <- unique(x)
  if (nrow(points) == 2) {
    return(as.numeric(apply(x, 1, function(row) all(row == points[2, ]))))
  }
  gmm_result <- mclust::Mclust(
    x, G = 2, verbose = FALSE,
    control = mclust::emControl(itmax = c(.Machine$integer.max, 1000))
  )
  if (is.null(gmm_result)) {
    stop("Mclust could not fit any covariance model.", call. = FALSE)
  }
  gmm_result$classification - 1
}

# Helper function to warn before a run with too many combinations or a clusters
# array too large for memory, which holds one 4-byte integer per sample,
# combination and method
warn_run_size <- function(n_samples, n_combinations, n_methods) {
  if (n_combinations > 1000000) {
    cli::cli_alert_warning("The number of combinations is > 1M. This may take a long time to complete. Consider reducing the number of genes.")
  }
  clusters_gb <- n_samples * n_combinations * n_methods * 4 / 1024^3
  if (clusters_gb > 1) {
    cli::cli_alert_warning("Storing the predicted clusters will take about {round(clusters_gb, 1)} GB of memory.")
  }
  invisible(clusters_gb)
}

# Helper function to split the scaled data of a combination into two clusters
# (0/1) with one of the valConsensus methods
cluster_samples <- function(x, method) {
  kmeans2 <- function(y) stats::kmeans(y, centers = 2, iter.max = 100, algorithm = "MacQueen")$cluster - 1
  switch(method,
    KMeans = kmeans2(x),
    GMM = gmm_clusters(x),
    HC = stats::cutree(stats::hclust(dist(x), method = "complete"), k = 2) - 1,
    PCA = kmeans2(stats::prcomp(x, scale = FALSE)$x[, 1:2]),
    tSNE = kmeans2(Rtsne::Rtsne(x, dims = 2, perplexity = 10, check_duplicates = FALSE)$Y),
    UMAP = kmeans2(umap::umap(x)$layout)
  )
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
#' columns. A matrix is not supported. The genes in \code{gene.list} must have no missing or infinite
#' values and must not be constant across samples, since each gene is scaled before clustering.
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
#' @param methods A character vector with the clustering methods to run: any of \code{"KMeans"},
#' \code{"GMM"}, \code{"HC"}, \code{"PCA"}, \code{"tSNE"} and \code{"UMAP"} (see Details). All six by
#' default. They always run in this order, whatever the order given.
#'
#' @return A list containing:
#' \describe{
#'   \item{allResults}{A list with the results of all tested combinations:
#'   \describe{
#'     \item{trueClass}{An integer vector with the true class (0/1) of each sample, named by the
#'     row names of \code{df.count}.}
#'     \item{clusters}{An integer array with dimensions \code{[sample, combination, method]} holding
#'     the cluster (0/1) predicted for each sample, already assigned to the classes as described in
#'     Details. Combinations are named by their genes joined by \code{"_"} and methods are those in
#'     \code{methods}, in the order \code{"KMeans"}, \code{"GMM"}, \code{"HC"}, \code{"PCA"},
#'     \code{"tSNE"}, \code{"UMAP"}.
#'     For example, \code{cbind(trueClass = trueClass, clusters[, "GENE1_GENE2", ])} gives the true
#'     and predicted classes of the combination of GENE1 and GENE2.}
#'     \item{metrics}{A data frame with one row per combination, in the order of the
#'     \code{combination} dimension of \code{clusters}. \code{Gene_Combination} holds the genes of the
#'     combination joined by \code{"_"}, \code{Mean_Metric} the mean of \code{metric} across the
#'     selected methods, and the remaining columns, named \code{<Method>_<Metric>} (e.g.
#'     \code{KMeans_FScore}), the four metrics of each method, computed from \code{trueClass} and
#'     \code{clusters}.}
#'     \item{failures}{Only present when a method failed on some combination (see Details): a
#'     data frame with one row per failure and columns \code{Gene_Combination}, \code{Method} and
#'     \code{Error}, the error message.}
#'   }}
#'   \item{topCombinations}{The rows of \code{allResults$metrics} for the top \code{N} combinations,
#'   sorted by decreasing \code{Mean_Metric}. Combinations where a method failed are left out, so
#'   there may be fewer than \code{N}.}
#'   \item{geneFrequency}{\code{NULL} if no combination ranks. Otherwise, a data frame with one row
#'   per gene in the top combinations, sorted by
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
#' standard deviation) and the samples are split into two clusters by each method in
#' \code{methods}:
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
#' Combinations are ranked by the mean of \code{metric} across the selected methods, so rankings
#' from runs with different \code{methods} are not comparable. The top
#' \code{N} are returned. For these, the function counts how many combinations contain each gene
#' and prints two plots: the metrics heatmap (\code{\link{plotTopMetrics}}) and the gene presence
#' heatmap (\code{\link{plotGeneCombinations}}).
#'
#' A method that fails on a combination does not stop the run: its clusters and metrics for that
#' combination are \code{NA}, the failure is recorded in \code{allResults$failures}, the
#' combination's \code{Mean_Metric} is \code{NA}, so it is left out of the ranking, and a warning
#' at the end reports the failures. To keep these combinations, run again without the failing
#' method. GMM fails when \code{mclust::Mclust()} cannot fit any covariance model; when the samples
#' of a combination take exactly two distinct values, they are split by value instead.
#'
#' The random seed is fixed at 123 inside the function, so results are reproducible but
#' cannot be varied. The methods that use random numbers (all but GMM and HC) draw them from a
#' single stream across all combinations, so their results change when \code{methods} changes.
#' t-SNE with perplexity 10 requires at
#' least 31 samples; with fewer, leave \code{"tSNE"} out of \code{methods}. The number of
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
valConsensus <- function(df.count, gene.list, class, N = 10, metric = "FScore",
                         methods = c("KMeans", "GMM", "HC", "PCA", "tSNE", "UMAP")) {

  cli::cli_h1("Validating Consensus Genes")

  # Fail before testing any combination, naming the missing features
  missing_genes <- setdiff(gene.list, colnames(df.count))
  if (length(missing_genes) > 0) {
    cli::cli_abort(c(
      "{length(missing_genes)} feature{?s} of {.arg gene.list} {?is/are} not in {.arg df.count}:",
      "x" = "{.val {missing_genes}}"
    ))
  }
  # Scaling turns a non-finite or constant feature into NaN, which breaks five of
  # the six methods; hierarchical clustering would silently ignore the feature
  gene_values <- lapply(df.count[unique(gene.list)], as.numeric)
  nonfinite_genes <- names(gene_values)[vapply(gene_values, function(x) any(!is.finite(x)), logical(1))]
  constant_genes <- names(gene_values)[vapply(gene_values, function(x) length(unique(x[is.finite(x)])) < 2, logical(1))]
  if (length(nonfinite_genes) > 0 || length(constant_genes) > 0) {
    cli::cli_abort(c(
      "Features of {.arg gene.list} must have finite values that vary across samples.",
      "x" = if (length(nonfinite_genes) > 0) "Missing or infinite values in {.val {nonfinite_genes}}.",
      "x" = if (length(constant_genes) > 0) "Same value in all samples: {.val {constant_genes}}."
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

  all_methods <- c("KMeans", "GMM", "HC", "PCA", "tSNE", "UMAP")
  if (!is.character(methods) || length(methods) == 0) {
    cli::cli_abort("{.arg methods} must be a character vector with at least one of {.val {all_methods}}.")
  }
  unknown_methods <- setdiff(methods, all_methods)
  if (length(unknown_methods) > 0) {
    cli::cli_abort(c(
      "{.arg methods} must be one of {.val {all_methods}}.",
      "x" = "Unknown method{?s}: {.val {unknown_methods}}."
    ))
  }
  # A fixed order keeps the random draws, and so the results, independent of
  # the order given
  methods <- intersect(all_methods, methods)
  # Rtsne requires nrow - 1 >= 3 * perplexity; checked here so the error comes
  # before testing any combination
  if ("tSNE" %in% methods && nrow(df.count) < 31) {
    cli::cli_abort(c(
      "t-SNE needs at least 31 samples, but {.arg df.count} has {nrow(df.count)}.",
      "i" = "Leave {.val tSNE} out of {.arg methods}."
    ))
  }

  gene_list <- unique(gene.list)
  if (length(gene_list) < 2) {
    cli::cli_abort("{.arg gene.list} must contain at least two distinct genes.")
  }
  labels <- class

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
  warn_run_size(nrow(df.count), total_combinations, length(methods))

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
  # Error message of each method that fails on a combination, NA otherwise
  error_messages <- matrix(NA_character_, nrow = total_combinations, ncol = length(methods),
                           dimnames = list(NULL, methods))

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

        # Cluster with each method, orient the clusters, compute the metrics and
        # save them for this combination. The methods run in a fixed order, so
        # the random draws, and the results, are reproducible. A failure leaves
        # NA results, so a long run is not lost to one combination
        for (m in methods) {
          evaluation <- tryCatch(
            evaluate_one_side(cluster_samples(count_subset, m), labels, metric),
            error = function(e) e
          )
          if (inherits(evaluation, "error")) {
            error_messages[i, m] <- conditionMessage(evaluation)
            next
          }
          clusters[, i, m] <- evaluation$clusters
          metric_values[i, paste0(m, "_", allowed_metrics)] <- unlist(evaluation[allowed_metrics])
        }
      }
    }
  )

  # End Progress bar
  cli::cli_progress_done()

  failed <- which(!is.na(error_messages), arr.ind = TRUE)
  failed <- failed[order(failed[, 1], failed[, 2]), , drop = FALSE]
  failures <- data.frame(
    Gene_Combination = combo_keys[failed[, 1]],
    Method = methods[failed[, 2]],
    Error = error_messages[failed],
    row.names = NULL
  )
  if (nrow(failures) == 0) {
    cli::cli_alert_success("All combinations tested successfully!")
  } else {
    n_failed_combos <- length(unique(failures$Gene_Combination))
    failed_methods <- unique(failures$Method)
    cli::cli_alert_warning("All combinations tested, but {nrow(failures)} clustering run{?s} failed.")
    cli::cli_warn(c(
      "{nrow(failures)} clustering run{?s} failed in {n_failed_combos} combination{?s}, left out of the ranking.",
      "i" = "Failed method{?s}: {.val {failed_methods}}.",
      "i" = "First error ({failures$Method[1]} on {failures$Gene_Combination[1]}): {failures$Error[1]}",
      "i" = "See {.code allResults$failures}, or leave the failing {qty(failed_methods)}method{?s} out of {.arg methods}."
    ))
  }

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
  if (length(top_index) < N) {
    cli::cli_alert_warning(
      "Only {length(top_index)} combination{?s} {?has/have} results from all methods, fewer than {.arg N} = {N}."
    )
    N <- length(top_index)
  }

  geneFrequency <- NULL
  if (N == 0) {
    cli::cli_alert_warning("No combination ranks, so nothing is plotted.")
  } else {
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
  }

  # Store the results in a list; failures only when a method failed
  all_results <- list(
    trueClass = stats::setNames(as.integer(as.character(class)), rownames(df.count)),
    clusters = clusters,
    metrics = all_metrics
  )
  if (nrow(failures) > 0) {
    all_results$failures <- failures
  }
  results <- list(
    allResults = all_results,
    topCombinations = topCombinations,
    geneFrequency = geneFrequency
  )

  return(results)
}
