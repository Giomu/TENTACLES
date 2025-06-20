#' @title Plot Heatmap of Model Performance Metrics
#'
#' @description
#' Creates a heatmap using `ggplot2` to visualize the performance of different
#' clustering and dimensionality reduction models across multiple evaluation metrics,
#' for each gene combination. Each combination is labeled with a concise and readable
#' gene summary on the y-axis, prefixed by a combination number (C1, C2, ...).
#'
#' @param top_results A data frame returned by `selectTopCombinations()`, containing
#' performance metrics of multiple models for each gene combination.
#' @param title A character string specifying the plot title. Default is "Top Model Metrics".
#' @param low_color Color for the lowest values in the heatmap gradient. Default: "white".
#' @param high_color Color for the highest values in the heatmap gradient. Default: "steelblue".
#' @param limits Numeric vector of length 2, specifying the minimum and maximum values for the fill scale. Default: c(0, 1).
#'
#' @return A `ggplot2` object representing the heatmap.
#'
#' @examples
#' \dontrun{
#' vc <- valConsensus(df.count = count_table, gene.list = gene_list,
#'                    labels = labels, N = 10, metric = "FScore")
#' p <- plotTopMetrics(
#'   vc$topCombinations,
#'   title = "Top Model Metrics",
#'   low_color = "white",
#'   high_color = "steelblue",
#'   limits = c(0, 1)
#' )
#' print(p)
#' }
#'
#' @export
#' @import ggplot2

plotTopMetrics <- function(top_results,
                           title = "Top Model Metrics",
                           low_color = "white",
                           high_color = "steelblue",
                           limits = c(0,1)) {

  if (!"Gene_Combination" %in% names(top_results)) {
    cli::cli_abort("Input data frame must contain a 'Gene_Combination' column.")
  }

  # Extract unique gene combinations
  unique_combos <- unique(top_results$Gene_Combination)

  # Generate pretty gene labels with combination number (C1, C2, ...)
  pretty_names <- vapply(
    seq_along(unique_combos),
    function(i) {
      genes <- strsplit(unique_combos[i], "_")[[1]]
      if (length(genes) > 3) {
        label <- paste(c(genes[1:3], "..."), collapse = ", ")
      } else {
        label <- paste(genes, collapse = ", ")
      }

      full_label <- paste0("C", i, ": ", label, " (", length(genes), " genes)")

      # Adjust if the label exceeds 30 characters
      if (nchar(full_label) > 30 && length(genes) > 2) {
        label <- paste(c(genes[1:2], "..."), collapse = ", ")
        full_label <- paste0("C", i, ": ", label, " (", length(genes), " genes)")
      }

      return(full_label)
    },
    character(1)
  )

  # Create mapping table
  gene_mapping <- data.frame(
    Gene_Combination = unique_combos,
    Pretty_Label = pretty_names,
    stringsAsFactors = FALSE
  )

  # Join labels into top_results using base R
  top_results <- merge(top_results, gene_mapping, by = "Gene_Combination", all.x = TRUE)

  # Convert data to long format
  reshape_cols <- grep("^(KMeans|GMM|HC|PCA|tSNE|UMAP)_", names(top_results), value = TRUE)
  long_results <- data.frame(
    Gene_Combination = rep(top_results$Gene_Combination, each = length(reshape_cols)),
    Pretty_Label = rep(top_results$Pretty_Label, each = length(reshape_cols)),
    Model = sub("_.*", "", reshape_cols),
    Metric = sub(".*_", "", reshape_cols),
    Value = as.vector(t(top_results[reshape_cols]))
  )

  # Define order of gene combinations for plotting
  level_order <- rev(unique(gene_mapping$Pretty_Label))

  # Create heatmap with ordered axis
  heatmap <- ggplot(long_results, aes(x = factor(Pretty_Label, levels = level_order), y = Model, fill = Value)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradient2(low = low_color,
                         high = high_color,
                         limits = limits) +
    facet_wrap(~Metric) +
    geom_text(aes(label = round(Value, 2)), color = "black", size = 3) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = title, fill = "Score") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "#7c7b7b"),
      axis.text.y = element_text(size = 7, color = "black"),
      panel.grid = element_blank(),
      strip.text = element_text(size = 13)
    ) +
    coord_flip()

  return(heatmap)
}
