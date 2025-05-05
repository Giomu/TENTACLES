#' @title Plot Heatmap of Model Performance Metrics
#'
#' @description
#' Creates a heatmap using `ggplot2` to visualize the performance of different
#' clustering and dimensionality reduction models across multiple evaluation metrics,
#' for each gene combination. Each combination is labeled with a concise and readable
#' gene summary on the y-axis.
#'
#' @param top_results A data frame returned by `selectTopCombinations()`, containing
#' performance metrics of multiple models for each gene combination.
#' @param title A character string specifying the plot title. Default is "Top Model Metrics".
#'
#' @return A `ggplot2` object representing the heatmap.
#'
#' @examples
#' \dontrun{
#' vc <- valConsensus(df.count = count_table, gene.list = gene_list,
#'                    labels = labels, N = 10, metric = "FScore")
#' p <- plotTopMetrics2(vc$topCombinations)
#' print(p)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#' @importFrom ggplot2 facet_wrap theme_minimal labs theme element_text
#' @importFrom ggplot2 coord_flip scale_x_discrete scale_x_continuous scale_y_discrete
#' @importFrom dplyr left_join select rename
#' @importFrom tidyr pivot_longer

plotTopMetrics2 <- function(top_results, title = "Top Model Metrics") {
  # Extract unique gene combinations
  unique_combos <- unique(top_results$Gene_Combination)

  # Generate pretty gene labels (first 3 genes + ... + count)
  pretty_names <- vapply(
    strsplit(unique_combos, "_"),
    function(x) {
      label <- if (length(x) <= 3) paste(x, collapse = ", ") else paste(c(x[1:3], "..."), collapse = ", ")
      paste0(label, " (", length(x), " genes)")
    },
    character(1)
  )

  # Create short labels with leading zeros using R base
  padded_labels <- paste0("comb", sprintf("%0*d", nchar(length(unique_combos)), seq_along(unique_combos)))

  # Map: short label <-> pretty gene name
  gene_mapping <- data.frame(
    Gene_Combination = unique_combos,
    Short_Label = padded_labels,
    Pretty_Label = pretty_names,
    stringsAsFactors = FALSE
  )

  # Join short labels into top_results
  top_results <- dplyr::left_join(top_results, gene_mapping, by = "Gene_Combination")

  # Convert data to long format for ggplot
  long_results <- top_results %>%
    tidyr::pivot_longer(
      tidyselect::matches("^(KMeans|GMM|HC|PCA|tSNE|UMAP)_"),
      names_to = c("Model", "Metric"),
      names_sep = "_",
      values_to = "Value"
    )

  # Create heatmap
  heatmap <- ggplot(long_results, aes(x = Pretty_Label, y = Model, fill = Value)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradient2(low = "white", high = "steelblue") +
    facet_wrap(~Metric) +
    geom_text(aes(label = round(Value, 2)), color = "black", size = 3) +
    theme_minimal() +
    labs(
      x = NULL,
      y = NULL,
      title = title,
      fill = "Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "#7c7b7b"),
      axis.text.y = element_text(size = 7, color = "black"),
      panel.grid = element_blank(),
      strip.text = element_text(size = 13)
    ) +
    coord_flip()

  return(heatmap)

}
