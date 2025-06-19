#' Plot PCA Scores and Top Loadings for Consensus Genes
#'
#' This function generates two plots based on precomputed PCA results:
#' (1) a scatterplot of sample scores on PC1 and PC2, colored by class labels,
#' and (2) a biplot of the top gene loadings for PC1 and PC2.
#'
#' @param pca_scores Data frame with columns `PC1`, `PC2`, and optionally sample labels.
#' @param pca_top_loadings Data frame with columns `PC1`, `PC2`, and `gene` for top loadings.
#' @param labels Vector or factor indicating class/group of each sample (used for coloring).
#'
#' @return A patchwork object combining the two ggplot2 plots: sample scores and top gene loadings.
#'
#' @details
#' Use this function after running PCA on your gene expression data and extracting the relevant
#' scores and top gene loadings. This function does not perform PCA; it only visualizes the results.
#'
#' @examples
#' \dontrun{
#' # After running PCA and extracting scores/loadings:
#' p <- pca.plot(pca_scores = my_scores, pca_top_loadings = my_top_loadings, labels = my_labels)
#' print(p)
#' }
#'
#' @import patchwork
#' @export
pca.plot <- function(pca_scores, pca_top_loadings, labels,
                     color_class_0 = "#535965",
                     color_class_1 = "#96CDCF",
                     arrow_color = "#595959") {

  # create the plot
  pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = factor(labels))) +
    scale_color_manual(values = c(color_class_0, color_class_1)) +
    geom_point(size = 3) +
    stat_ellipse(linewidth = 0.7, linetype = 2, type = "norm") +
    stat_ellipse(type = "t") +
    theme_minimal() +
    coord_fixed(ratio = 1) +
    labs(colour = "Class") +
    theme(
      legend.box.just = "right", legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 10, color = "#8f8f8f"),
      axis.title = element_text(size = 13, color = "#8f8f8f"),
    )

  # Plot of top loadings for each of PC1 and PC2
  loadings_plot <- ggplot(data = pca_top_loadings) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.08, "in")),
                 linewidth = 0.8,
                 colour = arrow_color
    ) +
    ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
                             nudge_y = 0.001, size = 2.5,
                             max.overlaps=Inf
    ) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    theme_minimal() +
    coord_fixed(ratio = 1) +
    labs(x = "PC1", y = "PC2") +
    theme(
      axis.text = element_text(size = 10, color = "#8f8f8f"),
      axis.title = element_text(size = 13, color = "#8f8f8f"),
      panel.grid.minor = element_blank()
    )

  return(list(pca_plot = pca_plot, loadings_plot = loadings_plot))
}

#' Plot AUROC and Fold Change for Genes
#'
#' Generates a ggplot2 plot displaying the AUROC (Area Under the ROC Curve) values for each gene,
#' with confidence intervals and fold change (FC) as a color scale. Useful for visualizing
#' the classification power and direction of regulation for consensus or selected genes.
#'
#' @param results A data frame with at least the following columns:
#'   - `gene`: Character, gene name.
#'   - `auroc`: Numeric, AUROC value for each gene.
#'   - `auroc_lower`: Numeric, lower confidence interval for AUROC.
#'   - `auroc_upper`: Numeric, upper confidence interval for AUROC.
#'   - `FC`: Numeric, fold change value for each gene.
#'
#' @return A ggplot2 object visualizing the AUROC (with confidence interval bars) and fold change for each gene.
#'
#' @details
#' Genes are ordered by AUROC. The color gradient represents the fold change, from blue (downregulated) through grey (neutral) to red (upregulated).
#' The horizontal dashed line at 0.5 indicates random classification performance.
#'
#' @examples
#' \dontrun{
#' # Example results data frame (results)
#' # results <- data.frame(
#' #   gene = c("GeneA", "GeneB"),
#' #   auroc = c(0.92, 0.85),
#' #   auroc_lower = c(0.88, 0.82),
#' #   auroc_upper = c(0.96, 0.88),
#' #   FC = c(1.5, -1.2)
#' # )
#' p <- auroc.fc.plot(results)
#' print(p)
#' }
#'
#' @export
auroc.fc.plot <- function(results,
                          low_color = "#187498",
                          middle_color = "#E5E5E5",
                          high_color = "#C62E2E"
) {

  p <- ggplot(results, aes(x = reorder(gene, auroc), y = auroc)) +
    geom_hline(yintercept = 0.5, linetype = "solid", color = "#7e7e7e", linewidth = 0.4) +
    geom_errorbar(aes(ymin = auroc_lower, ymax = auroc_upper, color = FC), width = 0, linewidth = 0.4, position = position_dodge(0.5)) +
    geom_point(aes(color = FC), size = 5, position = position_dodge(0.5)) +
    scale_color_gradient2(low = low_color, mid = middle_color, high = high_color, midpoint = 0, limits = c(-2.1, 2.1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "", y = "AUROC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1, color = "#7e7e7e"),
      axis.text.y = element_text(hjust = 1, size = 13, color = "#7e7e7e"),
      axis.title.y = element_text(size = 16, color = "#7e7e7e"),
      legend.title = element_text(size = 12, vjust = 0.5, hjust = 0),
      legend.text = element_text(size = 11, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks = element_line(color = "#7e7e7e", size = 0.3),
      axis.line.x = element_line(color = "#7e7e7e", linewidth = 0.3),
      axis.line.y = element_line(color = "#7e7e7e", linewidth = 0.3)
    )

  return(p)

}

#' Plot Hierarchical Clustering Heatmap for Consensus Genes
#'
#' Generates an interactive or static heatmap with hierarchical clustering for a matrix of gene expression values,
#' typically representing consensus genes across samples. This function wraps `heatmaply::heatmaply` with customization options.
#'
#' @param heatmap_data A list with:
#'   - `matrix`: Numeric matrix of expression values (genes x samples or vice-versa).
#'   - `hc_rows`: Hierarchical clustering object for rows (e.g., output of `hclust`).
#'   - `hc_cols`: Hierarchical clustering object for columns (e.g., output of `hclust`).
#'   - `side_colors` (optional): Vector or matrix for coloring sample columns.
#' @param colors Color palette (vector of colors or function) for the heatmap.
#' @param dendrogram Which dendrograms to display. One of `"both"`, `"row"`, `"column"`, `"none"`. Default: `"both"`.
#' @param show_dendrogram Logical vector of length 2; whether to show row and column dendrograms (e.g., `c(FALSE, TRUE)`).
#' @param scale How to scale the data: `"none"`, `"row"`, or `"column"`. Default: `"row"`.
#' @param custom_colors Named vector for coloring groups in side bar (e.g., `c("1"="#2e1457", "0"="#66a182")`).
#' @param margins Numeric vector of margins (top, right, bottom, left).
#' @param grid_color Color of the heatmap grid lines.
#' @param grid_width Width of the grid lines.
#' @param branches_lwd Line width for dendrogram branches.
#' @param fontsize_row Font size for row labels.
#' @param fontsize_col Font size for column labels.
#' @param scale_fill_gradient_fun Custom function for color scaling (advanced, optional).
#' @param limits Numeric vector of length 2, setting min and max of color scale (optional).
#' @param na.value Color for NA values.
#' @param cellnote Optional matrix of text to display in each cell.
#' @param cellnote_size Font size for cell notes (if provided).
#' @param cellnote_textposition Position for cell note text (if provided).
#' @param key.title Title for the heatmap color key.
#' @param key.xlab X-axis label for the heatmap color key.
#' @param key.ylab Y-axis label for the heatmap color key.
#'
#' @return An interactive heatmaply htmlwidget visualizing hierarchical clustering and expression values.
#'
#' @details
#' This function returns an interactive heatmap, allowing zoom, tooltip and dynamic exploration of gene/sample clustering.
#'
#' @examples
#' \dontrun{
#' p <- heatmap.plot(heatmap_data)
#' # In RStudio, p will display interactively in the Viewer pane.
#' }
#'
#' @import heatmaply
#' @export
heatmap.plot <- function(heatmap_data,
                         colors = grDevices::colorRampPalette(c("#2E4057", "#66A182", "#EDAE49"))(64),
                         dendrogram = "both",
                         show_dendrogram = c(FALSE, TRUE),
                         scale = "row",
                         custom_colors = c("1" = "#535965", "0" = "#96CDCF"),
                         margins = c(60, 100, 40, 20),
                         grid_color = "white",
                         grid_width = 0.00001,
                         branches_lwd = 0.4,
                         fontsize_row = 12,
                         fontsize_col = 5,
                         scale_fill_gradient_fun = NULL,
                         limits = NULL,
                         na.value = "grey50",
                         cellnote = NULL,
                         cellnote_size = NULL,
                         cellnote_textposition = "middle center",
                         key.title = "",
                         key.xlab = "Value",
                         key.ylab = "Frequency") {

  # Creating the heatmap with heatmaply
  p <- heatmaply::heatmaply(
    heatmap_data$matrix,
    Rowv = stats::as.dendrogram(heatmap_data$hc_rows),
    Colv = stats::as.dendrogram(heatmap_data$hc_cols),
    dendrogram = dendrogram,
    show_dendrogram = show_dendrogram,
    xlab = "", ylab = "", main = "",
    scale = scale, margins = margins,
    grid_color = grid_color, grid_width = grid_width,
    titleX = FALSE,
    hide_colorbar = FALSE,
    branches_lwd = branches_lwd,
    fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    showticklabels = c(FALSE, TRUE),
    labRow = rownames(heatmap_data$matrix),
    plot_method = "ggplot",
    colors = colors, col_side_colors = heatmap_data$side_colors, col_side_palette = custom_colors,
    scale_fill_gradient_fun = scale_fill_gradient_fun,
    limits = limits,
    na.value = na.value,
    cellnote = cellnote, cellnote_size = cellnote_size, cellnote_textposition = cellnote_textposition,
    key.title = key.title, key.xlab = key.xlab, key.ylab = key.ylab
  )

  return(p)
}

#' Plot Variable Importance for MLP Model with Performance Metrics
#'
#' Generates a bar plot visualizing the scaled importance of variables (genes) for a trained multilayer perceptron (MLP) model.
#' The plot also annotates the key model performance metrics: Accuracy, AUROC, and Brier Score.
#'
#' @param importances A data frame with variable importance results, typically with columns:
#'   - `Variable`: Character, gene/feature name.
#'   - `Importance`: Numeric, scaled variable importance (typically from -1 to 1).
#'   - `type`, `type_hjust`: Numeric, used internally for text positioning (optional, if precomputed).
#' @param test_performance A data frame with at least three performance metrics for the MLP model, expected columns:
#'   - `.estimate`: Numeric, containing values for Accuracy, AUROC, and Brier Score (in this order).
#'
#' @return A `ggplot2` object showing variable importance bars for the MLP model, with an annotation of the main performance metrics.
#'
#' @details
#' The bar plot displays the scaled importance of each variable, colored by sign (negative/positive).
#' The bottom annotation block displays MLP performance: Accuracy, AUROC, and Brier Score as calculated on the test set.
#'
#' @examples
#' \dontrun{
#' # Example importances data frame:
#' # importances <- data.frame(
#' #   Variable = c("Gene1", "Gene2"),
#' #   Importance = c(0.8, -0.6),
#' #   type = c(-0.026, 0.026),
#' #   type_hjust = c(1, 0)
#' # )
#' # test_performance <- data.frame(
#' #   .estimate = c(0.89, 0.95, 0.12),
#' #   .metric = c("accuracy", "roc_auc", "brier")
#' # )
#' p <- mlp.model.plot(importances, test_performance)
#' print(p)
#' }
#'
#' @import ggplot2
#' @export
mlp.model.plot <- function(importances,
                           test_performance,
                           colors = grDevices::colorRampPalette(c("#2E4057", "#66A182", "#EDAE49"))(64)
){

  # Create the plot
  cli::cli_alert_info("Creating plot ...")
  p <- ggplot(
    importances,
    aes(x = stats::reorder(Variable, Importance), y = Importance, label = Variable)
  ) +
    geom_segment(aes(y = 0, x = Variable, yend = Importance, xend = Variable), linewidth = 0.8, color = "gray") +
    geom_point(
      aes(color = Importance),
      stat = "identity",
      size = 4.5
    ) +
    geom_text(aes(label = Variable, y = type, hjust = type_hjust), color = "black", size = 3.5) +
    # geom_text(aes(label = round(Importance, 2), vjust = -2.1), color = "#7e7e7e", size = 3.5) +
    annotate("text",
             x = length(unique(importances$Variable)) - length(unique(importances$Variable)) * 0.1,
             y = -0.98, label = paste0(
               "MLP metrics: \n",
               "Accuracy: ", round(test_performance$.estimate[[1]], 3), "\n",
               "AUROC: ", round(test_performance$.estimate[[2]], 3), "\n",
               "Brier Score: ", round(test_performance$.estimate[[3]], 3)
             ),
             hjust = 0, vjust = 0.5, size = 5, color = "#7e7e7e"
    ) +
    labs(x = "", y = "Scaled importance") +
    theme_minimal() +
    guides(y = "none") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_text(size = 14, color = "#7e7e7e"),
      axis.text.x = element_text(size = 12, color = "#7e7e7e")
    ) +
    scale_color_gradientn(colors = colors) +
    coord_flip()
  cli::cli_alert_success("MLP Plot created successfully!")

  return(p)

}

