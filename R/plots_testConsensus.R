#' Perform PCA on Consensus Genes and Generate Plots
#'
#' This function performs Principal Component Analysis (PCA) on a given dataframe of gene expression data, specifically focusing on a subset of consensus genes. It generates two plots: one for PCA scores (PC1 vs. PC2) and another showing the loadings of the top 15 genes for the first two principal components.
#'
#' @param df A data frame containing gene expression data, where rows represent samples and columns represent genes. The column names must match the gene identifiers used in the `cons_genes` vector.
#' @param cons_genes A character vector of gene names to be included in the PCA analysis. These are the consensus genes used for the PCA computation.
#' @param class A factor or vector representing the class labels for the samples. This is used to color the points in the PCA plot and differentiate the sample groups.
#'
#' @return A combined plot containing two ggplot objects: one showing the PCA scores for PC1 and PC2, and another displaying the top 15 genes with the highest loadings for both principal components.
#'
#' @details The function first filters the gene expression data to include only the consensus genes provided.
#' It then performs PCA on this subset of genes, extracting eigenvalues, scores, and loadings.
#' The PCA scores are plotted to show how samples group based on their principal component values, colored by their class labels.
#' The top 15 genes with the highest loadings on PC1 and PC2 are also plotted as vectors in a second plot.
#' The function utilizes the `ggplot2` and `cli` packages to generate the plots.
#'
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming `df` is a data frame of gene expression,
#' # `cons_genes` is a vector of consensus gene names,
#' # and `class_labels` is a factor of class labels for each sample.
#' pca_plot <- pca.consensus(
#'     df = gene_expression_data,
#'     cons_genes = consensus_gene_list, class = class_labels
#' )
#' print(pca_plot)
#' }
#'
#' @import patchwork
#'
#' @export

pca_plot <- function(pca_scores, pca_top_loadings, labels) {

  # create the plot
  pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = factor(labels))) +
    scale_color_manual(values = c("skyblue3", "indianred")) +
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

  # Assign colors based on quadrants (4 colors)
  pca_top_loadings$Col <- ifelse(pca_top_loadings$PC1 > 0 & pca_top_loadings$PC2 > 0, "#66C2A5",
                                 ifelse(pca_top_loadings$PC1 < 0 & pca_top_loadings$PC2 < 0, "#E5C494",
                                        ifelse(pca_top_loadings$PC1 > 0 & pca_top_loadings$PC2 < 0, "#8DA0CB",
                                               "#FC8D62")))

  # Plot of top loadings for each of PC1 and PC2
  loadings_plot <- ggplot(data = pca_top_loadings) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.08, "in")),
                 linewidth = 0.8,
                 colour = pca_top_loadings$Col
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

  p <- (pca_plot | loadings_plot)

  return(p)
}

# Function for creating AUROC Plot
auroc_fc_plot <- function(results) {

  p <- ggplot(results, aes(x = reorder(gene, auroc), y = auroc)) +
    geom_hline(yintercept = 0.5, linetype = "solid", color = "#7e7e7e", linewidth = 0.4) +
    geom_errorbar(aes(ymin = auroc_lower, ymax = auroc_upper, color = FC), width = 0, linewidth = 0.4, position = position_dodge(0.5)) +
    geom_point(aes(color = FC), size = 5, position = position_dodge(0.5)) +
    scale_color_gradient2(low = "#187498", mid = "gray90", high = "#C62E2E", midpoint = 0, limits = c(-2.1, 2.1)) +
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

}

#####################################################################################################################

#' Heatmap Calculation with Hierarchical Clustering
#'
#' @description
#' Computes a hierarchical clustering heatmap matrix for a set of consensus genes in a gene expression data frame.

heatmap_plot_heatmaply <- function(heatmap_data,
                                   colors = grDevices::colorRampPalette(c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))(64),
                                   dendrogram = "both",
                                   show_dendrogram = c(FALSE, TRUE),
                                   scale = "row",
                                   custom_colors = c("1" = "#2e1457", "0" = "#66a182"),
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

mlp_model_plot <- function(importances, test_performance){

  # Create the plot
  cli::cli_alert_info("Creating plot ...")
  p <- ggplot(
    importances,
    aes(x = stats::reorder(Variable, Importance), y = Importance, label = Variable)
  ) +
    geom_segment(aes(y = 0, x = Variable, yend = Importance, xend = Variable), linewidth = 0.8, color = "gray") +
    geom_point(
      stat = "identity",
      size = 4.5,
      color = ifelse(importances$Importance < 0, "#995caf", "#e8ab1b")
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
    coord_flip()
  cli::cli_alert_success("MLP Plot created successfully!")

  return(p)

}

