# Helper function to pick black or white text for each heatmap tile, whichever
# has the higher WCAG contrast with the tile color. The tile color is
# interpolated from colors over limits, as the fill gradient does; values
# outside limits or missing get black text
label_colors <- function(values, colors, limits) {
  pos <- (values - limits[1]) / diff(limits)
  out <- is.na(pos) | pos < 0 | pos > 1
  pos[out] <- 0
  rgb <- grDevices::colorRamp(colors, space = "Lab")(pos) / 255
  # Relative luminance (WCAG 2)
  linear <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  luminance <- drop(linear %*% c(0.2126, 0.7152, 0.0722))
  contrast_black <- (luminance + 0.05) / 0.05
  contrast_white <- 1.05 / (luminance + 0.05)
  ifelse(out | contrast_black >= contrast_white, "black", "white")
}

#' @title Plot Heatmap of Model Performance Metrics
#'
#' @description
#' Creates a heatmap using \code{ggplot2} to visualize the performance of different
#' clustering and dimensionality reduction models across multiple evaluation metrics,
#' for each gene combination. Each combination is labeled with a concise and readable
#' gene summary on the y-axis, prefixed by a combination number (C1, C2, ...).
#'
#' @param top_results The \code{topCombinations} data frame returned by \code{\link{valConsensus}},
#' containing the performance metrics of the clustering methods for each gene combination.
#' @param title A character string specifying the plot title. Default is "Performance metrics of top combinations".
#' @param colors A character vector of at least two colors, from the lowest to the highest values,
#' used as the heatmap gradient. Default: \code{c("#b7efc5", "#10451d")}.
#' @param limits Numeric vector of length 2, specifying the minimum and maximum values for the fill scale. Default: c(0, 1).
#' @param low_color,high_color Deprecated, use \code{colors} instead.
#' If any of them is given, the gradient goes from \code{low_color} to \code{high_color}, and the missing one is
#' taken from the ends of \code{colors}.
#'
#' @details
#' Values are printed in black or white, whichever contrasts more with the tile color.
#'
#' @return A \code{ggplot2} object representing the heatmap.
#'
#' @examples
#' \dontrun{
#' vc <- valConsensus(df.count = count_table, gene.list = gene_list,
#'                    class = class, N = 10, metric = "FScore")
#' p <- plotTopMetrics(
#'   vc$topCombinations,
#'   title = "Performance metrics of top combinations",
#'   colors = c("white", "steelblue"),
#'   limits = c(0, 1)
#' )
#' print(p)
#' }
#'
#' @export
#' @import ggplot2

plotTopMetrics <- function(top_results,
                           title = "Performance metrics of top combinations",
                           colors = c("#b7efc5", "#10451d"),
                           limits = c(0, 1),
                           low_color = NULL,
                           high_color = NULL) {

  if (!"Gene_Combination" %in% names(top_results)) {
    cli::cli_abort("Input data frame must contain a 'Gene_Combination' column.")
  }

  # low_color and high_color were replaced by colors; keep old calls working
  if (!is.null(low_color) || !is.null(high_color)) {
    cli::cli_warn(c(
      "{.arg low_color} and {.arg high_color} are deprecated.",
      "i" = "Use {.arg colors} instead, e.g. {.code colors = c(\"white\", \"steelblue\")}."
    ))
    colors <- c(
      if (is.null(low_color)) colors[1] else low_color,
      if (is.null(high_color)) colors[length(colors)] else high_color
    )
  }
  if (!is.character(colors) || length(colors) < 2) {
    cli::cli_abort("{.arg colors} must be a character vector of at least two colors.")
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
  long_results$Text_Color <- label_colors(long_results$Value, colors, limits)

  # Define order of gene combinations for plotting
  level_order <- rev(unique(gene_mapping$Pretty_Label))

  # Create heatmap with ordered axis
  heatmap <- ggplot(long_results, aes(x = factor(Pretty_Label, levels = level_order), y = Model, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradientn(colors = colors, limits = limits) +
    facet_wrap(~Metric) +
    geom_text(aes(label = round(Value, 2), color = Text_Color), size = 3) +
    scale_color_identity() +
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

#' @title Plot Genes in Top Combinations
#'
#' @description
#' Creates a binary heatmap using \code{ggplot2} showing which genes (rows) are part of
#' each top gene combination (columns). A last column shows, for each gene, the number
#' of combinations containing it out of the total (e.g. "9/10").
#'
#' @param gene_freq The \code{geneFrequency} data frame returned by \code{\link{valConsensus}}: a
#' \code{Gene} column and one binary column per combination (\code{C1}, \code{C2}, ...). Other
#' columns, such as \code{Count} and \code{Combinations}, are ignored. Genes are plotted in the
#' order of the rows.
#' @param title A character string specifying the plot title. Default is "Gene frequencies in top combinations".
#' @param present_color Color of the tiles where the gene is in the combination. Default: "#1b4965".
#' @param absent_color Color of the tiles where the gene is not in the combination. Default: "#e1e5f2".
#'
#' @return A \code{ggplot2} object representing the heatmap.
#'
#' @examples
#' \dontrun{
#' vc <- valConsensus(df.count = count_table, gene.list = gene_list,
#'                    class = class, N = 10, metric = "FScore")
#' p <- plotGeneCombinations(vc$geneFrequency)
#' print(p)
#' }
#'
#' @export
#' @import ggplot2

plotGeneCombinations <- function(gene_freq,
                                 title = "Gene frequencies in top combinations",
                                 present_color = "#1b4965",
                                 absent_color = "#e1e5f2") {

  combo_cols <- grep("^C[0-9]+$", names(gene_freq), value = TRUE)
  if (!is.data.frame(gene_freq) || !"Gene" %in% names(gene_freq) || length(combo_cols) == 0) {
    cli::cli_abort(
      "Input must be a data frame with a {.field Gene} column and combination columns {.field C1}, {.field C2}, ..."
    )
  }

  # Binary gene x combination matrix, in the row and column order of the input
  gene_combinations <- as.matrix(gene_freq[combo_cols])
  rownames(gene_combinations) <- gene_freq$Gene

  # Convert matrix to long format, keeping the row and column order
  n_combos <- ncol(gene_combinations)
  long_results <- data.frame(
    Gene = rep(rownames(gene_combinations), times = n_combos),
    Position = rep(seq_len(n_combos), each = nrow(gene_combinations)),
    Present = factor(as.vector(gene_combinations), levels = c(0, 1), labels = c("No", "Yes"))
  )

  # Number of combinations containing each gene, shown as an extra column
  counts <- data.frame(
    Gene = rownames(gene_combinations),
    Label = paste0(rowSums(gene_combinations), "/", n_combos)
  )

  # Combinations sit at x = 1..n_combos, so the count column goes at n_combos + 1
  ggplot(long_results, aes(x = Position, y = factor(Gene, levels = rev(rownames(gene_combinations))))) +
    geom_tile(aes(fill = Present), color = "white", linewidth = 0.5) +
    geom_text(data = counts, aes(x = n_combos + 1, label = Label), size = 3) +
    scale_x_continuous(
      breaks = seq_len(n_combos + 1),
      labels = c(colnames(gene_combinations), "Count"),
      # Room on the right so labels like "10/10" are not clipped
      expand = expansion(add = c(0.1, 0.5))
    ) +
    scale_fill_manual(values = c(No = absent_color, Yes = present_color), drop = FALSE) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = title, fill = "In combination") +
    theme(
      axis.text.x = element_text(size = 9, color = "#7c7b7b"),
      axis.text.y = element_text(size = 9, color = "black"),
      panel.grid = element_blank()
    )
}
