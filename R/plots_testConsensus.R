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
#' @details The function first filters the gene expression data to include only the consensus genes provided. It then performs PCA on this subset of genes, extracting eigenvalues, scores, and loadings. The PCA scores are plotted to show how samples group based on their principal component values, colored by their class labels. The top 15 genes with the highest loadings on PC1 and PC2 are also plotted as vectors in a second plot. The function utilizes the `ggplot2`, `cli`, and `reshape2` packages to generate the plots.
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
#' @export
pca.consensus <- function(df, cons_genes, class) {
    # Create a matrix from our table of counts
    # rows: samples, cols: genes
    pca_matrix <- as.matrix(df[, colnames(df) %in% cons_genes])

    # Perform the PCA
    cli::cli_alert_info("Generating PCA plots for consensus genes...")
    sample_pca <- stats::prcomp(pca_matrix, scale = TRUE)
    # extract the eigenvalues
    pc_eigenvalues <- sample_pca$sdev^2

    # create a data frame with the eigenvalues
    pc_eigenvalues <- data.frame(
        PC = factor(1:length(pc_eigenvalues)),
        variance = pc_eigenvalues
    )
    # add a new column with the percent variance
    pc_eigenvalues$pct <- pc_eigenvalues$variance / sum(pc_eigenvalues$variance) * 100
    # add another column with the cumulative variance explained
    pc_eigenvalues$pct_cum <- cumsum(pc_eigenvalues$pct)

    # extract the PC scores
    pc_scores <- sample_pca$x
    pc_scores <- data.frame(sample = rownames(pc_scores), pc_scores)

    # extract the PC loadings
    pc_loadings <- sample_pca$rotation
    pc_loadings <- data.frame(gene = rownames(pc_loadings), pc_loadings)

    # reshape the loadings to long format
    pc_loadings_long <- reshape2::melt(pc_loadings, id.vars = "gene", measure.vars = c("PC1", "PC2"), variable.name = "PC", value.name = "loading")
    pc_loadings_long <- pc_loadings_long[order(-abs(pc_loadings_long$loading)), ]
    # get the top 15 genes for each PC
    top_genes <- unique(pc_loadings_long$gene[1:15]) # TODO: make this a parameter
    top_loadings <- pc_loadings[pc_loadings$gene %in% top_genes, ]

    # get the PC scores from prcomp object
    pc_scores_p <- data.frame(sample = rownames(sample_pca$x), sample_pca$x)
    # create the plot
    pca_plot <- ggplot(pc_scores_p, aes(x = PC1, y = PC2, colour = factor(class))) +
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


    # Create Color Column named 'Col'
    top_loadings$Col <- NA
    top_loadings[top_loadings$PC1 > 0 & top_loadings$PC2 > 0, "Col"] <- "#66C2A5"
    top_loadings[top_loadings$PC1 < 0 & top_loadings$PC2 < 0, "Col"] <- "#E5C494"
    top_loadings[top_loadings$PC1 > 0 & top_loadings$PC2 < 0, "Col"] <- "#8DA0CB"
    top_loadings[top_loadings$PC1 < 0 & top_loadings$PC2 > 0, "Col"] <- "#FC8D62"

    # plot dei top_loadings per ognuna delle PC1 e PC2
    loadings_plot <- ggplot(data = top_loadings) +
        geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
            arrow = arrow(length = unit(0.08, "in")),
            linewidth = 0.8,
            colour = top_loadings$Col
        ) +
        ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.001, size = 2.5
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
    cli::cli_alert_success("Plots created!")

    return(p)
}

#' Compute AUROC and Fold Change for Consensus Genes
#'
#' This function computes the Area Under the Receiver Operating Characteristic curve (AUROC)
#' and the Fold Change (FC) for each gene in a set of consensus genes. It returns a plot of
#' AUROC values for each gene, with error bars representing the confidence interval for each AUROC
#' and colors indicating the Fold Change values.
#'
#' @param df A data frame containing gene expression data with rows as samples and columns as genes.
#' @param cons_genes A character vector of selected consensus genes to be analyzed.
#' @param class A vector of class labels (0 or 1) corresponding to the samples in `df`.
#'
#' @return A ggplot object displaying the AUROC values for each gene, with fold change color-coding
#' and error bars for the confidence interval of AUROC.
#'
#' @details
#' The function uses the \code{\link[pROC]{roc}} function to compute the AUROC for each gene in the
#' consensus set and calculates the Fold Change manually as the difference in the mean expression
#' between the two classes (class 1 and class 0).
#'
#' The output plot orders genes based on their mean AUROC and color-codes the points according to
#' their corresponding Fold Change. The error bars represent the 95% confidence intervals for the
#' AUROC values.
#'
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming 'expression_data' is a data frame with
#' # gene expression data and 'labels' is a vector of class labels
#' # cons_genes <- c("Gene1", "Gene2", "Gene3") # Example of consensus genes
#' # result_plot <- auroc.FC.consensus(expression_data, cons_genes, labels)
#' # print(result_plot)
#' }
#'
#' @export
auroc.FC.consensus <- function(df, cons_genes, class) {
    df <- df[, colnames(df) %in% cons_genes]
    results <- tibble::tibble(
        gene = character(),
        auroc = numeric(),
        auroc_upper = numeric(),
        auroc_lower = numeric(),
        FC = numeric()
    )

    cli::cli_alert_info("Computing AUROC and Fold Change ...")
    for (gene in colnames(df)) {
        predictor <- df[[gene]]
        auc_obj <- pROC::roc(class, predictor, auc = TRUE, ci = TRUE, direction = "<", levels = (c(0, 1)))

        # Calculate Fold Change manually
        pred_with_class <- data.frame(predictor, class)
        FC <- mean(pred_with_class[pred_with_class$class == 1, "predictor"]) -
            mean(pred_with_class[pred_with_class$class == 0, "predictor"])

        results <- results %>%
            dplyr::add_row(
                gene = gene,
                auroc = auc_obj$auc[1],
                auroc_upper = auc_obj$ci[1],
                auroc_lower = auc_obj$ci[3],
                FC = FC
            )
    }
    cli::cli_alert_success("AUROC and Fold Change computed successfully!")

    # Auroc Plot
    cli::cli_alert_info("Creating AUROC plot ...")
    auroc <- NULL
    mean_auroc <- NULL
    auroc_lower <- NULL
    auroc_upper <- NULL
    p <- results %>%
        group_by(gene) %>%
        dplyr::mutate(mean_auroc = mean(auroc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gene = forcats::fct_reorder(gene, mean_auroc)) %>% # Order by AUC
        ggplot(aes(x = gene, y = auroc)) +
        geom_hline(yintercept = 0.5, linetype = "solid", color = "#7e7e7e", linewidth = 0.4) +
        geom_errorbar(aes(ymin = auroc_lower, ymax = auroc_upper, color = FC), width = 0, linewidth = 0.4, position = position_dodge(0.5)) +
        geom_point(aes(color = FC), size = 5, position = position_dodge(0.5)) +
        scale_color_gradient2(low = "#187498", mid = "gray90", high = "#C62E2E", midpoint = 0, limits = c(-2.1, 2.1)) +
        # geom_line(aes(y = mean_auc, group = 1), color = "green", linetype = "solid", linewidth = 0.2) +
        scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        scale_shape_manual(values = c(16, 17, 15, 3, 4, 18, 8)) +
        labs(x = "", y = "AUROC") +
        # coord_flip() +
        theme_minimal() +
        theme(
            # legend.position = "right",
            # legend.justification = "center",
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
    cli::cli_alert_success("AUROC plot created successfully!")

    return(p)
}

#' Create a Plot for MLP Model Variable Importance
#'
#' This function tunes and evaluates a multilayer perceptron (MLP) model using a dataset of consensus genes and creates a plot to visualize the variable importance. The plot shows the scaled importance of each variable in the model, and includes key performance metrics (accuracy, AUROC, and Brier score).
#'
#' @param df A data frame containing the gene expression data. The dataframe must include a column for the class variable (`class`), which is used for stratified sampling.
#' @param cons_genes A character vector containing the names of consensus genes to be used for model training.
#' @param class A factor or numeric vector representing the class labels. This is the outcome variable used for classification.
#'
#' @return A `ggplot` object showing the scaled variable importance for the MLP model, with an annotation of model performance metrics.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' plot <- mlp.model.plot(df = gene_data, cons_genes = consensus_genes, class = class_labels)
#' print(plot)
#' }
#'
#' @export
mlp.model.plot <- function(df, cons_genes, class) {
    # future::plan(future::multisession, workers = parallel::detectCores() - 1)

    withr::with_seed(
        seed = 123,
        code = {
            # Keep from df only genes of consensus
            df <- df[, colnames(df) %in% cons_genes]
            cli::cli_alert_info("Setting up MLP on {length(cons_genes)} genes ...")

            # Put class column at the end
            df$class <- class

            # Create a new train and test dataset
            cli::cli_alert_info("Splitting data into train and test ...")
            train_test_split <- rsample::initial_split(df, prop = 0.7, strata = class)
            # train_test_split <- split.train.test(df, prop = 0.7, seed = 456)
            train <- rsample::training(train_test_split)
            test <- rsample::testing(train_test_split)
            train_resamples <- rsample::vfold_cv(train, v = 3, strata = class)
            cli::cli_alert_success("Data splitted successfully!")


            # define the recipe
            cli::cli_alert_info("Setting up the recipe ...")
            df_recipe <-
                # which consists of the formula (outcome ~ predictors)
                recipes::recipe(class ~ ., data = df) %>%
                # and some pre-processing steps
                recipes::step_zv(recipes::all_predictors()) %>%
                recipes::step_normalize(recipes::all_numeric()) %>%
                recipes::step_corr(recipes::all_predictors())
            cli::cli_alert_success("Recipe set up successfully!")

            # apply the recipe to the training data
            cli::cli_alert_info("Applying the recipe to the training data ...")
            df_train_preprocessed <- df_recipe %>%
                # apply recipe to training data
                recipes::prep(train) %>%
                # extract pre-processed data
                recipes::juice()
            cli::cli_alert_success("Recipe applied successfully!")

            ## Setup the models
            # mlp model
            cli::cli_alert_info("Setting up MLP model ...")
            mlp_model <-
                # specify that the model is a random forest
                parsnip::mlp() %>%
                # specify parameters that need to be tuned
                parsnip::set_args(
                    hidden_units = tune(),
                    penalty = tune(),
                    epochs = tune()
                ) %>%
                # select the engine/package that underlies the model
                parsnip::set_engine("nnet") %>%
                # choose either the continuous regression or binary classification mode
                parsnip::set_mode("classification")

            # set the workflow
            mlp_workflow <- workflows::workflow() %>%
                # add the recipe
                workflows::add_recipe(df_recipe) %>%
                # add the model
                workflows::add_model(mlp_model)
            cli::cli_alert_success("MLP model set up successfully!")

            # extract results
            cli::cli_alert_info("Tuning MLP model ...")
            mlp_tune_results <- mlp_workflow %>%
                tune::tune_grid(
                    resamples = train_resamples, # CV object
                    grid = 10
                )
            cli::cli_alert_success("MLP model tuned successfully!")

            # print results
            mlp_tune_results %>%
                workflowsets::collect_metrics()

            # finalize the workflow
            cli::cli_alert_info("Finalizing MLP model ...")
            param_final <- mlp_tune_results %>%
                tune::select_best(metric = "roc_auc")

            mlp_workflow <- mlp_workflow %>%
                tune::finalize_workflow(param_final)
            cli::cli_alert_success("MLP model finalized successfully!")

            # Evaluate on test set
            cli::cli_alert_info("Evaluating MLP model on test set ...")
            mlp_fit <- mlp_workflow %>%
                # fit on the training set and evaluate on test set
                tune::last_fit(train_test_split)

            # Inspect performances on test set
            test_performance <- mlp_fit %>% workflowsets::collect_metrics()

            # Fitting and using final model
            final_model <- parsnip::fit(mlp_workflow, df)
            cli::cli_alert_success("MLP model evaluated successfully!")

            # Assess variable importance
            cli::cli_alert_info("Assessing variable importance ...")
            mlp_obj <- workflowsets::extract_fit_parsnip(final_model)

            importances <- as.data.frame(vip::vi(mlp_obj))
            cli::cli_alert_success("Variable importance assessed successfully!")
        }
    )

    # future::plan(future::sequential)

    # rescale Importance column in importances df in range -1 and 1
    importances$Importance <- scales::rescale(importances$Importance, to = c(-1, 1))
    importances$type <- ifelse(importances$Importance < 0, 0.026, -0.026)
    importances$type_hjust <- ifelse(importances$Importance < 0, 0, 1)
    # importances$score_hjust <- ifelse(importances$Importance < 0, 2, -1)

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
            color = if_else(importances$Importance < 0, "#995caf", "#e8ab1b")
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
            ), hjust = 0, vjust = 0.5, size = 5, color = "#7e7e7e",
            bg = "white",
            box.padding = unit(0.5, "lines")
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
    cli::cli_alert_success("Plot created successfully!")

    return(p)
}

#' Create a Heatmap with Hierarchical Clustering
#'
#' This function generates a heatmap from gene expression data, with hierarchical clustering of both rows and columns. The heatmap can be customized with various parameters, such as dendrogram visibility, scaling method, color schemes, and additional annotations.
#'
#' @param df A data frame containing gene expression data. The dataframe must include the genes of interest for the heatmap (specified in `cons_genes`).
#' @param cons_genes A character vector containing the names of genes to be used in the heatmap.
#' @param class A factor or numeric vector representing the class labels for the samples. This is used to assign colors to the side of the heatmap.
#' @param dendrogram A character string specifying which dendrograms to display. Possible values are "both" (both row and column), "row" (only row dendrogram), and "column" (only column dendrogram).
#' @param show_dendrogram A logical vector specifying whether to show the dendrogram for rows and/or columns.
#' @param scale A character string specifying how to scale the data for the heatmap. Options are "row", "column", or "none".
#' @param custom_colors A named vector specifying custom colors for the class labels. The names should match the values in `class`.
#' @param margins A numeric vector specifying the margins for the heatmap plot. The format is `c(bottom, left, top, right)`.
#' @param grid_color A color for the grid lines of the heatmap.
#' @param grid_width The width of the grid lines.
#' @param branches_lwd The line width for the branches in the dendrogram.
#' @param fontsize_row The font size for row labels.
#' @param fontsize_col The font size for column labels.
#' @param colors A vector of colors for the heatmap cells. Defaults to `viridis::viridis(n = 64)` if not specified.
#' @param hclust.method A character string specifying the hierarchical clustering method to use. Default is "complete".
#' @param distance.method A character string specifying the distance method for clustering. Default is "euclidean".
#' @param k_row The number of clusters for rows (optional).
#' @param k_col The number of clusters for columns (optional).
#' @param Rowv A logical value specifying whether to reorder rows based on hierarchical clustering. Default is `TRUE`.
#' @param Colv A logical value specifying whether to reorder columns based on hierarchical clustering. Default is `TRUE`.
#' @param scale_fill_gradient_fun A function to apply a custom gradient scale for the heatmap. Default is `NULL`.
#' @param limits A numeric vector specifying the color limits for the heatmap. Default is `NULL`.
#' @param na.value A color to use for missing values in the heatmap. Default is "grey50".
#' @param cellnote A matrix or data frame providing cell annotations. Default is `NULL`.
#' @param cellnote_size The font size for the cell annotations.
#' @param cellnote_textposition A string specifying the position of the text within the cell. Default is "middle center".
#' @param key.title A string specifying the title for the color key.
#' @param key.xlab A string specifying the x-axis label for the color key.
#' @param key.ylab A string specifying the y-axis label for the color key.
#'
#' @return A `ggplot` object representing the heatmap with hierarchical clustering.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' heatmap_plot <- heatmap.plot(df = gene_data, cons_genes = consensus_genes, class = class_labels)
#' print(heatmap_plot)
#' }
#'
#' @export
heatmap.plot <- function(
    df, cons_genes,
    class,
    dendrogram = "both",
    show_dendrogram = c(FALSE, TRUE),
    scale = "row",
    custom_colors = c("1" = "#2e1457", "0" = "#66a182"),
    # custom_colors = c("1" = "#9c3f2d", "0" = "#3f7ab9"),
    margins = c(60, 100, 40, 20),
    grid_color = "white",
    grid_width = 0.00001,
    branches_lwd = 0.4,
    fontsize_row = 12,
    fontsize_col = 5,
    colors = NULL,
    hclust.method = "complete",
    distance.method = "euclidean",
    k_row = NULL,
    k_col = NULL,
    Rowv = TRUE,
    Colv = TRUE,
    scale_fill_gradient_fun = NULL,
    limits = NULL,
    na.value = "grey50",
    cellnote = NULL,
    cellnote_size = NULL,
    cellnote_textposition = "middle center",
    key.title = "",
    key.xlab = "Value",
    key.ylab = "Frequency") {
    # Crea i colori laterali per la heatmap
    side_colors_df <- data.frame(class = as.factor(class))
    # Select genes of interest
    df <- df[, colnames(df) %in% cons_genes]
    # Traspose df
    df <- as.data.frame(t(df))


    # Assicurati che siano forniti i colori personalizzati
    if (is.null(colors)) {
        colors <- viridis::viridis(n = 64)
    }

    # Clustering per righe e colonne separatamente
    cli::cli_alert_info("Clustering rows and columns using {hclust.method} ...")
    hc_rows <- stats::hclust(stats::dist(as.matrix(df), method = distance.method), method = hclust.method)
    hc_cols <- stats::hclust(stats::dist(t(df), method = distance.method), method = hclust.method)
    cli::cli_alert_success("Rows and columns clustered successfully!")

    # Crea la heatmap usando heatmaply
    cli::cli_alert_info("Creating heatmap ...")
    p <- heatmaply::heatmaply(
        df,
        Rowv = stats::as.dendrogram(hc_rows),
        Colv = stats::as.dendrogram(hc_cols),
        dendrogram = dendrogram,
        show_dendrogram = show_dendrogram,
        xlab = "",
        ylab = "",
        main = "",
        scale = scale,
        margins = margins,
        grid_color = grid_color,
        grid_width = grid_width,
        titleX = FALSE,
        hide_colorbar = FALSE,
        branches_lwd = branches_lwd,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        showticklabels = c(FALSE, TRUE),
        labRow = rownames(df),
        plot_method = "ggplot",
        colors = colors,
        col_side_colors = side_colors_df,
        col_side_palette = custom_colors,
        scale_fill_gradient_fun = scale_fill_gradient_fun,
        limits = limits,
        na.value = na.value,
        cellnote = cellnote,
        cellnote_size = cellnote_size,
        cellnote_textposition = cellnote_textposition,
        key.title = key.title,
        key.xlab = key.xlab,
        key.ylab = key.ylab
    )
    cli::cli_alert_success("Heatmap created successfully!")

    # Stampa il plot
    return(p)
}
