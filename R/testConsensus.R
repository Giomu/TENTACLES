# Helper function to testConsensus
# This function computes the AUROC and Fold Change for the consensus genes
#' @export
compute.auroc.FC <- function(df, cons_genes, class) {
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
    p <- results %>%
        group_by(gene) %>%
        dplyr::mutate(mean_auroc = mean(auroc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gene = forcats::fct_reorder(gene, mean_auroc)) %>% # Order by AUC
        ggplot(aes(x = gene, y = auroc)) +
        geom_hline(yintercept = 0.5, linetype = "solid", color = "gray", linewidth = 0.2) +
        geom_errorbar(aes(ymin = auroc_lower, ymax = auroc_upper, color = FC), width = 0, linewidth = 0.4, position = position_dodge(0.5)) +
        geom_point(aes(color = FC), size = 5, position = position_dodge(0.5)) +
        scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0, limits = c(-2.6, 2.6)) +
        # geom_line(aes(y = mean_auc, group = 1), color = "green", linetype = "solid", linewidth = 0.2) +
        scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        scale_shape_manual(values = c(16, 17, 15, 3, 4, 18, 8)) +
        labs(x = "", y = "AUROC") +
        # coord_flip() +
        theme_minimal() +
        theme(
            # legend.position = "right",
            # legend.justification = "center",
            axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(hjust = 1, size = 19),
            axis.title.x = element_text(size = 18),
            legend.key.size = unit(0.5, "cm"),
            legend.title = element_text(size = 15, vjust = 0.5, hjust = 0),
            legend.text = element_text(size = 15, hjust = 0),
            panel.grid = element_blank(),
            axis.ticks = element_line(color = "black", size = 0.3),
            axis.line.x = element_line(color = "black", linewidth = 0.3),
            axis.line.y = element_line(color = "black", linewidth = 0.3)
        )
    cli::cli_alert_success("AUROC plot created successfully!")

    return(p)
}

# TODO fix the legend
# Helper function to testConsensus
# This function creates a heatmap plot. It is simply a wrapper around
# the heatmaply package
heatmap.plot <- function(
    df, cons_genes,
    class,
    dendrogram = "both",
    show_dendrogram = c(FALSE, TRUE),
    scale = "row",
    custom_colors = c("1" = "#440154FF", "0" = "#FDE725FF"),
    margins = c(60, 100, 40, 20),
    grid_color = "white",
    grid_width = 0.00001,
    branches_lwd = 0.4,
    fontsize_row = 5,
    fontsize_col = 7,
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
    key.title = "Legend",
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
        colors <- viridis::mako(n = 256)
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
        label_names = c("row", "column", "value"),
        labCol = colnames(df),
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


# Helper function to testConsensus
# This function performs the validation of the consensus genes using a MLP model
validation.model <- function(df, cons_genes, class) {
    #future::plan(future::multisession, workers = parallel::detectCores() - 1)

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
            #train_test_split <- split.train.test(df, prop = 0.7, seed = 456)
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
                    epochs = tune()) %>%
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

    #future::plan(future::sequential)

    # rescale Importance column in importances df in range -1 and 1
    importances$Importance <- scales::rescale(importances$Importance, to = c(-1, 1))
    importances$type <- ifelse(importances$Importance < 0, 0.09, -0.09)
    importances$type_score <- ifelse(importances$Importance < 0, 1.5, -0.6)

    # Create the plot
    cli::cli_alert_info("Creating plot ...")
    p <- ggplot(
        importances,
        aes(x = reorder(Variable, Importance), y = Importance, label = Variable)
    ) +
        geom_point(stat = "identity", aes(fill = Importance), size = 5, shape = 21, stroke = NA) +
        scale_fill_gradient2(low = "#ffee8c", high = "#218882") +
        geom_segment(aes(y = 0, x = Variable, yend = Importance, xend = Variable), color = "black") +
        geom_text(aes(label = Variable, y = type), color = "black", size = 2.5) +
        geom_text(aes(label = round(Importance, 2), hjust = type_score), color = "darkgrey", size = 2) +
        labs(
            title = "MLP Variables Importances",
            subtitle = paste0(
                "Accuracy: ", round(test_performance$.estimate[[1]], 3), "\n",
                "Roc Auc: ", round(test_performance$.estimate[[2]], 3), "\n",
                "Brier class: ", round(test_performance$.estimate[[3]], 3)
            ),
            x = "", y = ""
        ) +
        guides(y = "none") +
        coord_flip() +
        theme_minimal()
    cli::cli_alert_success("Plot created successfully!")

    return(p)
}

# Helper function to testConsensus
# This function performs PCA on the consensus genes
pca.validation <- function(df, cons_genes, class) {
    # Create a matrix from our table of counts
    # rows: samples, cols: genes
    pca_matrix <- df[, colnames(df) %in% cons_genes] %>%
        as.matrix()

    # Perform the PCA
    cli::cli_alert_info("Performing PCA on consensus genes ...")
    sample_pca <- stats::prcomp(pca_matrix, scale = TRUE)
    cli::cli_alert_success("PCA done!")

    # as_tibble(pca_matrix)
    # as_tibble(pca_matrix, rownames = "Sample_ID")


    # tiro fuori gli autovalori
    cli::cli_alert_info("Extracting eigenvalues ...")
    pc_eigenvalues <- sample_pca$sdev^2

    # create a "tibble" manually with
    # a variable indicating the PC number
    # and a variable with the variances
    pc_eigenvalues <- tibble::tibble(
        PC = factor(1:length(pc_eigenvalues)),
        variance = pc_eigenvalues
    ) %>%
        # add a new column with the percent variance
        dplyr::mutate(pct = variance / sum(variance) * 100) %>%
        # add another column with the cumulative variance explained
        dplyr::mutate(pct_cum = cumsum(pct))


    pc_scores <- sample_pca$x
    pc_scores <- pc_scores %>%
        # convert to a tibble retaining the sample names as a new column
        tibble::as_tibble(rownames = "sample") # TODO review tibble use
    cli::cli_alert_success("Eigenvalues extracted!")

    # tiro fuori autivettori
    cli::cli_alert_info("Extracting eigenvectors ...")
    pc_loadings <- sample_pca$rotation
    pc_loadings <- pc_loadings %>%
        tibble::as_tibble(rownames = "gene")
    cli::cli_alert_success("Eigenvectors extracted!")

    cli::cli_alert_info("Extracting top genes ...")
    top_genes <- pc_loadings %>%
        # select only the PCs we are interested in
        dplyr::select(gene, PC1, PC2) %>%
        # convert to a "long" format
        tidyr::pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>%
        # for each PC
        dplyr::group_by(PC) %>%
        # arrange by descending order of loading
        dplyr::arrange(desc(abs(loading))) %>%
        # take the 10 top rows
        dplyr::slice(1:15) %>% # TODO parameter asking the number of genes
        # pull the gene column as a vector
        dplyr::pull(gene) %>%
        # ensure only unique genes are retained
        unique()

    top_loadings <- pc_loadings %>%
        dplyr::filter(gene %in% top_genes)
    cli::cli_alert_success("Top genes extracted!")

    cli::cli_alert_info("Creating plots ...")
    # get the PC scores from prcomp object
    pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
        # convert to a tibble retaining the sample names as a new column
        tibble::as_tibble(rownames = "sample") %>%
        # create the plot
        ggplot(aes(x = PC1, y = PC2, colour = factor(class))) +
        scale_color_manual(values = c("lightblue", "indianred")) +
        geom_point(size = 3) +
        # suppressMessages(ggside::geom_xsidedensity(aes(color = class), linewidth = 0.65, show.legend = FALSE)) +
        stat_ellipse(linewidth = 0.7, linetype = 2, type = "norm") +
        stat_ellipse(type = "t") +
        theme_minimal() +
        coord_fixed(ratio = 1) +
        labs(title = "PC scores", colour = "Class") +
        theme(
            legend.box.just = "right", legend.position = "bottom",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 10)
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
            arrow = arrow(length = unit(0.05, "in")),
            colour = top_loadings$Col
        ) +
        geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 2
        ) +
        scale_x_continuous(expand = c(0.02, 0.02)) +
        theme_minimal() +
        coord_fixed(ratio = 1) +
        labs(title = "PC loadings", x = "", y = "")

    # ...e li metto insieme!
    p <- (pca_plot | loadings_plot) + patchwork::plot_annotation(tag_levels = "A")
    cli::cli_alert_success("Plots created!")

    return(p)
}




#' @title Test Consensus
#' @description This function tests a list of consensus genes on a df.count performing PCA,
#' AUROC, Heatmap, and tuning and fitting a MLP Model.
#'
#' @param df.count A data frame containing the counts data. Gene names on columns and samples on rows.
#' @param gene.list A vector list containing consensus genes to test.
#' @param class A vector containing the true class labels.
#' @return The PCA plot, AUROC plot, Heatmap plot, and Validation MLP Model plot
#'
#' @details
#' This function tests the consensus genes by performing PCA, AUROC, Heatmap,
#' and Validation MLP Model on the provided dataset. PCA is performed on the consensus genes,
#' AUROC is computed for each gene, Heatmap is created using the consensus genes,
#' and a MLP model is tuned and fitted on the filtered dataset.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom cli cli_alert_info cli_alert_success cli_abort cli_h1 cli_h2
#' @importFrom stats prcomp as.dendrogram hclust dist
#' @importFrom tibble tibble as_tibble add_row
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom viridis mako
#' @importFrom vip vi
#' @importFrom workflowsets collect_metrics extract_fit_parsnip
#' @importFrom workflows add_model add_recipe workflow
#' @importFrom parsnip fit mlp set_args set_engine set_mode
#' @importFrom patchwork plot_annotation
#' @importFrom scales rescale
#' @importFrom heatmaply heatmaply
#' @importFrom pROC roc
#' @importFrom sva ComBat
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom future plan multisession sequential
#' @importFrom rsample vfold_cv training testing
#' @importFrom recipes recipe step_zv step_normalize step_corr prep juice all_predictors all_numeric
#' @importFrom withr with_seed
#' @importFrom patchwork plot_annotation
#' @importFrom tune select_best finalize_workflow last_fit tune_grid tune
#'
#' @examples
#' /dontrun{
#' testConsensus(
#'    df.count, cons_genes, class)}
#'
#' @export
testConsensus <- function(df.count, gene.list, class) {

    cons_genes <- gene.list
    cli::cli_h2("Consensus Genes Testing")

    # Run the PCA anlysis
    cli::cli_h3("Principal Component Analysis")
    pca_plot <- pca.validation(df.count, cons_genes, class)
    print(pca_plot)
    cli::cli_alert_success("PCA analysis completed successfully!")

    cli::cli_h3("AUROC Analysis")
    auroc_plot <- compute.auroc.FC(df.count, cons_genes, class)
    print(auroc_plot)
    cli::cli_alert_success("AUROC analysis completed successfully!")

    cli::cli_h3("Heatmap Analysis")
    heatmap <- heatmap.plot(df.count, cons_genes, class)
    print(heatmap)
    cli::cli_alert_success("Heatmap analysis completed successfully!")

    # Run the validation model
    cli::cli_h3("Validation MLP Model")
    model_plot <- validation.model(df.count, cons_genes, class)
    print(model_plot)
    cli::cli_alert_success("Validation MLP Model completed successfully!")
}
