#' Class for the testConsensus Object
#'
#' @description
#' The `testConsensus.obj` class is an S4 class designed to store and organize the results of various analyses (PCA, AUROC, Heatmap, MLP) performed on a set of consensus genes.
#'
#' @slot data A list containing the input data:
#'   - `df.count`: A data frame of gene expression counts (samples x genes).
#'   - `gene.list`: A character vector of consensus genes to be tested.
#'   - `class`: A binary vector of class labels for each sample.
#'
#' @slot pca A list containing the PCA results:
#'   - `scores`: A data frame of PCA scores (samples x components).
#'   - `loadings`: A data frame of PCA loadings (genes x components).
#'   - `explained_variance`: A numeric vector of variance explained by each component.
#'   - `top_loadings`: A character vector of the top N genes with the highest absolute loadings.
#'   - `plot`: A ggplot object of the PCA visualization.
#'
#' @slot auroc_fc A list containing the AUROC and Fold Change results:
#'   - `data`: A data frame with AUROC values and Fold Change for each gene.
#'   - `plot`: A ggplot object of the AUROC/FC visualization.
#'
#' @slot heatmap A list containing the Heatmap and Hierarchical Clustering results:
#'   - `data`: A matrix of normalized gene expression values.
#'   - `plot`: A heatmap plot object.
#'
#' @slot mlp A list containing the Multi-Layer Perceptron (MLP) model results:
#'   - `importances`: A data frame of feature importances.
#'   - `test_performance`: A list with model performance metrics (accuracy, precision, recall).
#'   - `plot`: A ggplot object of the MLP performance visualization.
#'
#' @export
methods::setClass(
  "testConsensus.obj",
  slots = list(
    data = "list",       # Input data (df.count, gene.list, class)
    pca = "list",        # PCA results (scores, loadings, explained variance, plot)
    auroc_fc = "list",   # AUROC and FC results (data, plot)
    heatmap = "list",    # Heatmap with HC results (data, plot)
    mlp = "list"         # MLP model results (importances, test_performance, plot)
  )
)


#' Display Method for testConsensus.obj
#'
#' @description
#' Custom method for displaying the content of a `testConsensus.obj` object.
#' Provides a clear summary of the stored data, including the input data, PCA results, and AUROC/FC analysis.
#'
#' @param object An object of class `testConsensus.obj`.
#'
#' @details
#' The method dynamically adapts to the content of each slot, ensuring a clear display of large data frames (only the first 5 columns and rows) and a compact view of long vectors.
#'
#' @examples
#' \dontrun{
#' # Assuming `consensus_obj` is an existing testConsensus.obj object:
#' consensus_obj
#' }
#'
# Custom 'show' method for the 'testConsensus.obj' class with enhanced display
methods::setMethod(
  "show", "testConsensus.obj",
  function(object) {
    cat("Object of class 'testConsensus.obj'\n")

    # Iterating over each slot dynamically
    for (slot_name in slotNames(object)) {
      cat("\nSlot '", slot_name, "':\n", sep = "")
      slot_content <- slot(object, slot_name)

      if (length(slot_content) == 0) {
        cat("Empty\n")
      } else {
        # Specific treatment for each slot type
        if (is.list(slot_content)) {
          for (sub_name in names(slot_content)) {
            cat("\n- ", sub_name, ":\n", sep = "")
            sub_content <- slot_content[[sub_name]]

            if (is.data.frame(sub_content)) {
              cat(sprintf("* Data frame (%d rows, %d columns)\n",
                          nrow(sub_content), ncol(sub_content)))

              # Displaying only first 5 columns and 5 rows if too large
              if (ncol(sub_content) > 10) {
                cat("  Columns (first 5 shown):", paste(colnames(sub_content)[1:5], collapse = ", "), "...\n")
                print(sub_content[1:5, 1:5, drop = FALSE])
              } else {
                print(head(sub_content, 5))
              }
            } else if (is.vector(sub_content)) {
              # Compact display for large vectors
              if (length(sub_content) > 10) {
                cat("  Values (first 10 shown): ", paste(sub_content[1:10], collapse = ", "),
                    sprintf("... (%d total)\n", length(sub_content)))
              } else {
                cat("  Values: ", paste(sub_content, collapse = ", "), "\n")
              }
            } else if (inherits(sub_content, "ggplot")) {
              cat("* Plot (ggplot object)\n")
            } else {
              print(sub_content)
            }
          }
        } else {
          print(slot_content)
        }
      }
    }
  }
)

### ------------------------ Helper Functions ------------------------ ###

# Principal Component Analysis (PCA) Calculation
#
# @description
# Computes the PCA on a given data frame containing gene expression counts, selecting the top N genes based on absolute loadings.
#
# @param df A data frame where rows represent samples and columns represent genes.
# @param cons_genes A character vector containing the names of the consensus genes to be used in the PCA.
# @param top_n An integer specifying the number of top genes to select based on the absolute loading values. Default is 15.
#
# @return A list with the following components:
#   - `scores`: A data frame with PCA scores for each sample.
#   - `loadings`: A data frame with PCA loadings for each gene.
#   - `explained_variance`: A data frame containing the variance explained by each principal component.
#   - `top_loadings`: A data frame with the top N genes with the highest absolute loadings.
#
# @details
# This function performs PCA using the prcomp function from the stats package. The input data frame is filtered to keep only the specified consensus genes.
# The output includes PCA scores, loadings, explained variance, and the top N genes with the highest loadings.
#
# @examples
# \dontrun{
# pca_results <- pca_calculation(df = expression_data, cons_genes = gene_list, top_n = 20)
# }
#
# --- PCA Calculation Function ---
pca_calculation <- function(df, cons_genes, top_loadings = 15) {

  if (length(cons_genes) < 2) {
    cli::cli_abort("PCA requires at least two genes.")
  }

  # Create a matrix from our table of counts
  pca_matrix <- as.matrix(df[, colnames(df) %in% cons_genes])

  # Perform the PCA
  sample_pca <- stats::prcomp(pca_matrix, scale = TRUE)

  # Calculate eigenvalues and variance
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues_df <- data.frame(
    PC = seq_along(pc_eigenvalues),
    variance = pc_eigenvalues,
    pct = (pc_eigenvalues / sum(pc_eigenvalues)) * 100,
    pct_cum = cumsum((pc_eigenvalues / sum(pc_eigenvalues)) * 100)
  )

  # Extract PCA scores and loadings
  pc_scores <- data.frame(sample = rownames(sample_pca$x), sample_pca$x)
  pc_loadings <- data.frame(gene = rownames(sample_pca$rotation), sample_pca$rotation)

  # Reshape loadings to long format (Base R)
  pc_loadings_long <- data.frame(
    gene = rep(rownames(pc_loadings), 2),
    PC = rep(c("PC1", "PC2"), each = nrow(pc_loadings)),
    loading = c(pc_loadings$PC1, pc_loadings$PC2)
  )

  # Sort by absolute loading values and select top genes
  pc_loadings_long <- pc_loadings_long[order(-abs(pc_loadings_long$loading)), ]

  # Select the top N genes with the highest absolute loadings (long format)
  top_genes <- unique(pc_loadings_long$gene[1:top_loadings])

  # Convert to wide format, keeping only the top N genes
  top_loadings <- pc_loadings[pc_loadings$gene %in% top_genes, c("gene", "PC1", "PC2")]

  results <- list(
    scores = pc_scores,
    loadings = pc_loadings,
    explained_variance = pc_eigenvalues_df,
    top_loadings = top_loadings
  )

  return(results)
}

# AUROC and Fold Change Calculation
#
# @description
# Computes the Area Under the Receiver Operating Characteristic (AUROC) curve and Fold Change for each gene in a given set of consensus genes.
#
# @param df A data frame where rows represent samples and columns represent genes.
# @param cons_genes A character vector containing the names of the consensus genes to be tested.
# @param class A vector of class labels for each sample (binary, typically 0 and 1).
#
# @return A data frame with the following columns:
#   - `gene`: The name of each consensus gene.
#   - `auroc`: The AUROC value for each gene.
#   - `auroc_upper`: The upper limit of the AUROC confidence interval.
#   - `auroc_lower`: The lower limit of the AUROC confidence interval.
#   - `FC`: The Fold Change for each gene, calculated as the difference in mean expression between the two classes.
#
# @details
# This function calculates the AUROC for each gene using the pROC package. The Fold Change is computed as the difference in mean expression between the two classes.
#
# @examples
# \dontrun{
# auroc_results <- auroc_fc_calculation(df = expression_data, cons_genes = gene_list, class = labels)
# }
#
# --- AUROC and Fold Change Calculation Function ---
auroc_fc_calculation <- function(df, cons_genes, class) {
  # Filter dataframe to keep only consensus genes
  df <- df[, colnames(df) %in% cons_genes]

  # Empty dataframe to fill
  results <- data.frame(
    gene = character(),
    auroc = numeric(),
    auroc_upper = numeric(),
    auroc_lower = numeric(),
    FC = numeric(),
    stringsAsFactors = FALSE
  )

  # Looping in the gene
  for (gene in colnames(df)) {
    predictor <- df[[gene]]

    # Calculate AUROC
    roc_obj <- pROC::roc(class, predictor, auc = TRUE, ci = TRUE, direction = "<", levels = c(0, 1))

    # Calculate Fold Change
    FC <- mean(predictor[class == 1]) - mean(predictor[class == 0])

    # Append results
    results <- rbind(results, data.frame(
      gene = gene,
      auroc = roc_obj$auc[1],
      auroc_upper = roc_obj$ci[1],
      auroc_lower = roc_obj$ci[3],
      FC = FC
    ))
  }

  return(results)

}

# Heatmap Calculation with Hierarchical Clustering
#
# @description
# Computes a hierarchical clustering heatmap matrix for a set of consensus genes in a gene expression data frame.
#
# @param df A data frame where rows represent samples and columns represent genes.
# @param cons_genes A character vector containing the names of genes to be used in the heatmap.
# @param class A factor or numeric vector representing the class labels for the samples.
# @param hclust.method A character string specifying the hierarchical clustering method (default is "complete").
# @param distance.method A character string specifying the distance method for clustering (default is "euclidean").
#
# @return A list containing the following elements:
#   - `matrix`: The transposed data frame (genes x samples) filtered by consensus genes.
#   - `hc_rows`: The hierarchical cluster object for rows (genes).
#   - `hc_cols`: The hierarchical cluster object for columns (samples).
#   - `side_colors`: A data frame of class colors for the samples.
#
# @examples
# \dontrun{
# heatmap_data <- heatmap_calculation(df = expression_data, cons_genes = gene_list, class = labels)
# }
heatmap_calculation <- function(df, cons_genes, class,
                                hclust.method = "complete", distance.method = "euclidean") {
  # Filter the dataframe for the consensus genes
  df <- df[, colnames(df) %in% cons_genes, drop = FALSE]

  # Transpose the matrix (genes x samples)
  df_t <- as.data.frame(t(df))

  # Verificando se o número de classes é igual ao número de colunas
  if (length(class) != ncol(df_t)) {
    cli::cli_alert_danger("The number of class labels ({length(class)}) does not match the number of samples (columns) in the heatmap matrix ({ncol(df_t)}).")
    stop("The length of 'class' must match the number of columns in the heatmap matrix.")
  }

  # Create side colors based on the class (matching columns)
  side_colors_df <- data.frame(class = as.factor(class))
  rownames(side_colors_df) <- colnames(df_t)

  # Perform hierarchical clustering
  cli::cli_alert_info("Clustering rows and columns using {hclust.method} ...")
  hc_rows <- stats::hclust(stats::dist(df_t, method = distance.method), method = hclust.method)
  hc_cols <- stats::hclust(stats::dist(t(df_t), method = distance.method), method = hclust.method)
  cli::cli_alert_success("Rows and columns clustered successfully!")

  # Return the results
  return(list(
    matrix = df_t,
    hc_rows = hc_rows,
    hc_cols = hc_cols,
    side_colors = side_colors_df
  ))
}

# MLP Model Training and Variable Importance Calculation
#
# This helper function tunes and evaluates a multilayer perceptron (MLP) model
# using a dataset of consensus genes. It returns the variable importance
# and key performance metrics (accuracy, AUROC, Brier score) on the test set.
#
# Args:
#   df: Data frame with gene expression data (samples x genes).
#   cons_genes: Character vector of consensus gene names.
#   class: Factor or numeric vector of class labels (binary outcome).
#
# Returns:
#   List with:
#     - importances: Data frame with variable importances (scaled from -1 to 1).
#     - test_performance: Data frame with performance metrics on the test set.
#
# Example:
#   res <- mlp_model_calculation(df = gene_data, cons_genes = consensus_genes, class = class_labels)
#   head(res$importances)
#   res$test_performance
mlp_model_calculation <- function(df, cons_genes, class) {
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
      #TODO: NOT WORKING, IS THIS NECESSARY? THERE'S ALREADY A LOT OF THINGS PRINTING
      # mlp_tune_results %>%
      #   workflowsets::collect_metrics()

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

  # Return the results
  return(list(
    importances = importances,
    test_performance = test_performance
  ))

}

#' @title Consensus Gene Analysis with PCA, AUROC, Heatmap, and MLP
#' Test Consensus Genes with PCA, AUROC, and Other Analyses
#'
#' @description
#' This function performs a comprehensive analysis of consensus genes in a gene expression dataset, including:
#' - Principal Component Analysis (PCA) for dimensionality reduction and feature selection.
#' - AUROC and Fold Change analysis to assess the predictive power of each gene.
#' - Heatmap generation with hierarchical clustering to visualize gene expression patterns.
#' - Multi-Layer Perceptron (MLP) model for further classification performance evaluation.
#'
#' The results are returned as an S4 object of class `testConsensus.obj` for easy access and further analysis.
#'
#' @param df.count A numeric data frame where rows represent samples and columns represent genes (gene expression counts).
#' Each value should represent the expression level of a gene for a given sample.
#'
#' @param gene.list A character vector containing the names of consensus genes to be tested.
#' These gene names must match the column names in `df.count`.
#'
#' @param class A binary vector (0 and 1) indicating the class labels for each sample.
#' The length of this vector must match the number of rows in `df.count`.
#'
#' @param top_loadings An integer specifying the number of top genes to select based on absolute loadings in PCA.
#' Default is 15. Must be greater than zero and less than or equal to the total number of genes in `gene.list`.
#'
#' @param display_plots A logical value indicating whether to display the generated plots (PCA, AUROC, Heatmap, and MLP).
#' Default is `TRUE`. Set to `FALSE` to suppress the automatic display of plots.
#' Regardless of this setting, all plots will still be saved in the returned S4 object.
#'
#' @return An S4 object of class `testConsensus.obj`, containing:
#'
#' - `data`: A list with the input data (`df.count`, `gene.list`, `class`).
#' - `pca`: A list with the PCA results:
#'   - `scores`: The PCA-transformed coordinates of each sample.
#'   - `loadings`: The PCA loadings of each gene.
#'   - `top_loadings`: The selected top genes by absolute loading values.
#'   - `explained_variance`: The proportion of variance explained by each principal component.
#'   - `plot`: A ggplot object showing the PCA plot.
#'
#' - `auroc_fc`: A list with the AUROC and Fold Change results:
#'   - `data`: A data frame with AUROC and Fold Change values for each gene.
#'   - `plot`: A ggplot object showing the AUROC and FC plot.
#'
#' - `heatmap`: A list with the Heatmap and Hierarchical Clustering (HC) results:
#'   - `data`: A matrix of normalized gene expression values.
#'   - `plot`: A heatmap plot object for gene expression.
#'
#' - `mlp`: A list with the Multi-Layer Perceptron (MLP) model results:
#'   - `importances`: A data frame with the feature importances.
#'   - `test_performance`: A summary of model performance on the test data.
#'   - `plot`: A ggplot object showing the MLP model performance.
#'
#' @details
#' The `testConsensus` function performs a comprehensive analysis of consensus genes in a gene expression dataset.
#' It follows the steps below:
#'
#' 1. **Principal Component Analysis (PCA):**
#'    - Computes the PCA on the specified consensus genes.
#'    - Selects the top `top_loadings` genes based on absolute loading values.
#'    - Returns the PCA scores, loadings, explained variance, and a PCA plot.
#'
#' 2. **AUROC and Fold Change Analysis:**
#'    - Calculates the Area Under the Receiver Operating Characteristic Curve (AUROC) for each gene.
#'    - Computes the Fold Change between the two classes (0 and 1).
#'    - Provides a data frame with AUROC and FC values and a corresponding plot.
#'
#' 3. **Heatmap with Hierarchical Clustering (HC):**
#'    - Generates a heatmap of gene expression values for the consensus genes.
#'    - Applies hierarchical clustering to visualize gene expression patterns.
#'
#' 4. **Multi-Layer Perceptron (MLP) Model:**
#'    - Trains an MLP model using the input data.
#'    - Calculates feature importance scores and evaluates model performance.
#'    - Provides a model performance summary and a plot.
#'
#' The function organizes all results in an S4 object of class `testConsensus.obj` for easy access and further analysis.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' set.seed(123)
#' expression_data <- matrix(rnorm(2000), nrow = 100, ncol = 20)
#' colnames(expression_data) <- paste0("Gene", 1:20)
#' rownames(expression_data) <- paste0("Sample", 1:100)
#' gene_list <- paste0("Gene", 1:10)
#' class <- sample(c(0, 1), 100, replace = TRUE)
#'
#' # Running the analysis
#' test_obj <- testConsensus(
#'     df.count = expression_data,
#'     gene.list = gene_list,
#'     class = class,
#'     top_loadings = 5,
#'     display_plots = TRUE
#' )
#'
#' # Displaying the full results
#' print(test_obj)
#'
#' # Accessing specific results
#' test_obj@pca$scores       # PCA Scores
#' test_obj@auroc_fc$data    # AUROC and Fold Change Data
#' test_obj@heatmap$plot     # Heatmap Plot
#' test_obj@mlp$plot         # MLP Performance Plot
#' }
#'
#' @export
testConsensus <- function(df.count, gene.list, class, top_loadings = 15, display_plots = TRUE) {

  cli::cli_h2("Consensus Genes Testing")

  # Storing input data in the 'data' slot
  data_list <- list(
    df.count = df.count,
    gene.list = gene.list,
    class = class
  )

  # -------- PCA Analysis --------
  cli::cli_h3("Principal Component Analysis")
  pca_results <- tryCatch(
    pca_calculation(df.count, gene.list, top_loadings),
    error = function(e) cli::cli_abort(paste("Error in PCA analysis:", e$message))
  )

  cli::cli_alert_success("PCA analysis completed successfully!")

  # ----- PCA Plotting -----
  pca_plot <- pca.plot(pca_results$scores, pca_results$top_loadings, class)

  # Adding PCA Plot to Results
  pca_results$plot <- pca_plot

  cli::cli_alert_success("PCA plot created!")

  if (display_plots) {
    print(pca_results$plots$pca_plot)
    print(pca_results$plots$loadings_plot)
  }

  # -------- AUROC and FC Analysis --------

  cli::cli_h3("AUROC and Fold Change Analysis")

  auroc_fc.calc <- tryCatch(
    auroc_fc_calculation(df.count, gene.list, class),
    error = function(e) cli::cli_abort(paste("Error in AUROC and FC analysis:", e$message))
  )

  cli::cli_alert_success("AUROC and FC analysis completed successfully!")

  # ----- AUROC and FC Plotting -----
  auroc_fc.plot <- auroc.fc.plot(auroc_fc.calc)

  cli::cli_alert_success("AUROC and FC plot created!")

  if (display_plots) plot(auroc_fc.plot)

  #Gathering the results in one place
  auroc_fc_results <- list(
    data = auroc_fc.calc,
    plot = auroc_fc.plot
  )

  # -------- Heatmap with HC Analysis --------
  cli::cli_h3("Heatmap with HC Analysis")
  heatmap_results <- tryCatch(
    heatmap_calculation(df.count, gene.list, class),
    error = function(e) cli::cli_abort(paste("Error in Heatmap analysis:", e$message))
  )

  cli::cli_alert_success("Heatmap analysis completed successfully!")

  # ----- Heatmap with HC Plotting -----

  heatmap_plot <- heatmap.plot(heatmap_results)

  # Adding PCA Plot to Results
  heatmap_results$plot <- heatmap_plot

  cli::cli_alert_success("Heatmap plot created!")

  #TODO: Fix the Print plot in this heatmap
  if (display_plots) print(heatmap_plot)

  # -------- MLP Model Analysis --------

  cli::cli_h3("Validation MLP Model")

  mlp_results <- tryCatch(
    mlp_model_calculation(df.count, gene.list, class),
    error = function(e) cli::cli_abort(paste("Error in MLP Model:", e$message))
  )

  cli::cli_alert_success("MLP Model completed successfully!")

  # ----- MLP Model Plotting -----
  mlp.plot <- mlp.model.plot(mlp_results$importances, mlp_results$test_performance)

  # Adding PCA Plot to Results
  mlp_results$plot <- mlp.plot

  cli::cli_alert_success("MLP plot created!")

  if (display_plots) plot(mlp.plot)

  # -------- Finishing --------

  # Creating the S4 object for testConsensus
  testConsensus_obj <- methods::new(
    "testConsensus.obj",
    data = data_list,
    pca = pca_results,
    auroc_fc = auroc_fc_results,
    heatmap = heatmap_results,
    mlp = mlp_results
  )

  # Return the S4 object
  return(testConsensus_obj)

}
