#' Class for the runClassifiers object
#' @description
#' Create a class for the runClassifiers object.
#'
#' @slot data A list containing the adjusted data and metadata.
#' @slot models.info A list containing the finalized workflows for each model.
#' @slot model.features A list containing the variable importances for each model.
#' @slot performances A list containing the tuning and final metrics for each model.
#' @slot predictions A data frame containing the predictions for each model.
#' @export
methods::setClass("runClassifiers.obj",
  slots = list(
    data = "list",
    models.info = "list",
    model.features = "list",
    performances = "list",
    predictions = "data.frame"
  )
)

### ------------------------ Helper Functions ------------------------ ###
#' Check if required packages are installed
#'
#' @param pkgs Character vector of package names to check.
#' @return Invisibly TRUE if all packages installed, otherwise stops with a CLI alert.
#' @keywords internal
check_required_pkgs <- function(pkgs) {
  not_installed <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(not_installed) > 0) {
    cli::cli_abort(
      "The following packages are required but not installed: {paste(not_installed, collapse = ', ')}.
      Please install them with {.code install.packages()} before running this function."
    )
  }
  invisible(TRUE)
}



# Create model specifications
create_model_specs <- function() {
  list(
    xgboost = parsnip::boost_tree(
      mode = "classification", engine = "xgboost",
      mtry = tune(), trees = tune(),
      min_n = tune(), tree_depth = tune(),
      learn_rate = tune(), loss_reduction = tune(),
      sample_size = tune(), stop_iter = tune()
    ),
    bag_tree = parsnip::bag_tree(
      mode = "classification", engine = "rpart",
      tree_depth = tune(), min_n = tune()
    ),
    lightGBM = parsnip::boost_tree(
      mode = "classification", engine = "lightgbm",
      mtry = tune(), trees = tune(),
      min_n = tune(), tree_depth = tune(),
      learn_rate = tune(), loss_reduction = tune()
    ),
    pls = parsnip::pls(
      mode = "classification", num_comp = tune(),
      predictor_prop = tune(), engine = "mixOmics"
    ),
    logistic = parsnip::logistic_reg(
      mode = "classification", engine = "glm"
    ),
    C5_rules = parsnip::C5_rules(
      mode = "classification", engine = "C5.0",
      trees = tune(), min_n = tune()
    ),
    mars = parsnip::set_engine(
      parsnip::mars(
        mode = "classification",
        num_terms = tune(), prod_degree = tune()
      ),
      "earth"
    ),
    bag_mars = parsnip::set_engine(
      parsnip::bag_mars(
        mode = "classification", num_terms = tune(),
        prod_degree = tune()
      ),
      "earth"
    ),
    mlp = parsnip::mlp(
      mode = "classification", engine = "nnet",
      hidden_units = tune(), penalty = tune(),
      epochs = tune()
    ),
    bag_mlp = parsnip::bag_mlp(
      mode = "classification", engine = "nnet",
      hidden_units = tune(), penalty = tune(),
      epochs = tune()
    ),
    decision_tree = parsnip::decision_tree(
      mode = "classification", tree_depth = tune(),
      min_n = tune(), engine = "rpart"
    ),
    rand_forest = parsnip::set_engine(
      parsnip::rand_forest(
        mode = "classification", mtry = tune(),
        trees = 300, min_n = tune()
      ),
      "ranger",
      importance = "impurity"
    ),
    svm_linear = parsnip::svm_linear(
      mode = "classification", engine = "kernlab",
      cost = tune(), margin = tune()
    ),
    svm_poly = parsnip::svm_poly(
      mode = "classification", engine = "kernlab",
      cost = tune(), degree = tune(),
      scale_factor = tune(), margin = tune()
    ),
    svm_rbf = parsnip::svm_rbf(
      mode = "classification", engine = "kernlab",
      cost = tune(), rbf_sigma = tune(),
      margin = tune()
    )
  )
}

# Create dynamic workflow sets
create_dynamic_workflow_sets <- function(models, selector.recipes, data, downsample, boruta.maxRuns, selector.threshold) {
  model_specs <- create_model_specs()

  recipes_base <- recipes::recipe(class ~ ., data = data) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors())

  if (downsample == TRUE) {
    cli::cli_alert_info("Majority class will be downsampled ...")
    recipes_base <- recipes_base %>%
      themis::step_downsample(class)
  }

  # Recipe with Boruta Feature Selection
  recipe_boruta <- recipes_base %>%
    colino::step_select_boruta(recipes::all_predictors(), outcome = "class", options = list(maxRuns = boruta.maxRuns))

  # Recipe with ROC-based Feature Selection
  recipe_ROC <- recipes_base %>%
    colino::step_select_roc(recipes::all_predictors(), outcome = "class", threshold = selector.threshold)

  # Recipe with Information Gain Feature Selection
  recipe_INFGAIN <- recipes_base %>%
    colino::step_select_infgain(recipes::all_predictors(), outcome = "class", threshold = selector.threshold)

  # Recipe with Max Relevancy Min Redundancy Feature Selection
  recipe_MRMR <- recipes_base %>%
    colino::step_select_mrmr(recipes::all_predictors(), outcome = "class", threshold = selector.threshold)

  # Recipe correlation-based Feature Selection
  recipe_corr <- recipes_base %>%
    recipes::step_corr(recipes::all_predictors(), threshold = selector.threshold)

  # Filter to keep only models selected by User
  filtered_specs <- model_specs[models]

  # Map recipes names to respective recipe objects
  recipes_map <- list(
    base   = recipes_base,
    boruta = recipe_boruta,
    roc = recipe_ROC,
    infgain = recipe_INFGAIN,
    mrmr = recipe_MRMR,
    corr = recipe_corr
  )

  # Create workflow_set() for every model and recipe combination
  workflow_sets <- lapply(
    seq_along(filtered_specs),
    function(i) {
      model_name <- names(filtered_specs)[i]
      recipe_name <- selector.recipes[i]
      recipe_obj <- recipes_map[[recipe_name]]

      # Create the workflow_set with a custom wflow_id
      workflow_set <- workflowsets::workflow_set(
        preproc = list(recipe_obj),
        models = list(filtered_specs[[model_name]]),
        cross = TRUE
      )
      workflow_set$wflow_id <- paste0(model_name, "_", recipe_name)
      workflow_set
    }
  )

  # Combine all workflow_set() with rbind()
  combined_workflow_set <- do.call(rbind, workflow_sets)

  return(combined_workflow_set)
}

# Helper function to tune and fit models
tune_and_fit <- function(
    models, selector.recipes,
    tuning.method, n, metric, train_resamples, data, downsample, boruta.maxRuns, selector.threshold) {
  # Create dynamic workflow sets
  tune_workflows <- create_dynamic_workflow_sets(models, selector.recipes, data = data, downsample = downsample, boruta.maxRuns = boruta.maxRuns, selector.threshold = selector.threshold)

  cli::cli_h2("Tuning Model Parameters")
  # Tune the models
  if (tuning.method %in% c("tune_grid", "tune_race_anova", "tune_race_win_loss")) {
    tune_results <- tune_workflows %>%
      workflowsets::workflow_map(
        fn = tuning.method,
        resamples = train_resamples,
        grid = n,
        metrics = yardstick::metric_set(
          yardstick::accuracy, yardstick::precision, yardstick::recall, yardstick::f_meas
        ),
        verbose = TRUE
      )
  } else if (tuning.method %in% c("tune_bayes", "tune_sim_anneal")) {
    tune_results <- tune_workflows %>%
      workflowsets::workflow_map(
        fn = tuning.method,
        resamples = train_resamples,
        iter = n,
        metrics = yardstick::metric_set(
          yardstick::accuracy, yardstick::precision, yardstick::recall, yardstick::f_meas
        ),
        verbose = TRUE
      )
  } else {
    cli::cli_abort("Invalid tuning method. Please refer to the documentation for valid options.")
  }

  # Consider as successful those results where collect_metrics() returns a tibble with rows
  is_valid_tune_result <- function(res) {
    met <- tryCatch(tune::collect_metrics(res), error = function(e) NULL)
    !is.null(met) && nrow(met) > 0
  }

  successful_idx <- which(
    sapply(tune_results$result, is_valid_tune_result)
  )
  successful_ids <- tune_results$wflow_id[successful_idx]

  if (length(successful_ids) == 0) {
    cli::cli_abort("All models failed during tuning. Execution stopped.")
  }

  cli::cli_h3("Best Parameters Selection")
  cli::cli_alert_info("Selecting best parameters for each model ...")
  # Select the best parameters for each model
  # Select best parameters only for successful workflows
  best_params <- lapply(
    successful_ids,
    function(id) {
      tryCatch(
        tune::select_best(
          workflowsets::extract_workflow_set_result(tune_results, id = id),
          metric = metric
        ),
        error = function(e) {
          cli::cli_alert_warning("Model {id} failed in select_best: {conditionMessage(e)}")
          NULL
        }
      )
    }
  )
  names(best_params) <- successful_ids

  # Filter out any models where select_best failed (best_params is NULL)
  valid_idx <- which(!sapply(best_params, is.null))
  valid_ids <- successful_ids[valid_idx]
  best_params <- best_params[valid_idx]

  if (length(valid_ids) == 0) {
    cli::cli_abort("All models failed during parameter selection. Execution stopped.")
  }

  cli::cli_alert_success("Succesfully selected parameters!")

  cli::cli_h3("Workflows Finalization")
  cli::cli_alert_info("Updating workflows with best parameters ...")
  # Finalize workflows only for valid models
  final_workflows <- lapply(
    seq_along(valid_ids),
    function(i) {
      tryCatch(
        tune::finalize_workflow(
          workflowsets::extract_workflow(tune_workflows, id = valid_ids[i]),
          best_params[[i]]
        ),
        error = function(e) {
          cli::cli_alert_warning("Model {valid_ids[i]} failed in finalize_workflow: {conditionMessage(e)}")
          NULL
        }
      )
    }
  )
  names(final_workflows) <- valid_ids

  # Filter out any models where finalize_workflow failed (workflow is NULL)
  final_idx <- which(!sapply(final_workflows, is.null))
  final_ids <- valid_ids[final_idx]
  final_workflows <- final_workflows[final_idx]

  if (length(final_ids) == 0) {
    cli::cli_abort("All models failed during workflow finalization. Execution stopped.")
  }

  cli::cli_alert_success("Succesfully finalized the workflows!")

  cli::cli_h2("Model Fitting")
  cli::cli_alert_info("Fitting models on the entire dataset ...")

  # Fit only successful workflows
  last_fit_results <- lapply(
    final_workflows,
    function(workflow) {
      tryCatch(
        parsnip::fit(workflow, data = data),
        error = function(e) {
          cli::cli_alert_warning("Model failed during fitting: {conditionMessage(e)}")
          NULL
        }
      )
    }
  )
  names(last_fit_results) <- final_ids

  # Filter out failed fits
  fit_idx <- which(!sapply(last_fit_results, is.null))
  fit_ids <- final_ids[fit_idx]
  last_fit_results <- last_fit_results[fit_idx]
  final_workflows <- final_workflows[fit_idx]

  if (length(fit_ids) == 0) {
    cli::cli_abort("All models failed during final fitting. Execution stopped.")
  }

  cli::cli_alert_success("Models fitted succesfully.")

  cli::cli_alert_info("Computing metric performances ...")
  # Compute tuning metrics only for models that made it to final fitting
  tune_results_valid <- tune_results[tune_results$wflow_id %in% valid_ids, ]
  tuning_metrics <- as.data.frame(
    workflowsets::rank_results(
    tune_results_valid,
    rank_metric = metric, select_best = TRUE
  ))
  tuning_metrics <- tuning_metrics[tuning_metrics$wflow_id %in% fit_ids, ]

  # Define metrics to compute on test set
  multi_met <- yardstick::metric_set(
    yardstick::accuracy, yardstick::f_meas,
    yardstick::precision, yardstick::recall
  )
  # Compute predictions and metrics only for fitted models
  results <- lapply(seq_along(last_fit_results), function(i) {
    model <- last_fit_results[[i]]
    model_name <- names(last_fit_results)[i]
    tryCatch({
      # Prediction step
      preds <- stats::predict(model, new_data = data)
      predictions <- data.frame(
        preds,
        class = data$class,
        ID = row.names(data),
        model = model_name
      )
      # Calculate metrics using multi_met()
      metrics <- multi_met(predictions, truth = class, estimate = .pred_class)
      metrics$model <- model_name
      # Return a list with predictions and metrics
      list(predictions = predictions, metrics = metrics)
    }, error = function(e) {
      cli::cli_alert_warning("Model {model_name} failed during prediction or metric calculation: {conditionMessage(e)}")
      NULL
    })
  })
  names(results) <- names(last_fit_results)

  # Filter out models that failed in the prediction/metric step
  results <- Filter(Negate(is.null), results)

  # Extract predictions and metrics from results list
  predictions_df <- do.call(rbind, lapply(results, function(x) x$predictions))
  test_metrics <- do.call(rbind, lapply(results, function(x) x$metrics))

  cli::cli_alert_success("Metrics computed succesfully.")
  cli::cli_alert_success("Succesfully accomplished model fitting!")

  # Update the object of class runClassifiers.obj
  runClassifiers.obj <- methods::new(
    "runClassifiers.obj",
    data = list(),
    models.info = final_workflows,
    model.features = list(),
    performances = list(
      tuning_metrics = tuning_metrics,
      final_metrics = test_metrics
    ),
    predictions = predictions_df
  )

  return(list(runClassifiers.obj, last_fit_results, predictions_df))

}

# Helper function to calculate VIP for all models in last_fit_results
calculate_vip <- function(last_fit_results, test_x, test_y, n_sim, parallel = parallel) {
  # Wrapper function for predictions with models
  pfun <- function(model, new_data) {
    if ("ksvm" %in% class(model) || "_ksvm" %in% class(model)) {
      return(kernlab::predict(model, newdata = new_data, type = "probabilities")[, 2])
    } else if ("xrf" %in% class(model)) {
      pred <- stats::predict(model, new_data, type = "response")
      return(as.numeric(if (is.matrix(pred)) pred[, 1] else pred))
    } else if ("rda" %in% class(model) || "_rda" %in% class(model)) {
      pred <- stats::predict(model, new_data)
      return(as.numeric(pred$.pred_class))
    } else {
      probs <- stats::predict(model, new_data, type = "prob")
      if (".pred_1" %in% colnames(probs)) {
        return(probs$.pred_1)
      } else {
        cli::cli_alert_danger("The model did not return a valid probability.")
      }
    }
  }

  # Initialize an empty list to store results
  vip_list <- list()

  # Iterate only over non-NULL models
  for (name in names(last_fit_results)) {
    model_obj <- last_fit_results[[name]]
    if (is.null(model_obj)) {
      cli::cli_alert_warning("Skipping model {name} (fit object is NULL).")
      next
    }

    cli::cli_alert_info("Processing model {name} ...")

    # Try extracting the parsnip fit, catch errors
    model_fit <- tryCatch({
      workflows::extract_fit_parsnip(model_obj)$fit
    }, error = function(e) {
      cli::cli_alert_warning("Failed to extract parsnip fit for model {name}: {conditionMessage(e)}")
      return(NULL)
    })
    if (is.null(model_fit)) next

    # Try to compute variable importance
    vip_data <- tryCatch({
      if (inherits(model_fit, "bagger")) {
        cli::cli_alert_info("Model {name} is a bagger. Using baguette::var_imp() ...")
        vip_data <- baguette::var_imp(model_fit)
        colnames(vip_data) <- c("Variable", "Importance", "std.error", "used")
        vip_data[vip_data$Importance != 0, ]
      } else {
        num_features <- ncol(test_x)
        vip_result <- vip::vip(model_fit, num_features = num_features)$data
        vip_result[vip_result$Importance != 0, ]
      }
    }, error = function(e) {
      cli::cli_alert_info("Direct VIP failed for model {name}: {conditionMessage(e)}. Using permutation-based method...")
      # Try permutation-based VIP, catch error if even that fails
      tryCatch({
        num_features <- ncol(test_x)
        vip_result <- vip::vip(
          object = model_fit,
          method = "permute",
          parallel = parallel,
          nsim = n_sim,
          metric = "roc_auc",
          pred_wrapper = function(object, newdata) pfun(object, newdata),
          train = test_x,
          target = test_y,
          event_level = "second",
          num_features = num_features
        )$data
        vip_result[vip_result$Importance > 0, ]
      }, error = function(e2) {
        cli::cli_alert_warning("Permutation VIP failed for model {name}: {conditionMessage(e2)}. Skipping model.")
        return(NULL)
      })
    })

    # Only store VIP data if it exists and is not empty
    if (!is.null(vip_data) && nrow(vip_data) > 0) {
      vip_list[[name]] <- vip_data
      cli::cli_alert_success("Variable importance computed for model {name}!")
    } else {
      cli::cli_alert_warning("Variable importance not computed for model {name} (empty or failed).")
    }
  }

  return(vip_list)
}

# TODO include example in documentation using tables instead of the preProcess.obj
#' @title runClassifiers
#' @description This function run classifiers specified on an object of class preProcess.obj or
#' on the data provided in the arguments ...
#'
#' @param preProcess.obj An object of class preProcess.
#' @param models A character vector specifying the classifiers to be used. Supported models include:
#'
#'   - `"xgboost"`: Extreme Gradient Boosting, an ensemble method using boosted trees.
#'
#'   - `"bag_tree"`: Bagged decision trees, a bootstrapped ensemble of tree models.
#'
#'   - `"lightGBM"`: A fast gradient boosting method optimized for large datasets.
#'
#'   - `"pls"`: Partial Least Squares regression.
#'
#'   - `"logistic"`: Logistic regression, a simple linear classifier.
#'
#'   - `"C5_rules"`: Rule-based classifier using C5.0 decision trees.
#'   - `"mars"`: Multivariate Adaptive Regression Splines, a flexible non-linear regression method.
#'
#'   - `"bag_mars"`: Bagged version of MARS for increased stability.
#'
#'   - `"mlp"`: Multi-Layer Perceptron, a basic feedforward neural network.
#'
#'   - `"bag_mlp"`: Bagged version of MLP to reduce variance.
#'
#'   - `"decision_tree"`: A single decision tree.
#'
#'   - `"rand_forest"`: Random forest, an ensemble of decision trees with randomized splits.
#'
#'   - `"svm_linear"`: Support Vector Machine with a linear kernel.
#'
#'   - `"svm_poly"`: Support Vector Machine with a polynomial kernel.
#'
#'   - `"svm_rbf"`: Support Vector Machine with a radial basis function (RBF) kernel.
#' @param selector.recipes A character vector specifying the feature selection methods to be applied before classification.
#'
#' Supported selection strategies include:
#'
#'   - `"base"`: Uses only the default pre-processing steps (normalization, near-zero variance removal), with no feature selection.
#'
#'   - `"boruta"`: Wrapper method that iteratively removes unimportant features using random forest.
#'
#'   - `"roc"`: Selects features based on Receiver Operating Characteristic (ROC) AUC scores.
#'
#'   - `"infgain"`: Uses Information Gain to select the most informative features.
#'
#'   - `"mrmr"`: Minimum Redundancy Maximum Relevance (mRMR), selects features that maximize relevance and minimize redundancy.
#'
#'   - `"corr"`: Filters features based on correlation thresholds to reduce multicollinearity.
#'
#' `models` and `selector.recipes` will be paired in the order they are provided. If the number of recipes
#' does not match the number of models, the first recipe will be used for all models.
#' @param tuning.method A character string specifying the hyperparameter tuning strategy. Options include:
#'   - `"tune_grid"` (default): Grid search across a pre-defined set of hyperparameters.
#'   - `"tune_race_anova"`: Adaptive search using ANOVA-based pruning to speed up tuning.
#'   - `"tune_race_win_loss"`: Win-loss-based adaptive search that discards weak candidates early.
#'   - `"tune_bayes"`: Bayesian optimization to intelligently explore hyperparameter space.
#'   - `"tune_sim_anneal"`: Simulated annealing, a probabilistic approach to global optimization.
#' @param n An integer specifying the number of iterations for the tuning method.
#' @param v An integer specifying the number of folds for the cross-validation during the hyperparameters tuning.
#' @param boruta.maxRuns maxRuns parameter to pass to boruta feature selector when it is selected.
#' @param selector.threshold Threshold parameter used for feature selection (Applied on "roc", "infgain", "mrmr", "corr" selectors).
#' @param metric A character string specifying the metric to be used for tuning.
#' @param nsim An integer specifying the number of simulations for the permutation-based VIP.
#' @param seed An integer specifying the seed for reproducibility.
#' @param filter A logical specifying whether to filter genes not annotated in the GO and KEGG databases.
#' @param downsample A logical specifying whether to downsample the majority class.
#' @param plot A logical specifying whether to generate plots. Default is TRUE.
#' @param ... Arguments passed to the function. If preProcess.obj is not provided, the function looks
#' for df.count, df.clin and class in ... as specified in the documentation of preProcess function.
#' @param parallel A logical specifying whether to use parallel computation for model fitting and variable importance. Default is FALSE.
#'
#'   - TRUE: Use all available cores for parallel computation.
#'
#'   - FALSE: Run all computations sequentially.
#'
#'   This affects both model tuning and the permutation-based variable importance calculation.
#'
#' @return An object of class `runClassifiers.obj` containing:
#'   - `data`: A list containing the adjusted data and metadata.
#'   - `models.info`: A list containing the finalized workflows for each model.
#'   - `model.features`: A list containing the variable importance for each model.
#'   - `performances`: A list containing the tuning and final metrics for each model.
#'   - `predictions`: A data frame containing the predictions for each model.
#'
#' @details
#' The function runs specified classifiers paired with provided selector.recipes on pre-processed data.
#' If filter = TRUE, the function first filters the genes that are not annotated in the GO and KEGG databases.
#' Then, the function tunes and fits the models using the specified tuning method.
#' Finally it computes the variable importances for each model using the permutation-based VIP,
#' when the direct VIP computation fails.
#'
#' @import baguette
#' @import bonsai
#' @import plsmod
#' @import rules
#' @import nnet
#' @import NeuralNetTools
#' @import Boruta
#' @import praznik
#' @import FSelectorRcpp
#' @importFrom colino step_select_boruta step_select_roc step_select_infgain step_select_mrmr
#'
#' @examples
#' \dontrun{
#' rc <- runClassifiers(preProcess.obj,
#'   models = c("bag_mlp", "rand_forest", "svm_poly"),
#'   selector.recipes = "boruta", tuning.method = "tune_grid", n = 5,
#'   v = 3, metric = "accuracy", nsim = 2, seed = 123
#' )
#' }
#'
#' @export
runClassifiers <- function(
    preProcess.obj = NULL, ..., models = c("bag_mlp", "rand_forest", "svm_poly"),
    selector.recipes = c("boruta", "roc", "boruta"),
    tuning.method = "tune_grid", n = 5, v = 3,
    boruta.maxRuns = 100, selector.threshold = 0.95,
    metric = "accuracy",
    nsim = 2, filter = TRUE, seed = 123, downsample = FALSE, plot = TRUE,
    parallel = FALSE
    ) {
  #Message to show the start
  cli::cli_h1("runClassifiers")
    #Parallelization configuration for modeling and variable importance.
  if (isTRUE(parallel)) {
    cli::cli_alert_info("Parallelization enabled: running with multiple cores (multisession).")
    future::plan(future::multisession, workers = parallel::detectCores() - 1)
  } else {
    cli::cli_alert_info("Sequential mode: running on a single core.")
    future::plan(future::sequential)
  }

  # Take arguments passed in ...
  dots <- list(...)

  # Check if preProcess.obj is null
  if (is.null(preProcess.obj)) {
    cli::cli_alert_info("preProcess.obj not provided, looking for df.count, df.clin and class in ...")

    if ("df.count" %in% names(dots) && "df.clin" %in% names(dots) && "class" %in% names(dots)) {
      cli::cli_alert_success("Found df.count, df.clin and class in ...")
      # Create a preProcess object
      preProcess.obj <- preProcess(
        df.count = dots$df.count, df.clin = dots$df.clin, class = dots$class,
        is.normalized = TRUE, plot = FALSE
      )
      # cli::cli_alert_success("Successfully created preProcess object from data provided!")
    } else {
      cli::cli_abort("df.count, df.clin and class not found in ... Please refer to the documentation for valid options.")
    }
  }

  # Vector of available classifiers
  available_models <- c(
    "xgboost", "bag_tree", "lightGBM", "pls", "logistic",
    "C5_rules", "mars", "bag_mars", "mlp", "bag_mlp",
    "decision_tree", "rand_forest", "svm_linear", "svm_poly",
    "svm_rbf"
  )

  # Vector of available feature selectors
  available_selectors <- c("base", "boruta", "roc", "infgain", "mrmr", "corr")

  # Map each model to its required package (only non-base R packages)
  model_pkgs <- list(
    xgboost        = "xgboost",
    bag_tree       = "rpart",
    lightGBM       = "lightgbm",
    pls            = "mixOmics",
    logistic       = NULL,   # stats is base R
    C5_rules       = "C50",
    mars           = "earth",
    bag_mars       = "earth",
    mlp            = "nnet",
    bag_mlp        = "nnet",
    decision_tree  = "rpart",
    rand_forest    = "ranger",
    svm_linear     = "kernlab",
    svm_poly       = "kernlab",
    svm_rbf        = "kernlab"
  )
  # Identify which packages are required by the selected models
  required_pkgs <- unique(unlist(model_pkgs[models]))
  required_pkgs <- required_pkgs[!is.null(required_pkgs)]

  # Check if required packages are installed
  check_required_pkgs(required_pkgs)

  withr::with_seed(
    seed = seed,
    code = {
      if (filter == TRUE) {
        ## Keep only genes annotated in GO and KEGG databases
        cli::cli_alert_info("Filtering genes non-annotated in GO and KEGG db ...")
        # Upload annotated genes database
        data(annotated.genes, envir = environment())
        # Compute the number of columns before filtering. class col is excluded from computation
        ncol_pre <- ncol(preProcess.obj@processed$adjusted.data) - 1
        # Save class column to cbind it after filtering
        class <- preProcess.obj@processed$adjusted.data$class
        # Filter genes that are not present in the DB
        preProcess.obj@processed$adjusted.data <-
          preProcess.obj@processed$adjusted.data[, colnames(preProcess.obj@processed$adjusted.data) %in%
            annotated.genes$genes_to_retrieve]
        # cbind class column that has been filtered out
        preProcess.obj@processed$adjusted.data <- cbind(preProcess.obj@processed$adjusted.data, class)
        # Compute the number of columns after filter. class col is excluded from computation
        ncol_post <- ncol(preProcess.obj@processed$adjusted.data) - 1
        cli::cli_alert_success("Filtered {ncol_pre - ncol_post} genes. Total number of genes is now: {ncol_post}!")
      } else {
        preProcess.obj@processed$adjusted.data <- preProcess.obj@processed$adjusted.data
        cli::cli_alert_info("Skipping gene filtering using {ncol(preProcess.obj@processed$adjusted.data) - 1} genes ...")
      }

      # Split data into v resamples for tuning inside cross-validation
      data <- preProcess.obj@processed$adjusted.data
      train_resamples <- rsample::vfold_cv(data, v = v, strata = class)


      if (length(selector.recipes) >= 1 && length(selector.recipes) != length(models)) {
        cli::cli_alert_info("The number of recipes does not match the number of models. Using the first recipe for all models.")
        selector.recipes <- rep(selector.recipes[1], length(models))
      }

      if (!any(models %in% available_models)) {
        cli::cli_abort("Invalid model(s) selected. Please refer to the documentation for valid options.")
      }

      if (!any(selector.recipes %in% available_selectors)) {
        cli::cli_abort("Invalid selector recipe(s) selected. Please refer to the documentation for valid options.")
      }

      # Tune and fit the models
      t_and_f_output <- tune_and_fit(
        models = models,
        selector.recipes = selector.recipes,
        tuning.method = tuning.method,
        n = n,
        metric = metric,
        train_resamples = train_resamples,
        data = data,
        downsample = downsample,
        boruta.maxRuns = boruta.maxRuns,
        selector.threshold = selector.threshold
      )

      cli::cli_h2("Variable Importances")
      vip_results <- calculate_vip(
        last_fit_results = t_and_f_output[[2]],
        test_x = data[, -ncol(data)],
        test_y = data$class,
        n_sim = nsim,
        parallel = parallel
      )
    }
  )

  obj <- t_and_f_output[[1]]
  obj@model.features <- vip_results
  future::plan(future::sequential)

  if (plot == TRUE) {
    cli::cli_h2("Plots")

    # UpSet plot - needs at least 2 models/features
    if (!is.null(obj@model.features) && length(obj@model.features) >= 2) {
      cli::cli_alert_info("Generating UpSet plot ...")
      tryCatch({
        up <- upset.plot(obj)
        print(up)
        cli::cli_alert_success("Successfully generated UpSet plot!")
      }, error = function(e) {
        cli::cli_alert_warning("Failed to generate UpSet plot: {conditionMessage(e)}")
      })
    } else {
      cli::cli_alert_warning("UpSet plot requires at least 2 valid models with variable importance. Skipping.")
    }

    # Performances plot - needs at least 1 model with metrics
    has_perf <- !is.null(obj@performances$final_metrics) && nrow(obj@performances$final_metrics) > 0
    if (has_perf) {
      cli::cli_alert_info("Generating performances plot ...")
      tryCatch({
        perf <- performances.plot(obj)
        print(perf)
        cli::cli_alert_success("Successfully generated performances plot!")
      }, error = function(e) {
        cli::cli_alert_warning("Failed to generate performances plot: {conditionMessage(e)}")
      })
    } else {
      cli::cli_alert_warning("No valid model performance available for plotting. Skipping performances plot.")
    }

    # Wrong predictions heatmap - needs at least 1 prediction
    has_preds <- !is.null(obj@predictions) && nrow(obj@predictions) > 0
    if (has_preds) {
      cli::cli_alert_info("Generating predictions heatmap ...")
      tryCatch({
        wp <- wrong.preds.plot(obj@predictions)
        print(wp)
        cli::cli_alert_success("Successfully generated predictions heatmap plot!")
      }, error = function(e) {
        cli::cli_alert_warning("Failed to generate predictions heatmap plot: {conditionMessage(e)}")
      })
    } else {
      cli::cli_alert_warning("No predictions available to plot. Skipping predictions heatmap.")
    }
  }

  obj@data$adjusted.data <- preProcess.obj@processed$adjusted.data
  obj@data$metadata <- preProcess.obj@metadata


  cli::cli_alert_success("Successfully executed runClassifiers function!")

  return(obj)
}
