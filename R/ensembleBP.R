#' @export
# Create a class.
methods::setClass("ensBP.obj",
  slots = list(
    models.info = "list",
    model.features = "list",
    performances = "list"
  )
)

#' @importFrom magrittr %>%
#' @import baguette
#' @import bonsai
#' @import plsmod
#' @import rules

# Helper function to create model specifications
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
    mars = parsnip::mars(
      mode = "classification",
      num_terms = tune(), prod_degree = tune()
    ) %>% parsnip::set_engine("earth"),
    bag_mars = parsnip::bag_mars(
      mode = "classification", num_terms = tune(),
      prod_degree = tune()
    ) %>% parsnip::set_engine("earth"),
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
    rand_forest = parsnip::rand_forest(
      mode = "classification", mtry = tune(),
      trees = 300, min_n = tune()
    ) %>%
      parsnip::set_engine("ranger", importance = "impurity"),
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

# Helper function to create dynamic workflow sets
create_dynamic_workflow_sets <- function(models, selector.recipes, train) {
  model_specs <- create_model_specs()

  # Recipe with Boruta Feature Selection
  recipe_boruta <- recipes::recipe(class ~ ., data = train) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors()) %>%
    #recipes::step_corr(recipes::all_predictors(), threshold = 0.8) %>%
    colino::step_select_boruta(recipes::all_predictors(), outcome = "class")

  # Recipe with ROC-based Feature Selection
  recipe_ROC <- recipes::recipe(class ~ ., data = train) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors()) %>%
    #recipes::step_corr(recipes::all_predictors(), threshold = 0.8) %>%
    colino::step_select_roc(recipes::all_predictors(), outcome = "class", threshold = 0.95)

  # Recipe with Information Gain Feature Selection
  recipe_INFGAIN <- recipes::recipe(class ~ ., data = train) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors()) %>%
    #recipes::step_corr(recipes::all_predictors(), threshold = 0.8) %>%
    colino::step_select_infgain(recipes::all_predictors(), outcome = "class", threshold = 0.95)

  # Recipe with Max Relevancy Min Redundancy Feature Selection
  recipe_MRMR <- recipes::recipe(class ~ ., data = train) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors()) %>%
    #recipes::step_corr(recipes::all_predictors(), threshold = 0.8) %>%
    colino::step_select_mrmr(recipes::all_predictors(), outcome = "class", threshold = 0.95)

  # Filter to keep only models selected by User
  filtered_specs <- model_specs[models]

  #Map recipes names to respective recipe objects
  recipes_map <- list(
    boruta = recipe_boruta,
    roc = recipe_ROC,
    infgain = recipe_INFGAIN,
    mrmr = recipe_MRMR
  )

  # Create workflow_set() for every model and recipe combination
  workflow_sets <- purrr::map2(
    names(filtered_specs), selector.recipes,
    ~ {
      # Seleziona la ricetta in base al nome
      recipe_obj <- recipes_map[[.y]]
      model_name <- .x
      recipe_name <- .y

      # Crea il workflow_set con un wflow_id personalizzato
      workflowsets::workflow_set(
        preproc = list(recipe_obj),
        models = list(filtered_specs[[model_name]]),
        cross = TRUE
      ) %>%
        dplyr::mutate(wflow_id = paste0(model_name, "_", recipe_name))
    }
  )

  # Combine all workflow_set() with rbind()
  combined_workflow_set <- do.call(rbind, workflow_sets)

  return(combined_workflow_set)
}

# Helper function to tune and fit models
tune_and_fit <- function(
    models, selector.recipes,
    tuning.method, n, metric, train_resamples, train, train_test_split) {
  # Create dynamic workflow sets
  tune_workflows <- create_dynamic_workflow_sets(models, selector.recipes, train = train)

  cli::cli_h2("Tuning Model Parameters")
  # Tune the models
  # tictoc::tic("Model tuning and fitting")
  if (tuning.method %in% c("tune_grid", "tune_race_anova", "tune_race_win_loss")) {
    tune_results <- tune_workflows %>%
      workflowsets::workflow_map(
        fn = tuning.method,
        resamples = train_resamples,
        grid = n,
        verbose = TRUE
      )
  } else if (tuning.method %in% c("tune_bayes", "tune_sim_anneal")) {
    tune_results <- tune_workflows %>%
      workflowsets::workflow_map(
        fn = tuning.method,
        resamples = train_resamples,
        iter = n,
        verbose = TRUE
      )
  } else {
    cli::cli_abort("Invalid tuning method. Please refer to the documentation for valid options.")
  }
  # tictoc::toc()

  cli::cli_h3("Best Parameters Selection")
  cli::cli_alert_info("Selecting best parameters for each model ...")
  # Select the best parameters for each model
  best_params <- purrr::map(
    tune_results$wflow_id,
    ~ tune::select_best(workflowsets::extract_workflow_set_result(tune_results, id = .x),
      metric = metric
    )
  )
  cli::cli_alert_success("Succesfully selected parameters!")

  cli::cli_h3("Workflows Finalization")
  cli::cli_alert_info("Updating workflows with best parameters ...")
  # Using the best parameters to finalize the workflows
  final_workflows <- purrr::map2(
    tune_results$wflow_id,
    best_params,
    ~ tune::finalize_workflow(
      workflowsets::extract_workflow(tune_workflows, id = .x), .y
    )
  )
  names(final_workflows) <- tune_results$wflow_id
  cli::cli_alert_success("Succesfully finalized the workflows!")

  cli::cli_h2("Model Fitting")
  cli::cli_alert_info("Fitting models on the entire split ...")
  # Fit the finalized workflows on the full training set and test on the test set
  last_fit_results <- purrr::map(
    final_workflows,
    ~ tune::last_fit(.x, split = train_test_split)
  )
  cli::cli_alert_success("Models fitted succesfully.")

  cli::cli_alert_info("Computing performances on Test set ...")
  #cli::cli_h3("Model Performances")
  # Collect the tuning and test metrics
  tuning_metrics <- as.data.frame(workflowsets::rank_results(tune_results, select_best = TRUE))
  test_metrics <- last_fit_results %>%
    purrr::imap_dfr(~ workflowsets::collect_metrics(.x) %>% mutate(model = .y))
  test_metrics <- as.data.frame(test_metrics)
  cli::cli_alert_success("Metrics computed succesfully.")
  cli::cli_alert_success("Succesfully accomplished model fitting!")


  ensBP.obj <- methods::new("ensBP.obj",
    models.info = final_workflows,
    model.features = list(),
    performances = list(tuning_metrics = tuning_metrics, test_metrics = test_metrics)
  )

  return(list(ensBP.obj, last_fit_results))
}

# TODO: remove negative Importances from VIP computed with permutation-based method
# Helper function to calculate VIP for all models in last_fit_results
calculate_vip <- function(last_fit_results, test_x, test_y, n_sim) {
  # Wrapper function for predictions with models
  pfun <- function(model, new_data) {
    if ("ksvm" %in% class(model) || "_ksvm" %in% class(model)) {
      return(kernlab::predict(model, newdata = new_data, type = "probabilities")[, 2])
    } else if ("xrf" %in% class(model)) {
      pred <- predict(model, new_data, type = "response")
      return(as.numeric(if (is.matrix(pred)) pred[, 1] else pred))
    } else if ("rda" %in% class(model) || "_rda" %in% class(model)) {
      pred <- predict(model, new_data)
      return(as.numeric(pred$.pred_class))
    } else {
      probs <- predict(model, new_data, type = "prob")
      if (".pred_1" %in% colnames(probs)) {
        return(probs$.pred_1)
      } else {
        cli::cli_alert_danger("The model did not return a valid probability.")
      }
    }
  }

  # Initialize an empty list to store results
  vip_list <- list()

  # Iteration over all Models contained in last_fit_results
  for (name in names(last_fit_results)) {
    cli::cli_alert_info("Processing model {name} ...")
    model_fit <- workflows::extract_fit_parsnip(last_fit_results[[name]])$fit

    # Check if the model is of class 'bagger'
    if (inherits(model_fit, "bagger")) {
      cli::cli_alert_info("Model {name} is a bagger. Using baguette::var_imp() ...")
      vip_data <- baguette::var_imp(model_fit)
      colnames(vip_data) <- c("Variable", "Importance", "std.error", "used")
      # Keep only non-zero importance values
      vip_data <- dplyr::filter(vip_data, Importance != 0)
    } else {
      # Use VIP or permutation-based VIP
      num_features <- ncol(test_x)

      vip_data <- tryCatch(
        {
          vip::vip(model_fit, num_features = num_features)$data %>%
            # Keep only non-zero importance values
            dplyr::filter(Importance != 0)
        },
        error = function(e) {
          cli::cli_alert_info("Direct VIP failed for model {name}. Using permutation-based method...")
          vip::vip(
            object = model_fit,
            method = "permute",
            parallel = TRUE,
            nsim = n_sim,
            metric = "roc_auc",
            pred_wrapper = function(object, newdata) pfun(object, newdata),
            train = test_x,
            target = test_y,
            event_level = "second",
            num_features = num_features
          )$data %>%
            # Keep only positive importance values
            dplyr::filter(Importance > 0)
        }
      )
    }

    # Store the VIP data if Importance > 0
    vip_list[[name]] <- vip_data
    cli::cli_alert_success("Variable importance computed!")
  }

  return(vip_list)
}



#' @title ensembleBP
#' @description This function performs the ensemble of models using the ensBP approach.
#'
#' @param preProcess.obj An object of class preProcess.
#' @param models A character vector specifying the models to be used.
#'               Possible values are 'xgboost', 'bag_tree', 'lightGBM', 'pls', 'logistic',
#'               'C5_rules', 'mars', 'bag_mars', 'mlp', 'bag_mlp', 'decision_tree',
#'               'rand_forest', 'svm_linear', 'svm_poly', 'svm_rbf'.
#' @param selector.recipes A character vector specifying the selector recipes to be used.
#'                         Possible values are 'boruta', 'roc', 'infgain', 'mrmr'.
#' @param tuning.method A character string specifying the tuning method to be used.
#'                      Possible values are 'tune_grid (default)', tune_race_anova',
#'                      tune_race_win_loss', 'tune_bayes', 'tune_sim_anneal'.
#' @param n An integer specifying the number of iterations for the tuning method.
#' @param v An integer specifying the number of folds for the cross-validation during the hyperparameters tuning.
#' @param metric A character string specifying the metric to be used for tuning.
#'              Possible values are 'accuracy', 'roc_auc', 'sensitivity', 'specificity'.
#' @param nsim An integer specifying the number of simulations for the permutation-based VIP.
#' @param seed An integer specifying the seed for reproducibility.
#'
#' @return An object of class ensBP.
#'
#' @details
#' The function performs the ensemble of models using the ensBP approach.
#' The function first filters the genes that are not annotated in the GO and KEGG databases.
#' Then, the function tunes and fits the models using the specified tuning method.
#' The function computes the variable importances for each model using the permutation-based VIP,
#' when the direct VIP computation fails.
#'
#' @import dplyr
#' @importFrom baguette var_imp
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success cli_alert_danger cli_abort
#' @importFrom future plan multisession sequential
#' @importFrom parallel detectCores
#' @importFrom recipes recipe step_nzv step_normalize
#' @importFrom colino step_select_boruta step_select_roc step_select_infgain step_select_mrmr
#' @importFrom purrr map2 map imap_dfr
#' @importFrom rsample vfold_cv training testing
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom workflows extract_fit_parsnip
#' @importFrom vip vip
#' @importFrom kernlab predict
#' @importFrom withr with_seed
#' @importFrom tune select_best finalize_workflow last_fit
#' @importFrom workflowsets workflow_set extract_workflow collect_metrics rank_results workflow_map
#' @importFrom parsnip boost_tree bag_tree pls logistic_reg C5_rules mars bag_mars mlp bag_mlp decision_tree rand_forest svm_linear svm_poly svm_rbf
#'
#' @examples
#' /dontrun{
#' enso <- ensembleBP(preProcess.obj, models = c("bag_mlp", "rand_forest", "svm_poly"),
#'                    selector.recipes = "boruta", tuning.method = "tune_grid", n = 5,
#'                    v = 3, metric = "accuracy", nsim = 2, seed = 123)}
#'
#' @export
ensembleBP <- function(
    preProcess.obj, models = c("bag_mlp", "rand_forest", "svm_poly"), selector.recipes = "boruta",
    tuning.method = "tune_grid", n = 5, v = 3, metric = "accuracy",
    nsim = 2, seed = 123) {
  future::plan(future::multisession, workers = parallel::detectCores() - 1)
  cli::cli_h1("ensembleBP")

  available_models <- c(
    "xgboost", "bag_tree", "lightGBM", "pls", "logistic",
    "C5_rules", "mars", "bag_mars", "mlp", "bag_mlp",
    "decision_tree", "rand_forest", "svm_linear", "svm_poly",
    "svm_rbf"
  )
  available_selectors <- c("boruta", "roc", "infgain", "mrmr")

  withr::with_seed(
    seed = seed,
    code = {

      ## KEEP ONLY GENES ANNOTATED IN GO AND KEGG DATABASE
      cli::cli_alert_info("Filtering genes non-annotated in GO and KEGG db ...")
      # Upload annotated genes database
      data(annotated.genes)
      # Compute the number of columns before filtering. class col is excluded from computation
      ncol_pre <- ncol(preProcess.obj@splits$adjusted.split$data) - 1
      # Save class column to cbind it after filtering
      class <- preProcess.obj@splits$adjusted.split$data$class
      # Filter genes that are not present in the DB
      preProcess.obj@splits$adjusted.split$data <-
        preProcess.obj@splits$adjusted.split$data[, colnames(preProcess.obj@splits$adjusted.split$data) %in%
                                                    annotated.genes$genes_to_retrieve]
      # cbind class column that has been filtered out
      preProcess.obj@splits$adjusted.split$data <- cbind(preProcess.obj@splits$adjusted.split$data, class)
      # Compute the number of columns after filter. class col is excluded from computation
      ncol_post <- ncol(preProcess.obj@splits$adjusted.split$data) - 1
      cli::cli_alert_success("Filtered {ncol_pre - ncol_post} genes. Total number of genes is now: {ncol_post}!")

      train_test_split <- preProcess.obj@splits$adjusted.split
      train <- rsample::training(train_test_split)
      test <- rsample::testing(train_test_split)

      train_resamples <- rsample::vfold_cv(train, v = v, strata = class)


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
        train = train,
        train_test_split = train_test_split
      )

      cli::cli_h2("Variable Importances")
      vip_results <- calculate_vip(
        last_fit_results = t_and_f_output[[2]],
        test_x = test[, -ncol(test)],
        test_y = test$class,
        n_sim = nsim
      )
    }
  )

  obj <- t_and_f_output[[1]]
  obj@model.features <- vip_results
  future::plan(future::sequential)

  cli::cli_alert_success("Successfully executed ensembleBP function!")

  return(obj)
}
