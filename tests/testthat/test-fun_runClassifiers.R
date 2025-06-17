test_that("runClassifiers results are consistent with snapshot", {
  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]

  # Run preProcess
  pp <- preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )


  withr::with_options(
    new = list(future.globals.maxSize = 2 * 1024^3),
    code = {
      test_results <- runClassifiers(pp,
        models = c("bag_mlp", "rand_forest"),
        selector.recipes = c("boruta", "roc"),
        filter = TRUE, downsample = TRUE, plot = FALSE
      )
    }
  )

  expect_snapshot(test_results)
})

test_that("runClassifiers results with all options enabled are consistent with snapshot", {
  # Skip heavy tests on CRAN
  skip_on_cran()
  # Skip heavy tests if RUN_HEAVY_TESTS is not set to true
  skip_if(Sys.getenv("RUN_HEAVY_TESTS") != "true")

  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:1500]

  # Run preProcess
  pp <- preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )

  all_models <- c(
    "xgboost", "bag_tree", "lightGBM",
    "pls", "logistic", "C5_rules", "mars", "bag_mars", "mlp",
    "bag_mlp", "decision_tree", "rand_forest", "svm_linear",
    "svm_poly", "svm_rbf"
  )

  all_selectors <- rep(c("boruta", "roc", "infgain", "mrmr", "corr"), each = 3)

  withr::with_options(
    new = list(future.globals.maxSize = 2 * 1024^3),
    code = {
      test_results <- runClassifiers(pp,
        models = all_models,
        selector.recipes = all_selectors,
        filter = TRUE, downsample = TRUE, plot = FALSE
      )
    }
  )

  expect_snapshot(test_results)
})
