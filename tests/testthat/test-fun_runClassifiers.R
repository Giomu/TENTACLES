test_that("runClassifiers is deterministic within the same OS", {
  # This test ensures that runClassifiers() produces identical results
  # when run multiple times under the same conditions, OS, and R environment.
  # This is important for reproducibility: any hidden or uncontrolled sources
  # of randomness would break this expectation.
  #
  # Note: Due to numeric library and platform differences, it is expected that
  # results may differ between operating systems (e.g., macOS vs. Windows/Linux).
  # However, repeated runs on the *same* OS with the same seed and package versions
  # must produce identical outputs.

  # Load test datasets into a temporary environment for reproducibility
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the count matrix and clinical data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep only the first 100 genes to reduce test runtime
  acc.count <- acc.count[, 1:100]

  # Preprocess the data (including normalization and batch correction)
  pp <- preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )

  # Ensure a high memory limit for future globals (for compatibility)
  withr::with_options(
    new = list(future.globals.maxSize = 2 * 1024^3),
    code = {
      # Run runClassifiers() twice in the exact same environment and configuration
      # This should produce identical results on the same OS, seed, and package versions
      res1 <- runClassifiers(pp,
        models = c("bag_mlp", "rand_forest"),
        selector.recipes = c("boruta", "roc"),
        filter = TRUE, downsample = TRUE, plot = FALSE,
        parallel = FALSE
      )
      res2 <- runClassifiers(pp,
        models = c("bag_mlp", "rand_forest"),
        selector.recipes = c("boruta", "roc"),
        filter = TRUE, downsample = TRUE, plot = FALSE,
        parallel = FALSE
      )
    }
  )

  # Assert that both results are strictly identical (deterministic within OS)
  expect_equal(res1, res2)
})


test_that("runClassifiers object is created even if a model fails to tune", {
  # This test verifies that the `runClassifiers` function successfully creates a results object
  # even if one of the specified models fails to tune.
  # It ensures that the function is robust to model tuning failures
  # and still returns a valid results object for further analysis.

  # Note: For this particular test, the random forest model is expected to fail
  # to tune due to the small number of cases in the clinical data.

  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep a few genes
  acc.count <- acc.count[, 1:100]

  # Break the response variable so only random forest will fail to tune
  acc.clin$class <- factor(c(rep("control", nrow(acc.clin) - 5), rep("case", 5))) # 5 cases, rest controls

  pp <- preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )

  withr::with_options(
    new = list(future.globals.maxSize = 2 * 1024^3),
    code = {
      suppressWarnings(
        test_results <- runClassifiers(pp,
          models = c("bag_mlp", "rand_forest"),
          selector.recipes = c("boruta"),
          filter = TRUE, downsample = TRUE, plot = FALSE
        )
      )
    }
  )

  expect_true(exists("test_results"))
  expect_false(is.null(test_results))
})

test_that("runClassifiers is deterministic with all options enabled (heavy test)", {
  # This test ensures that runClassifiers() returns identical results
  # when run twice with all models and selectors enabled, under the same
  # OS, seed, and R environment. This checks for hidden randomness or
  # non-determinism in complex pipelines.
  #
  # Note: Results may still differ across operating systems, but must be
  # identical on repeated runs within the same OS and environment.

  # Skip heavy tests on CRAN and unless explicitly requested
  skip_on_cran()
  skip_if(Sys.getenv("RUN_HEAVY_TESTS") != "true")

  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 1500 genes for a realistic heavy test
  acc.count <- acc.count[, 1:1500]

  # Preprocess the data
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
      # Run twice under identical conditions
      res1 <- runClassifiers(pp,
        models = all_models,
        selector.recipes = all_selectors,
        filter = TRUE, downsample = TRUE, plot = FALSE,
        parallel = FALSE
      )
      res2 <- runClassifiers(pp,
        models = all_models,
        selector.recipes = all_selectors,
        filter = TRUE, downsample = TRUE, plot = FALSE,
        parallel = FALSE
      )
    }
  )

  # Assert strict equality (deterministic within OS)
  expect_equal(res1, res2)
})
