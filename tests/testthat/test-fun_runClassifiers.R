test_that("runClassifiers results are reproducible", {
  # Load preProcess object
  pp <- readRDS(test_path("fixtures", "preProcess_rep_test.rds"))

  # Load saved runClassifiers object
  reference_results <- readRDS(test_path("fixtures", "runClassifiers_rep_test.rds"))

  # Run classifiers
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

  # Compare results
  expect_identical(unclass(reference_results@models.info[[1]]), unclass(test_results@models.info[[1]]))
  expect_identical(unclass(reference_results@models.info[[2]]), unclass(test_results@models.info[[2]]))
  expect_identical(unclass(reference_results@model.features[[1]]), unclass(test_results@model.features[[1]]))
  expect_identical(unclass(reference_results@model.features[[2]]), unclass(test_results@model.features[[2]]))
  expect_identical(unclass(reference_results@performances$tuning_metric), unclass(test_results@performances$tuning_metric))
  expect_identical(unclass(reference_results@performances$final_metrics), unclass(test_results@performances$final_metrics))
})
