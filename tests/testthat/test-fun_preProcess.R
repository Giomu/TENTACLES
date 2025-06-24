test_that("Importing with different class column name works", {
  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]

  # Modify "class" column to "labels"
  acc.clin$labels <- acc.clin$class
  acc.clin$class <- NULL

  # Test the functionality
  expect_no_error(preProcess(acc.count, acc.clin,
    class = "labels", is.normalized = FALSE,
    batch = "patient.gender", plot = FALSE
  ))
})

test_that("Support to different label formats and error to non-binary labels", {
  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]
  class_col <- acc.clin$class

  # Test the functionality
  # Use 0 and 1 as labels (default labels)
  expect_no_error(preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))

  # Change the labels 0 and 1 to "A" and "B"
  acc.clin$class <- ifelse(acc.clin$class == 0, "A", "B")
  expect_no_error(preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))

  # Add a number to the labels
  acc.clin$class <- class_col
  acc.clin$class <- acc.clin$class + 2
  expect_no_error(preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))

  # Create non-binary labels and expect an error
  acc.clin$class <- class_col
  acc.clin$class <- seq(1, nrow(acc.clin))
  expect_error(preProcess(acc.count, acc.clin,
    class = "class_non_binary", is.normalized = FALSE,
    batch = "patient.gender", plot = FALSE
  ))
})

test_that("Throw an error if all rows (samples) do not match between tables", {
  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]

  # Transpose the count data (samples will be columns)
  acc.count.2 <- t(acc.count)
  expect_error(preProcess(acc.count.2, acc.clin,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))

  # Transpose the clinical data (samples will be columns)
  acc.clin.3 <- t(acc.clin)
  expect_error(preProcess(acc.count, acc.clin.3,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))
})

test_that("No error when just some samples do not match between tables", {
  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]

  # Remove some samples from the count data
  acc.count.2 <- acc.count[-(1:10), ]
  expect_no_error(preProcess(acc.count.2, acc.clin,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))

  # Remove some samples from the clinical data
  acc.clin.2 <- acc.clin[-(1:10), ]
  expect_no_error(preProcess(acc.count, acc.clin.2,
    is.normalized = FALSE, batch = "patient.gender",
    plot = FALSE
  ))
})


test_that("preProcess is deterministic within the same OS", {
  # This test checks that preProcess() produces the same result
  # when run twice under identical conditions, OS, and R environment.
  # This guards against hidden randomness in the preprocessing pipeline.
  #
  # Note: Results may still differ across operating systems due to
  # underlying libraries or data.table/BLAS differences.

  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]

  # Run preProcess with different configurations, twice each
  res_batch_1 <- preProcess(acc.count, acc.clin,
                            is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )
  res_batch_2 <- preProcess(acc.count, acc.clin,
                            is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )

  res_covariate_1 <- preProcess(acc.count, acc.clin,
                                is.normalized = FALSE, batch = "patient.gender", plot = FALSE,
                                covar.mod = "patient.primary_pathology.laterality"
  )
  res_covariate_2 <- preProcess(acc.count, acc.clin,
                                is.normalized = FALSE, batch = "patient.gender", plot = FALSE,
                                covar.mod = "patient.primary_pathology.laterality"
  )

  res_no_batch_1 <- preProcess(acc.count, acc.clin,
                               is.normalized = FALSE, plot = FALSE
  )
  res_no_batch_2 <- preProcess(acc.count, acc.clin,
                               is.normalized = FALSE, plot = FALSE
  )

  # Assert that the results are strictly identical for each configuration
  expect_equal(res_batch_1, res_batch_2)
  expect_equal(res_covariate_1, res_covariate_2)
  expect_equal(res_no_batch_1, res_no_batch_2)
  # Note: Differences across OSs are expected, but not within a single environment.
})

