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

test_that("preProcess results are reproducible after normalization and batch correction", {
  # Load test datasets into a temporary environment
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")

  # Access the data
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  # Keep 100 genes only
  acc.count <- acc.count[, 1:100]

  reference_results <- readRDS(test_path("fixtures", "preProcess_rep_test.rds"))
  test_results <- preProcess(acc.count, acc.clin,
    is.normalized = FALSE, batch = "patient.gender", plot = FALSE
  )

  # Compare the results
  expect_equal(unclass(reference_results@raw), unclass(test_results@raw), tolerance = 1e-5)
  expect_equal(unclass(reference_results@processed$normalized), unclass(test_results@processed$normalized), tolerance = 1e-5)
  expect_equal(unclass(reference_results@processed$sbatched), unclass(test_results@processed$sbatched), tolerance = 1e-5)
  expect_equal(unclass(reference_results@processed$adjusted.data), unclass(test_results@processed$adjusted.data), tolerance = 1e-5)
})
