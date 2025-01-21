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

# test_that("Support to different labels formats and error to non-binary labels"){

# }
