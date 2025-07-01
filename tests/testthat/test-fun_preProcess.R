# test-preProcess.R

# --------------------------------------------------------#
#                  S4 class preProcess                    #
# --------------------------------------------------------#

test_that("preProcess.obj can be instantiated correctly", {
  # Dummy data
  df.count <- as.data.frame(matrix(1:12, nrow = 3))
  df.meta <- data.frame(class = factor(c(0, 1, 1)))
  rownames(df.count) <- rownames(df.meta) <- paste0("S", 1:3)
  data.info <- list(type = "rnaseq", normalized = FALSE)

  # Create the object
  obj <- new("preProcess.obj",
             raw = df.count,
             processed = list(),
             metadata = df.meta,
             data.info = data.info
  )

  expect_s4_class(obj, "preProcess.obj")
  expect_equal(dim(obj@raw), c(3, 4))
  expect_equal(slotNames(obj), c("raw", "processed", "metadata", "data.info"))
  expect_equal(obj@data.info$type, "rnaseq")
})

# --------------------------------------------------------#
#                       data.check                        #
# --------------------------------------------------------#

test_that("data.check correctly validates inputs", {
  # Test: valid input
  res <- data.check("array", TRUE)
  expect_type(res, "list")
  expect_equal(res$type, "array")
  expect_true(res$normalized)

  # Test: valid input
  res2 <- data.check("array", FALSE)
  expect_equal(res2$type, "array")
  expect_false(res2$normalized)

  # Test: wrong inputs
  expect_error(data.check("invalid_type", TRUE), "Invalid data type")
  expect_error(data.check("rnaseq", "nope"))
  expect_error(data.check(123, TRUE), "must be a character string of length 1")
  expect_error(data.check(c("rnaseq", "array"), TRUE), "must be a character string of length 1")
  expect_error(data.check("rnaseq", "yes"), "single logical value")
  expect_error(data.check("rnaseq", c(TRUE, FALSE)), "single logical value")
})

# --------------------------------------------------------#
#                      match.samples                      #
# --------------------------------------------------------#

test_that("match.samples works and filters correctly", {
  df.count <- as.data.frame(matrix(1:20, nrow = 4))
  rownames(df.count) <- c("S1", "S2", "S3", "S4")

  df.meta <- data.frame(class = c(0, 1, 0, 1), age = c(10, 20, 30, 40))
  rownames(df.meta) <- c("S2", "S4", "S1", "S3")

  matched <- match.samples(df.count, df.meta)
  expect_type(matched, "character")
  expect_true(identical(matched, c("S1", "S2", "S3", "S4")))

  # No matches
  df.meta2 <- df.meta
  rownames(df.meta2) <- c("X1", "X2", "X3", "X4")
  expect_error(match.samples(df.count, df.meta2))

  # Some unmatched in meta
  df.meta3 <- df.meta
  rownames(df.meta3)[1] <- "missing_sample"
  matched_partial <- match.samples(df.count, df.meta3)
  expect_equal(matched_partial, c("S1", "S3", "S4"))

  # Some unmatched in counts
  df.count2 <- df.count[c(1,3), ]
  matched_partial <- match.samples(df.count2, df.meta)
  expect_equal(matched_partial, c("S1", "S3"))
})

# --------------------------------------------------------#
#                  data.import                            #
# --------------------------------------------------------#

test_that("data.import returns a valid preProcess.obj object", {
  df.count <- as.data.frame(matrix(rnbinom(40, size = 10, mu = 100), nrow = 4,
                                   dimnames = list(paste0("S", 1:4), paste0("G", 1:10))))
  df.clin <- data.frame(
    class = c(0, 1, 1, 0),
    age = c(20, 25, 30, 22),
    sex = c("M", "F", "F", "M")
  )
  rownames(df.clin) <- rownames(df.count)

  obj <- data.import(df.count, df.clin, class = "class", is.normalized = FALSE)
  expect_s4_class(obj, "preProcess.obj")
  expect_true(is.data.frame(obj@raw))
  expect_equal(dim(obj@raw), c(4, 11))  # 10 genes + 1 class
  expect_true(is.data.frame(obj@metadata))
  expect_type(obj@data.info, "list")
  expect_equal(levels(obj@metadata$class), c("0", "1"))
})

test_that("data.import handles renamed class and enforces binary labels", {
  # Class renaming with case.label
  df.count <- as.data.frame(matrix(rnbinom(40, size = 10, mu = 100), nrow = 4,
                                   dimnames = list(paste0("S", 1:4), paste0("G", 1:10))))
  df.clin1 <- data.frame(
    group = c("A", "B", "A", "B"),
    age = c(20, 25, 30, 22)
  )
  rownames(df.clin1) <- rownames(df.count)

  obj <- data.import(df.count, df.clin1, class = "group", case.label = "A", is.normalized = FALSE)
  expect_s4_class(obj, "preProcess.obj")
  expect_equal(levels(obj@metadata$group), c("0", "1"))

  # Error: class not binary
  df.clin2 <- data.frame(
    class = c("A", "B", "C", "D"),
    age = c(20, 25, 30, 22)
  )
  rownames(df.clin2) <- rownames(df.count)
  expect_error(data.import(df.count, df.clin2, class = "class"), "2 unique values")
})

test_that("data.import throws error with mismatched sample IDs or invalid input types", {
  df.count <- as.data.frame(matrix(rnbinom(40, size = 10, mu = 100), nrow = 4,
                                   dimnames = list(paste0("S", 1:4), paste0("G", 1:10))))
  df.clin <- data.frame(
    class = c(0, 1, 1, 0),
    age = c(20, 25, 30, 22)
  )

  # Sample ID mismatch
  rownames(df.clin) <- paste0("X", 1:4)
  expect_error(data.import(df.count, df.clin, class = "class"), "No common samples")

  # Count matrix not a data.frame
  df.count.mat <- matrix(rnbinom(40, size = 10, mu = 100), nrow = 4,
                         dimnames = list(paste0("S", 1:4), paste0("G", 1:10)))
  rownames(df.clin) <- paste0("S", 1:4)
  expect_error(data.import(df.count.mat, df.clin, class = "class"), "data.frame")
})


# --------------------------------------------------------#
#                validate_batch_args tests                #
# --------------------------------------------------------#

test_that("validate_batch_args accepts valid input", {
  meta <- data.frame(batch = c("A", "B", "A", "B"),
                     cov1 = 1:4,
                     cov2 = letters[1:4])
  expect_no_error(validate_batch_args(meta, batch = "batch", covar.mod = c("cov1", "cov2")))
  expect_no_error(validate_batch_args(meta, batch = "batch", covar.mod = NULL))
  expect_no_error(validate_batch_args(meta, batch = NULL, covar.mod = c("cov1")))
})

test_that("validate_batch_args catches batch-related issues", {
  meta1 <- data.frame(batch = c("A", "B", "A", "B"))
  meta2 <- data.frame(group = c("A", "B", "A", "B"))
  meta3 <- data.frame(batch = c("A", NA, "A", "B"))

  expect_error(validate_batch_args(meta1, batch = 123), "must be a character string")
  expect_error(validate_batch_args(meta1, batch = c("a", "b")), "must be a character string")
  expect_error(validate_batch_args(meta2, batch = "batch"), "not found in metadata")
  expect_error(validate_batch_args(meta3, batch = "batch"), "Missing values detected")
})

test_that("validate_batch_args catches covar.mod-related issues", {
  meta1 <- data.frame(c1 = 1:4)
  meta2 <- data.frame(c1 = c(1, 2, NA, 4))
  meta3 <- data.frame(batch = c("A", "B", "A", "B"), other = 1:4)

  expect_error(validate_batch_args(meta1, covar.mod = 123), "must be a character vector")
  expect_error(validate_batch_args(meta1, covar.mod = "c2"), "not found in metadata")
  expect_error(validate_batch_args(meta2, covar.mod = "c1"), "Missing values detected")
  expect_error(validate_batch_args(meta3, batch = "batch", covar.mod = c("batch", "other")),
               "must not be included in covar.mod")
})

# --------------------------------------------------------#
#                    normalization tests                  #
# --------------------------------------------------------#

test_that("normalization removes lowly expressed genes based on minfraction", {
  df <- data.frame(
    gene1 = c(100, 120, 110, 130, 115, 125),
    gene2 = c(0, 5, 3, 0, 0, 4),
    gene3 = c(0, 0, 0, 3, 0, 5),
    gene4 = c(0, 0, 0, 0, 0, 0),
    class = c(0, 1, 0, 1, 0, 1)
  )

  # Case 1: minfraction = 0.5 → keep gene1 e gene2
  norm1 <- normalization(df, class = "class", mincpm = 1, minfraction = 0.5)
  expect_equal(sort(colnames(norm1)), sort(c("gene1", "gene2", "class")))

  # Case 2: minfraction = 0.3 → keep gene1, gene2 e gene3
  norm2 <- normalization(df, class = "class", mincpm = 1, minfraction = 0.3)
  expect_equal(sort(colnames(norm2)), sort(c("gene1", "gene2", "gene3", "class")))

  # Case 3: minfraction = 0.9 → keep solo gene1
  norm3 <- normalization(df, class = "class", mincpm = 1, minfraction = 0.9)
  expect_equal(sort(colnames(norm3)), sort(c("gene1", "class")))
})

test_that("normalization removes lowly expressed genes based on mincpm", {
  df <- data.frame(
    gene1 = c(1, 1, 1, 1, 1, 1),     # very low CPM
    gene2 = c(5, 6, 7, 8, 9, 10),    # intermediate CPM
    gene3 = c(500, 600, 550, 520, 580, 610),  # very high CPM
    class = c(0, 1, 0, 1, 0, 1)
  )

  # Case 1: mincpm = 1 → keep all genes (1 is the minimum threshold)
  norm1 <- normalization(df, class = "class", mincpm = 1, minfraction = 1)
  expect_equal(sort(colnames(norm1)), sort(c("gene1", "gene2", "gene3", "class")))

  # Case 2: mincpm = 2000 → remove gene1, keep gene2 and gene3
  norm2 <- normalization(df, class = "class", mincpm = 2000, minfraction = 1)
  expect_equal(sort(colnames(norm2)), sort(c("gene2", "gene3", "class")))

  # Case 3: mincpm = 15000 → remove gene1 and gene2, keep only gene3
  norm3 <- normalization(df, class = "class", mincpm = 15000, minfraction = 1)
  expect_equal(sort(colnames(norm3)), sort(c("gene3", "class")))

  # Case 4: mincpm = 2000000 → remove all genes
  expect_error(normalization(df, class = "class", mincpm = 2000000, minfraction = 1),
               "No genes passed the expression filtering")
})

test_that("normalization returns expected log2 CPM values", {
  df <- data.frame(
    gene1 = c(1, 1, 1, 1, 1, 1),
    gene2 = c(5, 6, 7, 8, 9, 10),
    gene3 = c(500, 600, 550, 520, 580, 610),
    class = c(0, 1, 0, 1, 0, 1)
  )

  # Output da confrontare
  expected <- data.frame(
    gene2 = c(13.26881, 13.26881, 13.61666, 13.89015, 13.90253, 13.98175),
    gene3 = c(19.91252, 19.91252, 19.91248, 19.91242, 19.91242, 19.91240),
    class = c(0, 1, 0, 1, 0, 1),
    row.names = paste0("Sample", 1:6)
  )

  # Calcolo
  norm2 <- normalization(df, class = "class", mincpm = 2000, minfraction = 1)

  # Confronto con tolleranza
  expect_equal(norm2, expected, tolerance = 1e-5)
})














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

