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
  # Suggested packages needed by some models and selectors
  for (pkg in c("xgboost", "bonsai", "lightgbm", "plsmod", "mixOmics", "earth",
                "FSelectorRcpp", "praznik")) {
    skip_if_not_installed(pkg)
  }

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

test_that("make_tricky_pp keeps the non-syntactic names through preProcess", {
  pp <- make_tricky_pp()
  expect_true(all(tricky_names %in% colnames(pp@processed$adjusted.data)))
})

test_that("runClassifiers fits every workflow when feature names are non-syntactic", {
  skip_if_not_installed("FSelectorRcpp")
  skip_if_not_installed("earth")
  # Formula-based components (colino's infgain selector, baguette's bag_tree
  # and bag_mars) parse "HLA-B" as "HLA - B" and fail during tuning.
  pp <- make_tricky_pp()

  rc <- runClassifiers(pp,
                       models = c("decision_tree", "bag_tree", "bag_mars", "mlp"),
                       selector.recipes = c("infgain", "base", "base", "base"),
                       n = 2, v = 3, nsim = 1,
                       filter = FALSE, plot = FALSE, parallel = FALSE
  )

  expect_setequal(
    names(rc@models.info),
    c("decision_tree_infgain", "bag_tree_base", "bag_mars_base", "mlp_base")
  )
})

test_that("runClassifiers fits engines that may reject non-ASCII column names", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("bonsai")
  skip_if_not_installed("lightgbm")
  pp <- make_tricky_pp()

  rc <- runClassifiers(pp,
                       models = c("xgboost", "C5_rules", "lightGBM"),
                       selector.recipes = c("base", "base", "base"),
                       n = 2, v = 3, nsim = 1,
                       filter = FALSE, plot = FALSE, parallel = FALSE
  )

  expect_setequal(
    names(rc@models.info),
    c("xgboost_base", "C5_rules_base", "lightGBM_base")
  )
})

test_that("runClassifiers returns variable importances with the original feature names", {
  # Some engines escape non-syntactic names in their importance output
  # (e.g. mlp returns `HLA-B` with backticks), so features no longer match
  # the input column names.
  pp <- make_tricky_pp()

  rc <- runClassifiers(pp,
                       models = c("mlp", "rand_forest"),
                       selector.recipes = c("base", "base"),
                       n = 2, v = 3, nsim = 1,
                       filter = FALSE, plot = FALSE, parallel = FALSE
  )

  all_features <- unique(unlist(lapply(rc@model.features, function(x) x$Variable)))
  expect_true(all(all_features %in% colnames(rc@data$adjusted.data)))
  expect_true("HLA-B" %in% rc@model.features$mlp_base$Variable)
})

test_that("runClassifiers stores a one-to-one map from internal to original names", {
  pp <- make_tricky_pp()

  rc <- runClassifiers(pp,
                       models = "rand_forest",
                       selector.recipes = "base",
                       n = 2, v = 3, nsim = 1,
                       filter = FALSE, plot = FALSE, parallel = FALSE
  )

  map <- rc@data$feature.map
  features <- setdiff(colnames(pp@processed$adjusted.data), "class")
  expect_type(map, "character")
  expect_setequal(unname(map), features)
  expect_false(anyDuplicated(names(map)) > 0)
  expect_equal(make.names(names(map)), names(map))
  # Internal names must be ASCII: C5.0 cannot parse its names file otherwise
  # (C5imp() fails with "cannot parse names file" for "\u03b2.alanine").
  expect_false(any(grepl("[^ -~]", names(map))))
  # Names that are already syntactic are kept as they are
  expect_equal(unname(map[c("HLA.B", "C18.1", "A2ML1")]), c("HLA.B", "C18.1", "A2ML1"))
})

test_that("runClassifiers aborts on missing or empty feature names", {
  pp <- make_tricky_pp()

  expect_error(
    runClassifiers(rename_pp_features(pp, "A2ML1", NA),
                   models = "rand_forest", selector.recipes = "base",
                   filter = FALSE, plot = FALSE, parallel = FALSE),
    class = "tentacles_invalid_feature_names"
  )
  expect_error(
    runClassifiers(rename_pp_features(pp, "A2ML1", ""),
                   models = "rand_forest", selector.recipes = "base",
                   filter = FALSE, plot = FALSE, parallel = FALSE),
    class = "tentacles_invalid_feature_names"
  )
})

test_that("runClassifiers aborts on duplicated feature names, naming them", {
  pp <- rename_pp_features(make_tricky_pp(), "A2ML1", "HLA-B")

  expect_error(
    runClassifiers(pp,
                   models = "rand_forest", selector.recipes = "base",
                   filter = FALSE, plot = FALSE, parallel = FALSE),
    regexp = "HLA-B",
    class = "tentacles_invalid_feature_names"
  )
})

test_that("runClassifiers aborts when a feature is named 'class'", {
  # It would collide with the outcome column.
  pp <- rename_pp_features(make_tricky_pp(), "A2ML1", "class")

  expect_error(
    runClassifiers(pp,
                   models = "rand_forest", selector.recipes = "base",
                   filter = FALSE, plot = FALSE, parallel = FALSE),
    class = "tentacles_invalid_feature_names"
  )
})

test_that("required_model_pkgs lists engines, parsnip extensions and selector backends", {
  expect_setequal(
    required_model_pkgs(c("lightGBM", "C5_rules", "pls"), c("base", "roc", "corr")),
    c("bonsai", "lightgbm", "rules", "C50", "plsmod", "mixOmics")
  )
  expect_setequal(
    required_model_pkgs("logistic", c("infgain", "mrmr", "boruta")),
    c("FSelectorRcpp", "praznik", "Boruta")
  )
  expect_null(required_model_pkgs("logistic", "base"))
})

test_that("check_required_pkgs names the model or selector behind each missing package", {
  local_mocked_bindings(is_installed = function(pkg) !pkg %in% c("plsmod", "mixOmics", "FSelectorRcpp"))

  err <- expect_error(
    check_required_pkgs(c("pls", "rand_forest"), c("infgain", "base")),
    class = "tentacles_missing_packages"
  )
  msg <- conditionMessage(err)
  expect_match(msg, 'Model "pls" needs plsmod and mixOmics', fixed = TRUE)
  expect_match(msg, 'Selector "infgain" needs FSelectorRcpp', fixed = TRUE)
  expect_no_match(msg, "rand_forest", fixed = TRUE)
})

test_that("check_required_pkgs gives install commands for CRAN and Bioconductor", {
  local_mocked_bindings(is_installed = function(pkg) !pkg %in% c("plsmod", "mixOmics", "earth"))

  err <- expect_error(check_required_pkgs(c("pls", "mars"), "base"),
                      class = "tentacles_missing_packages")
  msg <- conditionMessage(err)
  expect_match(msg, 'install.packages(c("plsmod", "earth"))', fixed = TRUE)
  expect_match(msg, 'BiocManager::install("mixOmics")', fixed = TRUE)
})

test_that("check_required_pkgs adds BiocManager when a Bioconductor package needs it", {
  local_mocked_bindings(is_installed = function(pkg) !pkg %in% c("mixOmics", "BiocManager"))

  err <- expect_error(check_required_pkgs("pls", "base"),
                      class = "tentacles_missing_packages")
  expect_match(conditionMessage(err), 'install.packages("BiocManager")', fixed = TRUE)
})

test_that("check_required_pkgs passes when every package is installed", {
  local_mocked_bindings(is_installed = function(pkg) TRUE)

  expect_invisible(check_required_pkgs(c("pls", "mars"), "infgain"))
})
