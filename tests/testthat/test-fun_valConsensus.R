test_that("valConsensus aborts naming features missing from df.count", {
  # Without the check, the subset fails with "undefined columns selected".
  data <- make_tricky_pp()@processed$adjusted.data
  df.count <- data[, setdiff(colnames(data), "class")]

  expect_error(
    valConsensus(df.count, c("HLA-B", "C18:1", "NOT_A_FEATURE"), data$class, N = 2),
    regexp = "NOT_A_FEATURE"
  )
})

test_that("valConsensus aborts when class is not coded 0/1", {
  data <- make_tricky_pp()@processed$adjusted.data
  genes <- c("A1CF", "A2BP1")
  class <- as.numeric(as.character(data$class)) + 1

  expect_error(valConsensus(data[, genes], genes, class, N = 1), "0 \\(control\\)")
})

test_that("valConsensus aborts when class does not match the samples", {
  data <- make_tricky_pp()@processed$adjusted.data
  genes <- c("A1CF", "A2BP1")

  expect_error(valConsensus(data[, genes], genes, data$class[-1], N = 1), "samples")
  expect_error(valConsensus(data[, genes], genes, rep(0, nrow(data)), N = 1), "both classes")
})

test_that("evaluate_one_side handles a method that puts all samples in one cluster", {
  truth <- c(1, 0, 1, 0, 1, 1)
  pred <- rep(0, 6)

  # Predicting no 1s leaves FScore undefined, so the inverted assignment is kept
  expect_no_warning(res <- evaluate_one_side(pred, truth, "FScore"))
  expect_equal(res$clusters, rep(1L, 6))
  expect_equal(res$Accuracy, 4 / 6)
  expect_equal(res$Precision, 4 / 6)
  expect_equal(res$Recall, 1)
  expect_equal(res$FScore, 0.8)

  # Accuracy is defined for both; with mostly controls, all 0 wins and Precision is NA
  truth <- c(0, 0, 0, 0, 1, 1)
  expect_no_warning(res <- evaluate_one_side(pred, truth, "Accuracy"))
  expect_equal(res$clusters, rep(0L, 6))
  expect_equal(res$Accuracy, 4 / 6)
  expect_equal(res$Recall, 0)
  expect_true(is.na(res$Precision))
  expect_true(is.na(res$FScore))
})

test_that("gene_combination_matrix marks genes per combination, most frequent first", {
  combos <- list(c("B", "GENE_X"), c("A", "GENE_X"), c("A", "B", "GENE_X"))
  m <- gene_combination_matrix(combos)

  expect_equal(colnames(m), c("C1", "C2", "C3"))
  # GENE_X is in all three; A and B tie at two and are sorted by name
  expect_equal(rownames(m), c("GENE_X", "A", "B"))
  expect_equal(unname(m["A", ]), c(0L, 1L, 1L))
  expect_equal(unname(m["B", ]), c(1L, 0L, 1L))
})

test_that("gene_frequency keeps the binary columns and adds counts and combinations", {
  m <- gene_combination_matrix(list(c("A", "B"), c("A", "C")))
  freq <- gene_frequency(m)

  expect_equal(names(freq), c("Gene", "C1", "C2", "Count", "Combinations"))
  expect_equal(freq$Gene, c("A", "B", "C"))
  expect_equal(freq$C1, c(1L, 1L, 0L))
  expect_equal(freq$Count, c(2L, 1L, 1L))
  expect_equal(freq$Combinations, c("C1, C2", "C1", "C2"))
})

test_that("plotGeneCombinations returns a ggplot object", {
  freq <- gene_frequency(gene_combination_matrix(list(c("A", "B"), c("A", "C"))))

  expect_s3_class(plotGeneCombinations(freq), "ggplot")
  expect_error(plotGeneCombinations(data.frame(A = 1)), "Gene")
  expect_error(plotGeneCombinations(as.matrix(freq)), "data frame")
})

test_that("plotGeneCombinations shows how many combinations contain each gene", {
  freq <- gene_frequency(gene_combination_matrix(list(c("A", "B"), c("A", "C"))))
  p <- plotGeneCombinations(freq)

  counts <- p$layers[[2]]$data
  expect_equal(counts$Gene, c("A", "B", "C"))
  expect_equal(counts$Label, c("2/2", "1/2", "1/2"))
})

test_that("valConsensus returns gene frequencies of the top combinations", {
  data <- make_tricky_pp()@processed$adjusted.data
  genes <- c("A1CF", "A2BP1", "A2LD1")
  df.count <- data[, genes]
  # Genes are taken from the combinations, not split on "_"
  colnames(df.count)[3] <- genes[3] <- "GENE_X"

  withr::local_pdf(tempfile(fileext = ".pdf"))
  vc <- suppressMessages(valConsensus(df.count, genes, data$class, N = 4))

  # N = 4 keeps all combinations: three pairs and the triple, so each gene is in three
  expect_setequal(vc$geneFrequency$Gene, genes)
  expect_equal(vc$geneFrequency$Count, c(3L, 3L, 3L))
  expect_equal(names(vc$geneFrequency), c("Gene", paste0("C", 1:4), "Count", "Combinations"))
  expect_equal(sort(unname(colSums(vc$geneFrequency[paste0("C", 1:4)]))), c(2, 2, 2, 3))
  expect_null(vc$geneCombinations)
})

test_that("valConsensus warns and caps N at the number of combinations", {
  data <- make_tricky_pp()@processed$adjusted.data
  genes <- c("A1CF", "A2BP1")

  withr::local_pdf(tempfile(fileext = ".pdf"))
  # Two genes give a single combination
  expect_message(
    vc <- valConsensus(data[, genes], genes, data$class, N = 5),
    "larger than the 1 possible combination. Using"
  )
  expect_equal(nrow(vc$topCombinations), 1)
})

test_that("valConsensus keeps the clusters and metrics of all combinations", {
  data <- make_tricky_pp()@processed$adjusted.data
  genes <- c("A1CF", "A2BP1", "A2LD1")

  withr::local_pdf(tempfile(fileext = ".pdf"))
  vc <- suppressMessages(valConsensus(data[, genes], genes, data$class, N = 2, metric = "Accuracy"))
  res <- vc$allResults
  methods <- c("KMeans", "GMM", "HC", "PCA", "tSNE", "UMAP")
  combos <- c("A1CF_A2BP1", "A1CF_A2LD1", "A2BP1_A2LD1", "A1CF_A2BP1_A2LD1")

  expect_equal(unname(res$trueClass), as.integer(as.character(data$class)))
  expect_equal(names(res$trueClass), rownames(data))
  expect_equal(dim(res$clusters), c(nrow(data), 4L, 6L))
  expect_equal(dimnames(res$clusters)$combination, combos)
  expect_equal(dimnames(res$clusters)$method, methods)
  expect_type(res$clusters, "integer")
  expect_equal(res$metrics$Gene_Combination, combos)
  expect_equal(res$metrics$Mean_Metric, rowMeans(res$metrics[paste0(methods, "_Accuracy")]))

  # Every metric can be recomputed from trueClass and clusters
  truth <- factor(res$trueClass, levels = c(0, 1))
  for (combo in combos) {
    for (m in methods) {
      estimate <- factor(res$clusters[, combo, m], levels = c(0, 1))
      row <- res$metrics[res$metrics$Gene_Combination == combo, ]
      expect_equal(yardstick::accuracy_vec(truth, estimate, event_level = "second"), row[[paste0(m, "_Accuracy")]])
      expect_equal(yardstick::recall_vec(truth, estimate, event_level = "second"), row[[paste0(m, "_Recall")]])
      expect_equal(
        suppressWarnings(yardstick::f_meas_vec(truth, estimate, event_level = "second")),
        row[[paste0(m, "_FScore")]]
      )
    }
  }

  # topCombinations is the top of metrics
  top <- res$metrics[order(res$metrics$Mean_Metric, decreasing = TRUE)[1:2], ]
  rownames(top) <- NULL
  expect_equal(vc$topCombinations, top)
})

test_that("selectTopCombinations ranks by mean metric, keeping ties in order and NA last", {
  expect_equal(selectTopCombinations(c(0.5, NA, 0.9, 0.5), N = 4), c(3L, 1L, 4L, 2L))
  expect_equal(selectTopCombinations(c(0.5, NA, 0.9, 0.5), N = 2), c(3L, 1L))
})

test_that("valConsensus aborts on an unknown metric before testing combinations", {
  data <- make_tricky_pp()@processed$adjusted.data
  genes <- c("A1CF", "A2BP1")

  expect_error(valConsensus(data[, genes], genes, data$class, metric = "AUC"), "must be one of")
})

# Two combinations scored by two clustering methods, as in topCombinations
make_top_results <- function() {
  data.frame(
    Gene_Combination = c("GeneA_GeneB", "GeneC_GeneD"),
    Mean_Metric = c(0.855, 0.54),
    KMeans_Accuracy = c(0.85, 0.10), KMeans_Precision = c(0.82, 0.15),
    KMeans_Recall = c(0.80, 0.20), KMeans_FScore = c(0.81, 0.17),
    GMM_Accuracy = c(0.86, 0.98), GMM_Precision = c(0.83, NA),
    GMM_Recall = c(0.82, 0.97), GMM_FScore = c(0.82, 0.95)
  )
}

test_that("label_colors picks the text color with the higher contrast", {
  pal <- c("#EEE2DF", "#EED7C5", "#C89F9C", "#C97C5D", "#8B4B41")

  expect_equal(label_colors(c(0, 0.5, 1), pal, c(0, 1)), c("black", "black", "white"))
  expect_equal(label_colors(c(0, 1), c("black", "white"), c(0, 1)), c("white", "black"))
  # Missing and out-of-range values fall back to black
  expect_equal(label_colors(c(NA, 1.5), pal, c(0, 1)), c("black", "black"))
})

test_that("plotTopMetrics uses colors as the fill gradient", {
  p <- plotTopMetrics(make_top_results(), colors = c("white", "red"))
  fill_scale <- p$scales$get_scales("fill")

  expect_equal(toupper(fill_scale$palette(c(0, 1))), c("#FFFFFF", "#FF0000"))
  expect_equal(fill_scale$limits, c(0, 1))
  expect_error(plotTopMetrics(make_top_results(), colors = "red"), "at least two colors")
})

test_that("plotTopMetrics still accepts the deprecated low_color and high_color", {
  expect_warning(
    p <- plotTopMetrics(make_top_results(), low_color = "white", high_color = "steelblue"),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  # Only one given: the other end comes from colors
  expect_warning(plotTopMetrics(make_top_results(), high_color = "steelblue"), "deprecated")
})
