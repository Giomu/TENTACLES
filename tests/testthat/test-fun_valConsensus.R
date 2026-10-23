test_that("valConsensus aborts naming features missing from df.count", {
  # Without the check, the subset fails with "undefined columns selected".
  data <- make_tricky_pp()@processed$adjusted.data
  df.count <- data[, setdiff(colnames(data), "class")]

  expect_error(
    valConsensus(df.count, c("HLA-B", "C18:1", "NOT_A_FEATURE"), data$class, N = 2),
    regexp = "NOT_A_FEATURE"
  )
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

# Two combinations scored by two clustering methods, as selectTopCombinations() returns
make_top_results <- function() {
  results <- list(
    GeneA_GeneB = list(
      KMeans = list(Accuracy = 0.85, Precision = 0.82, Recall = 0.80, FScore = 0.81),
      GMM    = list(Accuracy = 0.86, Precision = 0.83, Recall = 0.82, FScore = 0.82)
    ),
    GeneC_GeneD = list(
      KMeans = list(Accuracy = 0.10, Precision = 0.15, Recall = 0.20, FScore = 0.17),
      GMM    = list(Accuracy = 0.98, Precision = NA, Recall = 0.97, FScore = 0.95)
    )
  )
  selectTopCombinations(results, N = 2, metric = "Accuracy")
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
