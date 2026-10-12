#Consensus features with non-syntactic names, plus the expression table and
# class labels they are tested on.
make_tricky_consensus <- function() {
  data <- make_tricky_pp()@processed$adjusted.data
  list(
    df.count = data[, setdiff(colnames(data), "class")],
    class = data$class,
    gene.list = c("HLA-B", "HLA.B", "3-hydroxybutyrate", "C18:1", "PC(34:1)")
  )
}

test_that("testConsensus returns MLP importances with the original feature names", {
  # nnet's importances come back escaped (`HLA-B`), so they no longer match
  # the consensus list.
  cons <- make_tricky_consensus()

  tc <- testConsensus(cons$df.count, cons$gene.list, cons$class,
                      top_loadings = 3, plot = FALSE)

  variables <- tc@mlp$importances$Variable
  expect_gt(length(variables), 0)
  expect_true(all(variables %in% cons$gene.list))
})

test_that("testConsensus warns about features missing from df.count", {
  cons <- make_tricky_consensus()

  expect_warning(
    testConsensus(cons$df.count, c(cons$gene.list, "NOT_A_FEATURE"), cons$class,
                  top_loadings = 3, plot = FALSE),
    regexp = "NOT_A_FEATURE"
  )
})
