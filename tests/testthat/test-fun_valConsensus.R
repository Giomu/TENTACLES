test_that("valConsensus aborts naming features missing from df.count", {
  # Without the check, the subset fails with "undefined columns selected".
  data <- make_tricky_pp()@processed$adjusted.data
  df.count <- data[, setdiff(colnames(data), "class")]

  expect_error(
    valConsensus(df.count, c("HLA-B", "C18:1", "NOT_A_FEATURE"), data$class, N = 2),
    regexp = "NOT_A_FEATURE"
  )
})
