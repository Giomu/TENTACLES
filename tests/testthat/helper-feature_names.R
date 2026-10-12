# Non-syntactic feature names found in real data: gene symbols, metabolites,
# lipids and module names. "HLA-B"/"HLA.B" and "C18:1"/"C18.1" collide once
# made syntactic. The beta is written as an escape to keep the file ASCII.
tricky_names <- c(
  HLA.B = "HLA-B",
  A2M   = "HLA.B",
  AARS  = "3-hydroxybutyrate",
  ABHD2 = "C18:1",
  ABCB1 = "C18.1",
  ABAT  = "PC(34:1)",
  ABCF3 = "T cell activation (I) (M7.0)",
  A1BG  = "\u03b2-alanine"
)

# Builds a small preProcess object whose features include tricky_names.
# The renamed genes are well expressed, so they survive the CPM filter.
make_tricky_pp <- function() {
  test_env <- load_test_data("acc.count", "acc.clin", package = "TENTACLES")
  acc.count <- test_env$acc.count
  acc.clin <- test_env$acc.clin

  acc.count <- acc.count[, c(1:100, which(colnames(acc.count) == "HLA.B"))]
  idx <- match(names(tricky_names), colnames(acc.count))
  colnames(acc.count)[idx] <- tricky_names

  preProcess(acc.count, acc.clin, is.normalized = FALSE, plot = FALSE)
}

# Renames feature columns of a preProcess object directly, bypassing
# preProcess, to simulate invalid names reaching runClassifiers.
rename_pp_features <- function(pp, from, to) {
  cols <- colnames(pp@processed$adjusted.data)
  cols[match(from, cols)] <- to
  colnames(pp@processed$adjusted.data) <- cols
  pp
}
