#' DataFrame containing genes annotated in GO and KEGG databases
#'
#' This dataframe contains all genes annotated in GO and KEGG databases and will be used
#' to filter non-annotated genes with a double goal:
#' 1. Speed up the entire classification process.
#' 2. Remove genes that are not annotated and will not contribute to an enrichment.
#'
#' @format A dataframe with 10760 rows and 2 columns.
#' @usage data(annotated.genes)
"annotated.genes"
