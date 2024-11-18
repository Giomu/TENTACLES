# Helper function to create binary dataframe from ensBP object
create_binary_df <- function(b) {
  # Check if b@model.features is a list
  if (!is.list(b@model.features)) {
    cli::cli_abort("ensBP.obj@model.features must be a list containing algorithm results.")
  }

  # Initialize a list to store all unique genes
  all_genes <- unique(unlist(lapply(b@model.features, function(x) x$Variable)))

  # Create a binary matrix
  binary_matrix <- matrix(0,
    nrow = length(b@model.features), ncol = length(all_genes),
    dimnames = list(names(b@model.features), all_genes)
  )

  # Populate the binary matrix
  for (alg in names(b@model.features)) {
    genes <- b@model.features[[alg]]$Variable
    binary_matrix[alg, genes] <- 1
  }

  # Convert the binary matrix to a data frame
  df_binary <- as.data.frame(binary_matrix)

  return(df_binary)
}


#' @title getConsensus
#' @description This function computes the consensus genes from an ensBP object.
#'
#' @param ensBP.obj An object of class ensBP.
#' @param n.min An integer specifying the minimum number of algorithms in which a gene must be present to be considered as a consensus gene.
#' @param group1 A character vector specifying the algorithms to be considered in the first group.
#' @param group2 A character vector specifying the algorithms to be considered in the second group.
#' @param meth1 A character string specifying the method to be used to compute the consensus genes in group1. Possible values are 'intersect' and 'union'.
#' @param meth2 A character string specifying the method to be used to compute the consensus genes in group2. Possible values are 'intersect' and 'union'.
#' @param meth.comb A character string specifying the method to be used to combine the consensus genes from group1 and group2. Possible values are 'intersect' and 'union'.
#' @param exclude A character vector specifying the algorithms to be excluded from the analysis.
#'
#' @return A list containing the consensus genes, the binary data frame,
#'         and the input parameters.
#'
#' @details
#' The function computes the consensus genes from an ensBP object.
#' The consensus genes are the genes that are present in at least n.min algorithms,
#' or the genes that are present in the algorithms specified in group1,
#' or the genes that are present in the algorithms specified in group1 and group2.
#' The consensus genes can be computed using the 'intersect' or 'union' method.
#'
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_danger cli_abort cli_warn
#'
#' @examples
#' /dontrun{
#' cons.1 <- getConsensus(ensBP.obj, n.min = 3)
#' cons.2 <- getConsensus(ensBP.obj, group1 = c("alg1", "alg2"), meth1 = "intersect")
#' cons.3 <- getConsensus(ensBP.obj, group1 = c("alg1", "alg2"), group2 = c("alg3", "alg4"),
#'                        meth1 = "intersect", meth2 = "intersect", meth.comb = "intersect")}
#'
#' @export
getConsensus <- function(ensBP.obj, n.min = NULL,
                         group1 = NULL, group2 = NULL,
                         meth1 = NULL, meth2 = NULL, meth.comb = NULL,
                         exclude = NULL) {
  # Store input parameters for reference
  inputParams <- list(
    n.min = n.min,
    group1 = group1, group2 = group2,
    meth1 = meth1, meth2 = meth2, meth.comb = meth.comb,
    exclude = exclude
  )

  # Preliminary print
  cli::cli_h1("Computing Consensus Genes")

  # Check class of input object. If binary data.frame simply assign it, else apply create_binary_df()
  if (class(ensBP.obj) == "data.frame") {
    # Assign ensBP.obj to df and print message
    df <- ensBP.obj
    cli::cli_alert_success("Informations successfully extracted from your data frame.")
  } else {
    # Print message to advise
    cli::cli_alert_info("Getting variables from ensBP.obj")
    df <- create_binary_df(ensBP.obj)
    cli::cli_alert_success("Variables information successfully extracted from ensBP.obj.")
  }

  # Preliminary check: verify if specified columns exist in the dataframe
  all_algorithms <- rownames(df)

  if (!is.null(exclude) && any(!exclude %in% rownames(df))) {
    cli::cli_abort("Some algorithms to exclude do not exist: {exclude[!exclude %in% rownames(df)]}.")
  }

  if (!is.null(group1) && any(!group1 %in% all_algorithms)) {
    cli::cli_abort("Some algorithms in 'group1' do not exist: {group1[!group1 %in% all_algorithms]}.")
  }

  if (!is.null(group2) && any(!group2 %in% all_algorithms)) {
    cli::cli_abort("Some algorithms in 'group2' do not exist: {group2[!group2 %in% all_algorithms]}.")
  }

  # Warn: n.min cannot be used with group1 or group2
  if (!is.null(n.min) && (!is.null(group1) || !is.null(group2))) {
    cli::cli_alert_info("You cannot use 'n.min' together with 'group1' or 'group2'. Selecting genes that appear in n.min algorithms")
  }

  # Remove specified columns from the dataframe
  if (!is.null(exclude)) {
    cli::cli_alert_info("Excluding the following algorithms: {paste(exclude, collapse = ', ')}")
    df <- df[!rownames(df) %in% exclude, ]
  }

  # Generate consensus based on n_min
  if (!is.null(n.min)) {
    selected_genes <- colnames(df)[colSums(df) >= n.min]
    if (length(selected_genes) == 0) {
      cli::cli_warn("No genes found in at least {n.min} algorithms.")
    } else {
      cli::cli_alert_success("{length(selected_genes)} genes found in at least {n.min} algorithms.")
    }
    consensus.list <- list(consensusGenes = selected_genes)
  }
  # Generate consensus using group1 only
  else if (!is.null(group1) & is.null(group2)) {
    subset1 <- df[group1, , drop = FALSE]

    if (meth1 == "intersect") {
      consensus.list <- list(colnames(subset1)[colSums(subset1) == nrow(subset1)])
      cli::cli_alert_success("Selecting genes using method 'intersect' for algorithms in group1.")
    } else if (meth1 == "union") {
      consensus.list <- list(colnames(subset1)[colSums(subset1) > 0])
      cli::cli_alert_success("Selecting genes using method 'union' for algorithms in group1.")
    }
  }
  # Generate consensus using both group1 and group2
  else if (!is.null(group1) & !is.null(group2)) {
    subset1 <- df[group1, , drop = FALSE]
    subset2 <- df[group2, , drop = FALSE]

    # Compute consensus for group1
    if (meth1 == "intersect") {
      consensus1 <- colnames(subset1)[colSums(subset1) == nrow(subset1)]
      cli::cli_alert_info("Found {length(consensus1)} genes with method 'intersect' in group1")
    } else if (meth1 == "union") {
      consensus1 <- colnames(subset1)[colSums(subset1) > 0]
      cli::cli_alert_info("Found {length(consensus1)} genes with method 'union' in group1")
    }

    # Compute consensus for group2
    if (meth2 == "intersect") {
      consensus2 <- colnames(subset2)[colSums(subset2) == nrow(subset2)]
      cli::cli_alert_info("Found {length(consensus2)} genes with method 'intersect' in group2")
    } else if (meth2 == "union") {
      consensus2 <- colnames(subset2)[colSums(subset2) > 0]
      cli::cli_alert_info("Found {length(consensus2)} genes with method 'union' in group2")
    }

    # Combine group1 and group2 using meth_comb
    if (meth.comb == "intersect") {
      final_consensus <- intersect(consensus1, consensus2)
      cli::cli_alert_success("Selecting genes using meth.comb 'intersect' between group1 and group2.")
    } else if (meth.comb == "union") {
      final_consensus <- union(consensus1, consensus2)
      cli::cli_alert_success("Selecting genes using meth.comb 'union' between group1 and group2.")
    }

    consensus.list <- list(consensusGenes = final_consensus)
  } else {
    cli::cli_abort("Function parameters uncorrectly specified. Please refer to the documentation!")
  }

  # Create the result object
  consensus <- list(
    consensusGenes = consensus.list[[1]],
    dataFrame = df,
    inputParams = inputParams
  )

  if (length(consensus$consensusGenes) == 0){
    # If no consensus genes are found print alert
    cli::cli_alert_danger("No features overlapping. Please try with different set up.")
  } else {
    # Final success message
    cli::cli_alert_info("Found {length(consensus$consensusGenes)} genes")
    cli::cli_alert_success("Consensus successfully generated!")
  }

  # Test consensus genes in the adjusted dataset contained in ensBP.obj@data
  testConsensus(df.count = ensBP.obj@data$adjusted.data,
                cons_genes = consensus.list[[1]],
                class = as.factor(ensBP.obj@data$metadata$class))
  cli::cli_alert_success("Consensus genes computed and tested!")

  return(consensus)
}
