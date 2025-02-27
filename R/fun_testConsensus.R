#' @title Test Consensus
#' @description This function tests a list of consensus genes on a df.count performing PCA,
#' AUROC, Heatmap, and tuning and fitting a MLP Model.
#'
#' @param df.count A data frame containing the counts data. Gene names on columns and samples on rows.
#' @param gene.list A vector list containing consensus genes to test.
#' @param class A vector containing the true class labels.
#' @return The PCA plot, AUROC plot, Heatmap plot, and Validation MLP Model plot
#'
#' @details
#' This function tests the consensus genes by performing PCA, AUROC, Heatmap,
#' and Validation MLP Model on the provided dataset. PCA is performed on the consensus genes,
#' AUROC is computed for each gene, Heatmap is created using the consensus genes,
#' and a MLP model is tuned and fitted on the filtered dataset.
#'
#'
#' @examples
#' \dontrun{
#' testConsensus(
#'     df.count, cons_genes, class
#' )
#' }
#'
#' @export
testConsensus <- function(df.count, gene.list, class) {
    cons_genes <- gene.list
    cli::cli_h2("Consensus Genes Testing")

    # Run the PCA anlysis
    cli::cli_h3("Principal Component Analysis")
    pca_plot <- pca.consensus(df.count, cons_genes, class)
    print(pca_plot)
    cli::cli_alert_success("PCA analysis completed successfully!")

    cli::cli_h3("AUROC Analysis")
    auroc_plot <- auroc.FC.consensus(df.count, cons_genes, class)
    print(auroc_plot)
    cli::cli_alert_success("AUROC analysis completed successfully!")

    cli::cli_h3("Heatmap Analysis")
    heatmap <- heatmap.plot(df.count, cons_genes, class)
    print(heatmap)
    cli::cli_alert_success("Heatmap analysis completed successfully!")

    # Run the validation model
    cli::cli_h3("Validation MLP Model")
    model_plot <- mlp.model.plot(df.count, cons_genes, class)
    print(model_plot)
    cli::cli_alert_success("Validation MLP Model completed successfully!")
}
