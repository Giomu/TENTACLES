# suppressMessages(library(clusterProfiler))
# suppressMessages(library(org.Hs.eg.db))
# suppressMessages(library(stringr))




#' @title enrich.GO
#' @description Perform Gene Ontology Enrichment Analysis. This function performs Gene Ontology (GO)
#' enrichment analysis using the enrichGO function from the clusterProfiler package.
#' It retrieves enriched GO terms based on a given gene list and specified parameters.
#'
#' @param gene_list Character vector of gene symbols or IDs for which GO enrichment analysis is performed.
#' @param keyType Character specifying the type of gene identifier used in gene_list (default: "SYMBOL").
#' @param pval Numeric specifying the p-value cutoff for enriched GO terms (default: 0.05).
#' @param ont Character specifying the ontology ("BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component) (default: "BP").
#' @param groupGO Logical indicating whether to group similar GO terms (default: TRUE).
#' @param signif Logical indicating whether to filter significantly enriched GO terms based on p-value cutoff (default: TRUE).
#' @param level Integer specifying the depth level to group GO terms if groupGO is TRUE (default: 4).
#'
#' @details
#' The function uses enrichGO to perform GO enrichment analysis based on the specified parameters.
#' It optionally groups similar GO terms if groupGO is TRUE and filters significantly enriched terms
#' based on the specified p-value cutoff if signif is TRUE. The function utilizes org.Hs.eg.db for
#' gene annotation in human species.
#'
#' @importFrom clusterProfiler enrichGO groupGO
#' @import org.Hs.eg.db
#'
#' @return A data frame containing enriched GO terms with associated statistics.
#'
#' @examples
#' \dontrun{
#' # Example usage of enrich.GO function
#' # Assuming gene_list is a vector of gene symbols or IDs
#' enrich.GO(gene_list)
#' go_grouped_signif <- enrich.GO(gene_list = gene_list, pval = 0.05, keyType = "SYMBOL")
#' go_grouped <- enrich.GO(gene_list = gene_list, pval = 0.05, signif = F)
#' go_results <- enrich.GO(gene_list = gene_list, pval = 0.05, signif = F, groupGO = F)
#' go_signif <- enrich.GO(gene_list = gene_list, pval = 0.05, signif = T, groupGO = F)
#' }
#'
#' @export
enrich.GO <- function(gene_list, keyType = "SYMBOL", pval = 0.05,
                      ont = "BP", groupGO = TRUE, signif = TRUE, level = 4) {
  go_results <- clusterProfiler::enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = keyType,
    ont = ont,
    pvalueCutoff = pval,
    pool = 1
  )

  if (groupGO == TRUE & signif == TRUE) {
    padj <- pval
    group_go <- clusterProfiler::groupGO(gene_list,
      OrgDb = "org.Hs.eg.db",
      keyType = keyType, ont = ont,
      level = level
    )
    go_level <- group_go@result$ID
    go_signif <- go_results@result[go_results@result$p.adjust < padj, ]
    go_signif <- go_signif[go_signif$ID %in% go_level, ]
    return(go_signif)
  } else if (groupGO == TRUE & signif == FALSE) {
    padj <- pval
    group_go <- clusterProfiler::groupGO(gene_list,
      OrgDb = "org.Hs.eg.db",
      keyType = keyType, ont = ont,
      level = level
    )
    go_level <- group_go@result$ID
    go_signif <- go_results@result
    go_signif <- go_signif[go_signif$ID %in% go_level, ]
    return(go_signif)
  } else if (groupGO == FALSE & signif == TRUE) {
    padj <- pval
    go_signif <- go_results@result[go_results@result$p.adjust < padj, ]
    return(go_signif)
  } else {
    return(go_results@result)
  }
}


#' @title enrich.kegg
#' @description Perform KEGG Pathway Enrichment Analysis. This function performs KEGG pathway enrichment
#' analysis using the enrichKEGG function from the clusterProfiler package. It retrieves enriched KEGG
#' pathways based on a given gene list and specified parameters.
#'
#' @param gene_list Character vector of gene symbols or IDs for which KEGG pathway enrichment analysis is performed.
#' @param keyType Character specifying the type of gene identifier used in gene_list (default: "SYMBOL").
#' @param pval Numeric specifying the p-value cutoff for enriched KEGG pathways (default: 0.05).
#' @param signif Logical indicating whether to filter significantly enriched KEGG pathways based on p-value cutoff (default: TRUE).
#'
#' @details
#' The function maps gene symbols to entrez IDs using org.Hs.eg.db, performs KEGG pathway enrichment analysis
#' using enrichKEGG, and optionally filters significantly enriched pathways based on the specified p-value cutoff.
#' The organism used is "hsa" (Homo sapiens).
#'
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @return A data frame containing enriched KEGG pathways with associated statistics.
#'
#' @examples
#' \dontrun{
#' # Example usage of enrich.kegg function
#' # Assuming gene_list is a vector of gene symbols or IDs
#' kegg_res_signif <- enrich.kegg(gene_list = gene_list, signif = TRUE)
#' }
#'
#' @export
enrich.kegg <- function(gene_list, keyType = "SYMBOL", pval = 0.05, signif = TRUE) {
  # Map gene SYMBOL into entrezID
  entrez_genes <- as.character(AnnotationDbi::mapIds(org.Hs.eg.db, gene_list, "ENTREZID", keyType))
  kegg_results <- clusterProfiler::enrichKEGG(
    gene = entrez_genes,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = pval
  )

  # Condition to return only significant results or not
  if (signif == TRUE) {
    padj <- pval
    kegg_results_filtered <- kegg_results@result[kegg_results@result$p.adjust < padj, ]
  } else {
    kegg_results_filtered <- kegg_results@result
  }

  # Convert KEGG IDs back to SYMBOL
  kegg_results_filtered$geneID <- sapply(strsplit(kegg_results_filtered$geneID, "/"), function(x) {
    symbols <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, x, "SYMBOL", "ENTREZID"))
    paste(symbols, collapse = "/")
  })

  return(kegg_results_filtered)
}
