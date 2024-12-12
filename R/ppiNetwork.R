setOldClass("igraph", igraph::make_empty_graph())

#' Class ppiNetwork.obj
#'
#' A class to store the results of the ppiNetwork function, including the network graph, initial interactions, new interactions, final genes, and enrichment results.
#'
#' @slot graph An object of class `igraph` representing the protein-protein interaction network.
#' @slot initialInteractions A `data.frame` containing the initial interactions among the input genes.
#' @slot newInteractions A `data.frame` containing the new interactions among the input genes and their neighbors.
#' @slot finalGenes A `data.frame` containing the final list of genes with their STRING identifiers and descriptions.
#' @slot enrichmentResults A `list` containing the enrichment analysis results for each cluster and the total network.
#' @slot inputParams A list of input parameters used to generate the ppiNetwork result.
#'
#' @export
setClass("ppiNetwork.obj",
  slots = list(
    graph = "igraph",
    initialInteractions = "data.frame",
    newInteractions = "data.frame",
    finalGenes = "data.frame",
    enrichmentResults = "list",
    inputParams = "list"
  )
)

#' @title ppiNetwork
#' @description Protein-Protein Interaction Network Analysis. This function builds a protein-protein
#' interaction (PPI) network from a list of genes, finds communities in the Network, and performs enrichment
#' analysis for each cluster.
#'
#' @param gene_list A character vector of gene identifiers.
#' @param version A character string specifying the STRING database version. Default is "12".
#' @param species An integer specifying the NCBI taxonomy ID of the species. Default is 9606 (Homo sapiens).
#' @param score_threshold A numeric value indicating the minimum interaction score for STRING. Default is 400.
#' @param max_interactors An integer specifying the maximum number of interactors to include. Default is 30.
#' @param frac A numeric value indicating the minimum degree threshold to filter nodes. Default is 8.
#' @param layout A character string specifying the layout for network visualization. Default is "stress".
#' @param circular A logical value indicating whether to plot the network in a circular layout. Default is FALSE.
#' @param palette A character string specifying the color palette for clusters. Default is 'Paired'.
#' @param max.overlaps An integer specifying the maximum number of overlaps for node labels. Default is 20.
#' @param cluster.method A character string specifying the method for community detection.
#'        Options include "leiden", "edge_betweenness", "infomap", "label_prop", "leading_eigen", "spinglass", "walktrap". Default is "louvain".
#' @param size_label Label size. Default is 2.
#' @param bundling Logical. Wether or not bundling edges. Default is TRUE.
#'
#'
#' @return An object of class `ppiNetwork` containing the PPI network, initial and
#' new interactions, final genes, and enrichment results (only significant ones).
#' @examples
#' \dontrun{
#' gene_list <- c("BRCA1", "TP53", "EGFR", "MYC")
#' result <- ppiNetwork(gene_list, cluster.method = "louvain")
#' }
#'
#' @export
ppiNetwork <- function(gene_list, version = "12", species = 9606, score_threshold = 400,
                       max_interactors = 30, frac = 8, layout = "stress", circular = FALSE,
                       palette = "Paired", max.overlaps = 20, cluster.method, bundling = TRUE,
                       size_label = 2) {
  # Save the input parameters in a list
  inputParams <- list(
    gene_list = gene_list, version = version,
    species = species, score_threshold = score_threshold,
    max_interactors = max_interactors, frac = frac,
    cluster.method = cluster.method,
    layout = layout, circular = circular,
    palette = palette, max.overlaps = max.overlaps
  )

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for this function.")
  }

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("ggraph package is required for this function.")
  }

  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("STRINGdb package is required for this function.")
  }

  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("S4Vectors package is required for this function.")
  }


  cli::cli_h1("ppi-Network")
  # Change type of gene_list from list to df
  gene_data <- as.data.frame(gene_list)
  # Set colnames properly
  colnames(gene_data) <- "Gene"

  cli::cli_h2("Genes Mapping")
  # Initialize STRING database
  cli::cli_alert_info("Uploading STRINGdb ...")
  string_db <- STRINGdb::STRINGdb$new(
    version = version, species = species,
    network_type = "full", input_directory = "",
    score_threshold = score_threshold
  )
  cli::cli_alert_success("STRINGdb uploaded!")

  # Map genes to STRING database and remove unmapped ones
  cli::cli_alert_info("Mapping input genes to STRINGdb ...")
  mapped_genes <- string_db$map(gene_data, "Gene", removeUnmappedRows = TRUE, takeFirst = F)
  # Retrieve STRINGid of my genes
  mapped_genes_ids <- mapped_genes$STRING_id
  cli::cli_alert_success("Genes mapped!")
  # Compute interactions of input genes
  cli::cli_alert_info("Computing interactions among input genes ...")
  initial_interactions <- string_db$get_interactions(mapped_genes_ids)
  if (nrow(initial_interactions) == 0) {
    cli::cli_alert_danger("No interactions found!")
  } else {
    cli::cli_alert_success("Interactions computed!")
  }


  # Find neighbors of input genes
  cli::cli_alert_info("Looking for the best {max_interactors} neighbors ...")
  neighbors <- string_db$get_neighbors(mapped_genes_ids)
  # Add input genes to the neighbors list
  all_genes <- unique(c(mapped_genes_ids, neighbors))
  # Compute all interactions among all_genes
  all_interactions <- string_db$get_interactions(all_genes)
  # Remove initial_interactions from all_interactions
  all_interactions <- all_interactions %>%
    dplyr::filter(!(S4Vectors::from %in% mapped_genes_ids & S4Vectors::to %in% mapped_genes_ids))
  # Take the top max_interactors interactions that involve at least one initial gene
  new_interactions <- all_interactions %>%
    filter(S4Vectors::from %in% mapped_genes_ids | S4Vectors::to %in% mapped_genes_ids) %>%
    arrange(desc(combined_score)) %>%
    utils::head(max_interactors)
  # create a vector containing new genes
  new_genes <- unique(c(new_interactions$from, new_interactions$to))
  # Combine new_genes with initial genes to retrieve a final list of genes
  final_genes <- unique(c(mapped_genes_ids, setdiff(new_genes, mapped_genes_ids)))
  names <- string_db$add_proteins_description(data.frame(STRING_id = final_genes))
  # Save string_ids of final_genes for later enrichment
  final_mapped_genes_ids <- final_genes
  cli::cli_alert_success("Neighbors added!")
  cli::cli_alert_success("Gene Mapping done!")

  # Build Network with final_genes
  cli::cli_h2("Network Building")
  g <- string_db$get_subnetwork(final_genes)

  # Change protein names with preferred_name
  igraph::V(g)$name <- names$preferred_name

  # remove unconnected nodes
  cli::cli_alert_info("Removing nodes with degree k < {frac} from Network visualization")
  isolated <- which(igraph::degree(g) < frac)
  g2 <- igraph::delete_vertices(g, isolated)
  isolated <- which(igraph::degree(g2) < frac)
  g2 <- igraph::delete_vertices(g2, isolated)
  cli::cli_alert_success("Nodes removed!")

  cli::cli_alert_info("Starting community detection using {cluster.method} algorithm ...")
  if (cluster.method == "leiden") {
    # Detect communities with leiden algorithm
    cat("\nFinding Network communities with leiden algorithm")
    clst2 <- igraph::cluster_leiden(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else if (cluster.method == "edge_betweenness") {
    # Detect communities with edge_betweenness algorithm
    cat("\nFinding Network communities with edge_betweenness algorithm")
    clst2 <- igraph::cluster_edge_betweenness(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else if (cluster.method == "infomap") {
    # Detect communities with infomap algorithm
    cat("\nFinding Network communities with infomap algorithm")
    clst2 <- igraph::cluster_infomap(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else if (cluster.method == "label_prop") {
    # Detect communities with label_prop algorithm
    cat("\nFinding Network communities with label_prop algorithm")
    clst2 <- igraph::cluster_label_prop(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else if (cluster.method == "leading_eigen") {
    # Detect communities with leading_eigen algorithm
    cat("\nFinding Network communities with leading_eigen algorithm")
    clst2 <- igraph::cluster_leading_eigen(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else if (cluster.method == "spinglass") {
    # Detect communities with spinglass algorithm
    cat("\nFinding Network communities with spinglass algorithm")
    clst2 <- igraph::cluster_spinglass(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else if (cluster.method == "walktrap") {
    # Detect communities with walktrap algorithm
    cat("\nFinding Network communities with walktrap algorithm")
    clst2 <- igraph::cluster_walktrap(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  } else {
    # Detect communities with louvain algorithm
    cat("\nFinding Network communities with louvain algorithm")
    clst2 <- igraph::cluster_louvain(graph = g2)
    # Assign color based on community membership
    igraph::V(g2)$color <- clst2$membership
  }

  igraph::V(g2)$size <- igraph::degree(graph = g2) # , mode = c("all"))#, normalized = TRUE)
  # Set edge weight based on edge betweenness
  # g2 <- set_edge_attr(g2, "weight", value = edge_betweenness(g2))
  cli::cli_alert_success("Communities detected!")

  cli::cli_alert_info("Generating Network visualization ...")
  # Function to create custom labels for communities
  create_cluster_labels <- function(cluster_vector) {
    unique_clusters <- unique(cluster_vector)
    labels <- paste("Cluster", unique_clusters)
    names(labels) <- unique_clusters
    return(labels)
  }

  # Generate Cluster labels
  cluster_labels <- create_cluster_labels(igraph::V(g2)$color)
  ## Generate color palette
  cols_f <- colorRampPalette(RColorBrewer::brewer.pal(8, "Paired"))

  if (bundling == TRUE) {
    # Generate Graph plot
    p <- ggraph::ggraph(g2, layout = layout, circular = circular) +
      # geom_edge_arc(strength = 0.1, width = 0.1, colour = "lightgrey") +
      ggraph::geom_edge_bundle_path(width = 0.1, colour = "lightgrey") +
      ggraph::scale_edge_size_continuous(range = c(0.1, 1)) +
      ggraph::geom_node_point(aes(size = (size / max(size)), colour = factor(color))) +
      scale_color_manual(values = cols_f(length(cluster_labels)), labels = cluster_labels) +
      ggraph::geom_node_label(aes(label = ifelse(size > stats::median(igraph::V(g2)$size), name, NA_character_), size = size),
        repel = TRUE,
        max.overlaps = max.overlaps,
        color = "black",
        fill = "white",
        # alpha = 0.5,
        size = size_label, # log10(V(g2)$size),
        label.r = unit(.1, "pt"),
        label.size = 0,
        label.padding = unit(.1, "lines")
      ) +
      scale_size_continuous(name = "Connectivity", range = c(0.1, 8)) +
      guides(colour = guide_legend(title = "Clusters")) +
      theme_void()
    print(p)
  } else {
    # Generate Graph plot
    p <- ggraph::ggraph(g2, layout = layout, circular = circular) +
      ggraph::geom_edge_arc(strength = 0.1, width = 0.1, colour = "lightgrey") +
      ggraph::geom_edge_link(width = 0.1, colour = "lightgrey") +
      # geom_edge_bundle_path(width = 0.1, colour = "lightgrey") +
      ggraph::scale_edge_size_continuous(range = c(0.1, 1)) +
      ggraph::geom_node_point(aes(size = (size / max(size)), colour = factor(color))) +
      scale_color_manual(values = cols_f(length(cluster_labels)), labels = cluster_labels) +
      ggraph::geom_node_label(aes(label = ifelse(size > stats::median(igraph::V(g2)$size), name, NA_character_), size = size),
        repel = TRUE,
        max.overlaps = max.overlaps,
        color = "black",
        fill = "white",
        # alpha = 0.5,
        size = size_label,
        label.r = unit(.1, "pt"),
        label.size = 0,
        label.padding = unit(.1, "lines")
      ) +
      scale_size_continuous(name = "Connectivity", range = c(0.1, 8)) +
      guides(colour = guide_legend(title = "Clusters")) +
      theme_void()
    print(p)
  }
  cli::cli_alert_success("Plot generated!")

  cli::cli_alert_info("Starting enrichment analysis ...")
  # Perform enrichment analysis for each cluster
  enrichment_results <- list()
  for (cluster in unique(igraph::V(g2)$color)) {
    cli::cli_alert_info("Computing enrichment of Cluster {cluster} ...")
    cluster_genes <- igraph::V(g2)$name[igraph::V(g2)$color == cluster]
    cluster_ids <- names$STRING_id[names$preferred_name %in% cluster_genes]

    # GO enrichment
    cli::cli_alert_info("Computing GO enrichment ...")
    enrichmentGO <- string_db$get_enrichment(cluster_ids, category = "Process")
    # Adjust enrichmentGO colnames to be plotted with enrichPlot.R
    rownames(enrichmentGO) <- enrichmentGO$term
    enrichmentGO$ID <- enrichmentGO$term
    enrichmentGO$term <- NULL
    enrichmentGO$geneID <- enrichmentGO$preferredNames
    enrichmentGO$geneID <- gsub(",", "/", enrichmentGO$geneID, fixed = TRUE)
    enrichmentGO$p.adjust <- stats::p.adjust(enrichmentGO$p_value, method = "BH")
    colnames(enrichmentGO)[which(names(enrichmentGO) == "p_value")] <- "pvalue"
    colnames(enrichmentGO)[which(names(enrichmentGO) == "description")] <- "Description"
    cli::cli_alert_success("GO enrichment computed!")

    # KEGG enrichment
    cli::cli_alert_info("Computing KEGG enrichment ...")
    enrichmentKEGG <- string_db$get_enrichment(cluster_ids, category = "KEGG")
    # Adjust enrichmentKEGG colnames to be plotted with enrichPlot.R
    rownames(enrichmentKEGG) <- enrichmentKEGG$term
    enrichmentKEGG$ID <- enrichmentKEGG$term
    enrichmentKEGG$term <- NULL
    enrichmentKEGG$geneID <- enrichmentKEGG$preferredNames
    enrichmentKEGG$geneID <- gsub(",", "/", enrichmentKEGG$geneID, fixed = TRUE)
    enrichmentKEGG$p.adjust <- stats::p.adjust(enrichmentKEGG$p_value, method = "BH")
    colnames(enrichmentKEGG)[which(names(enrichmentKEGG) == "p_value")] <- "pvalue"
    colnames(enrichmentKEGG)[which(names(enrichmentKEGG) == "description")] <- "Description"
    cli::cli_alert_success("GO enrichment computed!")

    enrichment_results[[paste("Cluster", cluster)]] <- list(
      GO = enrichmentGO,
      KEGG = enrichmentKEGG
    )
    cli::cli_alert_success("Enrichment of Cluster {cluster} done!")
  }

  # Perform enrichment analysis of total cluster
  cli::cli_alert_info("Computing enrichment of total ppi-Network ...")
  # Compute GO enrichment of the network genes in g
  cli::cli_alert_info("Computing GO enrichment ...")
  enrichmentGO <- string_db$get_enrichment(final_mapped_genes_ids, category = "Process")
  # Adjust enrichmentGO colnames to be plotted with enrichPlot.R
  rownames(enrichmentGO) <- enrichmentGO$term
  enrichmentGO$ID <- enrichmentGO$term
  enrichmentGO$term <- NULL
  enrichmentGO$geneID <- enrichmentGO$preferredNames
  enrichmentGO$geneID <- gsub(",", "/", enrichmentGO$geneID, fixed = TRUE)
  enrichmentGO$p.adjust <- stats::p.adjust(enrichmentGO$p_value, method = "BH")
  colnames(enrichmentGO)[which(names(enrichmentGO) == "p_value")] <- "pvalue"
  colnames(enrichmentGO)[which(names(enrichmentGO) == "description")] <- "Description"
  cli::cli_alert_success("GO enrichment computed!")

  # Compute KEGG enrichment of the network genes in g
  cli::cli_alert_info("Computing KEGG enrichment ...")
  enrichmentKEGG <- string_db$get_enrichment(final_mapped_genes_ids, category = "KEGG")
  # Adjust enrichmentKEGG colnames to be plotted with enrichPlot.R
  rownames(enrichmentKEGG) <- enrichmentKEGG$term
  enrichmentKEGG$ID <- enrichmentKEGG$term
  enrichmentKEGG$term <- NULL
  enrichmentKEGG$geneID <- enrichmentKEGG$preferredNames
  enrichmentKEGG$geneID <- gsub(",", "/", enrichmentKEGG$geneID, fixed = TRUE)
  enrichmentKEGG$p.adjust <- stats::p.adjust(enrichmentKEGG$p_value, method = "BH")
  colnames(enrichmentKEGG)[which(names(enrichmentKEGG) == "p_value")] <- "pvalue"
  colnames(enrichmentKEGG)[which(names(enrichmentKEGG) == "description")] <- "Description"
  cli::cli_alert_success("KEGG enrichment computed!")

  enrichment_results[["ppiNetwork"]] <- list(
    GO = enrichmentGO,
    KEGG = enrichmentKEGG
  )
  cli::cli_alert_success("Enrichment of total ppi-Network done!")

  result <- new("ppiNetwork.obj",
    graph = g2,
    initialInteractions = initial_interactions,
    newInteractions = new_interactions,
    finalGenes = names,
    enrichmentResults = enrichment_results,
    inputParams = inputParams
  )
  cli::cli_alert_success("ppiNetwork Analysis successfully accomplished!")
  return(result)
}
