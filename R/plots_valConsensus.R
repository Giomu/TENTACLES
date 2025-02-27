# Plot the performance metrics for the top gene combinations.
plotTopMetrics <- function(top_results) {
    # Mappatura delle combinazioni di geni con etichette brevi
    gene_mapping <- data.frame(
        Gene_Combination = unique(top_results$Gene_Combination),
        Short_Label = paste0("comb", seq_along(unique(top_results$Gene_Combination)))
    )

    # Uniamo la mappatura ai risultati
    top_results <- dplyr::left_join(top_results, gene_mapping, by = "Gene_Combination")

    # Preparazione dei dati: da wide a long usando tidyr::pivot_longer
    long_results <- top_results %>%
        tidyr::pivot_longer(
            # cols = starts_with("KMeans"):starts_with("UMAP"), # Seleziona tutte le colonne dei metodi
            tidyselect::matches("^(KMeans|GMM|HC|PCA|tSNE|UMAP)_"),
            names_to = c("Model", "Metric"),
            names_sep = "_", # Divide il nome in base al separatore "_"
            values_to = "Value"
        )

    # Definisci i colori manualmente per ogni modello
    # model_colors <- c(
    #     "KMeans" = "#212129", "GMM" = "#323949", "HC" = "#3d3e51",
    #     "PCA" = "#40445a", "tSNE" = "#4c5265", "UMAP" = "#5C637A"
    # )

    # Crea un'etichetta per la casella di testo che mostra le corrispondenze
    gene_correspondence_text <- paste(
        gene_mapping$Short_Label, ":", gene_mapping$Gene_Combination,
        collapse = "\n"
    )

    # Create a heatmap using ggplot2
    heatmap_plot <- ggplot(long_results, aes(x = Short_Label, y = Model, fill = Value)) +
        geom_tile(color = "white", size = 0.2) +
        scale_fill_gradient2(low = "white", high = "steelblue") +
        facet_wrap(~Metric, scales = "fixed") +
        geom_text(aes(label = round(Value, 2)), color = "black", size = 3) +
        theme_minimal() +
        labs(
            x = "",
            y = "",
            title = "NOT FINISHED YET",
            fill = "Score",
            caption = gene_correspondence_text
        ) +
        theme(
            axis.text.x = element_text(angle = 45, size = 12, hjust = 1, color = "#7c7b7b"),
            axis.text.y = element_text(size = 12, color = "#7c7b7b"),
            panel.grid = element_blank(),
            strip.text = element_text(size = 13, color = "black"),
            plot.caption = element_text(
                size = 5, hjust = 0, vjust = 1, margin = margin(t = 10),
                color = "darkgray", face = "italic"
            )
        ) +
        coord_flip()

    return(heatmap_plot)
}
