#' @import ggplot2

# ----------------- BATCH PLOTS ----------------------------

# Helper function to plot PCA before and after batch correction
batch.pca.plot <- function(data.before, data.after, batch, metadata) {
    # Helpler function to perform PCA
    perform.pca <- function(data) {
        # Perform PCA
        pca.data <- stats::prcomp(data[, -ncol(data)], scale = TRUE)
        pca.data <- as.data.frame(pca.data$x)
        pca.data <- pca.data[, 1:2] # Keep only first two PCs

        return(pca.data)
    }

    metadata <- match.samples(data.before, metadata)

    pca_before <- perform.pca(data.before)
    pca_before$correction <- "Before"
    pca_before$batch <- metadata[[batch]]

    # Perform PCA after batch correction (train and test)
    pca_after <- perform.pca(data.after)
    pca_after$correction <- "After"
    pca_after$batch <- metadata[[batch]]

    # Combine PCA data
    pca_data <- rbind(pca_before, pca_after)
    pca_data$correction <- factor(pca_data$correction, levels = c("Before", "After"))

    # Plot PCA
    PC1 <- NULL
    PC2 <- NULL
    colors <- colorRampPalette(brewer.pal(brewer.pal.info["Set2", 1], name = "Set2"))(length(unique(pca_data$batch)))
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = batch)) +
        geom_point( # size = 4,
            alpha = 0.75
        ) +
        suppressMessages(ggside::geom_ysidedensity(aes(color = batch), linewidth = 0.65, show.legend = FALSE)) +
        suppressMessages(ggside::geom_xsidedensity(aes(color = batch), linewidth = 0.65, show.legend = FALSE)) +
        labs(color = "Batch") +
        scale_color_manual(values = colors) +
        guides(fill = "none") +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text = element_blank(),
            axis.title = element_text(size = 14, color = "#696969"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            strip.text = element_text(size = 14, color = "#696969"),
            panel.grid.major = element_line(color = "#d3d3d355"),
            panel.border = element_rect(colour = "gray", fill = NA, size = 0.8),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            ggside.panel.grid = element_blank(),
            ggside.panel.border = element_blank(),
            ggside.axis.line = element_blank()
        ) +
        facet_wrap(~correction, scales = "fixed")

    return(p)
}

# Helper function to plot PVCA before and after batch correction
batch_pvca_plot <- function(data.before, data.after, metadata, class, batch, covar) {
    # Match data with metadata
    metadata <- match.samples(data.before, metadata)

    # Run PVCA
    run_PVCA <- function(data, metadata) {
        matrix_data <- t(as.matrix(data))
        metadata <- Biobase::AnnotatedDataFrame(metadata)
        eset <- Biobase::ExpressionSet(assayData = matrix_data, phenoData = metadata)

        if (!is.null(covar)) {
            batch_factors <- c(batch, class, covar)
        } else {
            batch_factors <- c(batch, class)
        }

        pvca_result <- suppressMessages(pvca::pvcaBatchAssess(eset, batch_factors, threshold = 0.6))
        pvca_values <- pvca_result$dat
        colnames(pvca_values) <- pvca_result$label
        pvca_df <- as.data.frame(t(pvca_values))

        return(pvca_df)
    }

    # PVCA before correction
    pvca_before <- run_PVCA(data.before, metadata)
    pvca_before$correction <- "Before"
    pvca_before$Effects <- rownames(pvca_before)

    # PVCA after correction
    pvca_after <- run_PVCA(data.after, metadata)
    pvca_after$correction <- "After"
    pvca_after$Effects <- rownames(pvca_after)

    pvca_data <- rbind(pvca_before, pvca_after)
    pvca_data$correction <- factor(pvca_data$correction, levels = c("Before", "After"))

    # Filter out the resid
    pvca_data <- pvca_data[!grepl("resid", pvca_data$Effects), ]

    # order Effects based on mean V1
    pvca_data$Effects <- factor(pvca_data$Effects, levels = unique(pvca_data[order(pvca_data$V1), "Effects"]))

    # Plot PVCA
    p <- ggplot(pvca_data, aes(x = V1, y = Effects, color = correction)) +
        geom_line(aes(group = Effects), size = 2.5, color = "lightgray", alpha = 0.75) +
        geom_point(size = 10, alpha = 0.95) +
        scale_color_brewer(palette = "Paired") +
        labs(x = "Weighted average proportion variance", color = "Batch correction") +
        scale_x_continuous(
            labels = scales::percent,
            limits = c(0, max(pvca_data$V1)),
            breaks = c(0, max(pvca_data$V1) / 2, max(pvca_data$V1))
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(size = 14, color = "#696969"),
            axis.text.y = element_text(size = 16, color = "#696969"),
            axis.title = element_text(size = 14, color = "#696969"),
            axis.title.y = element_blank(),
            axis.line.x = element_line(color = "gray", size = 0.2),
            axis.ticks.x = element_line(color = "gray", size = 0.8),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.5, "cm"),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(size = 16, color = "#696969"),
        )

    return(p)
}


# ----------------- CLASSIFICATION PLOTS ----------------------------

#' @title upset.plot
#' @description This function generates an UpSet plot to visualize the overlap of
#' features across different models using the ComplexUpSet library.
#'
#' @param data.obj An object of class ensBP.
#'
#' @return An UpSet plot visualizing the overlap of features across different models.
#'
#' @import ggplot2
#' @importFrom ComplexUpset intersection_size upset upset_modify_themes
#' @importFrom UpSetR fromList
#'
#' @examples
#' \dontrun{
#' # Example usage of upset.plot2 function
#' upset.plot(data.obj)}
#'
#' @export
upset.plot <- function(data.obj) {
    # Extract features
    model_features_list <- data.obj@model.features
    features <- lapply(model_features_list, function(x) unique(x$Variable))
    names(features) <- names(model_features_list)

    # Convert Data in right format
    upset_data <- UpSetR::fromList(features)

    # Build Upset plot with ComplexUpset library
    p <- ComplexUpset::upset(
        upset_data,
        intersect = names(features), # Nomi dei set
        width_ratio = 0.2,
        height_ratio = 1.5,
        name = "",
        set_sizes = ComplexUpset::upset_set_size(
            geom = ggplot2::geom_bar(
                fill = "#EDAE49", width = 0.4, alpha = 0.6),
            position = "right"
        ),
        stripes = c(alpha("grey90", 0.45), alpha("white", 0.3)),
        base_annotations = list(
            "Intersection size" = ComplexUpset::intersection_size(
                text = list(size = 2),
                fill = "#00798C",
                width = 0.8
            )
        ),
        themes = ComplexUpset::upset_modify_themes(
            list(
                "intersections_matrix" = ggplot2::theme(
                    axis.text = ggplot2::element_text(size = 10),
                    plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                    panel.grid.major = ggplot2::element_line(linewidth = 0.5),
                    panel.grid.minor = ggplot2::element_blank()
                )
            )
        )
    )

    # Modifica i pallini nella matrice delle intersezioni per renderli quadrati
    p$layers <- lapply(p$layers, function(layer) {
        if ("GeomPoint" %in% class(layer$geom)) {
            layer$geom$default_aes$shape <- 15  # Quadrati
            layer$geom$default_aes$size <- 2.5    # Dimensione
        }
        if ("GeomSegment" %in% class(layer$geom)) {
            layer$geom$default_aes$linetype <- "dotted"  # Linee tratteggiate
            layer$geom$default_aes$size <- 0.2           # Spessore della linea
            layer$geom$default_aes$colour <- "grey30"       # Colore grigio
        }
        layer
    })

    return(p)
}



#' @title performances.plot
#' @description This function generates a bar plot to visualize the performances of different models.
#'
#' @param performances A data frame containing the performances of different models.
#'
#' @return A bar plot visualizing the performances of different models.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # Example usage of performances.plot function
#' performances.plot(performances)
#' }
#'
#' @export
performances.plot <- function(performances) {
    tuning_performances <- performances$tuning_metrics
    test_performances <- performances$final_metrics

    tuning_performances$type <- "Tuning"
    test_performances$type <- "Test"

    colnames(tuning_performances)[colnames(tuning_performances) == ".metric"] <- "metric"
    colnames(tuning_performances)[colnames(tuning_performances) == "wflow_id"] <- "model"
    colnames(test_performances)[colnames(test_performances) == ".metric"] <- "metric"
    colnames(test_performances)[colnames(test_performances) == ".estimate"] <- "mean"

    tuning_performances <- tuning_performances[, c("model", "metric", "mean", "type")]
    test_performances <- test_performances[, c("model", "metric", "mean", "type")]

    performances <- rbind(tuning_performances, test_performances)

    # rename accuracy to Accuracy brier_class to Brier score and roc_auc to AUC
    performances$metric <- factor(performances$metric,
        levels = c("accuracy", "f_meas", "precision", "recall"),
        labels = c("Accuracy", "F-score", "Precision", "Recall")
    )
    performances$type <- factor(performances$type, levels = c("Tuning", "Test"))

    fill_colors <- c("Tuning" = "#8da0cb", "Test" = "#e78ac3")

    # Circular plot using circlize
    # To test, filter by type = Tuining
    # filt_perf <- performances[performances$type == "Tuning", ]
    # circos_colors <- c("firebrick", "dodgerblue", "forestgreen", "darkorange")

    # # Circos parameters
    # circos.par(cell.padding = c(0, 0, 0, 0), track.height = 0.4)

    # # Initialize
    # circos.initialize(factors = filt_perf$model, x = 1, y = filt_perf$mean, xlim = c(0, 1))

    # # Track
    # circos.trackPlotRegion(
    #     factors = filt_perf$model,
    #     y = filt_perf$mean,
    #     panel.fun = function(x, y) {
    #         circos.axis(h = "top", labels.cex = 0.5)
    #         circos.text(
    #             x = 1, y = 0.5, labels = filt_perf$metric, facing = "clockwise",
    #             niceFacing = TRUE, cex = 0.5
    #         )
    #     },
    #     bg.border = NA
    # )

    # circos.trackLines(
    #     factors = filt_perf$model,
    #     y = filt_perf$mean,
    #     col = circos_colors,
    #     lwd = 5
    # )


    model <- NULL
    mean <- NULL
    type <- NULL
    p <- ggplot(performances, aes(x = model, y = mean, fill = type)) +
        geom_bar(stat = "identity", position = "dodge", color = "darkgrey", width = 0.5) +
        labs(title = "Model performances", x = "", y = "", fill = "") +
        theme_minimal() +
        scale_fill_manual(values = fill_colors) +
        theme_minimal() +
        # theme(axis.text.x = element_text(angle = 270, vjust = 0.1, hjust=0.1)) +
        coord_flip() +
        # theme(
        #     axis.text.x = element_text(size = 12),
        #     axis.text.y = element_text(size = 12),
        #     axis.title.x = element_text(size = 14),
        #     axis.title.y = element_text(size = 14),
        #     legend.position = "bottom",
        #     legend.title = element_text(size = 14),
        #     legend.text = element_text(size = 12),
        #     strip.text = element_text(size = 16),
        #     panel.grid.minor = element_blank(),
        #     panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        #     plot.title = element_text(size = 18, hjust = 0.5)
        # ) +
        facet_wrap(~metric) # , scales = "free_y")

    return(p)
}

## TODO: make it more beautiful!!
#' @export
predheat.plot <- function(predictions_df) {
    predictions_df <- predictions_df %>%
        dplyr::mutate(correct = ifelse(class == .pred_class, "Correct", "Incorrect"))

    # Creiamo un plot con ggplot
    p <- ggplot(predictions_df, aes(x = model, y = ID, fill = correct)) +
        geom_tile(
            color = "white", lwd = 1, linetype = 1
        ) +
        scale_fill_manual(
            values = c("Correct" = "lightblue", "Incorrect" = "indianred")
        ) +
        theme_minimal() +
        labs(
            title = "Predictions Heatmap",
            x = "",
            y = "",
            fill = "Prediction Status"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
            axis.text.y = element_text(size = 7)) +
        coord_flip()

    return(p)
}


# ---------------------------------------------- ENRICHMENT PLOTS ----------------------------------------------

# TODO: truncate long y labels. Add a parameter to set the max length of the labels
#' Generate Bar Plot for Enrichment Results
#'
#' This function generates a bar plot visualizing enrichment results based on adjusted p-values.
#' Each bar represents a different enriched module or pathway, ordered by significance.
#'
#' @param enrichment_results A data frame containing enrichment results with columns 'Description' and 'p.adjust'.
#' @param pcutoff Numeric specifying the adjusted p-value cutoff for highlighting significance (default: 0.05).
#' @param low.col Character specifying the color for low p-values (default: "indianred").
#' @param high.col Character specifying the color for high p-values (default: "lightblue").
#'
#' @details
#' The function orders the enrichment results by adjusted p-values and creates a bar plot where each bar represents
#' a module or pathway. Bars are colored based on the adjusted p-value, using a gradient from 'low.col' to 'high.col'.
#' A vertical dashed line is added at -log10(pcutoff) for visualizing the significance cutoff.
#'
#' @import ggplot2
#' @importFrom stats reorder
#'
#' @return A bar plot visualizing enrichment results based on adjusted p-values.
#'
#' @examples
#' \dontrun{
#' # Example usage of bar.plot function
#' # Assuming enrichment_results is a data frame with columns 'Description' and 'p.adjust'
#' bar.plot(enrichment_results = enrichment_results, pcutoff = 0.05)
#' bar.plot(enrichment_results = go_grouped_signif)
#' bar.plot(enrichment_results = go_grouped_signif[go_grouped_signif$p.adjust < 0.08, ])
#' bar.plot(enrichment_results = kegg_res)
#' }
#'
#' @export
bar.plot <- function(enrichment_results, pcutoff = 0.05, low.col = "indianred",
                     high.col = "lightblue") {
    # suppressMessages(library(viridis))

    enrichment_results <- enrichment_results[order(enrichment_results$p.adjust), ]

    ggplot(data = enrichment_results, aes(
        x = -log10(enrichment_results$p.adjust),
        y = reorder(
            enrichment_results$Description,
            -enrichment_results$p.adjust
        ),
        fill = enrichment_results$p.adjust
    )) +
        geom_bar(stat = "identity") +
        # scale_fill_viridis(option = viridis.pal) +
        scale_fill_gradient(low = low.col, high = high.col) +
        geom_vline(
            xintercept = -log10(pcutoff), linetype = "dashed",
            color = "indianred"
        ) +
        theme_minimal() +
        ylab("") +
        xlab("-log10(adjusted p-value)") +
        labs(fill = "p.adj")
}

# TODO: Solve legend problem, it is never in a good position
#' Generate Circos Plot for Enrichment Results
#'
#' This function generates a Circos plot visualizing enrichment results using the chordDiagram function
#' from the circlize package. It plots genes and their enriched modules or pathways based on enrichment
#' results provided.
#'
#' @param enrichment_results A data frame containing enrichment results with columns 'ID' and 'Description'.
#' @param palette_genes Character specifying the color palette for genes (default: "Set2").
#' @param palette_modules Character specifying the color palette for modules or pathways (default: "Paired").
#' @param transparency Numeric specifying the transparency level for the plot elements (default: 0.5).
#' @param facing Character specifying the facing direction of text labels ("clockwise" or "counterclockwise", default: "clockwise").
#' @param cex Numeric specifying the character expansion factor for text labels (default: 0.7).
#' @param legend Logical indicating whether to include a legend in the plot (default: TRUE).
#' @param legend_title Character specifying the title for the legend (default: "circos_plot_legend").
#'
#' @details
#' The function prepares data for the Circos plot by converting enrichment results into a matrix format
#' suitable for chordDiagram. It then creates a Circos plot using chordDiagram from circlize, with genes
#' and modules/pathways represented by different colors. If legend is enabled, it adds a legend to the plot
#' showing the modules or pathways and their descriptions.
#'
#' @importFrom circlize chordDiagram circos.trackPlotRegion circos.axis get.cell.meta.data
#' @importFrom stringr str_trunc str_split
#' @import reshape2
#' @import RColorBrewer
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics par
#'
#' @return  A Circos plot visualizing enrichment results based on genes and modules or pathways.
#'
#' @examples
#' \dontrun{
#' # Example usage of circos.plot function
#' # Assuming enrichment_results is a data frame with columns 'ID' and 'Description'
#' circos.plot(
#'     enrichment_results = enrichment_results,
#'     palette_genes = "Set2", palette_modules = "Paired"
#' )
#' }
#'
#' @export
circos.plot <- function(enrichment_results,
                        palette_genes = "Set2", palette_modules = "Paired",
                        transparency = 0.5, facing = "clockwise",
                        cex = 0.7, legend = TRUE, legend_title = "circos_plot_legend") {
    geni_sel <- unique(unlist(str_split(enrichment_results$geneID, "/")))

    big_dat <- matrix(nrow = length(geni_sel), ncol = nrow(enrichment_results))
    rownames(big_dat) <- geni_sel
    colnames(big_dat) <- c(enrichment_results$ID)

    for (i in 1:ncol(big_dat)) {
        big_dat[, i] <- as.numeric(rownames(big_dat) %in%
            unlist(str_split(enrichment_results$geneID[i], "/")))
    }

    big_dat <- as.data.frame(big_dat)
    df <- big_dat
    df$ID <- rownames(df)
    df <- melt(df, by = df$ID)

    colori_geni <- colorRampPalette(brewer.pal(brewer.pal.info[palette_genes, 1],
        name = palette_genes
    ))(length(unique(df$ID)))
    colori_moduli <- colorRampPalette(brewer.pal(brewer.pal.info[palette_modules, 1],
        name = palette_modules
    ))(length(unique(df$variable)))

    circos_enriched <- function() {
        if (legend == T) {
            par(mar = c(0, 0, 0, 15))
        } else {
            par(mar = c(0, 0, 0, 0))
        }
        circlize::chordDiagram(df,
            transparency = 0.5,
            annotationTrack = c("grid", "axis"),
            preAllocateTracks = 1,
            big.gap = 10,
            grid.col = c(colori_geni, colori_moduli),
            col = c(colori_geni, colori_moduli)
        )

        circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
            xlim <- circlize::get.cell.meta.data("xlim")
            ylim <- circlize::get.cell.meta.data("ylim")
            sector.name <- circlize::get.cell.meta.data("sector.index")
            circlize::circos.text(mean(xlim),
                ylim[1] + .1,
                sector.name,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.7
            )
            circlize::circos.axis(
                h = "top",
                labels.cex = 0.2,
                sector.index = sector.name,
                track.index = 2
            )
        }, bg.border = NA)
        #
        if (legend == T) {
            # Aggiungi la legenda
            legend_df <- enrichment_results[, c("ID", "Description")]
            # Imposta le coordinate per la legenda
            par(xpd = TRUE)
            lgd.x <- par("usr")[2] * 0.9
            lgd.y <- par("usr")[3] * (-0.9)
            legend(
                x = lgd.x, y = lgd.y,
                title = legend_title,
                title.cex = 0.7,
                legend = str_trunc(legend_df$Description, 40),
                fill = colori_moduli,
                col = "black", cex = 0.6, bty = "n", inset = c(0, 0.1)
            )
        }
    }

    circos_enriched()
}



# ----------------- VALIDATION ----------------------------

#' Plot the performance metrics for the top gene combinations.
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
    model_colors <- c(
        "KMeans" = "#83B8C6", "GMM" = "#F3D17C", "HC" = "#C5DBC4",
        "PCA" = "#F49595", "tSNE" = "#C4C7E2", "UMAP" = "#EEC18E"
    )

    # Crea un'etichetta per la casella di testo che mostra le corrispondenze
    gene_correspondence_text <- paste(
        gene_mapping$Short_Label, ":", gene_mapping$Gene_Combination,
        collapse = "\n"
    )


    # Creazione del plot usando ggplot2
    p <- ggplot2::ggplot(long_results, ggplot2::aes(x = Short_Label, y = Value, fill = Model)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") + # Bar plot con i modelli affiancati
        ggplot2::facet_wrap(~Metric, scales = "fixed") + # Un riquadro per ogni metrica
        ggplot2::scale_fill_manual(values = model_colors) + # Imposta i colori personalizzati
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Performance Metrics by Model and Gene Combination",
            x = "Gene Combination",
            y = "Metric Value",
            fill = "Model",
            caption = gene_correspondence_text
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, face = "bold"),
            plot.caption = ggplot2::element_text(
                size = 5, hjust = 0, vjust = 1, margin = ggplot2::margin(t = 10),
                color = "darkgray", face = "italic"
            )
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

    return(p)
}
