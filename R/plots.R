#' @import ggplot2

# ----------------- PLOTS FOR preProcess() ----------------------------
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
        geom_point(
            size = 4,
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
    V1 <- NULL
    Effects <- NULL
    correction <- NULL
    p <- ggplot(pvca_data, aes(x = V1, y = Effects, color = correction)) +
        geom_line(aes(group = Effects), size = 2.5, color = "lightgray", alpha = 0.75) +
        geom_point(size = 12, alpha = 0.95) +
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
            # panel.grid.major.x = element_blank(),
            strip.text = element_text(size = 16, color = "#696969")
        )

    return(p)
}


# ----------------- PLOTS FOR runClassifiers() ----------------------------

# Helper function to plot the intersection of features between models.
upset.plot <- function(data.obj) {
    # Extract features
    model_features_list <- data.obj@model.features
    features <- lapply(model_features_list, function(x) unique(x$Variable))
    names(features) <- names(model_features_list)

    # Create membership matrix
    all_elements <- unique(unlist(features))
    membership_matrix <- sapply(features, function(set) as.integer(all_elements %in% set))
    rownames(membership_matrix) <- all_elements
    membership_matrix <- as.data.frame(membership_matrix)

    # Build Upset plot with ComplexUpset library
    p <- ComplexUpset::upset(
        membership_matrix,
        intersect = colnames(membership_matrix),
        width_ratio = 0.2,
        height_ratio = 1.5,
        name = "",
        set_sizes = ComplexUpset::upset_set_size(
            geom = geom_bar(
                fill = "#EDAE49", width = 0.4
            ),
            position = "right"
        ) + theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.title.x = element_text(size = 18, color = "darkgray"),
            axis.text.x = element_text(size = 13, color = "gray"),
        ),
        base_annotations = list(
            "Intersection size" = ComplexUpset::intersection_size(
                fill = "#00798C",
                width = 0.8,
                text = list(size = 4.5),
            ) + theme(
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = 15, color = "gray"),
                axis.title.y = element_text(size = 20, color = "darkgray"),
            )
        ),
        stripes = "white",
        matrix = (
            ComplexUpset::intersection_matrix(
                segment = geom_segment(linetype = "dotted", color = "black"),
                geom = geom_point(size = 5, shape = "square")
            )
        ),
        themes = ComplexUpset::upset_modify_themes(
            list(
                "intersections_matrix" = theme(
                    text = element_text(size = 20, color = "gray")
                )
            )
        )
    )
    return(p)
}

# Helper function to plot the performances of models during tunning and fitting.
performances.plot <- function(data.obj) {
    performances <- data.obj@performances
    tuning_performances <- performances$tuning_metrics
    fit_performances <- performances$final_metrics

    tuning_performances$type <- "Tuning"
    fit_performances$type <- "Fitting"

    colnames(tuning_performances)[colnames(tuning_performances) == ".metric"] <- "metric"
    colnames(tuning_performances)[colnames(tuning_performances) == "wflow_id"] <- "model"
    colnames(fit_performances)[colnames(fit_performances) == ".metric"] <- "metric"
    colnames(fit_performances)[colnames(fit_performances) == ".estimate"] <- "mean"

    tuning_performances <- tuning_performances[, c("model", "metric", "mean", "type")]
    fit_performances <- fit_performances[, c("model", "metric", "mean", "type")]

    performances <- rbind(tuning_performances, fit_performances)

    # rename accuracy to Accuracy brier_class to Brier score and roc_auc to AUC
    performances$metric <- factor(performances$metric,
        levels = c("accuracy", "f_meas", "precision", "recall"),
        labels = c("Accuracy", "F1-Score", "Precision", "Recall")
    )
    performances$type <- factor(performances$type, levels = c("Fitting", "Tuning"))

    fill_colors <- c("Tuning" = "#8d96a3", "Fitting" = "#2e4057")

    model <- NULL
    mean <- NULL
    type <- NULL
    p <- ggplot(performances, aes(x = model, y = mean, fill = type)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.9, color = "black", width = 0.5, linewidth = 0.3) +
        theme_minimal() +
        scale_fill_manual(values = fill_colors, breaks = c("Tuning", "Fitting")) +
        theme_minimal() +
        labs(y = "Score", x = "", fill = "") +
        coord_flip() +
        theme(
            axis.text.x = element_text(size = 11.5, color = "#7e7e7e"),
            axis.text.y = element_text(size = 13),
            axis.title.x = element_text(size = 14, color = "#7e7e7e"),
            #     legend.position = "bottom",
            #     legend.title = element_text(size = 14),
            #     legend.text = element_text(size = 12),
            strip.text = element_text(size = 13),
            #     panel.grid.minor = element_blank(),
            #     panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            #     plot.title = element_text(size = 18, hjust = 0.5)
        ) +
        facet_wrap(~metric) # , scales = "free_y")

    return(p)
}

wrong.preds.plot <- function(predictions_df) {
    predictions_df$result <- ifelse(predictions_df$.pred_class == predictions_df$class, "Correct", "Wrong")
    wrong_ids <- unique(predictions_df[predictions_df$result == "Wrong", "ID"])
    x_breaks <- ifelse(predictions_df$ID %in% wrong_ids, predictions_df$ID, "")

    p <- ggplot(predictions_df, aes(x = ID, y = model, color = result)) +
        geom_point(size = if_else(predictions_df$result == "Wrong", 4, 1.5), alpha = 0.9) +
        scale_color_manual(values = c("Correct" = "gray", "Wrong" = "indianred")) +
        labs(
            x = "",
            y = "",
            color = "Prediction"
        ) +
        # put breaks on the x-axis only on the Wrong samples
        scale_x_discrete(breaks = x_breaks) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2, margin = margin(t = 50), size = 9, color = "#7c7b7b"),
            axis.text.y = element_text(size = 13),
            panel.grid.y = element_blank(),
            panel.grid.major.x = element_line(color = "lightgray", size = 0.5, linetype = "solid"),
            legend.position = "right"
        )

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
    p <- ggplot(long_results, ggplot2::aes(x = Short_Label, y = Value, fill = Model)) +
        geom_bar(stat = "identity", position = "dodge") + # Bar plot con i modelli affiancati
        facet_wrap(~Metric, scales = "fixed") + # Un riquadro per ogni metrica
        scale_fill_manual(values = model_colors) + # Imposta i colori personalizzati
        theme_minimal() +
        labs(
            title = "Performance Metrics by Model and Gene Combination",
            x = "Gene Combination",
            y = "Metric Value",
            fill = "Model",
            caption = gene_correspondence_text
        ) +
        theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, face = "bold"),
            plot.caption = ggplot2::element_text(
                size = 5, hjust = 0, vjust = 1, margin = ggplot2::margin(t = 10),
                color = "darkgray", face = "italic"
            )
        ) +
        guides(fill = ggplot2::guide_legend(ncol = 1))

    return(p)
}


# ----------------- PLOTS FOR testConsensus() ----------------------------
# Helper function to plot the PCA for consensus genes
pca.consensus <- function(df, cons_genes, class) {
    # Create a matrix from our table of counts
    # rows: samples, cols: genes
    pca_matrix <- as.matrix(df[, colnames(df) %in% cons_genes])

    # Perform the PCA
    cli::cli_alert_info("Generating PCA plots for consensus genes...")
    sample_pca <- stats::prcomp(pca_matrix, scale = TRUE)
    # extract the eigenvalues
    pc_eigenvalues <- sample_pca$sdev^2

    # create a data frame with the eigenvalues
    pc_eigenvalues <- data.frame(
        PC = factor(1:length(pc_eigenvalues)),
        variance = pc_eigenvalues
    )
    # add a new column with the percent variance
    pc_eigenvalues$pct <- pc_eigenvalues$variance / sum(pc_eigenvalues$variance) * 100
    # add another column with the cumulative variance explained
    pc_eigenvalues$pct_cum <- cumsum(pc_eigenvalues$pct)

    # extract the PC scores
    pc_scores <- sample_pca$x
    pc_scores <- data.frame(sample = rownames(pc_scores), pc_scores)

    # extract the PC loadings
    pc_loadings <- sample_pca$rotation
    pc_loadings <- data.frame(gene = rownames(pc_loadings), pc_loadings)

    # reshape the loadings to long format
    pc_loadings_long <- reshape2::melt(pc_loadings, id.vars = "gene", measure.vars = c("PC1", "PC2"), variable.name = "PC", value.name = "loading")
    pc_loadings_long <- pc_loadings_long[order(-abs(pc_loadings_long$loading)), ]
    # get the top 15 genes for each PC
    top_genes <- unique(pc_loadings_long$gene[1:15]) # TODO: make this a parameter
    top_loadings <- pc_loadings[pc_loadings$gene %in% top_genes, ]

    # get the PC scores from prcomp object
    pc_scores_p <- data.frame(sample = rownames(sample_pca$x), sample_pca$x)
    # create the plot
    pca_plot <- ggplot(pc_scores_p, aes(x = PC1, y = PC2, colour = factor(class))) +
        scale_color_manual(values = c("skyblue3", "indianred")) +
        geom_point(size = 3) +
        stat_ellipse(linewidth = 0.7, linetype = 2, type = "norm") +
        stat_ellipse(type = "t") +
        theme_minimal() +
        coord_fixed(ratio = 1) +
        labs(colour = "Class") +
        theme(
            legend.box.just = "right", legend.position = "bottom",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            axis.text = element_text(size = 10, color = "#8f8f8f"),
            axis.title = element_text(size = 13, color = "#8f8f8f"),
        )


    # Create Color Column named 'Col'
    top_loadings$Col <- NA
    top_loadings[top_loadings$PC1 > 0 & top_loadings$PC2 > 0, "Col"] <- "#66C2A5"
    top_loadings[top_loadings$PC1 < 0 & top_loadings$PC2 < 0, "Col"] <- "#E5C494"
    top_loadings[top_loadings$PC1 > 0 & top_loadings$PC2 < 0, "Col"] <- "#8DA0CB"
    top_loadings[top_loadings$PC1 < 0 & top_loadings$PC2 > 0, "Col"] <- "#FC8D62"

    # plot dei top_loadings per ognuna delle PC1 e PC2
    loadings_plot <- ggplot(data = top_loadings) +
        geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
            arrow = arrow(length = unit(0.08, "in")),
            linewidth = 0.8,
            colour = top_loadings$Col
        ) +
        ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.001, size = 2.5
        ) +
        scale_x_continuous(expand = c(0.02, 0.02)) +
        theme_minimal() +
        coord_fixed(ratio = 1) +
        labs(x = "PC1", y = "PC2") +
        theme(
            axis.text = element_text(size = 10, color = "#8f8f8f"),
            axis.title = element_text(size = 13, color = "#8f8f8f"),
            panel.grid.minor = element_blank()
        )
    p <- (pca_plot | loadings_plot)
    cli::cli_alert_success("Plots created!")

    return(p)
}

auroc.FC.consensus <- function(df, cons_genes, class) {
    df <- df[, colnames(df) %in% cons_genes]
    results <- tibble::tibble(
        gene = character(),
        auroc = numeric(),
        auroc_upper = numeric(),
        auroc_lower = numeric(),
        FC = numeric()
    )

    cli::cli_alert_info("Computing AUROC and Fold Change ...")
    for (gene in colnames(df)) {
        predictor <- df[[gene]]
        auc_obj <- pROC::roc(class, predictor, auc = TRUE, ci = TRUE, direction = "<", levels = (c(0, 1)))

        # Calculate Fold Change manually
        pred_with_class <- data.frame(predictor, class)
        FC <- mean(pred_with_class[pred_with_class$class == 1, "predictor"]) -
            mean(pred_with_class[pred_with_class$class == 0, "predictor"])

        results <- results %>%
            dplyr::add_row(
                gene = gene,
                auroc = auc_obj$auc[1],
                auroc_upper = auc_obj$ci[1],
                auroc_lower = auc_obj$ci[3],
                FC = FC
            )
    }
    cli::cli_alert_success("AUROC and Fold Change computed successfully!")

    # Auroc Plot
    cli::cli_alert_info("Creating AUROC plot ...")
    p <- results %>%
        group_by(gene) %>%
        dplyr::mutate(mean_auroc = mean(auroc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gene = forcats::fct_reorder(gene, mean_auroc)) %>% # Order by AUC
        ggplot(aes(x = gene, y = auroc)) +
        geom_hline(yintercept = 0.5, linetype = "solid", color = "#7e7e7e", linewidth = 0.4) +
        geom_errorbar(aes(ymin = auroc_lower, ymax = auroc_upper, color = FC), width = 0, linewidth = 0.4, position = position_dodge(0.5)) +
        geom_point(aes(color = FC), size = 5, position = position_dodge(0.5)) +
        scale_color_gradient2(low = "#187498", mid = "gray90", high = "#C62E2E", midpoint = 0, limits = c(-2.1, 2.1)) +
        # geom_line(aes(y = mean_auc, group = 1), color = "green", linetype = "solid", linewidth = 0.2) +
        scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        scale_shape_manual(values = c(16, 17, 15, 3, 4, 18, 8)) +
        labs(x = "", y = "AUROC") +
        # coord_flip() +
        theme_minimal() +
        theme(
            # legend.position = "right",
            # legend.justification = "center",
            axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1, color = "#7e7e7e"),
            axis.text.y = element_text(hjust = 1, size = 13, color = "#7e7e7e"),
            axis.title.y = element_text(size = 16, color = "#7e7e7e"),
            legend.title = element_text(size = 12, vjust = 0.5, hjust = 0),
            legend.text = element_text(size = 11, hjust = 0),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.ticks = element_line(color = "#7e7e7e", size = 0.3),
            axis.line.x = element_line(color = "#7e7e7e", linewidth = 0.3),
            axis.line.y = element_line(color = "#7e7e7e", linewidth = 0.3)
        )
    cli::cli_alert_success("AUROC plot created successfully!")

    return(p)
}

# This function performs the test of the consensus genes using a MLP model
mlp.model.plot <- function(df, cons_genes, class) {
    # future::plan(future::multisession, workers = parallel::detectCores() - 1)

    withr::with_seed(
        seed = 123,
        code = {
            # Keep from df only genes of consensus
            df <- df[, colnames(df) %in% cons_genes]
            cli::cli_alert_info("Setting up MLP on {length(cons_genes)} genes ...")

            # Put class column at the end
            df$class <- class

            # Create a new train and test dataset
            cli::cli_alert_info("Splitting data into train and test ...")
            train_test_split <- rsample::initial_split(df, prop = 0.7, strata = class)
            # train_test_split <- split.train.test(df, prop = 0.7, seed = 456)
            train <- rsample::training(train_test_split)
            test <- rsample::testing(train_test_split)
            train_resamples <- rsample::vfold_cv(train, v = 3, strata = class)
            cli::cli_alert_success("Data splitted successfully!")


            # define the recipe
            cli::cli_alert_info("Setting up the recipe ...")
            df_recipe <-
                # which consists of the formula (outcome ~ predictors)
                recipes::recipe(class ~ ., data = df) %>%
                # and some pre-processing steps
                recipes::step_zv(recipes::all_predictors()) %>%
                recipes::step_normalize(recipes::all_numeric()) %>%
                recipes::step_corr(recipes::all_predictors())
            cli::cli_alert_success("Recipe set up successfully!")

            # apply the recipe to the training data
            cli::cli_alert_info("Applying the recipe to the training data ...")
            df_train_preprocessed <- df_recipe %>%
                # apply recipe to training data
                recipes::prep(train) %>%
                # extract pre-processed data
                recipes::juice()
            cli::cli_alert_success("Recipe applied successfully!")

            ## Setup the models
            # mlp model
            cli::cli_alert_info("Setting up MLP model ...")
            mlp_model <-
                # specify that the model is a random forest
                parsnip::mlp() %>%
                # specify parameters that need to be tuned
                parsnip::set_args(
                    hidden_units = tune(),
                    penalty = tune(),
                    epochs = tune()
                ) %>%
                # select the engine/package that underlies the model
                parsnip::set_engine("nnet") %>%
                # choose either the continuous regression or binary classification mode
                parsnip::set_mode("classification")

            # set the workflow
            mlp_workflow <- workflows::workflow() %>%
                # add the recipe
                workflows::add_recipe(df_recipe) %>%
                # add the model
                workflows::add_model(mlp_model)
            cli::cli_alert_success("MLP model set up successfully!")

            # extract results
            cli::cli_alert_info("Tuning MLP model ...")
            mlp_tune_results <- mlp_workflow %>%
                tune::tune_grid(
                    resamples = train_resamples, # CV object
                    grid = 10
                )
            cli::cli_alert_success("MLP model tuned successfully!")

            # print results
            mlp_tune_results %>%
                workflowsets::collect_metrics()

            # finalize the workflow
            cli::cli_alert_info("Finalizing MLP model ...")
            param_final <- mlp_tune_results %>%
                tune::select_best(metric = "roc_auc")

            mlp_workflow <- mlp_workflow %>%
                tune::finalize_workflow(param_final)
            cli::cli_alert_success("MLP model finalized successfully!")

            # Evaluate on test set
            cli::cli_alert_info("Evaluating MLP model on test set ...")
            mlp_fit <- mlp_workflow %>%
                # fit on the training set and evaluate on test set
                tune::last_fit(train_test_split)

            # Inspect performances on test set
            test_performance <- mlp_fit %>% workflowsets::collect_metrics()

            # Fitting and using final model
            final_model <- parsnip::fit(mlp_workflow, df)
            cli::cli_alert_success("MLP model evaluated successfully!")

            # Assess variable importance
            cli::cli_alert_info("Assessing variable importance ...")
            mlp_obj <- workflowsets::extract_fit_parsnip(final_model)

            importances <- as.data.frame(vip::vi(mlp_obj))
            cli::cli_alert_success("Variable importance assessed successfully!")
        }
    )

    # future::plan(future::sequential)

    # rescale Importance column in importances df in range -1 and 1
    importances$Importance <- scales::rescale(importances$Importance, to = c(-1, 1))
    importances$type <- ifelse(importances$Importance < 0, 0.026, -0.026)
    importances$type_hjust <- ifelse(importances$Importance < 0, 0, 1)
    # importances$score_hjust <- ifelse(importances$Importance < 0, 2, -1)

    # Create the plot
    cli::cli_alert_info("Creating plot ...")
    p <- ggplot(
        importances,
        aes(x = reorder(Variable, Importance), y = Importance, label = Variable)
    ) +
        geom_segment(aes(y = 0, x = Variable, yend = Importance, xend = Variable), linewidth = 0.8, color = "gray") +
        geom_point(
            stat = "identity",
            size = 4.5,
            color = if_else(importances$Importance < 0, "#995caf", "#e8ab1b")
        ) +
        geom_text(aes(label = Variable, y = type, hjust = type_hjust), color = "black", size = 3.5) +
        # geom_text(aes(label = round(Importance, 2), vjust = -2.1), color = "#7e7e7e", size = 3.5) +
        annotate("text",
            x = length(unique(importances$Variable)) - length(unique(importances$Variable)) * 0.1,
            y = -0.98, label = paste0(
                "MLP metrics: \n",
                "Accuracy: ", round(test_performance$.estimate[[1]], 3), "\n",
                "AUROC: ", round(test_performance$.estimate[[2]], 3), "\n",
                "Brier Score: ", round(test_performance$.estimate[[3]], 3)
            ), hjust = 0, vjust = 0.5, size = 5, color = "#7e7e7e",
            bg = "white",
            box.padding = unit(0.5, "lines")
        ) +
        labs(x = "", y = "Scaled importance") +
        theme_minimal() +
        guides(y = "none") +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.x = element_text(size = 14, color = "#7e7e7e"),
            axis.text.x = element_text(size = 12, color = "#7e7e7e")
        ) +
        coord_flip()
    cli::cli_alert_success("Plot created successfully!")

    return(p)
}

heatmap.plot <- function(
    df, cons_genes,
    class,
    dendrogram = "both",
    show_dendrogram = c(FALSE, TRUE),
    scale = "row",
    custom_colors = c("1" = "#2e1457", "0" = "#66a182"),
    # custom_colors = c("1" = "#9c3f2d", "0" = "#3f7ab9"),
    margins = c(60, 100, 40, 20),
    grid_color = "white",
    grid_width = 0.00001,
    branches_lwd = 0.4,
    fontsize_row = 12,
    fontsize_col = 5,
    colors = NULL,
    hclust.method = "complete",
    distance.method = "euclidean",
    k_row = NULL,
    k_col = NULL,
    Rowv = TRUE,
    Colv = TRUE,
    scale_fill_gradient_fun = NULL,
    limits = NULL,
    na.value = "grey50",
    cellnote = NULL,
    cellnote_size = NULL,
    cellnote_textposition = "middle center",
    key.title = "",
    key.xlab = "Value",
    key.ylab = "Frequency") {
    # Crea i colori laterali per la heatmap
    side_colors_df <- data.frame(class = as.factor(class))
    # Select genes of interest
    df <- df[, colnames(df) %in% cons_genes]
    # Traspose df
    df <- as.data.frame(t(df))


    # Assicurati che siano forniti i colori personalizzati
    if (is.null(colors)) {
        colors <- viridis::viridis(n = 64)
    }

    # Clustering per righe e colonne separatamente
    cli::cli_alert_info("Clustering rows and columns using {hclust.method} ...")
    hc_rows <- stats::hclust(stats::dist(as.matrix(df), method = distance.method), method = hclust.method)
    hc_cols <- stats::hclust(stats::dist(t(df), method = distance.method), method = hclust.method)
    cli::cli_alert_success("Rows and columns clustered successfully!")

    # Crea la heatmap usando heatmaply
    cli::cli_alert_info("Creating heatmap ...")
    p <- heatmaply::heatmaply(
        df,
        Rowv = stats::as.dendrogram(hc_rows),
        Colv = stats::as.dendrogram(hc_cols),
        dendrogram = dendrogram,
        show_dendrogram = show_dendrogram,
        xlab = "",
        ylab = "",
        main = "",
        scale = scale,
        margins = margins,
        grid_color = grid_color,
        grid_width = grid_width,
        titleX = FALSE,
        hide_colorbar = FALSE,
        branches_lwd = branches_lwd,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        showticklabels = c(FALSE, TRUE),
        labRow = rownames(df),
        plot_method = "ggplot",
        colors = colors,
        col_side_colors = side_colors_df,
        col_side_palette = custom_colors,
        scale_fill_gradient_fun = scale_fill_gradient_fun,
        limits = limits,
        na.value = na.value,
        cellnote = cellnote,
        cellnote_size = cellnote_size,
        cellnote_textposition = cellnote_textposition,
        key.title = key.title,
        key.xlab = key.xlab,
        key.ylab = key.ylab
    )
    cli::cli_alert_success("Heatmap created successfully!")

    # Stampa il plot
    return(p)
}
