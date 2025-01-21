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
    colors <- colorRampPalette(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info["Set2", 1], name = "Set2"))(length(unique(pca_data$batch)))
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
