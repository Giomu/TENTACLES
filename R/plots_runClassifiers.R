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
