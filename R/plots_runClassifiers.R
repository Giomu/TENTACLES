#' Generate an Upset Plot for Feature Overlaps Across Models
#'
#' This function creates an Upset plot using the `ComplexUpset` package to visualize the overlap
#' of selected features across multiple models stored in a `runClassifiers.obj` object.
#'
#' @param data.obj An object of class `runClassifiers.obj`, which contains the selected features for each model
#'        under the slot `@model.features`.
#'
#' @details
#' The function extracts features from `data.obj@model.features` and constructs a membership matrix
#' indicating the presence or absence of each feature across models. The resulting Upset plot
#' displays the intersection sizes and set sizes, highlighting the most frequently selected features.
#'
#' @return A `ggplot2` object representing the Upset plot.
#'
#' @examples
#' \dontrun{
#' # Example usage of upset.plot
#' pp <- preProcess(df.count = acc.count, df.clin = acc.clin,
#'                 batch = "batch", covar.mod = "covar")
#' rc <- runClassifiers(pp, models = c("bag_mlp", "rand_forest", "mlp", "C5_rules"),
#'                      selector.recipes = c("boruta", "roc", "boruta", "boruta"),
#'                      filter = TRUE, downsample = TRUE)
#' upset.plot(rc)
#' }
#'
#' @export
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

#' Generate a Bar Plot of Model Performances
#'
#' This function creates a bar plot using `ggplot2` to visualize the performance metrics
#' of different models during tuning and final fitting stages.
#'
#' @param data.obj An object of class `runClassifiers.obj`, which contains the performance metrics
#'        under the slot `@performances`.
#'
#' @details
#' The function extracts model performance metrics from `data.obj@performances`, separating
#' tuning (`tuning_metrics`) and final fitting (`final_metrics`) performances. Metrics such as
#' accuracy, F1-score, precision, and recall are visualized using grouped bar plots.
#'
#' The bars are colored differently for tuning and fitting phases.
#' The models are displayed on the y-axis, while their performance scores are shown on the x-axis.
#' The plot uses facets to display different metrics.
#'
#' @return A `ggplot2` object representing the bar plot of model performances.
#'
#' @examples
#' \dontrun{
#' # Example usage of performances.plot
#' pp <- preProcess(df.count = acc.count, df.clin = acc.clin,
#'                 batch = "batch", covar.mod = "covar")
#' rc <- runClassifiers(pp, models = c("bag_mlp", "rand_forest", "mlp", "C5_rules"),
#'                      selector.recipes = c("boruta", "roc", "boruta", "boruta"),
#'                      filter = TRUE, downsample = TRUE)
#' performances.plot(rc)
#' }
#'
#' @export
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

#' @title Plot of Model Prediction Errors
#' @description This function creates a scatter plot to visualize correct and incorrect predictions
#' across multiple models. Wrong predictions are highlighted, and only misclassified sample IDs
#' are displayed on the x-axis.
#'
#' @param predictions_df A data frame containing prediction results. It must include the following columns:
#' - `ID`: Unique identifier for each sample.
#' - `.pred_class`: The predicted class label.
#' - `class`: The true class label.
#' - `model`: The name of the model that made the prediction.
#'
#' @return A ggplot2 object displaying a scatter plot where:
#' - Each point represents a prediction made by a model.
#' - Correct predictions are shown in gray, and incorrect predictions in red.
#' - The x-axis only labels misclassified samples.
#'
#' @details
#' The function first compares the predicted class (`.pred_class`) with the true class (`class`)
#' to determine whether each prediction is correct or wrong. It then highlights incorrect predictions
#' with a larger point size and assigns unique IDs on the x-axis only for misclassified samples.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' pp <- preProcess(df.count = acc.count, df.clin = acc.clin,
#'                  batch = "batch", covar.mod = "covar")
#' rc <- runClassifiers(pp, models = c("bag_mlp", "rand_forest", "mlp", "C5_rules"),
#'                      selector.recipes = c("boruta", "roc", "boruta", "boruta"),
#'                      filter = TRUE, downsample = TRUE)
#' wrong.preds.plot(rc@predictions)
#' }
#'
#' @export
wrong.preds.plot <- function(predictions_df) {
    predictions_df$result <- ifelse(predictions_df$.pred_class == predictions_df$class, "Correct", "Wrong")
    wrong_ids <- unique(predictions_df[predictions_df$result == "Wrong", "ID"])
    x_breaks <- ifelse(predictions_df$ID %in% wrong_ids, predictions_df$ID, "")

    p <- ggplot(predictions_df, aes(x = ID, y = model, color = result)) +
      geom_point(size = ifelse(predictions_df$result == "Wrong", 4, 1.5), alpha = 0.9) +
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
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "lightgray", size = 0.5, linetype = "solid"),
          legend.position = "right"
      )

    return(p)
}
