#' class preProcess.obj
#' @description
#' Creates a class for data import that will be used to stock pre-processed data.
#'
#' @slot raw A data frame containing the raw data.
#' @slot processed A list containing the processed data.
#' @slot metadata A data frame containing the metadata.
#' @slot data.info A list containing the data type and normalization status.
#' @export
methods::setClass("preProcess.obj",
  slots = list(
    raw = "data.frame",
    processed = "list",
    metadata = "data.frame",
    data.info = "list"
  )
)

# Helper function to check data type and normalization status
data.check <- function(data.type = "rnaseq", is.normalized = FALSE) {
  if (!data.type %in% c("rnaseq", "array")) {
    cli::cli_abort("Invalid data type. Please choose from 'rnaseq' or 'array'.")
  }
  if (!is.logical(is.normalized)) {
    cli::cli_abort("is.normalized must be a logical value.")
  }
  output <- list(type = data.type, normalized = is.normalized)
  return(output)
}

# Helper function to match samples in both data tables
match.samples <- function(df.count, df.clin) {
  # Match samples in both data tables.
  m <- match(rownames(df.count), rownames(df.clin))
  # Check if there are common samples between count and clinical data
  if (all(is.na(m))) {
    cli::cli_abort("No common samples found between count and clinical data.
      Count data should have samples as rows and genes as columns.")
  } else if (any(is.na(m))) {
    cli::cli_alert_warning("One or more samples are not present in both tables. These samples will be removed.")
    m <- m[!is.na(m)]
  }

  # Filter out samples that are not in the clinical data
  df.clin <- df.clin[m, ]
  return(rownames(df.clin))
}

# Helper function to import data
# User should provide count table with samples on rows and genes on columns.
# Sample IDs should be the rownames of both count and clinical data.
data.import <- function(
    df.count, df.clin,
    class = "class",
    case.label = NULL,
    data.type = "rnaseq",
    is.normalized = FALSE) {
  # Check if the input data is a data frame
  if (!is.data.frame(df.count) || !is.data.frame(df.clin)) {
    cli::cli_abort("Input data must be in data frame format.")
  }

  # Check if the class column exists in the clinical data
  if (!class %in% colnames(df.clin)) {
    cli::cli_abort("Class column not found in clinical data.")
  }

  # Check if the class column has more or less than 2 unique values
  if (length(unique(df.clin[, class])) != 2) {
    cli::cli_abort("Class column must have exactly 2 unique values.")
  }

  # Get data information.
  data.info <- data.check(data.type, is.normalized)

  # Get matching samples in both data tables
  samples_in_common <- match.samples(df.count, df.clin)
  df.count <- df.count[samples_in_common, ]
  df.clin <- df.clin[samples_in_common, ]

  # Transform class labels to binary factors
  cli::cli_alert_info("Transforming class labels to binary factors...")
  if (all(df.clin[, class] %in% c(0, 1))) {
    cli::cli_alert_info("Binary (0 and 1) class labels found. Assuming '1' as the case label.")
    df.clin[, class] <- as.factor(df.clin[, class])
  } else if (!is.null(case.label) && case.label %in% unique(df.clin[, class])) {
    not.case <- unique(df.clin[, class])[unique(df.clin[, class]) != case.label]
    cli::cli_alert_info(paste("Transforming '{not.case}' to 0 and '{case.label}' to 1."))
    df.clin[, class] <- as.factor(ifelse(df.clin[, class] == case.label, 1, 0))
  } else {
    levels <- unique(df.clin[, class])
    not.case <- levels[levels != levels[1]]
    cli::cli_alert_info("Case label not found or not provided. Using {levels[1]} as the case label.")
    cli::cli_alert_info(paste("Transforming '{not.case}' to 0 and '{levels[1]}' to 1."))
    df.clin[, class] <- as.factor(ifelse(df.clin[, class] == levels[1], 1, 0))
  }

  # Join the count and clinical data class column
  df.count[, class] <- df.clin[, class]

  preProcess.obj <- methods::new("preProcess.obj",
    raw = df.count,
    processed = list(),
    metadata = df.clin,
    data.info = data.info
  )

  return(preProcess.obj)
}

# Helper function to normalize data in log2(CPM + 1) scale
normalization <- function(df.count, class, mincpm = 1, minfraction = 0.1) {
  # Filter Low CPM function
  filter.low.cpm <- function(normalized.counts) {
    keep <- rowSums(edgeR::cpm(normalized.counts) > mincpm) >= ncol(normalized.counts) * minfraction
    return(normalized.counts[keep, ])
  }

  # Normalize data using edgeR function
  edgeR.normalize <- function(data, filter = FALSE) {
    data <- edgeR::calcNormFactors(data)
    norm_data <- edgeR::cpm(data, normalized.lib.sizes = TRUE, log = FALSE)
    if (filter) {
      norm_data <- filter.low.cpm(norm_data)
    }
    norm_data <- log2(norm_data + 1)
    norm_data <- as.data.frame(t(norm_data))
  }

  # Remove the class column and create a DGEList object
  data_noclass <- df.count[, -which(colnames(df.count) == class), drop = FALSE]
  data_noclass <- t(data_noclass)
  data_noclass <- edgeR::DGEList(
    counts = data_noclass, gene = rownames(data_noclass),
    group = df.count[, class]
  )

  # Normalize and filter training data
  norm_data <- edgeR.normalize(data_noclass, filter = TRUE)

  # Add the class column back
  norm_data[, class] <- df.count[, class]

  return(norm_data)
}

# Helper function to perform batch correction using ComBat
correct.batches <- function(data, class, metadata,
                            batch,
                            covar.mod) {
  # Check input parameters
  if (!is.null(batch)) {
    if (!batch %in% names(metadata)) {
      cli::cli_abort("Batch variable {batch} not found in metadata.")
    }
  }

  if (!is.null(covar.mod)) {
    if (!all(covar.mod %in% names(metadata))) {
      cli::cli_abort("covar.mod variable {covar.mod} not found in metadata.")
    }

    # Ensure covar.mod is either a single string or a character vector
    if (length(covar.mod) == 1) {
      # When covar.mod is a single string
      covar_mod_matrix <- stats::model.matrix(~ as.factor(metadata[[covar.mod]]), data = metadata)
    } else {
      # When covar.mod is a character vector (multiple columns)
      # Combine the columns to create a single factor variable
      # Create a combined factor variable
      combined_factor <- do.call(paste, c(metadata[covar.mod], sep = "_"))
      covar_mod_matrix <- stats::model.matrix(~ as.factor(combined_factor), data = metadata)
    }
  } else {
    covar_mod_matrix <- covar.mod
  }

  # Format data for ComBat
  data_noclass <- data[, -which(colnames(data) == class), drop = FALSE]
  data_noclass <- t(data_noclass)
  batch_factor <- as.factor(metadata[, batch])

  # Perform batch correction using ComBat
  corrected_data <- suppressMessages(sva::ComBat(dat = data_noclass, batch = batch_factor, mod = covar_mod_matrix))

  corrected_data <- as.data.frame(t(corrected_data))
  corrected_data[, class] <- data[, class]

  return(corrected_data)
}

#' @title Pre-process data
#' @description This function performs data pre-processing steps including
#' data normalization, batch correction, and PCA and PVCA plots before and after
#' batch correction.
#'
#' @param df.count A data frame containing the count data. Rows are samples and columns are genes.
#' @param df.clin A data frame containing the clinical data. Rows are samples and columns are clinical variables.
#' @param class A character string specifying the column name in df.clin that contains the class labels.
#' @param case.label A character string specifying the case label in the class column. Default is NULL.
#' @param mincpm An integer specifying the minimum count per million (CPM) value for filtering genes. Default is 1.
#' @param minfraction A numeric value specifying the minimum fraction of samples a gene must be present in to be retained. Default is 0.1.
#' @param data.type A character string specifying the data type. Possible values are 'rnaseq' and 'array'. Default is 'rnaseq'.
#' @param is.normalized A logical value specifying if the data is already normalized. Default is FALSE.
#' @param batch A character string specifying the column name in df.clin that contains the batch variable. Default is NULL.
#' @param covar.mod A character string or a character vector specifying the column name(s) in df.clin that contains the covariate(s) for batch correction. Default is NULL.
#' @param plot A logical specifying whether to generate plots. Default is TRUE.
#'
#' @return An object of class preProcess.obj containing the pre-processed data.
#'
#' @details
#' The function performs the following steps:
#' 1. Import the data and match samples in both data tables.
#' 2. Normalize the data if the data type is RNA-seq and the data is not normalized.
#' 3. Perform batch correction if the batch variable is provided.
#' 4. Generate PCA and PVCA plots before and after batch correction.
#'
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom sva ComBat
#'
#' @examples
#' \dontrun{
#' preProcess(df.count, df.clin,
#'   class = "class",
#'   case.label = "case", mincpm = 1, minfraction = 0.1,
#'   data.type = "rnaseq", is.normalized = FALSE,
#'   batch = "batch", covar.mod = "covar", plot = TRUE
#' )
#' }
#'
#' @export
preProcess <- function(
    df.count, df.clin,
    class = "class", case.label = NULL,
    mincpm = 1, minfraction = 0.1,
    data.type = "rnaseq", is.normalized = FALSE,
    batch = NULL, covar.mod = NULL, plot = TRUE) {
  # Print starting Message
  cli::cli_h1("preProcess")

  # Import data
  cli::cli_alert_info("Importing data...")
  data.obj <- data.import(df.count, df.clin, class, case.label, data.type, is.normalized)
  cli::cli_alert_success("Data Imported!")

  # Normalize data if data type is RNA-seq and data is not normalized.
  # We assume data array data is already normalized.
  data.info <- data.obj@data.info
  if (data.info[["type"]] == "rnaseq" && data.info[["normalized"]] == FALSE) {
    cli::cli_alert_info("Normalizing and filtering data...")
    normalized_data <- normalization(data.obj@raw, class = class, mincpm = mincpm, minfraction = minfraction)
    data.obj@processed$normalized <- normalized_data
    cli::cli_alert_success("Data Normalized!")
  } else {
    if (data.info[["type"]] == "array") {
      cli::cli_alert_warning("Array data is assumed to be already normalized.")
    }
    data.obj@processed$normalized <- data.obj@raw
    normalized_data <- data.obj@raw
  }

  # Perform batch correction if batch variable is provided.
  if (!is.null(batch)) {
    cli::cli_alert_info("Performing batch correction...")
    corrected_data <- correct.batches(data.obj@processed$normalized, data.obj@metadata, class = class, batch = batch, covar.mod = covar.mod)
    data.obj@processed$sbatched <- corrected_data
    cli::cli_alert_success("Batch correction complete!")

    if (plot) {
      cli::cli_alert_info("Creating PCA plots...")
      p <- batch.pca.plot(normalized_data, corrected_data, batch = batch, metadata = data.obj@metadata)
      print(p)
      cli::cli_alert_success("PCA plot generated")

      cli::cli_alert_info("Creating PVCA plots...")
      pv <- batch_pvca_plot(normalized_data, corrected_data, class = class, batch = batch, covar = covar.mod, metadata = data.obj@metadata)
      print(pv)
      cli::cli_alert_success("PVCA plots generated")
    }

    adjusted_data <- corrected_data
  } else {
    adjusted_data <- data.obj@processed$normalized
  }

  data.obj@processed$adjusted.data <- adjusted_data

  cli::cli_alert_success("Data pre-processing complete!")
  return(data.obj)
}
