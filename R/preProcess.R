#' @export
# create a class for data import that will be used throughout the package
methods::setClass("preProcess.obj",
  slots = list(
    unsplit = "data.frame",
    splits = "list",
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
  }

  # Filter out samples that are not in the clinical data
  df.clin <- df.clin[m, ]
  return(df.clin)
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

  # Get data information.
  data.info <- data.check(data.type, is.normalized)

  # Use the new function in the data.import function
  df.clin <- match.samples(df.count, df.clin)

  # Transform class labels to binary factors.
  if (!is.null(case.label)) {
    df.clin$class <- as.factor(ifelse(df.clin$class == case.label, 1, 0))
  } else if (all(df.clin$class %in% c(0, 1))) {
    df.clin$class <- as.factor(df.clin$class)
  } else if (is.null(case.label) || !case.label %in% unique(df.clin$class)) {
    # Identify the class levels
    levels <- unique(df.clin$class)
    cli::cli_alert_info("Case label not found or not provided. Using {levels[1]} as the case label.")
    df.clin$class <- as.factor(ifelse(df.clin$class == levels[1], 1, 0))
  }

  # Join the count and clinical data class column
  df.count$class <- df.clin$class

  preProcess.obj <- methods::new("preProcess.obj",
    unsplit = df.count,
    splits = list(),
    metadata = df.clin,
    data.info = data.info
  )

  return(preProcess.obj)
}

# Function to split data into train and test
#' @export
split.train.test <- function(unsplit.df, prop = 0.7, seed = 123) {
  withr::with_seed(
    seed = seed,
    code = {
      # Split the data into training and test sets
      split <- rsample::initial_split(unsplit.df, prop = prop, strata = class)
    }
  )

  return(split)
}

# Helper function to create a split object
create.split.obj <- function(reference.split.obj, training.split, test.split) {
  reference.split.obj$data <- rbind(training.split, test.split)
  reference.split.obj$in_id <- seq_len(nrow(training.split))
  return(reference.split.obj)
}

# Helper function to normalize data
normalization <- function(split, mincpm = 1, minfraction = 0.1) {
  data_train <- rsample::training(split)
  data_test <- rsample::testing(split)

  filter.low.cpm <- function(normalized.counts) {
    keep <- rowSums(edgeR::cpm(normalized.counts) > mincpm) >= ncol(normalized.counts) * minfraction
    return(normalized.counts[keep, ])
  }

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
  data_noclass_train <- data_train[, -ncol(data_train)]
  data_noclass_train <- t(data_noclass_train)
  data_noclass_train <- edgeR::DGEList(counts = data_noclass_train, gene = rownames(data_noclass_train), group = data_train$class)

  # Normalize and filter training data
  norm_train <- edgeR.normalize(data_noclass_train, filter = TRUE)

  # Remove class column, filtered genes from training data and create a DGEList object of test data
  data_noclass_test <- data_test[, -ncol(data_test)]
  data_noclass_test <- t(data_noclass_test)
  data_noclass_test <- data_noclass_test[colnames(norm_train), ]
  data_noclass_test <- edgeR::DGEList(counts = data_noclass_test, gene = rownames(data_noclass_test), group = data_test$class)

  # Normalize test data
  norm_test <- edgeR.normalize(data_noclass_test, filter = FALSE)

  # Add the class column back
  norm_train$class <- data_train$class
  norm_test$class <- data_test$class

  norm_split <- create.split.obj(split, norm_train, norm_test)

  return(norm_split)
}

# Helper function to perform batch correction
correct.batches <- function(split, metadata,
                            batch,
                            covar.mod) {
  data_train <- rsample::training(split)
  data_test <- rsample::testing(split)
  metadata_train <- match.samples(data_train, metadata)
  metadata_test <- match.samples(data_test, metadata)

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
      covar_mod_train <- model.matrix(~ as.factor(metadata_train[[covar.mod]]), data = metadata_train)
      covar_mod_test <- model.matrix(~ as.factor(metadata_test[[covar.mod]]), data = metadata_test)
    } else {
      # When covar.mod is a character vector (multiple columns)
      # Combine the columns to create a single factor variable
      combined_factor_train <- do.call(paste, c(metadata_train[covar.mod], sep = "_")) # Create a combined factor variable
      combined_factor_train <- as.factor(combined_factor_train)
      covar_mod_train <- model.matrix(~combined_factor_train, data = metadata_train)

      combined_factor_test <- do.call(paste, c(metadata_test[covar.mod], sep = "_")) # Create a combined factor variable
      combined_factor_test <- as.factor(combined_factor_test)
      covar_mod_test <- model.matrix(~combined_factor_test, data = metadata_test)
    }
  } else {
    covar_mod_train <- covar_mod_test <- covar.mod
  }

  # Format data for ComBat
  data_noclass_train <- data_train[, -ncol(data_train)]
  data_noclass_test <- data_test[, -ncol(data_test)]
  data_noclass_train <- t(data_noclass_train)
  data_noclass_test <- t(data_noclass_test)
  batch_train <- as.factor(metadata_train[, batch]) # TODO possible bugs after dividing batch and covar.mod between splits. They could disappear from one split.
  batch_test <- as.factor(metadata_test[, batch])

  # Perform batch correction using ComBat
  corrected_train <- suppressMessages(sva::ComBat(dat = data_noclass_train, batch = batch_train, mod = covar_mod_train))
  corrected_test <- suppressMessages(sva::ComBat(dat = data_noclass_test, batch = batch_test, mod = covar_mod_test))

  corrected_train <- as.data.frame(t(corrected_train))
  corrected_test <- as.data.frame(t(corrected_test))
  corrected_train$class <- data_train$class
  corrected_test$class <- data_test$class

  corrected_split <- create.split.obj(split, corrected_train, corrected_test)

  return(corrected_split)
}

# Helper function to adjust test data using training data parameters
train.test.adjustment <- function(split) {
  split_df <- split$data
  split_indices <- split$in_id

  # Remove the class
  df_noclass <- split_df[, -ncol(split_df)]
  df_noclass <- t(df_noclass)

  # Train/Test labels based on the indices
  split_vector <- rep("Test", ncol(df_noclass))
  split_vector[split_indices] <- "Train"
  split_vector <- as.factor(split_vector)

  # Using ComBat to apply parameters from training data to test data.
  adjusted_split_df <- suppressMessages(sva::ComBat(dat = df_noclass, batch = split_vector, mod = NULL, ref.batch = "Train"))
  adjusted_split_df <- as.data.frame(t(adjusted_split_df))

  # Add the class column back
  adjusted_split_df$class <- split_df$class

  # Regenerate the split object with the adjusted data
  adjusted_split <- split
  adjusted_split$data <- adjusted_split_df

  return(adjusted_split)
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
#' @param plot A logical value specifying if PCA and PVCA plots should be generated. Default is TRUE.
#'
#' @return An object of class preProcess.obj containing the pre-processed data.
#'
#' @details
#' The function performs the following steps:
#' 1. Import the data.
#' 2. Split the data into training and test sets.
#' 3. Normalize the data if the data type is RNA-seq and the data is not normalized.
#' 4. Perform batch correction if the batch variable is provided.
#' 5. Generate PCA and PVCA plots before and after batch correction.
#' 6. Apply training data parameters to test data.
#'
#' @import ggplot2
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success cli_alert_danger cli_abort
#' @importFrom withr with_seed
#' @importFrom rsample initial_split training testing
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
#' @importFrom methods new
#'
#' @examples
#' /dontrun{
#' preProcess(df.count, df.clin, class = "class",
#'            case.label = "case", mincpm = 1, minfraction = 0.1,
#'            data.type = "rnaseq", is.normalized = FALSE,
#'            batch = "batch", covar.mod = "covar", plot = TRUE)}
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

  # Split data into training and test sets.
  cli::cli_alert_info("Splitting data into training and test sets...")
  data.obj@splits$raw.split <- split.train.test(data.obj@unsplit)
  cli::cli_alert_success("Data splitted!")

  # Normalize data if data type is RNA-seq and data is not normalized.
  # We assume data array data is already normalized.
  data.info <- data.obj@data.info
  if (data.info[["type"]] == "rnaseq" && data.info[["normalized"]] == FALSE) {
    cli::cli_alert_info("Normalizing and filtering data...")
    norm.split <- normalization(data.obj@splits$raw.split, mincpm = mincpm, minfraction = minfraction)
    data.obj@splits$norm.split <- norm.split
    cli::cli_alert_success("Normalized data!")
  } else {
    data.obj@splits$norm.split <- data.obj@splits$raw.split
  }

  # Perform batch correction if batch variable is provided.
  if (!is.null(batch)) {
    cli::cli_alert_info("Performing batch correction...")
    corrected.split <- correct.batches(data.obj@splits$norm.split, data.obj@metadata, batch = batch, covar.mod = covar.mod)
    data.obj@splits$bc.split <- corrected.split
    cli::cli_alert_success("Batch correction!")

    if (plot) {
      cli::cli_alert_info("Creating PCA plots...")
      p <- batch.pca.plot(norm.split, corrected.split, batch = batch, metadata = data.obj@metadata)
      print(p)
      cli::cli_alert_success("PCA plot generated")

      cli::cli_alert_info("Creating PVCA plots...")
      pv <- batch_pvca_plot(norm.split, corrected.split, class = class, batch = batch, covar = covar.mod, metadata = data.obj@metadata)
      print(pv)
      cli::cli_alert_success("PVCA plots generated")
    }

    cli::cli_alert_info("Applying training data parameters to test data...")
    adjusted.split <- train.test.adjustment(data.obj@splits$bc.split)
    cli::cli_alert_success("Training and Test successfully adjusted!")
  } else {
    cli::cli_alert_info("Applying training data parameters to test data...")
    adjusted.split <- train.test.adjustment(data.obj@splits$norm.split)
    cli::cli_alert_success("Training and Test successfully adjusted!")
  }

  data.obj@splits$adjusted.split <- adjusted.split

  cli::cli_alert_success("Data pre-processing complete!")
  return(data.obj)
}
