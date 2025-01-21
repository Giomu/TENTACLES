# Helper function to load test datasets
load_test_data <- function(..., package = "yourpkg") {
    # List of datasets to load
    dataset_names <- list(...)
    test_env <- new.env()

    # Load each dataset into the test environment
    for (dataset_name in dataset_names) {
        data(list = dataset_name, package = package, envir = test_env)
    }

    # Return the local environment containing the datasets
    return(test_env)
}
