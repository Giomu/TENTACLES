# DO NOT RUN UNLESS COMPLETELY SURE THE RESULTS WILL BE CORRECT

# generate_references <- function() {
#     data(acc.count, acc.clin)

#     # Keep 100 genes only
#     acc.count <- acc.count[, 1:100]

#     pp <- preProcess(acc.count, acc.clin,
#         is.normalized = FALSE,
#         batch = "patient.gender", plot = FALSE
#     )

#     withr::with_options(
#         new = list(future.globals.maxSize = 2 * 1024^3),
#         code = {
#             # Run classifiers
#             rc <- runClassifiers(pp,
#                 models = c("bag_mlp", "rand_forest"),
#                 selector.recipes = c("boruta", "roc"),
#                 filter = TRUE, downsample = TRUE, plot = FALSE
#             )
#         }
#     )

#     saveRDS(pp, "tests/testthat/fixtures/preProcess_rep_test.rds")
#     saveRDS(rc, "tests/testthat/fixtures/runClassifiers_rep_test.rds")
# }
