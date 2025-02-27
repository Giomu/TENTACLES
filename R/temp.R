# data(acc.count, acc.clin)

# # Keep 100 genes only
# acc.count <- acc.count[, 1:100]

# # acc.clin.2 <- acc.clin[-c(1:50), ]

# t1 <- preProcess(acc.count, acc.clin,
#     is.normalized = FALSE,
#     batch = "patient.gender", plot = TRUE
# )

# options(future.globals.maxSize = 2 * 1024^3)

# # # Run classifiers
# rc <- runClassifiers(t1,
#     models = c("bag_mlp", "rand_forest"),
#     selector.recipes = c("boruta", "roc"),
#     filter = TRUE, downsample = TRUE, plot = TRUE
# )
