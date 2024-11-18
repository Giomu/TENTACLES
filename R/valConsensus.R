
# ------------------- OPZIONE 1 -------------------

# # Carichiamo le librerie necessarie
# library(cluster)
# library(combinat)
# library(MASS)
# library(keras)
# library(caret)
# library(dplyr)
# library(ks)
#
# # Funzione per calcolare il Silhouette Score
# calculate_silhouette <- function(data, labels) {
#   labels <- as.numeric(as.factor(labels))
#   dist_matrix <- dist(data)
#
#   # Controllo per evitare silhouette su matrici di distanza non valide
#   if (any(is.na(dist_matrix)) || any(is.infinite(dist_matrix))) {
#     warning("La matrice delle distanze contiene NA/Inf, impostando il Silhouette Score a NA.")
#     return(NA)
#   }
#
#   sil <- silhouette(as.numeric(labels), dist_matrix)
#   mean(sil[, 3], na.rm = TRUE) # Restituisce il Silhouette Score medio
# }
#
# # Funzione per eseguire PCA e calcolare il Silhouette Score
# pca_silhouette <- function(data, labels, n_components = 2) {
#   pca_result <- prcomp(data, scale. = TRUE)
#
#   # Impostiamo il numero di componenti in base alla dimensione dei dati
#   actual_components <- min(n_components, ncol(pca_result$x))
#
#   # Riduciamo i dati al numero effettivo di componenti disponibili
#   reduced_data <- pca_result$x[, 1:actual_components, drop = FALSE]
#
#   calculate_silhouette(reduced_data, labels)
# }
#
# # Funzione per K-means clustering e valutazione con Adjusted Rand Index (ARI)
# kmeans_ari <- function(data, labels, k = 2) {
#   kmeans_result <- kmeans(data, centers = k, nstart = 25)
#   mclust::adjustedRandIndex(labels, kmeans_result$cluster)
# }
#
# # Funzione per Autoencoder e calcolo Silhouette Score sullo spazio latente
# autoencoder_silhouette <- function(data, labels) {
#   n_features <- ncol(data)
#
#   # Definizione dell'autoencoder
#   input_layer <- layer_input(shape = n_features)
#   encoded <- input_layer %>%
#     layer_dense(units = round(n_features / 2), activation = 'relu') %>%
#     layer_dense(units = round(n_features / 4), activation = 'relu')
#
#   decoded <- encoded %>%
#     layer_dense(units = round(n_features / 2), activation = 'relu') %>%
#     layer_dense(units = n_features, activation = 'linear')
#
#   autoencoder_model <- keras_model(inputs = input_layer, outputs = decoded)
#   encoder_model <- keras_model(inputs = input_layer, outputs = encoded)
#
#   autoencoder_model %>% compile(
#     loss = 'mean_squared_error',
#     optimizer = optimizer_adam()
#   )
#
#   # Addestramento dell'autoencoder
#   autoencoder_model %>% fit(as.matrix(data), as.matrix(data),
#                             epochs = 50, batch_size = 8, verbose = 0)
#
#   # Estrazione delle feature dallo spazio latente
#   latent_space <- encoder_model %>% predict(as.matrix(data))
#
#   # Controllo dei valori NA/NaN/Inf
#   if (any(is.na(latent_space)) || any(is.nan(latent_space)) || any(is.infinite(latent_space))) {
#     warning("Lo spazio latente contiene NA/NaN/Inf, impostando il Silhouette Score a NA.")
#     return(NA)
#   }
#
#   # Calcolo del Silhouette Score sullo spazio latente
#   calculate_silhouette(latent_space, labels)
# }
#
# # Funzione per Kernel Density Estimation (KDE)
# # kde_separation <- function(data, labels) {
# #   # Calcolo della densità usando la funzione kde() del pacchetto 'ks'
# #   density_case <- kde(data[labels == "case", ])
# #   density_control <- kde(data[labels == "control", ])
# #
# #   # Calcolo della distanza tra le due distribuzioni
# #   overlap <- sum(pmin(density_case$estimate, density_control$estimate))
# #   1 - overlap # Più basso è l'overlap, migliore è la separazione
# # }
# kde_separation <- function(subset_data, classes) {
#   # Separazione dei dati in due gruppi basati sulle classi
#   data_control <- subset_data[classes == 0, , drop = FALSE]
#   data_case <- subset_data[classes == 1, , drop = FALSE]
#
#   # Controlli aggiuntivi sui dati
#   cat("\nControlli aggiuntivi:\n")
#   cat("Range control:", range(data_control), "\n")
#   cat("Range case:", range(data_case), "\n")
#
#   # Controlla se i dati contengono solo un valore unico
#   if (length(unique(data_control)) == 1 || length(unique(data_case)) == 1) {
#     cat("Attenzione: dati con varianza nulla\n")
#     return(NA)
#   }
#
#   # Calcolo delle densità kernel
#   tryCatch({
#     kde_control <- kde(data_control)
#     kde_case <- kde(data_case)
#
#     # Verifica se la densità è finita
#     if (any(!is.finite(kde_control$estimate)) || any(!is.finite(kde_case$estimate))) {
#       cat("Densità non finite rilevate\n")
#       return(NA)
#     }
#
#     # Calcolo della separazione (ad esempio, distanza tra le medie)
#     separation_score <- abs(mean(kde_control$estimate) - mean(kde_case$estimate))
#
#     return(separation_score)
#   }, error = function(e) {
#     cat("Errore in kde_separation:", e$message, "\n")
#     return(NA)
#   })
# }
#
# # Funzione per Linear Discriminant Analysis (LDA)
# lda_accuracy <- function(data, labels) {
#   lda_model <- lda(labels ~ ., data = as.data.frame(data))
#   predictions <- predict(lda_model)$class
#   mean(predictions == labels) # Accuratezza della classificazione LDA
# }
#
# # Funzione principale per calcolare le metriche per ogni combinazione di geni
# evaluate_gene_combinations <- function(genes, count_table, classes) {
#   results <- data.frame()
#
#   # Genera tutte le combinazioni di geni (fino a 10)
#   for (i in 1:length(genes)) {
#     gene_combinations <- combn(genes, i, simplify = FALSE)
#
#     for (comb in gene_combinations) {
#       print(paste0("i = ", i, "\ncomb = ", comb))
#       subset_data <- count_table[, comb, drop = FALSE]
#
#       # Calcola tutte le metriche
#       silhouette <- calculate_silhouette(subset_data, classes)
#       pca_sil <- pca_silhouette(subset_data, classes)
#       kmeans_acc <- kmeans_ari(subset_data, classes)
#       autoencoder_sil <- autoencoder_silhouette(subset_data, classes)
#       kde_score <- kde_separation(subset_data, classes)
#       lda_acc <- lda_accuracy(subset_data, classes)
#
#       # Crea una riga con i risultati
#       result_row <- data.frame(
#         Genes = paste(comb, collapse = ", "),
#         Silhouette = silhouette,
#         PCA_Silhouette = pca_sil,
#         Kmeans_ARI = kmeans_acc,
#         Autoencoder_Silhouette = autoencoder_sil,
#         KDE_Score = kde_score,
#         LDA_Accuracy = lda_acc
#       )
#
#       # Aggiungi la riga ai risultati
#       results <- rbind(results, result_row)
#     }
#   }
#
#   return(results)
# }

# Esempio di utilizzo della funzione
# Supponiamo di avere:
# genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
# count_table <- matrix(rnorm(1000), nrow = 100, ncol = 5) # 100 campioni, 10 geni
# colnames(count_table) <- genes
# classes <- sample(c("case", "control"), 100, replace = TRUE)


# genes <- c$consensusGenes
# genes <- genes[1:4]
# count_table <- b@data$adjusted.data
# count_table <- count_table[, genes]
# count_table <- as.matrix(count_table)
# classes <- b@data$metadata$class


# # Calcola i risultati
# result_table <- evaluate_gene_combinations(genes, count_table, classes)
# print(result_table)

# ------------------- OPZIONE 2 -------------------

# Carichiamo i pacchetti necessari
# library(cluster)        # Per silhouette score e clustering
# library(factoextra)     # Per PCA e visualizzazione
# library(Rtsne)          # Per t-SNE
# library(umap)           # Per UMAP
# library(MASS)           # Per LDA
# library(e1071)          # Per SVM
# library(densityClust)   # Per KDE e analisi basate su densità
# library(fpc)            # Per Davies-Bouldin index
# library(DescTools)      # Per Dunn index
# library(EnvStats)       # Per Earth Mover's Distance
# library(ROCR)           # Per ROC AUC
# library(combinat)       # Per tutte le combinazioni possibili
#
# # Funzione aggiornata con controlli
# analyze_gene_combinations <- function(data, labels, genes) {
#
#   # Converti labels in numerico se sono fattori
#   if (is.factor(labels)) {
#     labels <- as.numeric(as.character(labels))
#   }
#
#   # Genera tutte le possibili combinazioni di geni
#   all_combinations <- unlist(lapply(1:length(genes), function(x) combn(genes, x, simplify = FALSE)), recursive = FALSE)
#
#   # Inizializziamo una lista per i risultati
#   results <- list()
#
#   # Funzione per calcolare silhouette score
#   calculate_silhouette <- function(data, labels) {
#     if (length(unique(labels)) > 1) {
#       dist_matrix <- dist(data)
#       sil <- silhouette(labels, dist_matrix)
#       mean(sil[, 3])  # Restituisce il valore medio del silhouette
#     } else {
#       NA
#     }
#   }
#
#   # Loop su tutte le combinazioni di geni
#   for (gene_set in all_combinations) {
#     # Estrai i dati relativi alla combinazione corrente
#     subset_data <- data[, gene_set, drop = FALSE]
#
#     # Verifica che il subset_data non sia vuoto
#     if (ncol(subset_data) == 0) {
#       print(paste("Combinazione vuota per:", paste(gene_set, collapse = ", ")))
#       next
#     }
#
#     # Stampa i dati per debuggare
#     print(paste("Gene set:", paste(gene_set, collapse = ", ")))
#     #print(subset_data)
#
#     # Inizializza una lista per i risultati della combinazione attuale
#     res <- list(
#       Gene_Set = paste(gene_set, collapse = ","),
#       Silhouette = calculate_silhouette(subset_data, labels),
#       Dunn = tryCatch(dunn(as.matrix(dist(subset_data)), labels), error = function(e) NA),
#       Davies_Bouldin = tryCatch(cluster.stats(as.matrix(dist(subset_data)), labels)$davies.bouldin, error = function(e) NA)
#     )
#
#     # PCA + k-means clustering
#     pca_res <- tryCatch(prcomp(subset_data, scale. = TRUE), error = function(e) NULL)
#     if (!is.null(pca_res) && ncol(pca_res$x) >= 2) {
#       kmeans_res <- kmeans(pca_res$x[, 1:2], centers = 2)
#       res$PCA_Separazione <- calculate_silhouette(pca_res$x[, 1:2], kmeans_res$cluster)
#     } else {
#       res$PCA_Separazione <- NA
#     }
#
#     # t-SNE con controllo sulla perplexity
#     num_samples <- nrow(subset_data)
#     tsne_perplexity <- min(30, max(1, floor(num_samples / 3)))
#     tsne_res <- tryCatch({
#       Rtsne(subset_data, perplexity = tsne_perplexity)$Y
#     }, error = function(e) matrix(NA, nrow = num_samples, ncol = 2))
#     if (!any(is.na(tsne_res))) {
#       tsne_kmeans <- kmeans(tsne_res, centers = 2)
#       res$tSNE_Separazione <- calculate_silhouette(tsne_res, tsne_kmeans$cluster)
#     } else {
#       res$tSNE_Separazione <- NA
#     }
#
#     # UMAP
#     umap_res <- tryCatch(umap(subset_data)$layout, error = function(e) matrix(NA, nrow = num_samples, ncol = 2))
#     if (!any(is.na(umap_res))) {
#       umap_kmeans <- kmeans(umap_res, centers = 2)
#       res$UMAP_Separazione <- calculate_silhouette(umap_res, umap_kmeans$cluster)
#     } else {
#       res$UMAP_Separazione <- NA
#     }
#
#     # LDA (Linear Discriminant Analysis)
#     lda_res <- tryCatch(lda(labels ~ ., data = as.data.frame(subset_data)), error = function(e) NULL)
#     if (!is.null(lda_res)) {
#       lda_pred <- predict(lda_res)$class
#       res$LDA_Accuracy <- mean(lda_pred == labels)
#     } else {
#       res$LDA_Accuracy <- NA
#     }
#
#     # ROC AUC (usando SVM con decision_function)
#     svm_model <- tryCatch(svm(as.factor(labels) ~ ., data = as.data.frame(subset_data), probability = TRUE), error = function(e) NULL)
#     if (!is.null(svm_model)) {
#       svm_preds <- predict(svm_model, as.data.frame(subset_data), decision.values = TRUE)
#       pred <- prediction(attributes(svm_preds)$decision.values, labels)
#       auc_perf <- performance(pred, measure = "auc")
#       res$ROC_AUC <- auc_perf@y.values[[1]]
#     } else {
#       res$ROC_AUC <- NA
#     }
#
#     # Wilcoxon Rank Sum Test per ogni gene
#     wilcox_pvalues <- sapply(gene_set, function(gene) {
#       tryCatch(wilcox.test(subset_data[, gene] ~ labels)$p.value, error = function(e) NA)
#     })
#     res$Wilcox_MeanP <- mean(wilcox_pvalues, na.rm = TRUE)
#
#     # Aggiungi i risultati alla lista generale
#     results[[length(results) + 1]] <- res
#   }
#
#   # Filtra i risultati vuoti e converti in una tabella
#   results <- Filter(function(x) !is.null(x), results)
#
#   results_df <- lapply(results, function(x) {
#     data.frame(
#       Gene_Set = x$Gene_Set,
#       Silhouette = x$Silhouette,
#       Dunn = x$Dunn,
#       Davies_Bouldin = ifelse(is.null(x$Davies_Bouldin), NA, x$Davies_Bouldin),  # Replace NULL with NA
#       PCA_Separazione = x$PCA_Separazione,
#       tSNE_Separazione = x$tSNE_Separazione,
#       UMAP_Separazione = x$UMAP_Separazione,
#       LDA_Accuracy = x$LDA_Accuracy,
#       ROC_AUC = x$ROC_AUC,
#       Wilcox_MeanP = x$Wilcox_MeanP
#     )
#   })
#
#   # Combine all rows into a single dataframe
#   results_df <- bind_rows(results_df)
#
#   # Converti i risultati in un data.frame
#   #results_df <- do.call(rbind, lapply(results, as.data.frame))
#   return(results)
# }
#
#
#
#
# # Esempio d'uso
# # dataset <- matrix(rnorm(1000), nrow = 50, ncol = 20) # Dataset di esempio
# # labels <- sample(c(0, 1), 50, replace = TRUE)        # True labels
# # cons_genes <- c("gene1", "gene2", "gene3", "gene4")  # Lista di geni
# genes <- c$consensusGenes
# genes <- genes[1:7]
# cons_genes <- genes
# count_table <- b@data$adjusted.data
# count_table <- count_table[, genes]
# count_table <- as.matrix(count_table)
# dataset <- count_table
# labels <- as.numeric(b@data$metadata$class)
# results_table <- analyze_gene_combinations(dataset, labels, cons_genes)
# print(results_table)
#
#
# library(dplyr)
#
# # Create a function to convert each list element into a dataframe row
# results_df <- lapply(results_table, function(x) {
#   data.frame(
#     Gene_Set = x$Gene_Set,
#     Silhouette = x$Silhouette,
#     Dunn = x$Dunn,
#     Davies_Bouldin = ifelse(is.null(x$Davies_Bouldin), NA, x$Davies_Bouldin),  # Replace NULL with NA
#     PCA_Separazione = x$PCA_Separazione,
#     tSNE_Separazione = x$tSNE_Separazione,
#     UMAP_Separazione = x$UMAP_Separazione,
#     LDA_Accuracy = x$LDA_Accuracy,
#     ROC_AUC = x$ROC_AUC,
#     Wilcox_MeanP = x$Wilcox_MeanP
#   )
# })
#
# # Combine all rows into a single dataframe
# results_df <- bind_rows(results_df)
#
# # View the results
# head(results_df)



# ------------------- OPZIONE 3 -------------------
# Required Libraries
library(cluster)        # For silhouette
library(fpc)            # For Davies-Bouldin
library(mclust)         # For Gaussian Mixture Model
library(factoextra)     # For clustering evaluation
library(umap)           # For UMAP
library(Rtsne)          # For t-SNE
library(ggplot2)        # For plotting

# Function for Clustering and Dimensionality Reduction
evaluate_clustering <- function(count_table, gene_list, labels) {

  # Subset the count table based on the given genes
  count_subset <- count_table[, gene_list, drop = FALSE]
  # Ensure the count data is numeric and remove non-numeric columns
  count_subset <- as.data.frame(lapply(count_subset, as.numeric))
  count_subset <- scale(count_subset) # Standardize the data
  rownames(count_subset) <- rownames(count_table)

  # 1. Clustering Methods
  clustering_results <- list()

  # a) K-Means Clustering
  set.seed(123) # for reproducibility
  kmeans_result <- stats::kmeans(count_subset, centers = 2, iter.max = 100, algorithm = "MacQueen")
  clustering_results$KMeans <- list(
    clusters = kmeans_result$cluster,
    silhouette = cluster::silhouette(kmeans_result$cluster, factoextra::get_dist(count_subset)),
    avg.silwidth = mean(as.data.frame(cluster::silhouette(kmeans_result$cluster, factoextra::get_dist(count_subset)))$sil_width)
  )

  # c) Gaussian Mixture Model (GMM)
  gmm_result <- Mclust(count_subset, G = 2)
  clustering_results$GMM <- list(
    clusters = gmm_result$classification,
    silhouette = cluster::silhouette(gmm_result$classification, factoextra::get_dist(count_subset)),
    avg.silwidth = mean(as.data.frame(cluster::silhouette(gmm_result$classification, factoextra::get_dist(count_subset)))$sil_width)
  )

  # d) Hierarchical Clustering
  hc_result <- hclust(dist(count_subset), method = "complete") # complete linkage
  hc_clusters <- cutree(hc_result, k = 2) # assuming 3 clusters
  clustering_results$Hierarchical <- list(
    clusters = hc_clusters,
    avg.silwidth = fpc::cluster.stats(dist(count_subset), hc_clusters)$avg.silwidth
  )

  # 2. Dimensionality Reduction and Separation Metrics
  reduction_results <- list()

  # a) PCA
  pca_result <- prcomp(count_subset, scale = F)
  pca_data <- data.frame(pca_result$x)
  pca_clusters <- kmeans(pca_data[, 1:2], centers = 2)$cluster
  reduction_results$PCA <- list(
    clusters = pca_clusters,
    avg.silwidth = fpc::cluster.stats(dist(pca_data[, 1:2]), pca_clusters)$avg.silwidth
  )

  # b) t-SNE
  tsne_result <- Rtsne(count_subset, dims = 2, perplexity = 10)
  tsne_data <- tsne_result$Y
  tsne_clusters <- kmeans(tsne_data, centers = 2)$cluster
  reduction_results$tSNE <- list(
    clusters = tsne_clusters,
    avg.silwidth = cluster.stats(dist(tsne_data), tsne_clusters)$avg.silwidth
  )

  # c) UMAP
  umap_result <- umap(count_subset)
  umap_data <- as.data.frame(umap_result$layout)
  umap_clusters <- kmeans(umap_data, centers = 3)$cluster
  reduction_results$UMAP <- list(
    clusters = umap_clusters,
    silhouette = cluster.stats(dist(umap_data), umap_clusters)$avg.silwidth
  )

  # 3. Return Results
  return(list(
    clustering = clustering_results,
    dimensionality_reduction = reduction_results
  ))
}

# Example of calling the function
# Example count_table (gene expression data) and labels (true labels for validation)
# Assuming 'count_table' is a data frame with genes as rows and samples as columns.
# Assuming 'gene_list' is a vector of genes you want to analyze.
# labels is a vector of known labels for validation purposes.



# count_table <- as.data.frame((b@data$adjusted.data))
# gene_list <- c$consensusGenes[1:5]
# labels <- as.numeric(b@data$metadata$class)
# result <- evaluate_clustering(count_table, gene_list, labels)







































