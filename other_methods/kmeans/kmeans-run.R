library(bluster)
library(scuttle)
library(scran)
library(scater)

samples <- c("2720", "6432", "6522", "8667")
clusters = c(7, 7, 7, 6)

# samples <- c("151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676")
# clusters <- c(7, 7, 7, 7, 5, 5, 5, 5, 7, 7, 7, 7)

for (i in 1:length(samples)) {
  sample = samples[i]
  n_cluster = clusters[i]
  spaData <- readRDS(paste0("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-DLPFC/", sample, "/IF", sample, ".rds"))
  # spaData <- readRDS(paste0("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/DLPFC-12/spaData-", sample, ".rds"))

  spaData <- spaData[, colData(spaData)$in_tissue == 1]
  spaData <- logNormCounts(spaData)
  dec <- modelGeneVar(spaData)
  hvgs <- getTopHVGs(dec, n=3000)

  set.seed(1000)
  spaData <- runPCA(spaData, ncomponents=20, subset_row=hvgs)
  spaData <- runUMAP(spaData, dimred="PCA")

  mat <- reducedDim(spaData, "PCA")
  dim(mat)

  set.seed(100)
  kmeans.out <- clusterRows(mat, KmeansParam(n_cluster))
  spaData$kmeans <- kmeans.out

  result = as.data.frame(kmeans.out)
  colnames(result) = c('kmeans')
  write.csv(result, file = paste0('./result/', sample, '_clusters.csv'), row.names = FALSE)
}
