library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)
library(spatialLIBD)

samples <- c(151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676)
clusters <- c(7, 7, 7, 7, 5, 5, 5, 5, 7, 7, 7, 7)
ehub <- ExperimentHub::ExperimentHub()
sce <- fetch_data(type = "sce", eh = ehub)

# samples <- c("2720", "6432", "6522", "8667")
# clusters = c(7, 7, 7, 6)

for (i in 1:length(samples)) {
  sample = samples[i]
  n_cluster = clusters[i]
  spaData <- sce[, colData(sce)$sample_name == sample]
  # spaData <- readVisium(paste0("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-DLPFC/", sample, "/"))
  set.seed(102)
  spaData <- spatialPreprocess(spaData, platform="Visium",
                                n.PCs=7, n.HVGs=2000, log.normalize=TRUE)
  # spaData <- qTune(spaData, qs = seq(2, 10), platform = "Visium", d = 7)
  # qPlot(spaData)
  set.seed(149)
  spaData <- spatialCluster(spaData, q=n_cluster, platform="Visium", d=7,
                             init.method="mclust", model="t", gamma=2,
                             nrep=10000, burn.in=100,
                             save.chain=TRUE)
  # result = as.data.frame(cbind(colData(spaData)$barcode, colData(spaData)$spatial.cluster))
  # colnames(result) = c('spot', 'bayesspace')
  # result <- result[order(result$spot), ]
  result = as.data.frame(colData(spaData)$spatial.cluster)
  colnames(result) = c('bayesspace')
  write.csv(result, file = paste0('./result/', sample, '_clusters.csv'), row.names = FALSE)
}


















