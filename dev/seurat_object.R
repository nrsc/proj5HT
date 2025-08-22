library(Seurat)

data("pbmc_small")
Seurat::DotPlot(object = pbmc_small, features = c("CD247", "CD3E", "CD9"))

pbmc_raw <- read.table(
  file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
  as.is = TRUE
)
pbmc_small <- CreateSeuratObject(counts = pbmc_raw)
