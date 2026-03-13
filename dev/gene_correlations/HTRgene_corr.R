brain = readRDS("~/proj5HT/data-raw/Seurat/offPipeline-Macaque_PatchSeq_excitatory.rds")

Seurat::DefaultAssay(brain) <- "RNA"

# create the "data" layer
brain <- Seurat::NormalizeData(brain, assay = "RNA")

# (optional but often helpful) ensure variable features exist
brain <- Seurat::FindVariableFeatures(brain, assay = "RNA")

ieg_genes <- c(
  "FOS","FOSB","JUN","JUNB","JUND",
  "EGR1","EGR2","EGR3",
  "ARC","NPAS4",
  "ATF3",
  "IER2","IER3",
  "DUSP1","DUSP5",
  "NR4A1","NR4A2",
  "BTG2","GADD45B"
)
htr_genes <- c(
  "HTR1A","HTR1B","HTR1D","HTR1E","HTR1F",
  "HTR2A","HTR2B","HTR2C",
  "HTR4","HTR5",
  "HTR7"
)


ieg_genes <- intersect(ieg_genes, rownames(brain[["RNA"]]))

brain <- Seurat::AddModuleScore(
  object   = brain,
  features = list(ieg_genes),
  assay    = "RNA",
  name     = "IEG",
  search   = FALSE
)

# score will be brain$IEG1

expr <- GetAssayData(brain, assay = "RNA", layer = "data")  # or logcounts
ieg_use <- intersect(ieg_genes, rownames(expr))

brain$IEG_score <- Matrix::colMeans(expr[ieg_use, , drop = FALSE])
summary(brain$IEG_score)

htr_use <- intersect(htr_genes, rownames(expr))
htr_mat <- t(expr[htr_use, , drop = FALSE])

colMeans(htr_mat > 0)

cors <- sapply(htr_use, function(g) {
  stats::cor(
    htr_mat[, g],
    brain$IEG_score,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

cors <- sort(cors, decreasing = TRUE)
cors

pvals <- sapply(htr_use, function(g) {
  stats::cor.test(
    htr_mat[, g],
    brain$IEG_score,
    method = "spearman"
  )$p.value
})

data.frame(
  gene = names(cors),
  rho  = cors,
  pval = pvals[names(cors)]
)



library(Seurat)
library(Matrix)

ieg_genes <- c("FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","ARC","NPAS4",
               "ATF3","IER2","IER3","DUSP1","DUSP5","NR4A1","NR4A2","BTG2","GADD45B")

htr_genes <- c("HTR1A","HTR1B","HTR1D","HTR1E","HTR1F","HTR2A","HTR2C","HTR4","HTR7")

DefaultAssay(brain) <- "RNA"

# pull expression matrix: prefer normalized "data", else do quick log-normalized counts
expr <- tryCatch(
  GetAssayData(brain, assay = "RNA", layer = "data"),
  error = function(e) NULL
)

if (is.null(expr)) {
  counts <- GetAssayData(brain, assay = "RNA", layer = "counts")
  lib <- Matrix::colSums(counts)
  expr <- log1p(t(t(counts) / pmax(lib, 1)) * 1e4)  # log1p(CP10K)
}

ieg_use <- intersect(ieg_genes, rownames(expr))
stopifnot(length(ieg_use) >= 5)

brain$IEG_score <- Matrix::colMeans(expr[ieg_use, , drop = FALSE])

DefaultAssay(brain) <- "RNA"



# If you only had counts earlier, ensure normalized data exists
if (!"data" %in% SeuratObject::Layers(brain[["RNA"]])) {
  brain <- NormalizeData(brain)
}

# Standard quick DR
brain <- FindVariableFeatures(brain, nfeatures = 2000)
brain <- ScaleData(brain, features = VariableFeatures(brain))
brain <- RunPCA(brain, features = VariableFeatures(brain))
brain <- RunUMAP(brain, dims = 1:30)

# Now plots work
FeaturePlot(brain, features = "IEG_score", reduction = "umap", order = TRUE)
FeaturePlot(brain, features = intersect(htr_genes, rownames(brain)), reduction = "umap",
            ncol = 3, order = TRUE)

FeaturePlot(brain, features = "IEG_score", order = TRUE)

htr_use <- intersect(htr_genes, rownames(expr))
FeaturePlot(brain, features = htr_use, ncol = 3, order = TRUE)



FeaturePlot(brain, features = c("IEG_score", "HTR1E", "HTR1F", "HTR2A", "HTR4", "HTR7"),
            reduction = "umap", ncol = 3, order = TRUE)

library(Matrix)

# pull expression matrix used earlier
expr <- tryCatch(
  GetAssayData(brain, assay = "RNA", layer = "data"),
  error = function(e) {
    counts <- GetAssayData(brain, assay = "RNA", layer = "counts")
    lib <- Matrix::colSums(counts)
    log1p(t(t(counts) / pmax(lib, 1)) * 1e4)
  }
)

genes_plot <- c("HTR1E","HTR1F","HTR2A","HTR4")
genes_plot <- intersect(genes_plot, rownames(expr))

df <- as.data.frame(t(expr[genes_plot, , drop = FALSE]))
df$IEG_score <- brain$IEG_score

par(mfrow = c(2,2))
for (g in genes_plot) {
  plot(df[[g]], df$IEG_score,
       pch = 16, cex = 0.5,
       xlab = g, ylab = "IEG_score",
       main = paste(g, "vs activity"))
  rho <- cor(df[[g]], df$IEG_score,
             method = "spearman", use = "pairwise.complete.obs")
  mtext(sprintf("Spearman rho = %.2f", rho),
        side = 3, line = -1.2, cex = 0.9)
}
par(mfrow = c(1,1))

library(dplyr)

subclass_var <- "assigned_subclass"  # change if needed

df_sub <- data.frame(
  cell = colnames(expr),
  subclass = brain@meta.data[["subclass_Corr"]],
  IEG_score = brain$IEG_score,
  t(expr[genes_plot, , drop = FALSE])
)

df_mean <- df_sub %>%
  group_by(subclass) %>%
  summarize(
    n_cells = n(),
    IEG = mean(IEG_score, na.rm = TRUE),
    across(all_of(genes_plot), mean, na.rm = TRUE),
    .groups = "drop"
  )

print(df_mean)

mat <- as.matrix(df_mean[, genes_plot])
rownames(mat) <- df_mean$subclass

heatmap(mat, scale = "row", Rowv = NA, Colv = NA,
        col = colorRampPalette(c("grey90","royalblue"))(50),
        margins = c(6,8))

g <- "HTR1F"
keep <- htr_mat[, g] > 0
sum(keep)

cor(htr_mat[keep, g], brain$IEG_score[keep],
    method = "spearman", use = "pairwise.complete.obs")

g <- "HTR4"
keep <- htr_mat[, g] > 0
sum(keep)

cor(htr_mat[keep, g], brain$IEG_score[keep],
    method = "spearman", use = "pairwise.complete.obs")


par(mfrow=c(1,2))
for (g in c("HTR4","HTR1F")) {
  keep <- htr_mat[, g] > 0
  boxplot(brain$IEG_score ~ keep,
          main = g,
          names = c("0", ">0"),
          xlab = paste0(g, " detected?"),
          ylab = "IEG_score")
}
par(mfrow=c(1,1))



