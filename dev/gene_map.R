#### Mapping Data ####
counts = feather::read_feather("den/GreatApes_Macaque_NCBI/counts.feather")
data = feather::read_feather("den/GreatApes_Macaque_NCBI/data.feather")
anno = feather::read_feather("den/GreatApes_Macaque_NCBI/anno.feather")
anno = anno[anno$subclass_label %in% unique(anno$subclass_label)[c(10,20,23)],]
counts = counts[,c(1,names(counts) %in% anno$sample_id)]
counts = cbind.data.frame(gene, counts)


colnames(counts)[1:10]
rownames(counts) = counts$gene
rownames(anno) = anno$sample_id

# feather::write_feather(counts, "data-raw/L5L23_counts.feather")
# feather::write_feather(data.frame(gene), "data-raw/genes_from_counts.feather")
# feather::write_feather(anno, "data-raw/L5L23_anno.feather")

#gene[grep("^MT", gene)]
#gene[grep("^MT", gene)][2]

counts = feather::read_feather("data-raw/L5L23_counts.feather")
gene = counts$gene
counts$gene = NULL

rownames(counts) = gene
colnames(counts)
dimnames(counts)[1]
counts = as.data.frame(counts, row.names = gene)
dimnames(counts)[2]


anno = feather::read_feather("data-raw/L5L23_anno.feather")

#colnames(counts) = data.frame(stringr::str_split(colnames(counts), "-", simplify = TRUE))[,1]
#rownames(counts) = gene

#rownames(anno) = data.frame(stringr::str_split(anno$sample_id, "-", simplify = TRUE))[,1]
anno = as.data.frame(anno, row.names = anno$sample_id)
#rownames(anno) = anno$sample_id

tst = Seurat::CreateSeuratObject(counts = counts, meta.data = anno)

tst[["RNA"]]$data <- NormalizeData(tst[["RNA"]]$counts)


tst[["RNA"]]$data <- NormalizeData(tst[["RNA"]]$counts)
tst[["percent.mt"]] <- PercentageFeatureSet(tst, pattern = "^MT-")

VlnPlot(tst, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(tst, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tst, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


VlnPlot(tst, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), pt.size = 0.2, ncol = 3)

VlnPlot(tst, features = HTRgenesINT, pt.size = 0.2, ncol = 3)

HTRgenes = c(grep("HTR1", counts$gene),
  grep("HTR2", counts$gene),
  grep("HTR3", counts$gene),
  grep("HTR4", counts$gene),
  grep("HTR5", counts$gene),
  grep("HTR6", counts$gene),
  grep("HTR7", counts$gene))


HTRgenesINT = c(
  "HTR1A",
  "HTR1B",
  "HTR1D",
  "HTR1E",
  "HTR1F",
  "HTR2A",
  "HTR2B",
  "HTR2C",
  "HTR3A",
  "HTR3B",
  "HTR3C",
  "HTR3D",
  "HTR3E",
  "HTR4",
  "HTR5A",
  "HTR6",
  "HTR7"
)

