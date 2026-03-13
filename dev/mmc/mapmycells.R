# Paths
cpm_file  <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/20251215_RSC-204-396_macaque_patchseq_star2.7_cpm.Rdata"
meta_file <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Macaque/patchseq/R_Object/20251215_RSC-204-396_macaque_patchseq_star2.7_samp.dat.Rdata"

# Load CPM
e1 <- new.env()
load(cpm_file, envir = e1)
cpmR <- e1$cpmR

colnames(cpmR)

# Load metadata (object name unknown; we'll detect it)
e2 <- new.env()
load(meta_file, envir = e2)
ls(e2)

# Replace 'METAOBJ' with whatever ls(e2) shows
meta <- e2$samp.dat

identical(meta$exp_component_name, colnames(cpmR))

# cpmR is genes x cells, with colnames = cell IDs
expr_cxg <- t(cpmR)  # cells x genes

# Make sure cell IDs are set as rownames
rownames(expr_cxg) <- meta$exp_component_name
rownames(meta)     <- meta$exp_component_name

# QC
dim(expr_cxg)                  # should be 3839 x 39670
head(rownames(expr_cxg))
head(colnames(expr_cxg))

out_dir <- "~/mmc_work"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(expr_cxg,
            file = file.path(out_dir, "query_cpm_cellsxgenes.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(meta,
            file = file.path(out_dir, "query_meta.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
