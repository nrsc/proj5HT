######################################################################################
#Patch-seq data analysis -  OVERVIEW OF ENTIRE PROJECT

#Code 1. Build the taxonomy using reference SEA-AD data subsampled to 100 cells per cluster

#(Internal Code 1.5) 1.5a. Map ALL data to neuron taxonomy using HANN mapping
#(Internal Code 1.5) 1.5b. Subset to include anything mapping to glutamatergic class (do NOT run QC)
#    ---- Save this starting file as csv with plan to share on GitHub: THIS IS THE STARTING FILE

#Code 2: 2. Map glutamatergic cells to neuron taxonomy using HANN mapping
#---- Be sure to report the top 3-5 mapped cell type for MET typing
#Code 3: 3a. Perform QC, including running Seurat to get Seurat QC pass/fail calls
#---- omit QC fail for downstream analysis
#Code 3: 3b. Create an integrated UMAP with the FACS and patch-seq data using Seurat
#---- share these coordinates for use in Cytosplore

#The issues with the missing samples have been resolved and I have appended them to the most recent R-objects here:
# \\allen\programs\celltypes\workgroups\rnaseqanalysis\SMARTer\STAR\Human\patchseq\R_Object 
# Use the files that begin with 250312, thanks again for the patience while we sorted through this!


######################################################################################
## Overview

# This script describes how to map all of the Allen Institute patch-seq data as of 14 March 2025 to the neuronal subset of the taxonomy for the normal reference data from SEA-AD (https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad), and returns gene expression and metadata tables corresponding to all cells mapping to any glutamatergic types using hierarchical (HANN) mapping regardless of QC parameters. These will be the input files included with the manuscript and the starting point for the next external code block.

# We strongly encourage running this code within the scrattch.taxonomy docker environment.  This code is run using docker://jeremyinseattle/scrattch:0.7.1 and may fail if run in any other environment.

######################################################################################
## Prep the session and open AIT32 for mapping

## Working directory
taxonomyDir = "/allen/programs/celltypes/workgroups/hct/HCT_RNAseq/Jeremy/patchseq_analysis/human_SEAAD_mapping_March2025/FOR_SUBMISSION/"
setwd(taxonomyDir)

## Load the libraries
suppressPackageStartupMessages({
  library(scrattch.taxonomy)
  library(scrattch.mapping)
  library(scrattch.patchseq)
  library(reticulate)
  library(data.table)
})
cell_type_mapper <- import("cell_type_mapper") # For hierarchical mapping
reticulate::use_python("/usr/bin/python3")
set.seed(42)

## Read in the taxonomy
AIT.anndata = loadTaxonomy(taxonomyDir, "AI_taxonomy.h5ad")


######################################################################################
## Read in the Patch-seq data

patchDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Human/patchseq/R_Object/"

load(file.path(patchDir,"20250312_RSC-122-384_human_patchseq_star2.7_intron.Rdata"))
load(file.path(patchDir,"20250312_RSC-122-384_human_patchseq_star2.7_exon.Rdata"))
load(file.path(patchDir,"20250312_RSC-122-384_human_patchseq_star2.7_samp.dat.Rdata"))

patch.counts <- as.matrix(intron)+as.matrix(exon)
patch.logCPM <- logCPM(patch.counts)

patch.anno   <- samp.dat
rownames(patch.anno) <- colnames(patch.counts)


######################################################################################
## Map patch-seq data to neuronal reference

## First set the mapping mode 
AIT.anndata = mappingMode(AIT.anndata, mode="neurons")

## Next run hierarchical mapping
mappingHierarchical <- hierarchicalMapMyCells(AIT.anndata, patch.logCPM)
	
	
######################################################################################
## Save the set of glutamatergic patch-seq cells 
## UPDATE: Also omit a set of cells provided by Rachel that are from non-human or juvenile cases

exclude.cells <- read.csv("human_glut_cells_remove.csv",row.names=1)
exclude <- is.element(rownames(patch.anno),rownames(exclude.cells))
 
is.glutamatergic <- mappingHierarchical$result$class_label_Hierarchical == "Neuronal: Glutamatergic"
keep <- is.glutamatergic&(!exclude)
fwrite(patch.counts[,keep],"patchseq_counts.csv.gz",row.names=TRUE)
fwrite(patch.anno[keep,],"patchseq_annotations.csv")

# End of code for Glutamatergic paper
# ------------------------------------------------------------------------------------------------------------


######################################################################################
## Save the above mapping results for a different set of cells for Femke

femke.cells     <- read.csv("Femke_cellIDs.csv")
mapping.results <- data.frame(cell_id = patch.anno$cell_id, as.data.frame(mappingHierarchical$result),mappingHierarchical$detail)
mapping.femke   <- mapping.results[match(femke.cells$cell_id,mapping.results$cell_id),c(1,11:52)]
mapping.femke$cell_id <- femke.cells$cell_id

fwrite(mapping.femke,"Femke_mapping_results.csv",row.names=FALSE)

# Add second set of Femke's cells from May
femke.cells2     <- read.csv("Femke_cellIDs2.csv")
mapping.femke2   <- mapping.results[match(femke.cells2$cell_id,mapping.results$cell_id),c(1,11:52)]
mapping.femke2$cell_id <- femke.cells2$cell_id

fwrite(mapping.femke2,"Femke_mapping_results2.csv",row.names=FALSE)


