######################################################################################
## Overview

# This script describes how to create an AIT version of the taxonomy for the normal reference data from SEA-AD (https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad), and then create a child taxonomy only include neuronal cells for patch-seq mapping.

# We strongly encourage running this code within the scrattch.taxonomy docker environment.  This code is run using docker://jeremyinseattle/scrattch:0.7.1 and may fail if run in any other environment.

######################################################################################
## Download and prepare reference data

# First let's download the existing QCed data set and associated dendrogram to our working directory and output some information about it.

# Working directory
taxonomyDir = "/allen/programs/celltypes/workgroups/hct/HCT_RNAseq/Jeremy/patchseq_analysis/human_SEAAD_mapping_March2025/FOR_SUBMISSION/"
if(!file.exists(taxonomyDir)) dir.create(taxonomyDir)
setwd(taxonomyDir)


######################################################################################
## Now let's prepare our environment. 

## Load the libraries
suppressPackageStartupMessages({
  library(scrattch.taxonomy)
  library(scrattch.mapping)
  library(scrattch.patchseq)
  library(reticulate)
})
cell_type_mapper <- import("cell_type_mapper") # For hierarchical mapping
set.seed(42)

######################################################################################
## Read in the reference taxonomy
  
#Read in and subset the data set to 100 cells per cluster (to reduce computational burden and more evenly sample cell types).  We note that "cluster_label", "subclass_label", and "class_label" correspond to SEA-AD supertype, subclass, and class, respectively, and are used for defining the hierarchy.  

## Download the reference data to the working directory and read it in
seaad_url  <- "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
dend_url   <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0f/37/0f3755cb-3acb-4b93-8a62-5d6adc74c673/dend.rds"
#download.file(seaad_url,"Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")  # NOTE: we recommend downloading via the web browser, as this command may fail
#download.file(dend_url,"Reference_MTG_dend.rds")
seaad_data <- read_h5ad("Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
seaad_dend <- readRDS("Reference_MTG_dend.rds")

## Subsample data (this can be done either here, or within buildTaxonomy)
keepCells <- subsampleCells(seaad_data$obs$cluster_label,100,seed=42)

## Get (subsampled) subset data and annotations
taxonomy.counts = (seaad_data$X)[keepCells,]
cn <- c("sample_name","cluster_label","cluster_confidence","subclass_label","class_label",
        "external_donor_name_label","age_label","donor_sex_label")
taxonomy.metadata = seaad_data$obs[keepCells,cn]

## Ensure count matrix and annotations are in the same order (this shouldn't be needed)
taxonomy.metadata = taxonomy.metadata[match(rownames(taxonomy.counts), taxonomy.metadata$sample_name),]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))

## Transpose the counts matrix (... for now; code in process to avoid transposing large matrices)
taxonomy.counts <- t(taxonomy.counts)
taxonomy.counts <- as(taxonomy.counts, "dgCMatrix")


######################################################################################
## Create the (parent) AIT Taxonomy

## First we need to compute some variable genes and a UMAP.  This are NOT used in this manuscript but are used for alternative mapping algorithms.

# Compute top 1000 binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.metadata$cluster, 1000)
# Compute UMAP coordinates
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout
# Set rownames to your annotation and UMAP data.frames (Required!)
rownames(umap.coords) = colnames(taxonomy.counts)


## Set up the levels of hierarchy for all mapping functions later
## -- This MUST be from broadest to most specific types, and NOT vice versa
## -- Note that "cluster" is the SEAAD supertypes
hierarchy = list("class_label", "subclass_label", "cluster_label")

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(
      counts        = taxonomy.counts,
      meta.data     = taxonomy.metadata,
      dend          = seaad_dend,  # If this is omitted buildTaxonomy will generate a dendrogram
      feature.set   = binary.genes,
      umap.coords   = umap.coords,
      taxonomyTitle = "AI_taxonomy",  # Determines the file name
      taxonomyDir   = taxonomyDir,
      subsample     = 100, # Not strictly necessary here, since we've already subsampled
      hierarchy     = hierarchy
)



######################################################################################
### Create child taxonomies for patch-seq mapping of all "neurons" and of "glutamatergic" neurons

# First let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. In this case we are going to exclude all the non-neuronal cells since we know that (1) patch-seq cells are all neurons and (2) the transcriptomics of patch-seq cells are often contaminated by non-neuronal transcripts that can "trick" algorithms into thinking the cells are non-neuronal. 

## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "neurons",
                                    subsample = 100, ## Subsampling not needed in this case
                                    subclass.column = "subclass_label", 
                                    class.column = "class_label", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal and Non-neural"), ## The off-target class.column labels for patchseqQC.
                                    subclass.subsample = 100, ## Subsampling is for PatchseqQC contamination calculation.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = taxonomyDir) # This will create a subfolder in the reference taxonomy directory

#The buildPatchseqTaxonomy function does the following, updating the anndata variable and file accordingly:
# - Creates a "child" taxonomy that only includes the cells and clusters requested (e.g., neurons for patch-seq mapping)
# - Defines marker genes for each node of the dendrogram for use with tree mapping
# - Creates a subsetted version of the parent dendrogram that also includes node marker genes
# - Creates the required marker and expression variables for 'QC_markers' for use with patchseqQC
# - Creates a table of cell to cluster probability (e.g., 'membership') values for calculation of KL divergence and creation of constellation diagrams Currently a subfolder for the child taxonomy is also created, but everything in that folder is also stored in the anndata, and so it can be safely ignored.

# At this point the reference taxonomy is created and ready for Patch-seq mapping!  It should be saved in your current directory in a file called "AI_taxonomy.h5ad"


## We also set up a taxonomy for only glutamatergic neurons. The reason for this is define a set of marker genes for patch-seq integration of glutamatergic neurons.  This may not be needed; it also may be useful to map cells to glutamatergic neuron types directly.

## Setup this taxonomy for patch-seq mapping and marker gene identification
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "glutamatergic",
                                    subsample = 100, ## Subsampling not needed in this case
                                    subclass.column = "subclass_label", 
                                    class.column = "class_label", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal and Non-neural","Neuronal: GABAergic"), ## The off-target class.column labels for patchseqQC.
                                    subclass.subsample = 100, ## Subsampling is for PatchseqQC contamination calculation.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = taxonomyDir) # This will create a subfolder in the reference taxonomy directory
# We note that the patch-seq QC parameters in this child taxonomy should NOT be used.




