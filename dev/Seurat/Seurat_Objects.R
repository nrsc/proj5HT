# Patch-Seq offpipeline Macaque ---------------------------------------------------------------
SeuratHCT_offPipelineMacaque = function(save_object = TRUE, #subclass = c('L5 IT','L5 ET',"L2/3 IT")
                                 subclass = c('L5 IT', 'L5 ET', "L5/6 NP", "L4 IT", "L2/3 IT", "L6 IT", "L6 CT", "L6b")
                                 ) {
  MD = projHCT::sheets$MD

  # Setup metadata from annotation file
  aPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Macaque_NCBI_RSC-204-386_map_full/", "anno.feather"))

  # Select only our cells
  tstAP = aPatch[match(MD$cell_name, aPatch$cell_name_label, nomatch = 0), ]
  #pipeAP = aPatch[!match(MD$cell_name, aPatch$cell_name_label, nomatch = 0),]

  # Pull those cells from the MD to merge
  tstMD = MD[match(tstAP$cell_name_label, MD$cell_name, nomatch = 0), ]

  # Check for matches
  tryCatch({
    identical(tstMD$cell_name, tstAP$cell_name_label)
  }, error = function(e) {
    warning("Metadata names and annotation names do not match", e$message)
    stop("Critical Failure: ", e$message)
  })

  tstMD$patched_cell_container_label = tstMD$patched_cell_container
  tstMD$patched_cell_container = NULL
  tstMD$sample_id = NULL

  # Metadata merging from our MD with the annotation data from the sequencing data.
  sMD = merge(tstAP, tstMD, "patched_cell_container_label")

  # Isolate subclass
  sMD = sMD[is.element(sMD$subclass_Corr_label, subclass), ]

  # Remove Squirrel Monkey Cells that are labeled with a Q instead of a SC
  sMD = sMD[-which(sMD$Species == "S.sciureus"),]

  # Load counts data. Data is already in RPM
  dPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Macaque_NCBI_RSC-204-386_map_full/", "data.feather"))

  # Select out metadata. Will select what's left in the sMD after cell and subclass filtering.
  dPatch = dPatch[match(sMD$sample_id, dPatch$sample_id), ]
  dPatch = as.data.frame(dPatch)

  #### Set rownames and colnames before transposing ####
  rownames(dPatch) = dPatch$sample_id
  dPatch$sample_id = NULL
  dPatch = as.data.frame(t(dPatch))
  #dPatch = as.data.frame(t(dPatch))
  #### Dp a log transformation of the data object
  dPatch = log2(dPatch + 1)

  sMD$set = rep("allMacaque-PatchSeq", dim(dPatch)[2])

  # dPatch.metadata =  data.frame(
  #   set = rep("Patch-seq", dim(dPatch)[2]),
  #   cell_name = sMD$cell_name,
  #   Cortical_area = sMD$Cortical_area,
  #   patched_cell_container = sMD$patched_cell_container_label,
  #   subclass = sMD$subclass_Corr_label,
  #   subclassTree = sMD$subclass_Tree_label,
  #   cluster = sMD$cluster_Corr_label,
  #   clusterTree = sMD$cluster_Tree_label,
  #   scoreCorr = sMD$score.Corr,
  #   scoreTree = sMD$score.Tree,
  #   praIntrons = sMD$percent_reads_aligned_to_introns,
  #   markerSumNorm = sMD$marker_sum_norm,
  #   roi = sMD$roi_label,
  #   species = sMD$Species,
  #   neighbourhoodCorr = sMD$neighborhood_Corr_label,
  #
  #   row.names = colnames(dPatch)
  # )


  dPatch <- Seurat::CreateSeuratObject(counts = dPatch,
                                       meta.data = sMD,
                                       project = "opMacaque-PatchSeq")

  dPatch[["percent.mt"]] <- PercentageFeatureSet(dPatch, pattern = "^MT*")

  if (save_object) {
    saveRDS(dPatch, file = "~/proj/proj5HT/data-raw/Seurat/offPipeline-Macaque_PatchSeq.rds")
  }

  ## Construct data set lists

  return(dPatch)

}

mopPS = SeuratHCT_offPipelineMacaque()

# Patch-Seq offpipeline Macaque ---------------------------------------------------------------
SeuratHCT_offPipelineHuman = function(save_object = TRUE, #subclass = c('L5 IT','L5 ET',"L2/3 IT")
                                      subclass = c('L5 IT', 'L5 ET', "L5/6 NP", "L4 IT", "L2/3 IT", "L6 IT", "L6 CT", "L6b")
                                      ) {
  MD = projHCT::sheets$MD

  # Setup metadata from annotation file
  aPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Human_RSC-122-380_map_full/", "anno.feather"))

  # Select only our cells
  tstAP = aPatch[match(MD$cell_name, aPatch$cell_name_label, nomatch = 0), ]
  #pipeAP = aPatch[!match(MD$cell_name, aPatch$cell_name_label, nomatch = 0),]

  # Pull those cells from the MD to merge
  tstMD = MD[match(tstAP$cell_name_label, MD$cell_name, nomatch = 0), ]
  #unique(tstMD$Cortical_area)
  tstMD[which(tstMD$Cortical_area == "MTG"),"Cortical_area"] = "TCx"


  # Check for matches
  tryCatch({
    identical(tstMD$cell_name, tstAP$cell_name_label)
  }, error = function(e) {
    warning("Metadata names and annotation names do not match", e$message)
    stop("Critical Failure: ", e$message)
  })

  tstMD$patched_cell_container_label = tstMD$patched_cell_container
  tstMD$patched_cell_container = NULL
  tstMD$sample_id = NULL

  # Metadata merging from our MD with the annotation data from the sequencing data.
  sMD = merge(tstAP, tstMD, "patched_cell_container_label")

  # Isolate subclass
  sMD = sMD[is.element(sMD$subclass_Corr_label, subclass), ]

  # Remove Squirrel Monkey Cells that are labeled with a Q instead of a SC
  #sMD = sMD[-which(sMD$Species == "S.sciureus"),]

  # Load counts data. Data is already in RPM
  dPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Human_RSC-122-380_map_full/", "data.feather"))

  # Select out metadata. Will select what's left in the sMD after cell and subclass filtering.
  dPatch = dPatch[match(sMD$sample_id, dPatch$sample_id), ]
  dPatch = as.data.frame(dPatch)

  #### Set rownames and colnames before transposing ####
  rownames(dPatch) = dPatch$sample_id
  dPatch$sample_id = NULL
  dPatch = as.data.frame(t(dPatch))
  #dPatch = as.data.frame(t(dPatch))
  #### Dp a log transformation of the data object
  dPatch = log2(dPatch + 1)

  sMD$set = rep("Human offp-Patch-seq", dim(dPatch)[2])

  dPatch <- Seurat::CreateSeuratObject(counts = dPatch,
                                       meta.data = sMD,
                                       project = "opHuman-PatchSeq")

  dPatch[["percent.mt"]] <- PercentageFeatureSet(dPatch, pattern = "^MT*")

  if (save_object) {
    saveRDS(dPatch, file = "~/proj/proj5HT/data-raw/Seurat/offPipeline-Human_PatchSeq.rds")
  }

  ## Construct data set lists

  return(dPatch)

}

hopPS = SeuratHCT_offPipelineHuman()

# all Patch-Seq Macaque ---------------------------------------------------------------
SeuratHCT_excPatchSeqMacaque = function(save_object = TRUE, #subclass = c('L5 IT','L5 ET',"L2/3 IT")
                              subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b")
                              ) {
  MD = projHCT::sheets$MD

  # Setup metadata from annotation file
  aPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Macaque_NCBI_RSC-204-386_map_full/", "anno.feather"))

  aPatch = aPatch[is.element(aPatch$subclass_Corr_label, subclass), ]

  #aPatch = aPatch[-which(aPatch$species_label == "S.sciureus"),]


  # Load counts data. Data is already in RPM
  dPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Macaque_NCBI_RSC-204-386_map_full/", "data.feather"))

  # Select out metadata. Will select what's left in the sMD after cell and subclass filtering.
  dPatch = dPatch[match(aPatch$sample_id, dPatch$sample_id), ]
  dPatch = as.data.frame(dPatch)

  #### Set rownames and colnames before transposing ####
  rownames(dPatch) = dPatch$sample_id
  dPatch$sample_id = NULL
  dPatch = as.data.frame(t(dPatch))
  #### Dp a log transformation of the data object
  dPatch = log2(dPatch + 1)

  aPatch$pipe_origin = dplyr::if_else(is.na(
    match(aPatch$cell_name_label, MD$cell_name, nomatch = NA)
  ), "on-pipeline", "off-pipeline")
  aPatch$set = "excMacaque-PatchSeq"


  dPatch <- Seurat::CreateSeuratObject(counts = dPatch,
                                       meta.data = aPatch,
                                       project = "excMacaque-PatchSeq")

  dPatch[["percent.mt"]] <- PercentageFeatureSet(dPatch, pattern = "^MT*")


  if (save_object) {
    saveRDS(dPatch, file = "~/proj/proj5HT/data-raw/Seurat/excMacaque_PatchSeq.rds")
  }

  ## Construct data set lists

  return(dPatch)

}
mePS = SeuratHCT_excPatchSeqMacaque()

# all Patch-Seq Human ---------------------------------------------------------------
SeuratHCT_excPatchSeqHuman = function(save_object = TRUE, #subclass = c('L5 IT','L5 ET',"L2/3 IT")
                              subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b")
                              ) {
  MD = projHCT::sheets$MD

  # Setup metadata from annotation file
  aPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Human_RSC-122-380_map_full/", "anno.feather"))

  aPatch = aPatch[is.element(aPatch$subclass_Corr_label, subclass), ]

  # Load counts data. Data is already in RPM
  dPatch <- feather::read_feather(paste0("~/proj/proj5HT/den/macaque/MTG/GreatApes_Human_RSC-122-380_map_full/", "data.feather"))

  # Select out metadata. Will select what's left in the sMD after cell and subclass filtering.
  dPatch = dPatch[match(aPatch$sample_id, dPatch$sample_id), ]
  dPatch = as.data.frame(dPatch)

  #### Set rownames and colnames before transposing ####
  rownames(dPatch) = dPatch$sample_id
  dPatch$sample_id = NULL
  dPatch = as.data.frame(t(dPatch))
  #### Dp a log transformation of the data object
  dPatch = log2(dPatch + 1)

  aPatch$pipe_origin = dplyr::if_else(is.na(
    match(aPatch$cell_name_label, MD$cell_name, nomatch = NA)
  ), "on-pipeline", "off-pipeline")

  aPatch$set = "excHuman-Patch-Seq"


  dPatch <- Seurat::CreateSeuratObject(counts = dPatch,
                                       meta.data = aPatch,
                                       project = "excHuman Patch-seq")

  dPatch[["percent.mt"]] <- PercentageFeatureSet(dPatch, pattern = "^MT*")


  if (save_object) {
    saveRDS(dPatch, file = "~/proj/proj5HT/data-raw/Seurat/excHuman_PatchSeq.rds")
  }

  ## Construct data set lists

  return(dPatch)

}

hePS = SeuratHCT_excPatchSeqHuman()


# 10xFacs and SmartSeq Macaque ---------------------------------------------------------------
SeuratHCT_10xSmart = function(pathSmart =  "~/proj/proj5HT/den/macaque/MTG/GreatApes_Macaque_NCBI_RSC-204-386_map_full/",
                              path_10x = "~/proj/proj5HT/den/GreatApes_Macaque_NCBI/",
                              save_object = TRUE,
                              save_file = "~/proj/proj5HT/data-raw/Seurat/Macaque_Smart10x.rds",
                              subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b","Chandelier","Pvalb")
                              ) {
  # 10x FACs ----------------------------------------------------------------

  # Annotation data
  aMTG <- feather::read_feather(paste(path_10x, "anno.feather", sep = ""))
  aMTG$Cortical_area = "TCx"

  # Select Subclass
  aMTG <- aMTG[is.element(aMTG$subclass_label, subclass), ]

  # Read in counts
  cMTG <-  feather::read_feather(paste(path_10x, "counts.feather", sep =""))
  # Feather loads in tibble. Change to dataframe
  cMTG = as.data.frame(cMTG)
  cMTG = e(cMTG)
  gene = cMTG$gene
  rownames(cMTG) = gene
  cMTG$gene = NULL

  # Select the cells, identified by column names, by those that are left in the annotation data
  cMTG = cMTG[, match(aMTG$sample_id, colnames(cMTG), nomatch = 0)]
  # Annotation data is filtered by subtype at this point.

  if (!identical(seq(1, length(names(cMTG))), match(colnames(cMTG), aMTG$sample_id))) {
    errorCondition(message = "Counts colnames and sample_id do not match")
  }

  cMTG = as.data.frame(cMTG)
  cMTG <- sweep(cMTG, 2, colSums(cMTG), FUN = "/") * 1e6
  #cMTG = cMTG * 10^6 / rep(colSums(cMTG), each = nrow(cMTG))
  cMTG = log2(cMTG + 1)

  cMTG.metadata =  data.frame(
    set = rep("FACs", dim(cMTG)[2]),
    pipe = rep("FACs", dim(cMTG)[2]),
    cell_name = rep("FACs", dim(cMTG)[2]),
    patched_cell_container = rep("FACs", dim(cMTG)[2]),
    subclassCorr = aMTG$subclass_label,
    subclassTree = aMTG$subclass_label,
    clusterCorr = aMTG$cluster_label,
    clusterTree = aMTG$cluster_label,
    roi = aMTG$roi_label,
    species = aMTG$species_label,
    neighbourhoodCorr = aMTG$neighborhood_label,
    row.names = colnames(cMTG)
  )

  gc()

  # Smart-Seq ---------------------------------------------------------------
  # Setup metadata from annotation file
  aPatch <- feather::read_feather(paste(pathSmart, "anno.feather", sep =""))
  # Isolate the subclass of interest
  aPatch = aPatch[is.element(aPatch$subclass_Tree_label, subclass), ]

  # Load counts data. Data is already in RPM
  dPatch <- feather::read_feather(paste0(pathSmart, "data.feather"))

  # Select out metadata. Will select what's left in the sMD after cell and subclass filtering.
  dPatch = dPatch[match(aPatch$sample_id, dPatch$sample_id), ]
  dPatch = as.data.frame(dPatch)

  #### Set rownames and colnames before transposing ####
  rownames(dPatch) = dPatch$sample_id
  dPatch$sample_id = NULL
  dPatch = as.data.frame(t(dPatch))
  #### Dp a log transformation of the data object
  dPatch = log2(dPatch + 1)

  pipe_origin = dplyr::if_else(is.na(match(
    aPatch$cell_name_label, MD$cell_name, nomatch = NA
  )), "on-pipeline", "off-pipeline")

  dPatch.metadata =  data.frame(
    set = rep("Smart-Seq", dim(dPatch)[2]),
    pipe = pipe_origin,
    cell_name = aPatch$cell_name_label,
    patched_cell_container = aPatch$patched_cell_container_label,
    subclassCorr = aPatch$subclass_Corr_label,
    subclassTree = aPatch$subclass_Tree_label,
    clusterCorr = aPatch$cluster_Corr_label,
    clusterTree = aPatch$cluster_Tree_label,
    roi = aPatch$roi_label,
    species = aPatch$species_label,
    neighbourhoodCorr = aPatch$neighborhood_Corr_label,
    row.names = colnames(dPatch)
  )


  # Combine the two datasets

  cMTG = cMTG[-which(is.na(match(rownames(cMTG), rownames(dPatch)))), ]
  dPatch = dPatch[-which(is.na(match(rownames(dPatch), rownames(cMTG)))), ]

  brain.data = cbind(cMTG, dPatch)

  brain.metadata = rbind(cMTG.metadata, dPatch.metadata)

  Smart10x <- Seurat::CreateSeuratObject(counts = brain.data,
                                         meta.data = brain.metadata,
                                         project = "Smart10x")

  Smart10x[["percent.mt"]] <- PercentageFeatureSet(Smart10x, pattern = "^MT*")

  if (save_object) {
    saveRDS(Smart10x, file = save_file)
  }

  return(Smart10x)

}


Smart10x = SeuratHCT_10xSmart()


# 10xFACs Macaque ---------------------------------------------------------------
SeuratHCT_10xMTG = function(path_10x = "~/proj/proj5HT/den/GreatApes_Macaque_NCBI/",
                            save_object = TRUE,
                            save_file = "~/proj/proj5HT/data-raw/MacaqueMTG_10xFACs.rds",
                            subclass = c('L5 IT','L5 ET',"L2/3 IT")
                            #subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b","Chandelier","Pvalb")
) {
  #### Feather for 10x data
  # Annotation data
  aMTG <- feather::read_feather(paste(path_10x, "anno.feather", sep = ""))

  # Select Subclass
  aMTG <- aMTG[is.element(aMTG$subclass_label, subclass), ]

  # Read in counts
  cMTG <-  feather::read_feather(paste(path_10x, "counts.feather", sep = ""))

  # Feather loads in tibble. Change to dataframe
  cMTG = as.data.frame(cMTG)
  gene = cMTG$gene
  rownames(cMTG) = gene
  cMTG$gene = NULL

  # Select the cells, identified by column names, by those that are left in the annotation data
  cMTG = cMTG[, match(aMTG$sample_id, colnames(cMTG), nomatch = 0)]
  # Annotation data is filtered by subtype at this point.

  if (!identical(seq(1, length(names(cMTG))), match(colnames(cMTG), aMTG$sample_id))) {
    errorCondition(message = "Counts colnames and sample_id do not match")
  }

  cMTG = as.data.frame(cMTG)

  # In counts so need to transform to Reads per million
  cMTG = cMTG * 10^6 / rep(colSums(cMTG), each = nrow(cMTG))

  cMTG = log2(cMTG + 1)

  ## Construct data set lists
  cMTG <- Seurat::CreateSeuratObject(counts = cMTG,
                                     meta.data = aMTG,
                                     project = "cMTG")

  cMTG[["percent.mt"]] <- PercentageFeatureSet(cMTG, pattern = "^MT*")

  if (save_object) {
    saveRDS(cMTG, file = save_file)
  }

  return(cMTG)

}

# 10x Human ---------------------------------------------------------------
SeuratHCT_10xHuman = function(human_10x = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Human/",
                              save_object = TRUE,
                              save_file = "~/proj/proj5HT/data-raw/Seurat/HumanMTG_10xFACs.rds",
                              subclass = c('L5 IT','L5 ET',"L2/3 IT")
                              #subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b")
){

  aMTG <- feather::read_feather(paste("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Human/", "anno.feather", sep = ""))
  aMTG$Cortical_area = "TCx"

  # Select Subclass
  aMTG <- aMTG[is.element(aMTG$subclass_label, subclass), ]

  # Read in counts
  cMTG <-  feather::read_feather(paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Human/", "counts.feather"))

  cMTG = as.data.frame(cMTG)
  gene = cMTG$gene
  rownames(cMTG) = gene
  cMTG$gene = NULL

  # Select the cells, identified by column names, by those that are left in the annotation data
  cMTG = cMTG[, match(aMTG$sample_id, colnames(cMTG), nomatch = 0)]
  # Annotation data is filtered by subtype at this point.

  #saveRDS(cMTG, "human_cMTG.rds")
  #cMTG = readRDS("human_cMTG.rds")

  if (!identical(seq(1, length(names(cMTG))), match(colnames(cMTG), aMTG$sample_id))) {
    errorCondition(message = "Counts colnames and sample_id do not match")
  }

  cMTG = as.data.frame(cMTG)
  cMTG = cMTG * 10^6 / rep(colSums(cMTG), each = nrow(cMTG))

  cMTG = log2(cMTG + 1)

  metadata.cMTG <- data.frame(
    set = rep("FACs-Human", dim(cMTG)[2]),
    Cortical_area = "MTG",
    subclass = aMTG$subclass_label,
    cluster = aMTG$cluster_label,
    crossSpeciesCluster = aMTG$cross_species_cluster_label,
    neighbourhood = aMTG$neighborhood_label,
    row.names = colnames(cMTG)
  )

  ## Construct data set lists
  cMTG <- Seurat::CreateSeuratObject(counts = cMTG,
                                     meta.data = aMTG,
                                     project = "cMTG")

  cMTG[["percent.mt"]] <- PercentageFeatureSet(cMTG, pattern = "^MT*")

  if (save_object) {
    saveRDS(cMTG, file = "~/proj/proj5HT/data-raw/Seurat/HumanMTG_10xFACs.rds")
  }

  return(cMTG)

}


SeuratHCT_10xChimp = function(path10x = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Chimp//",
                              save_object = TRUE,
                              save_file = "~/proj/proj5HT/data-raw/Seurat/ChimpMTG_10xFACs.rds",
                              subclass = c('L5 IT','L5 ET',"L2/3 IT")
                              #subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b")
){

  aMTG <- feather::read_feather(paste("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Chimp/", "anno.feather", sep = ""))
  aMTG$Cortical_area = "TCx"

  # Select Subclass
  aMTG <- aMTG[is.element(aMTG$subclass_label, subclass), ]

  # Read in counts
  cMTG <-  feather::read_feather(paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Chimp/", "counts.feather"))

  cMTG = as.data.frame(cMTG)
  gene = cMTG$gene
  rownames(cMTG) = gene
  cMTG$gene = NULL

  # Select the cells, identified by column names, by those that are left in the annotation data
  cMTG = cMTG[, match(aMTG$sample_id, colnames(cMTG), nomatch = 0)]
  # Annotation data is filtered by subtype at this point.

  #saveRDS(cMTG, "human_cMTG.rds")
  #cMTG = readRDS("human_cMTG.rds")

  if (!identical(seq(1, length(names(cMTG))), match(colnames(cMTG), aMTG$sample_id))) {
    errorCondition(message = "Counts colnames and sample_id do not match")
  }

  cMTG = as.data.frame(cMTG)
  cMTG = cMTG * 10^6 / rep(colSums(cMTG), each = nrow(cMTG))

  cMTG = log2(cMTG + 1)

  metadata.cMTG <- data.frame(
    set = rep("FACs-Human", dim(cMTG)[2]),
    Cortical_area = "MTG",
    subclass = aMTG$subclass_label,
    cluster = aMTG$cluster_label,
    crossSpeciesCluster = aMTG$cross_species_cluster_label,
    neighbourhood = aMTG$neighborhood_label,
    row.names = colnames(cMTG)
  )

  ## Construct data set lists
  cMTG <- Seurat::CreateSeuratObject(counts = cMTG,
                                     meta.data = aMTG,
                                     project = "chimp10x")

  cMTG[["percent.mt"]] <- PercentageFeatureSet(cMTG, pattern = "^MT*")

  if (save_object) {
    saveRDS(cMTG, file = "~/proj/proj5HT/data-raw/Seurat/ChimpMTG_10xFACs.rds")
  }

  return(cMTG)

}

SeuratHCT_10xGorilla = function(path10x = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Gorilla/",
                              save_object = TRUE,
                              save_file = "~/proj/proj5HT/data-raw/Seurat/GorillaMTG_10xFACs.rds",
                              subclass = c('L5 IT','L5 ET',"L2/3 IT")
                              #subclass = c('L5 IT','L5 ET',"L5/6 NP","L4 IT","L2/3 IT","L6 IT","L6 CT","L6b")
){

  aMTG <- feather::read_feather(paste(path10x, "anno.feather", sep = ""))
  aMTG$Cortical_area = "TCx"

  # Select Subclass
  aMTG <- aMTG[is.element(aMTG$subclass_label, subclass), ]

  # Read in counts
  gMTG <-  feather::read_feather(paste0(path10x, "counts.feather"))

  gMTG = as.data.frame(gMTG)
  gene = gMTG$gene
  rownames(gMTG) = gene
  gMTG$gene = NULL

  # Select the cells, identified by column names, by those that are left in the annotation data
  gMTG = gMTG[, match(aMTG$sample_id, colnames(gMTG), nomatch = 0)]
  # Annotation data is filtered by subtype at this point.

  #saveRDS(gMTG, "gorilla_gMTG.rds")
  gMTG = readRDS("gorilla_cMTG.rds")

  if (!identical(seq(1, length(names(gMTG))), match(colnames(gMTG), aMTG$sample_id))) {
    errorCondition(message = "Counts colnames and sample_id do not match")
  }

  gMTG = as.data.frame(gMTG)
  gMTG = gMTG * 10^6 / rep(colSums(gMTG), each = nrow(gMTG))

  gMTG = log2(gMTG + 1)

  ## Construct data set lists

  gMTG <- Seurat::CreateSeuratObject(counts = gMTG,
                                     meta.data = aMTG,
                                     project = "Gorilla10x")

  gMTG[["percent.mt"]] <- PercentageFeatureSet(gMTG, pattern = "^MT*")

  if (save_object) {
    saveRDS(gMTG, file = "~/proj/proj5HT/data-raw/Seurat/GorillaMTG_10xFACs.rds")
  }

  return(gMTG)

}
