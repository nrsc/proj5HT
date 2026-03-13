#subclass_col_brain <- "subclass_Corr"
subclass_col_brain <- "subclass"
subclass_col_cMTG  <- "subclass"

assay_use <- "RNA"
slot_detect <- "counts"   # for detection
slot_expr   <- "counts"   # for CPM

expr_thr <- 0


get_subclass_summary <- function(obj, genes, subclass_col,
                                 assay = "RNA",
                                 slot_counts = "counts",
                                 detect_thr = 0,
                                 dataset_name = "brain") {

  DefaultAssay(obj) <- assay

  genes_use <- intersect(genes, rownames(obj))

  counts <- GetAssayData(obj, slot = slot_counts)[genes_use, , drop = FALSE]
  subclass <- obj[[subclass_col]][,1]

  libsize <- Matrix::colSums(GetAssayData(obj, slot = slot_counts))

  df <- lapply(sort(unique(subclass)), function(sc) {

    cells <- colnames(obj)[subclass == sc]

    mat <- counts[, cells, drop = FALSE]

    # % expressing
    pct <- Matrix::rowMeans(mat > detect_thr) * 100

    # CPM (pseudobulk)
    pb_counts <- Matrix::rowSums(mat)
    pb_lib    <- sum(libsize[cells])
    cpm <- (pb_counts / pb_lib) * 1e6

    tibble::tibble(
      subclass = sc,
      gene     = genes_use,
      pct      = pct,
      cpm      = cpm,
      dataset  = dataset_name,
      n_cells  = length(cells)
    )
  }) %>% dplyr::bind_rows()

  df
}

df_brain <- get_subclass_summary(brain, HTRgenesINT, subclass_col_brain, dataset_name = "brain")
df_cMTG  <- get_subclass_summary(cMTG,  HTRgenesINT, subclass_col_cMTG,  dataset_name = "cMTG")

common_subclasses <- intersect(df_brain$subclass, df_cMTG$subclass)

df_brain <- dplyr::filter(df_brain, subclass %in% common_subclasses)
df_cMTG  <- dplyr::filter(df_cMTG,  subclass %in% common_subclasses)

df_diag <- df_brain %>%
  dplyr::select(subclass, gene, pct_brain = pct, cpm_brain = cpm, n_cells_brain = n_cells) %>%
  dplyr::left_join(
    df_cMTG %>%
      dplyr::select(subclass, gene, pct_cMTG = pct, cpm_cMTG = cpm, n_cells_cMTG = n_cells),
    by = c("subclass", "gene")
  )



### Circle map ###
ggplot2::ggplot(df_diag,
                ggplot2::aes(gene, subclass)) +
  ggplot2::geom_point(
    ggplot2::aes(size = pct_brain, color = log1p(cpm_brain)-log1p(cpm_cMTG))
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = pct_cMTG),
    shape = 1, stroke = 1
  ) +
  ggplot2::scale_color_viridis_c() +
  ggplot2::labs(
    title = "HTR expression: Patch-seq (filled) vs 10x (open)",
    size = "% expressing",
    color = "log1p(CPM)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


#### Delta circle map
df_diag <- df_diag %>%
  dplyr::mutate(delta_log1p_cpm = log1p(cpm_brain) - log1p(cpm_cMTG))

lim <- max(abs(df_diag$delta_log1p_cpm), na.rm = TRUE)  # 2.3100432

ggDOTs = ggplot2::ggplot(df_diag, ggplot2::aes(gene, subclass)) +
  ggplot2::geom_point(ggplot2::aes(size = pct_brain, color = delta_log1p_cpm)) +
  ggplot2::geom_point(ggplot2::aes(size = pct_cMTG), shape = 1, stroke = 1) +
  ggplot2::scale_color_gradient2(
    limits = c(-lim, lim),
    midpoint = 0,
    low = "yellow",
    mid = "sienna",
    high = "firebrick1",
    oob = scales::squish
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "HTR expression: Patch-seq (filled) vs 10x (open)",
    size = "% expressing",
    color = "Δ log1p(CPM)\n(pSeq - 10x)"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 axis.text = element_text(size = 16),
                 plot.title = element_text(size = 20),
                 axis.title = element_blank())


### Heat map
df_diag <- df_diag %>%
  dplyr::mutate(agreement = 1 - abs(pct_brain - pct_cMTG) / 100)

ggHM = ggplot2::ggplot(df_diag,
                ggplot2::aes(gene, subclass, fill = agreement)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    limits = c(0, 1),
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 axis.text = element_text(size = 16),
                 axis.title = element_blank())


grid.arrange(ggDOTs,ggHM, ncol = 1, layout_matrix = cbind(c(1,1,2)))






### % expressed
ggplot2::ggplot(df_diag,
                ggplot2::aes(pct_brain, pct_cMTG, label = gene)) +
  ggplot2::geom_abline(lty = 2) +
  ggplot2::geom_point() +
  ggplot2::geom_text(size = 3, vjust = -0.5) +
  ggplot2::facet_wrap(~subclass, nrow = 3) +
  ggplot2::labs(
    x = "% expressing (Patch-seq)",
    y = "% expressing (10x)"
  )# +
  #ggplot2::coord_equal()

df_diag %>%
  dplyr::group_by(subclass) %>%
  dplyr::summarise(r = cor(pct_brain, pct_cMTG, method = "spearman")) %>%
  ggplot2::ggplot(ggplot2::aes(reorder(subclass, r), r)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip()


#### diago ###
pal_HTR1  <- colorRampPalette(c("#08306B", "#4292C6"))(length(HTR1))
pal_HTR2  <- colorRampPalette(c("#00441B", "#41AB5D"))(length(HTR2))
pal_HTR0 <- colorRampPalette(c("#7F2704", "#F16913"))(length(HTR0))

gene_cols <- c(
  setNames(pal_HTR1,  HTR1),
  setNames(pal_HTR2,  HTR2),
  setNames(pal_HTR0, HTR0)
)



ggplot2::ggplot(
  df_diag,
  ggplot2::aes(pct_brain, pct_cMTG, label = gene, color = gene)
) +
  ggplot2::geom_abline(lty = 2) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = gene_cols) +
  ggplot2::scale_x_continuous(
    limits = c(-10, 100),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(-10, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  ggrepel::geom_text_repel(
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  ggplot2::facet_wrap(~subclass, nrow = 3) +
  ggplot2::labs(
    x = "% expressing (Patch-seq)",
    y = "% expressing (10x)",

  ) +
  theme_minimal()+
  ggplot2::coord_equal()+
  ggplot2::theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 18),
                 legend.text = element_text(size = 12))





