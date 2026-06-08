# `htr_scores.csv` / `htr_scores.rds`

Per-cell HTR receptor expression table. Built by
[`build_htr_scores.R`](build_htr_scores.R) from the macaque PatchSeq Seurat
object at `data-raw/Seurat/offPipeline-Macaque_PatchSeq.rds`.

One row per cell (`cell_name` key, joins to `compile_srtPuff.csv$cell_name`).

## Columns

**Metadata (subset of Seurat `@meta.data`):**
`cell_name`, `subclass_Corr`, `subclass_Tree`, `predicted_subclass`,
`Cortical_area`, `Species`, `LIMS_depth`, `Depth_from_pia`,
`cluster_Corr_label`, plus QC columns (`Genes.Detected`,
`percent_reads_aligned_to_introns`, `marker_sum_norm_label`,
`quality_score_label`, `score.Corr_label`).

**Receptor expression (counts-per-million, CPM, from Seurat `counts` layer):**
HTR1A, HTR1B, HTR1D, HTR1E, HTR1F, HTR2A, HTR2B, HTR2C, HTR4, HTR5A, HTR6, HTR7.

> The Seurat object's `counts` layer is curated as CPM — already per-cell
> normalized. Do **not** apply `NormalizeData()` before extraction.

**Provisional-source flags:**
`htr1a_status`, `htr1b_status` — currently both set to
`"provisional_seurat_primary"`. Scott will overlay corrected values from
accessory files in a follow-up pass; after patching, update these tags to
the accessory-file identifier (e.g. `"patched_<source>_<date>"`).

**Derived family scores (NA-tolerant rowSums):**
- `HTR1_sum` — Σ HTR1A/B/D/E/F (inhibitory family, Gi)
- `HTR2_sum` — Σ HTR2A/B/C (excitatory family, Gq)
- `HTR0_sum` — Σ HTR4/5A/6/7 (mixed Gs/Gi)
- `HTR_balance` — `HTR2_sum − HTR1_sum`. Positive = excitatory bias.

The derived columns are recomputed every run, so re-running the build
script after patching HTR1A/HTR1B will refresh `HTR1_sum` and `HTR_balance`
automatically.

## Status

**Provisional** — HTR1A and HTR1B values pending accessory-file patch.
All other receptors are final from the primary Seurat object.
