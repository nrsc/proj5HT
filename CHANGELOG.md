# proj5HT — Changelog

> Completed work log. Planned items live in [TODO.md](TODO.md).
> Newest first.

## Unreleased

### 2026-06-04 — hier metadata normalization + Seurat figure refactor

- **New `R/hier_metadata5HT.R` helpers** (all exported):
  `normalize_label_cols(df, cols)` — single source of truth for
  taxonomy/hier label cleanup (space → underscore, empty → NA, idempotent).
  `normalize_hier_assignments(df)` and `normalize_taxonomy_labels(df)`
  are thin wrappers with sensible defaults (4 `hier_*_assignment` cols
  and 8 AIBS taxonomy cols, respectively).
- **`add_hier_to_seurat(obj)`** now drops any pre-existing `hier_*`
  columns before merging, so re-runs are idempotent.
- **`build_seurat_from_feather()`** (`dev/Seurat/Seurat_Objects.R`)
  threads `normalize_taxonomy_labels()` unconditionally at build time
  and `normalize_hier_assignments()` after the optional `cleanMD` join.
  Dead `attach_hier` parameter removed.
- **New `rmd/HTR_violin_dotplots_hier.Rmd`** consolidates and replaces
  the legacy `Seurat-*.Rmd` plotting set. Introduces a
  `resolve_subclass()` fallback chain
  (`hier_subclass_assignment` → `subclass_Corr_label` → `subclass_label`)
  with per-cell `subclass_source` tagging and a per-cohort kable so
  panels declare exactly which call sourced each plot.
- **New vignette `vignettes/hier-metadata.Rmd`** documents the refactor
  and walks the `build_htr_scores` flow.
- **NAMESPACE**: added exports
  `add_hier_to_seurat`, `hier5HT`,
  `normalize_hier_assignments`,
  `normalize_label_cols`,
  `normalize_taxonomy_labels`.

### 2026-06-04 — Triage: 194 PatchSeq cells missing hier metadata

- Diagnosed root cause for 194 / 839 macaque PatchSeq cells being
  unmatched by `add_hier_to_seurat()`: all 194 are present in
  `headBCI::sheets$MD` but filtered out of `sheets$cleanMD` because
  `nwb_file == NA` (no NWB discovered under `den/`).
- User breakdown (filtered, n=194): Nikolai 119, Alejandro 23,
  Jeffrey 14, Brian 12, Meanhwan 8, Cristina 7, Kamilliam 3, ScottS 3,
  Lindsay 2, ScottO 2, Jonathan 1. The cleanMD set is dominated by
  MIES users (Brian, Cristina, ScottS, ScottO), confirming the gap
  tracks non-MIES acquisition workflows.
- Wrote tracking input
  [`dev/filtered_from_cleanMD.csv`](dev/filtered_from_cleanMD.csv)
  (194 rows × 11 cols, sorted User → Date → cell_name) for downstream
  cell-by-cell follow-up.
- Spawned two new headBCI/TODO items in section A:
  (1) non-MIES acquisition workflow support,
  (2) multipatch sort/categorize pipeline.

### 2026-06-02 — Deprecations

- **Archived `R/ketWashIn.R`** → [`.archive/2026-06-02_ketwash_spikepuff/`](.archive/2026-06-02_ketwash_spikepuff/).
  The newer `DrugWashIn()` (and the canonical pipeline path
  `DrugWashIn_from_asp()`) already handle both `Ket-wash-in` and
  `WAYwash` experiments (`R/DrugWashIn.R` line 47). No callers remained.
- **Archived `R/spikePuff.R`** → same dated archive folder. It duplicated
  the heavy NWB-extract + APanalysis work already performed by
  `build_asp()`; the canonical pipeline path is
  `spikePuff_from_asp()` (called in `run_asp_dfs_pipeline()`).
  `spikePuff()` was not exported and had no live callers.
- Removed corresponding `man/ketWashIn.Rd`; dropped `export(ketWashIn)`
  from `NAMESPACE`.
- README + `manu/README.md` updated to reflect the slimmer module set.

### 2026-06-02 — Project hygiene

- Added `CHANGELOG.md` (this file) and `TODO.md` as the canonical
  cross-session work log for proj5HT (mirroring the headBCI convention).

## Earlier (recovered from git)

- **2026-05-12** `9a38b11` — Fix classifier blind spots for silent periods
  in sparse traces.
- **2026-05-12** `0a91272` — Consolidate figure pipeline, add 4-way
  classifier, archive dead code.
- **2026-04-02** `5447d48` — Updates to infrastructure.
- **2026-03-25** `3aa6cfd` — README update.
- **2026-03-13** `0a59694` — First push in a long time (post-monolith
  decomposition).
