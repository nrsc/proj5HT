# proj5HT — Changelog

> Completed work log. Planned items live in [TODO.md](TODO.md).
> Newest first.

## Unreleased

### 2026-06-08 — migrate parent nnest from projHCT (-hct.rds) to headBCI (-bci.rds)

proj5HT no longer depends on projHCT for the parent-nnest lifecycle.
The child payload (`-srt.rds`) is now built off the headBCI parent
(`-bci.rds`) via `headBCI::bci_load()` / `headBCI::bci_nnest()`.

- **`nnest5HT()`** — parent loader switched from
  `projHCT::loadHCT(cell, "-hct.rds")` /
  `projHCT::nnestHCT(cell)` to
  `headBCI::bci_load(cell, parent_tag)` /
  `headBCI::bci_nnest(cell)`, where
  `parent_tag = headBCI::paramsBCI$tags$nnest` (`-bci.rds`). Local var
  renamed `hct` → `parent`.
- **`load5HT()`** — slow-rehydrate path now uses
  `headBCI::bci_load(cell, parent_tag)` and writes `parent_tag` into
  `child$parent$tag`. The no-child branch bootstraps via
  `nnest5HT(cell)` (matching the fast-path behavior), fixing a silent
  `NULL` return that blocked `update_5ht` for fresh cells.
  Default `overwrite` now includes both `"ver"` (headBCI) and
  `"version"` (legacy projHCT) header keys during the transition.
- **`hct_header_fields()`** — accepts both `version` and `ver`.
- **Path resolvers consolidated:**
  - `find_hct_path_fast()` deleted (duplicate of `find_child_path_fast`).
  - `resolve_hct_path()` renamed to `resolve_parent_path()`,
    parameterized on `parent_tag`, source defaults to
    `headBCI::paramsBCI$paths$rookery`.
  - `get_parent_path_from_child()` parameterized on `parent_tag`.
  - `find_child_path_fast()` source defaults to
    `headBCI::paramsBCI$paths$rookery` (the resolved local path, not the
    network token form stored in `paramsBCI$rookery`).
- **`get_idx_map()`** — switched to
  `headBCI::sheets$files$nnest_files` for `-bci.rds`. Child tags
  (`-srt.rds`, `-osc.rds`) fall back to disk scan (no projHCT index
  yet at the BCI level).
- **`params5HT$tags$parent`** — `"-hct.rds"` → `"-bci.rds"`
  (regenerated `data/params5HT.rda`).
- **`run_asp_pipeline()`** — added full Roxygen block and `@export`
  so the `update_5ht` startup driver can dispatch it.
- **`analysisNaAvail.R`** — `projHCT::loadHCT(x, tag = "-srt.rds")`
  → `load5HT(x, tag = "-srt.rds", rehydrate = FALSE)`.

Verified end-to-end: `~/proj5HT/dev/startup.sh --no-git` processes all
209 Master_UI cells (`ok=209 failed=0`); no `projHCT::` references
remain in the nnest lifecycle path. Remaining `projHCT::` references in
`R/` are scoped to NHP annotation tables in `data5HT.R`
(`sheets$map$mapCellsNHP`, `sheets$map$NHP_anno`) — separate migration
follow-up.

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
