# proj5HT — TODO

> Planned work. Completed items move to [CHANGELOG.md](CHANGELOG.md).
> Themed sections, roughly ordered by current priority.

---

## A. Today's focus (2026-06-02)

- [ ] **RMP–spike correlation work.** Cross-reference `build_rmp()` /
      `analyze_rmp_protocol_set()` output with spikePuff response
      classification to test whether baseline RMP predicts response
      direction / magnitude. Likely a new function in
      `R/response_metrics5HT.R` or a sibling
      `R/rmp_spike_concordance5HT.R`.
- [ ] **Baseline validation ideas** — extend `response_metrics5HT.R`
      with additional baseline sanity checks beyond drift / CV / Fano:
      - Stationarity test (e.g. ADF or KPSS on baseline rate)
      - Block-bootstrap CI on `blMean` to flag noisy baselines
      - Cross-protocol baseline reproducibility (compare blMean across
        sibling spikeTTL protocols in the same cell)
- [ ] **Seurat object hygiene** — exercise the refactored
      `dev/Seurat/Seurat_Objects.R` (see section D). Confirm each cohort
      builds, then attach `hier5HT()` via `add_hier_to_seurat()` and
      sanity-check the join coverage.

## B. Nnesting pipeline rerun

Goal: a clean, reproducible end-to-end nnest rebuild across the whole
rookery, validating the post-archive code surface.

- [ ] Decide the canonical entry point: `run_asp_pipeline(cell)` per
      cell, or `nnest5HT(cell) %>% build_asp() %>% build_rmp() %>%
      compile5HT()` chained in a per-cell driver.
- [ ] Inventory the tags currently written per cell:
      `-srt.rds`, `-srtPuff.csv`, `-srtDrugWash.csv`, `-rmp.csv`,
      `-asp.csv`. Confirm none reference the archived `ketWashIn`
      output (`-srtKetWash.csv`) — and sweep the rookery to delete any
      stale `-srtKetWash.csv` files.
- [ ] Add a `rerun_all_cells()` driver under `dev/` that:
      1. Lists candidate cells from `headBCI::sheets$MD`.
      2. Skips cells with no spikeTTL/TTLstack rows in the Master UI.
      3. Logs per-cell timing + success/fail to a manifest CSV.
      4. Optionally parallelises via `future.apply::future_lapply()`.
- [ ] Diff a sample of pre- vs post-rerun nnests (`-srt.rds`) to confirm
      the only changes are intended (drug-wash unification, no schema
      drift in `spike_puff_output`).
- [ ] Update `compile5HT()` if needed once outputs are uniform.

## C. DuckDB layer for proj5HT

Mirror the headBCI summary DB pattern (`bci_summary_*`, `R/summaryBCI.R`
in headBCI) but scoped to 5HT outputs.

- [ ] Create `R/summary5HT.R` exporting `srt_summary_build()`,
      `srt_summary_load(table)`, `srt_summary_update(srt)`, `srt_db()`.
- [ ] Tables (parquet, under `<rookery>/.summary5HT/`):
      - `cells` — one row per cell (subclass, hier_*, species, depth, …)
      - `spike_puff` — long-form rows from `-srtPuff.csv` outputs
      - `drug_wash` — long-form rows from `-srtDrugWash.csv`
      - `rmp` — rows from `-rmp.csv`
      - `protocol_sets` — per-protocol summary from `asp$smry`
      - `provenance` — last-modified, code version, params hash per cell
- [ ] Hook `srt_summary_update()` into `compile5HT()` / save5HT() so the
      store stays fresh without a full rebuild.
- [ ] Suggests: `arrow`, `duckdb`, `DBI`.
- [ ] Validate against the existing CSV compile (`compile5HT()` should
      produce the same `srtPuff` frame as `srt_db() %>% tbl("spike_puff")
      %>% collect()`).

## D. Seurat object cleanup

- [x] **2026-06-02** — drafted refactored builder at
      [`dev/Seurat/Seurat_Objects.R`](dev/Seurat/Seurat_Objects.R) that
      collapses the 8 cohort-specific functions into a single
      `build_seurat_from_feather()` core. Each cohort wrapper is now
      ~10 lines. Hier metadata is attached at build time via
      `add_hier_to_seurat()` when `cell_name_label` is present.
- [ ] Run the rebuild end-to-end and verify per-cohort row counts
      against the legacy outputs.
- [ ] Move the driver script (`mopPS = ...` etc.) out of source-time
      execution and into an explicit `build_all_seurat_objects()` entry
      point.
- [ ] Decide where the assembled Seurat .rds files live long-term — keep
      under `data-raw/Seurat/` or promote to lazy-loaded `data/`?
- [ ] Vignette: `vignettes/seurat-cohorts.Rmd` — show the cohort build
      + hier overlay + a DotPlot panel (mirrors Fig 2A in `manu/`).

## E. Manuscript figure scaffolding (manu/)

- [ ] Wire `manu/_shared/setup.R` against the rerun-validated outputs
      from section B (so figures don't silently inherit stale CSVs).
- [ ] Figure 2 (HTR landscape + hier classifier) — first real build
      target; uses `add_hier_to_seurat()` + the cleaned `offPipeline-*`
      Seurat objects.
- [ ] Figure 3 (response inventory) — depends on the response-metric
      additions from section A.

## F. Cross-project alignment

- [x] **2026-06-08** — Migrated proj5HT parent nnest from projHCT
      (`-hct.rds`) to headBCI (`-bci.rds`). `nnest5HT` / `load5HT` /
      path resolvers no longer reference `projHCT::`. See CHANGELOG
      `2026-06-08 — migrate parent nnest...`.
- [ ] Migrate the remaining `projHCT::sheets$map$mapCellsNHP` and
      `projHCT::sheets$map$NHP_anno` references in `R/data5HT.R` to
      `headBCI::sheets` equivalents (or absorb the NHP annotation
      mapping into headBCI). Track the 4-cell row delta noted in
      `manu/README.md` once consistent.
- [ ] Rename helpers that still carry "hct" in their names
      (`hct_header_fields`, `hct_header`, `apply_hct_header`) to
      `parent_header_*`. Mechanical rename; ripples to every call site
      in `child_payload5HT.R`.
- [ ] Index `-srt.rds` / `-osc.rds` files at the headBCI level
      (`headBCI::sheets$files$nnest_srt_files` etc.) so `get_idx_map`
      doesn't fall back to disk scan for child tags. headBCI-side
      change.
- [ ] Consider registering proj5HT as a `bci_register_project()` user
      and converting `build_asp` / `build_rmp` / `DrugWashIn` to
      `bci_register_step()` steps — final step of the proj5HT →
      framework integration noted in `headBCI/AGENTS.md`.

## G. Tests

- [ ] Stub `tests/testthat/test-DrugWashIn.R` — confirms one Ket-wash-in
      and one WAY-wash cell produce a non-empty `spike_puff_output`
      with the expected columns.
- [ ] Stub `tests/testthat/test-hier5HT.R` — `hier5HT()` returns one
      row per cell, no duplicate `cell_name`, "L2/3 IT" → "L2/3_IT"
      normalisation applied.
