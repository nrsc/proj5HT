---

## Repo cleanup & improvement roadmap

The sections below document concrete problems, their causes, and
recommended fixes. Items are ordered roughly from highest to lowest
impact.

---

### 1. The classifier pipeline has sprawled — consolidate it

The 5-HT response classifier was built iteratively with help from
another algorithm and now lives in **four overlapping locations**, each
containing its own copies of the core functions:

| Location | Key files | What it contains |
|---|---|---|
| `dev/classify_pct_change/` (20 files) | `df_prep_and_per_cell_build.R`, `df_prep_new.R`, `classify_pct_change.R`, plus 17 one-off visualization scripts | The most complete version: `prep_puff_df()`, `extract_bucket_features_one()`, `extract_bucket_features()`, `add_biphasic_features()`, `make_plot_df()`, `make_plot_df_ei()`, `plot_bucket_qc()`, etc. |
| `dev/build_response_classifier/` (4 files) | `puff_helpers.R`, `puff_processing.R`, `puff_plotting.R`, `run_puff_classifier.R` | A cleaner rewrite with `prepare_puff_df()`, `extract_cell_features_one()`, `build_per_cell_classifier()`, `plot_cell_qc()`, `plot_response_composition()`. |
| `dev/gpt_build_2/` | `df_build.R` | Yet another `prep_puff_df_trace()` variant. |
| `dev/updated_plotting_script/` (3 files) | `updated_plotting_classifier.R`, `dr_pre_square_plots.R` | Another `prep_puff_df()` copy, another `extract_bucket_features()`. |
| `claude/checking_consistency.r` (1688 lines) | Monolithic copy | A full dump of the `df_prep_and_per_cell_build.R` logic. |

**The actual pipeline** that matters is:

```
compile_srtPuff.csv (compiled per-cell CSVs)
        │
        ▼
   df (master dataframe — time, instRate, cell_name, assigned_subclass, bath, puff, …)
        │  filter by bath/subclass/puff
        ▼
   df_prepped (prep_puff_df: adds baseline_hz, pct_of_baseline per cell)
        │
        ├──► extract_bucket_features() ──► per_cell  (rapid/mid/late window features)
        │                                      │
        │                                      ├──► add_biphasic_features() ──► per_cell2
        │                                      │         (response_type, biphasic_order)
        │                                      │
        │                                      └──► plot_bucket_qc() / plot_bucket_qc_many()
        │
        ├──► make_plot_df(df_prepped, per_cell)  ──► plot output (trace panels by modality)
        │
        └──► make_plot_df_ei(df_prepped, per_cell2) ──► exc/inh trace + mean plots
```

**Recommended fix:**

1. **Promote the authoritative functions** into `R/` as proper exported
   functions. Create two new files:
   - `R/classify5HT.R` — `prep_puff_df()`, `extract_bucket_features()`,
     `add_biphasic_features()`, `classify_puff_response()` (a single
     top-level wrapper that runs the full chain and returns a list with
     `df_prepped`, `per_cell`, `per_cell2`).
   - `R/plot_classify5HT.R` — `plot_bucket_qc()`, `make_plot_df()`,
     `make_plot_df_ei()`, and the various summary visuals
     (`plot_response_composition`, `plot_biphasic_first_by_subclass`,
     donuts, split violins, etc.).
2. **Move the helpers** (`smooth_vec`, `has_persistent_excursion`,
   `find_runs`, `assert_has_cols`) into `R/classify5HT_helpers.R`.
3. **Archive the dev scripts** — move `dev/classify_pct_change/`,
   `dev/build_response_classifier/`, `dev/gpt_build_2/`,
   `dev/updated_plotting_script/`, and `claude/checking_consistency.r`
   into `dev/.archive/classifier_iterations/` so the history is
   preserved but the working tree is clean.
4. **Write a single driver** in `dev/run_classifier.R` that sources
   nothing, just calls the promoted package functions:
   ```r
   library(proj5HT)
   comp <- readRDS("data-raw/compiled5HT.rds")
   res  <- classify_puff_response(comp$srtPuff, subclass_keep = c("L23_IT","L5_ET","L5_IT"))
   # res$df_prepped, res$per_cell, res$per_cell2
   plot_bucket_qc(res$df_prepped, cell_id = res$per_cell2$cell_name[1])
   ```

---

### 2. Hardcoded paths — use package-level defaults

26+ lines across `R/` use `~/proj5HT/den/serotonin/Data/Master_UI_data.ods`
as a function default. This works on your machine but breaks for anyone
else (or when the project moves).

**Fix:** Add a central path resolver, similar to how projHCT uses
`projHCT::params`:

```r
# Already built: R/params5HT.R + data/params5HT.rda
# Access via:
proj5HT::params5HT$paths$master_ui
proj5HT::params5HT$paths$selected_md
proj5HT::params5HT$paths$rookery
proj5HT::params5HT$sheets$master_macaque
proj5HT::params5HT$classifier$rapid_window
proj5HT::params5HT$exclude
```

Then replace every hardcoded default with `params5HT$master_ui`. Same
for `sheets5HT.R` which hard-references `~/proj5HT/den/…` and
`~/proj5HT/rookery`.

---

### 3. `figures5HT.R` uses `source()` — convert to package functions

[R/figures5HT.R](R/figures5HT.R) calls `source("~/proj5HT/R/figures_df_build.R")`
and `source("~/proj5HT/figs/scripts/srtHumanData.R")`, etc. instead of
calling package functions. This:
- bypasses the package namespace
- means `devtools::check()` can't catch errors
- breaks portability

**Fix:** The figure-building scripts in `figs/scripts/` should either be
promoted into `R/` as proper functions, or converted to parameterised
Rmd reports in `rmd/`.

---

### 4. Duplicate function definitions

The function `prep_puff_df` / `prepare_puff_df` is defined in **5
separate files**. `extract_bucket_features`, `add_smoothed_trace`,
`has_persistent_excursion` each appear in 3+ files. This makes it
impossible to know which version is active.

**Rule going forward:** If a function is called by more than one script,
it belongs in `R/` with an `@export` tag. Dev scripts should only
contain one-off experiments and should call package functions, never
redefine them.

---

### 5. `R/` files that are not exported and may be dead

These files in `R/` define functions that are **not in NAMESPACE** (not
exported):

| File | Function | Status |
|---|---|---|
| `baseline_stim_5HTpuff.R` | `baseline_stim_5HTpuff()` | Superseded by `spikePuff_from_asp()` |
| `combine_sweeps.R` | `combine_sweeps()` | Stub (body references undefined `sweeps`) |
| `dfs5HT.R` | `dfs5HT()` | Stub (calls undefined `spikePuff5HT`) |
| `figures5HT.R` | `figures5HT()` | Source-based dispatcher, not a proper function |
| `figures_df_build.R` | `figure_df_build()` | Works, but only invoked via `source()` |
| `plotBasicPuff.R` | `plotBasicPuff()` | Legacy, may still be used interactively |
| `plotWashIn.R` | `plotWashIn()` | Legacy |
| `5HT_viewer.R` | `srt_viewer()` | Legacy viewer, references undefined `tst` |

**Fix:** Audit each one. If still used, add `@export` and fix. If
superseded, move to `dev/.archive/`.

---

### 6. The `claude/` directory

`claude/checking_consistency.r` is a 1,688-line monolithic dump that
duplicates `dev/classify_pct_change/df_prep_and_per_cell_build.R`
nearly verbatim. It appears to be a review/consistency-check artifact.

**Fix:** Move to `dev/.archive/claude_review/`. The useful bits should
already be captured in the promoted classifier functions.

---

### 7. Missing `Imports` in DESCRIPTION

The package uses `readODS`, `readr`, `stringr`, `jsonlite`, `zoo`,
`scales`, `forcats`, `rlang`, `tibble`, `reticulate`, and `Seurat` at
runtime, but none appear in the `DESCRIPTION` `Imports:` field.
`devtools::check()` will flag these.

**Fix:** Add runtime dependencies to `Imports:` and
suggest-only dependencies (Seurat, reticulate) to `Suggests:`:

```
Imports:
    dplyr,
    tidyr,
    ggplot2,
    readODS,
    readr,
    stringr,
    zoo,
    scales,
    tibble,
    rlang,
    nphys,
    projHCT
Suggests:
    Seurat,
    reticulate,
    jsonlite,
    forcats
```

---

### 8. `.Rhistory` files still tracked

Despite `.Rhistory` being in `.gitignore`, five `.Rhistory` files exist
on disk:
- `proj5HT/.Rhistory`
- `dev/.Rhistory`
- `dev/JeremyScripts/.Rhistory`
- `dev/Seurat/.Rhistory`
- `rmd/.Rhistory`

**Fix:**
```bash
find . -name '.Rhistory' -delete
git rm --cached $(git ls-files '*.Rhistory') 2>/dev/null
```

---

### 9. `tmp/` contains classification artifacts

`tmp/` holds `cell_classification.odt`, `classification_MD.csv`, and
subdirectories `basicPuff/` and `ketWash/`. These appear to be
intermediate exports that should either be moved into `data-raw/` (if
they're curated inputs) or deleted (if they're reproducible outputs).

`tmp/` is already in `.gitignore`, so this is a local housekeeping
issue.

---

### 10. Suggested new systems

#### a. A `params5HT` object (like projHCT `params`)

Centralise project paths, default windows, threshold values, and
excluded cells. Load once at package attach time via `.onLoad()` or as
a lazy-loaded data object.

#### b. A `sheets5HT` lazy-data upgrade

Currently `sheets5HT()` reads ODS files from disk every time it's
called. Consider:
1. Call it once during package build → save as `data/sh5HT.rda`
   (you already do this, but the rebuild logic in `sheets5HT()` is
   fragile).
2. Add a `refresh_sheets5HT()` that re-reads and updates the `.rda`,
   removing the `askYesNo` / `system("R CMD INSTALL …")` hack
   currently embedded in `sheets5HT()`.

#### c. A `dev/run_all.R` master driver

Replace the ~300-line `dev/proj5HT.R` scratch file with a small set
of purpose-named scripts:
- `dev/01_build_nnests.R` — loop over cells, build/update child payloads
- `dev/02_build_asp_rmp.R` — ASP + RMP pipelines
- `dev/03_classify.R` — run classifier on compiled data
- `dev/04_figures.R` — generate manuscript figures

#### d. Vignette for the classifier

Move the classifier workflow into a proper vignette
(`vignettes/response_classifier.Rmd`) that reads `compiled5HT.rds`,
runs the pipeline, and renders QC + summary plots. This makes the
methodology reproducible and reviewable.

---

### Quick-win checklist

- [ ] Archive `claude/checking_consistency.r`
- [ ] Archive `dev/gpt_build_2/`, `dev/updated_plotting_script/`
- [ ] Delete `.Rhistory` files
- [ ] Promote `prep_puff_df`, `extract_bucket_features`,
      `add_biphasic_features` into `R/classify5HT.R`
- [ ] Add `Imports` to DESCRIPTION
- [ ] Replace `~/proj5HT/…` defaults with `params5HT$…`
- [ ] Remove stubs: `combine_sweeps.R`, `dfs5HT.R`
- [ ] Fix `srt_viewer()` (references undefined `tst`)
