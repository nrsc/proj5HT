<!-- README.md — proj5HT -->

# proj5HT — Serotonin Modulation of Cortical Neurons

proj5HT is a development bundle for analysing serotonin (5-HT) effects
on cortical neuron spiking in macaque and human tissue. It is designed
as a **child project** that piggybacks on the
[projHCT](../projHCT/README.md) infrastructure: every cell in the 5-HT
study already lives inside the HCT ecosystem, so proj5HT re-uses HCT
metadata, file indexing, NWB extraction, and the nnest framework rather
than duplicating them.

---

## How proj5HT piggybacks on projHCT

### The parent → child nnest model

projHCT builds a canonical **header nnest** (`-hct.rds`) for every
patched cell. That object is the single source of truth for cell
identity, working/rookery directories, metadata, NWB paths, and
versioning.

proj5HT adds a **child payload** (`-srt.rds`) that sits beside the
parent in the same rookery directory. The child carries serotonin-specific
analysis results (`dfs`, `ana`, `smry`, `figs`, …) while
inheriting—*not copying*—the header truth from the parent at load time.

```
rookery/<cell>/
  ├── <cell>-hct.rds          ← projHCT header (source of truth)
  ├── <cell>-srt.rds          ← proj5HT child payload
  ├── <cell>-5ht_md_snapshot.csv
  ├── <cell>-srtPuff.csv      ← per-cell spike-puff export
  ├── <cell>-srtDrugWash.csv  ← per-cell drug-wash export
  ├── <cell>-rmp.csv          ← per-cell RMP export
  └── …
```

The key functions that manage this relationship live in
`R/child_payload5HT.R`:

| Function | Role |
|---|---|
| `nnest5HT(cell)` | Create / rebuild the child payload. Loads the HCT header via `projHCT::loadHCT()`, extracts "header truth" fields (`cell`, `wd`, `rd`, `md`, `files`, …), and writes `-srt.rds`. |
| `load5HT(cell)` | Load the child payload and (optionally) **rehydrate** it with current HCT header truth, so downstream code always sees up-to-date metadata without a full rebuild. |
| `hct_header()` / `apply_hct_header()` | Extract and merge the authoritative header fields from an HCT nnest into a child object. |

A fingerprint (`parent_state`) is stored inside each child so
`load5HT()` can detect when the parent HCT has changed on disk and
re-inject updated header fields in-memory without rewriting the child.

### Direct projHCT calls used throughout proj5HT

Almost every analysis function reaches into projHCT at run-time:

| projHCT reference | Where it appears | What it provides |
|---|---|---|
| `projHCT::sheets$MD` | `compile5HT`, `data5HT`, `DrugWashIn`, `build_asp`, `build_rmp`, `baseline_stim_5HTpuff` | The master metadata table for all patched cells (species, subclass, NWB path, depth, date, etc.). |
| `projHCT::sheets$map` | `data5HT` | NHP annotation / cell-mapping tables (used for cross-referencing cristina's L2/3 dataset). |
| `projHCT::sheets$files$nnest_*` | `child_payload5HT` (`get_idx_map`) | Fast file-to-cell index maps so `load5HT` / `nnest5HT` can locate nnests without disk scanning. |
| `projHCT::sheets$cleanMD` | `talk.Rmd` | A cleaned metadata frame used for presentations and Seurat integration. |
| `projHCT::loadHCT()` | `nnest5HT`, `load5HT`, `5HT_viewer`, `analysisNaAvail`, `vignettes` | Load or build the parent HCT nnest for a given cell. |
| `projHCT::nnestHCT()` | `nnest5HT` (fallback) | If no parent HCT nnest exists for a cell, build it on the fly before creating the child. |
| `projHCT::params$rookery` | `child_payload5HT` | The canonical rookery root path used for disk-scan fallbacks. |

### nphys — the shared foundation

Both projHCT and proj5HT depend on the [nphys](https://github.com/nrsc/nphys)
package for low-level electrophysiology utilities:

- `nphys::extractNWB()` — read sweeps from NWB files
- `nphys::APanalysis()` — spike detection & instantaneous rate
- `nphys::convertIndex()` — sample ↔ time conversion
- `nphys::fileD()` — extract a cell name from a file path

---

## Project layout

```
proj5HT/
├── R/                        ← Exported & internal functions
│   ├── child_payload5HT.R    ← nnest5HT / load5HT (parent–child wiring)
│   ├── build_asp.R           ← Aligned Spike Protocol builder
│   ├── build_rmp.R           ← Resting Membrane Potential builder
│   ├── spikePuff_from_asp.R  ← Spike-puff derived from ASP
│   ├── DrugWashIn.R          ← Drug wash-in analysis (covers Ket + WAY)
│   ├── DrugWashIn_from_asp.R ← Drug wash-in derived from ASP
│   ├── compile5HT.R          ← Compile per-cell CSVs into project-wide frames
│   ├── data5HT.R             ← Assemble dataframes (UI sheets + compiled data)
│   ├── sheets5HT.R           ← Load/cache 5-HT-specific ODS sheets
│   ├── figures5HT.R          ← Figure dispatch / sourcing
│   ├── analysisNaAvail.R     ← Sodium-availability modelling (R + Python)
│   ├── plot_asp_concat.R     ← Concatenated protocol heatmap + trace
│   ├── plot_asp_single_protocol.R ← Single-protocol trace + heat strip
│   ├── plot_sidebyside_heat_wash.R ← Side-by-side wash-in comparison
│   ├── run_asp_pipeline.R    ← Full pipeline wrapper (dfs + figs + washin)
│   ├── run_asp_dfs_pipeline.R    ← DFS-only sub-pipeline
│   ├── run_asp_figs_pipeline.R   ← Figures-only sub-pipeline
│   ├── run_asp_ana_pipeline.R    ← Analysis tables sub-pipeline
│   └── …
├── dev/
│   └── proj5HT.R             ← Main interactive loop / scratch driver
├── data/
│   └── sh5HT.rda             ← Cached sheets5HT() snapshot
├── data-raw/
│   ├── compile_srtPuff.csv   ← Compiled spike-puff data
│   ├── compile_srtDrugWash.csv
│   ├── compile_srtKetWash.csv
│   ├── compiled5HT.rds
│   ├── HCN23_dataset.csv
│   ├── hct_cluster_matched.csv
│   └── Seurat/               ← Seurat objects for transcriptomic analysis
├── den/serotonin/             ← Experiment-tracker ODS sheets & metadata
├── rookery                    ← Symlink → shared rookery (same tree as projHCT)
├── rmd/                       ← R Markdown reports & templates
│   ├── srtAP_analysis.Rmd
│   ├── AHP_types5HT.Rmd
│   ├── Seurat-*.Rmd           ← Transcriptomic figures
│   └── templates/
├── vignettes/
│   └── python_filtering.Rmd
├── figs/                      ← Output figures
│   ├── basicPuff/
│   ├── exp_concat_heat/
│   ├── exp_single_heat/
│   ├── pch5HT/
│   ├── washIn/
│   ├── traces/
│   ├── qc/
│   └── scripts/               ← Figure-building R scripts
├── py/                        ← Python helper modules
├── doc/                       ← Presentation templates
├── DESCRIPTION
├── NAMESPACE
└── Makefile
```

---

## Core workflow

### 1. Build or load child payload

```r
library(projHCT)
library(proj5HT)

# Create a 5-HT child nnest (reads the parent HCT header)
srt <- nnest5HT("QF25.26.023.19.06.03")

# Or load an existing one (rehydrates header truth from HCT)
srt <- load5HT("QF25.26.023.19.06.03", tag = "-srt.rds")
```

### 2. Build the Aligned Spike Protocol (ASP)

`build_asp()` reads the Master UI sheet, identifies `spikeTTL` /
`TTLstack` protocols, extracts NWB sweeps, runs `nphys::APanalysis()`,
and builds time-aligned spike-rate traces with baseline-normalised
percent-change metrics.

```r
mMD <- readODS::read_ods("den/serotonin/Data/Master_UI_data.ods",
                          sheet = "Master_Macaque")
mMD <- mMD[!is.na(mMD$cell_name), ]

srt$dfs$asp <- build_asp(srt, mMD = mMD)
```

### 3. Derive spike-puff and drug wash-in analyses from ASP

The `_from_asp` variants reuse the already-built ASP container so sweeps
are not re-extracted:

```r
srt$dfs$blSpike <- spikePuff_from_asp(srt, mMD = mMD)
srt$dfs$dwin    <- DrugWashIn_from_asp(srt, mMD = mMD)
```

### 4. Build RMP (resting membrane potential)

```r
srt$dfs$rmp <- build_rmp(srt, mMD = mMD)
```

### 5. Generate figures

```r
# Per-protocol trace + heat strip
plot_asp_single_protocol(srt$dfs$asp, protocol_key = srt$dfs$asp$protocol_keys[1])

# Concatenated multi-protocol view
plot_asp_concatenated(srt$dfs$asp)

# Side-by-side wash-in comparison
plot_sidebyside_heat_wash(srt)
```

### 6. Run the full pipeline (batch)

The main interactive loop in `dev/proj5HT.R` iterates over all cells
in the Master UI sheet:

```r
for (cell in cells) {
  srt <- load5HT(cell, tag = "-srt.rds", rehydrate = FALSE)
  if (is.null(srt)) srt <- nnest5HT(cell)

  srt$dfs$asp    <- build_asp(srt, mMD = mMD)
  srt$dfs$rmp    <- build_rmp(srt, mMD = mMD, save_srt = FALSE)
  srt$dfs$blSpike <- spikePuff_from_asp(srt, mMD = mMD)
  srt$dfs$dwin   <- DrugWashIn_from_asp(srt, mMD = mMD)

  saveRDS(srt, file = srt$child_files$nnest)
}
```

### 7. Compile project-wide dataframes

```r
comp5HT <- compile5HT()
# → writes data-raw/compile_srtPuff.csv
# → writes data-raw/compile_srtDrugWash.csv
# → writes data-raw/compiled5HT.rds
```

### 8. Transcriptomic integration

R Markdown documents in `rmd/` (e.g. `Seurat-PatchSeq.Rmd`,
`Seurat-Figures-depthHTR.Rmd`) load Patch-seq Seurat objects from
`data-raw/Seurat/` and cross-reference HTR gene expression
(HTR1A–HTR2C) with electrophysiology response phenotypes.

---

## Dependency summary

```
┌──────────┐       ┌──────────┐       ┌──────────┐
│  nphys   │◄──────│ projHCT  │◄──────│ proj5HT  │
│ (extract,│       │ (nnest,  │       │ (child   │
│  APana,  │       │  sheets, │       │  payload,│
│  fileD)  │       │  loadHCT,│       │  ASP,    │
└──────────┘       │  params) │       │  RMP,    │
                   └──────────┘       │  compile)│
                                      └──────────┘
```

- **nphys** — portable electrophysiology primitives (NWB I/O, spike
  detection, index conversion).
- **projHCT** — project-level infrastructure (nnest creation, metadata
  sheets, file indexing, rookery management). proj5HT calls into projHCT
  at runtime via `projHCT::` namespace calls; it does **not** copy or
  re-export HCT functions.
- **proj5HT** — serotonin-specific analysis layer. Builds child nnests
  that inherit the HCT header, adds spike-puff / drug-wash / RMP
  analysis, and produces figures and compiled dataframes.

### Other R dependencies

`dplyr`, `tidyr`, `ggplot2`, `readODS`, `readr`, `stringr`, `jsonlite`,
`Seurat`, `reticulate` (for Python interop).

---

## Getting started

1. Install [nphys](https://github.com/nrsc/nphys) and
   [projHCT](../projHCT/README.md) first (see projHCT README for full
   instructions).
2. Clone this repository and open `proj5HT.Rproj` in RStudio.
3. Build & install the package:
   ```bash
   R CMD INSTALL --preclean --no-multiarch --with-keep.source .
   ```
4. Open `dev/proj5HT.R`, set your working directory to the project
   root, and step through the interactive workflow.

