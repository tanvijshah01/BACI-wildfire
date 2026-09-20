# Wildfire Biomass Recovery — Causal Inference Study (BACI)

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**BACI-wildfire** is an academic research project (UCSB MESM) estimating the causal effect of
wildfire severity on forest biomass recovery using a Before-After-Control-Intervention (BACI) /
staggered difference-in-differences design:

- MTBS fire perimeters (2000–2023, Western US)
- eMapR (1990–2023) and ctrees (2000–2025) biomass (annual, Landsat-based)
- Callaway-Sant'Anna (2021) staggered DiD estimator, which handles the heterogeneous treatment
  timing inherent in wildfire occurrence data

**Current phase:** Exploratory data analysis (EDA) and West-wide data pipeline build-out

**Goal:** Panel-based causal inference paper for an ecology/fire science journal

---

## Docs Map — where to write what

Keep each kind of information in exactly one place; link instead of copying.

| File | Owns |
|---|---|
| `README.md` | Public one-screen overview and quick start |
| `CLAUDE.md` (this file) | Stable orientation: design, stack, directory tree, run order, style, current status |
| `DATA_DOWNLOAD_GUIDE.md` | Process reference: environment setup (GRIT + laptop), download → crop → extract commands |
| `NOTES.md` | Technical gotchas, design decisions + open questions, dated log of findings, literature notes |
| `EDA_PLAN.md` | Live EDA plan: questions, status, and outputs per analysis document |

---

## Directory Structure

```
BACI/
├── scripts/
│   ├── python/                                  # ctrees/arraylake extraction
│   │   ├── 02_explore_ctrees_zarr.py            # inspect arraylake zarr store structure
│   │   ├── 03_download_ctrees_ca.py             # download ctrees AGB: CA raster, fire-polygon
│   │   │                                        #   extraction, native ~100m GeoTIFFs (validated
│   │   │                                        #   baseline — do not modify)
│   │   └── 04_download_ctrees_west.py           # same, generalized to the 11-state West bbox
│   │                                             #   (~4.1x CA's pixel area); outputs "_west_"
│   │                                             #   tier shared across states by 08; memory-
│   │                                             #   hardened for GRIT, runs A (raw TIFF) → B → C
│   ├── r/                                       # data processing & analysis
│   │   ├── 00_crop_emapr_to_ca.R                # [retired] crop CONUS eMapR TIFs to CA
│   │   ├── 00_crop_emapr_to_west.R              # crop/mask each raw eMapR year to the union of
│   │   │                                        #   all 11 Western states (skip-safe, validity-checked)
│   │   ├── 01_create_emapr_100m_tifs.R          # downsample CA eMapR 30m -> ~90m ("100m")
│   │   ├── 01_create_emapr_300m_tifs.R          # downsample CA eMapR 30m -> 300m (fast extraction)
│   │   ├── 01_precompute_emapr_ca.R             # summary stats + display TIFs for EDA years
│   │   ├── 02_compare_biomass_ranges.R          # eMapR vs ctrees range comparison (2000/2001/2003)
│   │   ├── 02_extract_emapr_within_fires.R      # [retired] eMapR AGB within fires, CA-only mask
│   │   ├── 03_prepare_forest_mask.R             # [retired] CA-only NLCD 2004 mask (30m/100m, 1/NA)
│   │   ├── 04_extract_ctrees_within_fires.R     # [retired] ctrees AGB within fires, CA-only mask
│   │   ├── 05_prepare_forest_masks_west.R       # CURRENT: per-state NLCD 2004 forest masks (0/1)
│   │   ├── 06_extract_pct_forest_within_fires.R # % forest cover within each MTBS fire perimeter
│   │   ├── 07_extract_emapr_within_fires_new.R  # CURRENT: eMapR AGB within fires, per-state mask
│   │   ├── 08_extract_ctrees_within_fires_new.R # CURRENT: ctrees AGB within fires, per-state mask
│   │   ├── check_raw_emapr_files.R              # standalone validity check of raw eMapR composites
│   │   ├── validate_west_pipeline.R             # validation harness: 05–08 vs. retired CA baseline
│   │   └── tmp_mtbs_count.R                     # scratch — not part of any pipeline
│   └── run_extraction.bat                       # STALE: launches the archived GEE NBR script
├── analysis/                                    # Quarto exploratory & analysis documents
│   ├── 01_mtbs_exploration.qmd
│   ├── 02_ctrees_biomass_exploration.qmd
│   ├── 03_emapr_biomass_exploration.qmd
│   ├── 04_data_summary.qmd
│   ├── biomass_within_fires.qmd                 # CURRENT: ctrees vs eMapR within fires (Western pipeline)
│   ├── biomass_within_fires_old.qmd             # [retired] CA-only forest-mask version
│   ├── mtbs_assessment_comparison.qmd           # MTBS Initial vs Extended assessment bias
│   ├── west_pipeline_sanity_check.qmd           # regression-validation report for scripts 05–08
│   └── run_03_emapr.R                           # scratch: standalone copy of 03's setup code (1990/92/93)
├── manuscript/                                  # final manuscript (manuscript.qmd)
├── archive/nbr_landsat_approach/                # archived GEE/NBR biomass approach (see NOTES.md)
├── data/
│   ├── raw/                                     # untouched downloaded data — do not modify/delete
│   │   ├── mtbs/                                # MTBS + Burn Severity Program shapefiles
│   │   └── emapr_biomass/                       # CONUS-wide eMapR composites (~27.7 GB/year)
│   └── processed/
│       ├── ctrees/                               # ctrees TIFs + biomass_fire_polygons_ctrees*.csv
│       ├── emapr_biomass_ca/                     # CA-clipped eMapR TIFs, stats RDS, forested CSVs
│       ├── emapr_biomass_west/                   # West-wide clipped eMapR TIFs (~1 GB/year)
│       ├── forest_mask/                          # NLCD 2004 forest masks (per-state, multi-resolution)
│       ├── control_pixels/                       # never-burned control pixel locations & time series
│       └── validation/                           # outputs of validate_west_pipeline.R
├── figures/                                      # output plots and maps
├── output/                                       # panel/export CSVs
├── README.md  CLAUDE.md  NOTES.md  EDA_PLAN.md  DATA_DOWNLOAD_GUIDE.md
```

---

## Technology Stack

- **R** — primary language: all data processing, panel construction, visualization, and causal
  inference
- **Python** — biomass data acquisition only: ctrees AGB download/extraction via the arraylake
  zarr store (`scripts/python/`)
- **Quarto (`.qmd`)** — exploratory analysis documents and the final manuscript

### Key R Packages
- `tidyverse` — data wrangling
- `sf` — spatial vector data (MTBS fires)
- `terra` — raster data
- `did` — Callaway-Sant'Anna implementation
- `fixest` — fast fixed effects (alternative specifications)
- `modelsummary` — regression tables
- `ggplot2` + `tmap` — visualization
- `FedData` — download NLCD land cover data
- `here`, `glue` — path/string handling in scripts

### Key Python Packages
- `arraylake`, `zarr`, `xarray`, `netCDF4` — ctrees zarr store access and NetCDF output
- `geopandas`, `shapely`, `rasterio` — vector/polygon operations and GeoTIFF I/O
- `pandas`, `numpy` — tabular/array handling

### Environments
- **GRIT (primary for data work):** code at `~/BACI-wildfire`; `data/` is a symlink to
  `~/BACI-review/data` (shared project storage); Python via a venv at `~/BACI-wildfire/.venv`.
  Run long jobs inside `tmux`. Setup: `DATA_DOWNLOAD_GUIDE.md` Part 1.
- **Laptop (Windows):** the original development machine; too small for the ~1 TB raw eMapR archive.

---

## Data

### Data Sources

| Dataset | Source | Contents |
|---|---|---|
| `data/raw/mtbs/mtbs_fod_pts_data/` | USGS MTBS Program | Fire occurrence point locations (1984–2025, ~30,390 fires) |
| `data/raw/mtbs/mtbs_perimeter_data/` | USGS MTBS Program | Fire perimeter polygons — key fields: Fire_ID, Year, BurnBndAc, Severity1 |
| `data/raw/mtbs/burn_severity_fod_pts_data/` | USGS Burn Severity Program | BSP fire occurrence points |
| `data/raw/mtbs/burn_severity_perimeter_data/` | USGS Burn Severity Program | BSP fire perimeter polygons |
| `data/raw/emapr_biomass/` | eMapR lab FTP (`islay.ceoas.oregonstate.edu`) | 30m annual AGB, CONUS composites, 1990–2023 |
| `data/processed/ctrees/` | ctrees (arraylake zarr) | Annual ML-based AGB, ~100m native, 2000–2025 |

All raw data comes from federal/public sources — do not modify or delete files in `data/raw/`.
Shapefiles are standard ESRI format (`.shp`, `.dbf`, `.shx`, `.prj`, `.cpg`) with FGDC metadata
(`.xml`). Files are large (100MB+ each) — avoid loading entire datasets into memory; use spatial
filters or chunked reads.

**Coverage:** All fires >1000 acres; 11 contiguous Western US states (AZ, CA, CO, ID, MT, NV, NM,
OR, UT, WA, WY), conifer forests, fires 2000–2023 (biomass series extend earlier/later — see above).

**Control sites:** Never-burned sites within the same ecoregion — exact strategy still open (see
`NOTES.md` → Design decisions and `data/processed/control_pixels/`).

### Data Download
Full instructions — GRIT setup, the `arraylake auth login` step, eMapR FTP pull, cropping, and the
ctrees West download — are in **[`DATA_DOWNLOAD_GUIDE.md`](DATA_DOWNLOAD_GUIDE.md)**. Use the `rclone`
flow for new eMapR years, not the old manual `curl.exe --ftp-pasv` loop. Key facts:

- **Raw eMapR composites are kept permanently** (~27.7 GB/year, ~1 TB total; PI request) and must be
  cropped to the study region before use — loading them in Quarto causes multi-minute stalls.
  `00_crop_emapr_to_west.R` crops/masks each locally-available year to the union of the 11 Western states
  (`data/processed/emapr_biomass_west/`, ~1 GB/year). `00_crop_emapr_to_ca.R` is the retired CA-only
  predecessor. Both are skip-safe. Validate raw files with `check_raw_emapr_files.R` before trusting them.
- **Forest mask (NLCD):** `FedData::get_nlcd()` (no login). The current pipeline uses **NLCD 2004**,
  per-state, 0/1-encoded fraction-forest masks built by `05_prepare_forest_masks_west.R` (classes 41
  Deciduous, 42 Evergreen, 43 Mixed) — replaces the retired CA-only 1/NA mask.
- **Never trust `file.exists()`** as proof a raster is valid — see `NOTES.md` → Technical gotchas.
- **`terra::extract()` inside Quarto on Windows:** Quarto buffers chunk output until the chunk finishes,
  which makes terra's C++ threading appear frozen. Run extraction scripts from the R console or via
  `Rscript`, never inside a Quarto chunk.

---

## Script Pipeline

The pipeline builds per-state NLCD forest masks, then extracts biomass (eMapR and ctrees) within MTBS
fire perimeters, restricted to forested pixels. Ctrees acquisition (Python) and the R processing/
extraction pipeline run independently; their outputs are combined in the `analysis/` documents.

**Current run order (R, from project root, outside Quarto):**

```r
# 1 — per-state NLCD 2004 forest masks (0/1), multiple resolutions
Rscript scripts/r/05_prepare_forest_masks_west.R

# 2 — % forest cover within each MTBS fire perimeter
#     -> feeds analysis/mtbs_assessment_comparison.qmd
Rscript scripts/r/06_extract_pct_forest_within_fires.R

# 3 — eMapR AGB within fire perimeters, forested pixels only
#     -> data/processed/emapr_biomass_west/biomass_fire_polygons_emapr_west_<years>_100m_forested.csv
Rscript scripts/r/07_extract_emapr_within_fires_new.R

# 4 — ctrees AGB within fire perimeters, forested pixels only
#     -> data/processed/ctrees/biomass_fire_polygons_ctrees_west_forested.csv
Rscript scripts/r/08_extract_ctrees_within_fires_new.R
```

Scripts 05–08 feed `analysis/biomass_within_fires.qmd` and `analysis/mtbs_assessment_comparison.qmd`.
They are multi-state (driven by a `STATES_TO_RUN` vector; `WESTERN_STATES` lists all 11), skip-safe
(resume per state/year), and write one combined CSV with a `STUSPS` column. Delete the relevant output
CSV to force re-extraction. **`STATES_TO_RUN` is currently `c("CA", "WY")` in `05` and `c("CA")` in
`06`–`08`** — widen it as West-wide eMapR/ctrees rasters become available for more states/years.
Design rationale and validation results: `NOTES.md` (2026-08-30 entry).

**Retired pipeline (`00`–`04`, CA-only):** an earlier, 1/NA-encoded forest-mask version
(`03_prepare_forest_mask.R` → `02_extract_emapr_within_fires.R` / `04_extract_ctrees_within_fires.R`),
kept only because it feeds `analysis/biomass_within_fires_old.qmd` as a validation baseline. Do not
extend it for new work.

**Ctrees acquisition (Python, before the R pipeline needs ctrees TIFs):**

```
python scripts/python/02_explore_ctrees_zarr.py   # inspect zarr store structure (run first)
python scripts/python/03_download_ctrees_ca.py    # CA baseline — validated, do not modify
python scripts/python/04_download_ctrees_west.py  # West-wide (~6,800 fires); multi-hour, run in tmux
```

**Ancillary R scripts** (`01_create_emapr_100m_tifs.R`, `01_create_emapr_300m_tifs.R`,
`01_precompute_emapr_ca.R`, `02_compare_biomass_ranges.R`) build downsampled/display rasters and summary
stats for the EDA documents and are not required by the fire-extraction pipeline.

**Planned:** build the unit×year panel from these extraction outputs and run the Callaway-Sant'Anna
estimator.

---

## Analysis Approach

### Study Design
- **Type:** Staggered difference-in-differences (staircase design)
- **Treatment:** Fire occurrence (continuous treatment variable: severity)
- **Panel structure:** Sites (i) × Years (t)
- **Method:** Callaway-Sant'Anna (2021) staggered DiD estimator
- **Software:** R `did` package

### Key Assumptions
- Parallel trends (testable with pre-fire data)
- No anticipation (fires are unanticipated shocks)
- SUTVA (no interference between units)

### Identification Strategy
- Unit fixed effects control for time-invariant site characteristics
- Time fixed effects control for common shocks
- Requires only parallel trends (weaker than conditional independence)

Design decisions (fire-size threshold, severity classes, ecoregions, controls, fire complexes, buffers)
and open questions are tracked in `NOTES.md` → Design decisions & open questions. Literature notes are in
`NOTES.md` → Literature notes.

### Technical Notes
- **MTBS severity:** categorical (Unburned, Low, Moderate, High, Increased Greenness) or continuous
  dNBR/RdNBR
- **Panel balance:** unbalanced panels are fine (Callaway-Sant'Anna handles this)

### Potential Problems
- Not enough never-burned controls in fire-prone regions
- Fires too clustered spatially (violates SUTVA)
- eMapR biomass may have gaps/clouds
- Recent fires (2020–2023) have short recovery time

---

## Current Status

*Last updated 2026-09-20 (evening). Reflects the GRIT archive. Detailed history is in `NOTES.md`.*

### Done
- [x] MTBS fire perimeters downloaded
- [x] `07`/`08` rewritten for multi-state West reading shared West-wide rasters; ctrees matches the
      retired CA baseline exactly on fire selection (r = 1.000 in the Aug validation harness; eMapR
      re-validation pending)
- [x] Validation harness (`validate_west_pipeline.R` + `west_pipeline_sanity_check.qmd`)
- [x] West-wide eMapR crop script (`00_crop_emapr_to_west.R`) and raw-file validity checker
      (`check_raw_emapr_files.R`)
- [x] `04_download_ctrees_west.py` fully re-run on GRIT: 26/26 raw TIFs, valid `.nc`, 177,242-record fire
      CSV — **but see the corrupt-years item below before trusting 2000/2001 specifically**
- [x] `06`/`07`/`08`'s MTBS-loading OOM fixed (OGR SQL read-time filter) and verified on GRIT (304 CA
      fires, exact match to the validated baseline)
- [x] GRIT's actual memory cap identified: 4 GiB, enforced via cgroup v2 (not a laptop-era red herring —
      see `NOTES.md` "Technical gotchas")

### Incomplete — next actions
- [ ] **`ctrees_2000_west_100m.tif` / `ctrees_2001_west_100m.tif` are confirmed 100% NaN** — Part A's
      validity check is now fixed (commit `4686a2c`) but not yet exercised against these two files; delete
      them plus the downstream `.nc`/CSV and re-run `04` (exact commands in `NOTES.md`'s 2026-09-20 entry)
- [ ] **`05`'s `FedData::get_nlcd()` call OOMs on GRIT** downloading NLCD 2004 for CA — root cause
      suspected (see `NOTES.md`) but not fixed; forest masks for CA/WY are NOT currently confirmed present
      on GRIT despite earlier laptop-era status saying so — verify directly before assuming built
- [ ] **Cross-validate ctrees West CSV's CA rows against the `03` baseline** — `03_download_ctrees_ca.py`
      itself is separately dying with no traceback on GRIT (plausibly the same cgroup cap, not confirmed)
- [ ] **Raw eMapR archive:** 19/34 years confirmed complete on GRIT (2026-09-08); fetch the rest
      (`DATA_DOWNLOAD_GUIDE.md` §2.2), validate, then run `00_crop_emapr_to_west.R`
- [ ] **Widen `STATES_TO_RUN`** in `06`–`08` (and `05`) beyond CA/WY once West rasters exist for the needed
      years — blocked on the two items above first
- [ ] **Root-cause the residual eMapR bias vs. ctrees** and re-check the border-fire count gap on real
      multi-state output (`biomass_within_fires.qmd` §7; `NOTES.md`)
- [ ] **Resolve MTBS Initial vs. Extended assessment bias** (`mtbs_assessment_comparison.qmd`)
- [ ] **Expand EDA to all available years:** update `STUDY_YEARS` (currently 2005–2010 in `07`/`08`) and
      the year-range params in the analysis `.qmd` files once West crops are complete
- [ ] Clean up stale files (see "Stray files" below)

### Planned
- [ ] Build the unit × year panel (combine MTBS + eMapR + ctrees extraction outputs)
- [ ] Finalize control-site strategy
- [ ] Run Callaway-Sant'Anna analysis; event study plots; robustness checks
- [ ] Write manuscript

### Stray files (candidates for removal — confirm before deleting)
`scripts/run_extraction.bat` (launches an archived script), `scripts/r/tmp_mtbs_count.R` (scratch),
`analysis/run_03_emapr.R`, `analysis/04_data_summary.rmarkdown`, `analysis/01_exploration_files/`,
`analysis/mtbs_assessment_writeup_cache/`.

---

## File Naming Conventions

### Scripts
- Number prefix: `01_`, `02_`, `03_` (execution order within a pipeline; see Script Pipeline above)
- Descriptive name: `extract_biomass_gee`, `process_mtbs`
- Language suffix: `.py` for Python, `.R` for R, `.qmd` for Quarto

### Data Files
- `raw/` - untouched downloaded data
- `processed/` - cleaned, filtered data
- `final/` - analysis-ready datasets

### Figures
- Descriptive names: `fire_severity_map.png`, `event_study.png`
- High resolution: 300 dpi for publication

---

## Coding Style and Organization

- **Commented outline at the top of every script**: Each script should open with a block comment listing the major sections/steps in order (e.g., `# 1. Load data`, `# 2. Filter by severity`, `# 3. Export`). This acts as a table of contents so the logic is legible without reading every line.
- **Section headers throughout**: Divide scripts into clearly labeled sections that match the outline above.
- **Inline comments for non-obvious logic**: Explain *why*, not just *what* — especially for spatial operations, parameter choices, and DiD assumptions.
- **One concern per function**: Helper functions should do one thing and be named to reflect it.
- **All hardcoded values centralized**: Paths, CRS, date ranges, and filter thresholds belong at the top of each script (or a shared config), not scattered through the body.

---

## Notes for AI Assistant

When helping with this project:

1. **R is primary language** - use R for all analysis, Python only for ctrees data acquisition
2. **Quarto for reports** - use `.qmd` for exploratory analysis, `.R` for production scripts
3. **Spatial data:** Use `sf` package for vectors, `terra` for rasters
4. **Citations:** This is for academic publication, provide proper citations
5. **Causality:** Be precise about causal language vs. correlational
6. **Fire ecology:** Assume user knows fire ecology, focus on methods
7. **Current phase:** EDA plus West-wide pipeline build-out — understand the data before building the panel
8. **Pipeline currency:** Prefer scripts `05`–`08` over the retired `00`–`04` pipeline for any new
   extraction work
9. **Docs hygiene:** put new information in the file the Docs Map assigns it to; update "Current Status"
   here (not scattered copies) when pipeline state changes; record decisions/findings in `NOTES.md`

### Common Tasks
- Mapping fire perimeters
- Summary statistics by severity class
- Temporal/spatial distributions
- Identifying suitable control sites
- Checking data quality

### Avoid
- Suggesting cross-sectional methods (we're doing panel)
- Using TWFE without noting bias issues
- Mixing causal and correlational language
- Overcomplicated code (keep it readable)
- Whole-raster terra masking on the big ctrees/eMapR rasters (see `NOTES.md` → Technical gotchas)
