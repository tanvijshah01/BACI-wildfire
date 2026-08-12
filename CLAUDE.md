# Wildfire Biomass Recovery — Causal Inference Study (BACI)

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**BACI-wildfire** is an academic research project (UCSB MESM) estimating the causal effect of
wildfire severity on forest biomass recovery using a Before-After-Control-Intervention (BACI) /
staggered difference-in-differences design:

- MTBS fire perimeters (2000–2023, Western US)
- eMapR and ctrees biomass (annual, Landsat-based)
- Callaway-Sant'Anna (2021) staggered DiD estimator, which handles the heterogeneous treatment
  timing inherent in wildfire occurrence data

**Current phase:** Exploratory data analysis (EDA)

**Goal:** Panel-based causal inference paper for an ecology/fire science journal

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
│   │                                             #   tier, meant to be shared across states by 07/08
│   └── r/                                       # data processing & analysis
│       ├── 00_crop_emapr_to_ca.R                # [retired] crop CONUS eMapR TIFs to CA
│       ├── 00_crop_emapr_to_west.R              # crop/mask each raw eMapR year to the union of
│       │                                        #   all 11 Western states (skip-safe)
│       ├── 01_create_emapr_100m_tifs.R          # downsample CA eMapR 30m -> ~90m ("100m")
│       ├── 01_create_emapr_300m_tifs.R          # downsample CA eMapR 30m -> 300m (fast extraction)
│       ├── 01_precompute_emapr_ca.R             # summary stats + display TIFs for EDA years
│       ├── 02_compare_biomass_ranges.R          # eMapR vs ctrees range comparison (2000/2001/2003)
│       ├── 02_extract_emapr_within_fires.R      # [retired] eMapR AGB within fires, CA-only mask
│       ├── 03_prepare_forest_mask.R             # [retired] CA-only NLCD 2004 mask (30m/100m, 1/NA)
│       ├── 04_extract_ctrees_within_fires.R     # [retired] ctrees AGB within fires, CA-only mask
│       ├── 05_prepare_forest_masks_west.R       # CURRENT: per-state NLCD 2004 forest masks (0/1);
│       │                                        #   STATES_TO_RUN currently "CA" only, though
│       │                                        #   WESTERN_STATES already lists all 11
│       ├── 06_extract_pct_forest_within_fires.R # % forest cover within each MTBS fire perimeter
│       │                                        #   (STATES_TO_RUN currently "CA" only)
│       ├── 07_extract_emapr_within_fires_new.R  # CURRENT: eMapR AGB within fires, per-state mask
│       │                                        #   (STATE_FIPS currently hardcoded "CA")
│       └── 08_extract_ctrees_within_fires_new.R # CURRENT: ctrees AGB within fires, per-state mask
│                                                 #   (STATE_FIPS currently hardcoded "CA")
├── analysis/                                    # Quarto exploratory & analysis documents
│   ├── 01_mtbs_exploration.qmd
│   ├── 02_ctrees_biomass_exploration.qmd
│   ├── 03_emapr_biomass_exploration.qmd
│   ├── 04_data_summary.qmd
│   ├── 05_biomass_masking_slides.qmd            # slide deck on the forest-masking approach
│   ├── biomass_within_fires.qmd                 # CURRENT: ctrees vs eMapR within fires (Western pipeline)
│   ├── biomass_within_fires_old.qmd             # [retired] CA-only forest-mask version
│   └── mtbs_assessment_comparison.qmd           # MTBS Initial vs Extended assessment bias
├── paper/                                       # final manuscript (manuscript.qmd)
├── data/
│   ├── raw/                                     # untouched downloaded data — do not modify/delete
│   │   ├── mtbs/                                # MTBS + Burn Severity Program shapefiles
│   │   └── emapr_biomass/                       # CONUS-wide eMapR composites (~27.7 GB/year)
│   └── processed/
│       ├── ctrees/                               # ctrees TIFs + biomass_fire_polygons_ctrees*.csv
│       ├── emapr_biomass_ca/                     # CA-clipped eMapR TIFs, stats RDS, forested CSVs
│       ├── emapr_biomass_west/                   # West-wide clipped eMapR TIFs (~1 GB/year)
│       ├── forest_mask/                          # NLCD 2004 forest masks (per-state, multi-resolution)
│       └── control_pixels/                       # never-burned control pixel locations & time series
├── figures/                                      # output plots and maps
├── output/                                       # panel/export CSVs
└── DATA_DOWNLOAD_GUIDE.md                        # eMapR FTP/Nextcloud + ctrees arraylake download instructions
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
- `arraylake`, `zarr`, `xarray` — ctrees zarr store access
- `geopandas`, `shapely` — vector/polygon operations
- `pandas`, `numpy` — tabular/array handling

---

## Data

### Data Sources

| Dataset | Source | Contents |
|---|---|---|
| `data/raw/mtbs/mtbs_fod_pts_data/` | USGS MTBS Program | Fire occurrence point locations (1984–2025, ~30,390 fires) |
| `data/raw/mtbs/mtbs_perimeter_data/` | USGS MTBS Program | Fire perimeter polygons — key fields: Fire_ID, Year, BurnBndAc, Severity1 |
| `data/raw/mtbs/burn_severity_fod_pts_data/` | USGS Burn Severity Program | BSP fire occurrence points |
| `data/raw/mtbs/burn_severity_perimeter_data/` | USGS Burn Severity Program | BSP fire perimeter polygons |
| `data/raw/emapr_biomass/` | Google Earth Engine (eMapR/LandTrendR) | 30m annual AGB, CONUS composites |
| `data/processed/ctrees/` | ctrees (arraylake zarr) | Annual ML-based AGB, ~100m native |

All raw data comes from federal/public sources — do not modify or delete files in `data/raw/`.
Shapefiles are standard ESRI format (`.shp`, `.dbf`, `.shx`, `.prj`, `.cpg`) with FGDC metadata
(`.xml`). Files are large (100MB+ each) — avoid loading entire datasets into memory; use spatial
filters or chunked reads.

**Coverage:** All fires >1000 acres; 11 contiguous Western US states (AZ, CA, CO, ID, MT, NV, NM,
OR, UT, WA, WY), conifer forests, 2000–2023.

**Control sites:** Never-burned sites within the same ecoregion — exact strategy still being
determined in EDA (see `data/processed/control_pixels/`).

### Data Download

- **eMapR (FTP/Nextcloud) and ctrees (arraylake zarr)** — full download/access instructions,
  including the required `arraylake auth login` step and the Nextcloud (`rclone`) archive flow
  for raw eMapR composites, are in **[`DATA_DOWNLOAD_GUIDE.md`](DATA_DOWNLOAD_GUIDE.md)**. Use the
  Nextcloud/rclone flow for new years — do not go back to the old manual `curl.exe --ftp-pasv`
  loop.
- **Forest mask (NLCD)** — downloaded via `FedData::get_nlcd()` (no login required). The current
  pipeline uses **NLCD 2004**, per-state, 0/1-encoded fraction-forest masks built by
  `scripts/r/05_prepare_forest_masks_west.R` (forest classes 41 Deciduous, 42 Evergreen,
  43 Mixed). This replaces the retired CA-only, 1/NA-encoded mask from
  `scripts/r/03_prepare_forest_mask.R`.
- **CONUS eMapR composites** (~27.7 GB/year) must be cropped to the study region before use —
  loading them directly in Quarto causes multi-minute stalls. `scripts/r/00_crop_emapr_to_west.R`
  crops/masks each locally-available raw year to the union of all 11 Western states
  (`data/processed/emapr_biomass_west/`, ~1 GB/year); `00_crop_emapr_to_ca.R` is the retired
  CA-only predecessor. Both are skip-safe and only process missing years.

**Note on `terra::extract()` inside Quarto on Windows:** Quarto buffers chunk output until the
chunk finishes, which makes terra's C++ threading appear frozen. Always run extraction scripts
from the R console or via `Rscript`, never inside a Quarto chunk.

---

## Script Pipeline

The current pipeline builds per-state NLCD forest masks, then extracts biomass (eMapR and
ctrees) within MTBS fire perimeters, restricted to forested pixels. Ctrees data acquisition
(Python) and the R processing/extraction pipeline run independently; their outputs are combined
in the `analysis/` documents.

**Current run order (R, from project root, outside Quarto):**

```r
# 1 — per-state NLCD 2004 forest masks (0/1), multiple resolutions
#     STATES_TO_RUN currently c("CA") — edit to add more states as West-wide eMapR/ctrees
#     rasters become available for them
Rscript scripts/r/05_prepare_forest_masks_west.R

# 2 — % forest cover within each MTBS fire perimeter
#     -> feeds analysis/mtbs_assessment_comparison.qmd
Rscript scripts/r/06_extract_pct_forest_within_fires.R

# 3 — eMapR AGB within fire perimeters, forested pixels only
#     -> data/processed/emapr_biomass_ca/biomass_fire_polygons_emapr_..._newpipeline.csv
Rscript scripts/r/07_extract_emapr_within_fires_new.R

# 4 — ctrees AGB within fire perimeters, forested pixels only
#     -> data/processed/ctrees/biomass_fire_polygons_ctrees_forested_newpipeline.csv
Rscript scripts/r/08_extract_ctrees_within_fires_new.R
```

Scripts 05–08 feed `analysis/biomass_within_fires.qmd` and
`analysis/mtbs_assessment_comparison.qmd`, are skip-safe (append one state/year at a time), and
skip already-completed work. Delete the relevant output CSV to force a full re-extraction.

**Known gap:** scripts 07 and 08 still hardcode `STATE_FIPS <- "CA"`, and the `STATES_TO_RUN` in
scripts 05 and 06 is still `c("CA")` even though `WESTERN_STATES` already lists all 11 — these
need to be generalized once West-wide eMapR/ctrees rasters exist for the other states (see
`00_crop_emapr_to_west.R` above and Current Status below).

**Retired pipeline (`00`–`04`, CA-only):** an earlier, 1/NA-encoded forest-mask version of the
same idea (`03_prepare_forest_mask.R` → `02_extract_emapr_within_fires.R` /
`04_extract_ctrees_within_fires.R`), kept only because it feeds
`analysis/biomass_within_fires_old.qmd` as a validation baseline — the `05`–`08` pipeline was
cross-checked against it (ctrees matched exactly, r = 1.000; eMapR had a residual bias not yet
root-caused, see `biomass_within_fires.qmd` §7). Do not extend the retired pipeline for new work.

**Ctrees acquisition (Python, one-time, before the R pipeline needs ctrees TIFs):**

```
python scripts/python/02_explore_ctrees_zarr.py   # inspect zarr store structure (run first)
python scripts/python/03_download_ctrees_ca.py    # download CA raster + fire-polygon extraction
                                                    # + native ~100m annual GeoTIFFs
python scripts/python/04_download_ctrees_west.py  # same, generalized to the 11-state West bbox
                                                    # (~4.1x CA's pixel area, ~6,800 fires);
                                                    # multi-hour run — see DATA_DOWNLOAD_GUIDE.md
```

**Ancillary R scripts** (`01_create_emapr_100m_tifs.R`, `01_create_emapr_300m_tifs.R`,
`01_precompute_emapr_ca.R`, `02_compare_biomass_ranges.R`) build downsampled/display rasters and
summary stats for the EDA documents (`03_emapr_biomass_exploration.qmd`,
`04_data_summary.qmd`) and are not required by the fire-extraction pipeline above.

**Planned:** build the unit×year panel from these extraction outputs and run the Callaway-Sant'Anna
estimator (see Current Status below).

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

---

## Current Status

**Scope note:** EDA is expanding from a downscaled California-only pilot to the full 11-state
Western pipeline (scripts `05`–`08`). See `analysis/biomass_within_fires.qmd` and
`analysis/mtbs_assessment_comparison.qmd` for the current state of that expansion.

### Completed
- [x] Download MTBS fire perimeters
- [x] Build per-state NLCD 2004 forest masks for CA (`05_prepare_forest_masks_west.R`)
- [x] Extract eMapR + ctrees AGB within CA fire perimeters, forested pixels, current pipeline
      (`07`, `08`)
- [x] Cross-validate current pipeline against retired CA-only pipeline (ctrees r = 1.000; eMapR
      residual bias open, see `biomass_within_fires.qmd` §7)
- [x] West-wide crop script (`00_crop_emapr_to_west.R`) crops/masks each locally-available raw
      eMapR year to the union of all 11 Western states (skip-safe; run again as new raw years land)

### Incomplete — resume before expanding scope
- [ ] **Raw eMapR archive → Nextcloud:** years 1997–1999 and 2013–2023 are still not downloaded
      due to repeated local disk-space failures. Fix in progress: archive raw CONUS composites on
      Nextcloud instead of the laptop, fetched via `rclone` (installed at
      `C:\Users\shaht\bin\rclone.exe`, anonymous eMapR FTP remote configured, Nextcloud WebDAV
      remote still needs credentials). See "Nextcloud Archive (rclone)" in
      `DATA_DOWNLOAD_GUIDE.md` — do not use the old manual `curl.exe --ftp-pasv` loop for new years.
- [ ] **Generalize extraction scripts to all 11 states:** flip `STATE_FIPS`/`STATES_TO_RUN` in
      scripts `05`–`08` from `"CA"` to the full `WESTERN_STATES` list once West-wide eMapR/ctrees
      rasters exist for the needed years. `scripts/python/04_download_ctrees_west.py` (West-wide
      ctrees, ~6,800 fires) is written but not yet run — a multi-hour job, see
      `DATA_DOWNLOAD_GUIDE.md`.
- [ ] **Root-cause the residual eMapR bias vs. ctrees**
- [ ] **Resolve MTBS Initial vs. Extended assessment bias** (`mtbs_assessment_comparison.qmd`)
- [ ] **Expand EDA to all available years:** once West crops are complete for all years, update
      the relevant `STUDY_YEARS`/year-range params in the analysis `.qmd` files

### Planned
- [ ] Build the unit × year panel (combine MTBS + eMapR + ctrees extraction outputs)
- [ ] Finalize control-site strategy
- [ ] Run Callaway-Sant'Anna analysis
- [ ] Event study plots
- [ ] Robustness checks
- [ ] Write manuscript

---

## Exploratory Analysis Questions

### Data Quality
1. How many fires in MTBS 2000-2023, Western US?
2. What's the distribution of fire sizes?
3. What severity classes exist? How is it distributed?
4. Are there temporal trends in fire occurrence?
5. Spatial clustering of fires?

### Panel Structure
6. How many fires per year (treatment timing)?
7. What years have sufficient pre-fire data (for parallel trends testing)?

---

## Known Issues & Decisions

### Design Decisions (to be made in EDA)
- [ ] Which ecoregions to include?
- [ ] Minimum fire size threshold?
- [ ] Include moderate severity or only high severity?
- [ ] How to handle fire complexes (multiple fires same area)?
- [ ] Spatial buffer between sites?
- [ ] Time window (as long as possible - what the data allows)?

### Technical Notes
- **MTBS severity:** Categorical (Unburned, Low, Moderate, High, Increased Greenness)
- **Severity measurement:** Could use continuous dNBR instead of categories
- **Panel balance:** Don't need balanced panel (Callaway-Sant'Anna handles this)

### Potential Problems
- Not enough never-burned controls in fire-prone regions
- Fires too clustered spatially (violates SUTVA)
- eMapR biomass may have gaps/clouds
- Recent fires (2020-2023) have short recovery time

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

## Literature Context

### Existing Work
- **Descriptive:** Garcia 2017 (Rim Fire carbon), Reisch 2024, Stenzel 2019
- **Predictive:** Bright et al. 2019 (random forest, R²>0.7)
- **Causal attempts:** Ilangakoon et al. 2026 (GAM with space-for-time, lacks formal identification)

### Our Contribution
- First application of modern staggered DiD and continuous DiD to wildfire-biomass
- Relaxes untestable conditional independence assumption
- Leverages natural staggered timing of fires

### Methods Papers
- Callaway & Sant'Anna (2021) - Staggered DiD
- Goodman-Bacon (2021) - TWFE decomposition
- Sun & Abraham (2021) - Interaction-weighted estimator
- Liermann & Roni (2021) - Staircase design power analysis

### Domain Papers
- Ilangakoon et al. (2026) - Wildfire GAM study
- Bright et al. (2019) - Predictive modeling
- Garcia et al. (2017) - Rim Fire carbon

---

## Contact & Collaboration
- **PI:** [Name]
- **Collaborators:** [Names]
- **Code repository:** [GitHub URL if applicable]

---

## Notes for AI Assistant

When helping with this project:

1. **R is primary language** - use R for all analysis, Python only for ctrees/GEE data acquisition
2. **Quarto for reports** - use `.qmd` for exploratory analysis, `.R` for production scripts
3. **Spatial data:** Use `sf` package for vectors, `terra` for rasters
4. **Citations:** This is for academic publication, provide proper citations
5. **Causality:** Be precise about causal language vs. correlational
6. **Fire ecology:** Assume user knows fire ecology, focus on methods
7. **Current phase:** Focus on EDA - understanding data before extraction
8. **Pipeline currency:** Prefer scripts `05`–`08` over the retired `00`–`04` pipeline for any new
   extraction work (see Script Pipeline above)

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
