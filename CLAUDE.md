# Wildfire Biomass Recovery - Causal Inference Study

## Project Overview
Estimating causal effects of wildfire severity on forest biomass recovery using:
- MTBS fire perimeters (2000-2023, Western US)
- eMapR biomass data from Google Earth Engine
- Callaway-Sant'Anna staggered difference-in-differences

**Current phase:** Exploratory data analysis (EDA)

**Goal:** Panel-based causal inference paper for ecology/fire science journal

---

## Data Sources

### MTBS Fire Perimeters
- **Source:** https://www.mtbs.gov/direct-download
- **File:** `data/raw/mtbs_perims_DD.shp`
- **Key fields:** Fire_ID, Year, BurnBndAc, Severity1
- **Coverage:** All fires >1000 acres, 1984-2023
- **Focus:** 11 contiguous Western US states (AZ, CA, CO, ID, MT, NV, NM, OR, UT, WA, WY), conifer forests, 2000-2023

### eMapR Biomass
- **Source:** Google Earth Engine (requires Python extraction)
- **Asset:** `projects/eMapR/biomass` (check actual name in GEE catalog)
- **Resolution:** 30m Landsat-based
- **Temporal:** Annual, 2000-2023
- **Extraction status:** NOT YET DONE (planned for `scripts/python/01_extract_biomass_gee.py`)

### Control Sites
- **Approach:** Never-burned sites within same ecoregion
- **Status:** TO BE DETERMINED in EDA

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

## Literature Context

### Existing Work
- **Descriptive:** Garcia 2017 (Rim Fire carbon), Reisch 2024, Stenzel 2019
- **Predictive:** Bright et al. 2019 (random forest, R²>0.7)
- **Causal attempts:** Ilangakoon et al. 2026 (GAM with space-for-time, lacks formal identification)

### Our Contribution
- First application of modern staggered DiD and continuous DiD to wildfire-biomass
- Relaxes untestable conditional independence assumption
- Leverages natural staggered timing of fires

---

## Technical Stack

### Languages
- **Python:** Google Earth Engine data extraction only
- **R:** All analysis, visualization, manuscript

### Key R Packages
- `tidyverse` - data wrangling
- `sf` - spatial vector data (MTBS fires)
- `terra` - raster data
- `did` - Callaway-Sant'Anna implementation
- `fixest` - fast fixed effects (alternative specifications)
- `modelsummary` - regression tables
- `ggplot2` + `tmap` - visualization

### Python Packages (for GEE only)
- `earthengine-api`
- `geopandas`
- `pandas`

---

## Current Status

### Completed
- [ ] Download MTBS fire perimeters
- [ ] Set up project structure
- [ ] Install R packages

### In Progress
- [ ] **Exploratory Data Analysis (CURRENT PHASE)**
  - Understand MTBS data structure
  - Map fire locations
  - Examine severity distributions
  - Identify suitable study region
  - Determine control site strategy

### Planned
- [ ] Extract eMapR biomass from GEE (Python script)
- [ ] Create panel dataset (combine MTBS + biomass)
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
- Number prefix: `01_`, `02_`, `03_` (execution order)
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

## References

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

1. **R is primary language** - use R for all analysis, Python only for GEE
2. **Quarto for reports** - use `.qmd` for exploratory analysis, `.R` for production scripts
3. **Spatial data:** Use `sf` package for vectors, `terra` for rasters
4. **Citations:** This is for academic publication, provide proper citations
5. **Causality:** Be precise about causal language vs. correlational
6. **Fire ecology:** Assume user knows fire ecology, focus on methods
7. **Current phase:** Focus on EDA - understanding data before extraction

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


# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**BACI-wildfire** is an academic research project (UCSB MESM) analyzing the causal impact of wildfires on forest biomass using a Before-After-Control-Intervention (BACI) / Difference-in-Differences design. The primary causal inference method is the **Callaway-Sant'Anna staggered DiD estimator**, which handles the heterogeneous treatment timing inherent in wildfire occurrence data.

## Technology Stack

- **Python** — biomass data extraction via Google Earth Engine (GEE) API
- **R** — data processing, panel construction, and causal inference (primary analysis language)
- **Quarto (.qmd)** — reproducible analysis documents and final manuscript

## Script Pipelines

The two pipelines run independently and their outputs are merged in `03_create_panel.R`.

### R Pipeline — MTBS Investigation (`scripts/r/`)

The R pipeline investigates the MTBS dataset located in `data/mtbs/` to characterize fire events, assign treatment status, and build the analysis panel.

| Script | Purpose |
|---|---|
| `config.R` | Project-wide parameters: file paths, CRS, date ranges, filter thresholds |
| `04_helper_functions.R` | Shared utility functions — sourced at the top of every other R script |
| `02_process_mtbs.R` | Load MTBS shapefiles, filter by region/severity/size, standardize fields, export clean fire event table |
| `03_create_panel.R` | Join MTBS fire events with GEE biomass outputs; assign treatment/control; build unit-year panel for DiD |

**Execution order:** `config.R` → `04_helper_functions.R` → `02_process_mtbs.R` → `03_create_panel.R`

### Python Pipeline — Biomass Extraction via LandTrendR (`scripts/python/`)

The Python pipeline uses the **LandTrendR (LT-GEE)** spectral-temporal segmentation algorithm inside Google Earth Engine to extract annual biomass trajectories for fire-affected and control pixels. Reference: https://emapr.github.io/LT-GEE/

| Script | Purpose |
|---|---|
| `01_extract_biomass_gee.py` | Authenticate GEE, define AOI from MTBS perimeters, build annual SR image collection, run LandTrendR, export fitted biomass time series to `output/` |

**LandTrendR key inputs and parameters:**
- Annual surface reflectance image collection (one composite per year, cloud-masked, within a target season)
- Core parameters to configure in `config.R` / script header:
  - `maxSegments` — maximum number of fitted temporal segments
  - `spikeThreshold` (default 0.9) — dampens single-year spectral spikes
  - `recoveryThreshold` (default 0.25) — prevents implausibly fast post-fire recovery
  - `pvalThreshold` (default 0.1) — rejects poorly fitting segment models
  - `minObservationsNeeded` (default 6) — minimum annual observations required to fit

**LandTrendR key GEE API calls:** `buildSRcollection()` → `buildLTcollection()` → `runLT()` → `getSegmentData()`

Analysis documents in `analysis/` follow the same numeric ordering and depend on the processed outputs of both pipelines.

## Data

All raw data lives in `data/` and comes from federal sources — do not modify or delete these files.

| Dataset | Source | Contents |
|---|---|---|
| `data/mtbs_fod_pts_data/` | USGS MTBS Program | Fire occurrence point locations (1984–2025, ~30,390 fires) |
| `data/mtbs_perimeter_data/` | USGS MTBS Program | Fire perimeter polygons |
| `data/burn_severity_fod_pts_data/` | USGS Burn Severity Program | BSP fire occurrence points |
| `data/burn_severity_perimeter_data/` | USGS Burn Severity Program | BSP fire perimeter polygons |

All datasets are standard ESRI shapefiles (`.shp`, `.dbf`, `.shx`, `.prj`, `.cpg`) with FGDC metadata (`.xml`). Files are large (100MB+ each) — avoid loading entire datasets into memory; use spatial filters or chunked reads.

## Directory Structure

```
BACI/
├── scripts/
│   ├── python/        # GEE extraction scripts
│   └── r/             # Data processing & analysis scripts
├── analysis/          # Quarto exploratory & analysis documents
├── paper/             # Final manuscript (manuscript.qmd)
├── data/              # Raw shapefiles (MTBS + Burn Severity)
├── figures/           # Output plots and maps
└── output/            # Processed data outputs (e.g., panel CSVs)
```

## Coding Style and Organization

- **Commented outline at the top of every script**: Each script should open with a block comment listing the major sections/steps in order (e.g., `# 1. Load data`, `# 2. Filter by severity`, `# 3. Export`). This acts as a table of contents so the logic is legible without reading every line.
- **Section headers throughout**: Divide scripts into clearly labeled sections that match the outline above.
- **Inline comments for non-obvious logic**: Explain *why*, not just *what* — especially for spatial operations, parameter choices, and DiD assumptions.
- **One concern per function**: Helper functions in `04_helper_functions.R` should do one thing and be named to reflect it.
- **All hardcoded values in `config.R`**: Paths, CRS, date ranges, LandTrendR parameters, and filter thresholds belong in `config.R`, not scattered through scripts.

## Testing and Sanity Checks

Sanity checks and plot QA are not optional — embed them directly in each script/Quarto chunk immediately after the operation they validate. See `EDA_PLAN.md` for the full checklist. Summary of requirements:

### Data Checks (R and Python)
- Print row counts before and after every filter step to catch silent data loss
- Use `stopifnot()` (R) or `assert` (Python) to enforce expected field names, CRS, year range, and no duplicate IDs
- Check geometry validity with `st_is_valid()` before any spatial join; call `st_make_valid()` if needed
- After GEE export: confirm CSV has expected columns, year range 2000–2023, no all-NA fire records, and values within a plausible index range (e.g., NBR −1 to 1)
- Confirm all 11 target states (AZ, CA, CO, ID, MT, NV, NM, OR, UT, WA, WY) are present after spatial filtering

### Quick Visual Confirmation (run before saving final figures)
- After filtering MTBS: plot raw fire locations over state boundaries (`plot(st_geometry(...))`) to confirm spatial extent is correct
- After GEE export: plot a single fire's biomass time series to confirm the trajectory is non-flat and shows expected post-fire dip

### Plot and Map QA
Every saved figure must pass these checks before the script is considered complete:

**Legend**
- Title is human-readable (not a raw column name)
- Font size ≥ 10pt
- Colors are colorblind-safe — use `viridis`, `RColorBrewer "Set2"`, or `MetBrewer`
- Legend does not overlap data; reposition with `theme(legend.position = ...)` if needed

**Map scale** (`fire_locations_map.png` and any spatial figure)
- Scale bar present: `ggspatial::annotation_scale(location = "bl")` or `tmap::tm_scale_bar()`
- North arrow present: `ggspatial::annotation_north_arrow()` or `tmap::tm_compass()`
- Zoom covers all 11 Western states without excess whitespace; set bounding box explicitly if `ggplot2` auto-zoom crops states

**Plot 1 (biomass trend)**
- X-axis spans 2000–2023 with labeled ticks; Y-axis label includes index name and units
- Trend line and annual means are visually distinct; no flat line at zero (signals missing data)

**Plot 2 (before/after fire)**
- X-axis labeled "Years relative to fire"; year 0 marked with vertical dashed line
- Confidence intervals visible but not wider than the signal
- Severity strata (if shown) clearly labeled in legend

**Resolution check (R)**
```r
# Confirm saved PNG is publication-quality (300 dpi, ~10x7 in → ~3000x2100 px)
img <- png::readPNG("figures/biomass_trend_over_time.png")
cat("Image dimensions (px):", dim(img)[2], "x", dim(img)[1], "\n")
```

## Key Concepts

- **Treatment**: A pixel/unit being burned in a given fire event (identified by MTBS perimeter overlap)
- **Control**: Unburned pixels in the same region and time period
- **Running the Callaway-Sant'Anna estimator**: Use the `did` R package (`att_gt()` and `aggte()` functions)
- **CRS**: All spatial data should be projected to a consistent CRS before analysis — set this in `config.R`
