# EDA Plan — Wildfire Biomass Recovery

**Phase:** Pre-analysis exploration
**Goal:** Understand the MTBS and biomass data well enough to make key design decisions before building the analysis panel.

---

## Section 1 — Setup

- Install and load R packages: `tidyverse`, `sf`, `tigris`, `ggspatial`, `patchwork`, `scales`, `here`, `glue`, `png`
- Define study parameters in one place: `WESTERN_STATES`, `YEAR_MIN`, `YEAR_MAX`, `MIN_FIRE_ACRES`, `MTBS_PATH`
- Define `ASMT_COLORS` palette (colorblind-safe, one color per assessment type)
- Disable s2 with `sf_use_s2(FALSE)` — MTBS polygons have self-touching edges that s2 rejects

---

## Section 2 — MTBS Data Collection Methodology

Prose section (no code) documenting the MTBS methodology for reference:

- NBR and dNBR formulas (with LaTeX rendering)
- Per-fire threshold calibration using unburned reference areas
- Five burn severity classes (including Increased Greenness)
- Assessment types and image timing (Initial, Extended, Emergency, SS variants)
- Fire size thresholds (≥1,000 ac western US) and temporal coverage (1984–present)

**Citation:** Eidenshink et al. (2007); Dataset DOI: https://doi.org/10.5066/P9IED7RZ

---

## Section 3 — Load, Inspect, and Filter MTBS Data

**Script:** `analysis/01_mtbs_exploration.qmd`

- Load `data/raw/mtbs/mtbs_perimeter_data/mtbs_perims_DD.shp` with `sf::st_read()`
- Confirm key fields exist: `event_id`, `ig_date`, `incid_type`, `burnbndac`, `low_t`, `mod_t`, `high_t`, `dnbr_offst`, `asmnt_type`
- Parse `year` from `ig_date`; flag sentinel no-data threshold values (−9999, 9999)
- Check geometry validity; fix with `st_make_valid()` if needed

**Filter steps (print row count after each):**
- `incid_type == "Wildfire"` (excludes prescribed fires)
- `year >= 2000 & year <= 2023`
- Spatial join with `tigris::states()` for 11 Western states

**Severity sanity checks:**
- Threshold orientation (burned vs. greenness vs. mixed non-monotone)
- Scale check: values must NOT all fall in 0–100 (would indicate percentages, not dNBR)
- Variance check: non-zero variance across fires

---

### 3.1 Fire Locations Map

**What:** MTBS wildfire perimeters across the 11 Western states, colored by assessment type

**How:**
- Project to Albers Equal Area (EPSG:5070)
- Color by `asmnt_type` using `ASMT_COLORS` palette; `alpha = 0.65` for overlapping fires
- Scale bar (`ggspatial::annotation_scale()`) and north arrow

**Output:** `figures/fire_locations_map.png` ✓ Done

---

### 3.2 Fires Per Year — All Wildfires

**What:** Annual wildfire count (2000–2023) with a linear trend line

**How:**
- `count(year)` + `complete()` to fill zero-fire years
- `geom_col()` + `geom_smooth(method = "lm")` overlay

**Output:** `figures/fires_per_year.png` ✓ Done

---

### 3.3 Assessment Type Summary

**What:** Bar chart of fire count by MTBS assessment type category

**How:**
- `count(asmnt_type)` after `trimws()`; order by count descending
- `geom_col()` with count labels above each bar; `ASMT_COLORS` palette
- Caption explains SS limitation and why Extended is preferred

**Output:** `figures/asmnt_type_summary.png` ✓ Done

---

### 3.4 Extended Assessment Fires Only

**What:** Map restricted to `asmnt_type == "Extended"` with total count and % of all fires in subtitle

**How:**
- Filter `mtbs_west_alb` to `asmnt_type == "Extended"` (excludes Extended (SS))
- Also define `mtbs_extended_west` (non-projected) for use in 3.5 and 3.6
- Single fill color (#d6604d); subtitle auto-populates n and %

**Output:** `figures/fire_locations_extended_only.png` ✓ Done

---

### 3.5 Fires Per Year — All Fires and Extended Only

**What:** Side-by-side bar charts comparing annual wildfire counts for all assessments vs. Extended only

**How:**
- `count(year)` for `mtbs_extended_west`; `complete()` for zero-fire years
- Shared y-axis ceiling (`max(both series) * 1.15`) for fair visual comparison
- Combine with `patchwork`: `p35_all | p35_ext`

**Output:** `figures/fires_per_year_comparison.png` ✓ Done

---

### 3.6 Fire Size Distribution

**What:** Side-by-side log₁₀ histograms of fire size — all assessments vs. Extended only

**How:**
- `geom_histogram(bins = 60)` + `scale_x_log10()` for both panels
- Median annotated with dashed vertical line and text label on each panel
- Combine with `patchwork`

**Output:** `figures/fire_size_comparison.png` ✓ Done

---

### 3.7 Burn Severity — dNBR and RdNBR

**What:** Prose explanation of dNBR and RdNBR, then:
- (a) Map of all BSP CBI field plots colored by dNBR
- (b) dNBR histogram — all BSP CBI field plots, Western US
- (c) RdNBR histogram — all BSP CBI field plots, Western US

**Data source:** `data/raw/mtbs/burn_severity_perimeter_data/bsp_perims_DD.shp`, filtered to `map_prog == "CBI"`, excluding sentinel 9999 and `|preNBR| ≤ 10`

**Note:** No extended-only comparison — BSP CBI field-assessed fires and MTBS satellite-assessed fires are largely different programs covering different fires, so joining on `event_id` returns near-zero matches. The Extended/Initial distinction (satellite image timing) does not apply to CBI field data.

**RdNBR formula:** `rdnbr = dnbr_val / sqrt(abs(prenbr_val / 1000))` — ref: Miller & Thode (2007)

**Outputs:**
- `figures/severity_dnbr_map.png` ✓ Done
- `figures/severity_dnbr_histogram.png` ✓ Done
- `figures/severity_rdnbr_histogram.png` ✓ Done

---

## Section 4 — GEE Biomass Extraction

**Script:** `scripts/python/01_extract_biomass_gee.py`

**Overview:** Use the LandTrendR algorithm in Google Earth Engine to extract annual fitted biomass (or a biomass-correlated vegetation index) time series spatially constrained to MTBS fire perimeters. Reference: https://emapr.github.io/LT-GEE/

### 4a — Setup and Authentication
- Authenticate GEE: `ee.Authenticate()` then `ee.Initialize()`
- Import the LandTrendR module: `require('users/emaprlab/public:Modules/LandTrendr.js')` (JavaScript) or use the Python `earthengine-api` equivalent
- Verify access to the eMapR biomass asset (check exact asset path in GEE catalog: likely `projects/eMapR/biomass` or similar)

### 4b — Define Area of Interest from MTBS Perimeters
- Load the filtered MTBS fire perimeters (from Section 3) as a `ee.FeatureCollection`
- Use `.geometry()` or `.union()` to create an AOI geometry covering all fire perimeters in the study region
- This AOI is passed as the `aoi` parameter to all LandTrendR calls — it spatially constrains processing to fire-affected areas only

### 4c — Run LandTrendR and Extract Fitted Time Series
- Build the annual surface reflectance collection: `buildSRcollection(startYear, endYear, startDay, endDay, aoi)`
- Transform to the target vegetation index (NBR recommended for fire/biomass): `transformSRcollection()`
- Run LandTrendR segmentation: `runLT(startYear, endYear, startDay, endDay, aoi, index, ftvList, runParams)`
- Extract annual fitted values: `getFittedData(lt, startYear, endYear, index)`
  - Returns an image stack with one band per year — convert to tabular format with `arrayFlatten()` before export
- Key LandTrendR parameters (set in script header):
  - `maxSegments` — max number of temporal breakpoints
  - `spikeThreshold = 0.9`
  - `recoveryThreshold = 0.25`
  - `pvalThreshold = 0.1`
  - `minObservationsNeeded = 6`

### 4d — Export Results
- Export annual fitted values per fire perimeter as a CSV to `output/biomass_timeseries.csv`
- Columns should include: `event_id`, `fire_year`, `year`, `fitted_nbr`

---

### 4.1 Annual Mean NBR Over Time

**What:** Line plot of mean fitted NBR across all fire sites, by year (1984–2023)

**How:**
- Aggregate `output/biomass_timeseries.csv` by year (mean ± SE)
- Overlay a linear trend line (`geom_smooth(method = "lm")`)

**Output:** `figures/biomass_trend_over_time.png` — pending GEE extraction

---

### 4.2 Biomass Before and After Fire Events

**What:** Event-study plot showing mean NBR indexed to the fire year (year 0)

**How:**
- Compute `rel_year = year - fire_year`; subset to −5 to +10 window
- Average across fires by relative year (mean ± SE)
- Vertical dashed line at year 0

**Output:** `figures/biomass_before_after_fire.png` — pending GEE extraction

---

## Testing, Sanity Checks, and Visual QA

Run these checks inline within each script/Quarto chunk immediately after the operation they validates — do not batch them at the end.

---

### R — MTBS Data Checks

**Row counts at each filter step**
```r
cat("Total fires loaded:           ", nrow(mtbs_raw), "\n")
cat("After Wildfire filter:        ", nrow(mtbs_wf), "\n")
cat("After year filter (2000–23):  ", nrow(mtbs_wf), "\n")
cat("After Western US spatial join:", nrow(mtbs_west), "\n")
```

**Field and CRS checks**
```r
stopifnot(all(c("event_id", "ig_date", "incid_type", "burnbndac",
                "low_t", "mod_t", "high_t") %in% names(mtbs_west)))
stopifnot(!any(duplicated(mtbs_west$event_id)))
stopifnot(all(mtbs_west$year >= 2000 & mtbs_west$year <= 2023))
```

**Geometry validity**
```r
invalid <- sum(!st_is_valid(mtbs_raw))
if (invalid > 0) {
  warning(paste(invalid, "invalid geometries — run st_make_valid()"))
  mtbs_raw <- st_make_valid(mtbs_raw)
}
```

**State coverage check**
```r
target_states <- c("AZ","CA","CO","ID","MT","NV","NM","OR","UT","WA","WY")
missing <- setdiff(target_states, unique(mtbs_west$STUSPS))
if (length(missing) > 0) warning(paste("Missing states:", paste(missing, collapse=", ")))
```

---

### Python — GEE Biomass Export Checks (`01_extract_biomass_gee.py`)

**Post-export CSV validation**
```python
import pandas as pd

df = pd.read_csv("output/biomass_timeseries.csv")

print(f"Rows: {len(df)}, Columns: {list(df.columns)}")
assert {"event_id", "fire_year", "year", "fitted_nbr"}.issubset(df.columns), "Missing expected columns"
assert df["year"].between(1984, 2023).all(), "Years outside expected range"

na_by_fire = df.groupby("event_id")["fitted_nbr"].apply(lambda x: x.isna().all())
if na_by_fire.any():
    print(f"WARNING: {na_by_fire.sum()} fires have all-NA biomass values")

print(df["fitted_nbr"].describe())
assert df["fitted_nbr"].dropna().between(-1, 1).all(), "Fitted values outside expected NBR range"
```

---

### Plot QA Checklist

After generating each plot, verify the following before saving:

#### Legend (all plots)
- [ ] Legend title is descriptive (not a raw column name)
- [ ] Font size ≥ 10pt for all legend text
- [ ] Colors are colorblind-safe (`viridis`, `RColorBrewer "Set2"`, or `MetBrewer`)
- [ ] Legend does not overlap data

#### 3.1 — Fire Locations Map
- [ ] Scale bar and north arrow present
- [ ] All 11 Western states visible without excess whitespace
- [ ] State boundaries visible for geographic reference

#### 3.2 — Fires Per Year
- [ ] X-axis spans 2000–2023 with labeled ticks every 2–4 years
- [ ] Y-axis label reads "Number of wildfires"
- [ ] Trend line and bars are visually distinct

#### 3.3 — Assessment Type Summary
- [ ] Count labels visible above each bar
- [ ] Bars ordered by count descending

#### 3.4 — Extended Fires Map
- [ ] Scale bar and north arrow present
- [ ] Subtitle shows exact n and % of total

#### 3.5 — Fires Per Year Comparison
- [ ] Both panels share the same y-axis scale
- [ ] Panel titles clearly identify All vs. Extended

#### 3.6 — Fire Size Distribution
- [ ] X-axis is log₁₀ scale with readable tick labels
- [ ] Median annotated on each panel; both panels comparable

#### 3.7 — Burn Severity
- [ ] Map: color clipped to 2nd–98th percentile; scale bar present
- [ ] dNBR histogram: zero line (dotted) and median (dashed) visible
- [ ] RdNBR histogram: zero line (dotted) and median (dashed) visible

#### 4.1 — Biomass Trend Over Time
- [ ] X-axis spans 1984–2023
- [ ] Y-axis label includes index name (e.g., "Mean NBR (fitted)")
- [ ] No flat line at 0 (would indicate missing/failed data)

#### 4.2 — Biomass Before/After Fire
- [ ] X-axis labeled "Years relative to fire"
- [ ] Year 0 marked with a vertical dashed line
- [ ] Confidence intervals visible but not wider than the signal

**Resolution check (all plots)**
```r
img <- png::readPNG("figures/fire_locations_map.png")
cat("Dimensions (px):", dim(img)[2], "x", dim(img)[1], "\n")
# Expected: width/300 x height/300 inches at 300 dpi
```

---

## EDA Outputs

### Completed (MTBS only)
- `figures/fire_locations_map.png` — 3.1 (all fires, colored by assessment type)
- `figures/fires_per_year.png` — 3.2 (all wildfires, temporal trend)
- `figures/asmnt_type_summary.png` — 3.3 (bar chart: count by assessment type)
- `figures/fire_locations_extended_only.png` — 3.4 (Extended-only map with count)
- `figures/fires_per_year_comparison.png` — 3.5 (all vs. Extended comparison)
- `figures/fire_size_comparison.png` — 3.6 (size distribution, all vs. Extended)
- `figures/severity_dnbr_map.png` — 3.7a (CBI field plot locations)
- `figures/severity_dnbr_histogram.png` — 3.7b (dNBR distribution, all CBI plots)
- `figures/severity_rdnbr_histogram.png` — 3.7c (RdNBR distribution, all CBI plots)
- `analysis/01_mtbs_exploration.html` — rendered report

### Pending (requires GEE biomass extraction)
- `figures/biomass_trend_over_time.png` — 4.1
- `figures/biomass_before_after_fire.png` — 4.2

---

## References

### Primary MTBS Data Reference
Eidenshink, J., Schwind, B., Brewer, K., Zhu, Z. L., Quayle, B., & Howard, S. (2007). A project for monitoring trends in burn severity. *Fire Ecology*, 3(1), 3–21. https://www.mtbs.gov/sites/mtbs/files/inline-files/Eidenshink-final.pdf

**Dataset DOI:** https://doi.org/10.5066/P9IED7RZ

This paper describes the MTBS methodology: use of Landsat-derived NBR and dNBR to map burn severity, per-fire threshold calibration using unburned reference areas, the five severity classes (including Increased Greenness), fire size thresholds (≥1,000 ac western US), and temporal coverage since 1984.

### Burn Severity Index Reference
Miller, J. D., & Thode, A. E. (2007). Quantifying burn severity in a heterogeneous landscape with a relative version of the delta Normalized Burn Ratio (dNBR). *Remote Sensing of Environment*, 109(1), 66–80.

This paper introduces RdNBR (= dNBR / √|preNBR / 1000|), which normalizes dNBR for pre-fire vegetation density. RdNBR is preferred over raw dNBR for cross-fire comparison and as the continuous treatment variable in the DiD model.
