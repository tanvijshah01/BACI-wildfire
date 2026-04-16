# EDA Plan — Wildfire Biomass Recovery

**Phase:** Pre-analysis exploration
**Goal:** Understand the MTBS and biomass data well enough to make key design decisions before building the analysis panel.

---

## Step 1 — Load and Inspect MTBS Data

**Script:** `scripts/r/02_process_mtbs.R` / `analysis/01_exploration.qmd`

- Load `data/raw/mtbs/mtbs_perimeter_data/mtbs_perims_DD.shp` with `sf::st_read()`
- Print field names, classes, and a sample of rows (`glimpse()`, `head()`)
- Check CRS (NAD83, EPSG:4269) and geometry type
- Confirm key fields exist: `event_id`, `ig_date`, `incid_type`, `burnbndac`, `low_t`, `mod_t`, `high_t`
- Check for NAs, duplicates, malformed geometries (`st_is_valid()`)
- Note: disable s2 with `sf_use_s2(FALSE)` — MTBS polygons have self-touching edges that s2 rejects

**Output:** Console summary; note any data quality issues in `NOTES.md`

---

## Step 2 — Filter to Study Region and Period + MTBS Visualizations

- Filter `incid_type == "Wildfire"` (excludes prescribed fires)
- Subset to the 11 contiguous Western US states: Arizona, California, Colorado, Idaho, Montana, Nevada, New Mexico, Oregon, Utah, Washington, Wyoming — via spatial join with `tigris::states()`
- Subset to years 2000–2023 (parse year from `ig_date` string)
- Apply minimum fire size threshold (decision: see `NOTES.md`)
- Note how many fires remain after each filter

### Plot 1 — Fire Locations Map

**What:** MTBS wildfire perimeters across the 11 Western states, colored by dominant severity class

**How:**
- Project to Albers Equal Area (EPSG:5070) for display
- Color by dominant severity class (highest of `low_t`, `mod_t`, `high_t`)
- Add scale bar (`ggspatial::annotation_scale()`) and north arrow

**Script/doc:** `analysis/01_exploration.qmd`
**Output:** `figures/fire_locations_map.png` ✓ Done

---

### Plot 2 — Fires Per Year (Temporal Trend)

**What:** Annual wildfire count (2000–2023) with a linear trend line

**How:**
- Count fires by year (`count(year)`)
- Bar chart + `geom_smooth(method = "lm")` overlay
- Fill in any years with zero fires using `complete()`

**Script/doc:** `analysis/01_exploration.qmd`
**Output:** `figures/fires_per_year.png` ✓ Done

---

### Plot 3 — Fire Size Distribution

**What:** Histogram of fire sizes (`burnbndac`) on a log₁₀ scale

**How:**
- `geom_histogram()` with `scale_x_log10()`
- Annotate with median fire size

**Script/doc:** `analysis/01_exploration.qmd`
**Output:** `figures/fire_size_distribution.png` ✓ Done

---

### Plot 4 — Burn Severity Distribution

**What:** (4a) Bar chart of fire count by dominant severity class; (4b) histogram of high-severity area % per fire

**How:**
- Assign dominant severity: whichever of `low_t`, `mod_t`, `high_t` is highest
- 4a: `geom_col()` by dominant class with counts labeled
- 4b: `geom_histogram()` of `high_t` values across all fires

**Script/doc:** `analysis/01_exploration.qmd`
**Output:** `figures/severity_by_class.png` ✓ Done / `figures/severity_high_pct_distribution.png` ✓ Done

---

## Step 3 — Load and Inspect Biomass Time Series from eMapR LandTrendR (GEE)

**Script:** `scripts/python/01_extract_biomass_gee.py`

**Overview:** Use the LandTrendR algorithm in Google Earth Engine to extract annual fitted biomass (or a biomass-correlated vegetation index) time series spatially constrained to MTBS fire perimeters. Reference: https://emapr.github.io/LT-GEE/

### 3a — Setup and Authentication
- Authenticate GEE: `ee.Authenticate()` then `ee.Initialize()`
- Import the LandTrendR module: `require('users/emaprlab/public:Modules/LandTrendr.js')` (JavaScript) or use the Python `earthengine-api` equivalent
- Verify access to the eMapR biomass asset (check exact asset path in GEE catalog: likely `projects/eMapR/biomass` or similar)

### 3b — Define Area of Interest from MTBS Perimeters
- Load the filtered MTBS fire perimeters (from Step 2) as a `ee.FeatureCollection`
- Use `.geometry()` or `.union()` to create an AOI geometry covering all fire perimeters in the study region
- This AOI is passed as the `aoi` parameter to all LandTrendR calls — it spatially constrains processing to fire-affected areas only

### 3c — Run LandTrendR and Extract Fitted Time Series
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

### 3d — Export Results
- Export annual fitted values per fire perimeter as a CSV to `output/biomass_timeseries.csv`
- Columns should include: `event_id`, `year`, `fitted_index_value` (and optionally raw observed value)

---

### Plot 5 — Biomass Change Over Time

**What:** Line plot of mean fitted vegetation index (proxy for biomass) across all fire sites, by year (2000–2023)

**How:**
- Aggregate `output/biomass_timeseries.csv` by year (mean across all fire perimeters)
- Overlay a smoothed trend line
- Annotate years with anomalously high fire activity (cross-reference Plot 2)

**Script/doc:** `analysis/01_exploration.qmd`
**Output:** `figures/biomass_trend_over_time.png` — pending GEE extraction

---

### Plot 6 — Biomass Before and After Fire Events

**What:** Event-study style plot showing mean biomass indexed to the fire year (year 0), for N years before and after each fire event

**How:**
- For each fire in MTBS, compute relative year to fire: `rel_year = year - fire_year`
- Subset to a window (e.g., −5 to +10 years relative to fire)
- Average across fires by relative year, optionally stratified by dominant severity class
- Plot mean (± SE) biomass by relative year with a vertical line at year 0

**Script/doc:** `analysis/01_exploration.qmd`
**Output:** `figures/biomass_before_after_fire.png` — pending GEE extraction

---

## Testing, Sanity Checks, and Visual QA

Run these checks inline within each script/Quarto chunk immediately after the operation they validate — do not batch them at the end.

---

### R — MTBS Data Checks (`02_process_mtbs.R`)

**Row counts at each filter step**
```r
# Print before/after counts for every filter so you can catch silent drops
cat("Total fires loaded:          ", nrow(mtbs_raw), "\n")
cat("After Wildfire filter:       ", nrow(mtbs_wf), "\n")
cat("After year filter (2000-23): ", nrow(mtbs_years), "\n")
cat("After Western US spatial join:", nrow(mtbs_west), "\n")
```

**Field and CRS checks**
```r
stopifnot(all(c("event_id", "ig_date", "incid_type", "burnbndac",
                "low_t", "mod_t", "high_t") %in% names(mtbs_filtered)))
stopifnot(!any(duplicated(mtbs_filtered$event_id)))
stopifnot(all(mtbs_filtered$year >= 2000 & mtbs_filtered$year <= 2023))
```

**Geometry validity**
```r
invalid <- sum(!st_is_valid(mtbs_filtered))
if (invalid > 0) {
  warning(paste(invalid, "invalid geometries — run st_make_valid()"))
  mtbs_filtered <- st_make_valid(mtbs_filtered)
}
```

**State coverage check**
```r
target_states <- c("AZ","CA","CO","ID","MT","NV","NM","OR","UT","WA","WY")
missing <- setdiff(target_states, unique(mtbs_filtered$STUSPS))
if (length(missing) > 0) warning(paste("Missing states:", paste(missing, collapse=", ")))
```

---

### Python — GEE Biomass Export Checks (`01_extract_biomass_gee.py`)

**Post-export CSV validation**
```python
import pandas as pd

df = pd.read_csv("output/biomass_timeseries.csv")

print(f"Rows: {len(df)}, Columns: {list(df.columns)}")
assert {"event_id", "year", "fitted_index_value"}.issubset(df.columns), "Missing expected columns"
assert df["year"].between(2000, 2023).all(), "Years outside 2000-2023"

na_by_fire = df.groupby("event_id")["fitted_index_value"].apply(lambda x: x.isna().all())
if na_by_fire.any():
    print(f"WARNING: {na_by_fire.sum()} fires have all-NA biomass values")

print(df["fitted_index_value"].describe())
assert df["fitted_index_value"].dropna().between(-1, 1).all(), "Fitted values outside expected NBR range"
```

---

### Plot QA Checklist

After generating each plot, verify the following before saving:

#### Legend (all plots)
- [ ] Legend title is descriptive (not a raw column name)
- [ ] Font size ≥ 10pt for all legend text
- [ ] Colors are colorblind-safe (`viridis`, `RColorBrewer "Set2"`, or `MetBrewer`)
- [ ] Legend does not overlap data

#### Plot 1 — Fire Locations Map
- [ ] Scale bar present (`ggspatial::annotation_scale(location = "bl")`)
- [ ] North arrow present (`ggspatial::annotation_north_arrow()`)
- [ ] All 11 Western states visible without excess whitespace
- [ ] State boundaries visible for geographic reference

#### Plot 2 — Fires Per Year
- [ ] X-axis spans 2000–2023 with labeled ticks every 2–5 years
- [ ] Y-axis label reads "Number of wildfires"
- [ ] Trend line and bars are visually distinct

#### Plot 3 — Fire Size Distribution
- [ ] X-axis is log₁₀ scale with readable tick labels
- [ ] Median annotated with a vertical line and text label

#### Plot 4 — Severity Distribution
- [ ] (4a) Count labels visible above each bar
- [ ] (4b) X-axis labeled as percentage (0–100%)

#### Plot 5 — Biomass Trend Over Time
- [ ] X-axis spans 2000–2023
- [ ] Y-axis label includes index name (e.g., "Mean NBR (fitted)")
- [ ] No flat line at 0 (would indicate missing/failed data)

#### Plot 6 — Biomass Before/After Fire
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
- `figures/fire_locations_map.png` — Plot 1
- `figures/fires_per_year.png` — Plot 2
- `figures/fire_size_distribution.png` — Plot 3
- `figures/severity_by_class.png` — Plot 4a
- `figures/severity_high_pct_distribution.png` — Plot 4b
- `analysis/01_exploration.html` — rendered report

### Pending (requires GEE biomass extraction)
- `figures/biomass_trend_over_time.png` — Plot 5
- `figures/biomass_before_after_fire.png` — Plot 6
