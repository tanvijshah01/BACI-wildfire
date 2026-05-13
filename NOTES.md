# Project Notes — Wildfire Biomass Recovery

Working notes on decisions, findings, and open questions. Add entries in reverse chronological order (newest at top).

---

## 2026-05-13 — Fire polygon extraction: shapely thinning replaced by rasterio rasterization

Part B of `scripts/python/03_download_ctrees_ca.py` extracts mean Ctrees AGB within each MTBS CA fire polygon for every year. Two approaches were tried:

**Approach 1 — Shapely point-in-polygon with grid thinning (abandoned):**
For each fire × year combination, the script built a meshgrid of pixel centres inside the polygon's bounding box and called `shapely.within()` on each point. Large CA fires (e.g. the 2020 August Complex at ~1M acres) have bounding boxes containing millions of pixels. Even after thinning the grid to a maximum of 80,000 test points (`step = ceil(sqrt(n_bbox / 80_000))`), the 1,064-fire × 26-year loop ran at ~1 year per 5 minutes — a projected ~130-minute total runtime. The thinning also introduced a sampling approximation that under-sampled polygon boundaries.

**Approach 2 — Rasterio scanline rasterization + precomputed masks (current):**
`rasterio.features.rasterize()` burns each polygon onto a boolean grid aligned to the CA pixel coordinates using a C-level scanline algorithm — O(pixels in bounding box), no sampling, handles both Polygon and MultiPolygon. Crucially, masks are precomputed **once** for all 1,064 fires before the year loop. The 26-year extraction then uses only `raw[np.ix_(yi, xi)][mask]` — pure numpy array indexing with zero shapely/rasterio overhead per year. Falls back to `matplotlib.path.Path.contains_points()` (C extension) if rasterio is unavailable.

The rasterio approach is typically 10–50x faster for large polygons and eliminates the thinning approximation.

---

## 2026-05-12 — NBR/Landsat GEE extraction approach archived

The initial biomass extraction pipeline used raw annual Landsat composites from Google Earth Engine to compute **NBR (Normalized Burn Ratio)** as a biomass proxy. This approach was archived because **NBR is a unitless spectral index (−1 to +1), not a biomass quantity**. Using it as an outcome variable means DiD treatment effects are expressed in NBR units, which lack ecological interpretability and will not satisfy reviewers at ecology/fire science journals who expect biomass (Mg/ha) or carbon (MgC/ha).

Converting NBR to biomass requires an empirical calibration model — typically a Random Forest regression trained on co-located FIA field plot measurements, with climate and topographic variables added to handle spectral saturation above ~300 Mg/ha. This calibration pipeline is exactly what the eMapR lab built and published:

> Kennedy, R.E., Yang, Z., Gorelick, N., Braaten, J., Cavalcante, L., Cohen, W.B., & Healey, S. (2018). Implementation of the LandTrendr algorithm on Google Earth Engine. *Remote Sensing*, 10(5), 691. https://iopscience.iop.org/article/10.1088/1748-9326/aa9d9e

Replicating that calibration from scratch (FIA plot matching + spectral extraction + model training + validation) is a separate research project, not a task to embed in a causal inference paper.

**Archived files** (in `archive/nbr_landsat_approach/`):
- `01_extract_biomass_gee_nbr.py` — Landsat annual composite builder + NBR extraction via GEE
- `test_gee_debug_nbr.py` — GEE connectivity debug script
- `02_biomass_exploration_nbr.qmd` — EDA document for NBR time series (Plots 5 & 6)
- `02_biomass_exploration_nbr.html` — Rendered HTML of the above
- `biomass_timeseries_nbr.csv` — Output CSV from the pilot CA extraction

**Next step:** Exploring the Ctrees biomass dataset as an alternative outcome variable that provides pre-calibrated annual aboveground biomass in Mg/ha.

---

## 2026-04-15 — Technical Decisions: `01_exploration.qmd`

### `sf_use_s2(FALSE)` — use GEOS instead of s2 for spatial operations
MTBS fire perimeter polygons contain self-touching edges that the s2 spherical geometry engine rejects as invalid, causing `st_join()` to error. Switching to GEOS (planar geometry) resolves this. For coterminous US data the planar approximation is accurate enough and `st_make_valid()` alone is not sufficient to satisfy s2's stricter validity rules. **Apply this at the top of any R script that does spatial joins on MTBS data.**

### `here::i_am("analysis/01_exploration.qmd")` — anchor project root for Quarto docs
Quarto renders documents with the working directory set to the document's own folder (`analysis/`), so `here("data/...")` resolves to `analysis/data/...` — the wrong place. Calling `here::i_am()` with the document's path relative to the project root anchors `here()` correctly. **Every Quarto document in a subdirectory must include this call.**

### Auto-install block for `tigris` and `ggspatial`
These two packages were not pre-installed in the R environment. The QMD includes a `.install_if_missing()` helper at the top of the setup chunk so the document installs them automatically on first render. If the environment ever gets rebuilt, these two packages must be added.

---

## Decision Log

Use this section to record design decisions and their rationale once resolved.

| Date | Decision | Choice | Rationale |
|---|---|---|---|
| — | Minimum fire size | — | — |
| — | Severity inclusion | — | — |
| — | Treatment variable (categorical vs. dNBR) | — | — |
| — | Study region / ecoregions | — | — |
| — | Control site strategy | — | — |
| — | Fire complex handling | — | — |
| — | Spatial buffer size | — | — |

---

## Data Quality Findings

*Record any anomalies, missing fields, or unexpected values found during EDA.*

### MTBS Perimeters (`mtbs_perims_DD.shp`)
- [ ] Field names confirmed
- [ ] CRS confirmed
- [ ] NAs / duplicates checked
- Notes:

---

## Open Questions

Questions that have not yet been resolved and need investigation or a decision.

- Which ecoregion classification to use (EPA Level III vs. Bailey's)?
- Should fire complexes be treated as a single large fire or excluded?
- What is the actual GEE asset path for eMapR biomass? (`projects/eMapR/biomass` — verify in GEE catalog)
- Minimum pre-fire years needed for parallel trends test with Callaway-Sant'Anna?
- How to handle fires that re-burn the same site in a different year?

---

## Literature Notes

*Key findings and methodological details from papers relevant to design decisions.*

### Callaway & Sant'Anna (2021)
- Estimator: `att_gt()` for group-time ATTs, `aggte()` for aggregation
- Clean controls: "never treated" or "not yet treated" units
- Does not require balanced panel
- Pre-treatment parallel trends is testable with their placebo approach

### Bright et al. (2019)
- Predictive (not causal) — R² > 0.7 with random forest
- Useful benchmark for biomass signal magnitude

### Ilangakoon et al. (2026)
- GAM with space-for-time substitution — lacks formal causal identification
- Our study directly addresses this gap

---

## GEE / Python Notes

*Notes on Google Earth Engine setup and extraction issues.*

- GEE authentication: `earthengine authenticate` in terminal before running Python script
- eMapR biomass asset path needs verification in GEE catalog before running extraction
- LandTrendR reference: https://emapr.github.io/LT-GEE/
- LandTrendR key parameters (defaults): `maxSegments`, `spikeThreshold=0.9`, `recoveryThreshold=0.25`, `pvalThreshold=0.1`, `minObservationsNeeded=6`

---

## Meeting / Advisor Notes

*Notes from advisor meetings or collaborator discussions.*

---
