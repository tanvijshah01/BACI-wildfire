# Project Notes — Wildfire Biomass Recovery

Working notes on decisions, findings, and open questions. Add entries in reverse chronological order (newest at top).

---

## CLAUDE CODE MEMORY NOTES — Handoff for Positron Assistant (2026-08-30)

**What this section is:** Claude Code (the AI assistant previously used on this project, running
on a local laptop) keeps a persistent memory system *outside* this git repo, tied to that specific
machine. As the project moves to running on the GRIT server — data via Nextcloud/FileZilla, code
and assistant work happening on GRIT itself, likely via Positron Assistant instead of Claude Code —
that laptop-local memory doesn't travel automatically. This section consolidates everything from it
that's still relevant, so a new assistant (or a person) picking this up on GRIT has the context
without re-deriving it. Below is organized by how current/actionable each piece still is.

### Current data-pipeline status (as of 2026-08-30, verify before trusting — see gotcha below)

- **eMapR raw archive** (`data/raw/emapr_biomass/composite_YYYY_median.tif`): 20 of 34 years
  (1990–1996, 2000–2012) were present on the laptop. Years 1997–1999 and 2013–2023 were never
  fetched. The plan to archive these durably on Nextcloud was never completed (WebDAV credentials
  never obtained) — see "Nextcloud architecture — superseded" below. Whatever state this is in now
  depends entirely on what got uploaded to GRIT directly; it should not be assumed any of this
  eMapR raw/cropped work carried over — check `data/raw/emapr_biomass/` and
  `data/processed/emapr_biomass_west/` on GRIT directly.
- **eMapR West-crop** (`scripts/r/00_crop_emapr_to_west.R`, → `data/processed/emapr_biomass_west/composite_YYYY_west.tif`):
  repeatedly interrupted on the laptop by background-job kills (originally laptop lid-close sleep,
  fixed via a `powercfg` lid-action setting; later kills had no confirmed root cause — possibly
  session/environment teardown, not reproducible). Last confirmed local state: 1990 valid,
  1991 freshly rebuilt and valid, 1992 present but **unconfirmed** — sized suspiciously smaller than
  every other year and was mid-write when the job died, never verified before the conversation moved
  on to the GRIT migration. **Do not trust any `composite_YYYY_west.tif` file's mere presence** — see
  the corruption gotcha below.
- **ctrees**: CA-only download (`scripts/python/03_download_ctrees_ca.py`) is fully complete
  (26 years, 2000–2025). The West-wide version (`scripts/python/04_download_ctrees_west.py`,
  ~4.1x CA's pixel area, ~6,800 fires) had Parts A (1km NetCDF) and B (fire-polygon CSV) fully
  complete, and Part C (26 annual 100m GeoTIFFs) at 18/26 (years 2000–2017 done) the last time it
  was directly checked on the laptop — this was not rechecked again before the GRIT pivot, so
  treat it as possibly stale.
- **`07`/`08` extraction scripts**: rewritten to support multi-state `STATES_TO_RUN` instead of
  hardcoded `STATE_FIPS <- "CA"` (see "07/08 rewrite" below), and validated on CA+WY. Not yet run
  end-to-end for real production `STUDY_YEARS` output, since the West eMapR crop and West ctrees
  TIFs weren't both complete yet when work paused.

### Recurring technical gotcha: `file.exists()` ≠ "file is valid"

This machine (the laptop) had a repeated failure mode: a background raster-writing job gets killed
mid-`writeRaster()` (originally laptop-sleep, later unclear causes), leaving a GeoTIFF that's
**present on disk with a correct header but corrupt or truncated content** — either all-NA pixels,
or a file that just stops partway through. Scripts that only checked `file.exists()` to decide
whether to skip re-processing a year/state silently treated these as done. This caused real damage
three times: a corrupted WY forest mask (100% NA), a truncated ctrees `2018` GeoTIFF (15.6 MB
instead of ~730 MB), and the eMapR west-crop years 1991–1996/2000 disappearing between sessions
(rebuilt once, corruption status of 1992 still unconfirmed as of this handoff).

**Fix pattern applied** (in `05_prepare_forest_masks_west.R` and `00_crop_emapr_to_west.R`, and
worth applying anywhere else this pattern shows up): a `raster_is_valid()` helper —
```r
raster_is_valid <- function(path) {
  terra::global(terra::rast(path), "notNA")[[1]] > 0
}
```
— checked before trusting an existing output file, deleting and rebuilding if it fails. Note this
only catches *fully* corrupt (all-NA) files, not partial truncation with some real data still
readable — for that, compare file size against other years/states as a sanity check, or verify the
raster opens and has the expected extent/cell count.

**Whether this gotcha still applies on GRIT is unknown** — it may have been specific to the
laptop's sleep/background-task behavior. Worth being alert for it early on rather than assuming a
shared HPC server is immune.

### Other durable technical lessons

- **terra performance on large rasters** (from ctrees CA-scale, ~125M-cell TIFs): applying a
  forest mask via any "whole-raster" approach (`mask()`, `r * fm`, `extract()` on two full rasters)
  is extremely slow or OOMs, even though the file itself is manageable when file-backed and
  accessed only by cropped windows. **Only a polygon-by-polygon `crop()` → `mask()` → `mean()` loop
  works** — each `crop()` reads only that polygon's disk blocks. This is why `07`/`08` extract
  per-polygon rather than masking a whole state's raster up front. Don't try to "optimize" this into
  a single whole-raster operation — the bottleneck is random compressed disk I/O, not code
  structure, and this has been re-derived and confirmed multiple times.
- **Cache invalidation checks year coverage, not fire-set membership**: in the (now-retired)
  CA-only pipeline, changing which fires are included (e.g. a border-fire filter fix) did not
  invalidate an existing extraction CSV cache, because the cache validator only checked which years
  were present, not which `event_id`s. If a future cache-based script changes fire-selection logic,
  check whether existing caches need manual deletion, or whether the cache key/validator should
  include fire-set identity.

### Open, not-yet-confirmed item

**Border-fire eMapR/ctrees count gap**: in the old CA-only pipeline, ctrees matched all 272 study
fires (2005–2010 cohort) but eMapR only matched 234 — same ~38 fires missing consistently, not
random noise. Leading hypothesis: eMapR's old CA extraction was polygon-masked to the CA state
boundary while ctrees' wasn't, so a border fire whose forested area sits mostly outside CA could
lose its entire eMapR value to NA. The `07`/`08` rewrite's event_id-prefix MTBS dedup (see below)
was built partly to address the underlying mechanism, and a validation harness run confirmed old
vs. new fire *selection* matches exactly (0 disagreements) — but that's not the same as confirming
the eMapR *data availability* gap itself is closed, since `07` hadn't produced real multi-state
production output yet. Re-check the ctrees-vs-eMapR fire count gap directly once real output exists.

### `07`/`08` rewrite — design decisions

Rewritten from hardcoded `STATE_FIPS <- "CA"` to a `STATES_TO_RUN` vector, reading from shared
West-wide crops instead of per-state files:
- **eMapR (`07`)**: reads the native ~30m West-cropped TIF directly, does
  crop-to-polygon → `aggregate(fact=3)` → mask, per fire polygon — no separate whole-West
  downsampled precompute file (chosen deliberately to avoid an extra multi-GB/year file).
- **ctrees (`08`)**: same extraction algorithm as before, just pointed at the shared
  `ctrees_YYYY_west_100m.tif` (from `04_download_ctrees_west.py` Part C) instead of per-state files.
- **MTBS dedup**: both now use event_id-prefix dedup (spatial join to West states union, then
  `STUSPS = substr(event_id, 1, 2)` as authoritative) instead of single-state `st_filter()` +
  `startsWith()` — this is the fix for the border-fire mechanism above.
- **Output**: one combined CSV per script with a `STUSPS` column, resumable per (state, year), not
  one file per state.
- Both scripts wrap their per-polygon extraction in `tryCatch()` so one bad polygon logs a warning
  and returns `NA` instead of crashing a multi-hour, multi-state run.
- A validation harness (`scripts/r/validate_west_pipeline.R` +
  `analysis/west_pipeline_sanity_check.qmd`) was built specifically to gate this rewrite before
  trusting it at scale — resumable per (state, year), re-run with
  `Rscript scripts/r/validate_west_pipeline.R`. This is what caught the corrupted WY mask and a
  CRS bug (missing `st_transform(5070)`) in the harness's own comparison code — not in `07`/`08`
  themselves, which were confirmed to already transform correctly.

### Historical / superseded — context only, not current instructions

- **Nextcloud storage architecture**: the original plan was raw eMapR composites (~1 TB total)
  living durably on Nextcloud, fetched transiently via `rclone` to crop then delete locally. This
  was **never fully set up** (WebDAV credentials were never obtained) and is now superseded by the
  GRIT migration itself — GRIT is the new durable-storage answer, not Nextcloud-as-archive. Nextcloud
  (or FileZilla) may still be the *upload mechanism* to get data onto GRIT, per the current plan,
  but the "laptop ⇄ Nextcloud ⇄ crop-then-discard" architecture is not what's being built anymore.
- **Old NBR/Landsat GEE biomass approach** (`archive/nbr_landsat_approach/`, `START_YEAR = 1984`):
  archived — NBR is a unitless spectral index, not a biomass quantity, and calibrating it to
  biomass would be a separate research project. Current pipeline uses eMapR and ctrees, which
  provide pre-calibrated AGB directly. The underlying reason this came up (wanting long pre-fire
  baselines for parallel-trends testing) is naturally satisfied by eMapR's 1990 start year without
  needing this archived approach.
- **`reference_python_env`**: on the laptop, ctrees Python scripts needed the anaconda3 *base* env
  specifically (not the `wildfire` conda env, which was missing `arraylake`/`geopandas`/`rasterio`).
  This is laptop-specific and almost certainly **does not apply on GRIT** — GRIT will have its own
  Python/conda setup that needs to be checked independently; don't assume an env named `wildfire`
  or `base` means anything equivalent there.
- **Old CA-only retired pipeline** (`00_crop_emapr_to_ca.R` → `02_extract_emapr_within_fires.R` /
  `04_extract_ctrees_within_fires.R` → `analysis/biomass_within_fires_old.qmd`): kept only as a
  validation baseline the `05`–`08` pipeline was cross-checked against (ctrees matched exactly,
  r = 1.000). Do not extend this retired pipeline for new work.

---

## 2026-08-12 — ctrees West-wide download in progress, paused for the night

`scripts/python/04_download_ctrees_west.py` (generalizes `03_download_ctrees_ca.py` to
the 11-state West bbox, ~4.1x CA's pixel area, ~6,800 fires) is partway through. Progress
as of tonight:

- **Part A** (coarsened 1km NetCDF): ✅ complete — `ctrees_biomass_west_1km.nc` (315 MB)
- **Part B** (fire-polygon CSV): ✅ complete — `biomass_fire_polygons_ctrees_west.csv`
  (177,242 records = 6,817 fires × 26 years, checks out exactly)
- **Part C** (26 annual ~100m GeoTIFFs): **18/26 done** (2000–2017, ~730 MB each). Years
  **2018–2025 still need to run.**

**Resume command** (from project root, sleep-safe now — see below):
```powershell
& "C:\Users\shaht\anaconda3\python.exe" "scripts\python\04_download_ctrees_west.py" *>> "data\processed\ctrees\04_download_ctrees_west_log.txt"
```
Uses the `anaconda3` base env, not the `wildfire` conda env (which lacks `arraylake`/
`geopandas`/`rasterio`). Both Parts A and B will skip instantly (already done); Part C
resumes at 2018.

**Root cause of 3 interrupted runs tonight:** the laptop's **lid-close action** was
triggering Modern Standby regardless of the `standby-timeout-ac` idle setting (which was
already disabled) — closing the lid kills the background process outright, same failure
mode already documented for the eMapR west-crop. Fixed by unhiding and setting the lid
power setting directly:
```powershell
powercfg -attributes SUB_BUTTONS LIDACTION -ATTRIB_HIDE
powercfg /setacvalueindex SCHEME_CURRENT SUB_BUTTONS LIDACTION 0   # 0 = Do nothing
powercfg /setactive SCHEME_CURRENT
```
Confirmed via `powercfg /query SCHEME_CURRENT SUB_BUTTONS` → `Current AC Power Setting
Index: 0x00000000`. This should already be in place — verify it's still set before
walking away from a future long run, since `powercfg` state doesn't survive an OS
reinstall/reset.

**One more finding from tonight, already fixed in the script:** Part B originally had no
per-year checkpointing (only Part A did) — the second kill lost 16/26 years of Part B
progress because it held everything in memory and wrote the CSV once at the end. Part B
now checkpoints to `data/processed/ctrees/_west_fireagb_scratch/` per year, same pattern
as Part A's `_west_1km_scratch/`. Part C was already safe (checks each year's TIF file
before writing) — but note that check is existence-only, not validity: the 3rd kill left
a **truncated `ctrees_2018_west_100m.tif`** (15.6 MB instead of ~730 MB) that had to be
manually deleted before resuming, since Part C would have silently treated it as done. If
a run gets killed again mid-Part-C, check the last-written year's file size against the
others (~730 MB) before trusting it and resuming.

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
