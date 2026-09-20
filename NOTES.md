# Project Notes — Wildfire Biomass Recovery

Decisions, findings, and durable lessons for the project. Where things live:

| Section | What goes here |
|---|---|
| [Technical gotchas](#technical-gotchas) | Undated, durable lessons — check here before debugging |
| [Design decisions & open questions](#design-decisions--open-questions) | Single copy of what's decided and what isn't |
| [Dated log](#dated-log-newest-first) | What happened and why, newest first |
| [Literature notes](#literature-notes) | Papers and methods details |

Operating instructions (how to download, crop, extract, run on GRIT) live in
`DATA_DOWNLOAD_GUIDE.md`; project orientation and current status live in `CLAUDE.md`.

---

## Technical gotchas

### `file.exists()` ≠ "file is valid"
A background write killed mid-`writeRaster()` (laptop sleep, dropped connection, OOM) leaves a file
that is **present with a correct header but corrupt or truncated content**. Any skip-safe script that
only checks `file.exists()` silently treats it as done. This caused real damage several times: a
100%-NA WY forest mask, a truncated ctrees `2018` GeoTIFF (15.6 MB instead of ~730 MB), disappearing
eMapR west-crop years, and (on GRIT, 2026-09) a 7 KB header-only `ctrees_biomass_west_1km.nc`.
Validate before trusting an existing output, and delete + rebuild if it fails:

```r
raster_is_valid <- function(path) {
  terra::global(terra::rast(path), "notNA")[[1]] > 0
}
```

- Used in `05_prepare_forest_masks_west.R` and `00_crop_emapr_to_west.R`. It only catches *fully* corrupt
  (all-NA) files — for partial truncation, compare file size against sibling years/states or check
  extent/cell count.
- Python: `netcdf_is_valid()` in `04_download_ctrees_west.py` mirrors `raster_is_valid()` for the 1 km NetCDF.
  Part A's raw TIFs are still existence-checked only — check the last-written year's size against the
  others (~730 MB) after any interrupted run.
- Raw eMapR composites: `scripts/r/check_raw_emapr_files.R` checks exact expected byte count plus a small
  centered pixel-block read.

### terra on large rasters: only per-polygon crop → mask → mean works
On ~125M-cell ctrees TIFs, every whole-raster masking approach (`mask()`, `r * fm`, `extract()` on two
full rasters) is extremely slow or OOMs, even when the file is file-backed. Only a polygon-by-polygon
`crop()` → `mask()` → `mean()` loop works, because each `crop()` reads only that polygon's disk blocks.
This is why `07`/`08` extract per polygon instead of masking a whole state up front. The bottleneck is
random compressed disk I/O, not code structure — don't "optimize" this into a whole-raster operation
(re-derived and confirmed several times).

### Never materialize a full array (Python, GRIT memory cap)
GRIT sessions run under a Slurm cgroup memory limit. `04_download_ctrees_west.py` was OOM-killed three
times, each from holding a full ~2 GB year array (or several ~439 MB copies) in memory. The pattern that
fixed all three: stream — row-strip writes, windowed reads, one-year-at-a-time NetCDF writes with plain
`netCDF4`. See the 2026-09-20 log entry for the incident chronology and `DATA_DOWNLOAD_GUIDE.md` §3.3 for
what to do if it recurs. (An earlier "4 GiB `ulimit -m`" diagnosis was a red herring — `RLIMIT_RSS` isn't
enforced on modern Linux.)

### Cache validators must key on fire-set identity, not just years
In the (retired) CA-only pipeline, changing which fires are included (a border-fire filter fix) did not
invalidate an existing extraction CSV because the validator only checked which *years* were present, not
which `event_id`s. If fire-selection logic changes in any cache-based script, delete the cache manually or
extend the validator to include fire-set identity.

### R / Quarto conventions
- **`sf_use_s2(FALSE)`** at the top of any script that spatially joins MTBS data. MTBS perimeters have
  self-touching edges that the s2 engine rejects; `st_make_valid()` alone doesn't satisfy s2. Planar GEOS is
  accurate enough for coterminous-US data.
- **`here::i_am("analysis/<file>.qmd")`** in every Quarto document in a subdirectory — Quarto renders with the
  document's folder as working directory, so `here("data/...")` would otherwise resolve to `analysis/data/...`.
- **`tigris`, `ggspatial`** weren't pre-installed; `01_mtbs_exploration.qmd` auto-installs them via
  `.install_if_missing()` in its setup chunk. Add them if the R environment is ever rebuilt.
- **Extraction scripts run from `Rscript` or the R console, never a Quarto chunk** (Quarto buffers output, so
  terra's threading looks frozen).

---

## Design decisions & open questions

### Decided

| Date | Decision | Choice | Rationale |
|---|---|---|---|
| 2026-05-12 | Biomass outcome variable | eMapR + ctrees pre-calibrated AGB (Mg/ha); NBR approach archived | NBR is a unitless spectral index; calibrating it to biomass is a separate research project (see log) |
| 2026-08 | Forest mask | NLCD 2004, per-state, 0/1 fraction-forest (classes 41/42/43) | Replaces retired CA-only 1/NA mask; built by `05_prepare_forest_masks_west.R` |
| 2026-08-30 | Raw eMapR retention | Keep all raw CONUS composites permanently (~1 TB) | PI request; requires GRIT-scale storage, not a laptop |
| 2026-09 | Where data lives | GRIT is the durable store; Nextcloud plan dropped | See 2026-08-30 log entry |
| — | Time window | As long as the data allows (eMapR 1990–2023, ctrees 2000–2025, MTBS fires 2000–2023) | Long pre-fire baselines for parallel-trends testing |

### Leaning, not final
- **Treatment variable:** RdNBR (continuous; = dNBR / √|preNBR/1000|, Miller & Thode 2007) preferred over raw
  dNBR for cross-fire comparison — vs. categorical severity classes.

### Open
- Minimum fire size threshold?
- Include moderate severity, or only high?
- Which ecoregion classification (EPA Level III vs. Bailey's), and which ecoregions to include?
- Control-site strategy (never-burned, same ecoregion — details TBD; see `data/processed/control_pixels/`)?
- Fire complexes: treat as one large fire or exclude?
- Spatial buffer between sites (SUTVA)?
- Minimum pre-fire years needed for a parallel-trends test with Callaway-Sant'Anna?
- Re-burns: how to handle sites that burn again in a different year?
- **Residual eMapR-vs-ctrees bias** (`biomass_within_fires.qmd` §7) — root cause not yet found.
- **Border-fire eMapR/ctrees count gap** — see 2026-08-30 entry; re-check once real multi-state `07` output exists.
- **MTBS Initial vs. Extended assessment bias** — `mtbs_assessment_comparison.qmd`.

---

## Dated log (newest first)

### 2026-09-20 — GRIT migration: layout confirmed, ctrees West download hardened

**GRIT layout, as set up:** code at `~/BACI-wildfire`; data lives in a separate shared repo `~/BACI-review`,
joined by `ln -s ~/BACI-review/data ~/BACI-wildfire/data`, so every script's `here()`/`PROJ_ROOT`-relative
path works unchanged. Python env is a venv at `~/BACI-wildfire/.venv`, not the laptop's conda env. Setup steps:
`DATA_DOWNLOAD_GUIDE.md` Part 1.

**Raw eMapR validity checker added** (`scripts/r/check_raw_emapr_files.R`, commits `d54b32b`, `618aad1`).
19/34 years confirmed complete on GRIT as of 2026-09-08.

**`04_download_ctrees_west.py` — three OOM kills, three fixes** (commits `bcbc41b`, `8ad9dfe`, `4eeb88f`):
1. Coarsening's whole-array `reshape()` silently copied the full ~2 GB array because West's dimensions
   (25650 cols / 11) aren't an exact multiple of the coarsen factor.
2. Part A's plain raw download held one ~2 GB float32 year array on its own.
3. Part B's NetCDF assembly held up to three ~439 MB copies at once (list + `np.stack()` + `.astype()`).

The script was also reordered to run raw-download-first (A → B → C) so B/C read plain local GeoTIFFs instead
of re-querying arraylake, and `rasterio` became a hard requirement.

**Then the `file.exists()` gotcha hit GRIT for real.** A relaunch was started as a bare foreground command,
not under `tmux`; a dropped connection killed it mid-Part-B-write and left a 7 KB corrupted `.nc`. Diagnosed
live: 26/26 Part A TIFs present, `.nc` header-only, no `_west_fireagb_scratch/` (Part C never started), no
`tmux` session, `memory.events` showing no OOM.

**Fixes applied before relaunching** (commit `c931dcd`):
- `netcdf_is_valid()` added — Part B now deletes and rebuilds a corrupt `.nc` instead of trusting
  `OUT_NC.exists()`.
- Part C hardened pre-emptively (it had never run at West scale): per-fire `rasterio.Window` reads instead of
  a full ~2 GB array per year; mask cache shrunk with `np.packbits`; the previously dead `errors` list is now
  populated via per-fire `try/except`.
- End-of-script sanity check no longer loads the whole ~439 MB NetCDF (same bug as Part B's assembly).
- `log_peak_memory()` helper (`/proc/self/status` → `VmHWM`) after each part, so a future OOM is diagnosable
  from the log.

**Status:** code fixes committed; not yet re-run on GRIT. Next: delete the corrupt `.nc`, pull, relaunch under
`tmux` (command block in `DATA_DOWNLOAD_GUIDE.md` §3.3), and cross-validate the West CSV's CA rows against
the validated `03_download_ctrees_ca.py` baseline (expect correlation ≈ 1.000).

---

### 2026-08-30 — Laptop → GRIT handoff: decisions worth keeping

Consolidated when the project moved from the local laptop to GRIT. Pipeline *status* from that date is
superseded by `CLAUDE.md` "Current Status"; the decisions below still stand.

**`07`/`08` rewrite — design.** Rewritten from hardcoded `STATE_FIPS <- "CA"` to a `STATES_TO_RUN` vector,
reading shared West-wide crops instead of per-state files:
- **eMapR (`07`)** reads the native ~30m West-cropped TIF directly and does crop-to-polygon →
  `aggregate(fact=3)` → mask per fire — no separate whole-West downsampled file (avoids an extra multi-GB/year
  artifact).
- **ctrees (`08`)**: same algorithm as before, pointed at the shared `ctrees_YYYY_west_100m.tif`.
- **MTBS dedup:** both use event_id-prefix dedup (spatial join to the West states union, then
  `STUSPS = substr(event_id, 1, 2)` as authoritative) instead of single-state `st_filter()` +
  `startsWith()` — the fix for the border-fire mechanism below.
- **Output:** one combined CSV per script with a `STUSPS` column, resumable per (state, year).
- Per-polygon extraction is wrapped in `tryCatch()` so one bad polygon logs a warning and returns `NA`
  instead of crashing a multi-hour run.
- **Validation harness** (`scripts/r/validate_west_pipeline.R` + `analysis/west_pipeline_sanity_check.qmd`,
  re-run with `Rscript scripts/r/validate_west_pipeline.R`) gated the rewrite. It caught a corrupted WY mask
  and a CRS bug (missing `st_transform(5070)`) in the harness's own comparison code — not in `07`/`08`.
  Result: ctrees matches the retired baseline exactly (0% difference across 1,632 fire×year pairs); the CA
  dedup selects the identical 304-fire 2005–2010 set (0 disagreements); the WY smoke test gave 76 well-formed
  fires. **eMapR re-validation pending:** the first attempt returned all-NA because of the harness CRS bug
  (fixed). To finish it, delete `data/processed/validation/emapr_method_comparison_ca.csv` and `..._diff.csv`
  (they hold the bad run) and re-run the script — the ctrees output is skip-safe.

**Border-fire eMapR/ctrees count gap (open).** In the old CA-only pipeline, ctrees matched all 272 study fires
(2005–2010 cohort) but eMapR only 234 — the same ~38 fires missing consistently. Leading hypothesis: eMapR's
old CA extraction was polygon-masked to the CA state boundary while ctrees' wasn't, so a border fire with most
of its forested area outside CA lost its entire eMapR value to NA. The `07`/`08` dedup rewrite targets that
mechanism and the harness confirmed old-vs-new fire *selection* matches (0 disagreements), but that doesn't
confirm the eMapR *data availability* gap is closed. Re-check the ctrees-vs-eMapR fire count directly once
real multi-state output exists. A secondary hypothesis to rule out alongside it: eMapR's forest mask is coarser
(~90 m) than ctrees' (~100 m), so a small fire with only a sliver of forest could lose its only forest pixel
under the coarser mask regardless of the border issue. Diagnostic: take the `event_id`s present in ctrees but
not eMapR and check whether they cluster near the state border.

**Nextcloud storage architecture — dropped.** The plan was raw eMapR composites (~1 TB) living on Nextcloud,
fetched transiently via `rclone` to crop then discard. WebDAV credentials were never obtained and the
GRIT migration supersedes it: GRIT is the durable store (raw files are kept, per PI request).

**Retired CA-only pipeline** (`00_crop_emapr_to_ca.R` → `02_extract_emapr_within_fires.R` /
`04_extract_ctrees_within_fires.R` → `analysis/biomass_within_fires_old.qmd`): kept only as the validation
baseline the `05`–`08` pipeline was cross-checked against. Do not extend it for new work.

**eMapR West-crop interruptions (laptop).** `00_crop_emapr_to_west.R` was repeatedly killed on the laptop
(first by lid-close sleep, later with no confirmed cause), leaving unverified `composite_YYYY_west.tif` files.
Re-verify any West-crop file on GRIT rather than assuming earlier laptop output carried over.

---

### 2026-08-12 — ctrees West download: root causes from three interrupted laptop runs

- **Lid-close** triggered Modern Standby regardless of the idle-sleep setting, killing background jobs (same
  failure as the eMapR crop). Fixed on the laptop via `powercfg` `SUB_BUTTONS LIDACTION` = 0; only relevant if
  long runs move back to a laptop.
- **Part B had no per-year checkpointing** — a kill lost 16/26 years of progress. It now checkpoints to
  `data/processed/ctrees/_west_fireagb_scratch/` per year, like Part A's `_west_1km_scratch/`.
- **Part C's existence-only check** left a truncated `ctrees_2018_west_100m.tif` (15.6 MB vs ~730 MB) that had
  to be deleted by hand — the origin of the "validate, don't trust `file.exists()`" gotcha above.

---

### 2026-05-13 — Fire polygon extraction: shapely thinning replaced by rasterio rasterization

Part B of `03_download_ctrees_ca.py` extracts mean ctrees AGB within each MTBS CA fire polygon per year.

- **Abandoned — shapely point-in-polygon with grid thinning.** Large fires (e.g. the 2020 August Complex,
  ~1M acres) have bounding boxes with millions of pixels; even thinned to ≤80,000 test points the
  1,064-fire × 26-year loop ran ~1 year per 5 minutes (~130 min projected) and under-sampled boundaries.
- **Current — `rasterio.features.rasterize()` + precomputed masks.** A C-level scanline burn onto a boolean
  grid aligned to the pixel coordinates: no sampling, handles Polygon/MultiPolygon, and masks are computed
  once for all fires before the year loop, so the 26-year extraction is pure numpy indexing
  (`raw[np.ix_(yi, xi)][mask]`). ~10–50x faster for large polygons. (A `matplotlib.path` fallback for machines
  without rasterio existed at the time; `04_download_ctrees_west.py` has since made rasterio mandatory.)

---

### 2026-05-12 — NBR/Landsat GEE extraction approach archived

The first extraction pipeline built raw annual Landsat composites in Google Earth Engine and computed **NBR**
as a biomass proxy. Archived because **NBR is a unitless spectral index (−1 to +1), not a biomass
quantity** — DiD effects in NBR units lack ecological interpretability and won't satisfy ecology/fire-science
reviewers who expect biomass (Mg/ha) or carbon (MgC/ha).

Converting NBR to biomass needs an empirical calibration (typically random forest on co-located FIA plots,
plus climate/topography to handle saturation above ~300 Mg/ha) — exactly what the eMapR lab built:

> Kennedy, R.E., Yang, Z., Gorelick, N., Braaten, J., Cavalcante, L., Cohen, W.B., & Healey, S. (2018).
> Implementation of the LandTrendr algorithm on Google Earth Engine. *Remote Sensing*, 10(5), 691.
> https://iopscience.iop.org/article/10.1088/1748-9326/aa9d9e

Replicating that calibration is a separate research project. **Archived in `archive/nbr_landsat_approach/`:**
`01_extract_biomass_gee_nbr.py`, `test_gee_debug_nbr.py`, `02_biomass_exploration_nbr.qmd`/`.html`,
`biomass_timeseries_nbr.csv`. LandTrendR reference: https://emapr.github.io/LT-GEE/ (defaults used:
`spikeThreshold=0.9`, `recoveryThreshold=0.25`, `pvalThreshold=0.1`, `minObservationsNeeded=6`).

Data-access notes from the same day (scratch note, 2026-05-12): eMapR's own download is a GUI (one grid cell,
one year — impractical); direct GEE access yields NBR, not biomass. Chosen approach: download all years for a
region/state from the eMapR FTP, and compare against ctrees filtered to the same state.

The underlying goal (long pre-fire baselines for parallel-trends testing) is met by eMapR's 1990 start year
without the archived approach.

---

## Literature notes

### Callaway & Sant'Anna (2021) — staggered DiD
- `att_gt()` for group-time ATTs, `aggte()` for aggregation.
- Clean controls: "never treated" or "not yet treated" units.
- Does not require a balanced panel.
- Pre-treatment parallel trends testable via their placebo approach.

### Other methods papers
Goodman-Bacon (2021) — TWFE decomposition; Sun & Abraham (2021) — interaction-weighted estimator;
Liermann & Roni (2021) — staircase design power analysis.

### Bright et al. (2019)
Predictive (not causal); random forest, R² > 0.7. Useful benchmark for biomass signal magnitude.

### Ilangakoon et al. (2026)
GAM with space-for-time substitution — lacks formal causal identification. Our study addresses this gap.

### Other domain papers
Garcia et al. (2017) — Rim Fire carbon; Reisch (2024); Stenzel (2019).

### MTBS methodology
Eidenshink, J., Schwind, B., Brewer, K., Zhu, Z. L., Quayle, B., & Howard, S. (2007). A project for monitoring
trends in burn severity. *Fire Ecology*, 3(1), 3–21. Dataset DOI: https://doi.org/10.5066/P9IED7RZ. Covers NBR/dNBR
mapping, per-fire threshold calibration from unburned reference areas, the five severity classes (including
Increased Greenness), the ≥1,000 ac western US size threshold, and coverage since 1984.

### RdNBR
Miller, J. D., & Thode, A. E. (2007). Quantifying burn severity in a heterogeneous landscape with a relative
version of the delta Normalized Burn Ratio (dNBR). *Remote Sensing of Environment*, 109(1), 66–80.
RdNBR = dNBR / √|preNBR / 1000|; normalizes for pre-fire vegetation density.

### Our contribution
First application of modern staggered and continuous DiD to wildfire–biomass; relaxes the untestable
conditional-independence assumption of prior work and exploits the natural staggered timing of fires.
