# EDA Plan — Wildfire Biomass Recovery

**Phase:** Exploration → West-wide expansion
**Last updated:** 2026-09-24
**Goal:** Understand the MTBS, eMapR, and ctrees data well enough to make the design decisions in
`NOTES.md` (→ Design decisions & open questions) before building the analysis panel.

**Status legend:** ✅ done · 🔄 in progress / partly done · ⬜ not started · ⛔ stale (superseded)

Where this fits: `CLAUDE.md` = project orientation and pipeline status · `DATA_DOWNLOAD_GUIDE.md` = how to
run the pipeline · `NOTES.md` = decisions, hurdles, history · **this file = what EDA we're doing and how far
along each piece is.** Analysis documents live in `analysis/`; figures land in `figures/` (gitignored).

---

## 1. Core EDA questions

| # | Question | Status | Where answered / next step |
|---|---|---|---|
| 1 | How many MTBS fires in 2000–2023, Western US? | ✅ | `01_mtbs_exploration.qmd` §3 (row count printed at each filter step) |
| 2 | Distribution of fire sizes? | ✅ | `01` §3.6 (log₁₀ histograms, all vs. Extended) |
| 3 | Severity classes and their distribution? | ✅ | `01` §3.7 (BSP CBI dNBR/RdNBR). Note: the MTBS perimeter shapefile holds dNBR *threshold values*, not per-class area, so % low/moderate/high per fire is not available from it |
| 4 | Temporal trends in fire occurrence? | ✅ | `01` §3.2, §3.5 |
| 5 | Spatial clustering of fires? | ⬜ | Needed for the SUTVA / spatial-buffer decision |
| 6 | Fires per year (treatment timing)? | ✅ | `01` §3.2; per-cohort view in `biomass_within_fires.qmd` §2 |
| 7 | Which years have enough pre-fire data for parallel-trends testing? | ⬜ | Data-availability facts: eMapR starts 1990 (≥10 pre-fire years for 2000 fires), ctrees starts 2000 (none for 2000 fires). Quantify per cohort once the panel exists |
| 8 | Do ctrees and eMapR agree within fires? | 🔄 | `biomass_within_fires.qmd` §5, §7 — residual eMapR bias not root-caused |
| 9 | Is MTBS Initial vs. Extended assessment type spatially/seasonally biased? | 🔄 | `mtbs_assessment_comparison.qmd` — unresolved |

---

## 2. Analysis documents

| Document | Purpose | Status | Notes / open items |
|---|---|---|---|
| `01_mtbs_exploration.qmd` | MTBS perimeters, assessment types, sizes, timing, CBI severity | ✅ §1–3 · ⛔ §4 | §4 (GEE NBR extraction, 4.1/4.2) belongs to the archived NBR approach — remove or replace with a pointer to `biomass_within_fires.qmd` §6 (see `NOTES.md` 2026-05-12) |
| `02_ctrees_biomass_exploration.qmd` | CA ctrees AGB: annual trend, maps (2000/2020/2025), before/after fire event study | ✅ (CA) | Uses `03_download_ctrees_ca.py` outputs; no West version |
| `03_emapr_biomass_exploration.qmd` | CA eMapR AGB for 2000/2001/2003: distributions, year comparison, maps | ✅ (CA) | Reads retired CA crops (`emapr_biomass_ca/`) |
| `04_data_summary.qmd` | Side-by-side summary of MTBS, ctrees, eMapR (CA); dataset comparison table; §7 design decisions (assessment types, outcome variable, control-unit definition & SUTVA) | ✅ (CA) | Depends on the retired `03_prepare_forest_mask.R` mask. §7 decisions should be reconciled with `NOTES.md` |
| `biomass_within_fires.qmd` | **Current.** ctrees vs. eMapR within forested fire pixels, CA cohorts; masking progression; before/after event study (§6); design notes (§7) | 🔄 | Fed by `07`/`08` (`STATES_TO_RUN` = CA, `STUDY_YEARS` 2005–2010). Open: eMapR residual bias (§7); extend to West + all years |
| `biomass_within_fires_old.qmd` | Retired CA-only forest-mask version | ⛔ | Kept only as the validation baseline |
| `mtbs_assessment_comparison.qmd` | Initial vs. Extended: fire counts by state, ignition month, % forest cover, forest-cover-threshold sweep | 🔄 | Forest-cover cache currently CA only (`06`); widen to West |
| `west_pipeline_sanity_check.qmd` | Regression report for the 05–08 multi-state rewrite (dedup, mask validity, eMapR/ctrees old-vs-new) | ✅ | Re-run after any change to 05–08 (`Rscript scripts/r/validate_west_pipeline.R`, then render) |

---

## 3. Data products feeding the EDA documents

The EDA documents read small per-fire CSVs produced by the extraction scripts, plus (for maps only) a few
large rasters. This section lists what each raster is, how it is built from raw data, and where processing
currently stands. Run commands and environment setup are in `DATA_DOWNLOAD_GUIDE.md`.

### 3.1 Raster files and what each is for

**Forest masks.** Each mask records whether a pixel is forest, from NLCD 2004 land cover (classes 41
Deciduous, 42 Evergreen, 43 Mixed → 1; everything else → 0). Several versions exist because each biomass
dataset sits on a different grid, and a mask is only usable if its pixels line up with the raster it is
applied to. All are built per state by `05_prepare_forest_masks_west.R` into `data/processed/forest_mask/`.

| File (`<st>` = state code, e.g. `ca`) | Grid | Purpose | Used by |
|---|---|---|---|
| `nlcd2004_forestfrac_30m_<st>.tif` | 30 m, EPSG:5070 (NLCD native) | Full-detail master mask | `06` (% forest per fire); `08` (reprojected per fire onto the ctrees grid) |
| `nlcd2004_forestfrac_90m_<st>.tif` | 90 m, EPSG:5070 | Matches eMapR after `07` aggregates it from 30 m to 90 m | `07` (eMapR AGB per fire) |
| `nlcd2004_forestfrac_100m_<st>.tif` | ~100 m, EPSG:4326 (ctrees grid) | Matches ctrees' native grid, so a whole-state map can be masked in one step | `biomass_within_fires.qmd` §6 only (maps, masking-progression figures) |

**Biomass rasters.**

| File | Contents | Used by |
|---|---|---|
| `data/raw/emapr_biomass/…` | eMapR AGB, 30 m, CONUS-wide (~27.7 GB/year) | `00_crop_emapr_to_west.R` only |
| `emapr_biomass_west/composite_<yr>_west.tif` | eMapR AGB, 30 m, clipped to the 11 Western states (~1 GB/year) | `07`; source for the CA display files below |
| `composite_<yr>_ca_100m.tif` | eMapR AGB aggregated to 90 m, CA only | `biomass_within_fires.qmd` §6 only |
| `ctrees/ctrees_<yr>_west_100m.tif` | ctrees AGB, native ~100 m, EPSG:4326, one shared raster for the whole West | `08`; `biomass_within_fires.qmd` §6 |

The 100 m mask and the CA-100m eMapR composites exist only for figures. None of the per-fire extractions
(`06`–`08`) need them.

### 3.2 Pipeline: raw data → EDA

```
RAW                            PROCESSING                                  EDA
MTBS fire perimeters ──────────────────────────────────────┐
                                                           │
NLCD 2004 (fetched by 05 ──► 05  forest masks 30/90/100 m ─┤
  from MRLC WCS)                                           ├─► 06  % forest per fire ───► mtbs_assessment_comparison.qmd
                                                           │
eMapR CONUS (~28 GB/yr) ──► 00_crop_emapr_to_west ─────────┼─► 07  eMapR AGB per fire ──┐
                              (composite_<yr>_west.tif)    │                            │
                                 └─► 01_create_emapr_100m (CA 90 m display) ─────────┐ ├─► biomass_within_fires.qmd
                                                           │                         │ │
ctrees zarr store ──► 04_download_ctrees_west.py ──────────┴─► 08  ctrees AGB per fire ┘ │
                        (ctrees_<yr>_west_100m.tif) ─────────────────────────────────────┘
```

In words:

1. **Acquire raw data.** MTBS perimeters are downloaded once. eMapR comes from the eMapR lab FTP server, one
   CONUS composite per year. ctrees is pulled from the arraylake zarr store by `04_download_ctrees_west.py`,
   which writes one West-wide ~100 m GeoTIFF per year.
2. **Reduce raster size.** `00_crop_emapr_to_west.R` clips each raw eMapR year to the 11 Western states.
   ctrees already arrives at West extent.
3. **Build forest masks.** `05` downloads NLCD 2004 per state and writes the 30 m, 90 m, and ~100 m masks
   described above.
4. **Extract per-fire values.** For each MTBS fire, `06` computes the fraction of forest inside the
   perimeter; `07` and `08` compute mean forest-only AGB per year from eMapR and ctrees. Each writes one
   CSV with a `STUSPS` column. All three work fire-by-fire (crop → mask → mean) to stay within memory.
5. **EDA.** `mtbs_assessment_comparison.qmd` reads the `06` CSV. `biomass_within_fires.qmd` reads the `07`/`08`
   CSVs for its cohort and event-study analyses, and reads the large rasters directly for the §6 maps.

### 3.3 Where processing stands (2026-09-24)

**Complete and verified on GRIT:**
- ctrees West download: 26/26 yearly rasters valid (2000/2001 rebuilt after being found 100% NaN).
- MTBS loading in `06`–`08` filtered at read time; `08` reproduces the validated CA baseline (304 fires).

**Blocked by GRIT's 4 GiB memory cap** (cgroup v2; see `NOTES.md` → Technical hurdles):

| Step | State | Most recent fix (not yet run on GRIT) |
|---|---|---|
| `05`: CA ~100 m mask | Killed three times at the reprojection step | `5d93d1d`: replace `terra::project()` with `sf::gdal_utils("warp")`, working memory capped at 256 MB |
| `01_create_emapr_100m_tifs.R`: CA display composites | First GRIT run killed on the first year | `b5a67ee`: `todisk = TRUE`, restricted to 2005–2010 |
| `06`: % forest per fire, CA | Reached fire 970 of 1,044 before being killed | `0398410`: checkpoints every 10 fires, so re-runs resume |
| `biomass_within_fires.qmd` render | Fails at §6 masking-progression chunks | `2badaad`: read the `_west_` ctrees file and crop to CA before masking |

**Why the 100 m mask is the sticking point.** The 30 m and 90 m masks stay in NLCD's own projection, so
building them is a simple read → reclassify → write in strips. The ~100 m mask has to be reprojected from
EPSG:5070 onto the ctrees lat/lon grid. Reprojection is handled by GDAL's warp engine, which reads the
source in scattered blocks and manages its own working memory. That memory is controlled by neither
`GDAL_CACHEMAX` nor terra's `todisk` option, which is why the first two fixes (`3d1c6ef`, `7d86111`) had no
effect. The current fix calls GDAL warp directly with an explicit `-wm` limit. If it still fails, the next
option is to warp the state in tiles and write each tile as it completes.

Because only the §6 figures need the 100 m mask and the CA-100m eMapR composites, this does not block the
per-fire CSVs. `06`–`08` can finish with just the 30 m and 90 m masks, and the rest of
`biomass_within_fires.qmd` can render with §6 temporarily skipped.

**Scope limits.** Everything is currently CA-only (`STATES_TO_RUN` = CA in `06`–`08`; CA + WY in `05`) and
2005–2010 (`STUDY_YEARS`). Going West-wide waits on the remaining 15 raw eMapR years (19/34 confirmed
complete as of 2026-09-08) and on the steps above running reliably under the memory cap.

**Next run order on GRIT:** `05` → `06` (re-run until CA completes) → `07` → `08` →
`01_create_emapr_100m_tifs.R` → render `biomass_within_fires.qmd`.

---

## 4. Next EDA steps

1. **West-wide within-fire extraction** — after the ctrees West download and West eMapR crops finish
   (`CLAUDE.md` → Current Status), widen `STATES_TO_RUN` / `STUDY_YEARS` and refresh `biomass_within_fires.qmd`
   and `mtbs_assessment_comparison.qmd` for all states and years.
2. **Explain the eMapR-vs-ctrees residual bias** and re-check the border-fire count gap on real multi-state
   output (`NOTES.md` 2026-08-30).
3. **Pre-fire data sufficiency** (Q7) and **spatial clustering** (Q5), feeding the min-pre-years, fire-complex,
   and buffer decisions.
4. **Control-site candidates** — never-burned pixels in the same ecoregion (`data/processed/control_pixels/`).
5. **Tidy `01` §4 and `04` §7** per the table above.

---

## 5. Record: completed MTBS EDA (`01_mtbs_exploration.qmd`)

**Filter steps** (row count printed after each): `incid_type == "Wildfire"` → `year` 2000–2023 → spatial
join to the 11 Western states (`tigris::states()`). **Sanity checks:** required fields present, no duplicate
`event_id`, valid geometries (`st_make_valid()`), dNBR thresholds not all in 0–100 (would mean percentages),
sentinel values −9999/9999 flagged, all 11 target states present. `sf_use_s2(FALSE)` is required — see
`NOTES.md` → Technical hurdles.

**Figures produced** (`figures/`, gitignored — regenerate by rendering the qmd):

| Section | Figure | Content |
|---|---|---|
| 3.1 | `fire_locations_map.png` | All fires, colored by assessment type (EPSG:5070) |
| 3.2 | `fires_per_year.png` | Annual wildfire count with linear trend |
| 3.3 | `asmnt_type_summary.png` | Fire count by assessment type |
| 3.4 | `fire_locations_extended_only.png` | Extended-assessment fires only, with n and % of total |
| 3.5 | `fires_per_year_comparison.png` | All vs. Extended, shared y-axis |
| 3.6 | `fire_size_comparison.png` | log₁₀ size histograms, all vs. Extended, medians marked |
| 3.7 | `severity_dnbr_map.png`, `severity_dnbr_histogram.png`, `severity_rdnbr_histogram.png` | BSP CBI field-plot dNBR/RdNBR (`map_prog == "CBI"`, sentinel 9999 and `|preNBR| ≤ 10` excluded) |

**Notes:** No Extended-only comparison for CBI severity — BSP CBI (field-assessed) and MTBS (satellite-assessed)
are largely different fires, so joining on `event_id` returns near-zero matches. RdNBR = dNBR / √|preNBR / 1000|
(Miller & Thode 2007) and is the leading candidate for the continuous treatment variable (`NOTES.md`).
Methodology and citations: `NOTES.md` → Literature notes.

---

## 6. Figure QA (apply to every plot before saving)

- 300 dpi; check pixel dimensions match the intended inches
- Legend has a descriptive title (not a raw column name), text ≥ 10 pt, no overlap with data
- Colorblind-safe palette (`viridis`, `RColorBrewer "Set2"`, or `MetBrewer`)
- Maps: scale bar, north arrow, state boundaries for reference
- Time axes span the full year range with readable ticks; event-study plots mark year 0 with a dashed line
  and label the x-axis "Years relative to fire"
- Paired panels (all vs. Extended, CA vs. West) share axis scales
