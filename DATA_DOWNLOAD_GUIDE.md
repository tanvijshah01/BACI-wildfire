# Data Download Guide

This guide walks a new user through getting both biomass datasets from their
original source to analysis-ready CSVs. There are three stages, always in
this order:

1. **Download** the raw data (eMapR via anonymous FTP, ctrees via the
   arraylake zarr store).
2. **Crop** it to the 11-state Western study region (eMapR only — ctrees is
   pulled pre-cropped).
3. **Extract** biomass within forested MTBS fire-perimeter pixels (scripts
   `05`–`08`), which is what feeds the `analysis/*.qmd` documents.

This guide is the *process reference* (how to run things). See `CLAUDE.md`
for the directory structure and current project status, and `NOTES.md` for
the reasoning behind decisions and known gotchas. If you just want the
command to run right now, jump to **Part 5: End-to-End Checklist**.

**The two datasets don't work the same way — this trips people up, so read
this table before anything else:**

| | eMapR | ctrees |
|---|---|---|
| **Raw state** | A real, discrete file per year: the full CONUS `composite_YYYY_median.tif` (~27.7 GB BigTIFF), served over plain FTP. This is genuinely "the raw file" — an unmodified copy of the source. | **No discrete raw file exists.** The source is a live, cloud-hosted array (`arraylake`/zarr) — not something distributed as a downloadable file at all. |
| **How it's accessed (the "query")** | A literal file copy: FTP pull of one named path, byte-for-byte, same every time. | A live index/slice query into the remote array via the `arraylake` Python client — request a bounding box + year range, get back whatever pixels are inside it. Nothing is "downloaded" in the traditional sense until a script chooses to save what it read. |
| **Then manipulated by...** | Crop+mask to the 11-state West region (`00_crop_emapr_to_west.R`, §2.3), then per-fire extraction + forest-masking (`07`, Part 4). | The download scripts (`03`/`04`, Part 3) *are* the extraction step — in one pass they save native-resolution regional GeoTIFFs (the closest thing to "raw" that gets kept, see Part 3) alongside already-aggregated NetCDF/CSV outputs. Forest-masked fire extraction (`08`, Part 4) then reads those saved GeoTIFFs. |
| **Raw data retention** | **Kept permanently** (PI request, 2026-08-30) — every downloaded `composite_YYYY_median.tif` stays in `data/raw/emapr_biomass/`, not deleted after cropping. See the storage callout in Part 2. | N/A — the native-resolution GeoTIFFs from Part 3 are the retained artifact. |
| **Analysis-ready output** | `biomass_fire_polygons_emapr_west_<years>_100m_forested.csv` | `biomass_fire_polygons_ctrees_west_forested.csv` |

## Processing at a glance

Raw biomass → fire-level biomass, in the order things run (details in Parts 2–4):

| # | Step | Script | What it does | Output |
|---|---|---|---|---|
| 1 | Fetch raw eMapR | FTP pull + `check_raw_emapr_files.R` | Copy CONUS annual composites (30 m, 1990–2023); validate size and pixel content | `data/raw/emapr_biomass/composite_YYYY_median.tif` |
| 2 | Crop eMapR | `00_crop_emapr_to_west.R` | Crop and mask each year to the union of the 11 Western states | `data/processed/emapr_biomass_west/composite_YYYY_west.tif` (~1 GB/yr) |
| 3 | Download ctrees | `04_download_ctrees_west.py` | Query the arraylake zarr store for the West bbox: (A) native ~100 m annual GeoTIFFs, (B) coarsened 1 km NetCDF, (C) mean AGB per fire × year | `data/processed/ctrees/` |
| 4 | Forest masks | `05_prepare_forest_masks_west.R` | Fetch NLCD 2004 per state; classes 41/42/43 → 1, everything else → 0 (0/1, not 1/NA); write 30 m, 90 m (modal aggregate for eMapR) and ~100 m (ctrees grid) versions | `data/processed/forest_mask/nlcd2004_forestfrac_*_<st>.tif` |
| 5 | % forest per fire | `06_extract_pct_forest_within_fires.R` | Mean of the 30 m mask within each MTBS perimeter | `pct_forest_by_fire_west.csv` |
| 6 | eMapR within fires | `07_extract_emapr_within_fires_new.R` | Per fire polygon: crop → aggregate 30 → ~90 m (mean) → mask non-forest → mean AGB, for each year | `biomass_fire_polygons_emapr_west_<years>_100m_forested.csv` |
| 7 | ctrees within fires | `08_extract_ctrees_within_fires_new.R` | Same per-polygon crop → mask → mean on the ~100 m ctrees TIFs, forest mask projected per polygon | `biomass_fire_polygons_ctrees_west_forested.csv` |
| 8 | Validate | `validate_west_pipeline.R` | Regression-check 05–08 against the retired CA-only baseline | `data/processed/validation/`, `west_pipeline_sanity_check.qmd` |

**Fire selection (steps 5–7):** MTBS `Wildfire` perimeters ≥ 1,000 acres whose ignition year is in
`STUDY_YEARS`; each fire's state comes from its `event_id` prefix (not a raw spatial join), so border fires
are counted once. **Then:** the `analysis/*.qmd` documents read the step 5–7 CSVs; the unit × year panel and
Callaway-Sant'Anna estimation are the planned next stages (`CLAUDE.md`).

**Conventions used throughout:** every extraction is a per-polygon crop → mask → mean (whole-raster masking
OOMs on these files); all scripts are skip-safe and resume per state/year; existing outputs are validated, not
just checked for existence (`NOTES.md` → Technical gotchas).

**Contents**
- [Part 1: One-Time Setup](#part-1-one-time-setup)
- [Part 2: eMapR Biomass](#part-2-emapr-biomass)
- [Part 3: ctrees Biomass](#part-3-ctrees-biomass)
- [Part 4: Forest Mask + Fire-Extraction Pipeline (scripts 05–08)](#part-4-forest-mask--fire-extraction-pipeline-scripts-0508)
- [Part 5: End-to-End Checklist](#part-5-end-to-end-checklist)

---

## Part 1: One-Time Setup

Do these once per machine before downloading anything. **GRIT is the primary
environment** (the raw eMapR archive is ~1 TB and does not fit on a laptop);
laptop-only notes are labelled as such.

### 1.1 GRIT layout

The code repo (`~/BACI-wildfire`, this repo) and the actual data live in two
separate places — raw/processed data lives in `~/BACI-review`, a separate
project-storage repo also shared with other collaborators. `~/BACI-wildfire`
has no `data/` directory of its own; instead `data` there is a symlink:

```bash
ln -s ~/BACI-review/data ~/BACI-wildfire/data
```

This works transparently with every script's existing `here()`/`PROJ_ROOT`-relative
paths — no code changes needed — and is safe because `data/raw/`,
`data/processed/`, `data/final/` are gitignored. If `~/BACI-wildfire/data` is
ever missing, re-create the symlink rather than downloading a second copy.

### 1.2 GRIT Python environment (ctrees downloads)

A venv (not conda) — the safe default that avoids polluting whatever
global/shared Python environment other collaborators use on GRIT:

```bash
cd ~/BACI-wildfire
python3 -m venv .venv
source .venv/bin/activate
pip install arraylake zarr xarray netCDF4 geopandas rasterio
arraylake auth login        # opens a browser prompt; scripts fail with a connection error until done
```

Re-run `source ~/BACI-wildfire/.venv/bin/activate` in any new shell/tmux pane
before running a ctrees script — a fresh pane starts outside the venv. (No
`requirements*.txt` is committed yet; `pip freeze > requirements-ctrees.txt`
would make the env reproducible for the next person.)

`rasterio` is optional (a slower matplotlib-path fallback) for
`03_download_ctrees_ca.py` but a **hard requirement** for
`04_download_ctrees_west.py` — it backs the raw-GeoTIFF write/read path Parts
A–C depend on, and the script exits immediately if it's missing.

### 1.3 Long-running jobs on GRIT: use `tmux`

A dropped browser/SSH connection kills a plain foreground process, and can
leave corrupt half-written outputs (see `NOTES.md` → "`file.exists()` ≠ valid").
Start any multi-hour job inside `tmux`:

```bash
tmux new -s <name>          # start
# Ctrl-b then d             # detach
tmux attach -t <name>       # reattach from any terminal on GRIT
```

### 1.4 `rclone` (eMapR anonymous FTP pull)

`rclone` is only used to pull raw eMapR composites from the public FTP
server. Check it exists (`rclone listremotes`) and create the read-only
anonymous remote if not:

```bash
rclone config create emapr-ftp ftp host=islay.ceoas.oregonstate.edu user=anonymous pass=
```

On the laptop `rclone` is installed at `C:\Users\shaht\bin\rclone.exe` (not on
PATH — call it with the full path). On GRIT, confirm it's available
(`rclone version`) or use the `curl` fallback in §2.4.

### 1.5 Laptop only: prevent sleep during long jobs

Closing the lid interrupts transfers and can corrupt output files mid-write:

```powershell
powercfg /change standby-timeout-ac 0   # disable sleep, before starting
powercfg /change standby-timeout-ac 30  # re-enable, after it finishes
# lid-close can trigger Modern Standby even with the idle timeout off; also set:
powercfg -attributes SUB_BUTTONS LIDACTION -ATTRIB_HIDE
powercfg /setacvalueindex SCHEME_CURRENT SUB_BUTTONS LIDACTION 0   # 0 = Do nothing
powercfg /setactive SCHEME_CURRENT
```

---

## Part 2: eMapR Biomass

| | |
|---|---|
| **Source** | `islay.ceoas.oregonstate.edu` (FTP, anonymous login) |
| **Remote path** | `STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_YYYY_median.tif` |
| **Coverage** | 1990–2023, one file per year, ~27.7 GB each (BigTIFF), CONUS-wide |
| **CRS** | EPSG:5070 (NAD83 / Conus Albers), 30 m, ~96,815 × 153,809 px |
| **Local raw destination** | `data/raw/emapr_biomass/composite_YYYY_median.tif` |

Raw composites are too large to work with directly (loading one in Quarto
causes multi-minute stalls), so they always get cropped before use — see §2.3.

> **Retention policy: raw files are kept, never deleted after cropping.** The
> PI wants the raw archive on hand. The **full 34-year archive is ~950 GB–1 TB**,
> which will not fit on a laptop (that's what caused the repeated disk-space
> failures under the earlier delete-after-crop pattern) — it needs GRIT-scale
> storage. Confirm available quota before assuming "keep everything" is free.

### 2.1 Storage tiers

| Tier | Where it lives | Contents | Size |
|---|---|---|---|
| Raw (kept permanently) | `data/raw/emapr_biomass/` | full CONUS `composite_YYYY_median.tif` | ~30 GB/yr, ~1 TB for 1990–2023 |
| Processed (West-cropped) | `data/processed/emapr_biomass_west/` | `composite_YYYY_west.tif` | ~1 GB/yr |
| Analysis-ready | local / git | fire-polygon extraction CSVs (Part 4) | KB–MB |

### 2.2 Getting a raw year

**Use `rclone` (§1.4), not the manual `curl.exe` loop** — that loop is a
fallback kept in §2.4 only for machines without rclone.

```bash
# pull one year straight from the eMapR FTP into the raw tier
rclone copy emapr-ftp:STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_<yr>_median.tif data/raw/emapr_biomass/ --progress
```

(On the laptop, replace `rclone` with `& "C:\Users\shaht\bin\rclone.exe"`.)
If a year already exists on another machine, copying it over (`scp`/FileZilla)
is faster than re-pulling ~28 GB. Run this inside `tmux` on GRIT.

**Then validate before trusting it.** An interrupted transfer leaves a file
that looks present but is truncated (one real case: a year that held only
67.6% of its pixels while its header still claimed full dimensions):

```bash
Rscript scripts/r/check_raw_emapr_files.R
```

This checks, per year 1990–2023: file exists, size within tolerance of the
known-exact byte count (every complete year is the same size), header opens,
and a small centered pixel block has real values. Console output only; it
never reads the whole file. Delete and re-fetch any year that fails. (As of
2026-09-08, 19/34 years were confirmed complete on GRIT; `CLAUDE.md` tracks
current status.)

### 2.3 Crop to the study region (required before use in Quarto)

Cropping/masking turns a ~28 GB CONUS file into a usable ~1 GB regional file.

```bash
Rscript scripts/r/00_crop_emapr_to_west.R
```

- Input: `data/raw/emapr_biomass/composite_YYYY_median.tif`
- Output: `data/processed/emapr_biomass_west/composite_YYYY_west.tif` (~1 GB, masked to the union of all 11 Western study states)
- **Skip-safe**: only processes years whose raw file is present and whose cropped output doesn't already exist. Years missing from `data/raw/emapr_biomass/` are reported as skipped, not errored — fetch them first (§2.2), then re-run.
- The skip check is existence + validity (`raster_is_valid()`), so a corrupt cropped file is rebuilt automatically.

`scripts/r/00_crop_emapr_to_ca.R` is the retired CA-only predecessor
(outputs to `data/processed/emapr_biomass_ca/`) — it still works but is not
needed for new work.

### 2.4 Fallback: direct FTP with `curl.exe` (no rclone)

Only if rclone truly isn't available. **Windows `ftp.exe` does not support
passive mode** and fails with `Connection closed by remote host`; use
`curl.exe` (built into Windows 10/11), not plain `curl` in PowerShell (which
is aliased to `Invoke-WebRequest` and fails differently).

```powershell
# from a PowerShell prompt opened in data/raw/emapr_biomass/
# one year:
curl.exe --ftp-pasv --user "anonymous:" -O ftp://islay.ceoas.oregonstate.edu/STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_2002_median.tif
# a range (change the loop start to resume):
for ($yr = 1990; $yr -le 2023; $yr++) {
    curl.exe --ftp-pasv --user "anonymous:" -O "ftp://islay.ceoas.oregonstate.edu/STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_${yr}_median.tif"
}
```

On Mac/Linux the standard `ftp` client supports passive mode
(`ftp islay.ceoas.oregonstate.edu`, user `anonymous`, password = your email,
`cd STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/`, `mget composite_<yr>_median.tif`).
Don't `mget *` — that's the full ~1 TB series.

---

## Part 3: ctrees Biomass

| | |
|---|---|
| **Source** | `ucsb-emlab/BACI-wildfires` repo on [arraylake](https://arraylake.com), branch `main` |
| **Access method** | Live zarr reads over the network via the `arraylake` Python client — there is no bulk file download; each script pulls exactly the spatial/temporal subset it needs |
| **Group / variable** | `aboveground_biomass/agb`, shape `(26, 202500, 405000)`, dtype `int16` |
| **Coordinates** | `time` (2000–2025, annual), `x` (lon, ascending), `y` (lat, descending) |
| **Scale / fill** | divide `int16` values by 10 for Mg ha⁻¹; fill value `-9999` |
| **CRS / resolution** | WGS84 / EPSG:4326, ~0.000889° ≈ 100 m |

Unlike eMapR, there's nothing to crop afterward — each script pulls data
already limited to its target bounding box.

**On "raw" ctrees data:** there is no equivalent of eMapR's downloadable
CONUS file — the zarr array *is* the raw dataset and lives on arraylake's
infrastructure. A script "downloading" ctrees means: connect, request a
bounding box + year slice, and choose what to save. Of the three outputs
below, **Part A (the native ~100 m GeoTIFFs) is the closest thing to a
retained raw copy** — full pixel resolution, just spatially clipped to the
study region, no aggregation. The 1 km NetCDF and fire-polygon CSV are
already-aggregated derivatives. So if anyone asks "do we have the raw ctrees
data": yes, as the regional GeoTIFFs, not as a literal copy of arraylake's
storage.

### 3.1 Explore the store (optional — run first if the schema is unfamiliar)

```bash
python scripts/python/02_explore_ctrees_zarr.py
```

Connects to the repo, walks all groups/arrays, and prints coordinate ranges,
variable attributes, and a CA-subset size estimate. Console output only.

### 3.2 Download the California subset (validated baseline)

```bash
python scripts/python/03_download_ctrees_ca.py
```

Pulls the CA bounding box (`lon -124.5 to -114.1`, `lat 32.5 to 42.0`):

| Part | Output | Purpose |
|---|---|---|
| A | `data/processed/ctrees/ctrees_YYYY_ca_100m.tif` | One native-resolution (~100 m) GeoTIFF per year, for comparison against eMapR |
| B | `data/processed/ctrees/ctrees_biomass_ca_1km.nc` | Coarsened (~1 km) CA raster, 26 years, for R mapping |
| C | `data/processed/ctrees/biomass_fire_polygons_ctrees.csv` | Long panel: `event_id` × `year` × mean AGB within each MTBS CA fire polygon |

Skip-safe per part — delete a specific output to force re-extraction of just
that part. Part C rasterizes each fire polygon's mask once and reuses it
across all 26 years (see `NOTES.md` 2026-05-13). Sanity checks run
automatically at the end.

**This CA output is the validated baseline** (`biomass_within_fires.qmd` §7) —
don't modify it; the West-wide pipeline below is cross-checked against it.
Part letters match §3.3 (A = raw TIFF, B = 1 km NetCDF, C = fire CSV), but
this script's execution order is untouched: B → C → A.

### 3.3 Download the full 11-state Western subset

```bash
python scripts/python/04_download_ctrees_west.py
```

Same as §3.2 over the union bbox of all 11 Western study states
(`lon -124.8 to -102.0`, `lat 31.3 to 49.0`, ~4.1x the CA pixel area, ~6,800
fires). Unlike §3.2, this script runs in A → B → C order:

| Part | Output | Purpose |
|---|---|---|
| A | `data/processed/ctrees/ctrees_YYYY_west_100m.tif` | Raw download: one native-resolution (~100 m) GeoTIFF per year, straight from arraylake — runs first |
| B | `data/processed/ctrees/ctrees_biomass_west_1km.nc` | Coarsened (~1 km) West raster, 26 years, for R mapping — built from Part A's local TIFFs |
| C | `data/processed/ctrees/biomass_fire_polygons_ctrees_west.csv` | Long panel: `event_id` × `year` × mean AGB within each Western fire polygon — also reads Part A's local TIFFs |

Raw first means "raw ctrees data" is on disk as its own step before any
processing, and B/C read plain local files instead of re-querying arraylake.
Because B/C read Part A's TIFFs, **rasterio is a hard requirement**. Part A's
per-year TIFs are **shared across all 11 states** — each state's extraction
(`08`, Part 4) crops its own slice from the same file.

**Run it in `tmux`** (§1.3), from a real terminal, not a notebook. It's a
multi-hour job (~4x the reads, ~6.4x the fire polygons vs. CA); expect the 26
compressed GeoTIFFs to total ~8–20 GB, so check free space first (`df -h ~`).
On a laptop, disable sleep instead (§1.5).

```bash
tmux new -s ctrees_west
cd ~/BACI-wildfire && source .venv/bin/activate
python scripts/python/04_download_ctrees_west.py 2>&1 | tee -a data/processed/ctrees/04_download_ctrees_west_log.txt
```

**Resuming / validity.** Part B checkpoints each coarsened year to
`data/processed/ctrees/_west_1km_scratch/` and only assembles the NetCDF once
all 26 years exist; Part C checkpoints per year to `_west_fireagb_scratch/`.
On resume, Part B re-validates any existing `ctrees_biomass_west_1km.nc`
(`netcdf_is_valid()`) and deletes/rebuilds a corrupt one. Part A's TIFs are
still existence-checked only — after any interrupted run, compare the
last-written year's file size to the others (~730 MB) and delete it if it's
short. The CA outputs from §3.2 are never touched by this script.

**Memory on GRIT.** Interactive sessions run as Slurm jobs with a cgroup
memory cap, and this script was OOM-killed three times before being
restructured to *never materialize a full ~2 GB year array*: Part A writes in
row-strips (`STRIP_ROWS`), Part B coarsens from `factor`-row TIFF windows
(`coarsen_year_from_tif()`) and writes the NetCDF one year at a time with
plain `netCDF4`, and Part C reads per-fire `rasterio.Window`s with bit-packed
mask caching. The script logs peak memory (`VmHWM`) after each part — check
`04_download_ctrees_west_log.txt` first if a run dies. If OOM kills persist,
inspect the cap (`cat /sys/fs/cgroup/memory.max`) and ask whoever administers
GRIT for a larger allocation. Incident history: `NOTES.md` 2026-09-20.

**After the first successful full run**, cross-validate: the West CSV's CA
rows should match the `03_download_ctrees_ca.py` baseline (expect correlation
≈ 1.000, as the ctrees side of the `05`–`08` validation did).

---

## Part 4: Forest Mask + Fire-Extraction Pipeline (scripts 05–08)

Once eMapR is cropped (§2.3) and ctrees is downloaded (§3.3) for a state, four
R scripts turn those rasters into the CSVs the `analysis/*.qmd` documents
read — building a forest mask, then extracting biomass within fire
perimeters. Run them in order, from the project root, **outside Quarto**
(Quarto buffers chunk output, which makes terra's C++ threading look frozen).

```r
# 1 — per-state NLCD 2004 forest masks (0/1), three resolutions
Rscript scripts/r/05_prepare_forest_masks_west.R

# 2 — % forest cover within each MTBS fire perimeter
Rscript scripts/r/06_extract_pct_forest_within_fires.R

# 3 — eMapR AGB within fire perimeters, forested pixels only
Rscript scripts/r/07_extract_emapr_within_fires_new.R

# 4 — ctrees AGB within fire perimeters, forested pixels only
Rscript scripts/r/08_extract_ctrees_within_fires_new.R
```

| Script | Depends on | Output |
|---|---|---|
| `05` | NLCD 2004 (auto-downloaded from MRLC's WCS endpoint via `fetch_nlcd_landcover()`, no login needed); a ctrees `_100m.tif` template for the ~100 m mask variant | `data/processed/forest_mask/nlcd2004_forestfrac_{30m,90m,100m}_<state>.tif` |
| `06` | `05`'s 30 m masks | `data/processed/forest_mask/pct_forest_by_fire_west.csv` |
| `07` | `05`'s 90 m masks; §2.3's West-cropped eMapR TIFs | `data/processed/emapr_biomass_west/biomass_fire_polygons_emapr_west_<years>_100m_forested.csv` |
| `08` | `05`'s 30 m masks; §3.3's West-wide ctrees `_100m.tif` files | `data/processed/ctrees/biomass_fire_polygons_ctrees_west_forested.csv` |

**Each script is skip-safe**, resuming per state (and `07`/`08` per
state × year) — re-running only builds what's missing. Delete the relevant
output file/rows to force re-extraction.

**Adding a new state:** each script has a `STATES_TO_RUN` list near the top
(e.g. `STATES_TO_RUN <- c("CA", "WY")`) — add the state's two-letter code and
re-run `05`→`08` in order. `05` needs a ctrees template TIF from §3.3 for the
~100 m mask; `07`/`08` need `05`'s masks for that state first.

**Validation harness:** `Rscript scripts/r/validate_west_pipeline.R` (+
`analysis/west_pipeline_sanity_check.qmd`) re-checks the pipeline against the
retired CA baseline; resumable per (state, year).

---

## Part 5: End-to-End Checklist

Adding a brand-new state (or year) to the study, from nothing to
analysis-ready CSVs:

1. **One-time setup** (Part 1) done? GRIT symlink + venv, `arraylake auth login`, `emapr-ftp` remote.
2. **eMapR:** for each study year, get the raw composite (§2.2), validate it (`Rscript scripts/r/check_raw_emapr_files.R`), then crop (`Rscript scripts/r/00_crop_emapr_to_west.R`, §2.3).
3. **ctrees:** run `03_download_ctrees_ca.py` once ever (§3.2, baseline); run `04_download_ctrees_west.py` once ever for the West-wide TIFs (§3.3) — both are shared across all states, not re-run per state.
4. **Forest mask + extraction:** add the new state's code to `STATES_TO_RUN` in scripts `05`, `06`, `07`, `08`, then run them in that order (Part 4).
5. Confirm the new state's rows appear in the output CSVs in the Part 4 table, then point the relevant `analysis/*.qmd` at them.

For current progress against this checklist (which states/years are done,
what's still blocked), see "Current Status" in `CLAUDE.md` — that's the
living status log; this is the process reference.
