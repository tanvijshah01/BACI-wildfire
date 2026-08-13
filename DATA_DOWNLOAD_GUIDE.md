# Data Download Guide

This guide walks a new user through getting both biomass datasets from their
original source to analysis-ready CSVs. There are three stages, always in
this order:

1. **Download** the raw data (eMapR via FTP/Nextcloud, ctrees via the
   arraylake zarr store).
2. **Crop** it to the 11-state Western study region (eMapR only — ctrees is
   pulled pre-cropped).
3. **Extract** biomass within forested MTBS fire-perimeter pixels (scripts
   `05`–`08`), which is what feeds the `analysis/*.qmd` documents.

See `CLAUDE.md` for the full directory structure, pipeline diagram, and
current project status. If you just want the command to run right now, jump
to **Part 5: End-to-End Checklist**.

**Contents**
- [Part 1: One-Time Setup](#part-1-one-time-setup)
- [Part 2: eMapR Biomass](#part-2-emapr-biomass)
- [Part 3: ctrees Biomass](#part-3-ctrees-biomass)
- [Part 4: Forest Mask + Fire-Extraction Pipeline (scripts 05–08)](#part-4-forest-mask--fire-extraction-pipeline-scripts-0508)
- [Part 5: End-to-End Checklist](#part-5-end-to-end-checklist)

---

## Part 1: One-Time Setup

Do these once per machine before downloading anything.

| Tool | Needed for | Setup |
|---|---|---|
| `rclone` | eMapR raw downloads (Nextcloud + FTP) | Installed at `C:\Users\shaht\bin\rclone.exe` (not on PATH — call with the full path, or add the folder to PATH yourself) |
| `arraylake` Python client | ctrees downloads | `arraylake auth login` (opens a browser prompt; scripts fail with a connection error until this is done) |
| Python packages | ctrees downloads | `arraylake`, `zarr`, `xarray`, `netCDF4`, `geopandas`; `rasterio` optional but recommended (falls back to a slower rasterizer if missing) |

**rclone remotes** — check these exist before using rclone (`rclone listremotes`):

```powershell
# Anonymous eMapR FTP remote (read-only, no credentials needed)
& "C:\Users\shaht\bin\rclone.exe" config create emapr-ftp ftp host=islay.ceoas.oregonstate.edu user=anonymous pass=

# Nextcloud remote — needs a WebDAV URL + Nextcloud app password
# (Nextcloud -> Settings -> Security -> "Create new app password")
& "C:\Users\shaht\bin\rclone.exe" config create nextcloud webdav url=<WEBDAV_URL> vendor=nextcloud user=<USER> pass=<APP_PASSWORD>
```

(`rclone config create` accepts a plaintext `pass=` and obscures it in the stored config — no separate `rclone obscure` step needed.)

**Prevent sleep during any multi-hour download** — closing the laptop lid interrupts transfers and can corrupt output files mid-write:

```powershell
powercfg /change standby-timeout-ac 0   # disable sleep, before starting
powercfg /change standby-timeout-ac 30  # re-enable, after it finishes
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
causes multi-minute stalls) and too large to keep permanently on the laptop
(~1 TB for the full archive, which has already caused repeated disk-space
failures). The pattern is: **get a year's raw file onto the laptop just long
enough to crop it, then discard the raw file** — the small cropped output is
what actually gets used.

### 2.1 Storage tiers

| Tier | Where it lives | Contents | Size |
|---|---|---|---|
| Raw (archival) | Nextcloud: `BACI/raw/emapr_biomass/` | full CONUS `composite_YYYY_median.tif` | ~30 GB/yr |
| Processed (West-cropped) | `data/processed/emapr_biomass_west/`, also mirrored to Nextcloud `BACI/processed/emapr_biomass_west/` | `composite_YYYY_west.tif` | ~1 GB/yr |
| Analysis-ready | local / git | fire-polygon extraction CSVs (Part 4) | KB–MB |

### 2.2 Downloading a year (recommended: rclone)

Pick the branch that matches where the year currently lives. **Always use
`rclone`, not the manual `curl.exe --ftp-pasv` loop** — that loop is a
retired fallback kept in §2.4 only for machines without rclone available.

**A — Year's raw file is already local** (`data/raw/emapr_biomass/`): skip straight to §2.3 (crop).

**B — Year is archived on Nextcloud but not local:** pull it down, crop it, push the small result back up, then delete the scratch raw file so it doesn't accumulate on disk:

```powershell
& "C:\Users\shaht\bin\rclone.exe" copy nextcloud:BACI/raw/emapr_biomass/composite_<yr>_median.tif data/raw/emapr_biomass/ --progress
Rscript scripts/r/00_crop_emapr_to_west.R
& "C:\Users\shaht\bin\rclone.exe" copy data/processed/emapr_biomass_west/composite_<yr>_west.tif nextcloud:BACI/processed/emapr_biomass_west/
Remove-Item data/raw/emapr_biomass/composite_<yr>_median.tif
```

**C — Year isn't on Nextcloud or local yet (first-ever fetch):** stream it FTP → Nextcloud directly — bytes flow server-to-server, so it never lands as a full file on the laptop — then follow branch B:

```powershell
& "C:\Users\shaht\bin\rclone.exe" copy emapr-ftp:STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_<yr>_median.tif nextcloud:BACI/raw/emapr_biomass/ --progress
# then run branch B for the same year
```

**Archiving years you already have locally** (e.g. after using the legacy FTP method below), without deleting them:

```powershell
& "C:\Users\shaht\bin\rclone.exe" copy data/raw/emapr_biomass/ nextcloud:BACI/raw/emapr_biomass/ --progress
```

Deleting the local raw archive after an upload is confirmed is a manual,
user-triggered step — not automated by any script here.

### 2.3 Crop to the study region (required before use in Quarto)

Cropping/masking must happen before eMapR data is used anywhere — this is
what turns a ~28 GB CONUS file into a usable ~1 GB regional file.

```r
Rscript scripts/r/00_crop_emapr_to_west.R
```

- Input: `data/raw/emapr_biomass/composite_YYYY_median.tif`
- Output: `data/processed/emapr_biomass_west/composite_YYYY_west.tif` (~1 GB, masked to the union of all 11 Western study states)
- **Skip-safe**: only processes years whose raw file is present locally and whose cropped output doesn't already exist. Years missing from `data/raw/emapr_biomass/` are reported as skipped, not errored — fetch them first (§2.2), then re-run.

`scripts/r/00_crop_emapr_to_ca.R` is the retired CA-only predecessor
(outputs to `data/processed/emapr_biomass_ca/`) — it still works but is not
needed for any new work; the West-wide script above supersedes it.

### 2.4 Legacy fallback: direct FTP to the laptop (no rclone)

Only use this if rclone truly isn't available. It downloads the full raw
archive straight to `data/raw/emapr_biomass/` with no Nextcloud step, which
is exactly the disk-space-failure pattern §2 is designed to avoid — prefer
§2.2 whenever possible.

**Windows `ftp.exe` does not support passive mode** and will fail with
`Connection closed by remote host`. Use `curl.exe` instead (built into
Windows 10/11) — not plain `curl` in PowerShell, which is aliased to
`Invoke-WebRequest` and will fail silently in a different way.

```powershell
# from a PowerShell prompt opened in data/raw/emapr_biomass/
for ($yr = 1990; $yr -le 2023; $yr++) {
    $url = "ftp://islay.ceoas.oregonstate.edu/STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_${yr}_median.tif"
    Write-Host "Downloading $yr..."
    curl.exe --ftp-pasv --user "anonymous:" -O $url
}
```

To resume from a specific year, change the loop start (e.g. `$yr = 1992`).
Monitor progress from a second window:

```powershell
while ($true) {
    $files = Get-ChildItem . -Filter "*.tif"
    $tmp   = Get-ChildItem . -Filter "*.tmp"
    $totalGB = [math]::Round(($files + $tmp | Measure-Object Length -Sum).Sum / 1GB, 2)
    Write-Host "$(Get-Date -Format 'HH:mm:ss')  $($files.Count) tif complete | downloading: $($tmp.Name) | total on disk: $totalGB GB"
    ($files + $tmp) | ForEach-Object { $_.Refresh() }
    Start-Sleep 15
}
```

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

Unlike eMapR, there's nothing to crop afterward — each script below pulls
data already limited to its target bounding box.

### 3.1 Explore the store (optional — run first if the schema is unfamiliar)

```bash
python scripts/python/02_explore_ctrees_zarr.py
```

Connects to the repo, walks all groups/arrays, and prints coordinate
ranges, variable attributes, and a CA-subset size estimate. Console output
only — writes no files.

### 3.2 Download the California subset (validated baseline)

```bash
python scripts/python/03_download_ctrees_ca.py
```

Pulls the CA bounding box (`lon -124.5 to -114.1`, `lat 32.5 to 42.0`) and writes three outputs:

| Part | Output | Purpose |
|---|---|---|
| A | `data/processed/ctrees/ctrees_biomass_ca_1km.nc` | Coarsened (~1 km) CA raster, 26 years, for R mapping |
| B | `data/processed/ctrees/biomass_fire_polygons_ctrees.csv` | Long panel: `event_id` × `year` × mean AGB within each MTBS CA fire polygon |
| C | `data/processed/ctrees/ctrees_YYYY_ca_100m.tif` | One native-resolution (~100 m) GeoTIFF per year, for comparison against eMapR |

Skip-safe per part — delete a specific output file to force re-extraction
of just that part. Part B rasterizes each fire polygon's mask once
(`rasterio.features.rasterize()`) and reuses it across all 26 years —
10–50x faster than a shapely point-in-polygon approach. Sanity checks
(percent-valid-pixel counts) run automatically at the end.

**This CA output is the validated baseline** (`biomass_within_fires.qmd`
§7) — don't modify it; it's what the West-wide pipeline below is
cross-checked against.

### 3.3 Download the full 11-state Western subset

```bash
python scripts/python/04_download_ctrees_west.py
```

Same as §3.2 but over the union bbox of all 11 Western study states
(`lon -124.8 to -102.0`, `lat 31.3 to 49.0`, ~4.1x the CA pixel area,
~6,800 fires):

| Part | Output | Purpose |
|---|---|---|
| A | `data/processed/ctrees/ctrees_biomass_west_1km.nc` | Coarsened (~1 km) West raster, 26 years, for R mapping |
| B | `data/processed/ctrees/biomass_fire_polygons_ctrees_west.csv` | Long panel: `event_id` × `year` × mean AGB within each Western fire polygon |
| C | `data/processed/ctrees/ctrees_YYYY_west_100m.tif` | One native-resolution (~100 m) GeoTIFF per year, West-bbox extent |

Part C's per-year West TIFs are meant to be **shared across all 11
states** — each state's extraction (script `08`, Part 4 below) crops its
own slice from the same file rather than needing a separate download per
state.

**This is a multi-hour job** (4x the reads, ~6.4x the fire polygons vs.
CA). Before starting:
- Disable sleep (`powercfg /change standby-timeout-ac 0`, §1).
- Run from a real terminal, not a notebook.
- Expect Part C's 26 compressed GeoTIFFs to total ~8–20 GB on disk.

Part A checkpoints each coarsened year to
`data/processed/ctrees/_west_1km_scratch/` as it completes and only
assembles the final NetCDF once all 26 years are present — safe to resume
if interrupted. The CA-only outputs from §3.2 are left untouched by this
script.

---

## Part 4: Forest Mask + Fire-Extraction Pipeline (scripts 05–08)

Once eMapR is cropped (§2.3) and ctrees is downloaded (§3.2/3.3) for a
state, four R scripts turn those rasters into the CSVs the `analysis/*.qmd`
documents actually read — building a forest mask, then extracting biomass
within fire perimeters. Run them in order, from the project root, **outside
Quarto** (see the note below on why).

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

| Script | Downloads it depends on | Output |
|---|---|---|
| `05` | NLCD 2004 (auto-downloaded via `FedData::get_nlcd()`, no login needed); a ctrees `_100m.tif` template for the ~100 m mask variant | `data/processed/forest_mask/nlcd2004_forestfrac_{30m,90m,100m}_<state>.tif` |
| `06` | `05`'s 30 m masks | `data/processed/forest_mask/pct_forest_by_fire_west.csv` |
| `07` | `05`'s 90 m masks; §2.3's West-cropped eMapR TIFs | `data/processed/emapr_biomass_west/biomass_fire_polygons_emapr_west_<years>_100m_forested.csv` |
| `08` | `05`'s 30 m masks; §3.3's West-wide ctrees `_100m.tif` files | `data/processed/ctrees/biomass_fire_polygons_ctrees_west_forested.csv` |

**Each script is skip-safe**, resuming per state (and `07`/`08` per
state × year) — re-running only builds what's missing. Delete the relevant
output file/rows to force re-extraction.

**Adding a new state:** each script has a `STATES_TO_RUN` list near the
top (e.g. `STATES_TO_RUN <- c("CA", "WY")`) — add the state's two-letter
code and re-run scripts `05`→`08` in order. `05` requires that state's
NLCD download (automatic) and, for the ~100 m mask, a ctrees template TIF
from §3.3. `07`/`08` require `05`'s masks to exist first for that state.

**Windows/Quarto note:** always run these via `Rscript` or the R console,
never inside a Quarto chunk — Quarto buffers chunk output until the chunk
finishes, which makes terra's C++ threading look frozen even though it's
running.

---

## Part 5: End-to-End Checklist

Adding a brand-new state (or year) to the study, from nothing to
analysis-ready CSVs:

1. **One-time setup** (Part 1) done? `arraylake auth login`, rclone remotes configured.
2. **eMapR:** for each study year, get the raw composite local (Part 2.2), then crop it (`Rscript scripts/r/00_crop_emapr_to_west.R`, §2.3).
3. **ctrees:** run `03_download_ctrees_ca.py` once ever (§3.2, baseline); run `04_download_ctrees_west.py` once ever for the West-wide TIFs (§3.3) — both are shared across all states, not re-run per state.
4. **Forest mask + extraction:** add the new state's code to `STATES_TO_RUN` in scripts `05`, `06`, `07`, `08`, then run them in that order (Part 4).
5. Confirm the new state's rows appear in the output CSVs listed in the Part 4 table, then point the relevant `analysis/*.qmd` at them.

For current progress against this checklist (which states/years are done,
what's still blocked), see the "Current Status" section of `CLAUDE.md`
rather than this guide — that's a living status log, this is the
process reference.
