# Data Download Guide

Step-by-step instructions for pulling the raw biomass datasets used in this project. See `CLAUDE.md` for the overall data pipeline and directory structure.

---

## eMapR Biomass (islay.ceoas.oregonstate.edu)

**Server:** `islay.ceoas.oregonstate.edu` (FTP, anonymous login)
**Remote path:** `STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/`
**Local destination:** `data/raw/emapr_biomass/` (files confirmed here, not in `data/raw/` directly)
**Files:** `composite_YYYY_median.tif` — one per year, 1990–2023, ~27.7 GB each (BigTIFF format)
**CRS:** EPSG:5070 (NAD83 / Conus Albers), 30 m resolution, ~96,815 × 153,809 pixels CONUS-wide

**Important:** Windows `ftp.exe` does not support passive mode and will get `Connection closed by remote host` when running `mget *`. Use `curl.exe` instead (built into Windows 10/11). Do not use `curl` in PowerShell — it is an alias for `Invoke-WebRequest` and will fail; always use `curl.exe` explicitly.

Download all years from a PowerShell prompt opened in `data/raw/emapr_biomass/`:

```powershell
for ($yr = 1990; $yr -le 2023; $yr++) {
    $url = "ftp://islay.ceoas.oregonstate.edu/STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_${yr}_median.tif"
    Write-Host "Downloading $yr..."
    curl.exe --ftp-pasv --user "anonymous:" -O $url
}
```

To resume from a specific year (e.g., if 1990–1991 already downloaded), change the loop start:

```powershell
for ($yr = 1992; $yr -le 2023; $yr++) { ... }
```

Monitor download progress in a second PowerShell window:

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

Prevent the laptop from sleeping during download (required — closing the lid interrupts the transfer):

```powershell
powercfg /change standby-timeout-ac 0   # disable sleep on AC power
```

Re-enable after download completes:

```powershell
powercfg /change standby-timeout-ac 30
```

### Post-download preprocessing (required before rendering QMD)

The CONUS-wide composites (~27.7 GB each) must be cropped to the study region before use — loading them directly in Quarto causes multi-minute stalls.

- **California EDA:** Run `scripts/r/00_crop_emapr_to_ca.R` once after downloading. Outputs land in `data/processed/emapr_biomass_ca/composite_YYYY_ca.tif` (~100–300 MB each). `analysis/03_emapr_biomass_exploration.qmd` reads from these, not the raw CONUS files.
- **Western US analysis:** Run `scripts/r/00_crop_emapr_to_west.R`, which crops/masks to the union of all 11 Western study states (AZ, CA, CO, ID, MT, NV, NM, OR, UT, WA, WY). Outputs land in `data/processed/emapr_biomass_west/composite_YYYY_west.tif` (~1 GB per year). Only processes years whose raw file is already present in `data/raw/emapr_biomass/` — see the Nextcloud/rclone section below for fetching years that aren't local yet.

Both preprocessing scripts are skip-safe — they check for existing output files and only process missing years.

---

## Nextcloud Archive (rclone)

Raw eMapR CONUS composites are large (~30 GB/year, ~1 TB for the full
1990–2023 archive) and repeatedly caused local disk-space failures. The
long-term home for the raw archive is **Nextcloud** (accessed via its WebDAV
endpoint), fetched with **rclone** — not the laptop. Read pattern is
**crop-once, then work locally**: the raw Nextcloud tier is only touched by
the (occasional) crop step, never queried directly at analysis time.

### Storage tiers

| Tier | Nextcloud path | Contents | Size |
|---|---|---|---|
| Raw (archival) | `BACI/raw/emapr_biomass/` | full CONUS `composite_YYYY_median.tif` | ~30 GB/yr, ~1 TB total |
| Processed (also cached locally) | `BACI/processed/emapr_biomass_west/` | West-cropped `composite_YYYY_west.tif` | ~1 GB/yr |
| Analysis-ready | stays in git / local | fire-polygon extraction CSVs | KB–MB |

### rclone setup

`rclone` is installed at `C:\Users\shaht\bin\rclone.exe` (not on PATH —
call it with the full path, or add that folder to PATH yourself).

The anonymous eMapR FTP remote is already configured:
```powershell
& "C:\Users\shaht\bin\rclone.exe" config create emapr-ftp ftp host=islay.ceoas.oregonstate.edu user=anonymous pass=
```

The Nextcloud remote still needs to be created once WebDAV URL + an app
password are available (Nextcloud → Settings → Security → "Create new app
password"):
```powershell
& "C:\Users\shaht\bin\rclone.exe" config create nextcloud webdav url=<WEBDAV_URL> vendor=nextcloud user=<USER> pass=<APP_PASSWORD>
```
(`rclone config create` accepts a plaintext `pass=` and obscures it in the
stored config — no separate `rclone obscure` step needed.)

### A. Archive already-local raw years to Nextcloud

Doesn't touch or delete local files — just uploads what's already in
`data/raw/emapr_biomass/`:
```powershell
& "C:\Users\shaht\bin\rclone.exe" copy data/raw/emapr_biomass/ nextcloud:BACI/raw/emapr_biomass/ --progress
```

### B. Stream still-missing years directly FTP → Nextcloud

Never lands the full file permanently on the laptop — bytes flow
server-to-server through rclone:
```powershell
& "C:\Users\shaht\bin\rclone.exe" copy emapr-ftp:STEM_CONUS_BIOMASS/biomassfiaald-v1990-2023-1/composite_<yr>_median.tif nextcloud:BACI/raw/emapr_biomass/ --progress
```

### C. Crop a year that only exists on Nextcloud (not locally)

Pull it to a scratch path, run the crop script, upload the small result,
delete the scratch raw file — local disk only ever holds one transient
~30 GB file at a time, not an accumulating archive:
```powershell
& "C:\Users\shaht\bin\rclone.exe" copy nextcloud:BACI/raw/emapr_biomass/composite_<yr>_median.tif data/raw/emapr_biomass/ --progress
Rscript scripts/r/00_crop_emapr_to_west.R
& "C:\Users\shaht\bin\rclone.exe" copy data/processed/emapr_biomass_west/composite_<yr>_west.tif nextcloud:BACI/processed/emapr_biomass_west/
Remove-Item data/raw/emapr_biomass/composite_<yr>_median.tif
```

### D. Upload the West-cropped processed tier

```powershell
& "C:\Users\shaht\bin\rclone.exe" copy data/processed/emapr_biomass_west/ nextcloud:BACI/processed/emapr_biomass_west/ --progress
```

Deleting the local raw archive after upload is confirmed is a manual,
user-triggered step — not automated here.

---

## Ctrees Biomass (arraylake zarr store)

**Source:** `ucsb-emlab/BACI-wildfires` repo on [arraylake](https://arraylake.com), branch `main`
**Access method:** Direct zarr read over the network via the `arraylake` Python client — there is no bulk file download; data is pulled on demand for whatever spatial/temporal subset you request.
**Group / variable:** `aboveground_biomass/agb`, shape `(26, 202500, 405000)`, dtype `int16`
**Coordinates:** `time` (2000–2025, annual), `x` (lon, ascending), `y` (lat, descending, 90 to −90)
**Scale factor:** divide stored `int16` values by 10 to get Mg ha⁻¹
**Fill value:** `-9999`
**CRS:** WGS84 / EPSG:4326
**Resolution:** ~0.000889° ≈ 100 m

### Prerequisite — arraylake authentication

Before any script can connect to the zarr store, authenticate once per machine:

```bash
arraylake auth login
```

This opens a browser login prompt. Without it, `Client()` calls in the scripts below fail with a connection error (see `scripts/python/03_download_ctrees_ca.py`, which prints `Run arraylake auth login and retry` on failure).

Also ensure the Python dependencies are installed: `arraylake`, `zarr`, `xarray`, `netCDF4`, `geopandas`, `rasterio` (optional but recommended — falls back to a slower `matplotlib.path` rasterizer if missing).

### Step 1 — explore the store (optional, run first if the schema is unfamiliar)

```bash
python scripts/python/02_explore_ctrees_zarr.py
```

Connects to the repo, walks all groups/arrays, and prints coordinate ranges, variable attributes, and a CA-subset size estimate. Console output only — no files written.

### Step 2 — pull the California data

```bash
python scripts/python/03_download_ctrees_ca.py
```

This is the actual "download" step — it reads the CA bounding box (`lon -124.5 to -114.1`, `lat 32.5 to 42.0`) out of the zarr store and writes three local outputs:

| Part | Output | Purpose |
|---|---|---|
| A | `data/processed/ctrees/ctrees_biomass_ca_1km.nc` | Coarsened (~1 km) CA raster, 26 years, for R mapping |
| B | `data/processed/ctrees/biomass_fire_polygons_ctrees.csv` | Long panel: `event_id` × `year` × mean AGB within each MTBS CA fire polygon |
| C | `data/processed/ctrees/ctrees_YYYY_ca_100m.tif` | One native-resolution (~100 m) GeoTIFF per year, for comparison against eMapR |

The script is skip-safe: each part checks for its existing output and skips ahead if already present. Delete the relevant output file to force re-extraction of that part.

Runtime notes:
- Part B precomputes a rasterized mask per fire polygon once (via `rasterio.features.rasterize()`), then reuses it across all 26 years with pure numpy indexing — this replaced an earlier shapely point-in-polygon approach that took ~130 minutes; the rasterio approach is 10–50x faster. See `NOTES.md` (2026-05-13 entry) for the full comparison.
- Sanity checks run automatically at the end of the script (percent-valid pixel checks on all three outputs).

### Step 3 — pull the full Western US data (11 states)

```bash
python scripts/python/04_download_ctrees_west.py
```

Generalizes step 2 from the CA-only bounding box (`lon -124.5 to -114.1`, `lat 32.5 to 42.0`) to the union bbox of all 11 Western study states (`lon -124.8 to -102.0`, `lat 31.3 to 49.0`, ~4.1x the pixel area). Same three outputs, `_west_` instead of `_ca_`:

| Part | Output | Purpose |
|---|---|---|
| A | `data/processed/ctrees/ctrees_biomass_west_1km.nc` | Coarsened (~1 km) West raster, 26 years, for R mapping |
| B | `data/processed/ctrees/biomass_fire_polygons_ctrees_west.csv` | Long panel: `event_id` x `year` x mean AGB within each MTBS Western fire polygon (~6,800 fires) |
| C | `data/processed/ctrees/ctrees_YYYY_west_100m.tif` | One native-resolution (~100 m) GeoTIFF per year, West-bbox extent |

The CA-only outputs from step 2 are left untouched — they remain the validated baseline (`biomass_within_fires.qmd` §7). Part C's West-bbox TIFs are meant to be shared across all 11 states: once `scripts/r/07_extract_emapr_within_fires_new.R` and `08_extract_ctrees_within_fires_new.R` are generalized from `STATE_FIPS <- "CA"` to loop over `WESTERN_STATES`, each state crops its own slice from the same per-year West TIF rather than needing a separate download per state.

Runtime/resource notes:
- Per-year raw read is ~2 GB (float32) — fine in memory, but this is a multi-hour job overall (4x the reads, ~6.4x the fire polygons vs. the CA-only run). **Disable sleep before starting** (`powercfg /change standby-timeout-ac 0`, same as the eMapR download) and run from a real terminal, not a notebook.
- Part A checkpoints each coarsened year to `data/processed/ctrees/_west_1km_scratch/` as it's computed and only assembles the final NetCDF once all 26 years are present — resume-safe if the run is interrupted (the CA-only script's Part A does not do this, since a CA-scale run was short enough not to need it).
- Part C's 26 compressed GeoTIFFs total roughly 8–20 GB on disk.

