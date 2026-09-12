# =============================================================================
# 04_download_ctrees_west.py
# Download Ctrees aboveground biomass data for the full 11-state Western study
# region from the ucsb-emlab/BACI-wildfires arraylake zarr store.
#
# Generalizes 03_download_ctrees_ca.py from a California-only bounding box to
# the union bbox of all 11 Western study states (AZ, CA, CO, ID, MT, NV, NM,
# OR, UT, WA, WY) — same three-part structure, same extraction method. The
# CA-only outputs from 03 are left untouched as the validated baseline
# (see biomass_within_fires.qmd §7); this script writes a separate "_west_"
# output tier alongside them, mirroring how 00_crop_emapr_to_west.R produces
# a "_west_" tier alongside the retired CA-only eMapR crop.
#
# NOTE ON PART LETTERS vs. 03_download_ctrees_ca.py: this script's parts run
# in alphabetical/execution order (A = raw download, B = coarsened NetCDF,
# C = fire-polygon CSV) — 03's A/B/C instead match output type regardless of
# order (A = 1km NetCDF, B = fire CSV, C = raw TIFF) and still run A -> B -> C
# in that original order. The letter for a given output differs between the
# two scripts; DATA_DOWNLOAD_GUIDE.md documents each script's own mapping.
#
# Part A's per-year GeoTIFFs (ctrees_YYYY_west_100m.tif) are meant to be
# shared across all 11 states — scripts/r/07 and 08 will eventually loop
# STATE_FIPS over WESTERN_STATES and crop this same West-wide TIF per state
# per year, rather than needing one ctrees TIF per state.
#
# SCRIPT OUTLINE
# 1.  Setup — connect to arraylake, define parameters
# 2.  Open zarr store and resolve West bounding-box indices
# 3.  Part A — Native-resolution (~100 m) annual GeoTIFFs (raw download)
#       Runs FIRST: write one GeoTIFF per year at native 100 m resolution,
#       West-bbox extent, straight from arraylake. This is the retained
#       "raw" copy (see DATA_DOWNLOAD_GUIDE.md Part 3) and the only step in
#       this script that talks to arraylake for actual pixel data — running
#       it first means Part B/C below process a plain local raster instead
#       of a live remote read, isolating arraylake/icechunk overhead from
#       the heavier in-memory coarsening math.
# 4.  Part B — Coarsened West raster (for R mapping)
#       Reads each year's array back from Part A's local GeoTIFF (not from
#       arraylake again), coarsens to ~1 km, checkpoints each year's array
#       to scratch (resume-safe), then assembles into one NetCDF.
# 5.  Part C — Fire polygon extraction (for R event-study / DiD)
#       Precompute rasterized polygon masks once; then extract mean AGB
#       per fire x year using local GeoTIFF reads + numpy indexing (no
#       shapely per year, no arraylake reads).
# 6.  Sanity checks on all outputs
#
# OUTPUTS
#   data/processed/ctrees/ctrees_YYYY_west_100m.tif          — native 100 m annual GeoTIFFs (raw)
#   data/processed/ctrees/ctrees_biomass_west_1km.nc        — coarsened West raster, 26 years
#   data/processed/ctrees/biomass_fire_polygons_ctrees_west.csv — long panel: event_id x year x agb
#
# RESOURCE NOTES (vs. 03_download_ctrees_ca.py)
#   West bbox is ~4.1x the CA bbox by area (~511M px/year vs ~125M px/year at
#   ~100 m). Per-year raw read is ~2 GB (float32). Do NOT accumulate multiple
#   years of the raw (uncoarsened) array at once. Part A's 26 compressed
#   GeoTIFFs total on the order of 8-20 GB on disk. Runtime is dominated by
#   ~26 network reads at 4x the size, plus mask precompute over ~6.8k fires
#   (vs. 1.1k for CA) — budget for a multi-hour run. Disable sleep before
#   starting (`powercfg /change standby-timeout-ac 0`, matching the eMapR
#   download guidance in DATA_DOWNLOAD_GUIDE.md) and run from a real
#   terminal, not a notebook — Part A/B/C all checkpoint per-year/per-fire.
#
#   GRIT note: multiple runs were killed (OOM) against a 4 GiB per-session
#   memory cap — first during the coarsening step, then even during Part A's
#   plain raw download (no coarsening math at all). Root cause both times was
#   materializing a full ~2 GB year array: (1) West's dimensions aren't exact
#   multiples of COARSEN (25650 cols / 11 truncates to 25641), so a whole-
#   array reshape for coarsening was non-contiguous and numpy silently
#   copied the entire truncated array to satisfy it; (2) even without any
#   coarsening, holding one ~2 GB float32 year array (from an int16->float32
#   cast plus a same-shape boolean fill-mask) alongside library import
#   overhead (geopandas/rasterio/arraylake/xarray, easily several hundred
#   MB-1 GB) was enough to exceed the cap on its own. Both Part A (raw
#   download) and Part B (coarsening) now read/write in row-strips
#   (STRIP_ROWS for A, `factor`-row windows in coarsen_year_from_tif for B)
#   so no step ever holds more than one strip (~tens-hundreds of MB) instead
#   of a whole year array. Part C still loads a full raster per year (needed
#   for scattered polygon indexing across the whole extent) — if that OOMs
#   too, it would need the same per-fire windowed-read treatment. If OOM
#   kills persist after all this, ask GRIT admin for a higher memory cap.
#
# PART B CHECKPOINTING (new vs. 03)
#   03's coarsening step holds all 26 coarsened years in memory and writes
#   the NetCDF once at the end — a fine tradeoff at CA scale, but risky at
#   West scale given a multi-hour runtime (the eMapR West crop was already
#   interrupted twice by laptop sleep at a similar wall-clock scale). Each
#   coarsened year is now saved to a scratch .npy immediately after computing
#   it and skipped on re-run if already present; the final NetCDF assembly
#   step only runs once all years are checkpointed.
#
# EXTRACTION METHOD (Part C) — unchanged from 03; see that script for the
# shapely-vs-rasterio evaluation notes.
#
# DATASET NOTES (from 02_explore_ctrees_zarr.py)
#   - Group:        aboveground_biomass/
#   - Variable:     agb  shape=(26, 202500, 405000)  dtype=int16
#   - Coordinates:  time (2000-2025), x (lon), y (lat, descending 90 to -90)
#   - Scale factor: divide stored int16 by 10 to get Mg ha^-1
#   - Fill value:   -9999
#   - CRS:          WGS84 / EPSG:4326
#   - Resolution:   ~0.000889 degrees ~= 100 m
# =============================================================================

# --- 1. SETUP -----------------------------------------------------------------
import sys
sys.stdout.reconfigure(encoding="utf-8")

from arraylake import Client
import zarr
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import netCDF4          # required by xarray for NetCDF write
import shapely
from pathlib import Path

# Try rasterio for fast scanline polygon rasterization
try:
    import rasterio
    import rasterio.features
    from rasterio.transform import from_origin as _rio_from_origin
    from rasterio.windows import Window
    USE_RASTERIO = True
except ImportError:
    from matplotlib.path import Path as MplPath
    USE_RASTERIO = False

# Part A's raw GeoTIFFs are now a hard dependency for Part B/C (they read the
# raw array back from disk instead of arraylake) — rasterio is required to
# both write and read those TIFFs, so there's no fallback path here anymore.
if not USE_RASTERIO:
    print("ERROR: rasterio is required (it now backs the raw-GeoTIFF download\n"
          "that Part B/C read from). Install it with: pip install rasterio")
    sys.exit(1)

PROJ_ROOT     = Path(__file__).resolve().parent.parent.parent
OUT_TIFS_DIR  = PROJ_ROOT / "data" / "processed" / "ctrees"
OUT_NC        = OUT_TIFS_DIR / "ctrees_biomass_west_1km.nc"
OUT_CSV       = OUT_TIFS_DIR / "biomass_fire_polygons_ctrees_west.csv"
MTBS_PATH     = PROJ_ROOT / "data" / "raw" / "mtbs" / "mtbs_perimeter_data" / "mtbs_perims_DD.shp"

# Scratch dir for Part B per-year checkpoints (deleted after successful NetCDF assembly)
SCRATCH_DIR   = OUT_TIFS_DIR / "_west_1km_scratch"

# -- Ctrees zarr parameters (from exploration script) -------------------------
REPO_NAME    = "ucsb-emlab/BACI-wildfires"
BRANCH       = "main"
GROUP        = "aboveground_biomass"
VAR          = "agb"
SCALE_FACTOR = 10.0     # divide stored int16 by 10 -> Mg ha^-1
FILL_VALUE   = -9999

# -- Western study states (event_id prefix convention, matches WESTERN_STATES
#    in scripts/r/00_crop_emapr_to_west.R and the R extraction pipeline) -----
WESTERN_STATES = ["AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"]

# -- Western US bounding box (WGS84 degrees) — union of the 11 states' extents
WEST_LON = (-124.8, -102.0)
WEST_LAT = (31.3, 49.0)

# -- Coarsening factor: 11 x 0.000889 deg ~= 0.0098 deg ~= 1.1 km -----------
COARSEN = 11

# -- Row-strip size for Part A's read/write loop: bounds peak memory to one
#    strip (~200 MB at this size) instead of a whole ~2 GB year array. Even
#    the plain raw download (no coarsening math at all) was OOM-killed
#    against GRIT's 4 GiB cap, so Part A needs this too, not just Part B.
STRIP_ROWS = 2000

# -- MTBS filter --------------------------------------------------------------
MTBS_START = 2000
MTBS_END   = 2025    # match Ctrees temporal range

OUT_TIFS_DIR.mkdir(parents=True, exist_ok=True)
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)


def tif_path_for_year(yr):
    return OUT_TIFS_DIR / f"ctrees_{yr}_west_100m.tif"


def read_year_raster(yr):
    """
    Read one year's raw West-bbox array back from Part A's local GeoTIFF.

    The TIFF already has fill values converted to NaN and the int16->Mg ha^-1
    scale factor applied (done once, in Part A, before writing) — so this is
    a plain local raster read, no arraylake connection and no re-masking.
    Used by Part C, which needs the full array in memory anyway for
    scattered polygon indexing across the whole West extent.
    """
    with rasterio.open(tif_path_for_year(yr)) as src:
        return src.read(1).astype("float32")


def coarsen_year_from_tif(yr, factor):
    """
    Build one year's coarsened (~1 km) raster directly from Part A's local
    GeoTIFF, reading and averaging `factor` rows at a time via a rasterio
    window — never materializes the full ~2 GB year array at all.

    An earlier version loaded the whole raster into memory first
    (`read_year_raster`) and then block-averaged it with a numpy reshape.
    West's dimensions aren't exact multiples of `factor` (25650 cols / 11
    truncates to 25641), so that reshape was non-contiguous and numpy
    silently copied the entire truncated array to satisfy it — briefly
    doubling peak memory on top of the array already in memory, which is
    what OOM-killed runs against GRIT's 4 GiB cap. Reading window-by-window
    from disk sidesteps the whole-array reshape problem rather than just
    shrinking it.
    """
    with rasterio.open(tif_path_for_year(yr)) as src:
        nx_c = src.width // factor
        ny_c = src.height // factor
        nx_trim = nx_c * factor

        out = np.empty((ny_c, nx_c), dtype="float32")
        for i in range(ny_c):
            window = Window(0, i * factor, nx_trim, factor)
            strip = src.read(1, window=window)   # (factor, nx_trim)
            out[i] = strip.reshape(factor, nx_c, factor).mean(axis=(0, 2))
    return out


def precompute_mask(geom, x_arr, y_arr, res):
    """
    Rasterize a polygon onto the West pixel grid; return (yi, xi, mask_2d).

    Uses rasterio.features.rasterize() for fast C-level scanline rasterization
    (O(pixels), no sampling approximation). Handles MultiPolygon by unioning
    sub-polygon masks. Returns None when the polygon does not overlap the grid.

    Called once per fire before the year loop; the returned mask is reused
    across all 26 years using numpy array indexing (no shapely per year).

    Parameters
    ----------
    geom  : shapely geometry (Polygon or MultiPolygon), WGS84
    x_arr : 1D float64 array — x (lon) pixel centres, ascending
    y_arr : 1D float64 array — y (lat) pixel centres, descending (N first)
    res   : float — pixel spacing in degrees (~0.000889)
    """
    xmin, ymin, xmax, ymax = geom.bounds

    xi = np.where((x_arr >= xmin) & (x_arr <= xmax))[0]
    yi = np.where((y_arr >= ymin) & (y_arr <= ymax))[0]

    if xi.size == 0 or yi.size == 0:
        return None

    height, width = len(yi), len(xi)

    # from_origin(west, north, xsize, ysize): west/north are the top-left
    # corner of the top-left pixel. y_arr is descending so yi[0] is the
    # northernmost row.
    west  = float(x_arr[xi[0]]) - res / 2
    north = float(y_arr[yi[0]]) + res / 2
    transform = _rio_from_origin(west, north, res, res)
    mask = rasterio.features.rasterize(
        [(geom, 1)],
        out_shape=(height, width),
        transform=transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)

    if not mask.any():
        return None

    return yi, xi, mask


# --- 2. CONNECT TO ARRAYLAKE & RESOLVE WEST INDICES --------------------------
print("Connecting to arraylake...", flush=True)
try:
    client  = Client()
    repo    = client.get_repo(REPO_NAME)
    session = repo.readonly_session(branch=BRANCH)
    root    = zarr.open_group(session.store, zarr_format=3, mode="r")
    print("  Connected.", flush=True)
except Exception as e:
    print(f"ERROR: {e}\nRun `arraylake auth login` and retry.")
    sys.exit(1)

agb_zarr = root[f"{GROUP}/{VAR}"]     # (26, 202500, 405000), int16
time_raw = root[f"{GROUP}/time"][:]   # datetime64[D]
x_all    = root[f"{GROUP}/x"][:]      # ascending, -180 to +180
y_all    = root[f"{GROUP}/y"][:]      # descending, 90 to -90

n_years = agb_zarr.shape[0]
times   = pd.DatetimeIndex(time_raw.astype("datetime64[ns]"))
years   = times.year.tolist()

# West index slices (y is descending so north-side has smaller index)
x_mask = (x_all >= WEST_LON[0]) & (x_all <= WEST_LON[1])
y_mask = (y_all >= WEST_LAT[0]) & (y_all <= WEST_LAT[1])

x_idx  = np.where(x_mask)[0];  x_start, x_end = int(x_idx[0]), int(x_idx[-1]) + 1
y_idx  = np.where(y_mask)[0];  y_start, y_end = int(y_idx[0]), int(y_idx[-1]) + 1

x_west = x_all[x_idx]   # ascending
y_west = y_all[y_idx]   # descending (northernmost first)

print(f"  West subset: {y_end - y_start} rows x {x_end - x_start} cols x {n_years} years")
print(f"  Lat range: {y_west.min():.3f} - {y_west.max():.3f}")
print(f"  Lon range: {x_west.min():.3f} - {x_west.max():.3f}")


# --- 3. PART A — DOWNLOAD NATIVE-RESOLUTION (~100 m) RAW GeoTIFFs ------------
# Runs first: this is the only section that reads pixel data from arraylake.
# Writes one GeoTIFF per zarr year at the native ~100 m pixel grid, full West
# bbox. Meant to be shared across all 11 states once scripts/r/07 and 08 loop
# STATE_FIPS over WESTERN_STATES — each state crops its slice from the same
# per-year West TIF rather than needing a separate download per state.

RES = float(abs(x_west[1] - x_west[0]))   # ~0.000889 degrees, shared by A/B/C
west_c  = float(x_west[0]) - RES / 2
north_c = float(y_west[0]) + RES / 2
transform_c = _rio_from_origin(west_c, north_c, RES, RES)

tifs_needed = [yr for yr in years if not tif_path_for_year(yr).exists()]
if not tifs_needed:
    print(f"\nPart A: All {len(years)} raw 100 m TIFs already exist — skipping.")
else:
    print(f"\nPart A: Downloading {len(tifs_needed)} native-resolution (~100 m) raw GeoTIFFs"
          f" -> {OUT_TIFS_DIR.name}/")

    n_rows_total = y_end - y_start
    n_cols_total = x_end - x_start

    for t_idx, yr in enumerate(years):
        tif_path = tif_path_for_year(yr)
        if tif_path.exists():
            print(f"  {t_idx + 1}/{n_years} years — {yr} already downloaded, skipping", flush=True)
            continue

        # Written row-strip by row-strip (STRIP_ROWS at a time) rather than
        # building the whole ~2 GB year array in memory first — see
        # STRIP_ROWS' comment above for why even this plain download needs it.
        with rasterio.open(
            tif_path, "w",
            driver     = "GTiff",
            height     = n_rows_total,
            width      = n_cols_total,
            count      = 1,
            dtype      = "float32",
            crs        = "EPSG:4326",
            transform  = transform_c,
            compress   = "lzw",
            tiled      = True,
            blockxsize = 512,
            blockysize = 512,
            nodata     = float("nan"),
        ) as dst:
            for row_off in range(0, n_rows_total, STRIP_ROWS):
                row_end = min(row_off + STRIP_ROWS, n_rows_total)
                strip = agb_zarr[t_idx,
                                 y_start + row_off:y_start + row_end,
                                 x_start:x_end].astype("float32")
                strip[strip == FILL_VALUE] = np.nan
                strip /= SCALE_FACTOR   # -> Mg ha^-1
                dst.write(strip, 1, window=Window(0, row_off, n_cols_total, row_end - row_off))
                del strip

        size_mb = tif_path.stat().st_size / 1e6
        print(f"  {t_idx + 1}/{n_years} years downloaded ({yr}) — {size_mb:.1f} MB", flush=True)

    n_done = len(list(OUT_TIFS_DIR.glob("ctrees_*_west_100m.tif")))
    print(f"  Done. {n_done}/{len(years)} raw 100 m TIFs present in {OUT_TIFS_DIR.name}/")

assert all(tif_path_for_year(yr).exists() for yr in years), (
    "Part A did not produce a raw TIFF for every year — Part B/C below "
    "require all years present locally before they can proceed."
)


# --- 4. PART B — COARSENED WEST RASTER (reads Part A's local TIFFs) ---------
if OUT_NC.exists():
    print(f"\nPart B: {OUT_NC.name} already exists — skipping.")
else:
    print(f"\nPart B: Building coarsened (~1 km) West raster -> {OUT_NC.name}")

    for t_idx, yr in enumerate(years):
        scratch_path = SCRATCH_DIR / f"coarsened_{yr}.npy"
        if scratch_path.exists():
            print(f"  {t_idx + 1}/{n_years} years — {yr} already checkpointed, skipping", flush=True)
            continue

        coarsened = coarsen_year_from_tif(yr, COARSEN)

        np.save(scratch_path, coarsened)
        print(f"  {t_idx + 1}/{n_years} years processed ({yr}) — checkpointed", flush=True)

    # All years checkpointed — assemble into one NetCDF
    coarsened_layers = [np.load(SCRATCH_DIR / f"coarsened_{yr}.npy") for yr in years]

    ny_c, nx_c = coarsened_layers[0].shape
    x_c = x_west[: nx_c * COARSEN].reshape(nx_c, COARSEN).mean(axis=1)
    y_c = y_west[: ny_c * COARSEN].reshape(ny_c, COARSEN).mean(axis=1)

    agb_stack = np.stack(coarsened_layers, axis=0)   # (26, ny_c, nx_c)

    ds_out = xr.Dataset(
        {"agb": (["time", "y", "x"], agb_stack.astype("float32"))},
        coords={
            "time": times.values,
            "y":    ("y", y_c),
            "x":    ("x", x_c),
        },
        attrs={
            "title":        "Ctrees aboveground biomass — Western US ~1 km",
            "source":       "ucsb-emlab/BACI-wildfires (arraylake)",
            "units":        "Mg ha-1",
            "scale_note":   "Coarsened to ~1 km by block-averaging native 100 m pixels",
            "crs":          "EPSG:4326 (WGS84)",
        }
    )
    ds_out["agb"].attrs.update({"units": "Mg ha-1", "long_name": "Aboveground Biomass",
                                "_FillValue": -9999.0})

    encoding = {"agb": {"dtype": "float32", "zlib": True, "complevel": 4}}
    ds_out.to_netcdf(OUT_NC, encoding=encoding)
    print(f"  Saved: {OUT_NC}  ({OUT_NC.stat().st_size / 1e6:.1f} MB)")

    # Clean up scratch checkpoints now that the NetCDF is safely written
    for yr in years:
        (SCRATCH_DIR / f"coarsened_{yr}.npy").unlink(missing_ok=True)
    try:
        SCRATCH_DIR.rmdir()
    except OSError:
        pass   # leave it if anything unexpected remains


# --- 5. PART C — FIRE POLYGON EXTRACTION (reads Part A's local TIFFs) -------
if OUT_CSV.exists():
    print(f"\nPart C: {OUT_CSV.name} already exists — skipping.")
else:
    print(f"\nPart C: Extracting AGB within MTBS West fire polygons -> {OUT_CSV.name}")
    assert MTBS_PATH.exists(), f"MTBS shapefile not found: {MTBS_PATH}"

    # -- Load and filter MTBS Western US wildfires -----------------------------
    mtbs_raw = gpd.read_file(MTBS_PATH)
    mtbs_raw.columns = [c.lower() for c in mtbs_raw.columns]
    mtbs_raw["fire_year"] = mtbs_raw["ig_date"].str[:4].astype(int)

    mtbs_west = (mtbs_raw
                 .loc[mtbs_raw["incid_type"] == "Wildfire"]
                 .loc[mtbs_raw["fire_year"].between(MTBS_START, MTBS_END)]
                 .loc[mtbs_raw["event_id"].str[:2].isin(WESTERN_STATES)]
                 .to_crs("EPSG:4326")
                 .reset_index(drop=True))

    print(f"  West wildfires in MTBS ({MTBS_START}-{MTBS_END}): {len(mtbs_west)}", flush=True)
    assert len(mtbs_west) > 0, "No Western US wildfires found — check MTBS path and filters"

    # -- Precompute polygon masks (once per fire, reused across 26 years) -----
    # This is the key optimisation: rasterize each polygon onto the West pixel
    # grid a single time, storing (row_indices, col_indices, boolean_mask).
    # The year loop then uses only numpy indexing — no shapely or rasterio work.
    print(f"  Precomputing {len(mtbs_west)} polygon masks...", flush=True)
    fire_masks = []
    for _, fire in mtbs_west.iterrows():
        fire_masks.append((fire, precompute_mask(fire.geometry, x_west, y_west, RES)))

    n_with_mask = sum(1 for _, m in fire_masks if m is not None)
    print(f"  {n_with_mask}/{len(mtbs_west)} fire polygons overlap the West raster grid", flush=True)

    # -- Extract year by year using precomputed masks --------------------------
    # Checkpointed per year (mirrors Part B) — two prior runs of this script
    # were killed unexpectedly after several hours, once mid-extraction, losing
    # all in-memory progress since the original version only wrote the CSV
    # once at the end. Each year's records now land in SCRATCH_DIR_B
    # immediately and are skipped on re-run if already present.
    SCRATCH_DIR_B = OUT_TIFS_DIR / "_west_fireagb_scratch"
    SCRATCH_DIR_B.mkdir(parents=True, exist_ok=True)
    errors = []

    for t_idx, yr in enumerate(years):
        scratch_path = SCRATCH_DIR_B / f"records_{yr}.csv"
        if scratch_path.exists():
            print(f"  {t_idx + 1}/{n_years} years — {yr} already checkpointed, skipping", flush=True)
            continue

        # One local raster read per year (the only I/O in this loop)
        raw = read_year_raster(yr)

        year_records = []
        for fire, mask_result in fire_masks:
            if mask_result is None:
                mean_agb = np.nan
            else:
                yi, xi, mask = mask_result
                vals = raw[np.ix_(yi, xi)][mask]
                if vals.size > 0 and not np.all(np.isnan(vals)):
                    mean_agb = float(np.nanmean(vals))
                else:
                    mean_agb = np.nan

            year_records.append({
                "event_id":      fire["event_id"],
                "state":         fire["event_id"][:2],
                "fire_year":     fire["fire_year"],
                "year":          yr,
                "agb_mean_mgha": round(mean_agb, 2) if not np.isnan(mean_agb) else np.nan,
            })

        del raw
        pd.DataFrame(year_records).to_csv(scratch_path, index=False)
        print(f"  {t_idx + 1}/{n_years} years processed ({yr}) — checkpointed", flush=True)

    # All years checkpointed — assemble into one CSV
    df = pd.concat([pd.read_csv(SCRATCH_DIR_B / f"records_{yr}.csv") for yr in years],
                    ignore_index=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"  Saved: {OUT_CSV}  ({len(df):,} records)")

    # Clean up scratch checkpoints now that the CSV is safely written
    for yr in years:
        (SCRATCH_DIR_B / f"records_{yr}.csv").unlink(missing_ok=True)
    try:
        SCRATCH_DIR_B.rmdir()
    except OSError:
        pass   # leave it if anything unexpected remains

    if errors:
        print(f"  WARNING: {len(errors)} extraction errors (first 3):")
        for e in errors[:3]:
            print(f"    {e}")


# --- 6. SANITY CHECKS --------------------------------------------------------
print("\nSanity checks...")

# Part A
tif_count = len(list(OUT_TIFS_DIR.glob("ctrees_*_west_100m.tif")))
print(f"  Raw 100 m TIFs: {tif_count}/{n_years} year(s) in {OUT_TIFS_DIR.name}/")
if tif_count < n_years:
    print(f"  WARNING: Only {tif_count} of {n_years} raw 100 m TIFs present")

# Part B
ds_check = xr.open_dataset(OUT_NC)
assert "agb" in ds_check, "NetCDF missing 'agb' variable"
assert len(ds_check.time) == n_years, f"Expected {n_years} time steps"
agb_vals = ds_check["agb"].values
pct_valid = 100 * np.sum(~np.isnan(agb_vals)) / agb_vals.size
print(f"  NetCDF: {len(ds_check.time)} years, "
      f"{ds_check.dims['y']} x {ds_check.dims['x']} pixels, "
      f"{pct_valid:.1f}% valid")
assert pct_valid > 10, "Less than 10% valid pixels — check West bbox or fill masking"

# Part C
df_check = pd.read_csv(OUT_CSV)
assert {"event_id", "state", "fire_year", "year", "agb_mean_mgha"}.issubset(df_check.columns)
pct_valid_csv = 100 * df_check["agb_mean_mgha"].notna().mean()
print(f"  CSV: {df_check['event_id'].nunique()} fires x {df_check['year'].nunique()} years "
      f"= {len(df_check):,} records, {pct_valid_csv:.1f}% with valid AGB")
print("  Fires per state:")
print(df_check.drop_duplicates("event_id")["state"].value_counts().sort_index().to_string())
if pct_valid_csv < 50:
    print("  WARNING: <50% valid — check polygon alignment with Ctrees grid")

print("\nDone. Outputs ready for analysis/")
