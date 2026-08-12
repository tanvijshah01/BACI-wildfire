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
# Part C's per-year GeoTIFFs (ctrees_YYYY_west_100m.tif) are meant to be
# shared across all 11 states — scripts/r/07 and 08 will eventually loop
# STATE_FIPS over WESTERN_STATES and crop this same West-wide TIF per state
# per year, rather than needing one ctrees TIF per state.
#
# SCRIPT OUTLINE
# 1.  Setup — connect to arraylake, define parameters
# 2.  Open zarr store and resolve West bounding-box indices
# 3.  Part A — Coarsened West raster (for R mapping)
#       Load West data year-by-year, coarsen to ~1 km, checkpoint each year's
#       array to scratch (resume-safe — a multi-hour run at 4x the CA-script's
#       area is more exposed to interruption), then assemble into one NetCDF.
# 4.  Part B — Fire polygon extraction (for R event-study / DiD)
#       Precompute rasterized polygon masks once; then extract mean AGB
#       per fire x year using pure numpy indexing (no shapely per year).
# 5.  Part C — Native-resolution (~100 m) annual GeoTIFFs (for R comparison)
#       Write one GeoTIFF per year at native 100 m resolution, West-bbox
#       extent — shared input for scripts/r/07 and 08 once generalized to
#       loop over all 11 states.
# 6.  Sanity checks on all outputs
#
# OUTPUTS
#   data/processed/ctrees/ctrees_biomass_west_1km.nc        — coarsened West raster, 26 years
#   data/processed/ctrees/biomass_fire_polygons_ctrees_west.csv — long panel: event_id x year x agb
#   data/processed/ctrees/ctrees_YYYY_west_100m.tif          — native 100 m annual GeoTIFFs
#
# RESOURCE NOTES (vs. 03_download_ctrees_ca.py)
#   West bbox is ~4.1x the CA bbox by area (~511M px/year vs ~125M px/year at
#   ~100 m). Per-year raw read is ~2 GB (float32) — fine in memory, but do NOT
#   accumulate multiple years of the raw (uncoarsened) array at once. Part C's
#   26 compressed GeoTIFFs total on the order of 8-20 GB on disk. Runtime is
#   dominated by ~26 network reads at 4x the size, plus mask precompute over
#   ~6.8k fires (vs. 1.1k for CA) — budget for a multi-hour run. Disable sleep
#   before starting (`powercfg /change standby-timeout-ac 0`, matching the
#   eMapR download guidance in DATA_DOWNLOAD_GUIDE.md) and run from a real
#   terminal, not a notebook — Part B/C already checkpoint per-fire/per-year,
#   and Part A now checkpoints per-year too (see PART A CHECKPOINTING below).
#
# PART A CHECKPOINTING (new vs. 03)
#   03's Part A holds all 26 coarsened years in memory and writes the NetCDF
#   once at the end — a fine tradeoff at CA scale, but risky at West scale
#   given a multi-hour runtime (the eMapR West crop was already interrupted
#   twice by laptop sleep at a similar wall-clock scale). Each coarsened year
#   is now saved to a scratch .npy immediately after computing it and skipped
#   on re-run if already present; the final NetCDF assembly step only runs
#   once all years are checkpointed.
#
# EXTRACTION METHOD (Part B) — unchanged from 03; see that script for the
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
    USE_RASTERIO = True
except ImportError:
    from matplotlib.path import Path as MplPath
    USE_RASTERIO = False

PROJ_ROOT     = Path(__file__).resolve().parent.parent.parent
OUT_TIFS_DIR  = PROJ_ROOT / "data" / "processed" / "ctrees"
OUT_NC        = OUT_TIFS_DIR / "ctrees_biomass_west_1km.nc"
OUT_CSV       = OUT_TIFS_DIR / "biomass_fire_polygons_ctrees_west.csv"
MTBS_PATH     = PROJ_ROOT / "data" / "raw" / "mtbs" / "mtbs_perimeter_data" / "mtbs_perims_DD.shp"

# Scratch dir for Part A per-year checkpoints (deleted after successful NetCDF assembly)
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

# -- MTBS filter --------------------------------------------------------------
MTBS_START = 2000
MTBS_END   = 2025    # match Ctrees temporal range

OUT_TIFS_DIR.mkdir(parents=True, exist_ok=True)
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)


def coarsen_block(arr2d, factor):
    """Block-average a 2D numpy array by `factor` in each dimension."""
    ny, nx = arr2d.shape
    ny_c, nx_c = ny // factor, nx // factor
    arr2d = arr2d[: ny_c * factor, : nx_c * factor]
    return arr2d.reshape(ny_c, factor, nx_c, factor).mean(axis=(1, 3))


def precompute_mask(geom, x_arr, y_arr, res):
    """
    Rasterize a polygon onto the West pixel grid; return (yi, xi, mask_2d).

    Uses rasterio.features.rasterize() for fast C-level scanline rasterization
    (O(pixels), no sampling approximation). Falls back to matplotlib.path if
    rasterio is not available. Handles MultiPolygon by unioning sub-polygon
    masks. Returns None when the polygon does not overlap the grid.

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

    if USE_RASTERIO:
        # from_origin(west, north, xsize, ysize): west/north are the top-left
        # corner of the top-left pixel.  y_arr is descending so yi[0] is the
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
    else:
        # matplotlib.path fallback: vectorised C extension, O(pixels x edges).
        # Thin the pixel grid if the bounding box exceeds MAX_FALLBACK points —
        # masks are precomputed once so a generous cap (500k) keeps accuracy
        # high while preventing multi-minute hangs on 1M-acre fires.
        MAX_FALLBACK = 500_000
        n_bbox = height * width
        if n_bbox > MAX_FALLBACK:
            step = max(1, int(np.ceil(np.sqrt(n_bbox / MAX_FALLBACK))))
            xi = xi[::step]
            yi = yi[::step]
            height, width = len(yi), len(xi)

        X_grid, Y_grid = np.meshgrid(x_arr[xi], y_arr[yi])
        pts = np.column_stack([X_grid.ravel(), Y_grid.ravel()])

        if geom.geom_type == "MultiPolygon":
            in_poly = np.zeros(len(pts), dtype=bool)
            for part in geom.geoms:
                in_poly |= MplPath(np.array(part.exterior.coords)).contains_points(pts)
        else:
            in_poly = MplPath(np.array(geom.exterior.coords)).contains_points(pts)
        mask = in_poly.reshape(height, width)

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
print(f"  Extraction backend: {'rasterio.features' if USE_RASTERIO else 'matplotlib.path (fallback)'}")


# --- 3. PART A — COARSENED WEST RASTER ----------------------------------------
if OUT_NC.exists():
    print(f"\nPart A: {OUT_NC.name} already exists — skipping.")
else:
    print(f"\nPart A: Building coarsened (~1 km) West raster -> {OUT_NC.name}")

    for t_idx, (yr, ts) in enumerate(zip(years, times)):
        scratch_path = SCRATCH_DIR / f"coarsened_{yr}.npy"
        if scratch_path.exists():
            print(f"  {t_idx + 1}/{n_years} years — {yr} already checkpointed, skipping", flush=True)
            continue

        raw = agb_zarr[t_idx, y_start:y_end, x_start:x_end].astype("float32")
        raw[raw == FILL_VALUE] = np.nan
        raw /= SCALE_FACTOR                          # -> Mg ha^-1
        coarsened = coarsen_block(raw, COARSEN)
        del raw

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


# --- 4. PART B — FIRE POLYGON EXTRACTION -------------------------------------
if OUT_CSV.exists():
    print(f"\nPart B: {OUT_CSV.name} already exists — skipping.")
else:
    print(f"\nPart B: Extracting AGB within MTBS West fire polygons -> {OUT_CSV.name}")
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
    RES = float(abs(x_west[1] - x_west[0]))   # ~0.000889 degrees

    print(f"  Precomputing {len(mtbs_west)} polygon masks...", flush=True)
    fire_masks = []
    for _, fire in mtbs_west.iterrows():
        fire_masks.append((fire, precompute_mask(fire.geometry, x_west, y_west, RES)))

    n_with_mask = sum(1 for _, m in fire_masks if m is not None)
    print(f"  {n_with_mask}/{len(mtbs_west)} fire polygons overlap the West raster grid", flush=True)

    # -- Extract year by year using precomputed masks -------------------------
    records = []
    errors  = []

    for t_idx, (yr, ts) in enumerate(zip(years, times)):
        # One zarr read per year (the only I/O in this loop)
        raw = agb_zarr[t_idx, y_start:y_end, x_start:x_end].astype("float32")
        raw[raw == FILL_VALUE] = np.nan
        raw /= SCALE_FACTOR   # -> Mg ha^-1

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

            records.append({
                "event_id":      fire["event_id"],
                "state":         fire["event_id"][:2],
                "fire_year":     fire["fire_year"],
                "year":          yr,
                "agb_mean_mgha": round(mean_agb, 2) if not np.isnan(mean_agb) else np.nan,
            })

        del raw
        print(f"  {t_idx + 1}/{n_years} years processed ({yr})", flush=True)

    df = pd.DataFrame(records)
    df.to_csv(OUT_CSV, index=False)
    print(f"  Saved: {OUT_CSV}  ({len(df):,} records)")

    if errors:
        print(f"  WARNING: {len(errors)} extraction errors (first 3):")
        for e in errors[:3]:
            print(f"    {e}")


# --- 5. PART C — NATIVE-RESOLUTION (~100 m) ANNUAL GeoTIFFs -----------------
# Writes one GeoTIFF per zarr year at the native ~100 m pixel grid, full West
# bbox. Meant to be shared across all 11 states once scripts/r/07 and 08 loop
# STATE_FIPS over WESTERN_STATES — each state crops its slice from the same
# per-year West TIF rather than needing a separate download per state.
#
# Requires rasterio (USE_RASTERIO=True). If rasterio is not installed, Part C
# is skipped — install it with: pip install rasterio

if not USE_RASTERIO:
    print("\nPart C: rasterio not available — skipping 100 m GeoTIFF export.")
    print("  Install rasterio (pip install rasterio) and re-run to generate")
    print("  ctrees_YYYY_west_100m.tif files.")
else:
    tifs_needed = [yr for yr in years
                   if not (OUT_TIFS_DIR / f"ctrees_{yr}_west_100m.tif").exists()]
    if not tifs_needed:
        print(f"\nPart C: All {len(years)} 100 m TIFs already exist — skipping.")
    else:
        print(f"\nPart C: Writing {len(tifs_needed)} native-resolution (~100 m) GeoTIFFs"
              f" -> {OUT_TIFS_DIR.name}/")

        # Pixel spacing in degrees (~0.000889); all years share the same grid
        RES_C = float(abs(x_west[1] - x_west[0]))
        # rasterio transform: origin = upper-left corner of top-left pixel;
        # y_west is descending so y_west[0] is the northernmost (largest) latitude.
        west_c  = float(x_west[0]) - RES_C / 2
        north_c = float(y_west[0]) + RES_C / 2
        transform_c = _rio_from_origin(west_c, north_c, RES_C, RES_C)

        for t_idx, yr in enumerate(years):
            tif_path = OUT_TIFS_DIR / f"ctrees_{yr}_west_100m.tif"
            if tif_path.exists():
                print(f"  {yr} — already exists, skipping", flush=True)
                continue

            raw = agb_zarr[t_idx, y_start:y_end, x_start:x_end].astype("float32")
            raw[raw == FILL_VALUE] = np.nan
            raw /= SCALE_FACTOR   # -> Mg ha^-1

            with rasterio.open(
                tif_path, "w",
                driver     = "GTiff",
                height     = raw.shape[0],
                width      = raw.shape[1],
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
                dst.write(raw, 1)

            del raw
            size_mb = tif_path.stat().st_size / 1e6
            print(f"  {yr} — {size_mb:.1f} MB", flush=True)

        n_done = len(list(OUT_TIFS_DIR.glob("ctrees_*_west_100m.tif")))
        print(f"  Done. {n_done}/{len(years)} 100 m TIFs present in {OUT_TIFS_DIR.name}/")


# --- 6. SANITY CHECKS --------------------------------------------------------
print("\nSanity checks...")

# Part A
ds_check = xr.open_dataset(OUT_NC)
assert "agb" in ds_check, "NetCDF missing 'agb' variable"
assert len(ds_check.time) == n_years, f"Expected {n_years} time steps"
agb_vals = ds_check["agb"].values
pct_valid = 100 * np.sum(~np.isnan(agb_vals)) / agb_vals.size
print(f"  NetCDF: {len(ds_check.time)} years, "
      f"{ds_check.dims['y']} x {ds_check.dims['x']} pixels, "
      f"{pct_valid:.1f}% valid")
assert pct_valid > 10, "Less than 10% valid pixels — check West bbox or fill masking"

# Part B
df_check = pd.read_csv(OUT_CSV)
assert {"event_id", "state", "fire_year", "year", "agb_mean_mgha"}.issubset(df_check.columns)
pct_valid_csv = 100 * df_check["agb_mean_mgha"].notna().mean()
print(f"  CSV: {df_check['event_id'].nunique()} fires x {df_check['year'].nunique()} years "
      f"= {len(df_check):,} records, {pct_valid_csv:.1f}% with valid AGB")
print("  Fires per state:")
print(df_check.drop_duplicates("event_id")["state"].value_counts().sort_index().to_string())
if pct_valid_csv < 50:
    print("  WARNING: <50% valid — check polygon alignment with Ctrees grid")

# Part C
if USE_RASTERIO:
    tif_count = len(list(OUT_TIFS_DIR.glob("ctrees_*_west_100m.tif")))
    print(f"  100 m TIFs: {tif_count}/{n_years} year(s) in {OUT_TIFS_DIR.name}/")
    if tif_count < n_years:
        print(f"  WARNING: Only {tif_count} of {n_years} 100 m TIFs present")
else:
    print("  100 m TIFs: skipped (rasterio not installed)")

print("\nDone. Outputs ready for analysis/")
