# =============================================================================
# 01_extract_biomass_gee.py
# Extract annual fitted NBR (biomass proxy) from LandTrendR / Google Earth Engine
# for MTBS wildfire perimeters in the 11 Western US states, 2000-2023.
#
# SCRIPT OUTLINE
# 1.  Setup — imports, paths, study parameters
# 2.  GEE authentication and initialization
# 3.  Load and filter MTBS fire perimeters
# 4.  Sample fires for EDA (stratified by year)
# 5.  Landsat collection builder (harmonized across missions)
# 6.  LandTrendR extraction function (per fire)
# 7.  Run extraction loop and save to output/biomass_timeseries.csv
#
# BEFORE RUNNING:
#   earthengine authenticate        (one-time browser login)
#   Set GEE_PROJECT below to your Google Earth Engine project ID.
#
# OUTPUT:
#   output/biomass_timeseries.csv
#   Columns: event_id, fire_year, year, fitted_nbr
#   fitted_nbr is in INVERTED NBR space (loss = positive, as required by LT).
#   Negate when displaying: display_nbr = -fitted_nbr
# =============================================================================

# ─── 1. SETUP ─────────────────────────────────────────────────────────────────
import ee
import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# ── User must set their GEE project ID here ──────────────────────────────────
GEE_PROJECT = "your-gee-project-id"   # ← CHANGE THIS

# ── Resolve project root (two levels up from scripts/python/) ─────────────────
PROJ_ROOT = Path(__file__).resolve().parent.parent.parent

# ── Study parameters ──────────────────────────────────────────────────────────
START_YEAR  = 2000
END_YEAR    = 2023
START_DAY   = "07-01"   # July–September: peak summer growing season, Western US
END_DAY     = "09-30"
N_SAMPLE    = 200       # fires to sample for EDA; increase for full analysis
SCALE       = 90        # pixel resolution for reduceRegion (90m faster, 30m precise)

# LandTrendR algorithm parameters (see https://emapr.github.io/LT-GEE/)
RUN_PARAMS = {
    "maxSegments":            6,     # max temporal breakpoints per pixel
    "spikeThreshold":         0.9,   # dampens single-year spikes (1.0 = no dampening)
    "vertexCountOvershoot":   3,     # allows initial overshoot before pruning
    "preventOneYearRecovery": False, # allow recovery segments shorter than 1 year
    "recoveryThreshold":      0.25,  # max recovery rate (1/0.25 = 4-year minimum recovery)
    "pvalThreshold":          0.1,   # model fit p-value threshold
    "bestModelProportion":    1.25,  # model selection tolerance
    "minObservationsNeeded":  6,     # minimum annual obs required to run LT
}

# Western US approximate bounding box (pre-filter before spatial join in R)
WEST_LAT = (31.0, 49.0)
WEST_LON = (-124.5, -103.0)


# ─── 2. GEE AUTHENTICATION ────────────────────────────────────────────────────
print("Initializing Google Earth Engine...")
try:
    ee.Initialize(project=GEE_PROJECT)
    print("  GEE initialized successfully.")
except Exception as e:
    print(f"\nERROR: Could not initialize GEE.\n{e}")
    print("\nFix: Run `earthengine authenticate` in your terminal, then retry.")
    sys.exit(1)


# ─── 3. LOAD AND FILTER MTBS FIRE PERIMETERS ──────────────────────────────────
print("\nLoading MTBS fire perimeters...")
MTBS_PATH = PROJ_ROOT / "data/raw/mtbs/mtbs_perimeter_data/mtbs_perims_DD.shp"
assert MTBS_PATH.exists(), f"MTBS shapefile not found at: {MTBS_PATH}"

mtbs_raw = gpd.read_file(MTBS_PATH)
mtbs_raw["year"] = mtbs_raw["ig_date"].str[:4].astype(int)
mtbs_raw["burnbndlat"] = mtbs_raw["burnbndlat"].astype(float)
mtbs_raw["burnbndlon"] = mtbs_raw["burnbndlon"].astype(float)

mtbs = (mtbs_raw
    .loc[mtbs_raw["incid_type"] == "Wildfire"]
    .loc[mtbs_raw["year"].between(START_YEAR, END_YEAR)]
    .loc[mtbs_raw["burnbndlat"].between(*WEST_LAT)]
    .loc[mtbs_raw["burnbndlon"].between(*WEST_LON)]
    .reset_index(drop=True)
)
print(f"  Fires after filter: {len(mtbs)}")

# Sanity check: year range
assert mtbs["year"].between(START_YEAR, END_YEAR).all(), "Year range out of bounds"
assert (mtbs["incid_type"] == "Wildfire").all(), "Non-wildfire records present"


# ─── 4. SAMPLE FIRES (stratified by year for temporal coverage) ───────────────
n_per_year = max(1, N_SAMPLE // (END_YEAR - START_YEAR + 1))
sample = (mtbs
    .groupby("year", group_keys=False)
    .apply(lambda g: g.sample(min(len(g), n_per_year), random_state=42))
    .reset_index(drop=True)
)
print(f"  Sampled {len(sample)} fires across {sample['year'].nunique()} years for EDA")


# ─── 5. LANDSAT COLLECTION BUILDER ────────────────────────────────────────────
def mask_landsat_c2(image):
    """Mask cloud and cloud shadow using QA_PIXEL band (Landsat Collection 2)."""
    qa     = image.select("QA_PIXEL")
    cloud  = qa.bitwiseAnd(1 << 3)   # bit 3 = cloud
    shadow = qa.bitwiseAnd(1 << 4)   # bit 4 = cloud shadow
    return image.updateMask(cloud.eq(0).And(shadow.eq(0)))


def get_annual_nbr_image(year, aoi):
    """
    Build annual median NBR composite for one year, spatially clipped to aoi.

    Mission selection:
      2022+  → Landsat 9  (OLI-2, bands SR_B5=NIR, SR_B7=SWIR2)
      2013+  → Landsat 8  (OLI,   bands SR_B5=NIR, SR_B7=SWIR2)
      2000-12→ Landsat 5  (TM,    bands SR_B4=NIR, SR_B7=SWIR2)

    NBR is INVERTED (SWIR2 - NIR) so vegetation loss = positive delta,
    which is the convention required by LandTrendR.
    """
    date_start = f"{year}-{START_DAY}"
    date_end   = f"{year}-{END_DAY}"

    if year >= 2022:
        coll = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
                .filterDate(date_start, date_end)
                .filterBounds(aoi)
                .map(mask_landsat_c2)
                .select(["SR_B5", "SR_B7"], ["NIR", "SWIR2"]))
    elif year >= 2013:
        coll = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                .filterDate(date_start, date_end)
                .filterBounds(aoi)
                .map(mask_landsat_c2)
                .select(["SR_B5", "SR_B7"], ["NIR", "SWIR2"]))
    else:
        # Landsat 5 — preferred over L7 to avoid SLC-off data gaps (post-2003)
        coll = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
                .filterDate(date_start, date_end)
                .filterBounds(aoi)
                .map(mask_landsat_c2)
                .select(["SR_B4", "SR_B7"], ["NIR", "SWIR2"]))

    # Annual median composite → inverted NBR
    composite  = coll.median()
    nbr_inv    = composite.normalizedDifference(["SWIR2", "NIR"]).rename("NBR")

    return (nbr_inv
            .set("system:time_start", ee.Date.fromYMD(year, 7, 1).millis())
            .set("year", year))


def build_lt_collection(aoi):
    """
    Build annual NBR ImageCollection (START_YEAR to END_YEAR) for LandTrendR.
    One image per year; band ordering: NBR first (required by LT API).
    """
    images = [get_annual_nbr_image(y, aoi) for y in range(START_YEAR, END_YEAR + 1)]
    return ee.ImageCollection(images)


# ─── 6. LANDTRENDR EXTRACTION FUNCTION ────────────────────────────────────────
def extract_fire_nbr(row):
    """
    Run LandTrendR for one fire perimeter and return mean fitted NBR per year.

    LandTrendr band structure (axis 0):
      Row 0: year of observation
      Row 1: source (observed) spectral value
      Row 2: fitted (segmented) spectral value   ← we extract this row
      Row 3: is-vertex flag

    Parameters
    ----------
    row : GeoDataFrame row with fields event_id, ig_date, geometry

    Returns
    -------
    list of dicts: {event_id, fire_year, year, fitted_nbr}
    """
    fire_geom = ee.Geometry(row.geometry.__geo_interface__)
    fire_id   = row["event_id"]
    fire_year = int(row["ig_date"][:4])

    # Build annual collection constrained to this fire's extent
    lt_coll = build_lt_collection(fire_geom)

    # Run LandTrendR segmentation
    lt_result = ee.Algorithms.TemporalSegmentation.LandTrendr(
        **{**RUN_PARAMS, "timeSeries": lt_coll}
    )

    # Extract fitted values from LandTrendr band
    # arraySlice(0, 2, 3) → select row 2 (fitted values)
    # arrayProject([1])   → project to time axis (one value per year)
    # arrayFlatten        → convert to image with one band per year
    year_labels = [f"yr_{y}" for y in range(START_YEAR, END_YEAR + 1)]
    fitted_stack = (lt_result
                    .select("LandTrendr")
                    .arraySlice(0, 2, 3)
                    .arrayProject([1])
                    .arrayFlatten([year_labels]))

    # Mean fitted NBR over the fire polygon (synchronous call — suitable for EDA)
    mean_nbr = fitted_stack.reduceRegion(
        reducer   = ee.Reducer.mean(),
        geometry  = fire_geom,
        scale     = SCALE,
        maxPixels = int(1e8)
    ).getInfo()

    # Build per-year records
    records = []
    for y in range(START_YEAR, END_YEAR + 1):
        records.append({
            "event_id":   fire_id,
            "fire_year":  fire_year,
            "year":       y,
            "fitted_nbr": mean_nbr.get(f"yr_{y}", None)
        })
    return records


# ─── 7. RUN EXTRACTION LOOP ────────────────────────────────────────────────────
OUT_PATH = PROJ_ROOT / "output/biomass_timeseries.csv"

if OUT_PATH.exists():
    print(f"\nCSV already exists at {OUT_PATH} — delete it to re-run extraction.")
else:
    print(f"\nExtracting LandTrendR biomass for {len(sample)} fires...")
    print(f"  Scale: {SCALE}m | Years: {START_YEAR}–{END_YEAR} | Window: {START_DAY}–{END_DAY}")

    all_records = []
    errors      = []

    for i, (_, row) in enumerate(sample.iterrows(), 1):
        try:
            records = extract_fire_nbr(row)
            all_records.extend(records)
        except Exception as e:
            errors.append({"event_id": row["event_id"], "error": str(e)})

        if i % 25 == 0 or i == len(sample):
            print(f"  {i}/{len(sample)} fires processed...")

    df = pd.DataFrame(all_records)
    df.to_csv(OUT_PATH, index=False)

    print(f"\nSaved {len(df)} records ({len(sample)} fires × {END_YEAR - START_YEAR + 1} years)")
    print(f"Output: {OUT_PATH}")

    if errors:
        print(f"\nWARNING: {len(errors)} fires failed:")
        for e in errors[:5]:
            print(f"  {e}")

    # Sanity check on output
    df_check = pd.read_csv(OUT_PATH)
    assert set(["event_id", "fire_year", "year", "fitted_nbr"]).issubset(df_check.columns)
    assert df_check["year"].between(START_YEAR, END_YEAR).all()
    n_non_null = df_check["fitted_nbr"].notna().sum()
    pct_valid  = 100 * n_non_null / len(df_check)
    print(f"\nData quality: {n_non_null:,} / {len(df_check):,} records have valid NBR ({pct_valid:.1f}%)")
    if pct_valid < 50:
        print("WARNING: <50% valid — check date window, AOI, and Landsat availability.")
