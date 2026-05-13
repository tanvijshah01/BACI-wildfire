# =============================================================================
# 01_extract_biomass_gee.py
# Extract annual mean NBR (biomass proxy) from Landsat / Google Earth Engine
# for MTBS wildfire perimeters, 1984-2023.
#
# SCRIPT OUTLINE
# 1.  Setup — imports, paths, study parameters
# 2.  GEE authentication and initialization
# 3.  Load and filter MTBS fire perimeters
# 4.  Sample fires for EDA (stratified by year)
# 5.  Landsat annual composite builder (harmonized across missions)
# 6.  Build GEE FeatureCollection of sampled fire polygons
# 7.  Extract annual mean NBR (reduceRegions per year — one GEE call per year)
# 8.  Save to output/biomass_timeseries.csv
#
# BEFORE RUNNING:
#   earthengine authenticate        (one-time browser login — stores project in credentials)
#
# OUTPUT:
#   output/biomass_timeseries.csv
#   Columns: event_id, fire_year, year, fitted_nbr
#   fitted_nbr is raw (unsmoothed) inverted NBR: higher = more vegetation loss.
#   Negate when displaying: display_nbr = -fitted_nbr
#
# NOTE: This EDA script uses raw annual Landsat composites rather than
# LandTrendR fitted values. LandTrendR can be enabled for the final analysis
# to suppress inter-annual noise.
# =============================================================================

# ─── 1. SETUP ─────────────────────────────────────────────────────────────────
import ee
import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# ── Resolve project root (two levels up from scripts/python/) ─────────────────
PROJ_ROOT = Path(__file__).resolve().parent.parent.parent

# ── Study parameters ──────────────────────────────────────────────────────────
START_YEAR  = 1984       # Full Landsat 5 TM archive — gives pre-fire baseline for all MTBS fires
END_YEAR    = 2023
START_DAY   = "07-01"   # July–September: peak summer growing season, Western US
END_DAY     = "09-30"
N_SAMPLE    = 100       # ~2-3 fires per year for CA-only EDA; increase for full analysis
STATES      = ["CA"]               # California pilot; expand to all 11 states for full analysis
SCALE       = 90        # pixel resolution for reduceRegion (90m faster, 30m precise)

# Western US approximate bounding box (pre-filter before spatial join in R)
WEST_LAT = (31.0, 49.0)
WEST_LON = (-124.5, -103.0)


# ─── 2. GEE AUTHENTICATION ────────────────────────────────────────────────────
GEE_PROJECT = "emlab-gcp"   # Google Cloud project linked to GEE account

print("Initializing Google Earth Engine...")
try:
    ee.Initialize(project=GEE_PROJECT)
    print(f"  GEE initialized (project: {GEE_PROJECT}).")
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
    .loc[mtbs_raw["event_id"].str[:2].isin(STATES)]   # restrict to EDA states
    .reset_index(drop=True)
)
print(f"  Fires after filter ({', '.join(STATES)}): {len(mtbs)}")

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


# ─── 5. LANDSAT ANNUAL COMPOSITE BUILDER ──────────────────────────────────────
def mask_landsat_c2(image):
    """Mask cloud and cloud shadow using QA_PIXEL band (Landsat Collection 2)."""
    qa     = image.select("QA_PIXEL")
    cloud  = qa.bitwiseAnd(1 << 3)   # bit 3 = cloud
    shadow = qa.bitwiseAnd(1 << 4)   # bit 4 = cloud shadow
    return image.updateMask(cloud.eq(0).And(shadow.eq(0)))


def get_annual_nbr_image(year, aoi):
    """
    Build annual median NBR composite for one year over aoi.

    Mission selection:
      2022+  → Landsat 9  (OLI-2, SR_B5=NIR, SR_B7=SWIR2)
      2013+  → Landsat 8  (OLI,   SR_B5=NIR, SR_B7=SWIR2)
      ≤2012  → Landsat 5  (TM,    SR_B4=NIR, SR_B7=SWIR2)
               Landsat 5 preferred over L7 to avoid SLC-off gaps (post-2003)

    NBR is INVERTED: (SWIR2 − NIR) / (SWIR2 + NIR) so that vegetation loss
    is positive, matching the LandTrendR convention for future use.
    Returns a masked image when no cloud-free scenes exist for this year/aoi.
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
        coll = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
                .filterDate(date_start, date_end)
                .filterBounds(aoi)
                .map(mask_landsat_c2)
                .select(["SR_B4", "SR_B7"], ["NIR", "SWIR2"]))

    # Merge a fully-masked placeholder so the collection always has the right
    # band names even when no real scenes exist (prevents band-not-found errors).
    placeholder = (ee.Image.constant([0, 0])
                   .rename(["NIR", "SWIR2"])
                   .updateMask(ee.Image.constant(0)))
    composite = coll.merge(ee.ImageCollection([placeholder])).median()

    # Inverted NBR: higher = more vegetation loss
    return composite.normalizedDifference(["SWIR2", "NIR"]).rename("NBR")


# ─── 6. BUILD GEE FEATURE COLLECTION FROM SAMPLED FIRES ──────────────────────
print(f"\nBuilding GEE FeatureCollection from {len(sample)} sampled fires...")
features = [
    ee.Feature(
        ee.Geometry(row.geometry.__geo_interface__),
        {"event_id": row["event_id"], "fire_year": int(row["ig_date"][:4])}
    )
    for _, row in sample.iterrows()
]
fires_fc = ee.FeatureCollection(features)
fires_aoi = fires_fc.geometry().bounds()   # bounding box of fires — avoids self-touching edge errors in exact union


# ─── 7. EXTRACT ANNUAL MEAN NBR (reduceRegions per year) ──────────────────────
# One GEE call per year processes all fires at once — 40 calls total
# instead of 40 × N_fires calls in the per-fire loop approach.
OUT_PATH = PROJ_ROOT / "output/biomass_timeseries.csv"

if OUT_PATH.exists():
    print(f"\nCSV already exists at {OUT_PATH} — delete it to re-run extraction.")
else:
    n_years = END_YEAR - START_YEAR + 1
    print(f"\nExtracting annual NBR for {len(sample)} fires × {n_years} years...")
    print(f"  Scale: {SCALE}m | Years: {START_YEAR}–{END_YEAR} | Window: {START_DAY}–{END_DAY}")
    print(f"  Method: raw annual Landsat composites (EDA mode — no LandTrendR smoothing)")

    all_records = []
    errors      = []

    for y in range(START_YEAR, END_YEAR + 1):
        try:
            nbr_img = get_annual_nbr_image(y, fires_aoi)

            # reduceRegions: mean NBR over each fire polygon in one server-side call
            result = nbr_img.reduceRegions(
                collection = fires_fc,
                reducer    = ee.Reducer.mean(),
                scale      = SCALE
            ).getInfo()

            for feat in result["features"]:
                props = feat["properties"]
                all_records.append({
                    "event_id":   props.get("event_id"),
                    "fire_year":  props.get("fire_year"),
                    "year":       y,
                    "fitted_nbr": props.get("mean")  # reduceRegions with ee.Reducer.mean() uses "mean" as the key
                })
        except Exception as e:
            errors.append({"year": y, "error": str(e)})

        n_done = y - START_YEAR + 1
        if n_done % 5 == 0 or y == END_YEAR:
            print(f"  {n_done}/{n_years} years processed ({y})...")

    # ─── 8. SAVE AND VALIDATE ─────────────────────────────────────────────────
    df = pd.DataFrame(all_records)
    df.to_csv(OUT_PATH, index=False)

    print(f"\nSaved {len(df)} records ({len(sample)} fires × {n_years} years)")
    print(f"Output: {OUT_PATH}")

    if errors:
        print(f"\nWARNING: {len(errors)} years failed:")
        for e in errors[:5]:
            print(f"  {e}")

    df_check = pd.read_csv(OUT_PATH)
    assert set(["event_id", "fire_year", "year", "fitted_nbr"]).issubset(df_check.columns)
    assert df_check["year"].between(START_YEAR, END_YEAR).all(), \
        f"Years outside {START_YEAR}–{END_YEAR}"
    n_non_null = df_check["fitted_nbr"].notna().sum()
    pct_valid  = 100 * n_non_null / len(df_check)
    print(f"\nData quality: {n_non_null:,} / {len(df_check):,} records have valid NBR ({pct_valid:.1f}%)")
    if pct_valid < 50:
        print("WARNING: <50% valid — check date window, AOI, and Landsat availability.")
