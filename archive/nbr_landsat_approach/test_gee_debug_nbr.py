"""
Quick diagnostic: run reduceRegions for one year on one known CA fire polygon
and print ALL returned properties to see what key GEE uses for the NBR band.
"""
import ee
from pathlib import Path
import geopandas as gpd

PROJ_ROOT = Path(__file__).resolve().parent.parent.parent
GEE_PROJECT = "emlab-gcp"

ee.Initialize(project=GEE_PROJECT)
print("GEE initialized.")

# ── Load one California fire (first fire in MTBS after filtering) ──────────────
mtbs = gpd.read_file(PROJ_ROOT / "data/raw/mtbs/mtbs_perimeter_data/mtbs_perims_DD.shp")
mtbs["year"] = mtbs["ig_date"].str[:4].astype(int)
ca = (mtbs
      .loc[mtbs["incid_type"] == "Wildfire"]
      .loc[mtbs["year"] == 2005]
      .loc[mtbs["event_id"].str[:2] == "CA"]
      .head(1))

print(f"Test fire: {ca['event_id'].values[0]}, year 2005")
print(f"Geometry type: {ca.geometry.values[0].geom_type}")

# ── Build GEE feature ──────────────────────────────────────────────────────────
fire_geom = ee.Geometry(ca.geometry.values[0].__geo_interface__)
fire_feat = ee.Feature(fire_geom, {"event_id": ca["event_id"].values[0]})
fc = ee.FeatureCollection([fire_feat])

# ── Build NBR image for 2005 ───────────────────────────────────────────────────
aoi = fc.geometry().bounds()

coll = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
        .filterDate("2005-07-01", "2005-09-30")
        .filterBounds(aoi))

print(f"Landsat scenes found: {coll.size().getInfo()}")

# Cloud mask
def mask_c2(img):
    qa = img.select("QA_PIXEL")
    cloud  = qa.bitwiseAnd(1 << 3)
    shadow = qa.bitwiseAnd(1 << 4)
    return img.updateMask(cloud.eq(0).And(shadow.eq(0)))

composite = coll.map(mask_c2).median().select(["SR_B4", "SR_B7"], ["NIR", "SWIR2"])

# Check if composite has any valid pixels
pixel_count = composite.select("NIR").reduceRegion(
    reducer=ee.Reducer.count(),
    geometry=fire_geom,
    scale=90,
    maxPixels=1e8
).getInfo()
print(f"Valid NIR pixels in fire polygon: {pixel_count}")

# Compute NBR and check band name
nbr = composite.normalizedDifference(["SWIR2", "NIR"]).rename("NBR")
print(f"NBR band names: {nbr.bandNames().getInfo()}")

# Run reduceRegions and print ALL properties
result = nbr.reduceRegions(
    collection=fc,
    reducer=ee.Reducer.mean(),
    scale=90
).getInfo()

print(f"\nFeature properties returned:")
for feat in result["features"]:
    print(dict(feat["properties"]))
