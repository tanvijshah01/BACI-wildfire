# =============================================================================
# 02_explore_ctrees_zarr.py
# Explore the structure of the ucsb-emlab/BACI-wildfires arraylake zarr store.
# Run this FIRST to discover variable names, dimensions, resolution, and
# temporal range before writing the download script.
#
# SCRIPT OUTLINE
# 1. Connect to arraylake and open root zarr group
# 2. Walk all groups and arrays — print name, shape, dtype, chunks
# 3. Print coordinate arrays and their ranges
# 4. Print variable attributes (units, long_name, fill values)
# 5. Estimate CA subset data volume
# 6. Open with xarray and print Dataset repr
#
# BEFORE RUNNING:
#   pip install arraylake zarr xarray  (already done)
#   Log in once: python -c "from arraylake import Client; Client()"
#
# OUTPUT: Console only — no files written.
# =============================================================================

# ─── 1. CONNECT TO ARRAYLAKE ──────────────────────────────────────────────────
from arraylake import Client
import zarr
import numpy as np
import sys

REPO_NAME = "ucsb-emlab/BACI-wildfires"
BRANCH    = "main"

# CA bounding box for size estimation
CA_LAT = (32.5, 42.0)
CA_LON = (-124.5, -114.1)

print("=" * 60)
print(f"Connecting to arraylake repo: {REPO_NAME}")
print("=" * 60)

try:
    client  = Client()
    repo    = client.get_repo(REPO_NAME)
    session = repo.readonly_session(branch=BRANCH)
    root    = zarr.open_group(session.store, zarr_format=3, mode="r")
    print(f"Connected successfully (branch: {BRANCH})\n")
except Exception as e:
    print(f"\nERROR connecting to arraylake: {e}")
    print("Fix: make sure you are logged in — run `earthmover login` or follow")
    print("the browser prompt that appears when Client() is first called.")
    sys.exit(1)


# ─── 2. WALK ALL GROUPS AND ARRAYS ────────────────────────────────────────────
def walk_zarr(node, prefix=""):
    """Recursively print all groups and arrays with their metadata."""
    for key in node.keys():
        path = f"{prefix}/{key}" if prefix else key
        item = node[key]
        if isinstance(item, zarr.Group):
            print(f"  GROUP  {path}/")
            walk_zarr(item, prefix=path)
        else:
            # Array
            chunks = item.chunks if hasattr(item, "chunks") else "unknown"
            print(f"  ARRAY  {path}")
            print(f"           shape  = {item.shape}")
            print(f"           dtype  = {item.dtype}")
            print(f"           chunks = {chunks}")

print("─── Zarr store structure ───────────────────────────────────")
walk_zarr(root)
print()


# ─── 3. COORDINATE ARRAYS AND RANGES ─────────────────────────────────────────
# Look for likely coordinate names
COORD_CANDIDATES = ["time", "year", "lat", "latitude", "lon", "longitude",
                    "x", "y", "band", "level"]

print("─── Coordinate arrays ──────────────────────────────────────")
for name in COORD_CANDIDATES:
    # Search top-level and one level deep
    for node, prefix in [(root, ""), *[(root[g], g) for g in root.keys()
                                        if isinstance(root[g], zarr.Group)]]:
        if name in node:
            arr = node[name]
            if hasattr(arr, "shape") and arr.ndim == 1 and arr.size > 0:
                vals = arr[:]
                label = f"{prefix}/{name}" if prefix else name
                print(f"  {label}: shape={arr.shape}, dtype={arr.dtype}")
                print(f"    min={vals.min()}, max={vals.max()}")
                if arr.size <= 10:
                    print(f"    values={list(vals)}")
                else:
                    print(f"    first 5: {list(vals[:5])}")
                    print(f"    last  5: {list(vals[-5:])}")
print()


# ─── 4. VARIABLE ATTRIBUTES ───────────────────────────────────────────────────
print("─── Variable attributes ────────────────────────────────────")
def print_attrs(node, prefix=""):
    for key in node.keys():
        path = f"{prefix}/{key}" if prefix else key
        item = node[key]
        if isinstance(item, zarr.Group):
            print_attrs(item, prefix=path)
        else:
            attrs = dict(item.attrs) if hasattr(item, "attrs") else {}
            if attrs:
                print(f"  {path}:")
                for k, v in attrs.items():
                    print(f"    {k}: {v}")

print_attrs(root)
print()


# ─── 5. CA SUBSET SIZE ESTIMATE ───────────────────────────────────────────────
print("─── CA subset size estimate ────────────────────────────────")

# Find first data array that is 2D or 3D (likely the biomass variable)
def find_data_arrays(node, prefix=""):
    results = []
    for key in node.keys():
        path = f"{prefix}/{key}" if prefix else key
        item = node[key]
        if isinstance(item, zarr.Group):
            results.extend(find_data_arrays(item, prefix=path))
        elif hasattr(item, "ndim") and item.ndim >= 2:
            results.append((path, item))
    return results

data_arrays = find_data_arrays(root)
for path, arr in data_arrays[:3]:  # show first 3 candidates
    nbytes_total = arr.nbytes if hasattr(arr, "nbytes") else arr.size * arr.dtype.itemsize
    print(f"  {path}: total uncompressed = {nbytes_total / 1e9:.2f} GB")
    if arr.ndim == 3:
        # Assume (time, lat, lon) — estimate CA fraction
        lat_dim, lon_dim = arr.shape[-2], arr.shape[-1]
        ca_lat_frac = (CA_LAT[1] - CA_LAT[0]) / 180
        ca_lon_frac = (CA_LON[1] - CA_LON[0]) / 360
        ca_frac = ca_lat_frac * ca_lon_frac * 6  # rough multiplier for non-global
        ca_mb = (nbytes_total * ca_frac) / 1e6
        print(f"    CA subset estimate (rough): {ca_mb:.0f} MB")
print()


# ─── 6. XARRAY REPR ──────────────────────────────────────────────────────────
print("─── xarray Dataset repr ────────────────────────────────────")
try:
    import xarray as xr
    ds = xr.open_dataset(
        session.store,
        engine="zarr",
        zarr_format=3,
        consolidated=False
    )
    print(ds)
    print()
    print("Data variables:")
    for var in ds.data_vars:
        print(f"  {var}: {ds[var].dims} {ds[var].shape} | units={ds[var].attrs.get('units','?')}")
except Exception as e:
    print(f"xarray open failed: {e}")
    print("(This is OK — the zarr structure above is sufficient.)")

print()
print("=" * 60)
print("Exploration complete. Share this output to proceed with")
print("writing the download script (03_download_ctrees_ca.py).")
print("=" * 60)
