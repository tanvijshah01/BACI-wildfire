# =============================================================================
# 03_prepare_forest_mask.R
#
# Download NLCD 2001 for California and save two binary forest masks:
#   - 30 m  EPSG:5070  (matches eMapR native resolution and CRS)
#   - ~100 m EPSG:4326 (matches ctrees native resolution and CRS)
#
# Forest classes retained: 41 Deciduous, 42 Evergreen, 43 Mixed Forest
# All other classes set to NA.
#
# Run ONCE before running 02_extract_emapr_within_fires.R (forested) or
# 04_extract_ctrees_within_fires.R.
#
# Output:
#   data/processed/forest_mask/nlcd2001_forest_30m_ca.tif
#   data/processed/forest_mask/nlcd2001_forest_100m_ca.tif
#
# OUTLINE
# 1. Setup
# 2. Download NLCD 2001 clipped to CA
# 3. Reclassify to binary forest mask (30 m, EPSG:5070)
# 4. Project and resample to ctrees grid (~100 m, EPSG:4326)
# 5. Sanity checks and report
# =============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(here)
library(FedData)
library(tigris)
library(glue)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/03_prepare_forest_mask.R")

MASK_DIR   <- here("data", "processed", "forest_mask")
CTREES_DIR <- here("data", "processed", "ctrees")
EMAPR_DIR  <- here("data", "processed", "emapr_biomass_ca")

OUT_30M  <- file.path(MASK_DIR, "nlcd2001_forest_30m_ca.tif")
OUT_100M <- file.path(MASK_DIR, "nlcd2001_forest_100m_ca.tif")

dir.create(MASK_DIR, recursive = TRUE, showWarnings = FALSE)

FOREST_CLASSES <- c(41L, 42L, 43L)   # Deciduous, Evergreen, Mixed Forest

# ── 2. Download NLCD 2001 ─────────────────────────────────────────────────────
cat("Loading CA boundary...\n")
ca_sf   <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS == "CA")
ca_5070 <- sf::st_transform(ca_sf, 5070)

cat("Downloading NLCD 2001 for CA (this may take 1-2 minutes)...\n")
t0 <- proc.time()["elapsed"]

# FedData returns a SpatRaster; template must be an sf or SpatVector in any CRS
nlcd_raw <- FedData::get_nlcd(
  template  = ca_5070,
  label     = "CA",
  year      = 2001,
  dataset   = "landcover",
  extraction.dir = file.path(MASK_DIR, "nlcd_raw")
)

cat(glue("  Downloaded in {round(proc.time()['elapsed'] - t0)}s\n"))
cat(glue("  NLCD CRS:        {terra::crs(nlcd_raw, describe=TRUE)$code}\n"))
cat(glue("  NLCD resolution: {paste(round(terra::res(nlcd_raw)), collapse=' x ')} m\n"))
cat(glue("  NLCD extent:     {paste(round(as.vector(terra::ext(nlcd_raw))), collapse=' ')}\n"))

# ── 3. Binary forest mask — 30 m, EPSG:5070 ──────────────────────────────────
if (file.exists(OUT_30M)) {
  cat("[SKIP] 30 m mask already exists.\n")
} else {
  cat("Building 30 m forest mask...\n")

  # Reclassify: forest classes → 1, everything else → NA
  rcl <- matrix(
    c(FOREST_CLASSES, rep(1L, length(FOREST_CLASSES))),
    ncol = 2
  )
  mask_30m <- terra::classify(nlcd_raw, rcl, others = NA)
  mask_30m <- terra::clamp(mask_30m, lower = 1, upper = 1, values = TRUE)

  terra::writeRaster(mask_30m, OUT_30M, overwrite = FALSE,
                     datatype = "INT1U",
                     gdal = c("COMPRESS=LZW", "TILED=YES",
                               "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
  size_mb <- round(file.size(OUT_30M) / 1e6, 1)
  cat(glue("  Saved {OUT_30M} ({size_mb} MB)\n"))
}

mask_30m <- terra::rast(OUT_30M)

# ── 4. ~100 m mask — EPSG:4326, aligned to ctrees grid ───────────────────────
if (file.exists(OUT_100M)) {
  cat("[SKIP] ~100 m mask already exists.\n")
} else {
  cat("Building ~100 m forest mask (EPSG:4326)...\n")

  # Use a ctrees TIF as the target grid template
  ctrees_tif <- list.files(CTREES_DIR, pattern = "^ctrees_\\d{4}_ca_100m\\.tif$",
                            full.names = TRUE)[1]
  stopifnot("No ctrees 100m TIF found — run scripts/python/03_download_ctrees_ca.py first" =
              !is.null(ctrees_tif) && length(ctrees_tif) > 0)
  ctrees_template <- terra::rast(ctrees_tif)
  cat(glue("  ctrees template: {basename(ctrees_tif)}\n"))
  cat(glue("  ctrees CRS: {terra::crs(ctrees_template, describe=TRUE)$code}",
           "  res: {paste(round(terra::res(ctrees_template), 6), collapse=' x ')} deg\n"))

  # Project 30m mask to ctrees CRS/grid; use "near" to keep binary values clean
  mask_100m <- terra::project(mask_30m, ctrees_template, method = "near")

  terra::writeRaster(mask_100m, OUT_100M, overwrite = FALSE,
                     datatype = "INT1U",
                     gdal = c("COMPRESS=LZW", "TILED=YES",
                               "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
  size_mb <- round(file.size(OUT_100M) / 1e6, 1)
  cat(glue("  Saved {OUT_100M} ({size_mb} MB)\n"))
}

mask_100m <- terra::rast(OUT_100M)

# ── 5. Sanity checks ──────────────────────────────────────────────────────────
cat("\n=== Forest Mask Report ===\n")

for (lst in list(list(r = mask_30m,  label = "30 m  (EPSG:5070)"),
                 list(r = mask_100m, label = "~100 m (EPSG:4326)"))) {
  n_forest <- terra::global(lst$r, "notNA")[[1]]
  n_total  <- terra::ncell(lst$r)
  pct      <- round(100 * n_forest / n_total, 1)
  cat(glue("  {lst$label}: {scales::comma(n_forest)} forest px / {scales::comma(n_total)} total  ({pct}%)\n"))
}

cat("\nForest mask ready. Next steps:\n")
cat("  Rscript scripts/r/02_extract_emapr_within_fires.R   (forested eMapR extraction)\n")
cat("  Rscript scripts/r/04_extract_ctrees_within_fires.R  (forested ctrees extraction)\n")
