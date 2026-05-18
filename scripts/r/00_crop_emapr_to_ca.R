# =============================================================================
# 00_crop_emapr_to_ca.R
# One-time preprocessing: crop and mask each CONUS-wide eMapR composite to
# California and save as a compressed GeoTIFF.
#
# Input:  data/raw/emapr_biomass/composite_YYYY_median.tif  (~27.7 GB each)
# Output: data/processed/emapr_biomass_ca/composite_YYYY_ca.tif  (~100–300 MB)
#
# Run once; subsequent calls skip files that already exist.
# =============================================================================
# 1. Load packages
# 2. Build CA boundary (EPSG:5070)
# 3. Discover available years from raw directory
# 4. Loop: crop → mask → write (skip if output exists)
# 5. Report summary
# =============================================================================

# ── 1. Packages ───────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(here)
library(glue)
library(tigris)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("analysis/03_emapr_biomass_exploration.qmd")

# ── 2. California boundary in EPSG:5070 ───────────────────────────────────────
ca_wgs84    <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS == "CA")
stopifnot("CA boundary not found" = nrow(ca_wgs84) == 1)

# Project to NAD83 / Conus Albers — same CRS as eMapR composites
ca_sf   <- sf::st_transform(ca_wgs84, crs = 5070)
ca_vect <- terra::vect(ca_sf)
ca_bbox <- sf::st_bbox(ca_sf)
ca_ext  <- terra::ext(ca_bbox["xmin"], ca_bbox["xmax"],
                      ca_bbox["ymin"], ca_bbox["ymax"])

cat("CA extent (EPSG:5070):", paste(round(as.vector(ca_ext), 0), collapse = ", "), "\n")

# ── 3. Discover available years ───────────────────────────────────────────────
RAW_DIR  <- here("data", "raw", "emapr_biomass")
raw_files <- list.files(RAW_DIR, pattern = "^composite_\\d{4}_median\\.tif$", full.names = TRUE)
avail_years <- as.integer(sub(".*composite_(\\d{4})_median\\.tif", "\\1", raw_files))
avail_years <- sort(avail_years)

cat(length(avail_years), "raw composite files found:", paste(avail_years, collapse = ", "), "\n\n")

# ── 4. Crop → mask → write ────────────────────────────────────────────────────
OUT_DIR <- here("data", "processed", "emapr_biomass_ca")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- data.frame(year = integer(), status = character(), size_mb = numeric(),
                      stringsAsFactors = FALSE)

for (yr in avail_years) {
  in_path  <- file.path(RAW_DIR, glue("composite_{yr}_median.tif"))
  out_path <- file.path(OUT_DIR, glue("composite_{yr}_ca.tif"))

  if (file.exists(out_path)) {
    size_mb <- round(file.size(out_path) / 1e6, 1)
    cat(glue("[SKIP] {yr} — already exists ({size_mb} MB)\n"))
    results <- rbind(results, data.frame(year = yr, status = "skipped", size_mb = size_mb))
    next
  }

  cat(glue("[PROC] {yr} — cropping and masking..."))
  t0 <- proc.time()["elapsed"]

  r   <- terra::rast(in_path)
  r   <- terra::crop(r, ca_ext)
  r   <- terra::mask(r, ca_vect)
  terra::writeRaster(r, out_path, overwrite = FALSE, gdal = c("COMPRESS=LZW"))

  elapsed <- round(proc.time()["elapsed"] - t0, 0)
  size_mb <- round(file.size(out_path) / 1e6, 1)
  cat(glue(" done — {size_mb} MB in {elapsed}s\n"))

  results <- rbind(results, data.frame(year = yr, status = "saved", size_mb = size_mb))
}

# ── 5. Summary ────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat("Saved: ", sum(results$status == "saved"),  "files\n")
cat("Skipped:", sum(results$status == "skipped"), "files\n")
cat("Output directory:", OUT_DIR, "\n")
print(results)
