# =============================================================================
# 00_crop_emapr_to_west.R
# One-time preprocessing: crop and mask each CONUS-wide eMapR composite to
# the union of all 11 Western study states and save as a compressed GeoTIFF.
# Generalizes 00_crop_emapr_to_ca.R (kept as-is) from CA-only to the full
# study region.
#
# Input:  data/raw/emapr_biomass/composite_YYYY_median.tif  (~27.7 GB each)
# Output: data/processed/emapr_biomass_west/composite_YYYY_west.tif  (~1 GB)
#
# Only years whose raw file is already present locally are processed — years
# missing from data/raw/emapr_biomass/ are reported as skipped (see
# DATA_DOWNLOAD_GUIDE.md for the rclone flow that fetches a missing year to a
# scratch path, then this script picks it up on the next run).
#
# Run once; subsequent calls skip files that already exist.
# =============================================================================
# 1. Load packages
# 2. Build 11-state West boundary (EPSG:5070)
# 3. Discover available years from raw directory
# 4. Loop: crop → mask → write (skip if output exists or input missing)
# 5. Report summary
# =============================================================================

# ── 1. Packages ───────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(here)
library(glue)
library(tigris)
library(dplyr)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/00_crop_emapr_to_west.R")

WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

# ── 2. West boundary in EPSG:5070 (union of all 11 states) ────────────────────
west_wgs84 <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS %in% WESTERN_STATES)
stopifnot("Expected 11 Western states" = nrow(west_wgs84) == length(WESTERN_STATES))

# Project to NAD83 / Conus Albers — same CRS as eMapR composites
west_sf   <- sf::st_transform(west_wgs84, crs = 5070)
west_vect <- terra::vect(west_sf)
west_bbox <- sf::st_bbox(west_sf)
west_ext  <- terra::ext(west_bbox["xmin"], west_bbox["xmax"],
                        west_bbox["ymin"], west_bbox["ymax"])

cat("West extent (EPSG:5070):", paste(round(as.vector(west_ext), 0), collapse = ", "), "\n")

# A writeRaster() interrupted mid-write (laptop sleep, background-task kill —
# a recurring failure mode on this machine, see NOTES.md 2026-08-12/13 and
# scripts/r/05_prepare_forest_masks_west.R) can leave a GeoTIFF with a correct
# header/extent/CRS but no actual pixel data. file.exists() alone doesn't
# catch this. Mirrors 05's raster_is_valid() — an all-NA output here means the
# crop/mask/write was never really completed, so it isn't safe to skip.
raster_is_valid <- function(path) {
  terra::global(terra::rast(path), "notNA")[[1]] > 0
}

# ── 3. Discover available years ───────────────────────────────────────────────
RAW_DIR  <- here("data", "raw", "emapr_biomass")
raw_files <- list.files(RAW_DIR, pattern = "^composite_\\d{4}_median\\.tif$", full.names = TRUE)
avail_years <- as.integer(sub(".*composite_(\\d{4})_median\\.tif", "\\1", raw_files))
avail_years <- sort(avail_years)

cat(length(avail_years), "raw composite files found locally:", paste(avail_years, collapse = ", "), "\n\n")

# ── 4. Crop → mask → write ────────────────────────────────────────────────────
OUT_DIR <- here("data", "processed", "emapr_biomass_west")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- data.frame(year = integer(), status = character(), size_mb = numeric(),
                      stringsAsFactors = FALSE)

for (yr in avail_years) {
  in_path  <- file.path(RAW_DIR, glue("composite_{yr}_median.tif"))
  out_path <- file.path(OUT_DIR, glue("composite_{yr}_west.tif"))

  if (file.exists(out_path)) {
    if (raster_is_valid(out_path)) {
      size_mb <- round(file.size(out_path) / 1e6, 1)
      cat(glue("[SKIP] {yr} — already exists ({size_mb} MB)\n"))
      results <- rbind(results, data.frame(year = yr, status = "skipped", size_mb = size_mb))
      next
    } else {
      cat(glue("[REBUILD] {yr} — existing file is corrupt (all-NA, interrupted write), deleting and redoing...\n"))
      file.remove(out_path)
    }
  }

  cat(glue("[PROC] {yr} — cropping and masking..."))
  t0 <- proc.time()["elapsed"]

  r   <- terra::rast(in_path)
  r   <- terra::crop(r, west_ext)
  r   <- terra::mask(r, west_vect)
  terra::writeRaster(r, out_path, overwrite = FALSE, gdal = c("COMPRESS=LZW"))

  elapsed <- round(proc.time()["elapsed"] - t0, 0)
  size_mb <- round(file.size(out_path) / 1e6, 1)
  cat(glue(" done — {size_mb} MB in {elapsed}s\n"))

  results <- rbind(results, data.frame(year = yr, status = "saved", size_mb = size_mb))
}

# Report years missing raw input entirely (not fetched by this script)
missing_years <- setdiff(1990:2023, avail_years)
if (length(missing_years) > 0) {
  cat(glue("\n[MISSING RAW INPUT] {length(missing_years)} years not found in {RAW_DIR}: ",
           "{paste(missing_years, collapse = ', ')}\n"))
  cat("See DATA_DOWNLOAD_GUIDE.md for the rclone flow to fetch these without\n")
  cat("permanently storing the full CONUS file locally.\n")
}

# ── 5. Summary ────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat("Saved: ", sum(results$status == "saved"),  "files\n")
cat("Skipped:", sum(results$status == "skipped"), "files\n")
cat("Output directory:", OUT_DIR, "\n")
print(results)
