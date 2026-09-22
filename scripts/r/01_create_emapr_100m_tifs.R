# =============================================================================
# 01_create_emapr_100m_tifs.R
#
# One-time preprocessing: downsample CA eMapR composites from 30 m to ~90 m
# (labeled "~100 m" by convention, matching ctrees native ~100 m resolution).
#
# Why: Enables fair comparison with ctrees at matched ~100 m resolution.
# 30 m × 3 = 90 m is the closest clean block-average to ctrees ~100 m (~0.000889°).
# The resulting ~90 m files (~55 MB each) are also memory-safe for terra::extract().
#
# Source: prefers the retired CA-only 30 m crop (composite_YYYY_ca.tif, from
# the retired 00_crop_emapr_to_ca.R) when it's already on disk (laptop-era
# years); falls back to cropping + masking the CURRENT West-wide crop
# (composite_YYYY_west.tif, from 00_crop_emapr_to_west.R) down to CA when
# the CA-only file isn't available — e.g. on GRIT, which never ran the
# retired CA-only crop script, so years like 2005-2010 only exist as West
# crops there. Added 2026-09-21 so this script works from the current
# pipeline's output without depending on the retired script having run.
#
# How: terra::aggregate(fact = 3, fun = "mean") shrinks each 30 m pixel grid
# by 3× in each dimension → ~90 m output. Skip-safe: already-existing 100 m
# TIFs are never overwritten. Both eMapR crops share EPSG:5070 (confirmed in
# 00_crop_emapr_to_west.R), so the West-crop fallback needs crop()+mask()
# only, no reprojection.
#
# Run once from any working directory:
#   Rscript scripts/r/01_create_emapr_100m_tifs.R
#
# Outputs: data/processed/emapr_biomass_ca/composite_YYYY_ca_100m.tif
#          (one per available source year, from either source)
#
# OUTLINE
# 1. Setup + CA boundary (only needed for the West-crop fallback path)
# 2. Locate source years (CA-only 30 m crop, or West-wide crop as fallback)
# 3. Identify years missing 100 m versions
# 4. Aggregate and write
# 5. Report summary
# =============================================================================

library(terra)
library(sf)
library(tigris)
library(here)
library(glue)
library(dplyr)

here::i_am("scripts/r/01_create_emapr_100m_tifs.R")
sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)

# Must be set before any raster I/O touching the West-wide composites below
# — GDAL's block cache otherwise defaults to a % of the node's full system
# RAM, not the ~4 GiB cgroup cap this job actually runs under on GRIT (see
# scripts/r/05_prepare_forest_masks_west.R for the full writeup of this
# failure mode and why it's set this early, before any raster is touched).
terra::setGDALconfig("GDAL_CACHEMAX", "64")                        # MB
terra::setGDALconfig("GDAL_MAX_DATASET_POOL_RAM_USAGE", "64")      # MB

# terra decides whether to hold a raster in memory or chunk it through disk
# based on its own estimate of "available" memory, which reads the NODE's
# full system RAM, not this job's actual ~4 GiB cgroup cap — the same
# failure already diagnosed and fixed in 05_prepare_forest_masks_west.R.
# Missed here on the first attempt (2026-09-21): the West-crop fallback
# path's crop() of a ~954M-cell CA-sized extent out of the full West raster
# was OOM-killed immediately, exactly like 05's classify() was before this
# was added there. todisk = TRUE forces every terra operation for the rest
# of this script to chunk through disk instead of trusting that heuristic.
terra::terraOptions(todisk = TRUE)

# ── 1. Setup + CA boundary ────────────────────────────────────────────────────
EMAPR_CA_DIR   <- here("data", "processed", "emapr_biomass_ca")
EMAPR_WEST_DIR <- here("data", "processed", "emapr_biomass_west")
dir.create(EMAPR_CA_DIR, recursive = TRUE, showWarnings = FALSE)

# Restrict to the years biomass_within_fires.qmd's current params actually
# need (study_year_min/max = 2005/2010) rather than building all 34
# available West-crop years unattended — narrows both the OOM blast radius
# and the wasted time if something still goes wrong partway through.
# Set to NULL to build every available year instead.
YEARS_TO_BUILD <- 2005:2010

# Only used for the West-crop fallback path (crop + mask to CA); skipped
# entirely if every needed year already has a CA-only 30 m source on disk.
ca_5070 <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS == "CA") |>
  sf::st_transform(5070)
ca_vect_5070 <- terra::vect(ca_5070)

# ── 2. Locate source years (either source) ────────────────────────────────────
ca_files <- list.files(EMAPR_CA_DIR, pattern = "^composite_\\d{4}_ca\\.tif$",
                       full.names = TRUE)
ca_years <- as.integer(regmatches(basename(ca_files), regexpr("\\d{4}", basename(ca_files))))

west_files <- list.files(EMAPR_WEST_DIR, pattern = "^composite_\\d{4}_west\\.tif$",
                         full.names = TRUE)
west_years <- as.integer(regmatches(basename(west_files), regexpr("\\d{4}", basename(west_files))))

all_years <- sort(union(ca_years, west_years))
cat("CA-only 30 m crops found:  ", length(ca_years), "year(s) —",
    paste(ca_years, collapse = ", "), "\n")
cat("West-wide 30 m crops found:", length(west_years), "year(s) —",
    paste(west_years, collapse = ", "), "\n")
cat("Years available overall:   ", length(all_years), "\n\n")

if (!is.null(YEARS_TO_BUILD)) {
  all_years <- intersect(all_years, YEARS_TO_BUILD)
  cat("Restricting to YEARS_TO_BUILD:", paste(all_years, collapse = ", "), "\n\n")
}

# ── 3. Identify years missing 100 m TIFs ──────────────────────────────────────
out_files  <- file.path(EMAPR_CA_DIR, glue("composite_{all_years}_ca_100m.tif"))
need_build <- !file.exists(out_files)

cat(sum(!need_build), "100 m TIF(s) already exist — will skip.\n")
cat(sum(need_build),  "100 m TIF(s) to build:",
    paste(all_years[need_build], collapse = ", "), "\n\n")

if (!any(need_build)) {
  cat("Nothing to do — all 100 m TIFs are present.\n")
  quit(save = "no", status = 0)
}

# ── 4. Aggregate and write ────────────────────────────────────────────────────
# fact = 3: 30 m × 3 = 90 m (~100 m by convention); fun = "mean" preserves
# mean AGB within block. NAflag keeps existing nodata value from source raster.
t_start <- proc.time()
counter <- 0L
n_build <- sum(need_build)

for (i in which(need_build)) {
  counter <- counter + 1L
  yr      <- all_years[i]
  tif_out <- out_files[i]

  cat(glue("[{counter}/{n_build}] {yr} ... "))

  if (yr %in% ca_years) {
    # Preferred: already CA-shaped, no crop/mask needed.
    r_in <- terra::rast(ca_files[match(yr, ca_years)])
  } else {
    # Fallback: crop + mask the West-wide crop down to CA first. Masking
    # the full West raster directly (125M+ cells across 11 states) with no
    # crop() first is the whole-raster pattern that OOMs/takes 500s+
    # elsewhere in this project (see CLAUDE.md -> Avoid); cropping first
    # keeps this cheap by only reading CA's window.
    r_in <- terra::rast(west_files[match(yr, west_years)])
    r_in <- terra::crop(r_in, ca_vect_5070)
    r_in <- terra::mask(r_in, ca_vect_5070)
  }

  r_out <- terra::aggregate(r_in, fact = 3, fun = "mean", na.rm = TRUE)

  terra::writeRaster(r_out, tif_out,
                     filetype  = "GTiff",
                     overwrite = FALSE,
                     gdal      = c("COMPRESS=LZW", "TILED=YES",
                                   "BLOCKXSIZE=512", "BLOCKYSIZE=512"))

  rm(r_in, r_out)
  gc(verbose = FALSE)

  elapsed <- round((proc.time() - t_start)["elapsed"], 1)
  cat(glue("done  [{elapsed}s elapsed]\n"))
}

# ── 5. Summary ────────────────────────────────────────────────────────────────
all_100m <- sort(list.files(EMAPR_CA_DIR,
                            pattern = "^composite_\\d{4}_ca_100m\\.tif$"))
cat("\nDone. 100 m TIFs now available for",
    length(all_100m), "year(s):\n")
cat(paste(all_100m, collapse = "\n"), "\n")
