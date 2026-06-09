# =============================================================================
# 01_precompute_emapr_ca.R
# Pre-compute summary stats and 300 m display rasters for the EDA study years.
# Run once after 00_crop_emapr_to_ca.R has produced the CA-clipped files.
#
# Input:  data/processed/emapr_biomass_ca/composite_YYYY_ca.tif
# Output: data/processed/emapr_biomass_ca/stats_YYYY.rds    (scalar stats + sample)
#         data/processed/emapr_biomass_ca/composite_YYYY_ca_300m.tif
#         data/processed/emapr_biomass_ca/composite_YYYY_ca_100m.tif
#
# Skips any year whose output files already exist.
# =============================================================================
# 1. Setup
# 2. California boundary
# 3. Loop over study years: compute stats → save RDS; aggregate → save 300m and 100m TIFs
# 4. Report
# =============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(here)
library(glue)
library(tigris)
library(scales)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("analysis/03_emapr_biomass_exploration.qmd")

STUDY_YEARS <- c(2000, 2001, 2002, 2003)
AGG_FACTOR  <- 10
SAMPLE_N    <- 200000
CA_DIR      <- here("data", "processed", "emapr_biomass_ca")

# ── 2. California boundary for masking checks ─────────────────────────────────
ca_wgs84 <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS == "CA")
stopifnot("CA boundary not found" = nrow(ca_wgs84) == 1)
ca_vect <- terra::vect(sf::st_transform(ca_wgs84, crs = 5070))

# ── 3. Process each study year ────────────────────────────────────────────────
for (yr in STUDY_YEARS) {
  in_path      <- file.path(CA_DIR, glue("composite_{yr}_ca.tif"))
  rds_path     <- file.path(CA_DIR, glue("stats_{yr}.rds"))
  map_path_300 <- file.path(CA_DIR, glue("composite_{yr}_ca_300m.tif"))
  map_path_100 <- file.path(CA_DIR, glue("composite_{yr}_ca_100m.tif"))

  if (!file.exists(in_path)) stop(glue("CA file missing for {yr} — run 00_crop_emapr_to_ca.R first"))

  # ── Stats ──────────────────────────────────────────────────────────────────
  if (file.exists(rds_path)) {
    cat(glue("[SKIP] stats_{yr}.rds already exists\n"))
  } else {
    cat(glue("[STAT] {yr} — computing stats..."))
    t0 <- proc.time()["elapsed"]
    r  <- terra::rast(in_path)
    gs      <- terra::global(r, fun = c("min", "max", "mean", "sd"), na.rm = TRUE)
    n_total <- terra::ncell(r)
    n_valid <- terra::global(r, fun = "notNA")[[1]]
    vals    <- terra::spatSample(r, size = SAMPLE_N, na.rm = TRUE)[[1]]
    vals    <- vals[vals > 0]
    qs      <- quantile(vals, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    stats   <- list(year = yr, mean = gs$mean, sd = gs$sd,
                    n_valid = n_valid, n_total = n_total,
                    median = qs[3], q95 = qs[5], max = gs$max, vals = vals)
    saveRDS(stats, rds_path)
    cat(glue(" done in {round(proc.time()['elapsed'] - t0)}s\n"))
  }

  # ── 300 m aggregate ────────────────────────────────────────────────────────
  if (file.exists(map_path_300)) {
    cat(glue("[SKIP] composite_{yr}_ca_300m.tif already exists\n"))
  } else {
    cat(glue("[AGG]  {yr} — aggregating to 300 m..."))
    t0  <- proc.time()["elapsed"]
    r   <- terra::rast(in_path)
    lyr <- terra::aggregate(r, fact = AGG_FACTOR, fun = "mean", na.rm = TRUE)
    terra::writeRaster(lyr, map_path_300, overwrite = FALSE, gdal = c("COMPRESS=LZW"))
    size_mb <- round(file.size(map_path_300) / 1e6, 1)
    cat(glue(" done — {size_mb} MB in {round(proc.time()['elapsed'] - t0)}s\n"))
  }

  # ── ~100 m aggregate (fact = 3: 30 m × 3 = 90 m, labeled ~100 m) ──────────
  if (file.exists(map_path_100)) {
    cat(glue("[SKIP] composite_{yr}_ca_100m.tif already exists\n"))
  } else {
    cat(glue("[AGG]  {yr} — aggregating to ~100 m (fact = 3)..."))
    t0  <- proc.time()["elapsed"]
    r   <- terra::rast(in_path)
    lyr <- terra::aggregate(r, fact = 3, fun = "mean", na.rm = TRUE)
    terra::writeRaster(lyr, map_path_100, overwrite = FALSE,
                       gdal = c("COMPRESS=LZW", "TILED=YES",
                                "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
    size_mb <- round(file.size(map_path_100) / 1e6, 1)
    cat(glue(" done — {size_mb} MB in {round(proc.time()['elapsed'] - t0)}s\n"))
  }
}

# ── 4. Report ─────────────────────────────────────────────────────────────────
cat("\n=== Pre-compute complete ===\n")
for (yr in STUDY_YEARS) {
  rds  <- file.path(CA_DIR, glue("stats_{yr}.rds"))
  t300 <- file.path(CA_DIR, glue("composite_{yr}_ca_300m.tif"))
  t100 <- file.path(CA_DIR, glue("composite_{yr}_ca_100m.tif"))
  cat(glue("  {yr}: stats={if(file.exists(rds)) 'OK' else 'MISSING'}",
           "  300m={if(file.exists(t300)) 'OK' else 'MISSING'}",
           "  100m={if(file.exists(t100)) 'OK' else 'MISSING'}\n"))
}
