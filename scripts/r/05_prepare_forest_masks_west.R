# =============================================================================
# 05_prepare_forest_masks_west.R
#
# Download NLCD 2004 for each of the 11 Western study states and save binary
# (0/1, not 1/NA) forest masks per state, at three resolutions/grids —
# this is the full generalized replacement for 03_prepare_forest_mask.R
# (which was CA-only and 1/NA-encoded), now the SOLE forest-mask source for
# analysis/biomass_within_fires.qmd and the 06/07/08 extraction scripts.
#
# 0/1 (not 1/NA) rasters throughout: non-forest cells are 0, not NA, so a
# single crop -> mask(maskvalues=0) -> mean(na.rm=TRUE) per fire polygon (in
# 06/07/08) directly yields fraction-forest or forest-only AGB, with no
# separate "count total polygon pixels" step needed.
#
# Skip-safe per state AND per resolution: re-running only builds outputs that
# are missing. For a pilot run, edit STATES_TO_RUN below (e.g. c("WY", "CO"))
# before running the rest of WESTERN_STATES later — already-built states are
# skipped.
#
# Forest classes retained: 41 Deciduous, 42 Evergreen, 43 Mixed Forest.
#
# Output (st = lowercase state code):
#   data/processed/forest_mask/nlcd2004_forestfrac_30m_<st>.tif   30 m,  EPSG:5070 (native; 06 pct-forest extraction)
#   data/processed/forest_mask/nlcd2004_forestfrac_90m_<st>.tif   90 m,  EPSG:5070, modal-aggregated (07 eMapR extraction)
#   data/processed/forest_mask/nlcd2004_forestfrac_100m_<st>.tif  ~100 m, EPSG:4326, projected onto the ctrees
#                                                                  grid template (08 ctrees extraction + display maps)
#
# The ~100 m variant requires a ctrees_<year>_<st>_100m.tif to already exist
# as a template (currently only CA) — skipped with a message for states
# without one yet.
#
# OUTLINE
# 1. Setup
# 2. Per-state loop: download NLCD 2004, build 30 m / 90 m / ~100 m 0/1 masks
# 3. Sanity checks and report
# =============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(here)
library(FedData)
library(tigris)
library(glue)
library(dplyr)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/05_prepare_forest_masks_west.R")

WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

# Override for a pilot run, e.g.: STATES_TO_RUN <- c("WY", "CO")
STATES_TO_RUN <- c("CA", "WY")

MASK_DIR   <- here("data", "processed", "forest_mask")
CTREES_DIR <- here("data", "processed", "ctrees")
dir.create(MASK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(MASK_DIR, "nlcd_raw"), recursive = TRUE, showWarnings = FALSE)

FOREST_CLASSES <- c(41L, 42L, 43L)   # Deciduous, Evergreen, Mixed Forest
rcl <- matrix(c(FOREST_CLASSES, rep(1L, length(FOREST_CLASSES))), ncol = 2)

# A writeRaster() interrupted mid-write (e.g. laptop sleep — the same failure
# mode already seen with 00_crop_emapr_to_west.R) can leave a GeoTIFF with a
# correct header/extent/CRS but zero actual pixel data. file.exists() alone
# doesn't catch this, and these are 0/1 masks where every cell should be
# valid (never NA) — so any all-NA file is unambiguously corrupt. Caught this
# in practice on nlcd2004_forestfrac_30m_wy.tif (2026-08-12): correct 17066 x
# 20554 / EPSG:5070 header, 100% of cells NA.
raster_is_valid <- function(path) {
  terra::global(terra::rast(path), "notNA")[[1]] > 0
}

cat("States to build:", paste(STATES_TO_RUN, collapse = ", "), "\n\n")

# ── 2. Per-state loop ──────────────────────────────────────────────────────────
t_all <- proc.time()["elapsed"]

for (st in STATES_TO_RUN) {
  t0 <- proc.time()["elapsed"]
  out_30m  <- file.path(MASK_DIR, glue("nlcd2004_forestfrac_30m_{tolower(st)}.tif"))
  out_90m  <- file.path(MASK_DIR, glue("nlcd2004_forestfrac_90m_{tolower(st)}.tif"))
  out_100m <- file.path(MASK_DIR, glue("nlcd2004_forestfrac_100m_{tolower(st)}.tif"))

  ctrees_tif <- list.files(CTREES_DIR, pattern = glue("^ctrees_\\d{{4}}_{tolower(st)}_100m\\.tif$"),
                            full.names = TRUE)[1]
  need_100m  <- !file.exists(out_100m) && !is.na(ctrees_tif)

  ok_30m  <- file.exists(out_30m)  && raster_is_valid(out_30m)
  ok_90m  <- file.exists(out_90m)  && raster_is_valid(out_90m)
  ok_100m <- file.exists(out_100m) && raster_is_valid(out_100m)

  if (ok_30m && ok_90m && (ok_100m || is.na(ctrees_tif))) {
    cat(glue("[SKIP] {st} — 30 m/90 m/100 m masks already built (or 100 m unavailable, no ctrees template).\n"))
    next
  }

  # ── 30 m mask — download (or reuse cached NLCD) + reclassify ────────────────
  if (file.exists(out_30m) && !ok_30m) {
    cat(glue("[REBUILD] {st} — {basename(out_30m)} exists but is corrupt (all-NA); deleting and rebuilding.\n"))
    file.remove(out_30m)
  }
  if (ok_30m) {
    cat(glue("[SKIP] {st} — {basename(out_30m)} already exists.\n"))
    mask_30m <- terra::rast(out_30m)
  } else {
    cat(glue("[{st}] Loading state boundary...\n"))
    state_sf   <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
      dplyr::filter(STUSPS == st)
    stopifnot("State boundary not found" = nrow(state_sf) == 1)
    state_5070 <- sf::st_transform(state_sf, 5070)

    cat(glue("[{st}] Downloading NLCD 2004...\n"))
    dl_t0 <- proc.time()["elapsed"]

    nlcd_raw <- FedData::get_nlcd(
      template       = state_5070,
      label          = st,
      year           = 2004,
      dataset        = "landcover",
      extraction.dir = file.path(MASK_DIR, "nlcd_raw")
    )

    dl_elapsed <- round(proc.time()["elapsed"] - dl_t0)
    cat(glue("  Downloaded in {dl_elapsed}s | ",
             "res: {paste(round(terra::res(nlcd_raw)), collapse=' x ')} m | ",
             "cells: {scales::comma(terra::ncell(nlcd_raw))}\n"))

    cat(glue("[{st}] Reclassifying to 0/1 forest mask (30 m)...\n"))
    mask_30m <- terra::classify(nlcd_raw, rcl, others = 0L)

    terra::writeRaster(mask_30m, out_30m, overwrite = FALSE,
                       datatype = "INT1U",
                       gdal = c("COMPRESS=LZW", "TILED=YES",
                                 "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
    size_mb <- round(file.size(out_30m) / 1e6, 1)
    cat(glue("[{st}] Saved {basename(out_30m)} ({size_mb} MB)\n"))
    rm(nlcd_raw)
  }

  # ── 90 m mask — EPSG:5070, modal-aggregated (matches eMapR 100m TIF's ─────────
  # native resolution closely enough for 07's crop+resample step; mirrors old
  # 03_prepare_forest_mask.R §4, but 0/1-encoded like everything else here).
  if (file.exists(out_90m) && !raster_is_valid(out_90m)) {
    cat(glue("[REBUILD] {st} — {basename(out_90m)} exists but is corrupt (all-NA); deleting and rebuilding.\n"))
    file.remove(out_90m)
  }
  if (file.exists(out_90m)) {
    cat(glue("[SKIP] {st} — {basename(out_90m)} already exists.\n"))
  } else {
    cat(glue("[{st}] Building 90 m forest mask (EPSG:5070, fact=3, modal)...\n"))
    mask_90m <- terra::aggregate(mask_30m, fact = 3, fun = "modal")
    terra::writeRaster(mask_90m, out_90m, overwrite = FALSE,
                       datatype = "INT1U",
                       gdal = c("COMPRESS=LZW", "TILED=YES",
                                 "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
    size_mb <- round(file.size(out_90m) / 1e6, 1)
    cat(glue("[{st}] Saved {basename(out_90m)} ({size_mb} MB)\n"))
    rm(mask_90m)
  }

  # ── ~100 m mask — EPSG:4326, projected onto the ctrees grid template ─────────
  # Mirrors old 03_prepare_forest_mask.R §5. Requires a ctrees 100m TIF for
  # this state as the target grid — skip (not fail) if none exists yet, since
  # ctrees downloads currently only cover CA.
  if (file.exists(out_100m) && !raster_is_valid(out_100m)) {
    cat(glue("[REBUILD] {st} — {basename(out_100m)} exists but is corrupt (all-NA); deleting and rebuilding.\n"))
    file.remove(out_100m)
  }
  if (file.exists(out_100m)) {
    cat(glue("[SKIP] {st} — {basename(out_100m)} already exists.\n"))
  } else if (is.na(ctrees_tif)) {
    cat(glue("[{st}] No ctrees 100m TIF found — skipping ~100 m mask ",
             "(needed for ctrees extraction + display maps, not for 06's pct-forest sweep).\n"))
  } else {
    cat(glue("[{st}] Building ~100 m forest mask (EPSG:4326, projected onto ctrees grid)...\n"))
    ctrees_template <- terra::rast(ctrees_tif)
    mask_100m <- terra::project(mask_30m, ctrees_template, method = "near")
    terra::writeRaster(mask_100m, out_100m, overwrite = FALSE,
                       datatype = "INT1U",
                       gdal = c("COMPRESS=LZW", "TILED=YES",
                                 "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
    size_mb <- round(file.size(out_100m) / 1e6, 1)
    cat(glue("[{st}] Saved {basename(out_100m)} ({size_mb} MB)\n"))
    rm(mask_100m, ctrees_template)
  }

  total_elapsed <- round(proc.time()["elapsed"] - t0)
  cat(glue("[{st}] Done — {total_elapsed}s total\n\n"))

  rm(mask_30m)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE, full = TRUE)
}

# ── 3. Sanity checks and report ────────────────────────────────────────────────
report_res <- function(label, pattern) {
  cat(glue("\n=== Forest Mask Report (0/1, {label}) ===\n"))
  for (st in STATES_TO_RUN) {
    out_tif <- file.path(MASK_DIR, glue(pattern))
    if (!file.exists(out_tif)) {
      cat(glue("  {st}: MISSING\n"))
      next
    }
    r        <- terra::rast(out_tif)
    n_total  <- terra::ncell(r)
    n_forest <- terra::global(r, "sum", na.rm = TRUE)[[1]]
    pct      <- round(100 * n_forest / n_total, 1)
    cat(glue("  {st}: {scales::comma(n_forest)} forest px / {scales::comma(n_total)} total  ({pct}%)\n"))
  }
}

report_res("30 m, EPSG:5070",   "nlcd2004_forestfrac_30m_{tolower(st)}.tif")
report_res("90 m, EPSG:5070",   "nlcd2004_forestfrac_90m_{tolower(st)}.tif")
report_res("~100 m, EPSG:4326", "nlcd2004_forestfrac_100m_{tolower(st)}.tif")

cat(glue("\nTotal elapsed: {round((proc.time()['elapsed'] - t_all) / 60, 1)} min\n"))
cat("\nForest masks ready. Next steps:\n")
cat("  Rscript scripts/r/06_extract_pct_forest_within_fires.R\n")
cat("  Rscript scripts/r/07_extract_emapr_within_fires_new.R\n")
cat("  Rscript scripts/r/08_extract_ctrees_within_fires_new.R\n")
