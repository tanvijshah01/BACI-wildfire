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
library(tigris)
library(glue)
library(dplyr)
# FedData is no longer loaded — fetch_nlcd_landcover() below replaces its
# get_nlcd() call directly via httr/xml2 (both already FedData dependencies,
# so already installed) to skip a memory-heavy step; see that function.

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/05_prepare_forest_masks_west.R")

# terra decides whether to process a raster in memory or chunk it through
# disk based on its own estimate of "available" memory — which reads the
# NODE's full system RAM, not the ~4 GiB this job is actually capped to
# under GRIT's cgroup (confirmed 2026-09-20; see NOTES.md). That's why
# passing filename= to classify() alone didn't fix its OOM: filename= only
# says where the output goes, not how terra decides to compute it — terra
# still concluded a 954M-cell raster comfortably fits in memory and never
# switched to chunked processing. todisk = TRUE forces every terra
# operation for the rest of this script to always chunk through disk,
# overriding that heuristic outright rather than hoping any one call's
# arguments are enough to sidestep it.
terra::terraOptions(todisk = TRUE)

WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

# Override for a pilot run, e.g.: STATES_TO_RUN <- c("WY", "CO")
STATES_TO_RUN <- c("CA", "WY")

MASK_DIR   <- here("data", "processed", "forest_mask")
CTREES_DIR <- here("data", "processed", "ctrees")
dir.create(MASK_DIR, recursive = TRUE, showWarnings = FALSE)

FOREST_CLASSES <- c(41L, 42L, 43L)   # Deciduous, Evergreen, Mixed Forest
# (no reclassification matrix needed — the 30 m mask below reclassifies via
# a manual row-strip loop using FOREST_CLASSES directly, not terra::classify())

# A writeRaster() interrupted mid-write (e.g. laptop sleep — the same failure
# mode already seen with 00_crop_emapr_to_west.R) can leave a GeoTIFF with a
# correct header/extent/CRS but zero actual pixel data. file.exists() alone
# doesn't catch this, and these are 0/1 masks where every cell should be
# valid (never NA) — so any all-NA file is unambiguously corrupt. Caught this
# in practice on nlcd2004_forestfrac_30m_wy.tif (2026-08-12): correct 17066 x
# 20554 / EPSG:5070 header, 100% of cells NA.
raster_is_valid <- function(path) {
  # Wrapped in tryCatch() after a real crash on GRIT: a file killed mid-write
  # by the classify() OOM above (before the todisk=TRUE fix) was truncated
  # badly enough that GDAL couldn't even recognize it as a file at all
  # ("GDAL Error 4: ... not recognized as a supported file format"), which
  # makes terra::rast() throw a hard error rather than return an openable-
  # but-all-NA raster. Without this, that error propagated straight out of
  # raster_is_valid() and halted the whole script on its next run, instead
  # of being treated as "invalid, delete and rebuild" like this function is
  # supposed to do for any other kind of corruption.
  tryCatch(
    terra::global(terra::rast(path), "notNA")[[1]] > 0,
    error = function(e) FALSE
  )
}

# ── NLCD fetch, bypassing FedData::get_nlcd()'s factor/color-table step ──────
# FedData::get_nlcd() already crops server-side via MRLC's WCS endpoint (not a
# CONUS download — confirmed by reading FedData's own source,
# R/NLCD_FUNCTIONS.R on ropensci/FedData's GitHub), so the network fetch
# itself isn't what OOM-killed this step on GRIT (see NOTES.md's 2026-09-20
# entry). The step right after IS the suspect: for dataset == "Land_Cover",
# get_nlcd() converts the result to a categorical factor raster and attaches
# a full NLCD color table (terra::as.factor() + terra::coltab()) before
# writing it back out. We never use any of that — the very next thing this
# script does is reclassify raw class codes into a 0/1 mask via
# terra::classify() below, which reads identical underlying values whether or
# not a factor/color table is attached. So there's no behavior difference in
# skipping it, only a (hopefully) smaller peak memory footprint.
# Replicates get_nlcd()'s WCS request exactly (same URL pattern, same
# bbox-subset logic, same AEA projection for the subset coordinates) but
# returns the plain numeric raster straight from the WCS response instead.
fetch_nlcd_landcover <- function(template_sf, year = 2004, landmass = "L48") {
  coverage <- glue("NLCD_{year}_Land_Cover_{landmass}")
  source   <- glue("https://www.mrlc.gov/geoserver/mrlc_download/{coverage}/wcs")

  describe <- httr::GET(source, query = list(
    service = "WCS", version = "2.0.1",
    request = "DescribeCoverage", coverageid = coverage
  ))
  if (httr::status_code(describe) != 200L) {
    stop("No WCS coverage at ", source, " for NLCD ", year, " Land_Cover ", landmass)
  }

  xml_content <- describe |> httr::content(encoding = "UTF-8") |> xml2::as_list()
  envelope    <- xml_content$CoverageDescriptions$CoverageDescription$boundedBy$Envelope
  axis_labels <- envelope |> attr("axisLabels") |> strsplit(" ") |> unlist()

  # Same AEA projection FedData's WCS path uses, so the bbox subset lines up
  # with the service's own coordinate axes — template_sf can be any CRS.
  bbox <- template_sf |>
    sf::st_transform(
      "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    ) |>
    sf::st_bbox()

  tmp <- tempfile(fileext = ".tif")
  httr::GET(
    source,
    query = list(
      service    = "WCS", version = "2.0.1", request = "GetCoverage",
      coverageid = coverage,
      subset     = glue("{axis_labels[1]}({bbox['xmin']},{bbox['xmax']})"),
      subset     = glue("{axis_labels[2]}({bbox['ymin']},{bbox['ymax']})")
    ),
    httr::write_disk(tmp, overwrite = TRUE)
  )

  terra::rast(tmp)   # plain numeric raster, file-backed — no as.factor()/coltab()
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
    cat(glue("[REBUILD] {st} — {basename(out_30m)} exists but is corrupt (all-NA or unreadable); deleting and rebuilding.\n"))
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

    # fetch_nlcd_landcover() (defined above), not FedData::get_nlcd() directly
    # — see that function's comment for why: get_nlcd()'s post-download
    # as.factor()/coltab() step OOM-killed this on GRIT (NOTES.md 2026-09-20).
    nlcd_raw <- fetch_nlcd_landcover(state_5070, year = 2004)

    dl_elapsed <- round(proc.time()["elapsed"] - dl_t0)
    cat(glue("  Downloaded in {dl_elapsed}s | ",
             "res: {paste(round(terra::res(nlcd_raw)), collapse=' x ')} m | ",
             "cells: {scales::comma(terra::ncell(nlcd_raw))}\n"))

    cat(glue("[{st}] Reclassifying to 0/1 forest mask (30 m)...\n"))
    # terra::classify() OOM-killed this step on GRIT twice — once with no
    # filename= (holding the whole ~954M-cell result in memory), then again
    # with filename= AND terra::terraOptions(todisk = TRUE) both set. Neither
    # was enough: terra's own decision about how small to make its internal
    # processing chunks is driven by its estimate of "available" memory,
    # which reads the node's full system RAM rather than the ~4 GiB this job
    # is actually capped to (NOTES.md 2026-09-20) — so even in disk mode, the
    # chunks it chose were still too large. Rather than guess another terra
    # setting, this bypasses classify()'s automatic chunking entirely: a
    # manual row-strip read -> reclassify -> write loop via terra's own
    # low-level block I/O (writeStart/readValues/writeValues/writeStop),
    # sized explicitly ourselves instead of left to any heuristic — the same
    # technique already proven safe today in 04_download_ctrees_west.py's
    # Part A (STRIP_ROWS). `strip %in% FOREST_CLASSES` treats NA/unmatched
    # cells as 0, matching this script's documented 0/1-not-1/NA design
    # (see file header) exactly.
    n_rows_total <- terra::nrow(nlcd_raw)
    n_cols_total <- terra::ncol(nlcd_raw)
    strip_rows   <- max(1L, floor(5e6 / n_cols_total))  # ~5M cells/strip, any state width

    mask_30m <- terra::rast(nlcd_raw)   # blank template: same extent/res/CRS, no values
    terra::writeStart(mask_30m, out_30m, overwrite = FALSE, datatype = "INT1U",
                       gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
    for (row0 in seq(1L, n_rows_total, by = strip_rows)) {
      n_this <- min(strip_rows, n_rows_total - row0 + 1L)
      strip  <- terra::readValues(nlcd_raw, row = row0, nrows = n_this)
      strip  <- as.integer(strip %in% FOREST_CLASSES)
      terra::writeValues(mask_30m, strip, row0, n_this)
    }
    terra::writeStop(mask_30m)

    size_mb <- round(file.size(out_30m) / 1e6, 1)
    cat(glue("[{st}] Saved {basename(out_30m)} ({size_mb} MB)\n"))
    rm(nlcd_raw)
  }

  # ── 90 m mask — EPSG:5070, modal-aggregated (matches eMapR 100m TIF's ─────────
  # native resolution closely enough for 07's crop+resample step; mirrors old
  # 03_prepare_forest_mask.R §4, but 0/1-encoded like everything else here).
  if (file.exists(out_90m) && !raster_is_valid(out_90m)) {
    cat(glue("[REBUILD] {st} — {basename(out_90m)} exists but is corrupt (all-NA or unreadable); deleting and rebuilding.\n"))
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
    cat(glue("[REBUILD] {st} — {basename(out_100m)} exists but is corrupt (all-NA or unreadable); deleting and rebuilding.\n"))
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
