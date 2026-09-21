# =============================================================================
# 06_extract_pct_forest_within_fires.R
#
# Extract % forest cover within each MTBS fire perimeter (Western US, 2000-2023)
# using the per-state 0/1 NLCD 2004 forest masks built by
# 05_prepare_forest_masks_west.R.
#
# Run ONCE after 05_prepare_forest_masks_west.R has produced the masks for the
# states you want. Run from the project root before rendering
# analysis/mtbs_assessment_comparison.qmd.
#
# Why a separate script: same reason as 02_extract_emapr_within_fires.R and
# 04_extract_ctrees_within_fires.R — terra::extract()/heavy per-polygon loops
# can appear frozen inside Quarto on Windows (Quarto buffers chunk output while
# terra's C++ threading runs).
#
# Mirrors the state-attribution/dedup logic in analysis/mtbs_assessment_comparison.qmd
# §3 exactly (event_id prefix as authoritative state, not the raw spatial join)
# so this script's output lines up 1:1 with the qmd's own mtbs_study fires.
#
# Output:
#   data/processed/forest_mask/pct_forest_by_fire_west.csv
#   columns: event_id, STUSPS, fire_year, burnbndac, asmnt_binary, pct_forest, n_pixels
#
# OUTLINE
# 1. Setup
# 2. Load and filter MTBS to the Western US study cohort (matches the qmd)
# 3. Check cache — determine which states still need extraction
# 4. Per-state loop: crop + mask + mean per fire polygon
# 5. Report
# =============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(here)
library(dplyr)
library(readr)
library(glue)
library(tigris)
library(scales)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/06_extract_pct_forest_within_fires.R")

# Must be set before ANY raster I/O — GDAL fixes its block-cache size on
# first use, and changing this afterward has no effect. Confirmed on GRIT:
# peak RSS grew slowly but steadily across the per-fire loop below (from
# ~2098 MB at fire 770 to ~2150 MB at fire 930, then killed) even with
# periodic gc()/tmpFiles() every 10 fires — since plain R gc() can't touch
# it, this is GDAL's own C-level block cache, not R-managed memory. Its
# default size is a % of the NODE's full system RAM (the same wrong basis
# behind every other terra/GDAL surprise on GRIT today), which grows
# unbounded as each fire's crop() touches a new, essentially random region
# of the 954M-cell mask file — the opposite of 05's sequential row-strip
# reads, which never needed this because they touch each block exactly
# once, in order. Capped to an explicit small size instead.
terra::setGDALconfig("GDAL_CACHEMAX", "64")                        # MB
terra::setGDALconfig("GDAL_MAX_DATASET_POOL_RAM_USAGE", "64")      # MB

WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")
YEAR_MIN <- 2000
YEAR_MAX <- 2023

# Diagnostic only — periodic gc()/tmpFiles() (added after 06 was confirmed
# OOM-killed on GRIT partway through CA's 1044-fire loop) didn't fix it, and
# with no progress logging at all there's no way to tell whether it died on
# fire 1 or fire 1000 — meaning that fix might just never have gotten a
# chance to run. Mirrors log_peak_memory() in 04_download_ctrees_west.py
# (same /proc/self/status VmHWM technique) so the next failure is
# diagnosable from the log instead of another blind guess.
log_peak_memory <- function(label) {
  status_path <- "/proc/self/status"
  if (!file.exists(status_path)) return(invisible(NULL))
  hwm <- grep("^VmHWM:", readLines(status_path, warn = FALSE), value = TRUE)
  if (length(hwm) > 0) cat(glue("    [mem] {label}: peak RSS = {trimws(sub('VmHWM:', '', hwm))}\n"))
  invisible(NULL)
}

# Override for a pilot run, e.g.: STATES_TO_RUN <- c("WY", "CO")
STATES_TO_RUN <- c("CA")

MTBS_PATH   <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
MASK_DIR    <- here("data", "processed", "forest_mask")
OUT_CSV     <- here("data", "processed", "forest_mask", "pct_forest_by_fire_west.csv")

# ── 2. Load and filter MTBS to match analysis/mtbs_assessment_comparison.qmd §3 ─
# Filtered AT READ TIME via an OGR SQL query — see 08_extract_ctrees_within_
# fires_new.R's identical block for why (loading the full ~30k-fire national
# shapefile then subsetting in R hit GRIT's cgroup memory cap directly). No
# burnbndac threshold here (unlike 07/08) — this script's own filter chain
# below doesn't apply one, only incid_type + a year range, so only
# incid_type is safe to push down without changing behavior.
cat("Loading MTBS...\n")
mtbs_raw <- sf::st_read(
  MTBS_PATH, quiet = TRUE,
  query = "SELECT * FROM mtbs_perims_DD WHERE incid_type = 'Wildfire'"
)

n_invalid <- sum(!sf::st_is_valid(mtbs_raw))
if (n_invalid > 0) mtbs_raw <- sf::st_make_valid(mtbs_raw)

mtbs_parsed <- mtbs_raw |>
  dplyr::mutate(
    year = as.integer(substr(ig_date, 1, 4)),
    asmnt_binary = dplyr::case_when(
      asmnt_type %in% c("Extended", "Extended (SS)") ~ "Extended",
      TRUE                                            ~ "Initial"
    )
  )

mtbs_wf <- mtbs_parsed |>
  dplyr::filter(incid_type == "Wildfire", year >= YEAR_MIN, year <= YEAR_MAX)
cat("After Wildfire + year filter:", nrow(mtbs_wf), "\n")

western_states_sf <- tigris::states(cb = TRUE, progress_bar = FALSE) |>
  dplyr::filter(STUSPS %in% WESTERN_STATES) |>
  sf::st_transform(sf::st_crs(mtbs_raw))

mtbs_joined <- sf::st_join(mtbs_wf, western_states_sf["STUSPS"], join = sf::st_intersects, left = FALSE)

# Same dedup fix as the qmd: MTBS's own event_id state prefix is authoritative,
# not the (possibly duplicated) spatial-join match.
mtbs_study <- mtbs_joined |>
  dplyr::mutate(STUSPS = substr(event_id, 1, 2)) |>
  dplyr::filter(STUSPS %in% WESTERN_STATES) |>
  dplyr::distinct(event_id, .keep_all = TRUE)

stopifnot("Duplicate event_ids after dedup" = !any(duplicated(mtbs_study$event_id)))
cat("Study fires (deduped, all 11 states):", nrow(mtbs_study), "\n\n")

mtbs_5070 <- sf::st_transform(mtbs_study, 5070)

# ── 3. Cache check — which fires still need extraction (per-fire, not per-state) ──
# 06 was repeatedly OOM-killed partway through CA's 1044-fire loop (see
# NOTES.md 2026-09-20); with results only ever written once, at the very
# end of the whole per-state loop, every failed run lost ALL progress and
# started over from fire 1. Checkpointing per fire now (see the extraction
# loop below) needs the resume check to match: checking whether a STATE has
# any row at all (the original logic) would either wrongly treat a
# half-done state as fully complete, or force redoing already-checkpointed
# fires — checks individual event_ids instead.
done_fires <- character(0)
if (file.exists(OUT_CSV)) {
  partial    <- readr::read_csv(OUT_CSV, show_col_types = FALSE)
  done_fires <- unique(partial$event_id)
  cat("Partial cache found —", length(done_fires), "fire(s) already extracted.\n\n")
}

states_to_do <- STATES_TO_RUN[
  vapply(STATES_TO_RUN, function(st) {
    st_fires <- mtbs_5070$event_id[mtbs_5070$STUSPS == st]
    length(st_fires) > 0 && !all(st_fires %in% done_fires)
  }, logical(1))
]
if (length(states_to_do) == 0) {
  cat("Cache is up to date for all requested states. Nothing to do.\n")
  cat("Delete", basename(OUT_CSV), "to force re-extraction.\n")
  quit(save = "no")
}
cat("States to extract this run:", paste(states_to_do, collapse = ", "), "\n\n")

# ── 4. Per-state loop: crop + mask + mean per fire polygon ─────────────────────
terra::terraOptions(threads = 1, progress = 0)
t_all <- proc.time()["elapsed"]

for (st in states_to_do) {
  mask_tif <- file.path(MASK_DIR, glue("nlcd2004_forestfrac_30m_{tolower(st)}.tif"))
  if (!file.exists(mask_tif)) {
    cat(glue("[{st}] Mask not found ({basename(mask_tif)}) — ",
             "run scripts/r/05_prepare_forest_masks_west.R first. Skipping.\n\n"))
    next
  }

  fires_st_all <- mtbs_5070 |> dplyr::filter(STUSPS == st)
  fires_st     <- fires_st_all |> dplyr::filter(!event_id %in% done_fires)
  n_fires      <- nrow(fires_st)
  n_skipped    <- nrow(fires_st_all) - n_fires
  if (n_fires == 0) {
    cat(glue("[{st}] No study fires. Skipping.\n\n"))
    next
  }

  cat(glue("[{st}] Extracting % forest for {n_fires} fire polygon(s)",
           "{if (n_skipped > 0) glue(' ({n_skipped} already checkpointed, skipped)') else ''}...\n"))
  t0 <- proc.time()["elapsed"]

  forest_mask <- terra::rast(mask_tif)
  vect_st     <- terra::vect(fires_st)

  # Polygon-by-polygon: crop the mask to each polygon's extent, mask, mean.
  # Each crop reads only a few disk blocks (~KB RAM peak vs whole-raster ops) —
  # the same mitigation used in 02/04_extract_*_within_fires.R for Windows
  # terra memory limits. Values are 0/1, so mean(na.rm=TRUE) directly gives
  # fraction-forest — no separate denominator step needed.
  #
  # Checkpointed every 10 fires (written straight to OUT_CSV, not held until
  # the whole loop finishes) — this loop has been killed four times on GRIT
  # partway through, and every time, ALL progress was lost because nothing
  # was written until the very end. This turns "start over from fire 1
  # every failed attempt" into "resume from wherever it died" — the same
  # per-fire/per-year checkpointing pattern already used in
  # 04_download_ctrees_west.py and 07/08's per-(state,year) resumability,
  # just at per-fire granularity here since this script processes a whole
  # state's fires in one uninterrupted pass. Also still true regardless of
  # whether the underlying memory growth (NOTES.md 2026-09-20) ever gets
  # fully fixed — checkpointing alone guarantees this eventually completes
  # across enough re-runs, even if any single run can't get through all of
  # them.
  n_zero      <- 0L
  batch_start <- 1L
  flush_batch <- function(j) {
    idx <- batch_start:j
    df_batch <- data.frame(
      event_id     = fires_st$event_id[idx],
      STUSPS       = st,
      fire_year    = as.integer(fires_st$year[idx]),
      burnbndac    = fires_st$burnbndac[idx],
      asmnt_binary = fires_st$asmnt_binary[idx],
      pct_forest   = pct_forest[idx],
      n_pixels     = n_px[idx]
    )
    readr::write_csv(df_batch, OUT_CSV, append = file.exists(OUT_CSV))
    batch_start <<- j + 1L
    n_zero      <<- n_zero + sum(df_batch$n_pixels == 0L)
  }

  pct_forest <- numeric(n_fires)
  n_px       <- integer(n_fires)
  for (j in seq_len(n_fires)) {
    poly     <- vect_st[j]
    fm_c     <- terra::crop(forest_mask, poly, snap = "out")
    fm_m     <- terra::mask(fm_c, poly)
    vals     <- terra::values(fm_m, na.rm = TRUE)
    n_px[j]       <- length(vals)
    pct_forest[j] <- if (length(vals) > 0L) 100 * mean(vals) else NA_real_

    if (j %% 10 == 0 || j == n_fires) {
      flush_batch(j)
      terra::tmpFiles(remove = TRUE)
      gc(verbose = FALSE, full = TRUE)
      cat(glue("    ...{j}/{n_fires} fires processed (checkpointed)\n"))
      log_peak_memory(glue("after fire {j}"))
    }
  }

  elapsed <- round(proc.time()["elapsed"] - t0, 1)
  cat(glue("[{st}] Done — {elapsed}s ({round(elapsed/n_fires, 2)}s/fire) — ",
           "{n_zero} fire(s) with 0 masked pixels\n\n"))

  rm(forest_mask, vect_st)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE, full = TRUE)
}

# ── 5. Report ───────────────────────────────────────────────────────────────────
final         <- readr::read_csv(OUT_CSV, show_col_types = FALSE)
total_elapsed <- round((proc.time()["elapsed"] - t_all) / 60, 1)
cat(glue("\nDone. {nrow(final)} rows in {basename(OUT_CSV)}",
         " ({length(unique(final$STUSPS))} states) — {total_elapsed} min this run.\n"))
cat("You can now render analysis/mtbs_assessment_comparison.qmd.\n")
