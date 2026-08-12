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

WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")
YEAR_MIN <- 2000
YEAR_MAX <- 2023

# Override for a pilot run, e.g.: STATES_TO_RUN <- c("WY", "CO")
STATES_TO_RUN <- c("CA")

MTBS_PATH   <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
MASK_DIR    <- here("data", "processed", "forest_mask")
OUT_CSV     <- here("data", "processed", "forest_mask", "pct_forest_by_fire_west.csv")

# ── 2. Load and filter MTBS to match analysis/mtbs_assessment_comparison.qmd §3 ─
cat("Loading MTBS...\n")
mtbs_raw <- sf::st_read(MTBS_PATH, quiet = TRUE)

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

# ── 3. Cache check — which states still need extraction ────────────────────────
done_states <- character(0)
if (file.exists(OUT_CSV)) {
  partial     <- readr::read_csv(OUT_CSV, show_col_types = FALSE)
  done_states <- sort(unique(partial$STUSPS))
  cat("Partial cache found —", length(done_states), "state(s) already done:",
      paste(done_states, collapse = ", "), "\n\n")
}

states_to_do <- setdiff(STATES_TO_RUN, done_states)
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

  fires_st <- mtbs_5070 |> dplyr::filter(STUSPS == st)
  n_fires  <- nrow(fires_st)
  if (n_fires == 0) {
    cat(glue("[{st}] No study fires. Skipping.\n\n"))
    next
  }

  cat(glue("[{st}] Extracting % forest for {n_fires} fire polygon(s)...\n"))
  t0 <- proc.time()["elapsed"]

  forest_mask <- terra::rast(mask_tif)
  vect_st     <- terra::vect(fires_st)

  # Polygon-by-polygon: crop the mask to each polygon's extent, mask, mean.
  # Each crop reads only a few disk blocks (~KB RAM peak vs whole-raster ops) —
  # the same mitigation used in 02/04_extract_*_within_fires.R for Windows
  # terra memory limits. Values are 0/1, so mean(na.rm=TRUE) directly gives
  # fraction-forest — no separate denominator step needed.
  pct_forest <- numeric(n_fires)
  n_px       <- integer(n_fires)
  for (j in seq_len(n_fires)) {
    poly     <- vect_st[j]
    fm_c     <- terra::crop(forest_mask, poly, snap = "out")
    fm_m     <- terra::mask(fm_c, poly)
    vals     <- terra::values(fm_m, na.rm = TRUE)
    n_px[j]       <- length(vals)
    pct_forest[j] <- if (length(vals) > 0L) 100 * mean(vals) else NA_real_
  }

  df_st <- data.frame(
    event_id      = fires_st$event_id,
    STUSPS        = st,
    fire_year     = as.integer(fires_st$year),
    burnbndac     = fires_st$burnbndac,
    asmnt_binary  = fires_st$asmnt_binary,
    pct_forest    = pct_forest,
    n_pixels      = n_px
  )

  readr::write_csv(df_st, OUT_CSV, append = file.exists(OUT_CSV))

  elapsed <- round(proc.time()["elapsed"] - t0, 1)
  n_zero  <- sum(df_st$n_pixels == 0L)
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
