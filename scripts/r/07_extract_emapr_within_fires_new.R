# =============================================================================
# 07_extract_emapr_within_fires_new.R
#
# Extracts mean eMapR AGB within MTBS fire perimeters, forested pixels only,
# for a user-specified set of Western states (STATES_TO_RUN below). Feeds
# analysis/biomass_within_fires.qmd.
#
# Generalizes the CA-only version of this script: instead of reading a
# precomputed, CA-specific ~90 m TIF (composite_YYYY_ca_100m.tif, built by
# 01_create_emapr_100m_tifs.R), this script reads the shared native ~30 m
# West-wide crop (00_crop_emapr_to_west.R, data/processed/emapr_biomass_west/)
# and crops + aggregates to ~90 m PER FIRE POLYGON on the fly. No new
# precompute script or extra whole-West ~90 m file is needed — this mirrors
# the crop-first philosophy already used for the forest masks and for 08's
# ctrees extraction (see those scripts' comments), and keeps this script
# working incrementally as 00_crop_emapr_to_west.R fills in more years.
#
# Forest mask: per-state 90 m 0/1 mask from 05_prepare_forest_masks_west.R
# (data/processed/forest_mask/nlcd2004_forestfrac_90m_<st>.tif) — must already
# exist for every state in STATES_TO_RUN.
#
# MTBS state attribution: uses the event_id prefix (authoritative) after a
# spatial join to the West states union, not the raw spatial-join match alone
# — same border-fire-safe dedup as 06_extract_pct_forest_within_fires.R. This
# should resolve fires along state borders being dropped by the old CA-polygon
# `st_filter()` approach (see the eMapR-vs-ctrees fire-count gap flagged in
# biomass_within_fires.qmd §7).
#
# Output is a single combined CSV with a STUSPS column, resumable per
# state x year — delete rows (or the whole file) to force re-extraction.
#
# Run ONCE from the project root before rendering analysis/biomass_within_fires.qmd:
#   Rscript scripts/r/07_extract_emapr_within_fires_new.R
#
# Output: data/processed/emapr_biomass_west/biomass_fire_polygons_emapr_west_<years>_100m_forested.csv
#
# NOTE: analysis/biomass_within_fires.qmd currently reads the retired CA-only
# output path (data/processed/emapr_biomass_ca/..._newpipeline.csv). Update
# the qmd's EMAPR_FIRE_CSV path (and add a STUSPS filter/facet) before relying
# on this script's new output there — not done as part of this rewrite.
#
# OUTLINE
# 1. Load packages and paths
# 2. Load MTBS study fires (event_id-prefix dedup, all requested states)
# 3. Identify available West-wide ~30 m TIFs
# 4. Load per-state forest masks; check cache — resume per state x year
# 5. Extract mean AGB per fire polygon per year: crop -> aggregate(~90m) -> mask
# 6. Write CSV cache
# =============================================================================

library(terra)
library(sf)
library(here)
library(dplyr)
library(readr)
library(glue)
library(lubridate)
library(tigris)
library(stringr)
library(tibble)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/07_extract_emapr_within_fires_new.R")

# ── 1. Paths and parameters ───────────────────────────────────────────────────
STUDY_YEARS    <- 2005L:2010L   # matches biomass_within_fires.qmd params (fire ignition year filter)
WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

# Override for a pilot run, e.g.: STATES_TO_RUN <- c("CA", "WY")
# Each state must already have a 90 m forest mask (05_prepare_forest_masks_west.R).
STATES_TO_RUN <- c("CA")

MTBS_PATH <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
EMAPR_DIR <- here("data", "processed", "emapr_biomass_west")
MASK_DIR  <- here("data", "processed", "forest_mask")

years_hash     <- paste(min(STUDY_YEARS), max(STUDY_YEARS), sep = "_")
EMAPR_FIRE_CSV <- here("data", "processed", "emapr_biomass_west",
                       glue("biomass_fire_polygons_emapr_west_{years_hash}_100m_forested.csv"))

cat("Study years:   ", paste(STUDY_YEARS, collapse = ", "), "\n")
cat("States to run: ", paste(STATES_TO_RUN, collapse = ", "), "\n")
cat("Output cache:  ", basename(EMAPR_FIRE_CSV), "\n\n")

# ── 2. Load MTBS study fires — event_id-prefix dedup (mirrors 06) ────────────
cat("Loading MTBS...\n")
mtbs_raw <- sf::st_read(MTBS_PATH, quiet = TRUE)

invalid <- sum(!sf::st_is_valid(mtbs_raw))
if (invalid > 0) mtbs_raw <- sf::st_make_valid(mtbs_raw)

mtbs_wf <- mtbs_raw |>
  dplyr::filter(incid_type == "Wildfire") |>
  dplyr::mutate(year = lubridate::year(ig_date)) |>
  dplyr::filter(year %in% STUDY_YEARS, burnbndac >= 1000)

cat("Loading West states boundary...\n")
western_states_sf <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS %in% WESTERN_STATES) |>
  sf::st_transform(sf::st_crs(mtbs_raw))

mtbs_joined <- sf::st_join(mtbs_wf, western_states_sf["STUSPS"], join = sf::st_intersects, left = FALSE)

# event_id's own state prefix is authoritative, not the (possibly duplicated)
# spatial-join match — same fix as 06_extract_pct_forest_within_fires.R.
mtbs_study <- mtbs_joined |>
  dplyr::mutate(STUSPS = substr(event_id, 1, 2)) |>
  dplyr::filter(STUSPS %in% STATES_TO_RUN) |>
  dplyr::distinct(event_id, .keep_all = TRUE) |>
  sf::st_transform(5070)   # matches eMapR West crop CRS and forest mask CRS

stopifnot("Duplicate event_ids after dedup" = !any(duplicated(mtbs_study$event_id)))
cat("Study fires:   ", nrow(mtbs_study), "across", length(unique(mtbs_study$STUSPS)), "state(s)\n\n")

# ── 3. Inventory available West-wide ~30 m TIFs ───────────────────────────────
tifs_west_ready <- sort(as.integer(stringr::str_extract(
  list.files(EMAPR_DIR, pattern = "^composite_\\d{4}_west\\.tif$"),
  "\\d{4}"
)))
cat("West-wide TIFs available:", paste(tifs_west_ready, collapse = ", "), "\n\n")

if (length(tifs_west_ready) == 0) {
  stop("No West-wide eMapR TIFs found. Run scripts/r/00_crop_emapr_to_west.R first.")
}

# ── 4. Forest masks and resume check ──────────────────────────────────────────
missing_masks <- STATES_TO_RUN[!file.exists(file.path(
  MASK_DIR, glue("nlcd2004_forestfrac_90m_{tolower(STATES_TO_RUN)}.tif")
))]
if (length(missing_masks) > 0) {
  stop("Missing 90 m forest mask for: ", paste(missing_masks, collapse = ", "),
       "\nRun scripts/r/05_prepare_forest_masks_west.R first ",
       "(STATES_TO_RUN must include these states).")
}

forest_masks <- stats::setNames(
  lapply(STATES_TO_RUN, function(st) {
    terra::rast(file.path(MASK_DIR, glue("nlcd2004_forestfrac_90m_{tolower(st)}.tif")))
  }),
  STATES_TO_RUN
)
cat("Forest masks loaded for:", paste(STATES_TO_RUN, collapse = ", "), "\n\n")

needed <- expand.grid(STUSPS = STATES_TO_RUN, year = tifs_west_ready,
                      stringsAsFactors = FALSE) |> tibble::as_tibble()

done <- if (file.exists(EMAPR_FIRE_CSV)) {
  readr::read_csv(EMAPR_FIRE_CSV, show_col_types = FALSE) |>
    dplyr::distinct(STUSPS, year)
} else {
  tibble::tibble(STUSPS = character(0), year = integer(0))
}

missing <- dplyr::anti_join(needed, done, by = c("STUSPS", "year"))

if (nrow(missing) == 0) {
  cat("Cache is up to date for all requested state x year combinations.\n")
  cat("Nothing to do. Delete", basename(EMAPR_FIRE_CSV), "to force re-extraction.\n")
  quit(save = "no")
}

years_to_do <- sort(unique(missing$year))
cat("State x year combinations remaining:", nrow(missing),
    "across", length(years_to_do), "year(s)\n\n")

# ── 5. Extract mean AGB per fire x year: crop -> aggregate(~90m) -> mask ─────
terra::terraOptions(threads = 1, progress = 0)

n_total <- length(years_to_do)
t0 <- proc.time()

for (i in seq_along(years_to_do)) {
  yr  <- years_to_do[[i]]
  tif <- file.path(EMAPR_DIR, glue("composite_{yr}_west.tif"))

  if (!file.exists(tif)) {
    cat(glue("  [{i}/{n_total}] {yr} — West TIF not found (crop not yet done for ",
             "this year), skipping\n"))
    next
  }

  states_needed <- sort(unique(missing$STUSPS[missing$year == yr]))
  fires_yr <- mtbs_study |> dplyr::filter(STUSPS %in% states_needed)
  if (nrow(fires_yr) == 0) next

  t1 <- proc.time()
  r       <- terra::rast(tif)
  vect_yr <- terra::vect(fires_yr)

  # Polygon-by-polygon: crop the native ~30 m raster to each fire's extent
  # first (cheap, small read), THEN aggregate that small crop to ~90 m —
  # equivalent to the old whole-West precomputed ~90 m TIF for a single
  # polygon's mean, without ever materializing a full-West ~90 m file.
  #
  # Padding fix (same root cause as the retired precomputed-TIF version):
  # deriving the mask crop from r_agg's OWN extent (not independently from
  # the polygon) guarantees the mask crop can never be smaller than the
  # aggregated raster, avoiding spurious NA at the crop boundary that
  # maskvalues=0 would otherwise silently read as non-forest.
  agb_vec <- numeric(nrow(vect_yr))
  for (j in seq_len(nrow(vect_yr))) {
    poly  <- vect_yr[j]
    st_j  <- vect_yr$STUSPS[j]

    # A single problematic polygon (transient I/O hiccup, GDAL block-read
    # failure, a genuinely degenerate MTBS geometry, etc.) skips to NA with a
    # logged event_id instead of halting a multi-hour, many-state run.
    agb_vec[j] <- tryCatch({
      r_c     <- terra::crop(r, poly, snap = "out")
      r_agg   <- terra::aggregate(r_c, fact = 3, fun = "mean", na.rm = TRUE)
      fm_c    <- terra::crop(forest_masks[[st_j]], terra::ext(r_agg) + 100, snap = "out")
      fm_snap <- terra::resample(fm_c, r_agg, method = "near")
      r_masked <- terra::mask(r_agg, fm_snap, maskvalues = 0)
      vals     <- terra::values(r_masked, na.rm = TRUE)
      if (length(vals) > 0L) mean(vals) else NA_real_
    }, error = function(e) {
      cat(glue("    [WARN] {vect_yr$event_id[j]} ({yr}) failed: {conditionMessage(e)} — set NA\n"))
      NA_real_
    })
  }

  elapsed <- round((proc.time() - t1)[["elapsed"]], 1)

  df_yr <- data.frame(
    event_id      = vect_yr$event_id,
    STUSPS        = vect_yr$STUSPS,
    year          = yr,
    agb_mean_mgha = agb_vec
  )

  readr::write_csv(df_yr, EMAPR_FIRE_CSV, append = file.exists(EMAPR_FIRE_CSV))

  cat(glue("  [{i}/{n_total}] {yr} — {nrow(df_yr)} fire(s) across ",
           "{length(states_needed)} state(s) — {elapsed}s (saved)\n"))

  rm(r, df_yr)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE, full = TRUE)
}

# ── 6. Report ─────────────────────────────────────────────────────────────────
final <- if (file.exists(EMAPR_FIRE_CSV)) readr::read_csv(EMAPR_FIRE_CSV, show_col_types = FALSE) else NULL
total_elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)

if (!is.null(final)) {
  cat(glue("\nDone. {nrow(final)} rows in {basename(EMAPR_FIRE_CSV)}",
           " ({length(unique(final$year))} year(s), {length(unique(final$STUSPS))} state(s))",
           " — {total_elapsed} min this run.\n"))
} else {
  cat("\nNo TIFs were available to process this run.\n")
}
cat("You can now render analysis/biomass_within_fires.qmd",
    "(after updating its EMAPR_FIRE_CSV path — see header note above).\n")
