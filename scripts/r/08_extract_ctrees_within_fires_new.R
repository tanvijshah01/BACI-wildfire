# =============================================================================
# 08_extract_ctrees_within_fires_new.R
#
# Extracts mean ctrees AGB within MTBS fire perimeters, forested pixels only,
# for a user-specified set of Western states (STATES_TO_RUN below). Feeds
# analysis/biomass_within_fires.qmd.
#
# Generalizes the CA-only version of this script: instead of reading a
# per-state ctrees TIF (ctrees_YYYY_ca_100m.tif, from
# scripts/python/03_download_ctrees_ca.py), this script reads the single
# shared West-wide ctrees TIF (ctrees_YYYY_west_100m.tif, from
# scripts/python/04_download_ctrees_west.py) and crops it per state per fire
# polygon — the per-polygon crop-then-mask pattern here was already
# state-agnostic, so the only change is the source TIF and looping over
# multiple states' fires/masks instead of one hardcoded state.
#
# CRS note: the mask is EPSG:5070 (30 m native NLCD), while ctrees TIFs are
# EPSG:4326 (~100 m). This script projects the mask onto each cropped ctrees
# fragment PER POLYGON rather than precomputing a second whole-region aligned
# mask file — consistent with the crop-first philosophy used everywhere else
# in this pipeline for Windows memory safety (see 07's header for the same
# approach applied to eMapR).
#
# Forest mask: per-state 30 m 0/1 mask from 05_prepare_forest_masks_west.R
# (data/processed/forest_mask/nlcd2004_forestfrac_30m_<st>.tif) — must already
# exist for every state in STATES_TO_RUN.
#
# MTBS state attribution: uses the event_id prefix (authoritative) after a
# spatial join to the West states union, not the raw spatial-join match alone
# — same border-fire-safe dedup as 06_extract_pct_forest_within_fires.R and
# 07_extract_emapr_within_fires_new.R.
#
# Output is a single combined CSV with a STUSPS column, resumable per
# state x year — delete rows (or the whole file) to force re-extraction.
#
# Run ONCE from the project root before rendering analysis/biomass_within_fires.qmd:
#   Rscript scripts/r/08_extract_ctrees_within_fires_new.R
#
# Output: data/processed/ctrees/biomass_fire_polygons_ctrees_west_forested.csv
#
# NOTE: requires ctrees_YYYY_west_100m.tif files from
# scripts/python/04_download_ctrees_west.py (multi-hour run, West-wide) —
# these did not exist yet as of this rewrite. This script processes whatever
# years are present and resumes automatically as more land.
#
# NOTE: analysis/biomass_within_fires.qmd currently reads the retired CA-only
# output path (data/processed/ctrees/..._newpipeline.csv). Update the qmd's
# ctrees CSV path (and add a STUSPS filter/facet) before relying on this
# script's new output there — not done as part of this rewrite.
#
# OUTLINE
# 1. Setup
# 2. Load MTBS study fires (event_id-prefix dedup, all requested states)
# 3. Inventory available West-wide ctrees 100 m TIFs
# 4. Load per-state forest masks; check cache — resume per state x year
# 5. Extract mean forest AGB per fire polygon per year, using the 0/1 mask
# 6. Write CSV cache and report
# =============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────
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
here::i_am("scripts/r/08_extract_ctrees_within_fires_new.R")

STUDY_YEARS    <- 2005L:2010L   # matches biomass_within_fires.qmd params (fire ignition year filter)
WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

# Override for a pilot run, e.g.: STATES_TO_RUN <- c("CA", "WY")
# Each state must already have a 30 m forest mask (05_prepare_forest_masks_west.R).
STATES_TO_RUN <- c("CA")

MTBS_PATH  <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
CTREES_DIR <- here("data", "processed", "ctrees")
MASK_DIR   <- here("data", "processed", "forest_mask")
OUT_CSV    <- here("data", "processed", "ctrees", "biomass_fire_polygons_ctrees_west_forested.csv")

cat("Study years:   ", paste(STUDY_YEARS, collapse = ", "), "\n")
cat("States to run: ", paste(STATES_TO_RUN, collapse = ", "), "\n")
cat("Output cache:  ", basename(OUT_CSV), "\n\n")

# ── 2. Load MTBS study fires — event_id-prefix dedup (mirrors 06/07) ─────────
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
# spatial-join match — same fix as 06/07.
mtbs_study <- mtbs_joined |>
  dplyr::mutate(STUSPS = substr(event_id, 1, 2)) |>
  dplyr::filter(STUSPS %in% STATES_TO_RUN) |>
  dplyr::distinct(event_id, .keep_all = TRUE) |>
  sf::st_transform(5070)   # matches forest mask CRS; reprojected to 4326 below for ctrees

stopifnot("Duplicate event_ids after dedup" = !any(duplicated(mtbs_study$event_id)))

# 4326 to crop the ctrees raster (its native CRS). The mask crop is derived
# per-polygon from r_c's own (reprojected) extent — see §5 below.
mtbs_study_4326 <- sf::st_transform(mtbs_study, 4326)
cat("Study fires:   ", nrow(mtbs_study), "across", length(unique(mtbs_study$STUSPS)), "state(s)\n\n")

# ── 3. Inventory available West-wide ctrees 100 m TIFs ────────────────────────
tifs_available <- sort(as.integer(stringr::str_extract(
  list.files(CTREES_DIR, pattern = "^ctrees_\\d{4}_west_100m\\.tif$"),
  "\\d{4}"
)))
cat("West-wide ctrees 100 m TIFs available:", paste(tifs_available, collapse = ", "), "\n\n")

if (length(tifs_available) == 0) {
  stop("No West-wide ctrees 100m TIFs found. Run ",
       "scripts/python/04_download_ctrees_west.py first (Part C).")
}

# ── 4. Forest masks and resume check ──────────────────────────────────────────
missing_masks <- STATES_TO_RUN[!file.exists(file.path(
  MASK_DIR, glue("nlcd2004_forestfrac_30m_{tolower(STATES_TO_RUN)}.tif")
))]
if (length(missing_masks) > 0) {
  stop("Missing 30 m forest mask for: ", paste(missing_masks, collapse = ", "),
       "\nRun scripts/r/05_prepare_forest_masks_west.R first ",
       "(STATES_TO_RUN must include these states).")
}

forest_masks <- stats::setNames(
  lapply(STATES_TO_RUN, function(st) {
    terra::rast(file.path(MASK_DIR, glue("nlcd2004_forestfrac_30m_{tolower(st)}.tif")))
  }),
  STATES_TO_RUN
)
cat("Forest masks loaded for:", paste(STATES_TO_RUN, collapse = ", "), "\n\n")

needed <- expand.grid(STUSPS = STATES_TO_RUN, year = tifs_available,
                      stringsAsFactors = FALSE) |> tibble::as_tibble()

done <- if (file.exists(OUT_CSV)) {
  readr::read_csv(OUT_CSV, show_col_types = FALSE) |>
    dplyr::distinct(STUSPS, year)
} else {
  tibble::tibble(STUSPS = character(0), year = integer(0))
}

missing <- dplyr::anti_join(needed, done, by = c("STUSPS", "year"))

if (nrow(missing) == 0) {
  cat("Cache is up to date for all requested state x year combinations.\n")
  cat("Nothing to do. Delete", basename(OUT_CSV), "to force re-extraction.\n")
  quit(save = "no")
}

years_to_do <- sort(unique(missing$year))
cat("State x year combinations remaining:", nrow(missing),
    "across", length(years_to_do), "year(s)\n\n")

# ── 5. Extract mean forest AGB per fire x year, using the new 0/1 mask ───────
terra::terraOptions(threads = 1, progress = 0)

n_total <- length(years_to_do)
t0 <- proc.time()

for (i in seq_along(years_to_do)) {
  yr  <- years_to_do[[i]]
  tif <- file.path(CTREES_DIR, glue("ctrees_{yr}_west_100m.tif"))

  if (!file.exists(tif)) {
    cat(glue("  [{i}/{n_total}] {yr} — West TIF not found, skipping\n"))
    next
  }

  states_needed <- sort(unique(missing$STUSPS[missing$year == yr]))
  fires_yr_4326 <- mtbs_study_4326 |> dplyr::filter(STUSPS %in% states_needed)
  if (nrow(fires_yr_4326) == 0) next

  t1 <- proc.time()
  r       <- terra::rast(tif)
  vect_yr <- terra::vect(fires_yr_4326)

  # Polygon-by-polygon: crop the ctrees raster (4326) to each polygon's extent
  # first (cheap, small read), then derive the mask crop FROM r_c's own extent
  # (reprojected into the mask's CRS) rather than independently from the
  # polygon, before projecting the mask fragment onto r_c's exact grid.
  #
  # Padding fix: an axis-aligned rectangle in 4326 (r_c's extent) does not map
  # to an axis-aligned rectangle in 5070 — Albers Conic rotates/skews relative
  # to lat/lon depending on longitude. Cropping the mask independently from
  # the polygon's own (axis-aligned-in-5070) bounding box can therefore miss
  # the rotated corners of r_c's true footprint, leaving project() with no
  # source pixel there — it returns NA, and maskvalues=0 silently treats that
  # as non-forest (see biomass_within_fires.qmd §7). Fix: project r_c's own
  # extent into 5070 first, so the mask crop is guaranteed to cover r_c's true
  # rotated footprint, then pad with a small safety margin.
  agb_vec <- numeric(nrow(vect_yr))
  for (j in seq_len(nrow(vect_yr))) {
    poly_4326 <- vect_yr[j]
    st_j      <- vect_yr$STUSPS[j]

    r_c        <- terra::crop(r, poly_4326, snap = "out")
    r_ext_5070 <- terra::project(terra::as.polygons(r_c, extent = TRUE), "EPSG:5070")
    fm_c_5070  <- terra::crop(forest_masks[[st_j]], terra::ext(r_ext_5070) + 200, snap = "out")
    fm_c_4326  <- terra::project(fm_c_5070, r_c, method = "near")
    r_masked   <- terra::mask(r_c, fm_c_4326, maskvalues = 0)
    vals       <- terra::values(r_masked, na.rm = TRUE)
    agb_vec[j] <- if (length(vals) > 0L) mean(vals) else NA_real_
  }

  elapsed <- round((proc.time() - t1)[["elapsed"]], 1)

  df_yr <- data.frame(
    event_id      = vect_yr$event_id,
    STUSPS        = vect_yr$STUSPS,
    fire_year     = as.integer(vect_yr$year),
    year          = yr,
    agb_mean_mgha = agb_vec
  )

  readr::write_csv(df_yr, OUT_CSV, append = file.exists(OUT_CSV))

  cat(glue("  [{i}/{n_total}] {yr} — {nrow(df_yr)} fire(s) across ",
           "{length(states_needed)} state(s) — {elapsed}s (saved)\n"))

  rm(r, df_yr)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE, full = TRUE)
}

# ── 6. Report ─────────────────────────────────────────────────────────────────
final <- if (file.exists(OUT_CSV)) readr::read_csv(OUT_CSV, show_col_types = FALSE) else NULL
total_elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)

if (!is.null(final)) {
  cat(glue("\nDone. {nrow(final)} rows in {basename(OUT_CSV)}",
           " ({length(unique(final$year))} year(s), {length(unique(final$STUSPS))} state(s))",
           " — {total_elapsed} min this run.\n"))
} else {
  cat("\nNo TIFs were available to process this run.\n")
}
cat("You can now render analysis/biomass_within_fires.qmd",
    "(after updating its ctrees CSV path — see header note above).\n")
