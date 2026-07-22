# =============================================================================
# 02_extract_emapr_within_fires.R
#
# One-time preprocessing: extract mean eMapR AGB within MTBS fire perimeters
# for all available ~100 m CA TIFs and cache results to CSV.
#
# Run ONCE from the project root before rendering 05_biomass_within_fires.qmd:
#   Rscript scripts/r/02_extract_emapr_within_fires.R
#
# Why a separate script: terra::extract() on Windows can appear frozen inside
# Quarto (output is buffered until the chunk finishes). Running here gives
# real-time progress and avoids Quarto timeouts.
#
# Output: data/processed/emapr_biomass_ca/biomass_fire_polygons_emapr_<years>_100m.csv
#
# OUTLINE
# 1. Load packages and paths
# 2. Load MTBS study fires (same filters as QMD)
# 3. Identify available ~100 m TIFs
# 4. Check cache — skip if up to date
# 5. Extract mean AGB per fire polygon per year
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

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/02_extract_emapr_within_fires.R")

# ── 1. Paths and parameters ───────────────────────────────────────────────────
STUDY_YEARS <- 2005L:2010L   # must match QMD params (study_year_min / study_year_max)
STATE_FIPS  <- "CA"

MTBS_PATH   <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
EMAPR_DIR   <- here("data", "processed", "emapr_biomass_ca")
FOREST_MASK <- here("data", "processed", "forest_mask", "nlcd2004_forest_90m_ca.tif")

years_hash     <- paste(min(STUDY_YEARS), max(STUDY_YEARS), sep = "_")
EMAPR_FIRE_CSV <- here("data", "processed", "emapr_biomass_ca",
                       glue("biomass_fire_polygons_emapr_{years_hash}_100m_forested.csv"))

cat("Study years:   ", paste(STUDY_YEARS, collapse = ", "), "\n")
cat("Output cache:  ", basename(EMAPR_FIRE_CSV), "\n\n")

# ── 2. Load MTBS study fires ──────────────────────────────────────────────────
cat("Loading CA boundary...\n")
ca_sf   <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS == STATE_FIPS)
ca_5070 <- sf::st_transform(ca_sf, 5070)

cat("Loading MTBS...\n")
mtbs_raw <- sf::st_read(MTBS_PATH, quiet = TRUE)

mtbs_study <- mtbs_raw |>
  dplyr::filter(incid_type == "Wildfire") |>
  dplyr::mutate(year = lubridate::year(ig_date)) |>
  dplyr::filter(year %in% STUDY_YEARS) |>
  sf::st_transform(5070) |>
  sf::st_filter(ca_5070) |>
  dplyr::filter(startsWith(event_id, STATE_FIPS)) |>
  dplyr::filter(burnbndac >= 1000)

invalid <- sum(!sf::st_is_valid(mtbs_study))
if (invalid > 0) mtbs_study <- sf::st_make_valid(mtbs_study)

cat("Study fires:   ", nrow(mtbs_study), "\n\n")

# ── 3. Inventory available ~100 m TIFs ───────────────────────────────────────
tifs_100m_ready <- sort(as.integer(stringr::str_extract(
  list.files(EMAPR_DIR, pattern = "^composite_\\d{4}_ca_100m\\.tif$"),
  "\\d{4}"
)))
cat("100 m TIFs available:", paste(tifs_100m_ready, collapse = ", "), "\n\n")

if (length(tifs_100m_ready) == 0) {
  stop("No ~100 m TIFs found. Run scripts/r/01_create_emapr_100m_tifs.R first.")
}

# ── 4. Resume check — skip years already in a partial cache ──────────────────
# Writes one year at a time so a crash mid-run loses nothing already extracted.
tifs_needed <- tifs_100m_ready  # extract all available years, not just fire cohort years

done_years <- integer(0)
if (file.exists(EMAPR_FIRE_CSV)) {
  partial    <- readr::read_csv(EMAPR_FIRE_CSV, show_col_types = FALSE)
  done_years <- sort(unique(partial$year))
  if (all(tifs_needed %in% done_years)) {
    cat("Cache is up to date — covers all requested years.\n")
    cat("Nothing to do. Delete", basename(EMAPR_FIRE_CSV), "to force re-extraction.\n")
    quit(save = "no")
  }
  cat("Partial cache found —", length(done_years), "year(s) already done:",
      paste(done_years, collapse = ", "), "\nResuming...\n\n")
}

tifs_to_do <- setdiff(tifs_needed, done_years)

# ── 5. Extract mean AGB per fire × year ──────────────────────────────────────
# terra single-threaded: avoids Windows threading deadlocks.
# Results written after each year so a std::bad_alloc crash is resumable.
terra::terraOptions(threads = 1, progress = 0)

if (!file.exists(FOREST_MASK)) {
  stop("Forest mask not found: ", FOREST_MASK,
       "\nRun scripts/r/03_prepare_forest_mask.R first.")
}
forest_mask <- terra::rast(FOREST_MASK)
cat("Forest mask loaded:", basename(FOREST_MASK), "\n\n")

mtbs_vect <- terra::vect(mtbs_study)   # EPSG:5070, matches eMapR TIF CRS
n_total   <- length(tifs_100m_ready)
n_todo    <- length(tifs_to_do)

cat(glue("Extracting {n_todo} TIF(s) for {nrow(mtbs_study)} fire polygons (polygon-by-polygon)...\n\n"))
t0 <- proc.time()

for (i in seq_along(tifs_to_do)) {
  yr  <- tifs_to_do[[i]]
  pos <- which(tifs_100m_ready == yr)
  tif <- file.path(EMAPR_DIR, glue("composite_{yr}_ca_100m.tif"))

  if (!file.exists(tif)) {
    cat(glue("  [{pos}/{n_total}] {yr} — TIF not found, skipping\n"))
    next
  }

  t1 <- proc.time()
  r  <- terra::rast(tif)

  # Polygon-by-polygon: crop raster + mask to each fire polygon's extent, then mean.
  # Full-raster resample/mask/extract segfaults on Windows with terra regardless
  # of batching — the only stable path is the same approach used for ctrees.
  # Resamples the 90m forest mask to the eMapR 100m grid at polygon scale (fast).
  agb_vec <- numeric(nrow(mtbs_vect))
  for (j in seq_len(nrow(mtbs_vect))) {
    poly     <- mtbs_vect[j]
    r_c      <- terra::crop(r,           poly, snap = "out")
    fm_c     <- terra::crop(forest_mask, poly, snap = "out")
    fm_snap  <- terra::resample(fm_c, r_c, method = "near")
    r_masked <- terra::mask(r_c, fm_snap)
    vals     <- terra::values(r_masked, na.rm = TRUE)
    agb_vec[j] <- if (length(vals) > 0L) mean(vals) else NA_real_
  }

  elapsed <- round((proc.time() - t1)[["elapsed"]], 1)

  df_yr <- data.frame(
    event_id      = mtbs_vect$event_id,
    year          = yr,
    agb_mean_mgha = agb_vec
  )

  readr::write_csv(df_yr, EMAPR_FIRE_CSV,
                   append = file.exists(EMAPR_FIRE_CSV))

  cat(glue("  [{pos}/{n_total}] {yr} — {elapsed}s  (saved)\n"))

  rm(r, df_yr)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE, full = TRUE)
}

# ── 6. Report ─────────────────────────────────────────────────────────────────
final      <- readr::read_csv(EMAPR_FIRE_CSV, show_col_types = FALSE)
total_elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
cat(glue("\nDone. {nrow(final)} rows in {basename(EMAPR_FIRE_CSV)}",
         " ({length(unique(final$year))} years) — {total_elapsed} min this run.\n"))
cat("You can now render 05_biomass_within_fires.qmd.\n")
