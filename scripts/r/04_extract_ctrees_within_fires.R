# =============================================================================
# 04_extract_ctrees_within_fires.R
#
# Extract mean ctrees AGB within MTBS fire perimeters using local 100 m TIFs,
# applying the NLCD 2004 forest mask before extraction.
#
# Retired: this feeds analysis/biomass_within_fires_old.qmd only, which used the
# original CA-only forest-mask pipeline (03_prepare_forest_mask.R). The current
# pipeline is 08_extract_ctrees_within_fires_new.R + 05_prepare_forest_masks_west.R,
# feeding analysis/biomass_within_fires.qmd.
#
# Run ONCE after 03_prepare_forest_mask.R has produced the masks.
# Run from the project root before rendering biomass_within_fires_old.qmd.
#
# Why a separate script: same reason as 02_extract_emapr_within_fires.R —
# terra::extract() on Windows can appear frozen inside Quarto.
#
# Output:
#   data/processed/ctrees/biomass_fire_polygons_ctrees_forested.csv
#
# OUTLINE
# 1. Setup
# 2. Load CA boundary and MTBS study fires
# 3. Inventory available ctrees 100 m TIFs
# 4. Check cache — skip if up to date
# 5. Extract mean forest AGB per fire polygon per year
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

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/04_extract_ctrees_within_fires.R")

STUDY_YEARS <- 2005L:2010L   # must match QMD params (study_year_min / study_year_max)
STATE_FIPS  <- "CA"

MTBS_PATH    <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
CTREES_DIR   <- here("data", "processed", "ctrees")
FOREST_MASK  <- here("data", "processed", "forest_mask", "nlcd2004_forest_100m_ca.tif")
OUT_CSV      <- here("data", "processed", "ctrees", "biomass_fire_polygons_ctrees_forested.csv")

cat("Study years:  ", paste(STUDY_YEARS, collapse = ", "), "\n")
cat("Output cache: ", basename(OUT_CSV), "\n\n")

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

# Transform to EPSG:4326 to match ctrees TIF CRS
mtbs_4326 <- sf::st_transform(mtbs_study, 4326)
cat("Study fires:  ", nrow(mtbs_study), "\n\n")

# ── 3. Inventory available ctrees 100 m TIFs ──────────────────────────────────
tifs_available <- sort(as.integer(stringr::str_extract(
  list.files(CTREES_DIR, pattern = "^ctrees_\\d{4}_ca_100m\\.tif$"),
  "\\d{4}"
)))
cat("ctrees 100 m TIFs available:", paste(tifs_available, collapse = ", "), "\n\n")

if (length(tifs_available) == 0) {
  stop("No ctrees 100m TIFs found. Run scripts/python/03_download_ctrees_ca.py first.")
}

# ── 4. Forest mask and resume check ───────────────────────────────────────────
if (!file.exists(FOREST_MASK)) {
  stop("Forest mask not found: ", FOREST_MASK,
       "\nRun scripts/r/03_prepare_forest_mask.R first.")
}
forest_mask <- terra::rast(FOREST_MASK)
cat("Forest mask loaded:", basename(FOREST_MASK), "\n\n")

tifs_needed <- tifs_available  # extract all available years, not just fire cohort years

done_years <- integer(0)
if (file.exists(OUT_CSV)) {
  partial    <- readr::read_csv(OUT_CSV, show_col_types = FALSE)
  done_years <- sort(unique(partial$year))
  if (all(tifs_needed %in% done_years)) {
    cat("Cache is up to date — covers all requested years.\n")
    cat("Nothing to do. Delete", basename(OUT_CSV), "to force re-extraction.\n")
    quit(save = "no")
  }
  cat("Partial cache found —", length(done_years), "year(s) already done:",
      paste(done_years, collapse = ", "), "\nResuming...\n\n")
}

tifs_to_do <- setdiff(tifs_needed, done_years)

# ── 5. Extract mean forest AGB per fire × year ────────────────────────────────
terra::terraOptions(threads = 1, progress = 0)

mtbs_vect <- terra::vect(mtbs_4326)
n_total   <- length(tifs_available)
n_todo    <- length(tifs_to_do)

cat(glue("Extracting {n_todo} TIF(s) for {nrow(mtbs_study)} fire polygons (forest pixels only)...\n\n"))
t0 <- proc.time()

for (i in seq_along(tifs_to_do)) {
  yr  <- tifs_to_do[[i]]
  pos <- which(tifs_available == yr)
  tif <- file.path(CTREES_DIR, glue("ctrees_{yr}_ca_100m.tif"))

  if (!file.exists(tif)) {
    cat(glue("  [{pos}/{n_total}] {yr} — TIF not found, skipping\n"))
    next
  }

  t1 <- proc.time()
  r  <- terra::rast(tif)

  # Polygon-by-polygon: crop both rasters to each polygon's extent, mask, mean.
  # Each crop reads only a few disk blocks (~KB RAM peak vs ~2 GB for full-raster ops).
  agb_vec <- numeric(nrow(mtbs_4326))
  for (j in seq_len(nrow(mtbs_4326))) {
    poly      <- mtbs_vect[j]
    r_c       <- terra::crop(r,           poly, snap = "out")
    fm_c      <- terra::crop(forest_mask, poly, snap = "out")
    r_masked  <- terra::mask(r_c, fm_c)
    vals      <- terra::values(r_masked, na.rm = TRUE)
    agb_vec[j] <- if (length(vals) > 0L) mean(vals) else NA_real_
  }

  elapsed <- round((proc.time() - t1)[["elapsed"]], 1)

  df_yr <- data.frame(
    event_id      = mtbs_vect$event_id,
    fire_year     = as.integer(mtbs_vect$year),
    year          = yr,
    agb_mean_mgha = agb_vec
  )

  readr::write_csv(df_yr, OUT_CSV, append = file.exists(OUT_CSV))

  cat(glue("  [{pos}/{n_total}] {yr} — {elapsed}s  (saved)\n"))

  rm(r, ex, df_yr)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE, full = TRUE)
}

# ── 6. Report ─────────────────────────────────────────────────────────────────
final         <- readr::read_csv(OUT_CSV, show_col_types = FALSE)
total_elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
cat(glue("\nDone. {nrow(final)} rows in {basename(OUT_CSV)}",
         " ({length(unique(final$year))} years) — {total_elapsed} min this run.\n"))
cat("You can now render analysis/biomass_within_fires_old.qmd.\n")
