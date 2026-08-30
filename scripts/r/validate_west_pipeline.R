# =============================================================================
# validate_west_pipeline.R
#
# Regression/validation harness for the West-wide, multi-state pipeline
# (05_prepare_forest_masks_west.R, 06_extract_pct_forest_within_fires.R,
# 07_extract_emapr_within_fires_new.R, 08_extract_ctrees_within_fires_new.R)
# rewritten 2026-08-12, run BEFORE spending more time cropping/downloading new
# years or scaling to more states. Not part of the numbered production
# pipeline — this is a one-off (re-run as needed) check, not a pipeline step.
#
# Answers four questions with real data, not code review:
#   1. Does the new event_id-prefix MTBS dedup (used by 06/07/08) select the
#      same CA fires as the old single-state st_filter()+startsWith() dedup
#      it replaced? (Differences are expected at state borders — that's the
#      bug the new dedup fixes — but should be explainable, not arbitrary.)
#   2. Does the new dedup produce a sane fire set for a SECOND state (WY),
#      now that 05 has real WY forest masks?
#   3. Are all forest masks actually readable (not silently corrupt, per the
#      WY 30m mask found dead on 2026-08-12)?
#   4. Does today's 07/08 rewrite reproduce the validated old-CA-pipeline's
#      per-fire AGB values? eMapR's extraction METHOD genuinely changed
#      (precomputed whole-CA ~90m TIF -> per-polygon native-30m crop +
#      aggregate), so some small difference is expected; ctrees' extraction
#      math is UNCHANGED (only the source file path changed), so it should
#      match almost exactly — any real drift there is a refactor bug, not a
#      methodology difference, and needs to be caught NOW before the West
#      ctrees download lands and this pipeline is trusted for real output.
#
# Uses ONLY data already on disk today (CA's existing raw 30m + precomputed
# ~100m eMapR/ctrees TIFs, and the old-pipeline CSVs they fed) — does not
# require the West eMapR crop or West ctrees download to be complete, so it
# is not blocked by either of those in-progress jobs.
#
# Run from the project root (NOT inside Quarto — see CLAUDE.md's terra/Windows
# note):
#   Rscript scripts/r/validate_west_pipeline.R
#
# Outputs (data/processed/validation/), consumed by
# analysis/west_pipeline_sanity_check.qmd:
#   mtbs_dedup_comparison_ca.csv   — per-fire old-vs-new inclusion
#   wy_mtbs_smoke_test.csv         — new-dedup fire list for WY
#   forest_mask_validity_report.csv
#   emapr_method_comparison_ca.csv — per fire x year, old vs new AGB
#   ctrees_method_comparison_ca.csv
#
# OUTLINE
# 1. Setup
# 2. MTBS dedup comparison — CA, old filter vs new filter
# 3. MTBS dedup smoke test — WY, new filter only
# 4. Forest mask validity report — CA + WY, all tiers
# 5. eMapR method comparison — CA, STUDY_YEARS, old cache vs new-method re-extraction
# 6. ctrees method comparison — CA, STUDY_YEARS, old cache vs new-method re-extraction
# 7. Console summary
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
library(tibble)

sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)
here::i_am("scripts/r/validate_west_pipeline.R")

STUDY_YEARS    <- 2005L:2010L
WESTERN_STATES <- c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

MTBS_PATH   <- here("data", "raw", "mtbs", "mtbs_perimeter_data", "mtbs_perims_DD.shp")
MASK_DIR    <- here("data", "processed", "forest_mask")
EMAPR_CA_DIR <- here("data", "processed", "emapr_biomass_ca")
CTREES_DIR   <- here("data", "processed", "ctrees")
OUT_DIR      <- here("data", "processed", "validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OLD_EMAPR_CSV  <- here("data", "processed", "emapr_biomass_ca",
                       "biomass_fire_polygons_emapr_2005_2010_100m_forested_newpipeline.csv")
OLD_CTREES_CSV <- here("data", "processed", "ctrees",
                       "biomass_fire_polygons_ctrees_forested_newpipeline.csv")

raster_is_valid <- function(path) {
  if (!file.exists(path)) return(NA)
  terra::global(terra::rast(path), "notNA")[[1]] > 0
}

cat("=== validate_west_pipeline.R ===\n")
cat("Study years:", paste(STUDY_YEARS, collapse = ", "), "\n\n")

mtbs_raw <- sf::st_read(MTBS_PATH, quiet = TRUE)
invalid <- sum(!sf::st_is_valid(mtbs_raw))
if (invalid > 0) mtbs_raw <- sf::st_make_valid(mtbs_raw)

# ── 2. MTBS dedup comparison — CA, old filter vs new filter ──────────────────
cat("── 2. MTBS dedup: old (st_filter+startsWith) vs new (event_id-prefix) — CA ──\n")

ca_sf   <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS == "CA")
ca_5070 <- sf::st_transform(ca_sf, 5070)

mtbs_wf_common <- mtbs_raw |>
  dplyr::filter(incid_type == "Wildfire") |>
  dplyr::mutate(year = lubridate::year(ig_date)) |>
  dplyr::filter(year %in% STUDY_YEARS, burnbndac >= 1000)

# OLD: matches the pre-2026-08-12 07/08 (single-state st_filter + startsWith)
old_ca <- mtbs_wf_common |>
  sf::st_transform(5070) |>
  sf::st_filter(ca_5070) |>
  dplyr::filter(startsWith(event_id, "CA"))

# NEW: matches 06/07/08 as rewritten 2026-08-12 (West-states spatial join,
# then event_id prefix as the authoritative state)
western_states_sf <- tigris::states(cb = TRUE, year = 2022, resolution = "5m") |>
  dplyr::filter(STUSPS %in% WESTERN_STATES) |>
  sf::st_transform(sf::st_crs(mtbs_raw))
mtbs_joined <- sf::st_join(mtbs_wf_common, western_states_sf["STUSPS"],
                           join = sf::st_intersects, left = FALSE)
new_ca <- mtbs_joined |>
  dplyr::mutate(STUSPS = substr(event_id, 1, 2)) |>
  dplyr::filter(STUSPS == "CA") |>
  dplyr::distinct(event_id, .keep_all = TRUE) |>
  sf::st_transform(5070)   # matches eMapR/forest-mask CRS — used downstream as mtbs_common_5070

old_ids <- sort(unique(old_ca$event_id))
new_ids <- sort(unique(new_ca$event_id))

dedup_compare <- dplyr::full_join(
  tibble::tibble(event_id = old_ids, in_old = TRUE),
  tibble::tibble(event_id = new_ids, in_new = TRUE),
  by = "event_id"
) |>
  dplyr::mutate(in_old = !is.na(in_old), in_new = !is.na(in_new)) |>
  dplyr::left_join(
    sf::st_drop_geometry(mtbs_raw) |> dplyr::select(event_id, burnbndac, ig_date),
    by = "event_id"
  ) |>
  dplyr::arrange(event_id)

readr::write_csv(dedup_compare, file.path(OUT_DIR, "mtbs_dedup_comparison_ca.csv"))

n_only_old <- sum(dedup_compare$in_old & !dedup_compare$in_new)
n_only_new <- sum(!dedup_compare$in_old & dedup_compare$in_new)
n_both     <- sum(dedup_compare$in_old & dedup_compare$in_new)
cat(glue("  Old: {length(old_ids)} fires | New: {length(new_ids)} fires | ",
         "Both: {n_both} | Old-only: {n_only_old} | New-only: {n_only_new}\n\n"))

# ── 3. MTBS dedup smoke test — WY, new filter only ────────────────────────────
cat("── 3. MTBS dedup smoke test — WY (new filter only) ──\n")

new_wy <- mtbs_joined |>
  dplyr::mutate(STUSPS = substr(event_id, 1, 2)) |>
  dplyr::filter(STUSPS == "WY") |>
  dplyr::distinct(event_id, .keep_all = TRUE)

wy_summary <- sf::st_drop_geometry(new_wy) |>
  dplyr::select(event_id, year, burnbndac)
readr::write_csv(wy_summary, file.path(OUT_DIR, "wy_mtbs_smoke_test.csv"))

bad_prefix <- sum(!startsWith(new_wy$event_id, "WY"))
dup_ids    <- sum(duplicated(new_wy$event_id))
cat(glue("  WY fires: {nrow(new_wy)} | wrong event_id prefix: {bad_prefix} | ",
         "duplicate event_ids: {dup_ids}\n"))
cat(glue("  {if (nrow(new_wy) > 0 && bad_prefix == 0 && dup_ids == 0) 'PASS' else 'FAIL'} — ",
         "WY fire set is well-formed\n\n"))

# ── 4. Forest mask validity report — CA + WY, all tiers ──────────────────────
cat("── 4. Forest mask validity — CA + WY ──\n")

mask_report <- tibble::tibble()
for (st in c("ca", "wy")) {
  for (tier in c("30m", "90m", "100m")) {
    path <- file.path(MASK_DIR, glue("nlcd2004_forestfrac_{tier}_{st}.tif"))
    exists <- file.exists(path)
    valid  <- if (exists) raster_is_valid(path) else NA
    pct_forest <- NA_real_
    if (isTRUE(valid)) {
      r <- terra::rast(path)
      s <- terra::global(r, c("sum", "notNA"), na.rm = TRUE)
      pct_forest <- round(100 * s$sum / s$notNA, 1)
    }
    mask_report <- dplyr::bind_rows(mask_report, tibble::tibble(
      state = toupper(st), tier = tier, exists = exists, valid = valid,
      pct_forest = pct_forest
    ))
    status <- if (!exists) "MISSING" else if (!isTRUE(valid)) "*** CORRUPT ***" else glue("OK ({pct_forest}% forest)")
    cat(glue("  {toupper(st)} {tier}: {status}\n"))
  }
}
readr::write_csv(mask_report, file.path(OUT_DIR, "forest_mask_validity_report.csv"))
cat("\n")

# ── 5. eMapR method comparison — CA, old cache vs new-method re-extraction ───
cat("── 5. eMapR: old (precomputed ~90m TIF) vs new (per-polygon crop+aggregate) ──\n")

# Compare only fires present in BOTH dedup passes (§2) — isolates the
# extraction-method question from the dedup-set question (§2 already covers
# that separately).
common_fires <- dplyr::inner_join(
  sf::st_drop_geometry(old_ca) |> dplyr::select(event_id),
  sf::st_drop_geometry(new_ca) |> dplyr::select(event_id),
  by = "event_id"
)
mtbs_common_5070 <- new_ca |> dplyr::filter(event_id %in% common_fires$event_id)
cat(glue("  Comparing {nrow(mtbs_common_5070)} fires common to both dedup passes\n"))

if (!file.exists(OLD_EMAPR_CSV)) {
  cat("  OLD eMapR cache not found — skipping. Run 07 (pre-rewrite output) first.\n\n")
} else {
  old_emapr <- readr::read_csv(OLD_EMAPR_CSV, show_col_types = FALSE) |>
    dplyr::filter(event_id %in% mtbs_common_5070$event_id)

  forest_mask_90m <- terra::rast(file.path(MASK_DIR, "nlcd2004_forestfrac_90m_ca.tif"))
  terra::terraOptions(threads = 1, progress = 0)

  EMAPR_COMPARE_CSV <- file.path(OUT_DIR, "emapr_method_comparison_ca.csv")
  done_years <- integer(0)
  if (file.exists(EMAPR_COMPARE_CSV)) {
    done_years <- sort(unique(readr::read_csv(EMAPR_COMPARE_CSV, show_col_types = FALSE)$year))
  }
  years_to_do <- setdiff(STUDY_YEARS, done_years)

  vect_common <- terra::vect(mtbs_common_5070)

  for (yr in years_to_do) {
    tif <- file.path(EMAPR_CA_DIR, glue("composite_{yr}_ca.tif"))   # native 30 m
    if (!file.exists(tif)) {
      cat(glue("  {yr} — native 30m TIF not found, skipping\n"))
      next
    }
    t0 <- proc.time()
    r  <- terra::rast(tif)

    agb_vec <- numeric(nrow(vect_common))
    for (j in seq_len(nrow(vect_common))) {
      poly <- vect_common[j]
      # A single problematic polygon (transient I/O hiccup, GDAL block-read
      # failure, etc. — none reproduced on retest for this dataset, but a
      # multi-hour 11-state run shouldn't die over one fire) skips to NA with
      # a logged event_id instead of halting the whole run.
      agb_vec[j] <- tryCatch({
        r_c     <- terra::crop(r, poly, snap = "out")
        r_agg   <- terra::aggregate(r_c, fact = 3, fun = "mean", na.rm = TRUE)
        fm_c    <- terra::crop(forest_mask_90m, terra::ext(r_agg) + 100, snap = "out")
        fm_snap <- terra::resample(fm_c, r_agg, method = "near")
        r_masked <- terra::mask(r_agg, fm_snap, maskvalues = 0)
        vals     <- terra::values(r_masked, na.rm = TRUE)
        if (length(vals) > 0L) mean(vals) else NA_real_
      }, error = function(e) {
        cat(glue("    [WARN] {vect_common$event_id[j]} ({yr}) failed: {conditionMessage(e)} — set NA\n"))
        NA_real_
      })
    }

    df_yr <- data.frame(event_id = vect_common$event_id, year = yr, agb_new = agb_vec)
    readr::write_csv(df_yr, EMAPR_COMPARE_CSV, append = file.exists(EMAPR_COMPARE_CSV))

    elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
    cat(glue("  {yr} — {nrow(df_yr)} fires — {elapsed}s (saved)\n"))

    rm(r, df_yr)
    terra::tmpFiles(remove = TRUE)
    gc(verbose = FALSE, full = TRUE)
  }

  # Join old vs new and summarize
  new_emapr <- readr::read_csv(EMAPR_COMPARE_CSV, show_col_types = FALSE)
  emapr_diff <- dplyr::inner_join(old_emapr, new_emapr, by = c("event_id", "year")) |>
    dplyr::mutate(
      diff     = agb_new - agb_mean_mgha,
      pct_diff = 100 * diff / agb_mean_mgha
    )
  readr::write_csv(emapr_diff, file.path(OUT_DIR, "emapr_method_comparison_ca_diff.csv"))

  cat(glue("\n  eMapR old-vs-new: n={nrow(emapr_diff)} pairs | ",
           "median % diff = {round(median(emapr_diff$pct_diff, na.rm=TRUE), 2)}% | ",
           "mean abs % diff = {round(mean(abs(emapr_diff$pct_diff), na.rm=TRUE), 2)}% | ",
           "max abs % diff = {round(max(abs(emapr_diff$pct_diff), na.rm=TRUE), 2)}%\n\n"))
}

# ── 6. ctrees method comparison — CA, old cache vs new-method re-extraction ──
cat("── 6. ctrees: old (per-state TIF) vs new (shared-file code path) — math unchanged ──\n")

if (!file.exists(OLD_CTREES_CSV)) {
  cat("  OLD ctrees cache not found — skipping. Run 08 (pre-rewrite output) first.\n\n")
} else {
  old_ctrees <- readr::read_csv(OLD_CTREES_CSV, show_col_types = FALSE) |>
    dplyr::filter(event_id %in% mtbs_common_5070$event_id, year %in% STUDY_YEARS)

  mtbs_common_4326 <- sf::st_transform(mtbs_common_5070, 4326)
  vect_common_4326 <- terra::vect(mtbs_common_4326)
  forest_mask_30m  <- terra::rast(file.path(MASK_DIR, "nlcd2004_forestfrac_30m_ca.tif"))
  terra::terraOptions(threads = 1, progress = 0)

  CTREES_COMPARE_CSV <- file.path(OUT_DIR, "ctrees_method_comparison_ca.csv")
  done_years <- integer(0)
  if (file.exists(CTREES_COMPARE_CSV)) {
    done_years <- sort(unique(readr::read_csv(CTREES_COMPARE_CSV, show_col_types = FALSE)$year))
  }
  years_to_do <- setdiff(STUDY_YEARS, done_years)

  for (yr in years_to_do) {
    tif <- file.path(CTREES_DIR, glue("ctrees_{yr}_ca_100m.tif"))
    if (!file.exists(tif)) {
      cat(glue("  {yr} — ctrees CA 100m TIF not found, skipping\n"))
      next
    }
    t0 <- proc.time()
    r  <- terra::rast(tif)

    agb_vec <- numeric(nrow(vect_common_4326))
    for (j in seq_len(nrow(vect_common_4326))) {
      poly_4326 <- vect_common_4326[j]
      agb_vec[j] <- tryCatch({
        r_c        <- terra::crop(r, poly_4326, snap = "out")
        r_ext_5070 <- terra::project(terra::as.polygons(r_c, extent = TRUE), "EPSG:5070")
        fm_c_5070  <- terra::crop(forest_mask_30m, terra::ext(r_ext_5070) + 200, snap = "out")
        fm_c_4326  <- terra::project(fm_c_5070, r_c, method = "near")
        r_masked   <- terra::mask(r_c, fm_c_4326, maskvalues = 0)
        vals       <- terra::values(r_masked, na.rm = TRUE)
        if (length(vals) > 0L) mean(vals) else NA_real_
      }, error = function(e) {
        cat(glue("    [WARN] {vect_common_4326$event_id[j]} ({yr}) failed: {conditionMessage(e)} — set NA\n"))
        NA_real_
      })
    }

    df_yr <- data.frame(event_id = vect_common_4326$event_id, year = yr, agb_new = agb_vec)
    readr::write_csv(df_yr, CTREES_COMPARE_CSV, append = file.exists(CTREES_COMPARE_CSV))

    elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
    cat(glue("  {yr} — {nrow(df_yr)} fires — {elapsed}s (saved)\n"))

    rm(r, df_yr)
    terra::tmpFiles(remove = TRUE)
    gc(verbose = FALSE, full = TRUE)
  }

  new_ctrees <- readr::read_csv(CTREES_COMPARE_CSV, show_col_types = FALSE)
  ctrees_diff <- dplyr::inner_join(old_ctrees, new_ctrees, by = c("event_id", "year")) |>
    dplyr::mutate(
      diff     = agb_new - agb_mean_mgha,
      pct_diff = 100 * diff / agb_mean_mgha
    )
  readr::write_csv(ctrees_diff, file.path(OUT_DIR, "ctrees_method_comparison_ca_diff.csv"))

  cat(glue("\n  ctrees old-vs-new: n={nrow(ctrees_diff)} pairs | ",
           "median % diff = {round(median(ctrees_diff$pct_diff, na.rm=TRUE), 4)}% | ",
           "mean abs % diff = {round(mean(abs(ctrees_diff$pct_diff), na.rm=TRUE), 4)}% | ",
           "max abs % diff = {round(max(abs(ctrees_diff$pct_diff), na.rm=TRUE), 4)}%\n\n"))
}

# ── 7. Console summary ────────────────────────────────────────────────────────
cat("=== Validation complete ===\n")
cat("Outputs written to:", OUT_DIR, "\n")
cat("Render analysis/west_pipeline_sanity_check.qmd to view the full report.\n")
