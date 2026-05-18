# =============================================================================
# 02_compare_biomass_ranges.R
# Side-by-side comparison of eMapR and ctrees biomass ranges for 2000, 2001, 2003.
# eMapR stats are loaded from pre-computed RDS files (fast).
# ctrees stats are computed on-the-fly from the 1.1 km NetCDF.
# =============================================================================

library(terra)
library(here)
library(glue)

here::i_am("analysis/03_emapr_biomass_exploration.qmd")

STUDY_YEARS <- c(2000, 2001, 2003)
CA_DIR      <- here("data", "processed", "emapr_biomass_ca")

# ── eMapR: load pre-computed stats ────────────────────────────────────────────
emapr_rows <- lapply(STUDY_YEARS, function(yr) {
  rds <- file.path(CA_DIR, glue("stats_{yr}.rds"))
  if (!file.exists(rds)) stop(glue("Missing: {rds}"))
  s  <- readRDS(rds)
  qs <- quantile(s$vals, probs = c(0.95, 0.99), na.rm = TRUE)
  data.frame(dataset = "eMapR", year = yr,
             mean   = round(s$mean,   1),
             median = round(s$median, 1),
             sd     = round(s$sd,     1),
             p95    = round(qs[1],    1),
             p99    = round(qs[2],    1),
             max    = round(s$max,    1))
})

# ── ctrees: read NetCDF, extract per-year stats ───────────────────────────────
nc_path <- here("output", "ctrees_biomass_ca_1km.nc")
if (!file.exists(nc_path)) stop(glue("Missing ctrees NetCDF: {nc_path}"))

agb <- terra::rast(nc_path)

agb_times  <- terra::time(agb)
rast_years <- if (!is.null(agb_times) && !all(is.na(agb_times))) {
  as.integer(format(agb_times, "%Y"))
} else {
  as.integer(sub(".*_(\\d{4}).*", "\\1", names(agb)))
}

ctrees_rows <- lapply(STUDY_YEARS, function(yr) {
  idx <- which(rast_years == yr)
  if (length(idx) == 0) stop(glue("Year {yr} not found in ctrees NetCDF"))
  lyr  <- agb[[idx]]
  vals <- terra::values(lyr, na.rm = TRUE)
  vals <- vals[vals > 0]
  qs   <- quantile(vals, probs = c(0.95, 0.99), na.rm = TRUE)
  data.frame(dataset = "ctrees", year = yr,
             mean   = round(mean(vals),   1),
             median = round(median(vals), 1),
             sd     = round(sd(vals),     1),
             p95    = round(qs[1],        1),
             p99    = round(qs[2],        1),
             max    = round(max(vals),    1))
})

# ── Print comparison table ────────────────────────────────────────────────────
result <- do.call(rbind, c(emapr_rows, ctrees_rows))
result  <- result[order(result$year, result$dataset), ]
rownames(result) <- NULL

cat("\n=== Biomass Range Comparison (CA pixels > 0, Mg ha⁻¹) ===\n\n")
print(result, row.names = FALSE)

cat("\nNotes:\n")
cat("  eMapR : LandTrendR on Landsat ARD, 30 m native, stats from 200k-pixel sample\n")
cat("  ctrees: ML model, ~100 m native downsampled to 1.1 km for this analysis\n")
