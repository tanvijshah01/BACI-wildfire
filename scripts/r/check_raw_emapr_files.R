# =============================================================================
# check_raw_emapr_files.R
# Sanity-check the raw eMapR CONUS composites after a full-archive download.
# Standalone, re-runnable diagnostic — not part of the numbered pipeline
# (mirrors validate_west_pipeline.R's role for the extraction pipeline).
#
# For each expected year (1990-2023), checks:
#   1. File exists
#   2. File size is in a plausible range (catches truncated/failed downloads)
#   3. Raster header opens cleanly (catches files GDAL can't even read)
#   4. A small block of pixels near the raster's own center has real,
#      non-fill biomass values (never touches the whole ~28GB file — see
#      NOTES.md's terra-large-raster lesson: whole-raster reads on files
#      this size are slow/OOM-prone; this reads one small indexed block
#      instead of cropping to hardcoded coordinates, so it can't fail from
#      a wrong CRS/extent guess)
#
# Run from the project root, outside Quarto (same reasoning as the rest of
# the pipeline — terra's C++ threading looks frozen inside a Quarto chunk):
#   Rscript scripts/r/check_raw_emapr_files.R
#
# Output: console summary table; nothing written to disk.
# =============================================================================
# 1. Setup
# 2. Per-year checks
# 3. Summary report
# =============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────
library(terra)
library(here)
library(glue)

here::i_am("scripts/r/check_raw_emapr_files.R")

RAW_DIR        <- here("data", "raw", "emapr_biomass")
EXPECTED_YEARS <- 1990:2023
# 96,815 x 153,809 px, int16, uncompressed -> every complete year is this exact
# size (confirmed against 19 verified-complete downloaded years on GRIT,
# 2026-09-08). A wider 20-35 GB band was tried first but let a truncated file
# through (2009 landed at 20.1 GB, just above a 20 GB floor, containing only
# 67.6% of the expected pixel count while its header still claimed full
# dimensions) - so this is now a tight tolerance around the known-exact size
# rather than a generous nominal range.
SIZE_EXACT_BYTES <- 29783198966
SIZE_TOLERANCE_GB <- 0.05
SAMPLE_HALFWIN <- 250  # read a (2*250+1)^2 pixel block centered on the raster

cat("Checking", length(EXPECTED_YEARS), "expected years in", RAW_DIR, "\n\n")

# ── 2. Per-year checks ──────────────────────────────────────────────────────────
results <- data.frame()

for (yr in EXPECTED_YEARS) {
  path <- file.path(RAW_DIR, glue("composite_{yr}_median.tif"))
  row  <- data.frame(year = yr, exists = FALSE, size_gb = NA_real_,
                      size_ok = FALSE, opens = FALSE, dims = NA_character_,
                      sample_has_data = FALSE, status = "MISSING")

  if (file.exists(path)) {
    row$exists  <- TRUE
    row$size_gb <- round(file.size(path) / 1e9, 1)
    row$size_ok <- abs(file.size(path) - SIZE_EXACT_BYTES) <= SIZE_TOLERANCE_GB * 1e9

    r <- tryCatch(terra::rast(path), error = function(e) NULL)

    if (!is.null(r)) {
      row$opens <- TRUE
      row$dims  <- glue("{nrow(r)} x {ncol(r)}")

      # Read a small block indexed off the raster's own reported center —
      # never a hardcoded coordinate, so this can't fail from a wrong
      # CRS/extent assumption. Cheap: terra only reads this block from disk.
      center_row <- max(1, round(nrow(r) / 2) - SAMPLE_HALFWIN)
      center_col <- max(1, round(ncol(r) / 2) - SAMPLE_HALFWIN)
      n_rows     <- min(2 * SAMPLE_HALFWIN + 1, nrow(r) - center_row + 1)
      n_cols     <- min(2 * SAMPLE_HALFWIN + 1, ncol(r) - center_col + 1)

      vals <- tryCatch(
        terra::values(r, row = center_row, nrows = n_rows,
                       col = center_col, ncols = n_cols),
        error = function(e) NULL
      )
      if (!is.null(vals)) {
        vals <- vals[!is.na(vals)]
        row$sample_has_data <- length(vals) > 0 && any(vals > 0)
      }
    }

    row$status <- if (row$size_ok && row$opens && row$sample_has_data) {
      "OK"
    } else {
      "SUSPECT"
    }
  }

  results <- rbind(results, row)
  cat(glue(
    "{yr}: {row$status}",
    "{if (row$exists) glue(' ({row$size_gb} GB, dims {row$dims})') else ''}\n"
  ))
}

# ── 3. Summary ────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat("OK:     ", sum(results$status == "OK"), "/", length(EXPECTED_YEARS), "\n")
cat("SUSPECT:", sum(results$status == "SUSPECT"), "\n")
cat("MISSING:", sum(results$status == "MISSING"), "\n\n")

if (any(results$status != "OK")) {
  cat("Years needing attention:\n")
  print(results[results$status != "OK",
                c("year", "status", "exists", "size_gb", "opens", "sample_has_data")])
  cat("\nA SUSPECT file (exists but fails a check) most likely means an\n")
  cat("interrupted/truncated download — delete it and re-run the download\n")
  cat("loop for that year; it will skip everything already valid.\n")
}
