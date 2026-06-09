# =============================================================================
# 01_create_emapr_100m_tifs.R
#
# One-time preprocessing: downsample CA-clipped eMapR TIFs from 30 m to ~90 m
# (labeled "~100 m" by convention, matching ctrees native ~100 m resolution).
#
# Why: Enables fair comparison with ctrees at matched ~100 m resolution.
# 30 m × 3 = 90 m is the closest clean block-average to ctrees ~100 m (~0.000889°).
# The resulting ~90 m files (~55 MB each) are also memory-safe for terra::extract().
#
# How: terra::aggregate(fact = 3, fun = "mean") shrinks each 30 m pixel grid
# by 3× in each dimension → ~90 m output. Skip-safe: already-existing 100 m
# TIFs are never overwritten.
#
# Run once from any working directory:
#   Rscript scripts/r/01_create_emapr_100m_tifs.R
#
# Outputs: data/processed/emapr_biomass_ca/composite_YYYY_ca_100m.tif
#          (one per available source year)
#
# OUTLINE
# 1. Locate source (30 m CA) TIFs
# 2. Identify years missing 100 m versions
# 3. Aggregate and write
# 4. Report summary
# =============================================================================

library(terra)
library(here)
library(glue)

here::i_am("scripts/r/01_create_emapr_100m_tifs.R")

EMAPR_DIR <- here("data", "processed", "emapr_biomass_ca")

# ── 1. Locate source TIFs ─────────────────────────────────────────────────────
src_files <- sort(list.files(EMAPR_DIR,
                             pattern    = "^composite_\\d{4}_ca\\.tif$",
                             full.names = TRUE))
src_years <- as.integer(regmatches(basename(src_files),
                                   regexpr("\\d{4}", basename(src_files))))

cat("Source 30 m CA TIFs found:", length(src_files), "years\n")
cat("Years:", paste(src_years, collapse = ", "), "\n\n")

# ── 2. Identify years missing 100 m TIFs ──────────────────────────────────────
out_files  <- file.path(EMAPR_DIR, glue("composite_{src_years}_ca_100m.tif"))
need_build <- !file.exists(out_files)

cat(sum(!need_build), "100 m TIF(s) already exist — will skip.\n")
cat(sum(need_build),  "100 m TIF(s) to build:",
    paste(src_years[need_build], collapse = ", "), "\n\n")

if (!any(need_build)) {
  cat("Nothing to do — all 100 m TIFs are present.\n")
  quit(save = "no", status = 0)
}

# ── 3. Aggregate and write ────────────────────────────────────────────────────
# fact = 3: 30 m × 3 = 90 m (~100 m by convention); fun = "mean" preserves
# mean AGB within block. NAflag keeps existing nodata value from source raster.
t_start  <- proc.time()
counter  <- 0L

for (i in which(need_build)) {
  counter  <- counter + 1L
  yr       <- src_years[i]
  tif_in   <- src_files[i]
  tif_out  <- out_files[i]

  cat(glue("[{counter}/{sum(need_build)}] {yr} ... "))

  r_in  <- terra::rast(tif_in)
  r_out <- terra::aggregate(r_in, fact = 3, fun = "mean", na.rm = TRUE)

  terra::writeRaster(r_out, tif_out,
                     filetype  = "GTiff",
                     overwrite = FALSE,
                     gdal      = c("COMPRESS=LZW", "TILED=YES",
                                   "BLOCKXSIZE=512", "BLOCKYSIZE=512"))

  rm(r_in, r_out)
  gc(verbose = FALSE)

  elapsed <- round((proc.time() - t_start)["elapsed"], 1)
  cat(glue("done  [{elapsed}s elapsed]\n"))
}

# ── 4. Summary ────────────────────────────────────────────────────────────────
all_100m <- sort(list.files(EMAPR_DIR,
                            pattern = "^composite_\\d{4}_ca_100m\\.tif$"))
cat("\nDone. 100 m TIFs now available for",
    length(all_100m), "year(s):\n")
cat(paste(all_100m, collapse = "\n"), "\n")
