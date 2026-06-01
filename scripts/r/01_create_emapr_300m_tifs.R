# =============================================================================
# 01_create_emapr_300m_tifs.R
#
# One-time preprocessing: downsample CA-clipped eMapR TIFs from 30 m to 300 m.
#
# Why: terra::extract() on 500 MB 30 m files causes std::bad_alloc when called
# in a loop inside Quarto. 300 m TIFs (~12 MB each) allow fast, memory-safe
# polygon extraction with negligible accuracy loss for fire polygons ≥1,000 ac.
#
# How: terra::aggregate(fact = 10, fun = "mean") shrinks each 30 m pixel grid
# by 10× in each dimension → 300 m output. Skip-safe: already-existing 300 m
# TIFs are never overwritten.
#
# Run once from any working directory:
#   Rscript scripts/r/01_create_emapr_300m_tifs.R
#
# Outputs: data/processed/emapr_biomass_ca/composite_YYYY_ca_300m.tif
#          (one per available source year)
#
# OUTLINE
# 1. Locate source (30 m CA) TIFs
# 2. Identify years missing 300 m versions
# 3. Aggregate and write
# 4. Report summary
# =============================================================================

library(terra)
library(here)
library(glue)

here::i_am("scripts/r/01_create_emapr_300m_tifs.R")

EMAPR_DIR <- here("data", "processed", "emapr_biomass_ca")

# ── 1. Locate source TIFs ─────────────────────────────────────────────────────
src_files <- sort(list.files(EMAPR_DIR,
                             pattern  = "^composite_\\d{4}_ca\\.tif$",
                             full.names = TRUE))
src_years <- as.integer(regmatches(basename(src_files),
                                   regexpr("\\d{4}", basename(src_files))))

cat("Source 30 m CA TIFs found:", length(src_files), "years\n")
cat("Years:", paste(src_years, collapse = ", "), "\n\n")

# ── 2. Identify years missing 300 m TIFs ──────────────────────────────────────
out_files  <- file.path(EMAPR_DIR, glue("composite_{src_years}_ca_300m.tif"))
need_build <- !file.exists(out_files)

cat(sum(!need_build), "300 m TIF(s) already exist — will skip.\n")
cat(sum(need_build),  "300 m TIF(s) to build:",
    paste(src_years[need_build], collapse = ", "), "\n\n")

if (!any(need_build)) {
  cat("Nothing to do — all 300 m TIFs are present.\n")
  quit(save = "no", status = 0)
}

# ── 3. Aggregate and write ────────────────────────────────────────────────────
# fact = 10: 30 m × 10 = 300 m; fun = "mean" preserves mean AGB within block.
# NAflag keeps existing nodata value from the source raster.
t_start <- proc.time()

for (i in which(need_build)) {
  yr       <- src_years[i]
  tif_in   <- src_files[i]
  tif_out  <- out_files[i]

  cat(glue("[{which(need_build) == i |> which()}",
           "/{sum(need_build)}] {yr} ... "))

  r_in  <- terra::rast(tif_in)
  r_out <- terra::aggregate(r_in, fact = 10, fun = "mean", na.rm = TRUE)

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
all_300m <- sort(list.files(EMAPR_DIR,
                            pattern = "^composite_\\d{4}_ca_300m\\.tif$"))
cat("\nDone. 300 m TIFs now available for",
    length(all_300m), "year(s):\n")
cat(paste(all_300m, collapse = "\n"), "\n")
