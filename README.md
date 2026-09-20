# Wildfire Biomass Recovery — Causal Inference Study (BACI)

UCSB MESM research project estimating the causal effect of wildfire severity on forest biomass
recovery, using a Before-After-Control-Intervention (BACI) / staggered difference-in-differences design
(Callaway & Sant'Anna 2021) on MTBS fire perimeters (2000–2023, 11 contiguous Western US states) and two
annual Landsat-based aboveground-biomass products, **eMapR** (1990–2023) and **ctrees** (2000–2025).

**Status:** exploratory data analysis and West-wide data pipeline build-out. See
[`CLAUDE.md`](CLAUDE.md) → *Current Status* for what's done and what's next.

## Getting started

Large data files are not in this repository (`data/raw/`, `data/processed/` are gitignored). Everything
needed to reproduce them is in **[`DATA_DOWNLOAD_GUIDE.md`](DATA_DOWNLOAD_GUIDE.md)**:

| Data | Source | Where to look |
|---|---|---|
| MTBS fire perimeters | [mtbs.gov/direct-download](https://www.mtbs.gov/direct-download) → extract into `data/raw/mtbs/` | expected file `mtbs_perimeter_data/mtbs_perims_DD.shp` |
| eMapR biomass | anonymous FTP (`islay.ceoas.oregonstate.edu`), ~27.7 GB/year | guide Part 2 |
| ctrees biomass | arraylake zarr store (needs `arraylake auth login`) | guide Part 3 |
| Forest mask | NLCD 2004, auto-downloaded via `FedData` | guide Part 4 |

## Pipeline (R, from the project root, outside Quarto)

```r
Rscript scripts/r/00_crop_emapr_to_west.R              # crop raw eMapR to the 11-state West
python  scripts/python/04_download_ctrees_west.py      # ctrees West download (multi-hour; run in tmux)
Rscript scripts/r/05_prepare_forest_masks_west.R       # per-state NLCD 2004 forest masks
Rscript scripts/r/06_extract_pct_forest_within_fires.R # % forest per fire perimeter
Rscript scripts/r/07_extract_emapr_within_fires_new.R  # eMapR AGB within forested fire pixels
Rscript scripts/r/08_extract_ctrees_within_fires_new.R # ctrees AGB within forested fire pixels
```

Outputs feed the Quarto documents in `analysis/`. Scripts are skip-safe and resume where they left off.

## Documentation map

| File | What it covers |
|---|---|
| [`CLAUDE.md`](CLAUDE.md) | Project orientation: study design, directory layout, pipeline, current status, coding conventions |
| [`DATA_DOWNLOAD_GUIDE.md`](DATA_DOWNLOAD_GUIDE.md) | How to set up an environment and run download → crop → extract |
| [`NOTES.md`](NOTES.md) | Technical gotchas, design decisions and open questions, dated findings log, literature notes |
| [`EDA_PLAN.md`](EDA_PLAN.md) | EDA questions and per-document status |
