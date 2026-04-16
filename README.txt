## Data

Large data files are not included in this repository. To reproduce:

### MTBS Fire Perimeters
1. Download from https://www.mtbs.gov/direct-download
2. Extract to `data/raw/`
3. Expected file: `mtbs_perims_DD.shp`

### eMapR Biomass
1. Run `scripts/python/01_extract_biomass_gee.py`
2. Download from Google Drive
3. Place in `data/processed/fire_biomass_annual.csv`

### Generated Data
Run scripts in order:
```bash
Rscript scripts/r/02_process_mtbs.R
Rscript scripts/r/03_create_panel.R
```