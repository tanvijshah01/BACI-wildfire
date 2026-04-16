# Project Notes — Wildfire Biomass Recovery

Working notes on decisions, findings, and open questions. Add entries in reverse chronological order (newest at top).

---

## 2026-04-15 — Technical Decisions: `01_exploration.qmd`

### `sf_use_s2(FALSE)` — use GEOS instead of s2 for spatial operations
MTBS fire perimeter polygons contain self-touching edges that the s2 spherical geometry engine rejects as invalid, causing `st_join()` to error. Switching to GEOS (planar geometry) resolves this. For coterminous US data the planar approximation is accurate enough and `st_make_valid()` alone is not sufficient to satisfy s2's stricter validity rules. **Apply this at the top of any R script that does spatial joins on MTBS data.**

### `here::i_am("analysis/01_exploration.qmd")` — anchor project root for Quarto docs
Quarto renders documents with the working directory set to the document's own folder (`analysis/`), so `here("data/...")` resolves to `analysis/data/...` — the wrong place. Calling `here::i_am()` with the document's path relative to the project root anchors `here()` correctly. **Every Quarto document in a subdirectory must include this call.**

### Auto-install block for `tigris` and `ggspatial`
These two packages were not pre-installed in the R environment. The QMD includes a `.install_if_missing()` helper at the top of the setup chunk so the document installs them automatically on first render. If the environment ever gets rebuilt, these two packages must be added.

---

## Decision Log

Use this section to record design decisions and their rationale once resolved.

| Date | Decision | Choice | Rationale |
|---|---|---|---|
| — | Minimum fire size | — | — |
| — | Severity inclusion | — | — |
| — | Treatment variable (categorical vs. dNBR) | — | — |
| — | Study region / ecoregions | — | — |
| — | Control site strategy | — | — |
| — | Fire complex handling | — | — |
| — | Spatial buffer size | — | — |

---

## Data Quality Findings

*Record any anomalies, missing fields, or unexpected values found during EDA.*

### MTBS Perimeters (`mtbs_perims_DD.shp`)
- [ ] Field names confirmed
- [ ] CRS confirmed
- [ ] NAs / duplicates checked
- Notes:

---

## Open Questions

Questions that have not yet been resolved and need investigation or a decision.

- Which ecoregion classification to use (EPA Level III vs. Bailey's)?
- Should fire complexes be treated as a single large fire or excluded?
- What is the actual GEE asset path for eMapR biomass? (`projects/eMapR/biomass` — verify in GEE catalog)
- Minimum pre-fire years needed for parallel trends test with Callaway-Sant'Anna?
- How to handle fires that re-burn the same site in a different year?

---

## Literature Notes

*Key findings and methodological details from papers relevant to design decisions.*

### Callaway & Sant'Anna (2021)
- Estimator: `att_gt()` for group-time ATTs, `aggte()` for aggregation
- Clean controls: "never treated" or "not yet treated" units
- Does not require balanced panel
- Pre-treatment parallel trends is testable with their placebo approach

### Bright et al. (2019)
- Predictive (not causal) — R² > 0.7 with random forest
- Useful benchmark for biomass signal magnitude

### Ilangakoon et al. (2026)
- GAM with space-for-time substitution — lacks formal causal identification
- Our study directly addresses this gap

---

## GEE / Python Notes

*Notes on Google Earth Engine setup and extraction issues.*

- GEE authentication: `earthengine authenticate` in terminal before running Python script
- eMapR biomass asset path needs verification in GEE catalog before running extraction
- LandTrendR reference: https://emapr.github.io/LT-GEE/
- LandTrendR key parameters (defaults): `maxSegments`, `spikeThreshold=0.9`, `recoveryThreshold=0.25`, `pvalThreshold=0.1`, `minObservationsNeeded=6`

---

## Meeting / Advisor Notes

*Notes from advisor meetings or collaborator discussions.*

---
