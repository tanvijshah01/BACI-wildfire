# EDA Plan — Wildfire Biomass Recovery

**Phase:** Exploration → West-wide expansion
**Last updated:** 2026-09-20
**Goal:** Understand the MTBS, eMapR, and ctrees data well enough to make the design decisions in
`NOTES.md` (→ Design decisions & open questions) before building the analysis panel.

**Status legend:** ✅ done · 🔄 in progress / partly done · ⬜ not started · ⛔ stale (superseded)

Where this fits: `CLAUDE.md` = project orientation and pipeline status · `DATA_DOWNLOAD_GUIDE.md` = how to
run the pipeline · `NOTES.md` = decisions, gotchas, history · **this file = what EDA we're doing and how far
along each piece is.** Analysis documents live in `analysis/`; figures land in `figures/` (gitignored).

---

## 1. Core EDA questions

| # | Question | Status | Where answered / next step |
|---|---|---|---|
| 1 | How many MTBS fires in 2000–2023, Western US? | ✅ | `01_mtbs_exploration.qmd` §3 (row count printed at each filter step) |
| 2 | Distribution of fire sizes? | ✅ | `01` §3.6 (log₁₀ histograms, all vs. Extended) |
| 3 | Severity classes and their distribution? | ✅ | `01` §3.7 (BSP CBI dNBR/RdNBR). Note: the MTBS perimeter shapefile holds dNBR *threshold values*, not per-class area, so % low/moderate/high per fire is not available from it |
| 4 | Temporal trends in fire occurrence? | ✅ | `01` §3.2, §3.5 |
| 5 | Spatial clustering of fires? | ⬜ | Needed for the SUTVA / spatial-buffer decision |
| 6 | Fires per year (treatment timing)? | ✅ | `01` §3.2; per-cohort view in `biomass_within_fires.qmd` §2 |
| 7 | Which years have enough pre-fire data for parallel-trends testing? | ⬜ | Data-availability facts: eMapR starts 1990 (≥10 pre-fire years for 2000 fires), ctrees starts 2000 (none for 2000 fires). Quantify per cohort once the panel exists |
| 8 | Do ctrees and eMapR agree within fires? | 🔄 | `biomass_within_fires.qmd` §5, §7 — residual eMapR bias not root-caused |
| 9 | Is MTBS Initial vs. Extended assessment type spatially/seasonally biased? | 🔄 | `mtbs_assessment_comparison.qmd` — unresolved |

---

## 2. Analysis documents

| Document | Purpose | Status | Notes / open items |
|---|---|---|---|
| `01_mtbs_exploration.qmd` | MTBS perimeters, assessment types, sizes, timing, CBI severity | ✅ §1–3 · ⛔ §4 | §4 (GEE NBR extraction, 4.1/4.2) belongs to the archived NBR approach — remove or replace with a pointer to `biomass_within_fires.qmd` §6 (see `NOTES.md` 2026-05-12) |
| `02_ctrees_biomass_exploration.qmd` | CA ctrees AGB: annual trend, maps (2000/2020/2025), before/after fire event study | ✅ (CA) | Uses `03_download_ctrees_ca.py` outputs; no West version |
| `03_emapr_biomass_exploration.qmd` | CA eMapR AGB for 2000/2001/2003: distributions, year comparison, maps | ✅ (CA) | Reads retired CA crops (`emapr_biomass_ca/`) |
| `04_data_summary.qmd` | Side-by-side summary of MTBS, ctrees, eMapR (CA); dataset comparison table; §7 design decisions (assessment types, outcome variable, control-unit definition & SUTVA) | ✅ (CA) | Depends on the retired `03_prepare_forest_mask.R` mask. §7 decisions should be reconciled with `NOTES.md` |
| `biomass_within_fires.qmd` | **Current.** ctrees vs. eMapR within forested fire pixels, CA cohorts; masking progression; before/after event study (§6); design notes (§7) | 🔄 | Fed by `07`/`08` (`STATES_TO_RUN` = CA, `STUDY_YEARS` 2005–2010). Open: eMapR residual bias (§7); extend to West + all years |
| `biomass_within_fires_old.qmd` | Retired CA-only forest-mask version | ⛔ | Kept only as the validation baseline |
| `mtbs_assessment_comparison.qmd` | Initial vs. Extended: fire counts by state, ignition month, % forest cover, forest-cover-threshold sweep | 🔄 | Forest-cover cache currently CA only (`06`); widen to West |
| `west_pipeline_sanity_check.qmd` | Regression report for the 05–08 multi-state rewrite (dedup, mask validity, eMapR/ctrees old-vs-new) | ✅ | Re-run after any change to 05–08 (`Rscript scripts/r/validate_west_pipeline.R`, then render) |

---

## 3. Next EDA steps

1. **West-wide within-fire extraction** — after the ctrees West download and West eMapR crops finish
   (`CLAUDE.md` → Current Status), widen `STATES_TO_RUN` / `STUDY_YEARS` and refresh `biomass_within_fires.qmd`
   and `mtbs_assessment_comparison.qmd` for all states and years.
2. **Explain the eMapR-vs-ctrees residual bias** and re-check the border-fire count gap on real multi-state
   output (`NOTES.md` 2026-08-30).
3. **Pre-fire data sufficiency** (Q7) and **spatial clustering** (Q5), feeding the min-pre-years, fire-complex,
   and buffer decisions.
4. **Control-site candidates** — never-burned pixels in the same ecoregion (`data/processed/control_pixels/`).
5. **Tidy `01` §4 and `04` §7** per the table above.

---

## 4. Record: completed MTBS EDA (`01_mtbs_exploration.qmd`)

**Filter steps** (row count printed after each): `incid_type == "Wildfire"` → `year` 2000–2023 → spatial
join to the 11 Western states (`tigris::states()`). **Sanity checks:** required fields present, no duplicate
`event_id`, valid geometries (`st_make_valid()`), dNBR thresholds not all in 0–100 (would mean percentages),
sentinel values −9999/9999 flagged, all 11 target states present. `sf_use_s2(FALSE)` is required — see
`NOTES.md` → Technical gotchas.

**Figures produced** (`figures/`, gitignored — regenerate by rendering the qmd):

| Section | Figure | Content |
|---|---|---|
| 3.1 | `fire_locations_map.png` | All fires, colored by assessment type (EPSG:5070) |
| 3.2 | `fires_per_year.png` | Annual wildfire count with linear trend |
| 3.3 | `asmnt_type_summary.png` | Fire count by assessment type |
| 3.4 | `fire_locations_extended_only.png` | Extended-assessment fires only, with n and % of total |
| 3.5 | `fires_per_year_comparison.png` | All vs. Extended, shared y-axis |
| 3.6 | `fire_size_comparison.png` | log₁₀ size histograms, all vs. Extended, medians marked |
| 3.7 | `severity_dnbr_map.png`, `severity_dnbr_histogram.png`, `severity_rdnbr_histogram.png` | BSP CBI field-plot dNBR/RdNBR (`map_prog == "CBI"`, sentinel 9999 and `|preNBR| ≤ 10` excluded) |

**Notes:** No Extended-only comparison for CBI severity — BSP CBI (field-assessed) and MTBS (satellite-assessed)
are largely different fires, so joining on `event_id` returns near-zero matches. RdNBR = dNBR / √|preNBR / 1000|
(Miller & Thode 2007) and is the leading candidate for the continuous treatment variable (`NOTES.md`).
Methodology and citations: `NOTES.md` → Literature notes.

---

## 5. Figure QA (apply to every plot before saving)

- 300 dpi; check pixel dimensions match the intended inches
- Legend has a descriptive title (not a raw column name), text ≥ 10 pt, no overlap with data
- Colorblind-safe palette (`viridis`, `RColorBrewer "Set2"`, or `MetBrewer`)
- Maps: scale bar, north arrow, state boundaries for reference
- Time axes span the full year range with readable ticks; event-study plots mark year 0 with a dashed line
  and label the x-axis "Years relative to fire"
- Paired panels (all vs. Extended, CA vs. West) share axis scales
