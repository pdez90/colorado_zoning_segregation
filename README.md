# Where Segregation Follows Workers to Work

**Commuting geography, zoning, and metropolitan spatial structure in the Denver region.**
Replication code and outputs for a manuscript prepared for the *Journal of the American Planning Association* (draft: `Paper4_JAPA_draft_v2.docx`).

Residential segregation is not necessarily the segregation workers experience during the day, because commuting intervenes. This project links a harmonized zoning layer for 51 Denver-region jurisdictions (October 2023) to census-tract segregation indices built from LEHD LODES block data and to LODES origin–destination commuting flows, constructing for each tract its **workplace-location segregation exposure** — the flow-weighted segregation of the locations where its residents' jobs are. The central result is a decomposition: estimated naively, exclusionary zoning appears to moderate the coupling between residential and workplace-location segregation powerfully (+0.33 SD), but allowing metropolitan position to moderate the same relationship reduces the estimated zoning moderation by 60–72% across all substantively interpretable zoning measures; a position-adjusted zoning association of roughly one-third the baseline remains, stable across alternative spatial parameterizations of the segregation measure and robust to spatial-HAC inference.

## Repository contents

| Path | Contents |
|---|---|
| `60`–`73_co_*.R` | The analysis pipeline (see run order below) |
| `output/models/` | Model coefficients and paper tables (CSV) |
| `output/figures/` | Main and SI figures (PNG) |
| `diagnostics/` | QC tables written by each script |
| `CO_ZONING_DESIGN.md` | Internal design history and decision log (not the entry point — start here instead) |
| `Paper4_JAPA_draft_v2.docx` | Current manuscript draft |

## Script order and what each produces

Scripts are checkpointed: completed steps skip themselves on re-run.

| Script | Purpose | Feeds |
|---|---|---|
| `60_co_setup.R` | Configuration, paths, corrected-convention guards (sourced by all others) | — |
| `61_co_zoning_tract.R` | Zoning polygons → tract measures (dissolved areas, coverage screen) | Fig 1 |
| `62_co_od_wexp.R` | LODES OD → workplace-location segregation exposure, 2011–2023 | Fig 2; coupling-by-year diagnostic |
| `63_co_build_panel.R` | Analysis panels (zoning + segregation + covariates) | all models |
| `64_co_models.R` | Cross-sectional model blocks and sensitivities | — |
| `65_co_figures.R` | Exploratory figures | Fig S4 (trends by tercile) |
| `66`–`67_co_group_*.R` | Earnings/industry/age group flows and models | "Who bears the geography" |
| `68_co_accessibility.R` | Metropolitan-position measures; the three-rung specification ladder | Fig 3; Table of ladder estimates |
| `69_co_robustness.R` | Levels ladder, ADU validity, effect translation, tercile descriptives | Table 1 |
| `70_co_paper_figures.R` | Publication figures | Figs 1–4, S1–S3 |
| `71_co_si_descriptives.R` | Where each worker group lives and works | Tables S1–S2 |
| `72_co_measure_robustness.R` | Segregation-measure sensitivity (decay β = 0.25/0.5/1.0; cutoffs 5/10/20 km) | Table S3 |
| `73_co_spatial_inference.R` | Moran's I on residuals; Conley spatial-HAC SEs (+ `fixest` cross-check) | Table S4 |

## Data requirements (not redistributed here)

The pipeline reads four inputs that cannot be committed for licensing or size reasons; paths are configured at the top of `60_co_setup.R`:

1. **Harmonized Denver-region zoning shapefile** (`ALL_DenverMSA_10.3.23.shp`; 100,837 district polygons, October 2023). Not redistributable; contact the authors regarding access.
2. **Upstream segregation pipeline outputs** — the tract-level segregation panel, block/tract centroids, and covariate caches produced by the companion LODES replication pipeline (Tan & deSouza 2026 lineage), including the corrected covariate cache `p2_tract_lodes_cov_panel_v3.rds`. Script 60 refuses to run against pre-correction caches.
3. **LEHD LODES 8.4** (OD, RAC, WAC; public). Scripts download Colorado files automatically via `lehdr` and cache them.
4. **TIGER/Line 2024 Colorado tracts** (public; auto-downloaded via `tigris` if not found locally).

Generated intermediates (`clean/`, `*.rds`) are `.gitignore`d: they are large and fully regenerable from the sources above.

## Software

R ≥ 4.3 with: `tidyverse`, `sf`, `data.table`, `seg`, `lehdr`, `fixest`, `patchwork`, `ggrepel` (optionally `tigris`). The spatial-inference script implements Moran's I and Conley SEs manually (no `spdep` dependency) and soft-checks the Conley estimator against `fixest::vcov_conley` when available.

## Reproducibility notes

All committed tables and figures were produced by this pipeline run end-to-end on the inputs above. Results are deterministic given the LODES vintage (8.4) and the October 2023 zoning snapshot. The segregation index is the Reardon–O'Sullivan spatial dissimilarity computed within tracts on block counts (exp(−0.5·d km) decay, 10-km truncation, EPSG:5070), with sensitivity to those choices reported in `output/models/p4_measure_sensitivity.csv`. Standard-error choices (jurisdiction clustering vs Conley spatial HAC) are compared in `output/models/p4_spatial_inference.csv`.

## Citation and license

Code is released under the MIT License (see `LICENSE`). LODES and TIGER data are public-domain U.S. Census Bureau products; the zoning compilation remains subject to its own terms. Until the paper is published, please cite this repository and: deSouza, P. (2026). *Where segregation follows workers to work: Commuting geography, zoning, and metropolitan spatial structure in the Denver region.* Manuscript in preparation.
