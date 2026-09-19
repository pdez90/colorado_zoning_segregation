# When Segregation Follows Workers to Work

**Zoning, job access, and metropolitan geography in Denver.**
Replication code and outputs for a manuscript prepared for the *Journal of the American Planning Association*.

Residential segregation describes where people live, but not necessarily the segregation people experience during the workday, because commuting intervenes. This project links a harmonized zoning layer for 51 Denver-region jurisdictions (October 2023) to census-tract segregation indices built from LEHD LODES block data and to LODES origin–destination commuting flows, constructing for each tract its **workplace-location segregation exposure** — the flow-weighted segregation of the locations where its residents' jobs are. Four findings organize the paper: residential and workplace-location segregation are strongly but imperfectly related (r = .60, stable since 2011), most strongly in exclusionary, peripheral neighborhoods with poorer access to jobs; the estimated zoning relationship falls by about 70% once metropolitan position is incorporated; in an exact variance decomposition, 92% of realized exposure is accounted for by the composition of the labor market accessible from each neighborhood rather than sorting within it; and network job accessibility, by auto or by transit, adds no moderation independent of metropolitan position — estimates that suggest otherwise depend on one or two low-access fringe tracts.

## Repository contents

| Path | Contents |
|---|---|
| `60`–`88_co_*.R` | The analysis pipeline (see run order below) |
| `paper_pipeline/` | Upstream inputs built from public data (geography, LODES, segregation panel, covariates, Smart Location Database) and `run_all.sh`; see its README |
| `output/models/` | Model coefficients and paper tables (CSV) |
| `output/figures/` | Main and SI figures (PNG) |
| `diagnostics/` | QC tables written by each script |
| `CO_ZONING_DESIGN.md` | Internal design history and decision log (not the entry point — start here instead) |

## Script order and what each produces

The whole analysis runs with one command from a clone of this repository:

    export CO_ZONING_SHP=/path/to/ALL_DenverMSA_10.3.23.shp   # not redistributable, see below
    export CENSUS_API_KEY=...                                  # for the ACS income pull
    bash paper_pipeline/run_all.sh everything

Data and caches are written under `LODES_ROOT` (default `~/Downloads/LODES`). Downloads and slow steps are cached; `everything` clears the derived case-study caches before it starts, and any script that fails stops the run before the number dump, so outputs are never a mixture of old and new. Every scalar quoted in the manuscript is written to a file in `output/models/` (`p4_stats_*.csv`, `p4_manuscript_statistics.csv`, `p4_variance_allocation.csv`), and `ALL_CURRENT_VALUES.csv` collects every numeric output in one place.

| Script | Purpose | Feeds |
|---|---|---|
| `60_co_setup.R` | Configuration, paths, corrected-convention guards (sourced by all others) | — |
| `61_co_zoning_tract.R` | Zoning polygons → tract measures (dissolved areas, coverage screen) | Figs S1–S2 (via 70); zoning-layer counts |
| `62_co_od_wexp.R` | LODES OD → workplace-location segregation exposure, 2011–2023 | all exposure measures; coupling-by-year diagnostic |
| `63_co_build_panel.R` | Analysis panels (zoning + segregation + covariates) | all models |
| `64_co_models.R` | Cross-sectional model blocks and sensitivities | zoning/industrial models; change and persistence models (SI) |
| `65_co_figures.R` | Exploratory figures | Fig S5 (trends by tercile) |
| `66`–`67_co_group_*.R` | Earnings/industry/age group flows and models | Table S13; group-specific fits |
| `68_co_accessibility.R` | Metropolitan-position measures; the three-rung specification ladder | Fig 3; ladder estimates; Fig S4 |
| `69_co_robustness.R` | Levels ladder, ADU validity, effect translation, tercile descriptives | Table 1; Table S7; Fig S3 |
| `70_co_paper_figures.R` | Publication figures | Fig 3; Figs S1–S4, S6, S8 |
| `71_co_si_descriptives.R` | Where each worker group lives and works | Tables S10–S11 |
| `72_co_measure_robustness.R` | Segregation-measure sensitivity (decay β = 0.25/0.5/1.0; cutoffs 5/10/20 km) | Table S4 |
| `73_co_spatial_inference.R` | Moran's I on residuals; Conley spatial-HAC SEs (+ `fixest` cross-check) | Table S5 |
| `74_co_flow_maps.R` | Commute-flow maps: decile sheds and the exemplar pair | Fig 2; Fig S7 |
| `75_co_network_access.R` | Opportunity-set decomposition; SLD network-accessibility ladder; marginal effects | Figs S9–S10; Table S9 (upper panel); 92/8 allocation |
| `75b_co_access_influence.R` | Network-accessibility and headline interactions across SLD constructions, functional forms and leave-one-out deletions | Table S9 (lower panel); SI influence paragraph |
| `76_co_concept_figure.R` | Conceptual chain + menu-result figure | Fig 1 |
| `77_co_selfcontainment.R` | Jobs–workers dependence and commuting self-containment (tercile + jurisdiction) | Table 1 rows; inputs to Tables S15, S17 |
| `78_co_mismatch_portfolio.R` | RAC–WAC jobs–workers mismatch; destination-portfolio composition | Tables S8, S16 |
| `79_co_jurisdiction_dependence.R` | Jurisdictional boundary accounting; housing exclusion × low-wage dependence; balance-vs-matching typology | Figs S12–S13; Table S17 |
| `80_co_decentralization_network.R` | Employment decentralization 2011–2023; labor-market catchments | Figs S14–S15 |
| `81_co_zoning_flows.R` | Zoning-to-zoning flow matrix; tract-pair gravity model (PPML) | Fig S16; Table S18 |
| `82_co_excess_commuting.R` | Excess commuting: optimal worker–job assignment vs actual flows (needs `transport`) | Fig S17; Table S19 |
| `83_co_matched_access.R` | Earnings-matched accessibility (Shen competition-adjusted, per earnings band) | Table S12 |
| `84_co_alt_dimensions.R` | Position-adjusted rung for Hispanic–non-Hispanic and educational segregation | Table S6 |
| `85_co_jhbalance_literature.R` | Jobs–housing balance at four scales + nonlinearity; ladder on commute-distance outcomes | Tables S1–S3 |
| `86_co_cervero_test.R` | Cervero conjunction test: is exclusionary zoning located in employment-rich jurisdictions? | Tables S14–S15; Fig S11 |
| `87_co_workbased.R` | The work-based direction: mirrored ladder, mirrored menu, workplace-scale Cervero link | Table S20; Fig S18 |
| `88_co_manuscript_statistics.R` | Tract diameters, tract-area correlations, common-sample baseline, retention correlation | statistics quoted in text |

## Data requirements (not redistributed here)

The pipeline reads four inputs that cannot be committed for licensing or size reasons; paths are configured at the top of `60_co_setup.R`:

1. **Harmonized Denver-region zoning shapefile** (`ALL_DenverMSA_10.3.23.shp`; 100,837 district polygons, of which 100,828 have valid geometry and are used; October 2023). Not redistributable; contact the authors regarding access.
2. **Upstream segregation pipeline outputs** — the tract-level segregation panel, block/tract centroids, covariate caches and the Smart Location Database tract table. These are built from public data by `paper_pipeline/` (`bash paper_pipeline/run_all.sh everything`). Script 60 refuses to run against caches that predate the White-alone/Black-alone convention.
3. **LEHD LODES 8** (OD, RAC, WAC; public). Scripts download Colorado files automatically via `lehdr` (`version = "LODES8"`) and cache them; the committed outputs were built from the release current in September 2026 (8.4). A later LODES 8 release could change the inputs.
4. **TIGER/Line 2024 Colorado tracts and 2020 blocks** (public; downloaded by `paper_pipeline/10_geography.R` via `tigris`).
5. **EPA Smart Location Database v3** (public; downloaded by `paper_pipeline/53_sld.R`) and **ACS 5-year per-capita income, table B19301** (public; `tidycensus`, needs a Census API key).

Generated intermediates (`clean/`, `*.rds`) are `.gitignore`d: they are large and fully regenerable from the sources above.

## Software

R ≥ 4.3 with: `tidyverse` (incl. `readr`, `purrr`), `sf`, `data.table`, `lehdr`, `tigris`, `tidycensus`, `fixest`, `scales`, `patchwork`, `ggrepel`, `transport`. The spatial-inference script implements Moran's I and Conley SEs manually (no `spdep` dependency) and saves the `fixest::vcov_conley` standard error alongside for comparison; the manuscript reports the manual estimate, which is the larger of the two.

## Reproducibility notes

All committed tables and figures were produced by this pipeline run end-to-end on the inputs above. Results are deterministic given the LODES release and the October 2023 zoning snapshot, with one caveat: the per-tercile minimum commutes in the excess-commuting table come from one optimal assignment, which need not be unique; the regional minimum is unique. The segregation index is the Reardon–O'Sullivan spatial dissimilarity computed within tracts on block counts (exp(−0.5·d km) decay, 10-km truncation, EPSG:5070), with sensitivity to those choices reported in `output/models/p4_measure_sensitivity.csv`. Standard-error choices (jurisdiction clustering vs Conley spatial HAC) are compared in `output/models/p4_spatial_inference.csv`.

## Citation and license

Code is released under the MIT License (see `LICENSE`). LODES and TIGER data are public-domain U.S. Census Bureau products; the zoning compilation remains subject to its own terms. Until the paper is published, please cite this repository and: deSouza, P. (2026). *When segregation follows workers to work: Zoning, job access, and metropolitan geography in Denver.* Manuscript in preparation.
