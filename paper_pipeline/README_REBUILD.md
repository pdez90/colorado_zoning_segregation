# paper_pipeline — upstream inputs for the Colorado zoning scripts

Builds, from public data, every input that `60_co_setup.R` – `88_*.R` read.

    bash paper_pipeline/run_all.sh all      # new machine: geography, LODES, segregation panel, covariates, SLD, then 61–64
    bash paper_pipeline/run_all.sh fresh    # rebuild the SLD table, refit 63–87, influence checks, number dump
    bash paper_pipeline/run_all.sh refit    # as `fresh` without rebuilding the SLD table

| script | produces |
|---|---|
| `00_setup_and_functions.R`, `30_p2_setup.R`, `50_p3_setup.R` | conventions, paths, the Reardon–O'Sullivan spatial dissimilarity index |
| `10_geography.R` | 2020 block and tract centroids (EPSG:5070, km), tract land area |
| `20_lodes_blocks.R` | LODES 8 RAC/WAC block counts and tract OD flows, 2011–2023 |
| `32_segregation_panel.R` | tract-year segregation panel (White–Black, Hispanic–non-Hispanic, education; β = 0.25/0.5/1, aspatial) |
| `35_covariates.R` | LODES covariates; ACS per-capita income (B19301), pre-2021 releases bridged to 2020 tracts |
| `53_sld.R` | EPA Smart Location Database v3 → 2020 tracts. Primary: each 2020 block is assigned to the 2010 block group containing its internal point; tract value = POP20-weighted mean; transit no-service code read as zero. Seven sensitivity constructions are written alongside |
| `46_land_area_robustness.R` | tract-size robustness of the preferred specification |
| `45_manuscript_number_audit.R` | every numeric output in one file, `output/models/ALL_CURRENT_VALUES.csv` |
| `90_`, `91_validate_*.R` | compare a run against a reference copy of the outputs |

Paths assume the project lives at `~/Downloads/LODES` with the case-study scripts in `~/Downloads/LODES/Colorado`; `patch_60_paths.R` points `60_co_setup.R` at the zoning shapefile and TIGER folder on a new machine. ACS pulls need a Census API key. Every step caches; delete an output to force its rebuild.
