# ==============================================================================
# 63_co_build_panel.R      [COLORADO ZONING case study, step 3]
# Assemble the analysis datasets:
#
#   HOME panel  (tract-year, 2011-2023): wexp + commute measures (62) +
#     residential/workplace segregation (32) + ZONING measures (61, constant
#     within tract) + Paper-2 tract covariates (35) + ids + scope flags.
#
#   WORK panel  (tract-year): wres (62) + workplace segregation + zoning of
#     the WORK tract + industry mix covariates.
#
#   CHANGE table (tract, 2011 -> 2023): first-vs-last-year changes in
#     residential segregation and wexp, for the "does 2023 zoning predict
#     2011->2023 trajectories" block (zoning changes slowly; stated
#     explicitly as descriptive, not causal).
#
# SLD transit measures join IF the national 53 output exists (optional --
# lets 64 pit zoning against transit access as competing moderators).
#
# Output: output/co_analysis_home_panel.rds / .csv.gz
#         output/co_analysis_work_panel.rds
#         output/co_analysis_change_2011_2023.rds
# ==============================================================================

source("60_co_setup.R")

## ---- 1. load pieces ----------------------------------------------------------
home   <- readRDS(file.path(DIR_CO_CLEAN, "co_wexp_home_panel.rds"))
work   <- readRDS(file.path(DIR_CO_CLEAN, "co_wres_work_panel.rds"))
zoning <- readRDS(file.path(DIR_CO_CLEAN, "co_tract_zoning.rds"))
seg    <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds"))
lodes_cov <- readRDS(P4_LODES_COV_FILE)   # v3: corrected White (CR01) + counts
stopifnot(all(c("pct_white_rac", "pct_white_wac", "n_white_rac",
                "pct_black_rac") %in% names(lodes_cov)))  # 36's guard, mirrored
income    <- readRDS(file.path(DIR_CLEAN, "p2_tract_income_panel.rds"))
aland     <- readRDS(file.path(DIR_CLEAN, "tract_aland_2020.rds"))

sld_file <- file.path(DIR_CLEAN, "p3_tract_sld.rds")
sld <- if (file.exists(sld_file)) readRDS(sld_file) else NULL
if (is.null(sld)) message("No SLD file yet (paper_pipeline/53_sld.R) -- ",
                          "zoning-vs-transit comparisons will be skipped.")

seg_rac <- seg |> select(tract_id, year, starts_with("d_") & ends_with(
  c("_rac_half", "_rac_aspatial")))
seg_wac <- seg |> select(tract_id, year, starts_with("d_") & ends_with(
  c("_wac_half", "_wac_aspatial")))

zon_cols <- c("cover_zoned", "zoned_km2", "res_zoned_km2",
              paste0("pct_", names(CO_ZONE_GROUPS)),
              "pct_reslow_of_res", "pct_adu_res", "pct_grouphome_res",
              "pct_job_zone", "zoning_entropy", "n_jurisd", "jurisd_main")
zoning_j <- zoning |> select(tract_id, all_of(zon_cols))

## ---- 2. HOME panel -----------------------------------------------------------
dat_home <- home |>
  rename(CBSAFP = cbsa) |>
  inner_join(zoning_j,  by = "tract_id") |>   # zoning-covered tracts only
  left_join(seg_rac,    by = c("tract_id", "year")) |>
  left_join(seg_wac,    by = c("tract_id", "year")) |>
  left_join(lodes_cov,  by = c("tract_id", "year")) |>
  left_join(income,     by = c("tract_id", "year")) |>
  left_join(aland,      by = "tract_id") |>
  mutate(
    county_fips = substr(tract_id, 1, 5),
    state_fips  = substr(tract_id, 1, 2),
    pandemic    = year %in% P2_PANDEMIC,
    worker_density_rac = ifelse(aland_km2 > 0, workers_rac / aland_km2,
                                NA_real_),
    log_worker_density_rac = log(pmax(worker_density_rac, 1e-6)),
    income_percapita_k     = income_percapita / 1000,
    income_percapita_k_sq  = income_percapita_k^2,
    # jobs available locally per resident worker (jobs-housing balance;
    # zoning that forbids jobs should show up here first)
    jobs_housing_ratio = ifelse(workers_rac > 0, workers_wac / workers_rac,
                                NA_real_),
    in_scope   = cover_zoned >= CO_MIN_ZONED_COVER,
    denver_msa = CBSAFP == CO_CBSA_MAIN)
if (!is.null(sld)) dat_home <- dat_home |> left_join(sld, by = "tract_id")
if ("sld_D5BR" %in% names(dat_home))
  dat_home <- dat_home |>
    mutate(log_sld_D5BR = ifelse(sld_D5BR > 0, log(sld_D5BR), NA_real_))

saveRDS(dat_home, file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
data.table::fwrite(dat_home,
                   file.path(DIR_CO_OUT, "co_analysis_home_panel.csv.gz"))
message("Saved HOME panel: ", nrow(dat_home), " tract-years, ",
        n_distinct(dat_home$tract_id), " tracts (",
        sum(dat_home$in_scope[dat_home$year == CO_ANCHOR_YEAR]),
        " in-scope in ", CO_ANCHOR_YEAR, ")")

## ---- 3. WORK panel -----------------------------------------------------------
dat_work <- work |>
  rename(CBSAFP = cbsa) |>
  inner_join(zoning_j, by = "tract_id") |>
  left_join(seg_wac,   by = c("tract_id", "year")) |>
  left_join(lodes_cov, by = c("tract_id", "year")) |>
  left_join(aland,     by = "tract_id") |>
  mutate(
    county_fips = substr(tract_id, 1, 5),
    state_fips  = substr(tract_id, 1, 2),
    pandemic    = year %in% P2_PANDEMIC,
    worker_density_wac = ifelse(aland_km2 > 0, workers_wac / aland_km2,
                                NA_real_),
    log_worker_density_wac = log(pmax(worker_density_wac, 1e-6)),
    in_scope   = cover_zoned >= CO_MIN_ZONED_COVER,
    denver_msa = CBSAFP == CO_CBSA_MAIN)

saveRDS(dat_work, file.path(DIR_CO_OUT, "co_analysis_work_panel.rds"))
message("Saved WORK panel: ", nrow(dat_work), " tract-years")

## ---- 4. CHANGE table 2011 -> 2023 --------------------------------------------
y0 <- min(CO_YEARS); y1 <- CO_ANCHOR_YEAR
chg <- dat_home |>
  filter(year %in% c(y0, y1),
         !is.na(d_whiteblack_rac_half), !is.na(wexp_whiteblack_wac_half)) |>
  select(tract_id, CBSAFP, county_fips, year, in_scope, denver_msa,
         all_of(setdiff(zon_cols, "jurisd_main")), jurisd_main,
         res_seg = d_whiteblack_rac_half,
         wexp    = wexp_whiteblack_wac_half,
         n_commuters, pct_black_rac, pct_lowincome_rac, income_percapita_k,
         log_worker_density_rac) |>
  pivot_wider(names_from = year,
              values_from = c(res_seg, wexp, n_commuters, pct_black_rac,
                              pct_lowincome_rac, income_percapita_k,
                              log_worker_density_rac),
              names_sep = "_") |>
  mutate(d_res_seg = .data[[paste0("res_seg_", y1)]] -
                     .data[[paste0("res_seg_", y0)]],
         d_wexp    = .data[[paste0("wexp_", y1)]] -
                     .data[[paste0("wexp_", y0)]]) |>
  filter(!is.na(d_res_seg), !is.na(d_wexp))
saveRDS(chg, file.path(DIR_CO_OUT, "co_analysis_change_2011_2023.rds"))
message("Saved CHANGE table: ", nrow(chg), " tracts with both endpoints")

## ---- DIAGNOSTICS -------------------------------------------------------------
write_codiag(
  dat_home |> filter(year == CO_ANCHOR_YEAR) |>
    group_by(CBSAFP, in_scope) |>
    summarise(n_tracts = n(),
              med_commuters = median(n_commuters),
              mean_res_low = round(mean(pct_res_low, na.rm = TRUE), 1),
              .groups = "drop"),
  "63_scope_by_cbsa")

key_cols <- intersect(
  c("wexp_whiteblack_wac_half", "d_whiteblack_rac_half", "pct_res_low",
    "pct_reslow_of_res", "pct_adu_res", "zoning_entropy", "pct_job_zone",
    "mean_dist_km", "pct_black_rac", "income_percapita"),
  names(dat_home))
write_codiag(
  dat_home |> filter(in_scope) |> group_by(year) |>
    summarise(n_tracts = n(),
              across(all_of(key_cols), ~ round(100 * mean(is.na(.x)), 1),
                     .names = "pctna_{.col}")),
  "63_home_panel_completeness_by_year")

nonfin <- sapply(dat_home |> select(where(is.numeric)),
                 function(x) sum(is.infinite(x) | is.nan(x)))
if (any(nonfin > 0)) { print(nonfin[nonfin > 0])
  stop("Non-finite values in the home panel -- fix before modeling.") }
message("63 complete.")
