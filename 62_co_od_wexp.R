# ==============================================================================
# 62_co_od_wexp.R      [COLORADO ZONING case study, step 2]
# Colorado-only OD-weighted workplace-segregation exposure, 2011-2023.
# Same construct as 52 (Eq. 1), computed standalone so this case study does
# NOT wait on the national 51/52 run:
#     wexp_jt = sum_i (flow_jit / sum_i flow_jit) * D_wac_it
# plus the work-tract converse (wres) and OD commute measures.
#
# OD source, in order of preference:
#   1. the national 51 caches raw/od_tract_{part}_{yr}_co.rds (reused as-is)
#   2. else downloaded here via lehdr (co main+aux, JT01, tract-agg;
#      ~1 min/year) and cached in Colorado/clean/
# aux for Colorado = workers living OUT of state at Colorado workplaces; the
# same-MSA restriction removes them anyway (their home tract has no CO CBSA),
# but keeping the parts identical to 51/52 makes the numbers comparable.
#
# Restriction (as in 52): home AND work tract in the SAME MSA, and the MSA in
# CO_CBSA_KEEP (19740 Denver + 14500 Boulder + 24540 Greeley -- the
# jurisdictions the zoning file spans). Drop shares logged per year.
#
# Needs: p2_tract_segregation_panel.rds (32), tract_centroids_km.rds (21/02).
# Output: clean/co_wexp_{year}.rds; clean/co_wexp_home_panel.rds;
#         clean/co_wres_work_panel.rds
# ==============================================================================

source("60_co_setup.R")
library(lehdr)

seg_panel <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds"))
tract_centroids <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds"))

cent <- tract_centroids |>
  filter(CBSA_Code %in% CO_CBSA_KEEP) |>
  select(GEOID, CBSA_Code, X_km, Y_km)

seg_cols_wac <- as.vector(outer(P3_SEG_MEASURES, P3_SEG_POWERS,
                                \(m, p) sprintf("d_%s_wac_%s", m, p)))
seg_cols_rac <- as.vector(outer(P3_SEG_MEASURES, P3_SEG_POWERS,
                                \(m, p) sprintf("d_%s_rac_%s", m, p)))
stopifnot(all(c(seg_cols_wac, seg_cols_rac) %in% names(seg_panel)))

wtd_median <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  o <- order(x[ok]); x <- x[ok][o]; w <- w[ok][o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}
wmean_na <- function(d, w) {          # renormalizes over non-missing D (52)
  ok <- !is.na(d) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(d[ok] * w[ok]) / sum(w[ok])
}

## ---- 1. Colorado OD per year (reuse 51 caches when present) ------------------
grab_co_od <- function(yr) {
  per_part <- map(P3_OD_PARTS, function(part) {
    nat_cache <- p3_od_cache(part, yr, "co")     # written by the national 51
    co_cache  <- co_od_cache(part, yr)
    if (file.exists(nat_cache)) return(readRDS(nat_cache))
    if (file.exists(co_cache))  return(readRDS(co_cache))
    message(sprintf("Downloading OD %s %s co", part, yr))
    df <- tryCatch(
      grab_lodes(state = "co", year = yr, version = "LODES8",
                 lodes_type = "od", job_type = "JT01",
                 state_part = part, agg_geo = "tract"),
      error = function(e)
        # never build exposures from an incomplete (e.g. main-only) OD table
        stop(sprintf("OD download failed (%s %s co): %s", part, yr,
                     conditionMessage(e))))
    if (is.null(df)) return(NULL)
    df <- df |>
      transmute(w_tract = as.character(w_tract),
                h_tract = as.character(h_tract),
                across(any_of(P3_OD_COLS), ~ as.numeric(.x)))
    saveRDS(df, co_cache)
    df
  }) |> compact() |> bind_rows()
  if (nrow(per_part) == 0) return(NULL)
  per_part |>                                    # belt-and-braces, as in 51
    group_by(w_tract, h_tract) |>
    summarise(across(any_of(P3_OD_COLS), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
}

## ---- 2. weighted exposure per year (52's aggregation, CO scope) --------------
drop_log <- list()

for (yr in CO_YEARS) {
  ck <- co_wexp_year(yr)
  if (file.exists(ck)) { message(basename(ck), " exists"); next }
  od <- grab_co_od(yr)
  if (is.null(od)) { message("No OD for ", yr, " -- skipped."); next }

  od <- od |>
    inner_join(cent, by = c("h_tract" = "GEOID")) |>
    rename(cbsa_h = CBSA_Code, xh = X_km, yh = Y_km) |>
    inner_join(cent, by = c("w_tract" = "GEOID")) |>
    rename(cbsa_w = CBSA_Code, xw = X_km, yw = Y_km)

  jobs_metro  <- sum(od$S000)
  od          <- od |> filter(cbsa_h == cbsa_w)
  jobs_within <- sum(od$S000)
  drop_log[[as.character(yr)]] <- tibble(
    year = yr, jobs_both_metro = jobs_metro, jobs_same_msa = jobs_within,
    pct_dropped_cross_msa = 100 * (1 - jobs_within / jobs_metro))

  od <- od |> mutate(dist_km = sqrt((xh - xw)^2 + (yh - yw)^2))

  seg_yr <- seg_panel |> filter(year == yr)
  od <- od |>
    left_join(seg_yr |> select(tract_id, all_of(seg_cols_wac)),
              by = c("w_tract" = "tract_id")) |>
    left_join(seg_yr |> select(tract_id, all_of(seg_cols_rac)),
              by = c("h_tract" = "tract_id"))

  home <- od |>
    group_by(tract_id = h_tract) |>
    summarise(
      cbsa          = first(cbsa_h),
      n_commuters   = sum(S000),
      across(all_of(seg_cols_wac),
             ~ wmean_na(.x, S000), .names = "wexp_{.col}"),
      wexp_se01_d_whiteblack_wac_half = wmean_na(d_whiteblack_wac_half, SE01),
      wexp_se03_d_whiteblack_wac_half = wmean_na(d_whiteblack_wac_half, SE03),
      n_commuters_se01 = sum(SE01), n_commuters_se03 = sum(SE03),
      mean_dist_km    = weighted.mean(dist_km, w = S000),
      median_dist_km  = wtd_median(dist_km, S000),
      pct_same_tract  = 100 * sum(S000[h_tract == w_tract]) / sum(S000),
      pct_lt5km       = 100 * sum(S000[dist_km < 5])  / sum(S000),
      pct_gt24km      = 100 * sum(S000[dist_km > 24]) / sum(S000),
      eff_n_dest      = 1 / sum((S000 / sum(S000))^2),
      .groups = "drop") |>
    mutate(year = yr)
  names(home) <- sub("^wexp_d_", "wexp_", names(home))
  names(home) <- sub("^(wexp_se0[13])_d_", "\\1_", names(home))

  work <- od |>
    group_by(tract_id = w_tract) |>
    summarise(
      cbsa         = first(cbsa_w),
      n_workers_od = sum(S000),
      across(all_of(seg_cols_rac),
             ~ wmean_na(.x, S000), .names = "wres_{.col}"),
      mean_dist_km_work = weighted.mean(dist_km, w = S000),
      eff_n_origin      = 1 / sum((S000 / sum(S000))^2),
      .groups = "drop") |>
    mutate(year = yr)
  names(work) <- sub("^wres_d_", "wres_", names(work))

  saveRDS(list(home = home, work = work), ck)
  message(sprintf("%d: %s home tracts, %s work tracts, %.1f%% cross-MSA dropped",
                  yr, nrow(home), nrow(work),
                  drop_log[[as.character(yr)]]$pct_dropped_cross_msa))
  rm(od, home, work); gc()
}

## ---- 3. assemble panels ------------------------------------------------------
pieces <- map(CO_YEARS, function(yr) {
  f <- co_wexp_year(yr)
  if (file.exists(f)) readRDS(f) else NULL
}) |> compact()

home_panel <- map(pieces, "home") |> bind_rows()
work_panel <- map(pieces, "work") |> bind_rows()
saveRDS(home_panel, file.path(DIR_CO_CLEAN, "co_wexp_home_panel.rds"))
saveRDS(work_panel, file.path(DIR_CO_CLEAN, "co_wres_work_panel.rds"))
message("Saved home panel (", nrow(home_panel), " tract-years) and work panel (",
        nrow(work_panel), " tract-years)")

## ---- DIAGNOSTICS -------------------------------------------------------------
if (length(drop_log) > 0)
  write_codiag(bind_rows(drop_log), "62_cross_msa_jobs_dropped_by_year")
# Eyeball: pct_dropped_cross_msa should be modest and SMOOTH; Denver-Boulder
# cross-commuting is real, so expect a larger share than the national 5-10%
# (both flows are dropped by the same-MSA rule, consistent with 52).

write_codiag(
  home_panel |> group_by(cbsa, year) |>
    summarise(n_tracts = n(),
              med_commuters = median(n_commuters),
              pct_na_wexp = 100 * mean(is.na(wexp_whiteblack_wac_half)),
              mean_wexp   = mean(wexp_whiteblack_wac_half, na.rm = TRUE),
              .groups = "drop"),
  "62_home_panel_by_cbsa_year")

seg_res <- seg_panel |> select(tract_id, year, d_whiteblack_rac_half)
qc <- home_panel |>
  inner_join(seg_res, by = c("tract_id", "year")) |>
  group_by(year) |>
  summarise(r_resseg_wexp = cor(d_whiteblack_rac_half,
                                wexp_whiteblack_wac_half,
                                use = "complete.obs"))
write_codiag(qc, "62_res_vs_wexp_correlation_by_year")
# Eyeball: r should sit near the national tract-level value once 52 exists;
# the MSA paper found r ~ 0.5-0.6 at MSA level.
message("62 complete.")
