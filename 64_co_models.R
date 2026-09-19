# ==============================================================================
# 64_co_models.R      [COLORADO ZONING case study, step 4]
# Descriptives + models (fixest throughout, 55's conventions).
#
#  BLOCK ZA  Zoning -> residential segregation (2023 cross-section):
#              feols(res_seg ~ zoning + covs | county_fips, cl ~jurisd_main)
#            Each zoning measure separately (multicollinearity discipline,
#            as paper 1), then a joint spec.
#  BLOCK ZB  HEADLINE -- zoning as moderator of the res -> workhood link:
#              feols(wexp ~ res_seg * zoning + covs | county_fips,
#                    cl ~jurisd_main)
#            Does living in a segregated tract predict commuting into a more
#            segregated workhood MORE where zoning is exclusionary, and LESS
#            where residential land allows ADUs / uses are mixed?
#  BLOCK ZC  Converse (work tracts): what kind of ZONED land hosts workhoods
#            drawing from segregated neighborhoods?
#              feols(wres ~ work_seg + zoning-of-work-tract + covs | county)
#  BLOCK ZD  Trajectories: d(res_seg) and d(wexp) 2011->2023 ~ 2023 zoning +
#            baseline covs; plus persistence spec res23 ~ res11 x zoning.
#            (Zoning changes slowly; still a snapshot -- descriptive only.)
#  BLOCK ZM  Mechanism: commute distance / destination portfolio / jobs-
#            housing balance ~ zoning (the job-search-radius story of P3).
#  Sensitivities: aspatial D; hisp_nonhisp + college_lesshs; coverage >= 0.50;
#            + Boulder & Greeley pooled (CBSA FE); min commuters 50;
#            pct_reslow_of_res as the alternative exclusionary denominator;
#            NO-DENSITY variants of ZA + ZB (worker density is plausibly a
#            MEDIATOR of zoning -- zoning caps density -- so conditioning on
#            it may over-control; compare with/without);
#            exclusionary-measure comparison: pct_res_low vs pct_reslow_of_res
#            side-by-side across every variant + a joint horse race
#            (r = 0.91 between them -- read the joint spec with that in mind;
#            output: co_exclusionary_measure_comparison.csv);
#            transit moderator head-to-head IF the SLD file exists (53_sld.R).
#
# Epistemic framing (same as papers 1-3): descriptive/moderation language, NO
# causal claims. Zoning is endogenous to who lives where -- exclusionary
# zoning both produces and is produced by segregation. Say so in the text.
#
# Output: output/models/co_model_coefficients.csv
#         output/models/co_descriptives_2023.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
datw <- readRDS(file.path(DIR_CO_OUT, "co_analysis_work_panel.rds"))
chg  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_change_2011_2023.rds"))

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
tidy_fx <- function(fit, model_id, keep = "z_|:") {
  if (is.null(fit)) return(NULL)
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(grepl(keep, term)) |>
    mutate(conf.low  = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           model_id = model_id, n_obs = fit$nobs)
}
safe_feols <- function(fml, data, cluster, model_id, ...) {
  fit <- tryCatch(feols(as.formula(fml), data = data, cluster = cluster, ...),
                  error = function(e) { message("  FAILED ", model_id, ": ",
                                                conditionMessage(e)); NULL })
  tidy_fx(fit, model_id)
}
res <- list()
add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x

## =============================================================================
## Model frames
## =============================================================================
ZONING_VARS <- c("pct_res_low", "pct_reslow_of_res", "pct_adu_res",
                 "pct_grouphome_res", "zoning_entropy", "pct_job_zone",
                 "pct_mixed_res", "pct_industrial", "pct_commercial")
MODEL_VARS <- c("wexp_whiteblack_wac_half", "wexp_whiteblack_wac_aspatial",
                "d_whiteblack_rac_half", "d_whiteblack_rac_aspatial",
                "d_whiteblack_wac_half",
                "d_hisp_nonhisp_rac_half", "wexp_hisp_nonhisp_wac_half",
                "d_college_lesshs_rac_half", "wexp_college_lesshs_wac_half",
                "mean_dist_km", "eff_n_dest", "pct_same_tract",
                "jobs_housing_ratio", "log_sld_D5BR", "sld_NatWalkInd",
                ZONING_VARS,
                "pct_black_rac", "pct_lowincome_rac",
                "log_worker_density_rac", "income_percapita_k",
                "income_percapita_k_sq")

# cross-section frame: anchor year, in-scope, Denver MSA (main sample)
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  mutate(across(any_of(MODEL_VARS), zscore, .names = "z_{.col}"))
# pooled 3-CBSA frame (z-scores within the pooled sample; CBSA FE in models)
xs_all <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope,
         n_commuters >= P3_MIN_COMMUTERS) |>
  mutate(across(any_of(MODEL_VARS), zscore, .names = "z_{.col}"))
# relaxed-coverage sensitivity frame
xs_cov50 <- dat |>
  filter(year == CO_ANCHOR_YEAR, cover_zoned >= 0.50, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  mutate(across(any_of(MODEL_VARS), zscore, .names = "z_{.col}"))

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
# no-density covariate set: worker density is downstream of zoning (zoning
# caps density), so the density-conditional estimates may over-control --
# the _nodens variants bound the zoning association from the other side
COVS_NODENS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                       "z_income_percapita_k", "z_income_percapita_k_sq"),
                     collapse = " + ")

message(sprintf("Model frames: Denver %d tracts | pooled %d | cov50 %d",
                nrow(xs), nrow(xs_all), nrow(xs_cov50)))

## =============================================================================
## Descriptives: 2023 means by exclusionary-zoning tercile
## =============================================================================
desc <- xs |>
  mutate(reslow_tercile = ntile(pct_reslow_of_res, 3)) |>
  group_by(reslow_tercile) |>
  summarise(n = n(),
            across(c(pct_res_low, pct_adu_res, zoning_entropy, pct_job_zone,
                     d_whiteblack_rac_half, wexp_whiteblack_wac_half,
                     mean_dist_km, eff_n_dest, jobs_housing_ratio,
                     pct_black_rac, income_percapita_k),
                   ~ round(mean(.x, na.rm = TRUE), 3)),
            .groups = "drop")
write.csv(desc, file.path(DIR_CO_MOD, "co_descriptives_2023.csv"),
          row.names = FALSE)

# zoning-measure correlation matrix (multicollinearity eyeball for ZA-joint)
write_codiag(
  round(cor(xs |> select(any_of(ZONING_VARS)), use = "pairwise"), 2) |>
    as.data.frame() |> tibble::rownames_to_column("measure"),
  "64_zoning_correlations")

## =============================================================================
## BLOCK ZA -- zoning -> residential segregation (2023 cross-section)
## =============================================================================
for (ms in c("whiteblack", "hisp_nonhisp", "college_lesshs")) {
  y <- sprintf("z_d_%s_rac_half", ms)
  if (!y %in% names(xs)) next
  for (zv in ZONING_VARS)
    add(safe_feols(sprintf("%s ~ z_%s + %s | county_fips", y, zv, COVS),
                   xs, ~jurisd_main, sprintf("ZA_%s_%s", ms, zv)))
  # joint spec: the three least collinear dimensions of the zoning regime
  add(safe_feols(sprintf(
    "%s ~ z_pct_res_low + z_pct_adu_res + z_pct_job_zone + %s | county_fips",
    y, COVS), xs, ~jurisd_main, sprintf("ZA_joint_%s", ms)))
}
# own-tract WORKPLACE segregation too (is WAC seg zoned into industrial land?)
add(safe_feols(sprintf(
  "z_d_whiteblack_wac_half ~ z_pct_industrial + z_pct_job_zone + %s | county_fips",
  COVS), xs, ~jurisd_main, "ZA_wac_industrial"))

# no-density sensitivity (whiteblack): each zoning measure without the
# density covariate -- the raw res-low seg gradient (tercile 3 D ~ 2x
# tercile 1) is absorbed almost entirely by density in the main spec
for (zv in ZONING_VARS)
  add(safe_feols(sprintf(
    "z_d_whiteblack_rac_half ~ z_%s + %s | county_fips", zv, COVS_NODENS),
    xs, ~jurisd_main, sprintf("ZA_nodens_whiteblack_%s", zv)))

## =============================================================================
## BLOCK ZB -- HEADLINE: zoning moderates the res -> workhood link
## =============================================================================
ZB_MODS <- c("pct_res_low", "pct_reslow_of_res", "pct_adu_res",
             "zoning_entropy", "pct_job_zone")
for (ms in c("whiteblack", "hisp_nonhisp", "college_lesshs")) {
  y <- sprintf("z_wexp_%s_wac_half", ms)
  x <- sprintf("z_d_%s_rac_half", ms)
  if (!all(c(y, x) %in% names(xs))) next
  add(safe_feols(sprintf("%s ~ %s + %s | county_fips", y, x, COVS),
                 xs, ~jurisd_main, sprintf("ZB_base_%s", ms)))
  for (m in ZB_MODS)
    add(safe_feols(
      sprintf("%s ~ %s * z_%s + %s | county_fips", y, x, m, COVS),
      xs, ~jurisd_main, sprintf("ZB_mod_%s_%s", ms, m)))
}
# sensitivities (whiteblack, headline moderator set)
for (m in ZB_MODS) {
  add(safe_feols(sprintf(
    "z_wexp_whiteblack_wac_aspatial ~ z_d_whiteblack_rac_aspatial * z_%s + %s | county_fips",
    m, COVS), xs, ~jurisd_main, sprintf("ZB_aspatial_%s", m)))
  add(safe_feols(sprintf(   # county FE nest inside CBSA here, so county FE
    "z_wexp_whiteblack_wac_half ~ z_d_whiteblack_rac_half * z_%s + %s | county_fips",
    m, COVS), xs_all, ~jurisd_main, sprintf("ZB_pooled3_%s", m)))
  add(safe_feols(sprintf(
    "z_wexp_whiteblack_wac_half ~ z_d_whiteblack_rac_half * z_%s + %s | county_fips",
    m, COVS), xs_cov50, ~jurisd_main, sprintf("ZB_cov50_%s", m)))
  add(safe_feols(sprintf(
    "z_wexp_whiteblack_wac_half ~ z_d_whiteblack_rac_half * z_%s + %s | county_fips",
    m, COVS), xs |> filter(n_commuters >= 50), ~jurisd_main,
    sprintf("ZB_mincomm50_%s", m)))
  add(safe_feols(sprintf(
    "z_wexp_whiteblack_wac_half ~ z_d_whiteblack_rac_half * z_%s + %s | county_fips",
    m, COVS_NODENS), xs, ~jurisd_main, sprintf("ZB_nodens_%s", m)))
}
# exclusionary-measure horse race: both denominators' interactions in one
# model. r(pct_res_low, pct_reslow_of_res) = 0.91 -- coefficients split a
# shared signal; the separate models + comparison table are the honest read,
# this spec only asks which one carries the residual variation
add(safe_feols(paste(
  "z_wexp_whiteblack_wac_half ~ z_d_whiteblack_rac_half * z_pct_res_low +",
  "z_d_whiteblack_rac_half * z_pct_reslow_of_res +", COVS, "| county_fips"),
  xs, ~jurisd_main, "ZB_horserace_exclusionary"))
# zoning vs transit head-to-head (only if the SLD join happened in 63)
if ("z_log_sld_D5BR" %in% names(xs)) {
  add(safe_feols(paste(
    "z_wexp_whiteblack_wac_half ~ z_d_whiteblack_rac_half * z_pct_res_low +",
    "z_d_whiteblack_rac_half * z_log_sld_D5BR +", COVS, "| county_fips"),
    xs, ~jurisd_main, "ZB_horse_race_reslow_vs_transit"))
}

## =============================================================================
## BLOCK ZC -- converse: zoning of the workhoods that draw from segregated
## neighborhoods
## =============================================================================
xw <- datw |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_workers_od >= P3_MIN_COMMUTERS) |>
  mutate(across(any_of(c("wres_whiteblack_rac_half", "d_whiteblack_wac_half",
                         ZONING_VARS, "pct_manuf_wac", "pct_black_wac",
                         "pct_lowincome_wac", "log_worker_density_wac")),
                zscore, .names = "z_{.col}"))
COVS_W <- paste(c("z_pct_black_wac", "z_pct_lowincome_wac", "z_pct_manuf_wac",
                  "z_log_worker_density_wac"), collapse = " + ")
add(safe_feols(paste(
  "z_wres_whiteblack_rac_half ~ z_d_whiteblack_wac_half +", COVS_W,
  "| county_fips"), xw, ~jurisd_main, "ZC_base"))
for (zv in c("pct_industrial", "pct_commercial", "pct_job_zone",
             "zoning_entropy"))
  add(safe_feols(sprintf(
    "z_wres_whiteblack_rac_half ~ z_d_whiteblack_wac_half + z_%s + %s | county_fips",
    zv, COVS_W), xw, ~jurisd_main, sprintf("ZC_%s", zv)))

## =============================================================================
## BLOCK ZD -- trajectories 2011 -> 2023 (2023 zoning, baseline covariates)
## =============================================================================
chg_z <- chg |>
  filter(in_scope, denver_msa) |>
  mutate(across(any_of(c("d_res_seg", "d_wexp", "res_seg_2011", "res_seg_2023",
                         "wexp_2011", ZONING_VARS,
                         "pct_black_rac_2011", "pct_lowincome_rac_2011",
                         "income_percapita_k_2011",
                         "log_worker_density_rac_2011")),
                zscore, .names = "z_{.col}"))
COVS_0 <- paste(c("z_pct_black_rac_2011", "z_pct_lowincome_rac_2011",
                  "z_income_percapita_k_2011",
                  "z_log_worker_density_rac_2011"), collapse = " + ")
for (out in c("z_d_res_seg", "z_d_wexp"))
  for (zv in c("pct_reslow_of_res", "pct_res_low", "pct_adu_res",
               "zoning_entropy", "pct_job_zone"))
    add(safe_feols(sprintf("%s ~ z_%s + %s | county_fips", out, zv, COVS_0),
                   chg_z, ~jurisd_main,
                   sprintf("ZD_%s_%s", sub("^z_d_", "", out), zv)))
# persistence: is 2011 segregation MORE predictive of 2023 segregation in
# exclusionary-zoned tracts? (lock-in)
add(safe_feols(paste(
  "z_res_seg_2023 ~ z_res_seg_2011 * z_pct_res_low +", COVS_0,
  "| county_fips"), chg_z, ~jurisd_main, "ZD_persistence_reslow"))
# the same with the paper's primary exclusionary measure (share of RESIDENTIAL
# land zoned low-density); the all-zoned-land version above is kept for comparison
add(safe_feols(paste(
  "z_res_seg_2023 ~ z_res_seg_2011 * z_pct_reslow_of_res +", COVS_0,
  "| county_fips"), chg_z, ~jurisd_main, "ZD_persistence_reslow_of_res"))
# functional-form checks on the primary-measure persistence model: residential
# segregation is strongly right-skewed, so refit with both years winsorized at
# their 99th percentile and with both years entered as ranks
wins99 <- function(x) pmin(x, stats::quantile(x, .99, na.rm = TRUE))
chg_w <- chg_z |> mutate(z_res_seg_2011 = zscore(wins99(res_seg_2011)),
                         z_res_seg_2023 = zscore(wins99(res_seg_2023)))
chg_r <- chg_z |> mutate(z_res_seg_2011 = zscore(rank(res_seg_2011, na.last = "keep")),
                         z_res_seg_2023 = zscore(rank(res_seg_2023, na.last = "keep")))
add(safe_feols(paste(
  "z_res_seg_2023 ~ z_res_seg_2011 * z_pct_reslow_of_res +", COVS_0,
  "| county_fips"), chg_w, ~jurisd_main, "ZD_persistence_reslow_of_res_wins99"))
add(safe_feols(paste(
  "z_res_seg_2023 ~ z_res_seg_2011 * z_pct_reslow_of_res +", COVS_0,
  "| county_fips"), chg_r, ~jurisd_main, "ZD_persistence_reslow_of_res_rank"))
add(safe_feols(paste(
  "z_res_seg_2023 ~ z_res_seg_2011 * z_pct_adu_res +", COVS_0,
  "| county_fips"), chg_z, ~jurisd_main, "ZD_persistence_adu"))

## =============================================================================
## BLOCK ZM -- mechanism: commute geography and jobs-housing balance
## =============================================================================
for (out in c("z_mean_dist_km", "z_eff_n_dest", "z_pct_same_tract",
              "z_jobs_housing_ratio")) {
  add(safe_feols(sprintf(
    "%s ~ z_pct_res_low + z_pct_job_zone + %s | county_fips", out, COVS),
    xs, ~jurisd_main, sprintf("ZM_zoning_%s", sub("^z_", "", out))))
  add(safe_feols(sprintf(
    "%s ~ z_d_whiteblack_rac_half * z_pct_res_low + %s | county_fips",
    out, COVS), xs, ~jurisd_main,
    sprintf("ZM_interact_%s", sub("^z_", "", out))))
}

## ---- write tables ------------------------------------------------------------
all_res <- bind_rows(res)
write.csv(all_res, file.path(DIR_CO_MOD, "co_model_coefficients.csv"),
          row.names = FALSE)
message("Fitted ", n_distinct(all_res$model_id), " models.")

# exclusionary-measure comparison: every pct_res_low / pct_reslow_of_res
# coefficient (levels + interactions) side by side across all variants
cmp <- all_res |>
  filter(grepl("z_pct_res_low|z_pct_reslow_of_res", term)) |>
  mutate(measure = ifelse(grepl("reslow_of_res", term),
                          "reslow_of_res", "res_low"),
         is_interaction = grepl(":", term)) |>
  select(model_id, measure, is_interaction, term,
         estimate, std.error, p.value, n_obs) |>
  arrange(desc(is_interaction), measure, model_id)
write.csv(cmp, file.path(DIR_CO_MOD, "co_exclusionary_measure_comparison.csv"),
          row.names = FALSE)
message("Exclusionary-measure comparison: ", nrow(cmp), " rows.")

## ---- DIAGNOSTICS -------------------------------------------------------------
write_codiag(
  all_res |>
    mutate(block = sub("_.*", "", model_id),
           sig = p.value < 0.05,
           direction = ifelse(estimate > 0, "pos", "neg")) |>
    count(block, direction, sig),
  "64_model_sign_summary")
# Eyeball before interpreting: n_obs (Denver-only 2023 in-scope frame is a
# few hundred tracts -- power is limited; that is WHY jurisd_main clustering
# with 27 clusters, not county with ~7, is the default).
message("64 complete. Tables in ", DIR_CO_MOD)
