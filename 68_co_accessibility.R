# ==============================================================================
# 68_co_accessibility.R      [PAPER 4, step 8]
# THE CONFOUND TEST. Is the zoning moderation of the residential -> workplace
# segregation link (64's Block ZB) actually ZONING, or is it metropolitan
# POSITION -- peripheral tracts happening to have both exclusionary zoning and
# tighter coupling?
#
# Builds accessibility/centrality measures for every regional tract from data
# already on disk (no new downloads):
#   dist_cbd_km    straight-line km to Denver CBD (Union Station)
#   dist_empctr_km km to the nearest major employment center (tract in the top
#                  2% of 2023 WAC jobs)
#   jobs_10/20km   cumulative jobs within 10 / 20 km (opportunity counts)
#   jobs_grav      gravity-weighted jobs, sum_i C000_i * exp(-0.10 * d_ik)
#                  NOTE: this 0.10/km decay is a JOB-ACCESS parameter (a
#                  ~7-km half-distance commute shed). It is deliberately NOT
#                  the segregation decay (beta = 0.5, 10-km cutoff, the
#                  corrected 2018 convention) -- the two measure different
#                  things and must not be conflated.
#
# THE THREE-RUNG LADDER (this is the paper's central argument):
#   A_total                  wexp ~ res_seg * zoning + covs        (naive; the
#                            specification the zoning-segregation literature
#                            typically reports)
#   B_plus_regional_position + dist_cbd, dist_empctr AND their interactions
#                            with res_seg. This is the PREFERRED specification
#                            for decomposing the association. NOTE ON LANGUAGE:
#                            position is NOT thereby a proven confounder --
#                            metropolitan development, annexation, infrastructure,
#                            employment decentralization and zoning co-evolved.
#                            The ladder decomposes a cross-sectional association
#                            into components tracked by regional geography vs
#                            jurisdictional regulation; it does not identify a
#                            causal zoning effect.
#   C_plus_local_jobsurface  + jobs_grav and its interaction. The local job
#                            surface is plausibly partly downstream of zoning
#                            (that is what pct_job_zone allows), so this rung
#                            is an OVER-CONTROL SENSITIVITY that may remove
#                            part of the pathway through which zoning operates
#                            -- never the preferred estimate, and not a formal
#                            bound (conditioning on a mediator does not bound
#                            a total association).
#
# Verified result (Aug 2026): adjustment for metropolitan position reduces the
# estimated zoning moderation by roughly two-thirds; ~30% of the baseline
# estimate remains at rung B (the position-adjusted zoning association),
# individually significant only for the exclusionary measure; little remains
# at rung C.
# NOTE: jobs_grav correlates about -0.93 with dist_cbd_km in Denver -- in a
# broadly monocentric region "job access" and "centrality" are the same thing.
# Report that correlation; do not present them as independent controls.
#
# Needs: 63's home panel, co_wac_tract_2023.rds (66), tract_centroids_km.rds.
# Output: output/models/p4_accessibility_ladder.csv
#         clean/co_accessibility_2023.rds
# ==============================================================================

source("60_co_setup.R")
library(fixest)

## ---- 1. accessibility surface ------------------------------------------------
acc_file <- file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds")
if (file.exists(acc_file)) {
  message(basename(acc_file), " exists")
  acc <- readRDS(acc_file)
} else {
  cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
    filter(CBSA_Code %in% CO_CBSA_KEEP) |>
    transmute(tract_id = as.character(GEOID), CBSA_Code, X_km, Y_km)

  wac <- readRDS(file.path(DIR_CO_CLEAN,
                           sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
    transmute(tract_id = as.character(tract_id), jobs = C000)

  opp <- cent |> left_join(wac, by = "tract_id") |>
    mutate(jobs = ifelse(is.na(jobs), 0, jobs))
  message(nrow(opp), " regional tracts; ", format(sum(opp$jobs), big.mark = ","),
          " jobs")

  # Denver CBD (Union Station) -> the same CRS the centroids use (5070, km)
  cbd <- st_sfc(st_point(c(-105.0000, 39.7531)), crs = 4326) |>
    st_transform(CRS_METERS) |> st_coordinates()
  cbd_x <- cbd[1] / 1000; cbd_y <- cbd[2] / 1000

  XY <- as.matrix(opp[, c("X_km", "Y_km")])
  D  <- as.matrix(dist(XY))          # km, tract-centroid to tract-centroid
  J  <- opp$jobs

  opp <- opp |>
    mutate(
      dist_cbd_km = sqrt((X_km - cbd_x)^2 + (Y_km - cbd_y)^2),
      jobs_10km   = as.numeric((D <= 10) %*% J),
      jobs_20km   = as.numeric((D <= 20) %*% J),
      jobs_grav   = as.numeric(exp(-0.10 * D) %*% J))
  centers <- J >= quantile(J, 0.98)
  message(sum(centers), " employment centers (top 2% of tracts by jobs)")
  opp$dist_empctr_km <- apply(D[, centers, drop = FALSE], 1, min)

  acc <- opp |>
    mutate(across(c(jobs_10km, jobs_20km, jobs_grav),
                  ~ log(pmax(.x, 1)), .names = "log_{.col}")) |>
    select(tract_id, dist_cbd_km, dist_empctr_km,
           jobs_10km, jobs_20km, jobs_grav,
           log_jobs_10km, log_jobs_20km, log_jobs_grav)
  saveRDS(acc, acc_file)
  message("Saved accessibility measures for ", nrow(acc), " tracts")
}

## ---- 2. model frame (identical to 64's Denver cross-section) -----------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }

MODEL_VARS <- c("wexp_whiteblack_wac_half", "d_whiteblack_rac_half",
                "pct_res_low", "pct_reslow_of_res", "pct_adu_res",
                "zoning_entropy", "pct_job_zone",
                "pct_black_rac", "pct_lowincome_rac", "log_worker_density_rac",
                "income_percapita_k", "income_percapita_k_sq",
                "dist_cbd_km", "dist_empctr_km", "log_jobs_grav",
                "log_jobs_10km", "jobs_housing_ratio", "mean_dist_km",
                "eff_n_dest")
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  mutate(across(any_of(MODEL_VARS), zscore, .names = "z_{.col}"))
message("Denver ", CO_ANCHOR_YEAR, " in-scope frame: ", nrow(xs), " tracts")

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
X   <- "z_d_whiteblack_rac_half"
POS <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                     "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)
LOC <- sprintf("z_log_jobs_grav + %s:z_log_jobs_grav", X)

## ---- 3. the ladder -----------------------------------------------------------
tidy_fx <- function(fit, model_id, keep) {
  if (is.null(fit)) return(NULL)
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |>
    mutate(conf.low = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           model_id = model_id, n_obs = fit$nobs)
}
safe <- function(fml, model_id, keep)
  tidy_fx(tryCatch(feols(as.formula(fml), data = xs, cluster = ~jurisd_main),
                   error = function(e) { message("  FAILED ", model_id, ": ",
                                                 conditionMessage(e)); NULL }),
          model_id, keep)

ZMODS <- c("pct_reslow_of_res", "pct_res_low", "pct_adu_res",
           "zoning_entropy", "pct_job_zone")
res <- list()
for (zv in ZMODS) {
  base <- sprintf("z_wexp_whiteblack_wac_half ~ %s * z_%s + %s | county_fips",
                  X, zv, COVS)
  keep <- c(sprintf("%s:z_%s", X, zv), sprintf("z_%s", zv), X,
            sprintf("%s:z_log_jobs_grav", X), sprintf("%s:z_dist_cbd_km", X))
  specs <- list(
    A_total = base,
    B_plus_regional_position = sub("\\| county_fips",
                                   paste("+", POS, "| county_fips"), base),
    C_plus_local_jobsurface  = sub("\\| county_fips",
                                   paste("+", POS, "+", LOC, "| county_fips"),
                                   base))
  for (s in names(specs))
    res[[length(res) + 1]] <- safe(specs[[s]], sprintf("%s__%s", zv, s), keep)
}
# position-control sensitivity: CBD only / employment-centre only / both
for (zv in c("pct_reslow_of_res", "pct_adu_res", "zoning_entropy",
             "pct_job_zone")) {
  base <- sprintf("z_wexp_whiteblack_wac_half ~ %s * z_%s + %s", X, zv, COVS)
  vars <- list(
    B1_cbd_only    = sprintf("z_dist_cbd_km + %s:z_dist_cbd_km", X),
    B2_empctr_only = sprintf("z_dist_empctr_km + %s:z_dist_empctr_km", X),
    B3_both        = POS)
  for (s in names(vars))
    res[[length(res) + 1]] <- safe(
      sprintf("%s + %s | county_fips", base, vars[[s]]),
      sprintf("%s__%s", zv, s), c(sprintf("%s:z_%s", X, zv), X))
}
all_res <- bind_rows(res)
write.csv(all_res, file.path(DIR_CO_MOD, "p4_accessibility_ladder.csv"),
          row.names = FALSE)
message("Fitted ", n_distinct(all_res$model_id), " ladder models.")

## ---- 4. share of the baseline estimate remaining at each rung ----------------
surv <- all_res |>
  filter(grepl("__(A_total|B_plus_regional_position|C_plus_local_jobsurface)$",
               model_id)) |>
  mutate(zoning = sub("__.*", "", model_id),
         spec   = sub(".*__", "", model_id)) |>
  # keep ONLY each model's own zoning interaction: rung C also carries the
  # job-surface and distance interactions, which would duplicate (zoning, spec)
  # and turn the pivot into list-columns
  filter(term == paste0(X, ":z_", zoning)) |>
  select(zoning, spec, estimate, p.value) |>
  pivot_wider(names_from = spec, values_from = c(estimate, p.value)) |>
  mutate(pct_of_baseline_remaining_position =
           100 * estimate_B_plus_regional_position / estimate_A_total,
         pct_of_baseline_remaining_all =
           100 * estimate_C_plus_local_jobsurface / estimate_A_total)
write.csv(surv, file.path(DIR_CO_MOD, "p4_moderation_survival.csv"),
          row.names = FALSE)
print(as.data.frame(surv))

## ---- DIAGNOSTICS -------------------------------------------------------------
# The collinearity that drives the whole result -- REPORT THIS IN THE PAPER.
cc <- c("pct_reslow_of_res", "pct_job_zone", "zoning_entropy", "pct_adu_res",
        "log_jobs_grav", "dist_cbd_km", "dist_empctr_km",
        "log_worker_density_rac")
write_codiag(
  round(cor(xs[cc], use = "pairwise"), 3) |> as.data.frame() |>
    tibble::rownames_to_column("measure"),
  "68_zoning_vs_accessibility_correlations")
# Verified: cor(log_jobs_grav, dist_cbd_km) ~ -0.93. In a monocentric metro
# these are one construct; say so rather than implying two controls.

write_codiag(
  xs |> summarise(across(all_of(c("dist_cbd_km", "dist_empctr_km",
                                  "log_jobs_grav", "jobs_10km")),
                         list(min = ~min(.x, na.rm = TRUE),
                              med = ~median(.x, na.rm = TRUE),
                              max = ~max(.x, na.rm = TRUE)))) |>
    pivot_longer(everything()),
  "68_accessibility_distributions")
message("68 complete.")
