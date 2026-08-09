# ==============================================================================
# 83_co_matched_access.R      [PAPER 4, step 23 -- SI extension: matched access]
# Earnings-matched accessibility: how many LOW-WAGE jobs are accessible to
# LOW-WAGE resident workers from each neighborhood? Aggregate job access can
# look reasonable while access to jobs that match the resident workforce is
# poor; this measures the match directly.
#
# Measure: Shen's (1998) competition-adjusted accessibility, per earnings band:
#     A^g_j = sum_i  J^g_i f(d_ij) / D^g_i ,   D^g_i = sum_k W^g_k f(d_ik)
# where J^g_i are band-g jobs at tract i (WAC CE01/CE03), W^g_k band-g resident
# workers (RAC CE01/CE03), f(d) = exp(-0.10 d_km) -- the job-access impedance
# used throughout this project (68/75), deliberately NOT the segregation decay.
# A^g_j is expected band-g jobs available per band-g worker at j after
# accounting for competing workers; its worker-weighted regional mean equals
# the regional band-g jobs-per-worker ratio, so values are directly readable
# as matched jobs per worker. Euclidean-gravity opportunity set, as elsewhere;
# the SLD transit measure counts total jobs only and cannot be split by
# earnings, so no transit variant is possible with these data (stated).
#
# Questions:
#   (1) descriptives: matched access for low- vs high-earnings workers by
#       exclusionary tercile, and the low/high ratio;
#   (2) does exclusionary zoning predict poor LOW-WAGE matched access even
#       conditional on OVERALL gravity job access (z_log_jobs_grav)? Naive and
#       overall-access-conditional rungs, county FE, jurisdiction clusters.
#
# Output: output/models/p4_matched_access_tercile.csv
#         output/models/p4_matched_access_models.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

## ---- frames ------------------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code, X_km, Y_km) |>
  filter(CBSA_Code == CO_CBSA_MAIN)
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR)))
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR)))

# regional computation frame: ALL Denver-MSA tracts with centroids (jobs and
# competing workers outside the zoning frame still shape access within it)
reg <- cent |>
  left_join(rac |> transmute(tract_id, w_low = CE01, w_high = CE03,
                             w_all = C000), by = "tract_id") |>
  left_join(wac |> transmute(tract_id, j_low = CE01, j_high = CE03,
                             j_all = C000), by = "tract_id") |>
  mutate(across(c(w_low, w_high, w_all, j_low, j_high, j_all),
                ~ replace_na(.x, 0)))
D <- as.matrix(dist(reg[, c("X_km", "Y_km")]))
FD <- exp(-0.10 * D)                       # includes diagonal f(0) = 1

shen <- function(J, W) {
  demand <- as.numeric(FD %*% W)           # D_i: competing workers around i
  as.numeric(FD %*% (J / pmax(demand, 1e-9)))
}
reg$acc_low  <- shen(reg$j_low,  reg$w_low)
reg$acc_high <- shen(reg$j_high, reg$w_high)
reg$acc_all  <- shen(reg$j_all,  reg$w_all)

## ---- descriptives by tercile --------------------------------------------------
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc |> select(tract_id, log_jobs_grav), by = "tract_id") |>
  left_join(reg |> select(tract_id, acc_low, acc_high, w_low),
            by = "tract_id") |>
  mutate(tercile = ntile(pct_reslow_of_res, 3),
         match_ratio = acc_low / acc_high)
by_terc <- xs |>
  group_by(tercile) |>
  summarise(n = n(),
            matched_access_low  = weighted.mean(acc_low, pmax(w_low, 1),
                                                na.rm = TRUE),
            matched_access_high = mean(acc_high, na.rm = TRUE),
            low_over_high = weighted.mean(match_ratio, pmax(w_low, 1),
                                          na.rm = TRUE), .groups = "drop")
write.csv(by_terc, file.path(DIR_CO_MOD, "p4_matched_access_tercile.csv"),
          row.names = FALSE)
print(as.data.frame(by_terc), digits = 3)

## ---- models -------------------------------------------------------------------
ma <- xs |>
  mutate(across(c(acc_low, acc_high, match_ratio, pct_reslow_of_res,
                  pct_black_rac, pct_lowincome_rac, log_worker_density_rac,
                  income_percapita_k, income_percapita_k_sq, log_jobs_grav),
                zscore, .names = "z_{.col}"))
COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}
res <- list()
for (y in c("acc_low", "acc_high", "match_ratio")) {
  res[[length(res) + 1]] <- tidy1(feols(as.formula(sprintf(
    "z_%s ~ z_pct_reslow_of_res + %s | county_fips", y, COVS)),
    data = ma, cluster = ~jurisd_main),
    sprintf("MA_%s_naive", y), "z_pct_reslow_of_res")
  res[[length(res) + 1]] <- tidy1(feols(as.formula(sprintf(
    "z_%s ~ z_pct_reslow_of_res + z_log_jobs_grav + %s | county_fips",
    y, COVS)), data = ma, cluster = ~jurisd_main),
    sprintf("MA_%s_cond_totalaccess", y),
    c("z_pct_reslow_of_res", "z_log_jobs_grav"))
}
out <- bind_rows(res)
write.csv(out, file.path(DIR_CO_MOD, "p4_matched_access_models.csv"),
          row.names = FALSE)
print(as.data.frame(out), digits = 3)
message("83 complete.")
