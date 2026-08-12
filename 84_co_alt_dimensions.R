# ==============================================================================
# 84_co_alt_dimensions.R      [PAPER 4, step 24 -- SI: alternative dimensions]
# Does the paper's headline decomposition replicate across segregation
# dimensions? The naive exclusionary moderation is nearly identical for
# White-Black (+0.333), Hispanic-non-Hispanic (+0.306), and college vs
# less-than-high-school (+0.345) segregation (64's ZB block). This script runs
# the missing piece: the PREFERRED, position-moderated rung for each
# dimension, so the SI can report naive vs preferred side by side.
#
# Dimensions (all already in the corrected Paper 2 panel; conventions
# asserted in 60): whiteblack = CR01/CR02; hisp_nonhisp = CT02/CT01;
# college_lesshs = CD education margins (workers 30+). White-non-White is
# deliberately NOT computed: it is not in the corrected measure set and the
# residual non-White aggregate pools heterogeneous groups.
#
# Output: output/models/p4_alt_dimensions.csv   (SI Table S14)
# ==============================================================================

source("60_co_setup.R")
library(fixest)

dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }

DIMS <- c("whiteblack", "hisp_nonhisp", "college_lesshs")
need <- c(sprintf("d_%s_rac_half", DIMS), sprintf("wexp_%s_wac_half", DIMS))
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id")
stopifnot(all(need %in% names(xs)))
xs <- xs |>
  mutate(across(all_of(c(need, "pct_reslow_of_res",
                         "pct_black_rac", "pct_lowincome_rac",
                         "log_worker_density_rac", "income_percapita_k",
                         "income_percapita_k_sq",
                         "dist_cbd_km", "dist_empctr_km")),
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
for (ms in DIMS) {
  y <- sprintf("z_wexp_%s_wac_half", ms)
  x <- sprintf("z_d_%s_rac_half", ms)
  keep <- sprintf("%s:z_pct_reslow_of_res", x)
  f_naive <- sprintf(
    "%s ~ %s * z_pct_reslow_of_res + %s | county_fips", y, x, COVS)
  f_pos <- sprintf(paste0(
    "%s ~ %s * z_pct_reslow_of_res + %s + z_dist_cbd_km + z_dist_empctr_km",
    " + %s:z_dist_cbd_km + %s:z_dist_empctr_km | county_fips"),
    y, x, COVS, x, x)
  # naive rung (replicates 64's ZB_mod for reference in the same file)
  res[[length(res) + 1]] <- tidy1(
    feols(as.formula(f_naive), data = xs, cluster = ~jurisd_main),
    sprintf("ALT_%s_naive", ms), keep)
  # preferred rung: metropolitan position moderates too (68's rung B)
  res[[length(res) + 1]] <- tidy1(
    feols(as.formula(f_pos), data = xs, cluster = ~jurisd_main),
    sprintf("ALT_%s_position", ms), keep)
}
out <- bind_rows(res) |>
  mutate(pct_of_naive = NA_real_)
for (ms in DIMS) {
  nv <- out$estimate[out$model_id == sprintf("ALT_%s_naive", ms)]
  out$pct_of_naive[out$model_id == sprintf("ALT_%s_position", ms)] <-
    100 * out$estimate[out$model_id == sprintf("ALT_%s_position", ms)] / nv
}
write.csv(out, file.path(DIR_CO_MOD, "p4_alt_dimensions.csv"),
          row.names = FALSE)
print(as.data.frame(out), digits = 3)
message("84 complete.")
